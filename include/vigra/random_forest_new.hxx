#ifndef VIGRA_RANDOM_FOREST_NEW_HXX
#define VIGRA_RANDOM_FOREST_NEW_HXX



#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>

#include "multi_array.hxx"
#include "sampling.hxx"
#include "threadpool.hxx"
#include "random_forest_new/random_forest.hxx"
#include "random_forest_new/random_forest_common.hxx"



namespace vigra
{



namespace detail
{
    template <typename ACC>
    struct RFMapUpdater
    {
        template <typename A, typename B>
        void operator()(A & a, B const & b) const
        {
            a = b;
        }
    };

    template <>
    struct RFMapUpdater<ArgMaxAcc>
    {
        template <typename A, typename B>
        void operator()(A & a, B const & b) const
        {
            auto it = std::max_element(b.begin(), b.end());
            a = std::distance(b.begin(), it);
        }
    };

    template <typename FEATURES, typename LABELS, typename SAMPLER, typename SCORER>
    void split_score(
            FEATURES const & features,
            LABELS const & labels,
            std::vector<size_t> const & instance_weights,
            std::vector<size_t> const & instances,
            SAMPLER const & dim_sampler,
            SCORER & score
    ){
        typedef typename FEATURES::value_type FeatureType;

        auto feats = std::vector<FeatureType>(instances.size());
        auto sorted_indices = std::vector<size_t>(feats.size());

        for (size_t i = 0; i < dim_sampler.sampleSize(); ++i)
        {
            size_t const d = dim_sampler[i];

            // Copy the features to a vector with the correct size (so the sort is faster because of data locality).
            for (size_t kk = 0; kk < instances.size(); ++kk)
                feats[kk] = features(instances[kk], d);

            // Sort the features.
            indexSort(feats.begin(), feats.end(), sorted_indices.begin());
            auto tosort_instances = std::vector<size_t>(instances);
            applyPermutation(sorted_indices.begin(), sorted_indices.end(), instances.begin(), tosort_instances.begin());

            // Get the score of the splits.
            score(features, labels, instance_weights, tosort_instances.begin(), tosort_instances.end(), d);
        }
    }
}



template <typename FEATURES, typename LABELS>
struct DefaultRF
{
    typedef RandomForest<FEATURES,
                         LABELS,
                         LessEqualSplitTest<typename FEATURES::value_type>,
                         ArgMaxVectorAcc<size_t>,
                         VectorTag> type;
};



/**
 * @brief Train a single randomized decision tree.
 */
template <typename RF, typename SCORER, typename STOP, typename RANDENGINE>
void random_forest_single_tree(
        typename RF::Features const & features,
        typename RF::Labels const & labels,
        RandomForestNewOptions const & options,
        STOP stop,
        RF & tree,
        RANDENGINE const & randengine
){
    typedef typename RF::Features Features;
    typedef typename Features::value_type FeatureType;
    typedef LessEqualSplitTest<FeatureType> SplitTests;
    typedef typename RF::Node Node;
    typedef typename RF::ACC ACC;
    typedef typename ACC::input_type ACCInputType;

    static_assert(std::is_same<SplitTests, typename RF::SplitTests>::value,
                  "random_forest_single_tree(): Wrong Random Forest class.");

    auto const num_instances = features.shape()[0];
    auto const num_features = features.shape()[1];
    vigra_precondition(num_instances == labels.size(),
                       "random_forest_single_tree(): Shape mismatch between features and labels.");
    vigra_precondition(num_features == tree.problem_spec_.num_features_,
                       "random_forest_single_tree(): Wrong number of features.");

    // Create the index vector for bookkeeping.
    auto instance_indices = std::vector<size_t>(num_instances);
    std::iota(instance_indices.begin(), instance_indices.end(), 0);
    typedef std::vector<size_t>::iterator InstanceIter;

    // Create the weights for the bootstrap sample.
    auto instance_weights = std::vector<size_t>(num_instances, 1);
    if (options.bootstrap_sampling_)
    {
        std::fill(instance_weights.begin(), instance_weights.end(), 0);
        Sampler<MersenneTwister> sampler(num_instances,
                                         SamplerOptions().withReplacement().stratified(options.use_stratification_),
                                         &randengine);
        sampler.sample();
        for (size_t i = 0; i < sampler.sampleSize(); ++i)
        {
            int const index = sampler[i];
            ++instance_weights[index];
        }
    }

    // Create the sampler for the split dimensions.
    auto const mtry = tree.problem_spec_.actual_mtry_;
    Sampler<MersenneTwister> dim_sampler(num_features, SamplerOptions().withoutReplacement().sampleSize(mtry), &randengine);

    // Create the node stack and place the root node inside.
    std::stack<Node> node_stack;
    typedef std::pair<InstanceIter, InstanceIter> IterPair;
    PropertyMap<Node, IterPair> instance_range;
    PropertyMap<Node, std::vector<size_t> > node_distributions;
    PropertyMap<Node, size_t> node_depths;
    {
        auto const rootnode = tree.graph_.addNode();
        node_stack.push(rootnode);
        instance_range.insert(rootnode, IterPair(instance_indices.begin(), instance_indices.end()));

        auto priors = std::vector<size_t>(tree.problem_spec_.num_classes_, 0);
        for (auto i : instance_indices)
            priors[labels(i)] += instance_weights[i];

        node_distributions.insert(rootnode, priors);
        node_depths.insert(rootnode, 0);
    }

    // Split the nodes.
    detail::RFMapUpdater<ACC> node_map_updater;
    while (!node_stack.empty())
    {
        auto const node = node_stack.top();
        node_stack.pop();
        auto const begin = instance_range.at(node).first;
        auto const end = instance_range.at(node).second;
        auto const & priors = node_distributions.at(node);
        auto const depth = node_depths.at(node);

        // Get the instances with weight > 0.
        std::vector<size_t> used_instances;
        for (auto it = begin; it != end; ++it)
            if (instance_weights[*it] > 0)
                used_instances.push_back(*it);
 
        // Find the best split.
        dim_sampler.sample();
        SCORER score(priors);
        if (options.resample_count_ == 0 || used_instances.size() <= options.resample_count_)
        {
            detail::split_score(
                features,
                labels,
                instance_weights,
                used_instances,
                dim_sampler,
                score
            );
        }
        else
        {
            // Generate a random subset of the node instances.
            Sampler<MersenneTwister> resampler(used_instances.begin(), used_instances.end(), SamplerOptions().withoutReplacement().sampleSize(options.resample_count_), &randengine);
            resampler.sample();
            auto indices = std::vector<size_t>(options.resample_count_);
            for (size_t i = 0; i < options.resample_count_; ++i)
                indices[i] = used_instances[resampler[i]];

            detail::split_score(
                features,
                labels,
                instance_weights,
                indices,
                dim_sampler,
                score
            );
        }

        // If no split was found, the node is terminal.
        if (!score.split_found_)
        {
            tree.node_responses_.insert(node, ACCInputType());
            node_map_updater(tree.node_responses_.at(node), node_distributions.at(node));
            continue;
        }

        // Do the split.
        auto const n_left = tree.graph_.addNode();
        auto const n_right = tree.graph_.addNode();
        tree.graph_.addArc(node, n_left);
        tree.graph_.addArc(node, n_right);
        auto const best_split = score.best_split_;
        auto const best_dim = score.best_dim_;
        auto const split_iter = std::partition(begin, end,
            [&](size_t i)
            {
                return features(i, best_dim) <= best_split;
            }
        );
        instance_range.insert(n_left, IterPair(begin, split_iter));
        instance_range.insert(n_right, IterPair(split_iter, end));
        tree.split_tests_.insert(node, SplitTests(best_dim, best_split));
        node_depths.insert(n_left, depth+1);
        node_depths.insert(n_right, depth+1);

        // Update the prior list for the left child and check if the node is terminal.
        auto priors_left = std::vector<size_t>(tree.problem_spec_.num_classes_, 0);
        for (auto it = begin; it != split_iter; ++it)
            priors_left[labels(*it)] += instance_weights[*it];
        node_distributions.insert(n_left, priors_left);

        // Check if the left child is terminal.
        if (stop(labels, RFNodeDescription<decltype(priors_left)>(depth+1, priors_left)))
        {
            tree.node_responses_.insert(n_left, ACCInputType());
            node_map_updater(tree.node_responses_.at(n_left), node_distributions.at(n_left));
        }
        else
        {
            node_stack.push(n_left);
        }

        // Update the prior list for the right child.
        auto priors_right = std::vector<size_t>(tree.problem_spec_.num_classes_, 0);
        for (auto it = split_iter; it != end; ++it)
            priors_right[labels(*it)] += instance_weights[*it];
        node_distributions.insert(n_right, priors_right);

        // Check if the right child is terminal.
        if (stop(labels, RFNodeDescription<decltype(priors_right)>(depth+1, priors_right)))
        {
            tree.node_responses_.insert(n_right, ACCInputType());
            node_map_updater(tree.node_responses_.at(n_right), node_distributions.at(n_right));
        }
        else
        {
            node_stack.push(n_right);
        }
    }
}



/// \brief Preprocess the labels and call the train functions on the single trees.
template <typename FEATURES,
          typename LABELS,
          typename SCORER,
          typename STOP>
typename DefaultRF<FEATURES, LABELS>::type
random_forest_impl(
        FEATURES const & features,
        LABELS const & labels,
        RandomForestNewOptions const & options,
        int n_threads,
        STOP const & stop
){
    // typedef FEATURES Features;
    typedef LABELS Labels;
    // typedef typename Features::value_type FeatureType;
    typedef typename Labels::value_type LabelType;
    typedef typename DefaultRF<FEATURES, LABELS>::type RF;

    // Check the number of trees.
    size_t const tree_count = options.tree_count_;
    vigra_precondition(tree_count > 0, "random_forest(): tree_count must not be zero.");
    std::vector<RF> trees(tree_count);

    // Transform the labels to 0, 1, 2, ...
    std::set<LabelType> const dlabels(labels.begin(), labels.end());
    std::vector<LabelType> const distinct_labels(dlabels.begin(), dlabels.end());
    std::map<LabelType, size_t> label_map;
    for (size_t i = 0; i < distinct_labels.size(); ++i)
    {
        label_map[distinct_labels[i]] = i;
    }
    MultiArray<1, size_t> transformed_labels(Shape1(labels.size()));
    for (size_t i = 0; i < labels.size(); ++i)
    {
        transformed_labels(i) = label_map[labels(i)];
    }

    // Write the problem specification into the trees.
    auto const pspec = ProblemSpecNew<LabelType>()
                           .num_features(features.shape()[1])
                           .num_instances(features.shape()[0])
                           .distinct_classes(distinct_labels)
                           .actual_mtry(options.get_features_per_node(features.shape()[1]));
    for (auto & t : trees)
        t.problem_spec_ = pspec;

    // Find the correct number of threads.
    if (n_threads == -1)
        n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1)
        n_threads = 1;

    // Use the global random engine to create seeds for the random engines that run in the threads.
    UniformIntRandomFunctor<MersenneTwister> rand_functor;
    std::set<UInt32> seeds;
    while (seeds.size() < n_threads)
    {
        seeds.insert(rand_functor());
    }
    vigra_assert(seeds.size() == n_threads, "random_forest(): Could not create random seeds.");

    // Create the random engines that run in the threads.
    std::vector<MersenneTwister> rand_engines;
    for (auto seed : seeds)
    {
        rand_engines.push_back(MersenneTwister(seed));
    }

    // Train the trees.
    ThreadPool pool((size_t)n_threads);
    for (size_t i = 0; i < tree_count; ++i)
    {
        pool.enqueue([&features, &transformed_labels, &options, &stop, &trees, i, &rand_engines](size_t thread_id)
            {
                random_forest_single_tree<RF, SCORER, STOP>(features, transformed_labels, options, stop, trees[i], rand_engines[thread_id]);
            }
        );
    }
    pool.waitFinished();

    // Merge the trees together and return the result.
    RF rf = trees[0];
    for (size_t i = 1; i < trees.size(); ++i)
    {
        rf.merge(trees[i]);
    }
    return rf;
}



/// \brief Get the stop criterion from the option object and pass it as template argument.
template <typename FEATURES, typename LABELS, typename SCORER>
typename DefaultRF<FEATURES, LABELS>::type random_forest_impl0(
        FEATURES const & features,
        LABELS const & labels,
        RandomForestNewOptions const & options,
        int n_threads
){
    if (options.max_depth_ > 0)
        return random_forest_impl<FEATURES, LABELS, SCORER, DepthStop>(features, labels, options, n_threads, DepthStop(options.max_depth_));
    else if (options.min_num_instances_ > 1)
        return random_forest_impl<FEATURES, LABELS, SCORER, NumInstancesStop>(features, labels, options, n_threads, NumInstancesStop(options.min_num_instances_));
    else if (options.node_complexity_tau_ > 0)
        return random_forest_impl<FEATURES, LABELS, SCORER, NodeComplexityStop>(features, labels, options, n_threads, NodeComplexityStop(options.node_complexity_tau_));
    else
        return random_forest_impl<FEATURES, LABELS, SCORER, PurityStop>(features, labels, options, n_threads, PurityStop());
}



/// \brief Get the scorer from the option object and pass it as template argument.
template <typename FEATURES, typename LABELS>
typename DefaultRF<FEATURES, LABELS>::type random_forest(
        FEATURES const & features,
        LABELS const & labels,
        RandomForestNewOptions const & options = RandomForestNewOptions(),
        int n_threads = -1
){
    if (options.split_ == RF_GINI)
        return random_forest_impl0<FEATURES, LABELS, GiniScorer>(features, labels, options, n_threads);
    else if (options.split_ == RF_ENTROPY)
        return random_forest_impl0<FEATURES, LABELS, EntropyScorer>(features, labels, options, n_threads);
    else if (options.split_ == RF_KSD)
        return random_forest_impl0<FEATURES, LABELS, KSDScorer>(features, labels, options, n_threads);
    else
        throw std::runtime_error("random_forest(): Unknown split.");
}



} // namespace vigra

#endif
