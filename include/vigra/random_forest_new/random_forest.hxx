#ifndef VIGRA_RF_RANDOM_FOREST_HXX
#define VIGRA_RF_RANDOM_FOREST_HXX

#include <type_traits>
#include <thread>

#include "../multi_shape.hxx"
#include "../binary_forest.hxx"
#include "../threadpool.hxx"



namespace vigra
{




template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACCTYPE, ContainerTag CTag = MapTag>
class RandomForest
{
public:

    typedef FEATURES Features;
    typedef typename Features::value_type FeatureType;
    typedef LABELS Labels;
    typedef typename Labels::value_type LabelType;
    typedef SPLITTESTS SplitTests;
    typedef ACCTYPE ACC;
    typedef typename ACC::input_type AccInputType;
    typedef BinaryForest Graph;
    typedef Graph::Node Node;
    typedef std::vector<size_t> DistributionType;
    
    static ContainerTag const container_tag = CTag;

    // Default (empty) constructor.
    RandomForest();

    // Default constructor (copy all of the given stuff).
    RandomForest(
        Graph const & graph,
        PropertyMap<Node, SplitTests, CTag> const & split_tests,
        PropertyMap<Node, AccInputType, CTag> const & node_responses,
        std::vector<LabelType> const & distinct_labels,
        size_t num_features
    );

    /// \brief Grow this forest by incorporating the other.
    void merge(
        RandomForest const & other
    );

    /// \brief Predict the given data and return the average number of split comparisons.
    /// \note labels should have the shape (features.shape()[0],).
    double predict(
        FEATURES const & features,
        LABELS & labels,
        int n_threads = -1
    ) const;

    /// \brief Predict the probabilities of the given data and return the average number of split comparisons.
    /// \note probs should have the shape (features.shape()[0], num_trees).
    template <typename PROBS>
    double predict_proba(
        FEATURES const & features,
        PROBS & probs,
        int n_threads = -1
    ) const;

    /// \brief For each data point in features, compute the corresponding leaf ids and return the average number of split comparisons.
    /// \note ids should have the shape (features.shape()[0], num_trees).
    template <typename IDS>
    double leaf_ids(
        FEATURES const & features,
        IDS & ids,
        int n_threads = -1
    ) const;

    /// \brief Return the number of nodes.
    size_t num_nodes() const
    {
        return graph_.numNodes();
    }

    /// \brief Return the number of trees.
    size_t num_trees() const
    {
        return graph_.numRoots();
    }

    /// \brief The graph structure.
    Graph graph_;

    /// \brief Contains a test for each internal node, that is used to determine whether given data goes to the left or the right child.
    PropertyMap<Node, SplitTests, CTag> split_tests_;

    /// \brief Contains the responses of each node (for example the most frequent label).
    PropertyMap<Node, AccInputType, CTag> node_responses_;

    /// \brief The distinct labels that were found in training.
    std::vector<LabelType> distinct_labels_;

    /// \brief The number of features.
    size_t num_features_;

private:

    /// \brief Compute the leaf ids of the instances in [from, to).
    template <typename IDS>
    double leaf_ids_impl(
        FEATURES const & features,
        IDS & ids,
        size_t from,
        size_t to
    ) const;

};

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::RandomForest()
    :
    graph_(),
    split_tests_(),
    node_responses_(),
    distinct_labels_(),
    num_features_(0)
{}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::RandomForest(
    Graph const & graph,
    PropertyMap<Node, SplitTests, CTag> const & split_tests,
    PropertyMap<Node, AccInputType, CTag> const & node_responses,
    std::vector<LabelType> const & distinct_labels,
    size_t num_features
)   :
    graph_(graph),
    split_tests_(split_tests),
    node_responses_(node_responses),
    distinct_labels_(distinct_labels),
    num_features_(num_features)
{}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
void RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::merge(
    RandomForest const & other
){
    vigra_precondition(num_features_ == other.num_features_,
                       "RandomForest::merge(): Number of features must not be different.");
    vigra_precondition(distinct_labels_ == other.distinct_labels_,
                       "RandomForest::merge(): The distinct labels must not be different.");

    size_t const offset = num_nodes();
    graph_.merge(other.graph_);
    for (auto const & p : other.split_tests_)
    {
        split_tests_.insert(Node(p.first.id()+offset), p.second);
    }
    for (auto const & p : other.node_responses_)
    {
        node_responses_.insert(Node(p.first.id()+offset), p.second);
    }
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::predict(
    FEATURES const & features,
    LABELS & labels,
    int n_threads
) const {
    vigra_precondition(features.shape()[0] == labels.shape()[0],
                       "RandomForest::predict(): Shape mismatch between features and labels.");
    vigra_precondition(features.shape()[1] == num_features_,
                       "RandomForest::predict(): Number of features in prediction differs from training.");

    MultiArray<2, double> probs(Shape2(features.shape()[0], distinct_labels_.size()));
    double const average_split_counts = predict_proba(features, probs, n_threads);
    for (size_t i = 0; i < features.shape()[0]; ++i)
    {
        auto const sub_probs = probs.template bind<0>(i);
        auto it = std::max_element(sub_probs.begin(), sub_probs.end());
        size_t const label = std::distance(sub_probs.begin(), it);
        labels(i) = distinct_labels_[label];
    }
    return average_split_counts;
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
template <typename PROBS>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::predict_proba(
    FEATURES const & features,
    PROBS & probs,
    int n_threads
) const {
    vigra_precondition(features.shape()[0] == probs.shape()[0],
                       "RandomForest::predict_proba(): Shape mismatch between features and probabilities.");
    vigra_precondition(features.shape()[1] == num_features_,
                       "RandomForest::predict_proba(): Number of features in prediction differs from training.");
    vigra_precondition(probs.shape()[1] == distinct_labels_.size(),
                       "RandomForest::predict_proba(): Number of labels in probabilities differs from training.");

    size_t const num_roots = graph_.numRoots();
    MultiArray<2, size_t> ids(Shape2(features.shape()[0], num_roots));
    double const average_split_counts = leaf_ids(features, ids, n_threads);
    ACC acc;
    for (size_t i = 0; i < features.shape()[0]; ++i)
    {
        std::vector<AccInputType> tree_results(num_roots);
        for (size_t k = 0; k < num_roots; ++k)
        {
            tree_results[k] = node_responses_.at(Node(ids(i, k)));
        }
        auto sub_probs = probs.template bind<0>(i);
        acc(tree_results.begin(), tree_results.end(), sub_probs.begin());
    }
    return average_split_counts;
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
template <typename IDS>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::leaf_ids(
    FEATURES const & features,
    IDS & ids,
    int n_threads
) const {
    vigra_precondition(features.shape()[0] == ids.shape()[0],
                       "RandomForest::leaf_ids(): Shape mismatch between features and probabilities.");
    vigra_precondition(features.shape()[1] == num_features_,
                       "RandomForest::leaf_ids(): Number of features in prediction differs from training.");
    vigra_precondition(ids.shape()[1] == graph_.numRoots(),
                       "RandomForest::leaf_ids(): Leaf array has wrong shape.");

    size_t const num_instances = features.shape()[0];
    if (n_threads == -1)
        n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1)
        n_threads = 1;
    std::vector<double> split_comparisons(n_threads, 0.0);
    std::vector<size_t> indices(num_instances);
    std::iota(indices.begin(), indices.end(), 0);
    parallel_foreach(
        n_threads,
        num_instances,
        indices.begin(),
        indices.end(),
        [this, &features, &ids, &split_comparisons](size_t thread_id, size_t i) {
            split_comparisons[thread_id] += leaf_ids_impl(features, ids, i, i+1);
        }
    );

    double const sum_split_comparisons = std::accumulate(split_comparisons.begin(), split_comparisons.end(), 0.0);
    return sum_split_comparisons / features.shape()[0];
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC, ContainerTag CTag>
template <typename IDS>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC, CTag>::leaf_ids_impl(
    FEATURES const & features,
    IDS & ids,
    size_t from,
    size_t to
) const {
    vigra_precondition(features.shape()[0] == ids.shape()[0],
                       "RandomForest::leaf_ids_impl(): Shape mismatch between features and labels.");
    vigra_precondition(features.shape()[1] == num_features_,
                       "RandomForest::leaf_ids_impl(): Number of Features in prediction differs from training.");
    vigra_precondition(from >= 0 && from <= to && to <= features.shape()[0],
                       "RandomForest::leaf_ids_impl(): Indices out of range.");
    vigra_precondition(ids.shape()[1] == graph_.numRoots(),
                       "RandomForest::leaf_ids_impl(): Leaf array has wrong shape.");

    double split_comparisons = 0.0;
    for (size_t i = from; i < to; ++i)
    {
        auto const sub_features = features.template bind<0>(i);
        for (size_t k = 0; k < graph_.numRoots(); ++k)
        {
            Node node = graph_.getRoot(k);
            while (graph_.outDegree(node) > 0)
            {
                size_t const child_index = split_tests_.at(node)(sub_features);
                node = graph_.getChild(node, child_index);
                split_comparisons += 1.0;
            }
            ids(i, k) = node.id();
        }
    }
    return split_comparisons;
}



} // namespace vigra

#endif
