/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/
#ifndef VIGRA_RF3_RANDOM_FOREST_HXX
#define VIGRA_RF3_RANDOM_FOREST_HXX

#include <type_traits>
#include <thread>

#include "../multi_shape.hxx"
#include "../binary_forest.hxx"
#include "../threadpool.hxx"
#include "random_forest_common.hxx"



namespace vigra
{

namespace rf3
{

/********************************************************/
/*                                                      */
/*                    rf3::RandomForest                 */
/*                                                      */
/********************************************************/

/** \brief Random forest version 3.

    vigra::rf3::RandomForest is typicall constructed via the factory function \ref vigra::rf3::random_forest().
*/
template <typename FEATURES,
          typename LABELS,
          typename SPLITTESTS = LessEqualSplitTest<typename FEATURES::value_type>,
          typename ACCTYPE = ArgMaxVectorAcc<double>>
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

    static ContainerTag const container_tag = VectorTag;

    // FIXME:
    // Once the support for Visual Studio 2012 is dropped, replace this struct with
    // template <typename T>
    // using NodeMap = PropertyMap<Node, T, container_tag>;
    // Then the verbose typename NodeMap<T>::type, which typically shows up on NodeMap usages,
    // can be replace with NodeMap<T>.
    template <typename T>
    struct NodeMap
    {
        typedef PropertyMap<Node, T, container_tag> type;
    };

    // Default (empty) constructor.
    RandomForest();

    // Default constructor (copy all of the given stuff).
    RandomForest(
        Graph const & graph,
        typename NodeMap<SplitTests>::type const & split_tests,
        typename NodeMap<AccInputType>::type const & node_responses,
        ProblemSpec<LabelType> const & problem_spec
    );

    /// \brief Grow this forest by incorporating the other.
    void merge(
        RandomForest const & other
    );

    /// \brief Predict the given data and return the average number of split comparisons.
    /// \note labels must be a 1-D array with size <tt>features.shape(0)</tt>.
    void predict(
        FEATURES const & features,
        LABELS & labels,
        int n_threads = -1,
        const std::vector<size_t> & tree_indices = std::vector<size_t>()
    ) const;

    /// \brief Predict the probabilities of the given data and return the average number of split comparisons.
    /// \note probs should have the shape (features.shape()[0], num_classes).
    template <typename PROBS>
    void predict_probabilities(
        FEATURES const & features,
        PROBS & probs,
        int n_threads = -1,
        const std::vector<size_t> & tree_indices = std::vector<size_t>()
    ) const;

    /// \brief For each data point in features, compute the corresponding leaf ids and return the average number of split comparisons.
    /// \note ids should have the shape (features.shape()[0], num_trees).
    template <typename IDS>
    double leaf_ids(
        FEATURES const & features,
        IDS & ids,
        int n_threads = -1,
        const std::vector<size_t> tree_indices = std::vector<size_t>()
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

    /// \brief Return the number of classes.
    size_t num_classes() const
    {
        return problem_spec_.num_classes_;
    }
    
    /// \brief Return the number of classes.
    size_t num_features() const
    {
        return problem_spec_.num_features_;
    }

    /// \brief The graph structure.
    Graph graph_;

    /// \brief Contains a test for each internal node, that is used to determine whether given data goes to the left or the right child.
    typename NodeMap<SplitTests>::type split_tests_;

    /// \brief Contains the responses of each node (for example the most frequent label).
    typename NodeMap<AccInputType>::type node_responses_;

    /// \brief The specifications.
    ProblemSpec<LabelType> problem_spec_;

    /// \brief The options that were used for training.
    RandomForestOptions options_;

private:

    /// \brief Compute the leaf ids of the instances in [from, to).
    template <typename IDS, typename INDICES>
    double leaf_ids_impl(
        FEATURES const & features,
        IDS & ids,
        size_t from,
        size_t to,
        INDICES const & tree_indices
    ) const;

    template<typename PROBS>
    void predict_probabilities_impl(
        FEATURES const & features,
        PROBS & probs,
        const size_t i,
        const std::vector<size_t> & tree_indices) const;

};

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::RandomForest()
    :
    graph_(),
    split_tests_(),
    node_responses_(),
    problem_spec_()
{}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::RandomForest(
    Graph const & graph,
    typename NodeMap<SplitTests>::type const & split_tests,
    typename NodeMap<AccInputType>::type const & node_responses,
    ProblemSpec<LabelType> const & problem_spec
)   :
    graph_(graph),
    split_tests_(split_tests),
    node_responses_(node_responses),
    problem_spec_(problem_spec)
{}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
void RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::merge(
    RandomForest const & other
){
    vigra_precondition(problem_spec_ == other.problem_spec_,
                       "RandomForest::merge(): You cannot merge with different problem specs.");

    // FIXME: Eventually compare the options and only fix if the forests are compatible.

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

// FIXME TODO we don't support the selection of tree indices any more in predict_probabilities, might be a good idea
// to re-enable this.
template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
void RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::predict(
    FEATURES const & features,
    LABELS & labels,
    int n_threads,
    const std::vector<size_t> & tree_indices
) const {
    vigra_precondition(features.shape()[0] == labels.shape()[0],
                       "RandomForest::predict(): Shape mismatch between features and labels.");
    vigra_precondition((size_t)features.shape()[1] == problem_spec_.num_features_,
                       "RandomForest::predict(): Number of features in prediction differs from training.");

    MultiArray<2, double> probs(Shape2(features.shape()[0], problem_spec_.num_classes_));
    predict_probabilities(features, probs, n_threads, tree_indices);
    for (size_t i = 0; i < (size_t)features.shape()[0]; ++i)
    {
        auto const sub_probs = probs.template bind<0>(i);
        auto it = std::max_element(sub_probs.begin(), sub_probs.end());
        size_t const label = std::distance(sub_probs.begin(), it);
        labels(i) = problem_spec_.distinct_classes_[label];
    }
}


// FIXME TODO we don't support the selection of tree indices any more in predict_probabilities, might be a good idea
// to re-enable this.
template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
template <typename PROBS>
void RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::predict_probabilities(
    FEATURES const & features,
    PROBS & probs,
    int n_threads,
    const std::vector<size_t> & tree_indices
) const {
    vigra_precondition(features.shape()[0] == probs.shape()[0],
                       "RandomForest::predict_probabilities(): Shape mismatch between features and probabilities.");
    vigra_precondition((size_t)features.shape()[1] == problem_spec_.num_features_,
                       "RandomForest::predict_probabilities(): Number of features in prediction differs from training.");
    vigra_precondition((size_t)probs.shape()[1] == problem_spec_.num_classes_,
                       "RandomForest::predict_probabilities(): Number of labels in probabilities differs from training.");

    // By default, actual_tree_indices is empty. In that case we want to use all trees. 
    // We need to make a copy. I really don't know how the old code did compile...
    std::vector<size_t> tree_indices_cpy(tree_indices);
    if (tree_indices_cpy.size() == 0)
    {
        tree_indices_cpy.resize(graph_.numRoots());
        std::iota(tree_indices_cpy.begin(), tree_indices_cpy.end(), 0);
    }
    else {
        // Check the tree indices.
        std::sort(tree_indices_cpy.begin(), tree_indices_cpy.end());
        tree_indices_cpy.erase(std::unique(tree_indices_cpy.begin(), tree_indices_cpy.end()), tree_indices_cpy.end());
        for (auto i : tree_indices_cpy)
            vigra_precondition(i < graph_.numRoots(), "RandomForest::leaf_ids(): Tree index out of range.");
    }
    
    size_t const num_instances = features.shape()[0];
    
    if (n_threads == -1)
        n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1)
        n_threads = 1;
    
    parallel_foreach(
        n_threads,
        num_instances,
        [&features,&probs,&tree_indices_cpy,this](size_t, size_t i) {
            this->predict_probabilities_impl(features, probs, i, tree_indices_cpy);
        }
    );
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
template <typename PROBS>
void RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::predict_probabilities_impl(
    FEATURES const & features,
    PROBS & probs,
    const size_t i,
    const std::vector<size_t> & tree_indices
) const {

    // instantiate the accumulation function and the vector to store the tree node results
    ACC acc;
    std::vector<AccInputType> tree_results;
    tree_results.reserve(tree_indices.size());
    auto const sub_features = features.template bind<0>(i);
    
    // loop over the trees
    for (auto k : tree_indices)
    {
        Node node = graph_.getRoot(k);
        while (graph_.outDegree(node) > 0)
        {
            size_t const child_index = split_tests_.at(node)(sub_features);
            node = graph_.getChild(node, child_index);
        }
        tree_results.emplace_back(node_responses_.at(node));
    }

    // write the tree results into the probabilities
    auto sub_probs = probs.template bind<0>(i);
    acc(tree_results.begin(), tree_results.end(), sub_probs.begin());    
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
template <typename IDS>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::leaf_ids(
    FEATURES const & features,
    IDS & ids,
    int n_threads,
    std::vector<size_t> tree_indices
) const {
    vigra_precondition(features.shape()[0] == ids.shape()[0],
                       "RandomForest::leaf_ids(): Shape mismatch between features and probabilities.");
    vigra_precondition((size_t)features.shape()[1] == problem_spec_.num_features_,
                       "RandomForest::leaf_ids(): Number of features in prediction differs from training.");
    vigra_precondition(ids.shape()[1] == graph_.numRoots(),
                       "RandomForest::leaf_ids(): Leaf array has wrong shape.");

    // Check the tree indices.
    std::sort(tree_indices.begin(), tree_indices.end());
    tree_indices.erase(std::unique(tree_indices.begin(), tree_indices.end()), tree_indices.end());
    for (auto i : tree_indices)
        vigra_precondition(i < graph_.numRoots(), "RandomForest::leaf_ids(): Tree index out of range.");

    // By default, actual_tree_indices is empty. In that case we want to use all trees.
    if (tree_indices.size() == 0)
    {
        tree_indices.resize(graph_.numRoots());
        std::iota(tree_indices.begin(), tree_indices.end(), 0);
    }

    size_t const num_instances = features.shape()[0];
    if (n_threads == -1)
        n_threads = std::thread::hardware_concurrency();
    if (n_threads < 1)
        n_threads = 1;
    std::vector<double> split_comparisons(n_threads, 0.0);
    std::vector<size_t> indices(num_instances);
    std::iota(indices.begin(), indices.end(), 0);
    std::fill(ids.begin(), ids.end(), -1);
    parallel_foreach(
        n_threads,
        indices.begin(),
        indices.end(),
        [this, &features, &ids, &split_comparisons, &tree_indices](size_t thread_id, size_t i) {
            split_comparisons[thread_id] += this->leaf_ids_impl(features, ids, i, i+1, tree_indices);
        }
    );

    double const sum_split_comparisons = std::accumulate(split_comparisons.begin(), split_comparisons.end(), 0.0);
    return sum_split_comparisons / features.shape()[0];
}

template <typename FEATURES, typename LABELS, typename SPLITTESTS, typename ACC>
template <typename IDS, typename INDICES>
double RandomForest<FEATURES, LABELS, SPLITTESTS, ACC>::leaf_ids_impl(
    FEATURES const & features,
    IDS & ids,
    size_t from,
    size_t to,
    INDICES const & tree_indices
) const {
    vigra_precondition(features.shape()[0] == ids.shape()[0],
                       "RandomForest::leaf_ids_impl(): Shape mismatch between features and labels.");
    vigra_precondition(features.shape()[1] == problem_spec_.num_features_,
                       "RandomForest::leaf_ids_impl(): Number of Features in prediction differs from training.");
    vigra_precondition(from >= 0 && from <= to && to <= (size_t)features.shape()[0],
                       "RandomForest::leaf_ids_impl(): Indices out of range.");
    vigra_precondition(ids.shape()[1] == graph_.numRoots(),
                       "RandomForest::leaf_ids_impl(): Leaf array has wrong shape.");

    double split_comparisons = 0.0;
    for (size_t i = from; i < to; ++i)
    {
        auto const sub_features = features.template bind<0>(i);
        for (auto k : tree_indices)
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



} // namespace rf3
} // namespace vigra

#endif
