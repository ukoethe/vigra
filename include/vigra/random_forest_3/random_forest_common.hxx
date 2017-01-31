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
#ifndef VIGRA_RF3_COMMON_HXX
#define VIGRA_RF3_COMMON_HXX

#include <iterator>
#include <type_traits>
#include <cmath>
#include <numeric>

#include "../multi_array.hxx"
#include "../mathutil.hxx"

namespace vigra
{

namespace rf3
{

/** \addtogroup MachineLearning
**/
//@{

template <typename T>
struct LessEqualSplitTest
{
public:
    LessEqualSplitTest(size_t dim = 0, T const & val = 0)
        :
        dim_(dim),
        val_(val)
    {}
    
    template<typename FEATURES>
    size_t operator()(FEATURES const & features) const
    {
        return features(dim_) <= val_ ? 0 : 1;
    }
    
    size_t dim_;
    T val_;
};



struct ArgMaxAcc
{
public:
    typedef size_t input_type;

    template <typename ITER, typename OUTITER>
    void operator()(ITER begin, ITER end, OUTITER out)
    {
        std::fill(buffer_.begin(), buffer_.end(), 0);
        size_t max_v = 0;
        size_t n = 0;
        for (ITER it = begin; it != end; ++it)
        {
            size_t const v = *it;
            if (v >= buffer_.size())
            {
                buffer_.resize(v+1, 0);
            }
            ++buffer_[v];
            ++n;
            max_v = std::max(max_v, v);
        }
        for (size_t i = 0; i <= max_v; ++i)
        {
            *out = buffer_[i] / static_cast<double>(n);
            ++out;
        }
    }
private:
    std::vector<size_t> buffer_;
};



template <typename VALUETYPE>
struct ArgMaxVectorAcc
{
public:
    typedef VALUETYPE value_type;
    typedef std::vector<value_type> input_type;
    template <typename ITER, typename OUTITER>
    void operator()(ITER begin, ITER end, OUTITER out)
    {
        std::fill(buffer_.begin(), buffer_.end(), 0);
        size_t max_v = 0;
        for (ITER it = begin; it != end; ++it)
        {
            input_type const & vec = *it;
            if (vec.size() >= buffer_.size())
            {
                buffer_.resize(vec.size(), 0);
            }
            value_type const n = std::accumulate(vec.begin(), vec.end(), static_cast<value_type>(0));
            for (size_t i = 0; i < vec.size(); ++i)
            {
                buffer_[i] += vec[i] / static_cast<double>(n);
            }
            max_v = std::max(vec.size()-1, max_v);
        }
        for (size_t i = 0; i <= max_v; ++i)
        {
            *out = buffer_[i];
            ++out;
        }
    }
    private:
        std::vector<double> buffer_;
};



// struct LargestSumAcc
// {
// public:
//     typedef std::vector<size_t> input_type;
//     template <typename ITER>
//     size_t operator()(ITER begin, ITER end)
//     {
//         std::fill(buffer_.begin(), buffer_.end(), 0);
//         for (ITER it = begin; it != end; ++it)
//         {
//             auto const & v = *it;
//             if (v.size() > buffer_.size())
//             {
//                 buffer_.resize(v.size(), 0);
//             }
//             for (size_t i = 0; i < v.size(); ++i)
//             {
//                 buffer_[i] += v[i];
//             }
//         }
//         size_t max_label = 0;
//         size_t max_count = 0;
//         for (size_t i = 0; i < buffer_.size(); ++i)
//         {
//             if (buffer_[i] > max_count)
//             {
//                 max_count = buffer_[i];
//                 max_label = i;
//             }
//         }
//         return max_label;
//     }
// private:
//     std::vector<size_t> buffer_;
// };



// struct ForestGarroteAcc
// {
// public:
//     typedef double input_type;
//     template <typename ITER, typename OUTITER>
//     void operator()(ITER begin, ITER end, OUTITER out)
//     {
//         double s = 0.0;
//         for (ITER it = begin; it != end; ++it)
//         {
//             s += *it;
//         }
//         if (s < 0.0)
//             s = 0.0;
//         else if (s > 1.0)
//             s = 1.0;
//         *out = 1.0-s;
//         ++out;
//         *out = s;
//     }
// };



namespace detail
{

    /// Abstract scorer that iterates over all split candidates, uses FUNCTOR to compute a score,
    /// and saves the split with the minimum score.
    template <typename FUNCTOR>
    class GeneralScorer
    {
    public:

        typedef FUNCTOR Functor;

        GeneralScorer(std::vector<double> const & priors)
            :
            split_found_(false),
            best_split_(0),
            best_dim_(0),
            best_score_(std::numeric_limits<double>::max()),
            priors_(priors),
            n_total_(std::accumulate(priors.begin(), priors.end(), 0.0))
        {}

        template <typename FEATURES, typename LABELS, typename WEIGHTS, typename ITER>
        void operator()(
            FEATURES const & features,
            LABELS const & labels,
            WEIGHTS const & weights,
            ITER begin,
            ITER end,
            size_t dim
        ){
            if (begin == end)
                return;

            Functor score;

            std::vector<double> counts(priors_.size(), 0.0);
            double n_left = 0;
            ITER next = begin;
            ++next;
            for (; next != end; ++begin, ++next)
            {
                // Move the label from the right side to the left side.
                size_t const left_index = *begin;
                size_t const right_index = *next;
                size_t const label = static_cast<size_t>(labels(left_index));
                counts[label] += weights[left_index];
                n_left += weights[left_index];

                // Skip if there is no new split.
                auto const left = features(left_index, dim);
                auto const right = features(right_index, dim);
                if (left == right)
                    continue;

                // Update the score.
                split_found_ = true;
                double const s = score(priors_, counts, n_total_, n_left);
                bool const better_score = s < best_score_;
                if (better_score)
                {
                    best_score_ = s;
                    best_split_ = 0.5*(left+right);
                    best_dim_ = dim;
                }
            }
        }

        bool split_found_; // whether a split was found at all
        double best_split_; // the threshold of the best split
        size_t best_dim_; // the dimension of the best split
        double best_score_; // the score of the best split

    private:

        std::vector<double> const priors_; // the weighted number of datapoints per class
        double const n_total_; // the weighted number of datapoints
    };

} // namespace detail

/// \brief Functor that computes the gini score.
///
/// This functor is typically selected indirectly by passing the value <tt>RF_GINI</tt>
/// to vigra::rf3::RandomForestOptions::split().
class GiniScore
{
public:
    double operator()(std::vector<double> const & priors,
                      std::vector<double> const & counts, double n_total, double n_left) const
    {
        double const n_right = n_total - n_left;
        double gini_left = 1.0;
        double gini_right = 1.0;
        for (size_t i = 0; i < counts.size(); ++i)
        {
            double const p_left = counts[i] / n_left;
            double const p_right = (priors[i] - counts[i]) / n_right;
            gini_left -= (p_left*p_left);
            gini_right -= (p_right*p_right);
        }
        return n_left*gini_left + n_right*gini_right;
    }

        // needed for Gini-based variable importance calculation
    template <typename LABELS, typename WEIGHTS, typename ITER>
    static double region_score(LABELS const & labels, WEIGHTS const & weights, ITER begin, ITER end)
    {
        // Count the occurences.
        std::vector<double> counts;
        double total = 0.0;
        for (auto it = begin; it != end; ++it)
        {
            auto const d = *it;
            auto const lbl = labels[d];
            if (counts.size() <= lbl)
            {
                counts.resize(lbl+1, 0.0);
            }
            counts[lbl] += weights[d];
            total += weights[d];
        }

        // Compute the gini.
        double gini = total;
        for (auto x : counts)
        {
            gini -= x*x/total;
        }
        return gini;
    }
};

/// \brief Functor that computes the entropy score.
///
/// This functor is typically selected indirectly by passing the value <tt>RF_ENTROPY</tt>
/// to vigra::rf3::RandomForestOptions::split().
class EntropyScore
{
public:
    double operator()(std::vector<double> const & priors, std::vector<double> const & counts, double n_total, double n_left) const
    {
        double const n_right = n_total - n_left;
        double ig = 0;
        for (size_t i = 0; i < counts.size(); ++i)
        {
            double c = counts[i];
            if (c != 0)
                ig -= c * std::log(c / n_left);

            c = priors[i] - c;
            if (c != 0)
                ig -= c * std::log(c / n_right);
        }
        return ig;
    }

    template <typename LABELS, typename WEIGHTS, typename ITER>
    double region_score(LABELS const & /*labels*/, WEIGHTS const & /*weights*/, ITER /*begin*/, ITER /*end*/) const
    {
        vigra_fail("EntropyScore::region_score(): Not implemented yet.");
        return 0.0; // FIXME
    }
};

/// \brief Functor that computes the Kolmogorov-Smirnov score.
///
/// Actually, it reutrns the negated KSD score, because we want to minimize.
///
/// This functor is typically selected indirectly by passing the value <tt>RF_KSD</tt>
/// to vigra::rf3::RandomForestOptions::split().
class KolmogorovSmirnovScore
{
public:
    double operator()(std::vector<double> const & priors, std::vector<double> const & counts, double /*n_total*/, double /*n_left*/) const // Fix unused parameter warning, but leave in to not break compatibility with overall API
    {
        double const eps = 1e-10;
        double nnz = 0;
        std::vector<double> norm_counts(counts.size(), 0.0);
        for (size_t i = 0; i < counts.size(); ++i)
        {
            if (priors[i] > eps)
            {
                norm_counts[i] = counts[i] / priors[i];
                ++nnz;
            }
        }
        if (nnz < eps)
            return 0.0;

        // NOTE to future self:
        // In std::accumulate, it makes a huge difference whether you use 0 or 0.0 as init. Think about that before making changes.
        double const mean = std::accumulate(norm_counts.begin(), norm_counts.end(), 0.0) / nnz;

        // Compute the sum of the squared distances.
        double ksd = 0.0;
        for (size_t i = 0; i < norm_counts.size(); ++i)
        {
            if (priors[i] != 0)
            {
                double const v = (mean-norm_counts[i]);
                ksd += v*v;
            }
        }
        return -ksd;
    }

    template <typename LABELS, typename WEIGHTS, typename ITER>
    double region_score(LABELS const & /*labels*/, WEIGHTS const & /*weights*/, ITER /*begin*/, ITER /*end*/) const
    {
        vigra_fail("KolmogorovSmirnovScore::region_score(): Region score not available for the Kolmogorov-Smirnov split.");
        return 0.0;
    }
};

// This struct holds the depth and the weighted number of datapoints per class of a single node.
template <typename ARR>
struct RFNodeDescription
{
public:
    RFNodeDescription(size_t depth, ARR const & priors)
        :
        depth_(depth),
        priors_(priors)
    {}
    size_t depth_;
    ARR const & priors_;
};



// Return true if the given node is pure.
template <typename LABELS, typename ITER>
bool is_pure(LABELS const & /*labels*/, RFNodeDescription<ITER> const & desc)
{
    bool found = false;
    for (auto n : desc.priors_)
    {
        if (n > 0)
        {
            if (found)
                return false;
            else
                found = true;
        }
    }
    return true;
}

/// @brief Random forest 'node purity' stop criterion.
///
/// Stop splitting a node when it contains only instanes of a single class.
class PurityStop
{
public:
    template <typename LABELS, typename ITER>
    bool operator()(LABELS const & labels, RFNodeDescription<ITER> const & desc) const
    {
        return is_pure(labels, desc);
    }
};

/// @brief Random forest 'maximum depth' stop criterion.
///
/// Stop splitting a node when the its depth reaches a given value or when it is pure.
class DepthStop
{
public:
    /// @brief Constructor: terminate tree construction at \a max_depth.
    DepthStop(size_t max_depth)
        :
        max_depth_(max_depth)
    {}

    template <typename LABELS, typename ITER>
    bool operator()(LABELS const & labels, RFNodeDescription<ITER> const & desc) const
    {
        if (desc.depth_ >= max_depth_)
            return true;
        else
            return is_pure(labels, desc);
    }
    size_t max_depth_;
};

/// @brief Random forest 'number of datapoints' stop criterion.
///
/// Stop splitting a node when it contains too few instances or when it is pure.
class NumInstancesStop
{
public:
    /// @brief Constructor: terminate tree construction when node contains less than \a min_n instances.
    NumInstancesStop(size_t min_n)
        :
        min_n_(min_n)
    {}

    template <typename LABELS, typename ARR>
    bool operator()(LABELS const & labels, RFNodeDescription<ARR> const & desc) const
    {
        typedef typename ARR::value_type value_type;
        if (std::accumulate(desc.priors_.begin(), desc.priors_.end(), static_cast<value_type>(0)) <= min_n_)
            return true;
        else
            return is_pure(labels, desc);
    }
    size_t min_n_;
};

/// @brief Random forest 'node complexity' stop criterion.
///
/// Stop splitting a node when it allows for too few different data arrangements.
/// This includes purity, which offers only a sinlge data arrangement.
class NodeComplexityStop
{
public:
    /// @brief Constructor: stop when fewer than <tt>1/tau</tt> label arrangements are possible.
    NodeComplexityStop(double tau = 0.001)
        :
        logtau_(std::log(tau))
    {
        vigra_precondition(tau > 0 && tau < 1, "NodeComplexityStop(): Tau must be in the open interval (0, 1).");
    }

    template <typename LABELS, typename ARR>
    bool operator()(LABELS const & /*labels*/, RFNodeDescription<ARR> const & desc) // Fix unused parameter, but leave in for API compatability
    {
        typedef typename ARR::value_type value_type;

        // Count the labels.
        size_t const total = std::accumulate(desc.priors_.begin(), desc.priors_.end(), static_cast<value_type>(0));

        // Compute log(prod_k(n_k!)).
        size_t nnz = 0;
        double lg = 0.0;
        for (auto v : desc.priors_)
        {
            if (v > 0)
            {
                ++nnz;
                lg += loggamma(static_cast<double>(v+1));
            }
        }
        lg += loggamma(static_cast<double>(nnz+1));
        lg -= loggamma(static_cast<double>(total+1));
        if (nnz <= 1)
            return true;

        return lg >= logtau_;
    }

    double logtau_;
};

enum RandomForestOptionTags
{
    RF_SQRT,
    RF_LOG,
    RF_CONST,
    RF_ALL,
    RF_GINI,
    RF_ENTROPY,
    RF_KSD
};


/** \brief  Options class for \ref vigra::rf3::RandomForest version 3.

  <b>\#include</b> \<vigra/random_forest_3.hxx\><br/>
  Namespace: vigra::rf3
*/
class RandomForestOptions
{
public:

    RandomForestOptions()
        :
        tree_count_(255),
        features_per_node_(0),
        features_per_node_switch_(RF_SQRT),
        bootstrap_sampling_(true),
        resample_count_(0),
        split_(RF_GINI),
        max_depth_(0),
        node_complexity_tau_(-1),
        min_num_instances_(1),
        use_stratification_(false),
        n_threads_(-1),
        class_weights_()
    {}

    /**
     * @brief The number of trees.
     *
     * Default: 255
     */
    RandomForestOptions & tree_count(int p_tree_count)
    {
        tree_count_ = p_tree_count;
        return *this;
    }

    /**
     * @brief The number of features that are considered when computing the split.
     *
     * @param p_features_per_node the number of features
     *
     * Default: use sqrt of the total number of features.
     */
    RandomForestOptions & features_per_node(int p_features_per_node)
    {
        features_per_node_switch_ = RF_CONST;
        features_per_node_ = p_features_per_node;
        return *this;
    }

    /**
     * @brief The number of features that are considered when computing the split.
     *
     * @param p_features_per_node_switch possible values: <br/>
                    <tt>vigra::rf3::RF_SQRT</tt> (use square root of total number of features, recommended for classification), <br/>
                    <tt>vigra::rf3::RF_LOG</tt> (use logarithm of total number of features, recommended for regression), <br/>
                    <tt>vigra::rf3::RF_ALL</tt> (use all features).
     *
     * Default: <tt>vigra::rf3::RF_SQRT</tt>
     */
    RandomForestOptions & features_per_node(RandomForestOptionTags p_features_per_node_switch)
    {
        vigra_precondition(p_features_per_node_switch == RF_SQRT ||
                           p_features_per_node_switch == RF_LOG ||
                           p_features_per_node_switch == RF_ALL,
                           "RandomForestOptions::features_per_node(): Input must be RF_SQRT, RF_LOG or RF_ALL.");
        features_per_node_switch_ = p_features_per_node_switch;
        return *this;
    }

    /**
     * @brief Use bootstrap sampling.
     *
     * Default: true
     */
    RandomForestOptions & bootstrap_sampling(bool b)
    {
        bootstrap_sampling_ = b;
        return *this;
    }

    /**
     * @brief If resample_count is greater than zero, the split in each node is computed using only resample_count data points.
     *
     * Default: \a n = 0 (don't resample in every node)
     */
    RandomForestOptions & resample_count(size_t n)
    {
        resample_count_ = n;
        bootstrap_sampling_ = false;
        return *this;
    }

    /**
     * @brief The split criterion.
     *
     * @param p_split possible values: <br/>
                  <tt>vigra::rf3::RF_GINI</tt> (use Gini criterion, \ref vigra::rf3::GiniScorer), <br/>
                  <tt>vigra::rf3::RF_ENTROPY</tt> (use entropy criterion, \ref vigra::rf3::EntropyScorer),  <br/>
                  <tt>vigra::rf3::RF_KSD</tt> (use Kolmogorov-Smirnov criterion, \ref vigra::rf3::KSDScorer).
     *
     * Default: <tt>vigra::rf3::RF_GINI</tt>
     */
    RandomForestOptions & split(RandomForestOptionTags p_split)
    {
        vigra_precondition(p_split == RF_GINI ||
                           p_split == RF_ENTROPY ||
                           p_split == RF_KSD,
                           "RandomForestOptions::split(): Input must be RF_GINI, RF_ENTROPY or RF_KSD.");
        split_ = p_split;
        return *this;
    }

    /**
     * @brief Do not split a node if its depth is greater or equal to max_depth.
     *
     * Default: \a d = 0 (don't use depth as a termination criterion)
     */
    RandomForestOptions & max_depth(size_t d)
    {
        max_depth_ = d;
        return *this;
    }

    /**
     * @brief Value of the node complexity termination criterion.
     *
     * Default: \a tau = -1 (don't use complexity as a termination criterion)
     */
    RandomForestOptions & node_complexity_tau(double tau)
    {
        node_complexity_tau_ = tau;
        return *this;
    }

    /**
     * @brief Do not split a node if it contains less than min_num_instances data points.
     *
     * Default: \a n = 1 (don't use instance count as a termination criterion)
     */
    RandomForestOptions & min_num_instances(size_t n)
    {
        min_num_instances_ = n;
        return *this;
    }

    /**
     * @brief Use stratification when creating the bootstrap samples.
     *
     * That is, preserve the proportion between the number of class instances exactly
     * rather than on average.
     *
     * Default: false
     */
    RandomForestOptions & use_stratification(bool b)
    {
        use_stratification_ = b;
        return *this;
    }

    /**
     * @brief The number of threads that are used in training.
     *
     * \a n = -1 means use number of cores, \a n = 0 means single-threaded training.
     *
     * Default: \a n = -1 (use as many threads as there are cores in the machine).
     */
    RandomForestOptions & n_threads(int n)
    {
        n_threads_ = n;
        return *this;
    }

    /**
     * @brief Each datapoint is weighted by its class weight. By default, each class has weight 1.
     * @details
     * The classes in the random forest training have to follow a strict ordering. The weights must be given in that order.
     * Example:
     * You have the classes 3, 8 and 5 and use the vector {0.2, 0.3, 0.4} for the class weights.
     * The ordering of the classes is 3, 5, 8, so class 3 will get weight 0.2, class 5 will get weight 0.3
     * and class 8 will get weight 0.4.
     */
    RandomForestOptions & class_weights(std::vector<double> const & v)
    {
        class_weights_ = v;
        return *this;
    }

    /**
     * @brief Get the actual number of features per node.
     *
     * @param total the total number of features
     *
     * This function is normally only called internally before training is started.
     */
    size_t get_features_per_node(size_t total) const
    {
        if (features_per_node_switch_ == RF_SQRT)
            return std::ceil(std::sqrt(total));
        else if (features_per_node_switch_ == RF_LOG)
            return std::ceil(std::log(total));
        else if (features_per_node_switch_ == RF_CONST)
            return features_per_node_;
        else if (features_per_node_switch_ == RF_ALL)
            return total;
        vigra_fail("RandomForestOptions::get_features_per_node(): Unknown switch.");
        return 0;
    }

    int tree_count_;
    int features_per_node_;
    RandomForestOptionTags features_per_node_switch_;
    bool bootstrap_sampling_;
    size_t resample_count_;
    RandomForestOptionTags split_;
    size_t max_depth_;
    double node_complexity_tau_;
    size_t min_num_instances_;
    bool use_stratification_;
    int n_threads_;
    std::vector<double> class_weights_;

};



template <typename LabelType>
class ProblemSpec
{
public:

    ProblemSpec()
        :
        num_features_(0),
        num_instances_(0),
        num_classes_(0),
        distinct_classes_(),
        actual_mtry_(0),
        actual_msample_(0)
    {}

    ProblemSpec & num_features(size_t n)
    {
        num_features_ = n;
        return *this;
    }

    ProblemSpec & num_instances(size_t n)
    {
        num_instances_ = n;
        return *this;
    }

    ProblemSpec & num_classes(size_t n)
    {
        num_classes_ = n;
        return *this;
    }

    ProblemSpec & distinct_classes(std::vector<LabelType> v)
    {
        distinct_classes_ = v;
        num_classes_ = v.size();
        return *this;
    }

    ProblemSpec & actual_mtry(size_t m)
    {
        actual_mtry_ = m;
        return *this;
    }

    ProblemSpec & actual_msample(size_t m)
    {
        actual_msample_ = m;
        return *this;
    }

    bool operator==(ProblemSpec const & other) const
    {
        #define COMPARE(field) if (field != other.field) return false;
        COMPARE(num_features_);
        COMPARE(num_instances_);
        COMPARE(num_classes_);
        COMPARE(distinct_classes_);
        COMPARE(actual_mtry_);
        COMPARE(actual_msample_);
        #undef COMPARE
        return true;
    }

    size_t num_features_;
    size_t num_instances_;
    size_t num_classes_;
    std::vector<LabelType> distinct_classes_;
    size_t actual_mtry_;
    size_t actual_msample_;

};

//@}

} // namespace rf3

} // namespace vigra

#endif

