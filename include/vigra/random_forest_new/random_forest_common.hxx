#ifndef VIGRA_RF_RANDOM_FOREST_COMMON_HXX
#define VIGRA_RF_RANDOM_FOREST_COMMON_HXX

#include <iterator>
#include <type_traits>
#include <cmath>

#include "../multi_array.hxx"
#include "../mathutil.hxx"

namespace vigra
{



template <typename T>
struct LessEqualSplitTest
{
public:
    LessEqualSplitTest(size_t dim = 0, T const & val = 0)
        :
        dim_(dim),
        val_(val)
    {}
    size_t operator()(MultiArrayView<1, T> const & features) const
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
            value_type const n = std::accumulate(vec.begin(), vec.end(), 0);
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



// namespace detail
// {
//     class GiniTracker
//     {
//     public:
//         GiniTracker(std::vector<size_t> const & counts)
//             :
//             counts_(counts),
//             n_(std::accumulate(counts.begin(), counts.end(), 0))
//         {
//             score_ = 1.0;
//             if (n_ != 0)
//             {
//                 for (auto c : counts)
//                 {
//                     double const p = c / static_cast<double>(n_);
//                     score_ -= p*p;
//                 }
//             }
//         }
//         double score() const
//         {
//             return n_ * score_;
//         }
//         double increment(size_t i)
//         {
//             score_ *= n_*n_;
//             score_ += 2 * (n_-counts_[i]);
//             score_ /= static_cast<double>((n_+1)*(n_+1));
//             ++counts_[i];
//             ++n_;
//             return n_ * score_;
//         }
//         double decrement(size_t i)
//         {
//             if (n_ == 1)
//             {
//                 score_ = 1.0;
//             }
//             else
//             {
//                 score_ *= n_*n_;
//                 score_ -= 2 * (n_-counts_[i]);
//                 score_ /= static_cast<double>((n_-1)*(n_-1));
//             }
//             --counts_[i];
//             --n_;
//             return n_ * score_;
//         }
//     private:
//         std::vector<size_t> counts_;
//         size_t n_;
//         double score_;
//     };
// }



// class MultiClassGiniScorer
// {
// public:

//     MultiClassGiniScorer(std::vector<size_t> const & priors)
//         :
//         split_found_(false),
//         best_split_(0),
//         best_dim_(0),
//         best_score_(std::numeric_limits<double>::max()),
//         priors_(priors)
//     {}

//     template <typename FEATURES, typename LABELS, typename WEIGHTS, typename ITER>
//     void operator()(
//         FEATURES const & features,
//         LABELS const & labels,
//         WEIGHTS const & weights,
//         ITER begin,
//         ITER end,
//         size_t dim
//     ){
//         if (begin == end)
//             return;

//         detail::GiniTracker left_(std::vector<size_t>(priors_.size(), 0));
//         detail::GiniTracker right_(priors_);

//         ITER next = begin;
//         ++next;
//         for (; next != end; ++begin, ++next)
//         {
//             size_t const label = static_cast<size_t>(labels(*begin));

//             double l = left_.increment(label);
//             double r = right_.decrement(label);
//             for (size_t i = 1; i < weights[*begin]; ++i)
//             {
//                 l = left_.increment(label);
//                 r = right_.decrement(label);
//             }

//             // Skip if there is no new split.
//             auto const left = features(*begin, dim);
//             auto const right = features(*next, dim);
//             if (left == right)
//                 continue;

//             // Update the score.
//             split_found_ = true;
//             double const s = l+r;
//             if (s < best_score_)
//             {
//                 best_score_ = s;
//                 best_split_ = 0.5*(left+right);
//                 best_dim_ = dim;
//             }
//         }
//     }

//     bool split_found_;
//     double best_split_;
//     size_t best_dim_;
//     double best_score_;

// private:

//     std::vector<size_t> priors_;
// };


class GiniScorer
{
public:

    GiniScorer(std::vector<size_t> const & priors)
        :
        split_found_(false),
        best_split_(0),
        best_dim_(0),
        best_score_(std::numeric_limits<double>::max()),
        priors_(priors),
        n_total_(std::accumulate(priors.begin(), priors.end(), 0))
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

        std::vector<size_t> counts(priors_.size(), 0);
        size_t n_left = 0;
        ITER next = begin;
        ++next;
        for (; next != end; ++begin, ++next)
        {
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
            double const s = score(counts, n_left);
            if (s < best_score_)
            {
                best_score_ = s;
                best_split_ = 0.5*(left+right);
                best_dim_ = dim;
            }
        }
    }

    bool split_found_;
    double best_split_;
    size_t best_dim_;
    double best_score_;

private:

    double score(std::vector<size_t> const & counts, double n_left) const
    {
        double const n_right = static_cast<double>(n_total_) - n_left;
        double gini_left = 1;
        double gini_right = 1;
        for (size_t i = 0; i < counts.size(); ++i)
        {
            double const p_left = counts[i] / n_left;
            double const p_right = (priors_[i] - counts[i]) / n_right;
            gini_left -= (p_left*p_left);
            gini_right -= (p_right*p_right);
        }
        return n_left*gini_left + n_right*gini_right;
    }

    std::vector<size_t> const priors_;
    size_t n_total_;
};



class EntropyScorer
{
public:

    EntropyScorer(std::vector<size_t> const & priors)
        :
        split_found_(false),
        best_split_(0),
        best_dim_(0),
        best_score_(std::numeric_limits<double>::max()),
        priors_(priors),
        n_total_(std::accumulate(priors.begin(), priors.end(), 0))
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

        std::vector<size_t> counts(priors_.size(), 0);
        size_t n_left = 0;
        ITER next = begin;
        ++next;
        for (; next != end; ++begin, ++next)
        {
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
            double const s = score(counts, n_left);
            if (s < best_score_)
            {
                best_score_ = s;
                best_split_ = 0.5*(left+right);
                best_dim_ = dim;
            }
        }
    }

    bool split_found_;
    double best_split_;
    size_t best_dim_;
    double best_score_;

private:

    double score(std::vector<size_t> const & counts, double n_left) const
    {
        double const n_right = static_cast<double>(n_total_) - n_left;

        double ig = 0;
        for (size_t i = 0; i < counts.size(); ++i)
        {
            auto c = counts[i];
            if (c != 0)
            {
                ig -= c * std::log(c / n_left);
            }
            c = priors_[i] - c;
            if (c != 0)
            {
                ig -= c * std::log(c / n_right);
            }
        }
        return ig;
    }

    std::vector<size_t> const priors_;
    size_t n_total_;
};



// class KSDScorer
// {
// public:

//     KSDScorer(std::vector<size_t> const & priors)
//         :
//         split_found_(false),
//         best_split_(0),
//         best_dim_(0),
//         best_score_(std::numeric_limits<double>::lowest()),
//         priors_(priors),
//         nnz_(priors.size() - std::count(priors.begin(), priors.end(), 0))
//     {}

//     template <typename FEATURES, typename LABELS, typename ITER>
//     void operator()(
//         FEATURES const & features,
//         LABELS const & labels,
//         ITER begin,
//         ITER end,
//         size_t dim
//     ){
//         if (begin == end)
//             return;

//         std::vector<double> increment(priors_.size(), 0);
//         for (size_t i = 0; i < priors_.size(); ++i)
//         {
//             if (priors_[i] != 0)
//             {
//                 increment[i] = 1.0 / priors_[i];
//             }
//         }

//         std::vector<double> counts(priors_.size(), 0);
//         ITER next = begin;
//         ++next;
//         for (; next != end; ++begin, ++next)
//         {
//             size_t const left_index = *begin;
//             size_t const right_index = *next;
//             size_t const label = static_cast<size_t>(labels(left_index));
//             counts[label] += increment[label];

//             // Skip if there is no new split.
//             auto const left = features(left_index, dim);
//             auto const right = features(right_index, dim);
//             if (left == right)
//                 continue;

//             // Update the score.
//             split_found_ = true;
//             double const s = score(counts);
//             if (s > best_score_)
//             {
//                 best_score_ = s;
//                 best_split_ = 0.5*(left+right);
//                 best_dim_ = dim;
//             }
//         }
//     }

//     bool split_found_;
//     double best_split_;
//     size_t best_dim_;
//     double best_score_;

// private:

//     double score(std::vector<double> const & counts) const
//     {
//         if (nnz_ == 0)
//             return 0.0;

//         // Compute the means.
//         double const mean = std::accumulate(counts.begin(), counts.end(), 0.0) / nnz_;

//         // Compute the sum of the squared distances.
//         double ksd = 0.0;
//         for (size_t i = 0; i < counts.size(); ++i)
//         {
//             if (priors_[i] != 0)
//             {
//                 ksd += (mean-counts[i])*(mean-counts[i]);
//             }
//         }
//         return ksd;
//     }

//     std::vector<size_t> const priors_;
//     size_t const nnz_;
// };



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



template <typename LABELS, typename ITER>
bool is_pure(LABELS const & labels, RFNodeDescription<ITER> const & desc)
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



class PurityStop
{
public:
    template <typename LABELS, typename ITER>
    bool operator()(LABELS const & labels, RFNodeDescription<ITER> const & desc) const
    {
        return is_pure(labels, desc);
    }
};



class DepthStop
{
public:
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



class NumInstancesStop
{
public:
    NumInstancesStop(size_t min_n)
        :
        min_n_(min_n)
    {}
    template <typename LABELS, typename ITER>
    bool operator()(LABELS const & labels, RFNodeDescription<ITER> const & desc) const
    {
        if (std::accumulate(desc.priors_.begin(), desc.priors_.end(), 0) <= min_n_)
            return true;
        else
            return is_pure(labels, desc);
    }
    size_t min_n_;
};



class NodeComplexityStop
{
public:
    NodeComplexityStop(double tau = 0.001)
        :
        logtau_(std::log(tau))
    {
        vigra_precondition(tau > 0 && tau < 1, "NodeComplexityStop(): Tau must be in the open interval (0, 1).");
    }

    template <typename LABELS, typename ITER>
    bool operator()(LABELS const & labels, RFNodeDescription<ITER> const & desc)
    {
        // Count the labels.
        size_t const total = std::accumulate(desc.priors_.begin(), desc.priors_.end(), 0);

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
    RF_ENTROPY
};



class RandomForestNewOptions
{
public:

    RandomForestNewOptions()
        :
        tree_count_(256),
        features_per_node_(0),
        features_per_node_switch_(RF_SQRT),
        bootstrap_sampling_(true),
        resample_count_(0),
        split_(RF_GINI),
        max_depth_(0),
        node_complexity_tau_(-1),
        min_num_instances_(1)
    {}

    RandomForestNewOptions & tree_count(int p_tree_count)
    {
        tree_count_ = p_tree_count;
        return *this;
    }

    RandomForestNewOptions & features_per_node(int p_features_per_node)
    {
        features_per_node_switch_ = RF_CONST;
        features_per_node_ = p_features_per_node;
        return *this;
    }

    RandomForestNewOptions & features_per_node(RandomForestOptionTags p_features_per_node_switch)
    {
        vigra_precondition(p_features_per_node_switch == RF_SQRT ||
                           p_features_per_node_switch == RF_LOG ||
                           p_features_per_node_switch == RF_ALL,
                           "RandomForestNewOptions::features_per_node(): Input must be RF_SQRT, RF_LOG or RF_ALL.");
        features_per_node_switch_ = p_features_per_node_switch;
        return *this;
    }

    RandomForestNewOptions & bootstrap_sampling(bool b)
    {
        bootstrap_sampling_ = b;
        return *this;
    }

    RandomForestNewOptions & resample_count(int n)
    {
        resample_count_ = n;
        bootstrap_sampling_ = false;
        return *this;
    }

    RandomForestNewOptions & split(RandomForestOptionTags p_split)
    {
        vigra_precondition(p_split == RF_GINI ||
                           p_split == RF_ENTROPY,
                           "RandomForestNewOptions::split(): Input must be RF_GINI or RF_ENTROPY.");
        split_ = p_split;
        return *this;
    }

    RandomForestNewOptions & max_depth(size_t d)
    {
        max_depth_ = d;
        return *this;
    }

    RandomForestNewOptions & node_complexity_tau(double tau)
    {
        node_complexity_tau_ = tau;
        return *this;
    }

    RandomForestNewOptions & min_num_instances(size_t n)
    {
        min_num_instances_ = n;
        return *this;
    }

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
        vigra_fail("RandomForestNewOptions::get_features_per_node(): Unknown switch.");
        return 0;
    }

    int tree_count_;
    int features_per_node_;
    RandomForestOptionTags features_per_node_switch_;
    bool bootstrap_sampling_;
    int resample_count_;
    RandomForestOptionTags split_;
    size_t max_depth_;
    double node_complexity_tau_;
    size_t min_num_instances_;

};



} // namespace vigra

#endif