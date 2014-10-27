#ifndef RF_EARLY_STOPPING_P_HXX
#define RF_EARLY_STOPPING_P_HXX
#include <cmath>
#include <stdexcept>
#include "rf_common.hxx"

namespace vigra
{

#if 0    
namespace es_detail
{
    template<class T>
    T power(T const & in, int n)
    {
        T result = NumericTraits<T>::one();
        for(int ii = 0; ii < n ;++ii)
            result *= in;
        return result;
    }
}
#endif

/**Base class from which all EarlyStopping Functors derive.
 */
class StopBase
{
protected:
    ProblemSpec<> ext_param_;
    int tree_count_ ;
    bool is_weighted_;

public:
    template<class T>
    void set_external_parameters(ProblemSpec<T> const  &prob, int tree_count = 0, bool is_weighted = false)
    {
        ext_param_ = prob; 
        is_weighted_ = is_weighted;
        tree_count_ = tree_count;
    }

#ifdef DOXYGEN
        /** called after the prediction of a tree was added to the total prediction
         * \param weightIter Iterator to the weights delivered by current tree.
         * \param k          after kth tree
         * \param prob       Total probability array
         * \param totalCt    sum of probability array. 
         */
    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter weightIter, int k, MultiArrayView<2, T, C> const &  prob , double totalCt)
#else
    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter,  int /* k */, MultiArrayView<2, T, C> const & /* prob */, double /* totalCt */)
    {return false;}
#endif //DOXYGEN
};


/**Stop predicting after a set number of trees
 */
class StopAfterTree : public StopBase
{
public:
    double max_tree_p;
    int max_tree_;
    typedef StopBase SB;
    
    ArrayVector<double> depths;
    
    /** Constructor
     * \param max_tree number of trees to be used for prediction
     */
    StopAfterTree(double max_tree)
    :
        max_tree_p(max_tree)
    {}

    template<class T>
    void set_external_parameters(ProblemSpec<T> const  &prob, int tree_count = 0, bool is_weighted = false)
    {
        max_tree_ = ceil(max_tree_p * tree_count);
        SB::set_external_parameters(prob, tree_count, is_weighted);
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter,  int k, MultiArrayView<2, T, C> const & /* prob */, double /* totalCt */)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }
        if(k < max_tree_)
           return false;
        depths.push_back(double(k+1)/double(SB::tree_count_));
        return true;  
    }
};

/** Stop predicting after a certain amount of votes exceed certain proportion.
 *  case unweighted voting: stop if the leading class exceeds proportion * SB::tree_count_ 
 *  case weighted voting: stop if the leading class exceeds proportion * msample_ * SB::tree_count_ ;
 *                          (maximal number of votes possible in both cases)
 */
class StopAfterVoteCount : public StopBase
{
public:
    double proportion_;
    typedef StopBase SB;
    ArrayVector<double> depths;

    /** Constructor
     * \param proportion specify proportion to be used.
     */
    StopAfterVoteCount(double proportion)
    :
        proportion_(proportion)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter,  int k, MultiArrayView<2, T, C> const & prob, double /* totalCt */)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }


        if(SB::is_weighted_)
        {
            if(prob[argMax(prob)] > proportion_ *SB::ext_param_.actual_msample_* SB::tree_count_)
            {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return true;
            }
        }
        else
        {
            if(prob[argMax(prob)] > proportion_ * SB::tree_count_)
            {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return true;
            }
        }
        return false;
    }

};


/** Stop predicting if the 2norm of the probabilities does not change*/
class StopIfConverging : public StopBase

{
public:
    double thresh_;
    int num_;
    MultiArray<2, double> last_;
    MultiArray<2, double> cur_;
    ArrayVector<double> depths;
    typedef StopBase SB;

    /** Constructor
     * \param thresh: If the two norm of the probabilities changes less then thresh then stop
     * \param num   : look at atleast num trees before stopping
     */
    StopIfConverging(double thresh, int num = 10)
    :
        thresh_(thresh), 
        num_(num)
    {}

    template<class T>
    void set_external_parameters(ProblemSpec<T> const  &prob, int tree_count = 0, bool is_weighted = false)
    {
        last_.reshape(MultiArrayShape<2>::type(1, prob.class_count_), 0);
        cur_.reshape(MultiArrayShape<2>::type(1, prob.class_count_), 0);
        SB::set_external_parameters(prob, tree_count, is_weighted);
    }
    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> const & prob, double totalCt)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }
        if(k <= num_)
        {
            last_ = prob;
            last_/= last_.norm(1);
            return false;
        }
        else 
        {
            cur_ = prob;
            cur_ /= cur_.norm(1);
            last_ -= cur_;
            double nrm = last_.norm(); 
            if(nrm < thresh_)
            {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return true;
            }
            else
            {
                last_ = cur_;
            }
        }
        return false;
    }
};


/** Stop predicting if the margin prob(leading class) - prob(second class) exceeds a proportion
 *  case unweighted voting: stop if margin exceeds proportion * SB::tree_count_ 
 *  case weighted voting: stop if margin exceeds proportion * msample_ * SB::tree_count_ ;
 *                          (maximal number of votes possible in both cases)
 */
class StopIfMargin : public StopBase  
{
public:
    double proportion_;
    typedef StopBase SB;
    ArrayVector<double> depths;

    /** Constructor
     * \param proportion specify proportion to be used.
     */
    StopIfMargin(double proportion)
    :
        proportion_(proportion)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter,  int k, MultiArrayView<2, T, C> prob, double /* totalCt */)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }
        int index = argMax(prob);
        double a = prob[argMax(prob)];
        prob[argMax(prob)] = 0;
        double b = prob[argMax(prob)];
        prob[index] = a; 
        double margin = a - b;
        if(SB::is_weighted_)
        {
            if(margin > proportion_ *SB::ext_param_.actual_msample_ * SB::tree_count_)
            {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return true;
            }
        }
        else
        {
            if(prob[argMax(prob)] > proportion_ * SB::tree_count_)
            {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return true;
            }
        }
        return false;
    }
};


/**Probabilistic Stopping criterion (binomial test)
 *
 * Can only be used in a two class setting
 *
 * Stop if the Parameters estimated for the underlying binomial distribution
 * can be estimated with certainty over 1-alpha.
 * (Thesis, Rahul Nair Page 80 onwards: called the "binomial" criterion
 */
class StopIfBinTest : public StopBase  
{
public:
    double alpha_;  
    MultiArrayView<2, double> n_choose_k;
    /** Constructor
     * \param alpha specify alpha (=proportion) value for binomial test.
     * \param nck_ Matrix with precomputed values for n choose k
     * nck_(n, k) is n choose k. 
     */
    StopIfBinTest(double alpha, MultiArrayView<2, double> nck_)
    :
        alpha_(alpha),
        n_choose_k(nck_)
    {}
    typedef StopBase SB;
    
    /**ArrayVector that will contain the fraction of trees that was visited before terminating
     */
    ArrayVector<double> depths;

    double binomial(int N, int k, double p)
    {
//        return n_choose_k(N, k) * es_detail::power(p, k) *es_detail::power(1 - p, N-k);
        return n_choose_k(N, k) * std::pow(p, k) * std::pow(1 - p, N-k);
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k,
            MultiArrayView<2, T, C> const &prob, double totalCt)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }
        if(k < 10)
        {
            return false;
        }
        int index = argMax(prob);
        int n_a  = prob[index];
        int n_b  = prob[(index+1)%2];
        int n_tilde = (SB::tree_count_ - n_a + n_b);
        double p_a = double(n_b - n_a + n_tilde)/double(2* n_tilde);
        vigra_precondition(p_a <= 1, "probability should be smaller than 1");
        double cum_val = 0;
        int c = 0; 
  //      std::cerr << "prob: " << p_a << std::endl;
        if(n_a <= 0)n_a = 0;
        if(n_b <= 0)n_b = 0;
        for(int ii = 0; ii <= n_b + n_a;++ii)
        {
//            std::cerr << "nb +ba " << n_b + n_a << " " << ii <<std::endl;
            cum_val += binomial(n_b + n_a, ii, p_a); 
            if(cum_val >= 1 -alpha_)
            {
                c = ii;
                break;
            }
        }
//        std::cerr << c << " " << n_a << " " << n_b << " " << p_a <<   alpha_ << std::endl;
        if(c < n_a)
        {
            depths.push_back(double(k+1)/double(SB::tree_count_));
            return true;
        }

        return false;
    }
};

/**Probabilistic Stopping criteria. (toChange)
 *
 * Can only be used in a two class setting
 *
 * Stop if the probability that the decision will change after seeing all trees falls under
 * a specified value alpha.
 * (Thesis, Rahul Nair Page 80 onwards: called the "toChange" criterion
 */
class StopIfProb : public StopBase  
{
public:
    double alpha_;  
    MultiArrayView<2, double> n_choose_k;
    
    
    /** Constructor
     * \param alpha specify alpha (=proportion) value
     * \param nck_ Matrix with precomputed values for n choose k
     * nck_(n, k) is n choose k. 
     */
    StopIfProb(double alpha, MultiArrayView<2, double> nck_)
    :
        alpha_(alpha),
        n_choose_k(nck_)
    {}
    typedef StopBase SB;
    /**ArrayVector that will contain the fraction of trees that was visited before terminating
     */
    ArrayVector<double> depths;

    double binomial(int N, int k, double p)
    {
//        return n_choose_k(N, k) * es_detail::power(p, k) *es_detail::power(1 - p, N-k);
        return n_choose_k(N, k) * std::pow(p, k) * std::pow(1 - p, N-k);
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> prob, double totalCt)
    {
        if(k == SB::tree_count_ -1)
        {
                depths.push_back(double(k+1)/double(SB::tree_count_));
                return false;
        }
        if(k <= 10)
        {
            return false;
        }
        int index = argMax(prob);
        int n_a  = prob[index];
        int n_b  = prob[(index+1)%2];
        int n_needed = ceil(double(SB::tree_count_)/2.0)-n_a;
        int n_tilde = SB::tree_count_ - (n_a +n_b);
        if(n_tilde <= 0) n_tilde = 0;
        if(n_needed <= 0) n_needed = 0;
        double p = 0;
        for(int ii = n_needed; ii < n_tilde; ++ii)
            p += binomial(n_tilde, ii, 0.5);
        
        if(p >= 1-alpha_)
        {
            depths.push_back(double(k+1)/double(SB::tree_count_));
            return true;
        }

        return false;
    }
};


class DepthAndSizeStopping: public StopBase
{
public:
    int max_depth_;
    int min_size_;

    int max_depth_reached; //for debug maximum reached depth

    DepthAndSizeStopping() : max_depth_(-1), min_size_(0)
    {}

    /** Constructor DepthAndSize Criterion
     * Stop growing the tree if a certain depth or size is reached or make a
     * leaf if the node is smaller than a certain size. Note this is checked
     * before the split so it is still possible that smaller leafs are created
     */

    DepthAndSizeStopping(int depth, int size) :
        max_depth_(depth), min_size_(size)
    {}

    template<class T>
    void set_external_parameters(ProblemSpec<T> const &, int
     tree_count = 0, bool /* is_weighted_ */= false)
    {}

    template<class Region>
    bool operator()(Region& region)
    {
        if (region.depth() > max_depth_ + 1)
           throw std::runtime_error("violation in the stopping criterion");

        return (region.depth() > max_depth_) || (region.size() < min_size_) ;
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter, int /* k */,
            MultiArrayView<2, T, C> const &/* prob */, double /* totalCt */)
    {
        return true;
    }
};

} //namespace vigra;
#endif //RF_EARLY_STOPPING_P_HXX
