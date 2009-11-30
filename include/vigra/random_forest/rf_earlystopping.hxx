#ifndef RF_EARLY_STOPPING_P_HXX
#define RF_EARLY_STOPPING_P_HXX

#include "rf_common.hxx"

namespace vigra
{

class StopBase
{
protected:
    ProblemSpec<> ext_param_;
    int tree_count_ ;
    bool is_weighted_;

public:
    template<class T>
    void set_external_parameters(ProblemSpec<T> &prob, int tree_count = 0, bool is_weighted = false)
    {
        ext_param_ = prob; 
        is_weighted_ = is_weighted;
        tree_count_ = tree_count;
    }
};


/**Stop predicting after a set number of trees
 */
class StopAfterTree : public StopBase
{
public:
    int max_tree_;
    typedef StopBase SB;
    
    StopAfterTree(int max_tree)
    :
        max_tree_(max_tree)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> prob, double totalCt)
    {
        if(k < max_tree_)
           return false;
        return true;  
    }
};

/** Stop predicting after a certain amount of votes exceed certain proportion.
 *  case unweighted voting: stop if the leading class exceeds proportion * SB::tree_count_ 
 *  case weighted votion: stop if the leading class exceeds proportion * msample_ * SB::tree_count_ ;
 *                          (maximal number of votes possible in both cases)
 */
class StopAfterVoteCount : public StopBase
{
public:
    double proportion_;
    typedef StopBase SB;

    StopAfterVoteCount(double proportion)
    :
        proportion_(proportion)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> prob, double totalCt)
    {
        if(SB::is_weighted_)
        {
            if(prob[argMax(prob)] > proportion_ *SB::ext_param_.actual_msample_* SB::tree_count_)
                return true;
        }
        else
        {
            if(prob[argMax(prob)] > proportion_ * SB::tree_count_)
                return true;
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
    ArrayVector<double> deviation;
    typedef StopBase SB;

    StopIfConverging(double thresh, int num = 10)
    :
        thresh_(thresh), 
        num_(num)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> prob, double totalCt)
    {

        deviation.push_back(0);
        double totalWeight = 0;
        for(int ii = 0; ii < prob.shape(1); ++ii)
            totalWeight += iter[ii]; 
        for(int ii = 0; ii < prob.shape(1); ++ii)
            deviation.back() += sq(iter[ii]/totalWeight - prob[ii]/totalCt);
        if(deviation.size() < num_)
            return false;
        else
        {
            double totaldev = std::accumulate(deviation.end() - num_, deviation.end(), 0.0);
            totaldev /= num_;
            if(totaldev < thresh_)
                return true;
        }
        return false;
    }
};


/** Stop predicting if the margin prob(leading class) - prob(second class) exceeds a proportion
 *  case unweighted voting: stop if margin exceeds proportion * SB::tree_count_ 
 *  case weighted votion: stop if margin exceeds proportion * msample_ * SB::tree_count_ ;
 *                          (maximal number of votes possible in both cases)
 */
class StopIfMargin : public StopBase  
{
public:
    double proportion_;
    typedef StopBase SB;

    StopIfMargin(double proportion)
    :
        proportion_(proportion)
    {}

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter iter,  int k, MultiArrayView<2, T, C> prob, double totalCt)
    {
        double a = prob[argMax(prob)];
        prob[argMax(prob)] = 0;
        double b = prob[argMax(prob)];
        prob[argMax(prob)] = a; 
        double margin = a - b;
        if(SB::is_weighted_)
        {
            if(margin > proportion_ *SB::ext_param_.actual_msample_ * SB::tree_count_)
                return true;
        }
        else
        {
            if(prob[argMax(prob)] > proportion_ * SB::tree_count_)
                return true;
        }
        return false;
    }
};

} //namespace vigra;
#endif //RF_EARLY_STOPPING_P_HXX
