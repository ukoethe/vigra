/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
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
#ifndef VIGRA_RANDOM_FOREST_SPLIT_HXX
#define VIGRA_RANDOM_FOREST_SPLIT_HXX
#include <algorithm>
#include <map>
#include <numeric>
#include "../mathutil.hxx"
#include "../array_vector.hxx"
#include "../sized_int.hxx"
#include "../matrix.hxx"
#include "../random.hxx"
#include "../functorexpression.hxx"
#include "rf_nodeproxy.hxx"
#include "rf_sampling.hxx"
//#include "../hokashyap.hxx"
//#include "vigra/rf_helpers.hxx"

namespace vigra
{

// Incomplete Class to ensure that findBestSplit is always implemented in
// the derived classes of SplitBase
class CompileTimeError;


namespace detail
{
    template<class Tag>
    class Normalise
    {
    public:
        template<class Iter>
        static void exec(Iter begin, Iter  end)
        {}
    };

    template<>
    class Normalise<ClassificationTag>
    {
    public:
        template<class Iter>
        static void exec (Iter begin, Iter end)
        {
            double bla = std::accumulate(begin, end, 0.0);
            for(int ii = 0; ii < end - begin; ++ii)
                begin[ii] = begin[ii]/bla ;
        }
    };
}


/** Base Class for all SplitFunctors used with the \ref RandomForest class
    defines the interface used while learning a tree.
**/
template<class Tag>
class SplitBase
{
  public:

    typedef Tag           RF_Tag;
    typedef DT_StackEntry<ArrayVectorView<Int32>::iterator>
                                        StackEntry_t;

    ProblemSpec<>           ext_param_;

    NodeBase::T_Container_type          t_data;
    NodeBase::P_Container_type          p_data;

    NodeBase                            node_;

    /** returns the DecisionTree Node created by
        \ref findBestSplit or \ref makeTerminalNode.
    **/

    template<class T>
    void set_external_parameters(ProblemSpec<T> const & in)
    {
        ext_param_ = in;
        t_data.push_back(in.column_count_);
        t_data.push_back(in.class_count_);
    }

    NodeBase & createNode() 
    {
        return node_;
    }

    int classCount() const
    {
        return int(t_data[1]);
    }

    int featureCount() const
    {
        return int(t_data[0]);
    }

    /** resets internal data. Should always be called before
        calling findBestSplit or makeTerminalNode
    **/
    void reset()
    {
        t_data.resize(2);
        p_data.resize(0);
    }


    /** findBestSplit has to be implemented in derived split functor.
        these functions only insures That a CompileTime error is issued
        if no such method was defined.
    **/

    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2> labels,
                      Region region,
                      ArrayVector<Region> childs,
                      Random randint)
    {
        CompileTimeError SplitFunctor__findBestSplit_member_was_not_defined;
        return 0;
    }

    /** default action for creating a terminal Node.
        sets the Class probability of the remaining region according to
        the class histogram
    **/
    template<class T, class C, class T2,class C2, class Region, class Random>
    int makeTerminalNode(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region &                   region,
                      Random                     randint)
    {
        Node<e_ConstProbNode> ret(t_data, p_data);
        node_ = ret;
        if(ext_param_.class_weights_.size() != region.classCounts().size())
        {
        std::copy(          region.classCounts().begin(),
                            region.classCounts().end(),
                            ret.prob_begin());
        }
        else
        {
        std::transform(     region.classCounts().begin(),
                            region.classCounts().end(),
                            ext_param_.class_weights_.begin(),
                            ret.prob_begin(), std::multiplies<double>());
        }
        detail::Normalise<RF_Tag>::exec(ret.prob_begin(), ret.prob_end());
        ret.weights() = region.size();  
        return e_ConstProbNode;
    }


};

/** Functor to sort the indices of a feature Matrix by a certain dimension
**/
template<class DataMatrix>
class SortSamplesByDimensions
{
    DataMatrix const & data_;
    MultiArrayIndex sortColumn_;
    double thresVal_;
  public:

    SortSamplesByDimensions(DataMatrix const & data, 
                            MultiArrayIndex sortColumn,
                            double thresVal = 0.0)
    : data_(data),
      sortColumn_(sortColumn),
      thresVal_(thresVal)
    {}

    void setColumn(MultiArrayIndex sortColumn)
    {
        sortColumn_ = sortColumn;
    }
    void setThreshold(double value)
    {
        thresVal_ = value; 
    }

    bool operator()(MultiArrayIndex l, MultiArrayIndex r) const
    {
        return data_(l, sortColumn_) < data_(r, sortColumn_);
    }
    bool operator()(MultiArrayIndex l) const
    {
        return data_(l, sortColumn_) < thresVal_;
    }
};

template<class DataMatrix>
class DimensionNotEqual
{
    DataMatrix const & data_;
    MultiArrayIndex sortColumn_;

  public:

    DimensionNotEqual(DataMatrix const & data, 
                            MultiArrayIndex sortColumn)
    : data_(data),
      sortColumn_(sortColumn)
    {}

    void setColumn(MultiArrayIndex sortColumn)
    {
        sortColumn_ = sortColumn;
    }

    bool operator()(MultiArrayIndex l, MultiArrayIndex r) const
    {
        return data_(l, sortColumn_) != data_(r, sortColumn_);
    }
};

template<class DataMatrix>
class SortSamplesByHyperplane
{
    DataMatrix const & data_;
    Node<i_HyperplaneNode> const & node_;

  public:

    SortSamplesByHyperplane(DataMatrix              const & data, 
                            Node<i_HyperplaneNode>  const & node)
    :       
            data_(data), 
            node_()
    {}

    /** calculate the distance of a sample point to a hyperplane
     */
    double operator[](MultiArrayIndex l) const
    {
        double result_l = -1 * node_.intercept();
        for(int ii = 0; ii < node_.columns_size(); ++ii)
        {
            result_l +=     rowVector(data_, l)[node_.columns_begin()[ii]] 
                        *   node_.weights()[ii];
        }
        return result_l;
    }

    bool operator()(MultiArrayIndex l, MultiArrayIndex r) const
    {
        return (*this)[l]  < (*this)[r];
    }

};

/** makes a Class Histogram given indices in a labels_ array
 *  usage: 
 *      MultiArrayView<2, T2, C2> labels = makeSomeLabels()
 *      ArrayVector<int> hist(numberOfLabels(labels), 0);
 *      RandomForestClassCounter<T2, C2, ArrayVector> counter(labels, hist);
 *
 *      Container<int> indices = getSomeIndices()
 *      std::for_each(indices, counter);
 */
template <class DataSource, class CountArray>
class RandomForestClassCounter
{
    DataSource  const &     labels_;
    CountArray        &     counts_;

  public:

    RandomForestClassCounter(DataSource  const & labels, 
                             CountArray & counts)
    : labels_(labels),
      counts_(counts)
    {
        reset();
    }

    void reset()
    {
        counts_.init(0);
    }

    void operator()(MultiArrayIndex l) const
    {
        counts_[labels_[l]] +=1;
    }
};


/** Functor To Calculate the Best possible Split Based on the Gini Index
    given Labels and Features along a given Axis
*/

namespace detail
{
    template<int N>
    class ConstArr
    {
    public:
        double operator[](size_t) const
        {
            return (double)N;
        }
    };


}





/** Functor to calculate the gini impurity
 */
class GiniCriterion
{
public:
    /**caculate the weighted gini impurity based on class histogram
     * and class weights
     */
    template<class Array, class Array2>
    double operator()        (Array     const & hist, 
                              Array2    const & weights, 
                              double            total = 1.0) const
    {
        return impurity(hist, weights, total);
    }
    
    /** calculate the gini based impurity based on class histogram
     */
    template<class Array>
    double operator()(Array const & hist, double total = 1.0) const
    {
        return impurity(hist, total);
    }
    
    /** static version of operator(hist total)
     */
    template<class Array>
    static double impurity(Array const & hist, double total)
    {
        return impurity(hist, detail::ConstArr<1>(), total);
    }

    /** static version of operator(hist, weights, total)
     */
    template<class Array, class Array2>
    static double impurity   (Array     const & hist, 
                              Array2    const & weights, 
                              double            total)
    {

        int     class_count     = hist.size();
        double  gini            = 0.0;
        if(class_count == 2)
        {
            double w            = weights[0] * weights[1];
            gini                = w * (hist[0] * hist[1] / total);
        }
        else
        {
            for(int ii = 0; ii < class_count; ++ii)
            {
                double w        = weights[ii];
                gini           += w*( hist[ii]*( 1.0 - w * hist[ii]/total ) );
            }
        }
        return gini; 
    }
};


template <class DataSource, class Impurity= GiniCriterion>
class ImpurityLoss
{

    DataSource  const &         labels_;
    ArrayVector<double>        counts_;
    ArrayVector<double> const  class_weights_;
    double                      total_counts_;
    Impurity                    impurity_;

  public:

    template<class T>
    ImpurityLoss(DataSource  const & labels, 
                                ProblemSpec<T> const & ext_)
    : labels_(labels),
      counts_(ext_.class_count_, 0.0),
      class_weights_(ext_.class_weights_),
      total_counts_(0.0)
    {}

    void reset()
    {
        counts_.init(0);
        total_counts_ = 0.0;
    }

    template<class Counts>
    double increment_histogram(Counts const & counts)
    {
        std::transform(counts.begin(), counts.end(),
                       counts_.begin(), counts_.begin(),
                       std::plus<double>());
        total_counts_ = std::accumulate( counts_.begin(), 
                                         counts_.end(),
                                         0.0);
        return impurity_(counts_, class_weights_, total_counts_);
    }

    template<class Counts>
    double decrement_histogram(Counts const & counts)
    {
        std::transform(counts.begin(), counts.end(),
                       counts_.begin(), counts_.begin(),
                       std::minus<double>());
        total_counts_ = std::accumulate( counts_.begin(), 
                                         counts_.end(),
                                         0.0);
        return impurity_(counts_, class_weights_, total_counts_);
    }

    template<class Iter>
    double increment(Iter begin, Iter end)
    {
        for(Iter iter = begin; iter != end; ++iter)
        {
            counts_[labels_(*iter, 0)] +=1.0;
            total_counts_ +=1.0;
        }
        return impurity_(counts_, class_weights_, total_counts_);
    }

    template<class Iter>
    double decrement(Iter const &  begin, Iter const & end)
    {
        for(Iter iter = begin; iter != end; ++iter)
        {
            counts_[labels_(*iter,0)] -=1.0;
            total_counts_ -=1.0;
        }
        return impurity_(counts_, class_weights_, total_counts_);
    }

    template<class Iter, class Resp_t>
    double init (Iter begin, Iter end, Resp_t resp)
    {
        reset();
        std::copy(resp.begin(), resp.end(), counts_.begin());
        total_counts_ = std::accumulate(counts_.begin(), counts_.end(), 0.0); 
        return impurity_(counts_,class_weights_, total_counts_);
    }
    
    ArrayVector<double> const & response()
    {
        return counts_;
    }
};

template <class DataSource>
class RegressionForestCounter
{
    typedef MultiArrayShape<2>::type Shp;
    DataSource const &      labels_;
    ArrayVector <double>    mean_;
    ArrayVector <double>    variance_;
    ArrayVector <double>    tmp_;
    size_t                  count_;

    template<class T>
    RegressionForestCounter(DataSource const & labels, 
                            ProblemSpec<T> const & ext_)
    :
        labels_(labels),
        mean_(ext_.response_size, 0.0),
        variance_(ext_.response_size, 0.0),
        tmp_(ext_.response_size),
        count_(0)
    {}
    
    //  west's alorithm for incremental variance
    // calculation
    template<class Iter>
    double increment (Iter begin, Iter end)
    {
        for(Iter iter = begin; iter != end; ++iter)
        {
            ++count_;
            for(int ii = 0; ii < mean_.size(); ++ii)
                tmp_[ii] = labels_(*iter, ii) - mean_[ii]; 
            double f  = 1.0 / count_,
                   f1 = 1.0 - f;
            for(int ii = 0; ii < mean_.size(); ++ii)
                mean_[ii] += f*tmp_[ii]; 
            for(int ii = 0; ii < mean_.size(); ++ii)
                variance_[ii] += f1*sq(tmp_[ii]);
        }
        return std::accumulate(variance_.begin(), 
                               variance_.end(),
                               0.0,
                               std::plus<double>())
                /(count_ -1);
    }

    template<class Iter>
    double decrement (Iter begin, Iter end)
    {
        for(Iter iter = begin; iter != end; ++iter)
        {
            --count_;
            for(int ii = 0; ii < mean_.size(); ++ii)
                tmp_[ii] = labels_(*iter, ii) - mean_[ii]; 
            double f  = 1.0 / count_,
                   f1 = 1.0 + f;
            for(int ii = 0; ii < mean_.size(); ++ii)
                mean_[ii] -= f*tmp_[ii]; 
            for(int ii = 0; ii < mean_.size(); ++ii)
                variance_[ii] -= f1*sq(tmp_[ii]);
        }
        return std::accumulate(variance_.begin(), 
                               variance_.end(),
                               0.0,
                               std::plus<double>())
                /(count_ -1);
    }

    template<class Iter, class Resp_t>
    double init (Iter begin, Iter end, Resp_t resp)
    {
        reset();
        return increment(begin, end);
    }
    

    ArrayVector<double> const & response()
    {
        return mean_;
    }

    void reset()
    {
        mean_.init(0.0);
        variance_.init(0.0);
        count_ = 0; 
    }
};

template<class Tag, class Datatyp>
struct LossTraits;

struct LSQLoss
{};

template<class Datatype>
struct LossTraits<GiniCriterion, Datatype>
{
    typedef ImpurityLoss<Datatype, GiniCriterion> type;
};

template<class Datatype>
struct LossTraits<LSQLoss, Datatype>
{
    typedef RegressionForestCounter<Datatype> type;
};

template<class LineSearchLossTag>
class BestGiniOfColumn
{
public:
    ArrayVector<double>     class_weights_;
    ArrayVector<double>     bestCurrentCounts[2];
    double                  min_gini_;
    ptrdiff_t               min_index_;
    double                  min_threshold_;
    ProblemSpec<>           ext_param_;

    BestGiniOfColumn()
    {}

    template<class T> 
    BestGiniOfColumn(ProblemSpec<T> const & ext)
    :
        class_weights_(ext.class_weights_),
        ext_param_(ext)
    {
        bestCurrentCounts[0].resize(ext.class_count_);
        bestCurrentCounts[1].resize(ext.class_count_);
    }
    template<class T> 
    void set_external_parameters(ProblemSpec<T> const & ext)
    {
        class_weights_ = ext.class_weights_; 
        ext_param_ = ext;
        bestCurrentCounts[0].resize(ext.class_count_);
        bestCurrentCounts[1].resize(ext.class_count_);
    }
    /** calculate the best gini split along a Feature Column
     * \param column, the feature vector - has to support the [] operator
     * \param labels, the label vector 
     * \param begin 
     * \param end     (in and out)
     *                begin and end iterators to the indices of the
     *                samples in the current region. 
     *                the range begin - end is sorted by the column supplied
     *                during function execution.
     * \param class_counts
     *                class histogram of the range. 
     *
     *  precondition: begin, end valid range, 
     *                class_counts positive integer valued array with the 
     *                class counts in the current range.
     *                labels.size() >= max(begin, end); 
     *  postcondition:
     *                begin, end sorted by column given. 
     *                min_gini_ contains the minimum gini found or 
     *                NumericTraits<double>::max if no split was found.
     *                min_index_ countains the splitting index in the range
     *                or invalid data if no split was found.
     *                BestCirremtcounts[0] and [1] contain the 
     *                class histogram of the left and right region of 
     *                the left and right regions. 
     */
    template<   class DataSourceF_t,
                class DataSource_t, 
                class I_Iter, 
                class Array>
    void operator()(DataSourceF_t   const & column,
                    int                     g,
                    DataSource_t    const & labels,
                    I_Iter                & begin, 
                    I_Iter                & end,
                    Array           const & region_response)
    {
        std::sort(begin, end, 
                  SortSamplesByDimensions<DataSourceF_t>(column, g));
        typedef typename 
            LossTraits<LineSearchLossTag, DataSource_t>::type LineSearchLoss;
        LineSearchLoss left(labels, ext_param_);
        LineSearchLoss right(labels, ext_param_);

        

        min_gini_ = right.init(begin, end, region_response);
        min_threshold_ = *begin;
        min_index_     = 0;
        DimensionNotEqual<DataSourceF_t> comp(column, g); 
        
        I_Iter iter = begin;
        I_Iter next = std::adjacent_find(iter, end, comp);
        while( next  != end)
        {

            double loss = right.decrement(iter, next + 1) 
                    +     left.increment(iter , next + 1);
#ifdef CLASSIFIER_TEST
            if(loss < min_gini_ && !closeAtTolerance(loss, min_gini_))
#else
            if(loss < min_gini_ )
#endif 
            {
                bestCurrentCounts[0] = left.response();
                bestCurrentCounts[1] = right.response();
#ifdef CLASSIFIER_TEST
                min_gini_       = loss < min_gini_? loss : min_gini_;
#else
                min_gini_       = loss; 
#endif
                min_index_      = next - begin +1 ;
            }
            iter = next +1 ;
            next = std::adjacent_find(iter, end, comp);
        }
        next = begin + min_index_ - 1; 
        min_threshold_  = (column(*next,g) + column(*(next +1), g))/2;
    }

    template<class DataSource_t, class Iter, class Array>
    double loss_of_region(DataSource_t const & labels,
                          Iter & begin, 
                          Iter & end, 
                          Array const & region_response) const
    {
        typedef typename 
            LossTraits<LineSearchLossTag, DataSource_t>::type LineSearchLoss;
        LineSearchLoss region_loss(labels, ext_param_);
        return 
            region_loss.init(begin, end, region_response);
    }

};

template<class ColumnDecisionFunctor, class Tag>
class ThresholdSplit: public SplitBase<Tag>
{
  public:


    typedef SplitBase<Tag> SB;
    
    ArrayVector<Int32>          splitColumns;
    ColumnDecisionFunctor       bgfunc;

    double                      region_gini_;
    ArrayVector<double>         min_gini_;
    ArrayVector<ptrdiff_t>      min_indices_;
    ArrayVector<double>         min_thresholds_;

    int                         bestSplitIndex;

    double minGini() const
    {
        return min_gini_[bestSplitIndex];
    }
    int bestSplitColumn() const
    {
        return splitColumns[bestSplitIndex];
    }
    double bestSplitThreshold() const
    {
        return min_thresholds_[bestSplitIndex];
    }

    template<class T>
    void set_external_parameters(ProblemSpec<T> const & in)
    {
        SB::set_external_parameters(in);        
        bgfunc.set_external_parameters( SB::ext_param_);
        int featureCount_ = SB::ext_param_.column_count_;
        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
        min_gini_.resize(featureCount_);
        min_indices_.resize(featureCount_);
        min_thresholds_.resize(featureCount_);
    }


    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

        typedef typename Region::IndexIterator IndexIterator;
        if(region.size() == 0)
        {
           std::cerr << "SplitFunctor::findBestSplit(): stackentry with 0 examples encountered\n"
                        "continuing learning process...."; 
        }
        // calculate things that haven't been calculated yet. 
        
        if(std::accumulate(region.classCounts().begin(),
                           region.classCounts().end(), 0) != region.size())
        {
            RandomForestClassCounter<   MultiArrayView<2,T2, C2>, 
                                        ArrayVector<double> >
                counter(labels, region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }

        // Is the region pure already?
        region_gini_ = bgfunc.loss_of_region(labels,
                                             region.begin(), 
                                             region.end(),
                                             region.classCounts());
        if(region_gini_ <= SB::ext_param_.precision_)
            return  makeTerminalNode(features, labels, region, randint);

        // select columns  to be tried.
        for(int ii = 0; ii < SB::ext_param_.actual_mtry_; ++ii)
            std::swap(splitColumns[ii], 
                      splitColumns[ii+ randint(features.shape(1) - ii)]);

        // find the best gini index
        bestSplitIndex              = 0;
        double  current_min_gini    = region_gini_;
        int     num2try             = features.shape(1);
        for(int k=0; k<num2try; ++k)
        {
            //this functor does all the work
            bgfunc(features,
                   splitColumns[k],
                   labels, 
                   region.begin(), region.end(), 
                   region.classCounts());
            min_gini_[k]            = bgfunc.min_gini_; 
            min_indices_[k]         = bgfunc.min_index_;
            min_thresholds_[k]      = bgfunc.min_threshold_;
#ifdef CLASSIFIER_TEST
            if(     bgfunc.min_gini_ < current_min_gini
               &&  !closeAtTolerance(bgfunc.min_gini_, current_min_gini))
#else
            if(bgfunc.min_gini_ < current_min_gini)
#endif
            {
                current_min_gini = bgfunc.min_gini_;
                childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
                childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];
                childRegions[0].classCountsIsValid = true;
                childRegions[1].classCountsIsValid = true;

                bestSplitIndex   = k;
                num2try = SB::ext_param_.actual_mtry_;
            }
        }

        // did not find any suitable split
        if(closeAtTolerance(current_min_gini, region_gini_))
            return  makeTerminalNode(features, labels, region, randint);
        
        //create a Node for output
        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        node.threshold()    = min_thresholds_[bestSplitIndex];
        node.column()       = splitColumns[bestSplitIndex];
        
        // partition the range according to the best dimension 
        SortSamplesByDimensions<MultiArrayView<2, T, C> > 
            sorter(features, node.column(), node.threshold());
        IndexIterator bestSplit =
            std::partition(region.begin(), region.end(), sorter);
        // Save the ranges of the child stack entries.
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[0].rule = region.rule;
        childRegions[0].rule.push_back(std::make_pair(1, 1.0));
        childRegions[1].setRange(   bestSplit       , region.end()    );
        childRegions[1].rule = region.rule;
        childRegions[1].rule.push_back(std::make_pair(1, 1.0));

        return i_ThresholdNode;
    }
};

typedef  ThresholdSplit<BestGiniOfColumn<GiniCriterion> >                 GiniSplit;
typedef  ThresholdSplit<BestGiniOfColumn<LSQLoss>, RegressionTag>         RegressionSplit;

} //namespace vigra
#endif // VIGRA_RANDOM_FOREST_SPLIT_HXX
