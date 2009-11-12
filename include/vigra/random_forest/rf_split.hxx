/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


/** Base Class for all SplitFunctors used with the \ref RandomForestn class
    defines the interface used while learning a tree.
**/

class SplitBase
{
  public:

    typedef ClassificationTag   		RF_Tag;
    typedef DT_StackEntry<ArrayVectorView<Int32>::iterator>
                                		StackEntry_t;

    ProblemSpec<> 			ext_param_;

    NodeBase::T_Container_type 			t_data;
    NodeBase::P_Container_type 			p_data;

    NodeBase 							node_;

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
                      Region & region,
                      Random randint)
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

        double bla = std::accumulate(ret.prob_begin(), ret.prob_end(), 0.0);
        for(int ii = 0; ii < ret.prob_size(); ++ii)
            ret.prob_begin()[ii] = ret.prob_begin()[ii]/bla ;
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

    SortSamplesByHyperplane(DataMatrix 				const & data, 
							Node<i_HyperplaneNode>	const & node)
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
            result_l += 	rowVector(data_, l)[node_.columns_begin()[ii]] 
						* 	node_.weights()[ii];
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
 *  	MultiArrayView<2, T2, C2> labels = makeSomeLabels()
 *  	ArrayVector<int> hist(numberOfLabels(labels), 0);
 *  	RandomForestClassCounter<T2, C2, ArrayVector> counter(labels, hist);
 *
 *  	Container<int> indices = getSomeIndices()
 *  	std::for_each(indices, counter);
 */
template <class DataSource, class CountArray>
class RandomForestClassCounter
{
    DataSource  const & 	labels_;
    CountArray        &  	counts_;

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
	double operator()  		 (Array 	const & hist, 
							  Array2 	const & weights, 
					  		  double 			total = 1.0) const
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
	static double impurity   (Array 	const & hist, 
							  Array2 	const & weights, 
					  		  double 			total)
	{

		int 	class_count 	= hist.size();
		double 	gini			= 0;
		if(class_count == 2)
		{
			double w 			= weights[0] * weights[1];
			gini 				= w * (hist[0] * hist[1] / total);
		}
		else
		{
			for(int ii = 0; ii < class_count; ++ii)
			{
				double w 		= weights[ii];
				gini    	   += w*( hist[ii]*( 1.0 - w * hist[ii]/total ) );
			}
		}
		return gini; 
	}
};

template<class DataSource, class ImpurityFunctor = GiniCriterion>
class Impurity
{
protected:
	typedef	ArrayVector<double>					Weight_t;
	typedef ArrayVector<double>					Hist_t;

	ImpurityFunctor 			const & 		impurity_;
	DataSource 					const & 		data_;
	Weight_t					const &	 		class_weights_;
public:
	Hist_t										class_histogram_;
protected:
	double 										total_counts_;
	RandomForestClassCounter<DataSource,
							 Hist_t		>		counter_;
	Hist_t										tmp_;
	RandomForestClassCounter<DataSource,
							 Hist_t		>		tmp_counter_;
public:
	
	/** construct with data, class count, weights, and impurity functor
	 * 	 
	 * \param data 	a Linear array containing the class correspondence of 
	 * 				each instance (The label vector) 
	 * 				- should support operator[]
	 *
	 * \param class_count
	 * 				number of classes.
	 * \param class_weights 
	 * 				(optional) how to weight the class histogram from which
	 * 				impurity is calculated
	 * \param impurity
	 * 				(optional) supply a Impurity function that was not default
	 * 				constucted.
	 */
	Impurity(DataSource 	 const & data,
			 Int32 				   	 class_count 	= 0, 
			 Weight_t		 const & class_weights 	= Weight_t(),
			 ImpurityFunctor const & impurity 		= ImpurityFunctor())
	: 
		impurity_(impurity), 
		data_(data),
		class_weights_(class_weights),
		class_histogram_(class_count, 0.0),
		counter_(data, class_histogram_),
		tmp_(class_count, 0.0),
		tmp_counter_(data, tmp_)
	{}

	/** calculate the impurity given Class histogram.
	 */
	template<class Counts>
	double calculate_impurity(Counts const & counts)
   	{
		std::copy(counts.begin(), counts.end(), class_histogram_.begin());
		total_counts_ = 0;
		total_counts_ = std::accumulate( class_histogram_.begin(), 
		   								 class_histogram_.end(),
										 total_counts_);
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}

	/** calculate the impurity given indices to the data array
	 */
	template<class Iter>
	double calculate_impurity(Iter const & begin, Iter const & end)
   	{
		std::for_each(begin, end, counter_);
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}
	
	
	
	/** remove instances and calculate impurity
	 * Note: Impurity does not remember the instances used while 
	 * 		 calculating the impurity. The class counts of the removed
	 * 		 instances are calculated and then the 
	 * 		 decrement_impurity(class_counts) is called.
	 *
	 * \param begin, end: begin and end iterators to the indices of data
	 * 					  elements that should be removed.
	 */
	template<class Counts>
	double decrement_impurity(Counts const & counts_)
   	{
		std::transform(class_histogram_.begin(), class_histogram_.end(),
					   counts_.begin(), 		 class_histogram_.begin(),
					   std::minus<double>());
		total_counts_ = 0;
		total_counts_ = std::accumulate( class_histogram_.begin(), 
		   								 class_histogram_.end(),
										 total_counts_);
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}
	
	template<class Iter>
	double decrement_impurity(Iter const & begin, Iter const & end)
   	{
		tmp_counter_.reset();
		std::for_each(begin, end, tmp_counter_);
		return decrement_impurity(tmp_);
	}





	/** add instances and calculate impurity given difference histogram
	 */
	template<class Counts>
	double increment_impurity(Counts const & counts_)
   	{
		std::transform(class_histogram_.begin(), class_histogram_.end(),
					   counts_.begin(), 		 class_histogram_.begin(),
					   std::plus<double>());
		total_counts_ = 0;
		total_counts_ = std::accumulate( class_histogram_.begin(), 
		   								 class_histogram_.end(),
										 total_counts_);
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}

	/** add instances and calculate impurity given range of samples 
	 */
	template<class Iter>
	double increment_impurity(Iter const & begin, Iter const & end)
   	{
		tmp_counter_.reset();
		std::for_each(begin, end, tmp_counter_);
		return increment_impurity(tmp_);
	}
};

class BestGiniOfColumn
{
public:
	ArrayVector<double>		class_weights_;
	ArrayVector<double>		bestCurrentCounts[2];
	double 					min_gini_;
	ptrdiff_t 				min_index_;
	double					min_threshold_;


    BestGiniOfColumn(int classCount = 0)
    :
        class_weights_(ArrayVector<double>(classCount, 1.0))
    {
        bestCurrentCounts[0].resize(classCount);
        bestCurrentCounts[1].resize(classCount);
    }

    BestGiniOfColumn(int 						class_count, 
					 ArrayVector<double>const & class_weights)
    :
        class_weights_(class_weights)
    {
        bestCurrentCounts[0].resize(class_count);
        bestCurrentCounts[1].resize(class_count);
    }
  
	/** calculate the best gini split along a Feature Column
	 * \param column, the feature vector - has to support the [] operator
	 * \param labels, the label vector 
	 * \param begin 
	 * \param end	  (in and out)
	 * 				  begin and end iterators to the indices of the
	 * 				  samples in the current region. 
	 * 				  the range begin - end is sorted by the column supplied
	 * 				  during function execution.
	 * \param class_counts
	 * 				  class histogram of the range. 
	 *
	 *  precondition: begin, end valid range, 
	 *  			  class_counts positive integer valued array with the 
	 *  			  class counts in the current range.
	 *  			  labels.size() >= max(begin, end); 
	 *  postcondition:
	 *  			  begin, end sorted by column given. 
	 *  			  min_gini_ contains the minimum gini found or 
	 *  			  NumericTraits<double>::max if no split was found.
	 *  			  min_index_ countains the splitting index in the range
	 *  			  or invalid data if no split was found.
	 *  			  BestCirremtcounts[0] and [1] contain the 
	 *  			  class histogram of the left and right region of 
	 *  			  the left and right regions. 
	 */
	template<	class DataSourceF_t,
				class DataSource_t, 
				class I_Iter, 
				class Array>
	void operator()(DataSourceF_t	const & column,
					DataSource_t	const & labels,
				  	I_Iter 				  & begin, 
				  	I_Iter 				  & end,
				  	Array 			const & class_counts)
	{
		std::sort(begin, end, 
				  SortSamplesByDimensions<DataSourceF_t>(column, 0));

		Impurity<DataSource_t> 
			left(labels, class_counts.size(), class_weights_);
		Impurity<DataSource_t> 
			right(labels, class_counts.size(), class_weights_);

		right.calculate_impurity(class_counts);

		min_gini_ = NumericTraits<double>::max();

		DimensionNotEqual<DataSourceF_t> comp(column, 0); 
		
		I_Iter iter = begin;
		I_Iter next = std::adjacent_find(iter, end, comp);
		while( next  != end)
		{

			double gini = right.decrement_impurity(iter, next + 1) 
					+ 	  left.increment_impurity(iter , next + 1);
			if(gini < min_gini_ && !closeAtTolerance(gini, min_gini_))
			{
				bestCurrentCounts[0] = left.class_histogram_;
				bestCurrentCounts[1] = right.class_histogram_;
				min_gini_ 		= gini < min_gini_? gini : min_gini_;
				min_index_ 		= next - begin +1 ;
				min_threshold_	= (column[*next] + column[*(next +1)])/2;
			}
			iter = next +1 ;
			next = std::adjacent_find(iter, end, comp);
		}

	}
};


template<class ColumnDecisionFunctor>
class ThresholdSplit: public SplitBase
{
  public:


    typedef SplitBase SB;
    
	ArrayVector<Int32> 			splitColumns;
	ColumnDecisionFunctor		bgfunc;

	double 						region_gini_;
    ArrayVector<double> 		min_gini_;
	ArrayVector<ptrdiff_t>		min_indices_;
	ArrayVector<double>			min_thresholds_;

    int 						bestSplitIndex;

	double minGini()
	{
		return min_gini_[bestSplitIndex];
	}
	int bestSplitColumn()
	{
		return splitColumns[bestSplitIndex];
	}

	template<class T>
    void set_external_parameters(ProblemSpec<T> const & in)
    {
        SB::set_external_parameters(in);        
		bgfunc = ColumnDecisionFunctor( ext_param_.class_count_, 
				 				  		ext_param_.class_weights_);
        int featureCount_ = ext_param_.column_count_;
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
		
		// calculate things that haven't been calculated yet. 
        if(!region.classCountsIsValid)
        {
            RandomForestClassCounter<	MultiArrayView<2,T2, C2>, 
										ArrayVector<Int32> >
          		counter(labels, region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }

		// Is the region pure already?
        region_gini_ = GiniCriterion::impurity(region.classCounts(),
											   region.size());
		if(region_gini_ == 0)
            return  makeTerminalNode(features, labels, region, randint);

		// select columns  to be tried.
		for(int ii = 0; ii < SB::ext_param_.actual_mtry_; ++ii)
            std::swap(splitColumns[ii], 
					  splitColumns[ii+ randint(features.shape(1) - ii)]);

		// find the best gini index
        bestSplitIndex 				= 0;
		double	current_min_gini 	= NumericTraits<double>::max();
        int 	num2try 			= features.shape(1);
        for(int k=0; k<num2try; ++k)
        {
			//this functor does all the work
            bgfunc(columnVector(features, splitColumns[k]),
				   labels, 
				   region.begin(), region.end(), 
				   region.classCounts());
			
			min_gini_[k]			= bgfunc.min_gini_; 
			min_indices_[k] 		= bgfunc.min_index_;
			min_thresholds_[k]		= bgfunc.min_threshold_;

            if(		bgfunc.min_gini_ < current_min_gini
			   &&  !closeAtTolerance(bgfunc.min_gini_, current_min_gini))
            {
				current_min_gini = bgfunc.min_gini_;
				childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
				childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];
                childRegions[0].classCountsIsValid = true;
                childRegions[1].classCountsIsValid = true;

				bestSplitIndex   = k;
                num2try = std::max(int(k), SB::ext_param_.actual_mtry_);
            }
        }

		// did not find any suitable split
        if(closeAtTolerance(current_min_gini, NumericTraits<double>::max()))
            return  makeTerminalNode(features, labels, region, randint);
		
		//create a Node for output
        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        node.threshold() 	= min_thresholds_[bestSplitIndex];
		node.column() 		= splitColumns[bestSplitIndex];
		
		// partition the range according to the best dimension 
        SortSamplesByDimensions<MultiArrayView<2, T, C> > 
			sorter(features, node.column(), node.threshold());
        std::partition(region.begin(), region.end(), sorter);

        // Save the ranges of the child stack entries.
        IndexIterator bestSplit = 	region.begin()  + 
									min_indices_[bestSplitIndex];
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[1].setRange(   bestSplit       , region.end()    );

        return i_ThresholdNode;
    }
};

typedef  ThresholdSplit<BestGiniOfColumn> GiniSplit;

} //namespace vigra
#endif // VIGRA_RANDOM_FOREST_SPLIT_HXX
