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

    RF_Traits::ProblemSpec_t 			ext_param_;

    NodeBase::T_Container_type 			t_data;
    NodeBase::P_Container_type 			p_data;

    NodeBase 							node_;

    /** returns the DecisionTree Node created by
        \ref findBestSplit or \ref makeTerminalNode.
    **/

    void set_external_parameters(RF_Traits::ProblemSpec_t const & in)
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

  public:

    SortSamplesByDimensions(DataMatrix const & data, 
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
        return data_(l, sortColumn_) < data_(r, sortColumn_);
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
class BestGiniFunctor
{
    public:

    ArrayVector<double>         classWeights_;
    int                      	classCount;
    ArrayVector<Int32>          currentCounts[2];
    ArrayVector<Int32>          bestCurrentCounts[2];
    double                      totalCounts[2];
    double                      minGini;
	double						totalGini; 
    BestGiniFunctor()
    {}
    /** Initialize with A Class Weight Array
    **/
    BestGiniFunctor(ArrayVector<double>const & classWeights)
    :
        classWeights_(classWeights),
        classCount(classWeights.size())
    {
        currentCounts[0].resize(classCount);
        currentCounts[1].resize(classCount);
        bestCurrentCounts[0].resize(classCount);
        bestCurrentCounts[1].resize(classCount);
    }

    /** Initialize with the number of classes in unweighted
        case
    **/
    BestGiniFunctor(int classCount_)
    :
        classWeights_(ArrayVector<double>(classCount_, 1.0)),
        classCount(classCount_)
    {
        currentCounts[0].resize(classCount);
        currentCounts[1].resize(classCount);
        bestCurrentCounts[0].resize(classCount);
        bestCurrentCounts[1].resize(classCount);
    }

    BestGiniFunctor(int classCount_, ArrayVector<double>const & classWeights)
    :
        classWeights_(ArrayVector<double>(classCount_, 1.0)),
        classCount(classCount_)
    {
        if(classWeights.size() == classCount)
            std::copy(classWeights.begin(), 
					  classWeights.end(), 
					  classWeights_.begin());
        currentCounts[0].resize(classCount);
        currentCounts[1].resize(classCount);
        bestCurrentCounts[0].resize(classCount);
        bestCurrentCounts[1].resize(classCount);
    }
    /** Calculates the best possible split of Region given the labels 
	 * in labels and the features in featureColumn
     */
    template<class T, class C, class T2,class C2, class Region>
    typename Region::IndexIterator 
			operator()(  MultiArrayView<2, T, C>     const & featureColumn,
                         MultiArrayView<2, T2, C2>   const & labels,
                         Region &                            region      )
    {
        typedef typename Region::IndexIterator Iterator;
        Iterator 
			bestSplit  = Iterator();

        currentCounts[0].init(0);
        std::copy(  region.classCounts().begin(),
                    region.classCounts().end(),
                    currentCounts[1].begin());

        totalCounts[0] = totalCounts[1] = 0;
        totalCounts[1] = std::inner_product(currentCounts[1].begin(),
                                            currentCounts[1].end(),
                                            classWeights_.begin(),
                                            0.0);
        minGini        = NumericTraits<double>::max();

		totalGini = 0.0;
		if(classCount == 2)
		{
			double w = classWeights_[0] *classWeights_[1];
			totalGini 	 =w* double(  currentCounts[1][0]
									* currentCounts[1][1]) 
						   / totalCounts[1];
		}
		else
		{
			for(int l=0; l<classCount; ++l)
			{
				double w    = classWeights_[l];
				totalGini += w * double(currentCounts[1][l])
							   * (1.0 - w 
									* double(currentCounts[1][l]) 
									/ totalCounts[1]);
			}
		}

        for(int m = 0; m < region.size()-1; ++m)
        {

            int label = labels[region[m]];
            currentCounts[0][label] += 1;
            totalCounts[0]          += 1;
            currentCounts[1][label] -= 1;
            totalCounts[1]          -= 1;

            if (featureColumn[region[m]] == featureColumn[region[m+1]])
                continue ;

            double gini = 0.0;
            if(classCount == 2)
            {
                double w = classWeights_[0] *classWeights_[1];
                gini 	 = w *  (double(  currentCounts[0][0]
								  	    * currentCounts[0][1]) 
								 / totalCounts[0] 
							   + double(  currentCounts[1][0]
										* currentCounts[1][1]) 
							     / totalCounts[1]);
            }
            else
            {
                for(int l=0; l<classCount; ++l)
                {
                    double w    = classWeights_[l];
                    gini += w * (  double(currentCounts[0][l]) 
								 * (1.0 - w
								        * double(currentCounts[0][l])
								   		/ totalCounts[0]) 

							    +  double(currentCounts[1][l])
								 * (1.0 - w 
										* double(currentCounts[1][l]) 
										/ totalCounts[1])
								 );
                }
            }

            if(gini < minGini && !closeAtTolerance(gini, minGini))
            {
                minGini         = gini;
                bestSplit       = region.begin()+m;
                for(int ii = 0; ii < 2; ++ii)
                    std::copy(  currentCounts[ii].begin(),
                                currentCounts[ii].end(),
                                bestCurrentCounts[ii].begin());
            }

        }
        return bestSplit;
    }
};

namespace detail
{
	template<int N>
	class ConstArr
	{
		static const double val = N;
	public:
		double operator[](size_t in)
		{
			return val;
		}
	};
}



class GiniCriterion
{
	template<class Array>
	double operator()(Array const & hist, Int32 total) const
	{
		return (*this)(hist, detail::ConstArr<1>(), total);
	}

	template<class Array, class Array2>
	double operator()(Array 	const & hist, 
					  Array2 	const & weights, 
					  double 			total = 1.0) const
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
	Hist_t										class_histogram_;
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
			 Int32 					 class_count, 
			 Weight_t		 const & class_weights 	= Weight_t(),
			 ImpurityFunctor const & impurity 		= ImpurityFunctor())
	: 
		impurity_(impurity), 
		data_(data),
		class_weights_(class_weights),
		class_histogram_(class_count, 0),
		counter_(data, class_histogram_),
		tmp_(class_count, 0),
		tmp_counter_(data, tmp_)
	{}

	/** calculate the impurity given Class histogram.
	 */
	template<class Counts>
	double calculate_impurity(Counts const & counts)
   	{
		std::copy(counts.begin(), counts.end(), class_histogram_.begin());
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
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}
	
	template<class Iter>
	double decrement_impurity(Iter const & begin, Iter const & end)
   	{
		tmp_counter_.reset();
		std::for_each(begin, end, tmp_counter_);
		return decrement_impurity(tmp_);
	}





	/** add instances and calculate impurity
	 */
	template<class Counts>
	double increment_impurity(Counts const & counts_)
   	{
		std::transform(class_histogram_.begin(), class_histogram_.end(),
					   counts_.begin(), 		 class_histogram_.begin(),
					   std::plus<double>());
		return impurity_(class_histogram_, class_weights_, total_counts_);
	}

	template<class Iter>
	double increment_impurity(Iter const & begin, Iter const & end)
   	{
		tmp_counter_.reset();
		std::for_each(begin, end, tmp_counter_);
		return increment_impurity(tmp_);
	}
};


/** calculate the gini of a region - Only use if you want to calculate 
 * the gini only once.
 */
template<class Region>
double calculate_gini(Region region, 
					 ArrayVector<double> classWeights_ 
					 				= ArrayVector<double>())
{
    int classCount
            = region.classCounts.size();
    if(classWeights_.size() == 0 )
        classWeights_.resize(classCount, 1.0);
    double gini
            = 0;
    double totalCounts
            = std::inner_product(region.classCounts().begin(),
                                 region.classCounts.end(),
                                 classWeights_.begin(),
                                 0.0);

    if(classCount == 2)
    {
        double w
            = classWeights_[0] *classWeights_[1];
        gini
            = w *(double(region.classCounts()[0]*region.classCounts()[1]) / totalCounts);
    }
    else
    {
        for(int l=0; l<classCount; ++l)
        {
            double w
                = region.classWeights_[l];
            gini
               += w*(   double(region.classCounts()[l]) *
                    (1.0 - w * double(region.classCounts()[l]) / totalCounts));
        }
    }
    return gini;
}



class GiniSplit: public SplitBase
{
  public:


    typedef SplitBase SB;

	double						giniDecrease;
    ArrayVector<Int32> splitColumns;
    BestGiniFunctor bgfunc;
    double minGini;
    int bestSplitColumn;

    void set_external_parameters(RF_Traits::ProblemSpec_t const & in)
    {
        SB::set_external_parameters(in);
        bgfunc = BestGiniFunctor( ext_param_.class_count_, ext_param_.class_weights_);
        int featureCount_ = ext_param_.column_count_;

        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
    }


    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

        typedef typename Region::IndexIterator IndexIterator;
            // if the region is pure than just make a terminal Node.
        if(!region.classCountsIsValid)
        {
            RandomForestClassCounter<MultiArrayView<2,T2, C2>, ArrayVector<Int32> >
                                    counter(labels,region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }
        if(region.isPure())
            return  makeTerminalNode(features, labels, region, randint);

            // select columns to be tried for split first mtry entries in
            // splitColumns
        Sampler<Random> sampler(    SB::ext_param_.actual_mtry_,
                                    columnCount(features),
                                    SamplingOptions()
                                        .sampleWithoutReplacement()
                                        .useExternalMemory(splitColumns),
                                    randint);


        sampler.sample();
            // This Thing sor:ts a vector by a given dimension in the multi array
            //features

        SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorter(features, 0);

            // find the best gini index
        bestSplitColumn = 0;
        minGini = NumericTraits<double>::max();
        IndexIterator           bestSplit = IndexIterator();
        //ArrayVector<double> &   classWeights = SB::opt_.classWeights;
        int num2try = features.shape(1);

        for(int k=0; k<num2try; ++k)
        {

            sorter.setColumn(splitColumns[k]);
            std::sort(region.begin(), region.end(), sorter);
            // The calculation of the gini was put into the functor because, other
            // split criteria need a sorting of the samples according to one dimension
            // as well.
            IndexIterator  curSplit = bgfunc.operator()(columnVector(features, splitColumns[k]), labels, region);

            if(bgfunc.minGini < minGini &&  !closeAtTolerance(bgfunc.minGini, minGini))
            {
                std::copy(bgfunc.bestCurrentCounts[0].begin(), bgfunc.bestCurrentCounts[0].end(),
                                                                childRegions[0].classCounts().begin());
                std::copy(bgfunc.bestCurrentCounts[1].begin(), bgfunc.bestCurrentCounts[1].end(),
                                                                childRegions[1].classCounts().begin());
//                childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
//                childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];
				giniDecrease 	= bgfunc.totalGini - bgfunc.minGini; 
				minGini         = bgfunc.minGini;
                bestSplit       = curSplit;
                bestSplitColumn   = splitColumns[k];

                childRegions[0].classCountsIsValid = true;
                childRegions[1].classCountsIsValid = true;
                num2try = std::max(int(k), SB::ext_param_.actual_mtry_);
            }
        }

        if(minGini == NumericTraits<double>::max())
            return  makeTerminalNode(features, labels, region, randint);

        // sort samples according to best split dimension.
        sorter.setColumn(bestSplitColumn);
        std::sort(region.begin(), region.end(), sorter);


        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        //fix threshold value and update childRegions
        node.threshold() = (    features(bestSplit[0], bestSplitColumn)
                            +   features(bestSplit[1], bestSplitColumn)) / 2.0;
        node.column() = bestSplitColumn;

        ++bestSplit;

        // Save the childRanges.
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[1].setRange(   bestSplit       , region.end()    );
        return i_ThresholdNode;
    }
};

/**
template<class Options>
class OObGiniSplit: public GiniSplit<Options>
{
  public:


    typedef GiniSplit<Options> SB; //BaseType

    OObGiniSplit(int featureCount_, int classCount_, Options opt)
        :   SB(featureCount_, classCount_, opt)
    {
    }

    template<class T, class C, class T2,class C2, class Region, class Random>
    int makeTerminalNode(   MultiArrayView<2, T, C> features,
                            MultiArrayView<2, T2, C2>  labels,
                            Region & region,
                            Random randint)
    {
        Node<e_ConstProbNode> ret(SB::t_data, SB::p_data);
        SB::node_ = ret;
        ArrayVector<Int32> oobClassCounts(region.classCounts().size(), 0);
            RandomForestClassCounter<T2, C2, ArrayVector<Int32> >
                                    counter(labels,oobClassCounts);
            std::for_each(  region.oob_begin(), region.oob_end(), counter);
        std::transform(oobClassCounts.begin(), oobClassCounts.end(), SB::opt_.classWeights.begin(),
                       ret.prob_begin(), std::multiplies<double>());
        double bla = std::accumulate(ret.prob_begin(), ret.prob_end(), 0.0);
        for(int ii = 0; ii < ret.prob_size(); ++ii)
            ret.prob_begin()[ii] = ret.prob_begin()[ii]/bla ;
        return e_ConstProbNode;
    }

    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {
        typedef typename Region::IndexIterator IndexIterator;
        int normalBestSplit = SB::findBestSplit(features,
                                                labels,
                                                region,
                                                childRegions,
                                                randint);
        if(normalBestSplit == e_ConstProbNode)
        {
            SB::reset();
            //use oob data instead of ib data
            return makeTerminalNode(features, labels, region, randint);
        }
        else
        {

            SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorter(features, SB::bestSplitColumn);
            std::sort(region.oob_begin(), region.oob_end(), sorter);
            double threshold = Node<i_ThresholdNode>(SB::t_data, SB::p_data, 2).threshold();

#if 0
            for(int ii = 0; ii < region.oob_size(); ++ii)
            {
                if(ii%10 == 0)std::cerr << std::endl << int(ii/10) << "    ";
                std::cerr << features(region.oob_begin()[ii], SB::bestSplitColumn) << " ";
            }
#endif
            IndexIterator oobBestSplit = region.oob_begin();

            for(int ii = 0; ii < region.oob_size(); ++ii)
                if(features(region.oob_begin()[ii], SB::bestSplitColumn) > threshold)
                {
                    oobBestSplit = region.oob_begin() + ii;
                    break;
                }

            //std::cerr <<std::endl<< threshold<<" " << oobBestSplit - region.oob_begin()<<  " " << region.oob_end() - oobBestSplit <<  std::endl << std::endl;
            childRegions[0].set_oob_range(   region.oob_begin()  , oobBestSplit       );
            childRegions[1].set_oob_range(   oobBestSplit       , region.oob_end()    );
            // if one of the oob Ranges is empty return a Terminal node;
            if(childRegions[0].oob_size() == 0 || childRegions[1].oob_size() == 0)
            {
                SB::reset();
                return makeTerminalNode(features, labels, region, randint);
            }
        }
        return i_ThresholdNode;
    }

};
/*
template<class Options>
class UninformedSplit: public SplitBase
{
  public:


    typedef SplitBase SB;


    ArrayVector<Int32> splitColumns;
    BestGiniFunctor bgfunc;
    double minGini;
    int bestSplitColumn;
    UninformedSplit(int featureCount_, int classCount_, Options opt)
        :   SplitBase<Options>(featureCount_, classCount_, opt),
            bgfunc(classCount_, SB::opt_.classWeights)
    {
        if(SB::opt_.mtry == 0)
            SB::opt_.mtry = int(std::floor(std::sqrt(double(featureCount_)) + 0.5));
        SB::opt_.classWeights.resize(classCount_, 1.0);
        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
    }


    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

        typedef typename Region::IndexIterator IndexIterator;
            // if the region is pure than just make a terminal Node.
        if(!region.classCountsIsValid)
        {
            RandomForestClassCounter<T2, C2, ArrayVector<Int32> >
                                    counter(labels,region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }
        if(region.isPure())
            return  makeTerminalNode(features, labels, region, randint);

            // select columns to be tried for split first mtry entries in
            // splitColumns
        Sampler<Random> sampler(    SB::opt_.mtry,
                                    columnCount(features),
                                    SamplingOptions()
                                        .sampleWithoutReplacement()
                                        .useExternalMemory(splitColumns),
                                    randint);


        sampler.sample();
            // This Thing sorts a vector by a given dimension in the multi array
            //features

        SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorter(features, 0);

            // find the best gini index
        bestSplitColumn = 0;
        minGini = NumericTraits<double>::max();
        IndexIterator           bestSplit = IndexIterator();
        //ArrayVector<double> &   classWeights = SB::opt_.classWeights;
        int num2try = features.shape(1);

        for(int k=0; k<num2try; ++k)
        {

            sorter.setColumn(splitColumns[k]);
            std::sort(region.begin(), region.end(), sorter);
            // The calculation of the gini was put into the functor because, other
            // split criteria need a sorting of the samples according to one dimension
            // as well.
            IndexIterator  curSplit = bgfunc.operator()(columnVector(features, splitColumns[k]), labels, region);

            if(bgfunc.minGini < minGini &&  !closeAtTolerance(bgfunc.minGini, minGini))
            {
                std::copy(bgfunc.bestCurrentCounts[0].begin(), bgfunc.bestCurrentCounts[0].end(),
                                                                childRegions[0].classCounts().begin());
                std::copy(bgfunc.bestCurrentCounts[1].begin(), bgfunc.bestCurrentCounts[1].end(),
                                                                childRegions[1].classCounts().begin());
//                childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
//                childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];

                minGini         = bgfunc.minGini;
                bestSplit       = curSplit;
                bestSplitColumn   = splitColumns[k];

                childRegions[0].classCountsIsValid = true;
                childRegions[1].classCountsIsValid = true;
                num2try = std::max(int(k), SB::opt_.mtry);
            }
        }

        if(minGini == NumericTraits<double>::max())
            return  makeTerminalNode(features, labels, region, randint);

        // sort samples according to best split dimension.
        sorter.setColumn(bestSplitColumn);
        std::sort(region.begin(), region.end(), sorter);

        // check informativeness of best Split
        double totalCountss
            = std::inner_product(region.classCounts().begin(),
                                 region.classCounts.end(),
                                 SB::opt_.classWeights.begin(),
                                 0.0);
        if((calculateGini(region, SB::opt_.classWeights) - minGini)/totalCountss < 0.05)
        {
            bestSplit = region.begin() + int(region.size()/2);
        }
        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        //fix threshold value and update childRegions
        node.threshold() = (    features(bestSplit[0], bestSplitColumn)
                            +   features(bestSplit[1], bestSplitColumn)) / 2.0;
        node.column() = bestSplitColumn;

        ++bestSplit;

        // Save the childRanges.
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[1].setRange(   bestSplit       , region.end()    );
        return i_ThresholdNode;
    }
};

template<class Options>
class TotallyUninformedSplit: public SplitBase<Options>
{
  public:


    typedef SplitBase<Options> SB;


    ArrayVector<Int32> splitColumns;
    BestGiniFunctor bgfunc;
    double minGini;
    int bestSplitColumn;
    TotallyUninformedSplit(int featureCount_, int classCount_, Options opt)
        :   SplitBase<Options>(featureCount_, classCount_, opt),
            bgfunc(classCount_, SB::opt_.classWeights)
    {
        if(SB::opt_.mtry == 0)
            SB::opt_.mtry = int(std::floor(std::sqrt(double(featureCount_)) + 0.5));
        SB::opt_.classWeights.resize(classCount_, 1.0);
        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
    }


    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

        typedef typename Region::IndexIterator IndexIterator;
            // if the region is pure than just make a terminal Node.
        if(!region.classCountsIsValid)
        {
            RandomForestClassCounter<T2, C2, ArrayVector<Int32> >
                                    counter(labels,region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }
        if(region.isPure())
            return  makeTerminalNode(features, labels, region, randint);

            // select columns to be tried for split first mtry entries in
            // splitColumns
        Sampler<Random> sampler(    SB::opt_.mtry,
                                    columnCount(features),
                                    SamplingOptions()
                                        .sampleWithoutReplacement()
                                        .useExternalMemory(splitColumns),
                                    randint);


        sampler.sample();
            // This Thing sorts a vector by a given dimension in the multi array
            //features

        SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorter(features, 0);

            // find the best gini index
        bestSplitColumn = 0;
        minGini = NumericTraits<double>::max();
        IndexIterator           bestSplit = IndexIterator();
        //ArrayVector<double> &   classWeights = SB::opt_.classWeights;
        int num2try = features.shape(1);

        for(int k=0; k<1; ++k)
        {

            sorter.setColumn(splitColumns[k]);
            bestSplitColumn = splitColumns[k];
            std::sort(region.begin(), region.end(), sorter);
            bestSplit = region.begin() + int(region.size()/2);
        }

        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        //fix threshold value and update childRegions
        node.threshold() = (    features(bestSplit[0], bestSplitColumn)
                            +   features(bestSplit[1], bestSplitColumn)) / 2.0;
        node.column() = bestSplitColumn;

        ++bestSplit;

        // Save the childRanges.
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[1].setRange(   bestSplit       , region.end()    );
        return i_ThresholdNode;
    }
};



template<class Options>
class HyperplaneSplit: public SplitBase<Options>
{
  public:


    typedef SplitBase<Options> SB;


    ArrayVector<Int32> splitColumns;
    BestGiniFunctor bgfunc;
    double minGini;
    int bestSplitColumn;
    HyperplaneSplit(int featureCount_, int classCount_, Options opt)
        :   SplitBase<Options>(featureCount_, classCount_, opt),
            bgfunc(classCount_, SB::opt_.classWeights)
    {
        if(SB::opt_.mtry == 0)
            SB::opt_.mtry = int(std::floor(std::sqrt(double(featureCount_)) + 0.5));
        SB::opt_.classWeights.resize(classCount_, 1.0);
        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
    }


    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  labels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

        typedef typename Region::IndexIterator IndexIterator;
            // if the region is pure than just make a terminal Node.
        if(!region.classCountsIsValid)
        {
            RandomForestClassCounter<T2, C2, ArrayVector<Int32> >
                                    counter(labels,region.classCounts());
            std::for_each(  region.begin(), region.end(), counter);
            region.classCountsIsValid = true;
        }
        if(region.isPure())
            return  makeTerminalNode(features, labels, region, randint);

            // select columns to be tried for split first mtry entries in
            // splitColumns
        Sampler<Random> sampler(    SB::opt_.mtry,
                                    columnCount(features),
                                    SamplingOptions()
                                        .sampleWithoutReplacement()
                                        .useExternalMemory(splitColumns),
                                    randint);


        sampler.sample();
            // This Thing sorts a vector by a given dimension in the multi array
            //features

        SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorter(features, 0);

            // find the best gini index
        bestSplitColumn = 0;
        minGini = NumericTraits<double>::max();
        IndexIterator           bestSplit = IndexIterator();
        //ArrayVector<double> &   classWeights = SB::opt_.classWeights;
        int num2try = features.shape(1);
        ArrayVectorView<Int32>  	labelView(region.size(), &region[0]);

        MultiArray<2, T> ww(MultiArrayShape<2>::type(features.shape(1)+1,1));
        int num = HoKashyap( features,
                                labels,
                                ww,
                                1,
                                detail::local_feats,
                                labelView);

        if(num < 1000)
        {
            Node<i_HyperplaneNode>   node(features.shape(1), SB::t_data, SB::p_data);
            std::copy(ww.begin(), ww.end(), node.parameters_begin()+1);
            SB::node_ = node;
            SortSamplesByDimensions<MultiArrayView<2, T, C> >
                                sorterHyp(features, node);
            std::sort(region.begin(), region.end(), sorterHyp);
            bestSplit = region.begin();
            int ii = 0;
            while(sorterHyp[region[ii]] < 0)
                ++ii;
            bestSplit += ii;
            childRegions[0].setRange(   region.begin()  , bestSplit       );
            childRegions[1].setRange(   bestSplit       , region.end()    );
            return i_HyperplaneNode;
        }

        for(int k=0; k<num2try; ++k)
        {

            sorter.setColumn(splitColumns[k]);
            std::sort(region.begin(), region.end(), sorter);
            // The calculation of the gini was put into the functor because, other
            // split criteria need a sorting of the samples according to one dimension
            // as well.
            IndexIterator  curSplit = bgfunc.operator()(columnVector(features, splitColumns[k]), labels, region);

            if(bgfunc.minGini < minGini &&  !closeAtTolerance(bgfunc.minGini, minGini))
            {
                std::copy(bgfunc.bestCurrentCounts[0].begin(), bgfunc.bestCurrentCounts[0].end(),
                                                                childRegions[0].classCounts().begin());
                std::copy(bgfunc.bestCurrentCounts[1].begin(), bgfunc.bestCurrentCounts[1].end(),
                                                                childRegions[1].classCounts().begin());
//                childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
//                childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];

                minGini         = bgfunc.minGini;
                bestSplit       = curSplit;
                bestSplitColumn   = splitColumns[k];

                childRegions[0].classCountsIsValid = true;
                childRegions[1].classCountsIsValid = true;
                num2try = std::max(int(k), SB::opt_.mtry);
            }
        }

        if(minGini == NumericTraits<double>::max())
            return  makeTerminalNode(features, labels, region, randint);

        // sort samples according to best split dimension.
        sorter.setColumn(bestSplitColumn);
        std::sort(region.begin(), region.end(), sorter);


        Node<i_ThresholdNode>   node(SB::t_data, SB::p_data);
        SB::node_ = node;
        //fix threshold value and update childRegions
        node.threshold() = (    features(bestSplit[0], bestSplitColumn)
                            +   features(bestSplit[1], bestSplitColumn)) / 2.0;
        node.column() = bestSplitColumn;

        ++bestSplit;

        // Save the childRanges.
        childRegions[0].setRange(   region.begin()  , bestSplit       );
        childRegions[1].setRange(   bestSplit       , region.end()    );
        return i_ThresholdNode;
    }
};
*/
} //namespace vigra
#endif // VIGRA_RANDOM_FOREST_SPLIT_HXX
