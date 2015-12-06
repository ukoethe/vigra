/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2012 by Ullrich Koethe                  */
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

#ifndef VIGRA_ACCUMULATOR_HXX
#define VIGRA_ACCUMULATOR_HXX

#ifdef _MSC_VER
#pragma warning (disable: 4503)
#endif

#include "accumulator-grammar.hxx"
#include "config.hxx"
#include "metaprogramming.hxx"
#include "bit_array.hxx"
#include "static_assert.hxx"
#include "mathutil.hxx"
#include "utilities.hxx"
#include "multi_iterator_coupled.hxx"
#include "matrix.hxx"
#include "multi_math.hxx"
#include "eigensystem.hxx"
#include "histogram.hxx"
#include "polygon.hxx"
#include "functorexpression.hxx"
#include "labelimage.hxx"
#include <algorithm>
#include <iostream>

namespace vigra {
  
/** \defgroup FeatureAccumulators Feature Accumulators

The namespace <tt>vigra::acc</tt> provides the function \ref vigra::acc::extractFeatures() along with associated statistics functors and accumulator classes. Together, they provide a framework for efficient compution of a wide variety of statistical features, both globally for an entire image, and locally for each region defined by a label array. Many different statistics can be composed out of a small number of fundamental statistics and suitable modifiers. The user simply selects the desired statistics by means of their <i>tags</i> (see below), and a template meta-program automatically generates an efficient functor that computes exactly those statistics.

The function \ref acc::extractFeatures() "extractFeatures()" scans the data in as few passes as the selected statstics permit (usually one or two passes are sufficient). Statistics are computed by accurate incremental algorithms, whose internal state is maintained by accumulator objects. The state is updated by passing data to the accumulator one sample at a time. Accumulators are grouped within an accumulator chain. Dependencies between accumulators in the accumulator chain are automatically resolved and missing dependencies are inserted. For example, to compute the mean, you also need to count the number of samples. This allows accumulators to offload some of their computations on other accumulators, making the algorithms more efficient. Each accumulator only sees data in the appropriate pass through the data, called its "working pass". 

<b>\#include</b> \<vigra/accumulator.hxx\>

    
<b>Basic statistics:</b>
    - PowerSum<N> (computes @f$ \sum_i x_i^N @f$)
    - AbsPowerSum<N> (computes @f$ \sum_i |x_i|^N @f$)
    - Skewness, UnbiasedSkewness
    - Kurtosis, UnbiasedKurtosis
    - Minimum, Maximum
    - FlatScatterMatrix (flattened upper-triangular part of scatter matrix)
    - 4 histogram classes (see \ref histogram "below")
    - StandardQuantiles (0%, 10%, 25%, 50%, 75%, 90%, 100%)
    - ArgMinWeight, ArgMaxWeight (store data or coordinate where weight assumes its minimal or maximal value)
    - CoordinateSystem (identity matrix of appropriate size)
    
    <b>Modifiers:</b> (S is the statistc to be modified)
    - Normalization
      <table border="0">
      <tr><td> DivideByCount<S>        </td><td>  S/Count           </td></tr>
      <tr><td> RootDivideByCount<S>    </td><td>  sqrt( S/Count )     </td></tr>
      <tr><td> DivideUnbiased<S>       </td><td>  S/(Count-1)       </td></tr>
      <tr><td> RootDivideUnbiased<S> &nbsp; &nbsp;  </td><td>  sqrt( S/(Count-1) ) </td></tr>
      </table>    
    - Data preparation:
      <table border="0">
      <tr><td>  Central<S>   </td><td> substract mean before computing S </td></tr>
      <tr><td>  Principal<S> </td><td> project onto PCA eigenvectors   </td></tr>
      <tr><td>  Whitened<S> &nbsp; &nbsp;  </td><td> scale to unit variance after PCA   </td></tr>
      <tr><td>  Coord<S>        </td><td> compute S from pixel coordinates rather than from pixel values    </td></tr>
      <tr><td>  Weighted<S>     </td><td> compute weighted version of S   </td></tr>
      <tr><td>  Global<S>       </td><td> compute S globally rather than per region (per region is default if labels are given)   </td></tr>
      </table>
      
    Aliases for many important features are implemented (mainly as <tt>typedef FullName Alias</tt>). The alias names are equivalent to full names. Below are some examples for supported alias names. A full list of all available statistics and alias names can be found in the namespace reference <tt>vigra::acc</tt>. These examples also show how to compose statistics from the fundamental statistics and modifiers:
    
    <table border="0">
    <tr><th> Alias           </th><th>   Full Name                 </th></tr>
    <tr><td> Count           </td><td>  PowerSum<0>                </td></tr>
    <tr><td> Sum             </td><td>  PowerSum<1>                </td></tr>
    <tr><td> SumOfSquares    </td><td>  PowerSum<2>                </td></tr>
    <tr><td> Mean            </td><td>  DivideByCount<PowerSum<1>> </td></tr>
    <tr><td> RootMeanSquares &nbsp; </td><td>  RootDivideByCount<PowerSum<2>> </td></tr>
    <tr><td> Moment<N>       </td><td>  DivideByCount<PowerSum<N>>  </td></tr>
    <tr><td> Variance        </td><td>  DivideByCount<Central<PowerSum<2>>>  </td></tr>
    <tr><td> StdDev          </td><td>  RootDivideByCount<Central<PowerSum<2>>>  </td></tr>
    <tr><td> Covariance      </td><td>  DivideByCount<FlatScatterMatrix> </td></tr>
    <tr><td> RegionCenter    </td><td>  Coord<Mean>                </td></tr>
    <tr><td> CenterOfMass    </td><td>  Weighted<Coord<Mean>>      </td></tr>
    </table>
    
    There are a few <b>rules for composing statistics</b>:
    - modifiers can be specified in any order, but are internally transformed to standard order: Global<Weighted<Coord<normalization<data preparation<basic statistic
    - only one normalization modifier and one data preparation modifier (Central or Principal or Whitened) is permitted 
    - Count ignores all modifiers except Global and Weighted
    - Sum ignores Central and Principal, because sum would be zero
    - ArgMinWeight and ArgMaxWeight are automatically Weighted


    Here is an example how to use \ref acc::AccumulatorChain to compute statistics. (To use Weighted<> or Coord<> modifiers, see below):

    \code
    #include <vigra/multi_array.hxx>
    #include <vigra/impex.hxx>
    #include <vigra/accumulator.hxx>
    using namespace vigra::acc;
    typedef double DataType;
    int size = 1000;
    vigra::MultiArray<2, DataType> data(vigra::Shape2(size, size));
   
    AccumulatorChain<DataType, 
        Select<Variance, Mean, StdDev, Minimum, Maximum, RootMeanSquares, Skewness, Covariance> >
        a;
        
    std::cout << "passes required: " << a.passesRequired() << std::endl;
    extractFeatures(data.begin(), data.end(), a); 
    
    std::cout << "Mean: " << get<Mean>(a) << std::endl;
    std::cout << "Variance: " << get<Variance>(a) << std::endl;
    \endcode
    
    The \ref acc::AccumulatorChain object contains the selected statistics and their dependencies. Statistics have to be wrapped with \ref acc::Select. The statistics are computed with the acc::extractFeatures function and the statistics can be accessed with acc::get . 

    Rules and notes:
    - order of statistics in Select<> is arbitrary
    - up to 20 statistics in Select<>, but Select<> can be nested
    - dependencies are automatically inserted
    - duplicates are automatically removed
    - extractFeatures() does as many passes through the data as necessary
    - each accumulator only sees data in the appropriate pass (its "working pass")

    The Accumulators can also be used with vector-valued data (vigra::RGBValue, vigra::TinyVector, vigra::MultiArray or vigra::MultiArrayView):
    
    \code
    typedef vigra::RGBValue<double> DataType;
    AccumulatorChain<DataType, Select<...> > a;
    ...
    \endcode

    To compute <b>weighted statistics</b> (Weighted<>) or <b>statistics over coordinates</b> (Coord<>), the accumulator chain can be used with several coupled arrays, one for the data and another for the weights and/or the labels. "Coupled" means that statistics are computed over the corresponding elements of the involved arrays. This is internally done by means of \ref CoupledScanOrderIterator and \ref vigra::CoupledHandle which provide simultaneous access to several arrays (e.g. weight and data) and corresponding coordinates. The types of the coupled arrays are best specified by means of the helper class \ref vigra::CoupledArrays :
    
    \code 
    vigra::MultiArray<3, RGBValue<unsigned char> > data(...);
    vigra::MultiArray<3, double>                   weights(...);
    
    AccumulatorChain<CoupledArrays<3, RGBValue<unsigned char>, double>,
                     Select<...> > a;
    \endcode
    
This works likewise for label images which are needed for region statistics (see below). The indxx of the array holding data, weights, or labels respectively can be specified inside the Select wrapper. These <b>index specifiers</b> are: (INDEX is of type int)
    - DataArg<INDEX>: data are in array 'INDEX' (default INDEX=1)
    - LabelArg<INDEX>: labels are in array 'INDEX' (default INDEX=2)
    - WeightArg<INDEX>: weights are in array 'INDEX' (default INDEX=rightmost index)

Pixel coordinates are always at index 0. To collect statistics, you simply pass all arrays to the <tt>extractFeatures()</tt> function:
    \code
    using namespace vigra::acc;
    vigra::MultiArray<3, double> data(...), weights(...);
    
    AccumulatorChain<CoupledArrays<3, double, double>, // two 3D arrays for data and weights
        Select<DataArg<1>, WeightArg<2>,           // in which array to look (coordinates are always arg 0)
               Mean, Variance,                     //statistics over values  
               Coord<Mean>, Coord<Variance>,       //statistics over coordinates,
               Weighted<Mean>, Weighted<Variance>, //weighted values,
               Weighted<Coord<Mean> > > >          //weighted coordinates.
        a;
     
    extractFeatures(data, weights, a);
    \endcode
    
    This even works for a single array, which is useful if you want to combine values with coordinates. For example, to find the location of the minimum element in an array, you interpret the data as weights and select the <tt>Coord<ArgMinWeight></tt> statistic (note that the version of <tt>extractFeatures()</tt> below only works in conjunction with <tt>CoupledArrays</tt>, despite the fact that there is only one array involved):
    \code 
    using namespace vigra::acc;
    vigra::MultiArray<3, double> data(...);
    
    AccumulatorChain<CoupledArrays<3, double>,
                     Select<WeightArg<1>,           // we interprete the data as weights
                            Coord<ArgMinWeight> > > // and look for the coordinate with minimal weight
        a;
        
    extractFeatures(data, a);
    std::cout << "minimum is at " << get<Coord<ArgMinWeight> >(a) << std::endl;
    \endcode
    
    To compute <b>region statistics</b>, you use \ref acc::AccumulatorChainArray. Regions are defined by means of a label array whose elements specify the region ID of the corresponding point. Therefore, you will always need at least two arrays here, which are again best specified using the <tt>CoupledArrays</tt> helper:
    
    \code
    using namespace vigra::acc;
    vigra::MultiArray<3, double> data(...);
    vigra::MultiArray<3, int> labels(...);

    AccumulatorChainArray<CoupledArrays<3, double, int>,
        Select<DataArg<1>, LabelArg<2>,       // in which array to look (coordinates are always arg 0)
               Mean, Variance,                    //per-region statistics over values
               Coord<Mean>, Coord<Variance>,      //per-region statistics over coordinates
               Global<Mean>, Global<Variance> > > //global statistics
    a;

    a.ignoreLabel(0); //statistics will not be computed for region 0 (e.g. background)

    extractFeatures(data, labels, a);

    int regionlabel = ...;
    std::cout << get<Mean>(a, regionlabel) << std::endl; //get Mean of region with label 'regionlabel'
    \endcode

   
    In some application it will be known only at run-time which statistics have to be computed. An Accumulator with <b>run-time activation</b> is provided by the \ref acc::DynamicAccumulatorChain class. One specifies a set of statistics at compile-time and from this set one can activate the needed statistics at run-time:
  
    \code
    using namespace vigra::acc;
    vigra::MultiArray<2, double> data(...);
    DynamicAccumulatorChain<double, 
        Select<Mean, Minimum, Maximum, Variance, StdDev> > a; // at compile-time
    activate<Mean>(a);      //at run-time
    a.activate("Minimum");  //same as activate<Minimum>(a) (alias names are not recognized)
    
    extractFeatures(data.begin(), data.end(), a);
    std::cout << "Mean: " << get<Mean>(a) << std::endl;       //ok
    //std::cout << "Maximum: " << get<Maximum>(a) << std::endl; // run-time error because Maximum not activated
    \endcode
      
    Likewise, for run-time activation of region statistics, use \ref acc::DynamicAccumulatorChainArray. 

    <b>Accumulator merging</b> (e.g. for parallelization or hierarchical segmentation) is possible for many accumulators:

    \code
    using namespace vigra::acc;
    vigra::MultiArray<2, double> data(...);
    AccumulatorChain<double, Select<Mean, Variance, Skewness> > a, a1, a2;

    extractFeatures(data.begin(), data.end(), a); //process entire data set at once
    extractFeatures(data.begin(), data.begin()+data.size()/2, a1); //process first half
    extractFeatures(data.begin()+data.size()/2, data.end(), a2); //process second half
    a1 += a2; // merge: a1 now equals a0 (within numerical tolerances)
    \endcode

    Not all statistics can be merged (e.g. Principal<A> usually cannot, except for some important specializations). A statistic can be merged if the "+=" operator is supported (see the documentation of that particular statistic). If the accumulator chain only requires one pass to collect the data, it is also possible to just apply the extractFeatures() function repeatedly:

    \code
    using namespace vigra::acc;
    vigra::MultiArray<2, double> data(...);
    AccumulatorChain<double, Select<Mean, Variance> > a;

    extractFeatures(data.begin(), data.begin()+data.size()/2, a); // this works because 
    extractFeatures(data.begin()+data.size()/2, data.end(), a);   // all statistics only need pass 1
    \endcode
    
    More care is needed to merge coordinate-based statistics. By default, all coordinate statistics are computed in the local coordinate system of the current region of interest. That is, the upper left corner of the ROI has the coordinate (0, 0) by default. This behavior is not desirable when you want to merge coordinate statistics from different ROIs: then, all accumulators should use the same coordinate system, usually the global system of the entire dataset. This can be achieved by the <tt>setCoordinateOffset()</tt> function. The following code demonstrates this for the <tt>RegionCenter</tt> statistic:

    \code
    using namespace vigra;
    using namespace vigra::acc;
    
    MultiArray<2, double> data(width, height);
    MultiArray<2, int>    labels(width, height);
    
    AccumulatorChainArray<CoupledArrays<2, double, int>,
                          Select<DataArg<1>, LabelArg<2>, 
                                 RegionCenter> >
    a1, a2;
    
    // a1 is responsible for the left half of the image. The local coordinate system of this ROI 
    // happens to be identical to the global coordinate system, so the offset is zero.
    Shape2 origin(0,0);
    a1.setCoordinateOffset(origin);
    extractFeatures(data.subarray(origin, Shape2(width/2, height)), 
                    labels.subarray(origin, Shape2(width/2, height)),
                    a1);
                    
    // a2 is responsible for the right half, so the offset of the local coordinate system is (width/2, 0)
    origin = Shape2(width/2, 0);
    a2.setCoordinateOffset(origin);
    extractFeatures(data.subarray(origin, Shape2(width, height)), 
                    labels.subarray(origin, Shape2(width, height)),
                    a2);
   
    // since both accumulators worked in the same global coordinate system, we can safely merge them
    a1.merge(a2);
    \endcode
    
    When you compute region statistics in ROIs, it is sometimes desirable to use a local region labeling in each ROI. In this way, the labels of each ROI cover a consecutive range of numbers starting with 0. This can save a lot of memory, because <tt>AccumulatorChainArray</tt> internally uses dense arrays -- accumulators will be allocated for all labels from 0 to the maxmimum label, even when many of them are unused. This is avoided by a local labeling. However, this means that label 1 (say) may refer to two different regions in different ROIs. To adjust for this mismatch, you can pass a label mapping to <tt>merge()</tt> that provides a global label for each label of the accumulator to be merged. Thus, each region on the right hand side will be merged into the left-hand-side accumulator with the given <i>global</i> label. For example, let us assume that the left and right half of the image contain just one region and background. Then, the accumulators of both ROIs have the label 0 (background) and 1 (the region). Upon merging, the region from the right ROI should be given the global label 2, whereas the background should keep its label 0. This is achieved like this:
    
    \code
    std::vector<int> labelMapping(2);
    labelMapping[0] = 0;   // background keeps label 0
    labelMapping[1] = 2;   // local region 1 becomes global region 2
    
    a1.merge(a2, labelMapping);
    \endcode

    \anchor histogram
    Four kinds of <b>histograms</b> are currently implemented:
    
    <table border="0">
      <tr><td> IntegerHistogram      </td><td>   Data values are equal to bin indices   </td></tr>
      <tr><td> UserRangeHistogram    </td><td>  User provides lower and upper bounds for linear range mapping from values to indices.    </td></tr>
      <tr><td> AutoRangeHistogram    </td><td>  Range mapping bounds are defiend by minimum and maximum of the data (2 passes needed!)    </td></tr>
      <tr><td> GlobalRangeHistogram &nbsp;  </td><td>  Likewise, but use global min/max rather than region min/max as AutoRangeHistogram will </td></tr>
      </table>    
  

       
    - The number of bins is specified at compile time (as template parameter int BinCount) or at run-time (if BinCount is zero at compile time). In the first case the return type of the accumulator is TinyVector<double, BinCount> (number of bins cannot be changed). In the second case, the return type is MultiArray<1, double> and the number of bins must be set before seeing data (see example below). 
    - If UserRangeHistogram is used, the bounds for the linear range mapping from values to indices must be set before seeing data (see below).
    - Options can be set by passing an instance of HistogramOptions to the accumulator chain (same options for all histograms in the chain) or by directly calling the appropriate member functions of the accumulators.
    - Merging is supported if the range mapping of the histograms is the same.
    - Histogram accumulators have two members for outliers (left_outliers, right_outliers).

    With the StandardQuantiles class, <b>histogram quantiles</b> (0%, 10%, 25%, 50%, 75%, 90%, 100%) are computed from a given histgram using linear interpolation. The return type is TinyVector<double, 7> .

    \anchor acc_hist_options Usage:
    \code
    using namespace vigra::acc;
    typedef double DataType;
    vigra::MultiArray<2, DataType> data(...);
    
    typedef UserRangeHistogram<40> SomeHistogram;   //binCount set at compile time
    typedef UserRangeHistogram<0> SomeHistogram2; // binCount must be set at run-time
    typedef AutoRangeHistogram<0> SomeHistogram3;
    typedef StandardQuantiles<SomeHistogram3> Quantiles3;
    
    AccumulatorChain<DataType, Select<SomeHistogram, SomeHistogram2, SomeHistogram3, Quantiles3> > a;
    
    //set options for all histograms in the accumulator chain:
    vigra::HistogramOptions histogram_opt;         
    histogram_opt = histogram_opt.setBinCount(50); 
    //histogram_opt = histogram_opt.setMinMax(0.1, 0.9); // this would set min/max for all three histograms, but range bounds 
                                                         // shall be set automatically by min/max of data for SomeHistogram3
    a.setHistogramOptions(histogram_opt);  

    // set options for a specific histogram in the accumulator chain:
    getAccumulator<SomeHistogram>(a).setMinMax(0.1, 0.9); // number of bins must be set before setting min/max
    getAccumulator<SomeHistogram2>(a).setMinMax(0.0, 1.0);

    extractFeatures(data.begin(), data.end(), a);

    vigra::TinyVector<double, 40> hist = get<SomeHistogram>(a);
    vigra::MultiArray<1, double> hist2 = get<SomeHistogram2>(a);
    vigra::TinyVector<double, 7> quant = get<Quantiles3>(a);
    double right_outliers = getAccumulator<SomeHistogram>(a).right_outliers;
    \endcode


    
*/


/** This namespace contains the accumulator classes, fundamental statistics and modifiers. See \ref FeatureAccumulators for examples of usage.
*/
namespace acc {

/****************************************************************************/
/*                                                                          */
/*                             infrastructure                               */
/*                                                                          */
/****************************************************************************/

  /// \brief Wrapper for MakeTypeList that additionally performs tag standardization.

template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct Select
: public MakeTypeList<
    typename StandardizeTag<T01>::type, typename StandardizeTag<T02>::type, typename StandardizeTag<T03>::type, 
    typename StandardizeTag<T04>::type, typename StandardizeTag<T05>::type, typename StandardizeTag<T06>::type, 
    typename StandardizeTag<T07>::type, typename StandardizeTag<T08>::type, typename StandardizeTag<T09>::type, 
    typename StandardizeTag<T10>::type, typename StandardizeTag<T11>::type, typename StandardizeTag<T12>::type, 
    typename StandardizeTag<T13>::type, typename StandardizeTag<T14>::type, typename StandardizeTag<T15>::type, 
    typename StandardizeTag<T16>::type, typename StandardizeTag<T17>::type, typename StandardizeTag<T18>::type, 
    typename StandardizeTag<T19>::type, typename StandardizeTag<T20>::type
    >
{};

    // enable nesting of Select<> expressions 
template <class T01, class T02, class T03, class T04, class T05,
          class T06, class T07, class T08, class T09, class T10,
          class T11, class T12, class T13, class T14, class T15,
          class T16, class T17, class T18, class T19, class T20>
struct StandardizeTag<Select<T01, T02, T03, T04, T05,
                             T06, T07, T08, T09, T10,
                             T11, T12, T13, T14, T15,
                             T16, T17, T18, T19, T20>, 
                      Select<T01, T02, T03, T04, T05,
                             T06, T07, T08, T09, T10,
                             T11, T12, T13, T14, T15,
                             T16, T17, T18, T19, T20> >
{
    typedef typename  Select<T01, T02, T03, T04, T05,
                             T06, T07, T08, T09, T10,
                             T11, T12, T13, T14, T15,
                             T16, T17, T18, T19, T20>::type type;
};

struct AccumulatorBegin
{
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "AccumulatorBegin (internal)";
       // static const std::string n("AccumulatorBegin (internal)");
       // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {};
};


struct AccumulatorEnd;
struct DataArgTag;
struct WeightArgTag;
struct LabelArgTag;
struct CoordArgTag;
struct LabelDispatchTag;

template <class T, class TAG, class CHAIN>
struct HandleArgSelector;  // find the correct handle in a CoupledHandle

struct Error__Global_statistics_are_only_defined_for_AccumulatorChainArray;

/** \brief Specifies index of labels in CoupledHandle. 

    LabelArg<INDEX> tells the acc::AccumulatorChainArray which index of the Handle contains the labels. (Note that coordinates are always index 0)
 */
template <int INDEX>
class LabelArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("LabelArg<") + asString(INDEX) + "> (internal)";
        // static const std::string n = std::string("LabelArg<") + asString(INDEX) + "> (internal)";
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef LabelArgTag Tag;
        typedef void value_type;
        typedef void result_type;

        static const int value = INDEX;
        static const unsigned int workInPass = 0;
    };
};

template <int INDEX>
class CoordArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("CoordArg<") + asString(INDEX) + "> (internal)";
        // static const std::string n = std::string("CoordArg<") + asString(INDEX) + "> (internal)";
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef CoordArgTag Tag;
        typedef void value_type;
        typedef void result_type;

        static const int value = INDEX;
        static const unsigned int workInPass = 0;
    };
};

template <class T, class TAG, class NEXT=AccumulatorEnd>
struct AccumulatorBase;

template <class Tag, class A>
struct LookupTag;

template <class Tag, class A, class TargetTag=typename A::Tag>
struct LookupDependency;

#ifndef _MSC_VER  // compiler bug? (causes 'ambiguous overload error')

template <class TAG, class A>
typename LookupTag<TAG, A>::reference
getAccumulator(A & a);

template <class TAG, class A>
typename LookupDependency<TAG, A>::result_type
getDependency(A const & a);

#endif

namespace acc_detail {

/****************************************************************************/
/*                                                                          */
/*                   internal tag handling meta-functions                   */
/*                                                                          */
/****************************************************************************/

    // we must make sure that Arg<INDEX> tags are at the end of the chain because 
    // all other tags potentially depend on them
template <class T>
struct PushArgTagToTail
{
    typedef T type;
};

#define VIGRA_PUSHARGTAG(TAG) \
template <int INDEX, class TAIL> \
struct PushArgTagToTail<TypeList<TAG<INDEX>, TAIL> > \
{ \
    typedef typename Push<TAIL, TypeList<TAG<INDEX> > >::type type; \
};

VIGRA_PUSHARGTAG(DataArg)
VIGRA_PUSHARGTAG(WeightArg)
VIGRA_PUSHARGTAG(CoordArg)
VIGRA_PUSHARGTAG(LabelArg)

#undef VIGRA_PUSHARGTAG

    // Insert the dependencies of the selected functors into the TypeList and sort
    // the list such that dependencies come after the functors using them. Make sure 
    // that each functor is contained only once.
template <class T>
struct AddDependencies;

template <class HEAD, class TAIL>
struct AddDependencies<TypeList<HEAD, TAIL> >
{
    typedef typename AddDependencies<TAIL>::type                                   TailWithDependencies;
    typedef typename StandardizeDependencies<HEAD>::type                           HeadDependencies;
    typedef typename AddDependencies<HeadDependencies>::type                       TransitiveHeadDependencies;
    typedef TypeList<HEAD, TransitiveHeadDependencies>                             HeadWithDependencies;
    typedef typename PushUnique<HeadWithDependencies, TailWithDependencies>::type  UnsortedDependencies;
    typedef typename PushArgTagToTail<UnsortedDependencies>::type                  type;
};

template <>
struct AddDependencies<void>
{
    typedef void type;
};

    // Helper class to activate dependencies at runtime (i.e. when activate<Tag>(accu) is called,
    // activate() must also be called for Tag's dependencies).
template <class Dependencies>
struct ActivateDependencies;

template <class HEAD, class TAIL>
struct ActivateDependencies<TypeList<HEAD, TAIL> >
{
    template <class Chain, class ActiveFlags>
    static void exec(ActiveFlags & flags)
    {
        LookupTag<HEAD, Chain>::type::activateImpl(flags);
        ActivateDependencies<TAIL>::template exec<Chain>(flags);
    }
    
    template <class Chain, class ActiveFlags, class GlobalFlags>
    static void exec(ActiveFlags & flags, GlobalFlags & gflags)
    {
        LookupTag<HEAD, Chain>::type::template activateImpl<Chain>(flags, gflags);
        ActivateDependencies<TAIL>::template exec<Chain>(flags, gflags);
    }
};

template <class HEAD, class TAIL>
struct ActivateDependencies<TypeList<Global<HEAD>, TAIL> >
{
    template <class Chain, class ActiveFlags, class GlobalFlags>
    static void exec(ActiveFlags & flags, GlobalFlags & gflags)
    {
        LookupTag<Global<HEAD>, Chain>::type::activateImpl(gflags);
        ActivateDependencies<TAIL>::template exec<Chain>(flags, gflags);
    }
};

template <>
struct ActivateDependencies<void>
{
    template <class Chain, class ActiveFlags>
    static void exec(ActiveFlags &)
    {}
    
    template <class Chain, class ActiveFlags, class GlobalFlags>
    static void exec(ActiveFlags &, GlobalFlags &)
    {}
};

template <class List>
struct SeparateGlobalAndRegionTags;

template <class HEAD, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<HEAD, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef TypeList<HEAD, typename Inner::RegionTags>  RegionTags;
    typedef typename Inner::GlobalTags                  GlobalTags;
};

template <class HEAD, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<Global<HEAD>, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef typename Inner::RegionTags                  RegionTags;
    typedef TypeList<HEAD, typename Inner::GlobalTags>  GlobalTags;
};

template <int INDEX, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<DataArg<INDEX>, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef TypeList<DataArg<INDEX>, typename Inner::RegionTags>  RegionTags;
    typedef TypeList<DataArg<INDEX>, typename Inner::GlobalTags>  GlobalTags;
};

template <int INDEX, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<LabelArg<INDEX>, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef TypeList<LabelArg<INDEX>, typename Inner::RegionTags>  RegionTags;
    typedef TypeList<LabelArg<INDEX>, typename Inner::GlobalTags>  GlobalTags;
};

template <int INDEX, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<WeightArg<INDEX>, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef TypeList<WeightArg<INDEX>, typename Inner::RegionTags>  RegionTags;
    typedef TypeList<WeightArg<INDEX>, typename Inner::GlobalTags>  GlobalTags;
};

template <int INDEX, class TAIL>
struct SeparateGlobalAndRegionTags<TypeList<CoordArg<INDEX>, TAIL> >
{
    typedef SeparateGlobalAndRegionTags<TAIL>           Inner;
    typedef TypeList<CoordArg<INDEX>, typename Inner::RegionTags>  RegionTags;
    typedef TypeList<CoordArg<INDEX>, typename Inner::GlobalTags>  GlobalTags;
};

template <>
struct SeparateGlobalAndRegionTags<void>
{
    typedef void RegionTags;
    typedef void GlobalTags;
};

/****************************************************************************/
/*                                                                          */
/*          helper classes to handle tags at runtime via strings            */
/*                                                                          */
/****************************************************************************/

template <class Accumulators>
struct CollectAccumulatorNames;

template <class HEAD, class TAIL>
struct CollectAccumulatorNames<TypeList<HEAD, TAIL> >
{
    template <class BackInsertable>
    static void exec(BackInsertable & a, bool skipInternals=true)
    {
        if(!skipInternals || HEAD::name().find("internal") == std::string::npos)
            a.push_back(HEAD::name());
        CollectAccumulatorNames<TAIL>::exec(a, skipInternals);
    }
};

template <>
struct CollectAccumulatorNames<void>
{
    template <class BackInsertable>
    static void exec(BackInsertable & a, bool skipInternals=true)
    {}
};

template <class T>
struct ApplyVisitorToTag;

template <class HEAD, class TAIL>
struct ApplyVisitorToTag<TypeList<HEAD, TAIL> >
{
    template <class Accu, class Visitor>
    static bool exec(Accu & a, std::string const & tag, Visitor const & v)
    {
        static std::string * name = VIGRA_SAFE_STATIC(name, new std::string(normalizeString(HEAD::name())));
        if(*name == tag)
        {
            v.template exec<HEAD>(a);
            return true;
        }
        else
        {
            return ApplyVisitorToTag<TAIL>::exec(a, tag, v);
        }
    }
};

template <>
struct ApplyVisitorToTag<void>
{
    template <class Accu, class Visitor>
    static bool exec(Accu & a, std::string const & tag, Visitor const & v)
    {
        return false;
    }
};

struct ActivateTag_Visitor
{
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        a.template activate<TAG>();
    }
};

struct TagIsActive_Visitor
{
    mutable bool result;
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        result = a.template isActive<TAG>();
    }
};

/****************************************************************************/
/*                                                                          */
/*                    histogram initialization functors                     */
/*                                                                          */
/****************************************************************************/

template <class TAG>
struct SetHistogramBincount
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {}
};

template <template <int> class Histogram>
struct SetHistogramBincount<Histogram<0> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {
        a.setBinCount(options.binCount);
    }
};

template <class TAG>
struct ApplyHistogramOptions
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {}
};

template <class TAG>
struct ApplyHistogramOptions<StandardQuantiles<TAG> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {}
};

template <class TAG, template <class> class MODIFIER>
struct ApplyHistogramOptions<MODIFIER<TAG> >
: public ApplyHistogramOptions<TAG>
{};

template <>
struct ApplyHistogramOptions<IntegerHistogram<0> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {
        SetHistogramBincount<IntegerHistogram<0> >::exec(a, options);
    }
};

template <int BinCount>
struct ApplyHistogramOptions<UserRangeHistogram<BinCount> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {
        SetHistogramBincount<UserRangeHistogram<BinCount> >::exec(a, options);
        if(a.scale_ == 0.0 && options.validMinMax())
            a.setMinMax(options.minimum, options.maximum);
    }
};

template <int BinCount>
struct ApplyHistogramOptions<AutoRangeHistogram<BinCount> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {
        SetHistogramBincount<AutoRangeHistogram<BinCount> >::exec(a, options);
        if(a.scale_ == 0.0 && options.validMinMax())
            a.setMinMax(options.minimum, options.maximum);
    }
};

template <int BinCount>
struct ApplyHistogramOptions<GlobalRangeHistogram<BinCount> >
{
    template <class Accu>
    static void exec(Accu & a, HistogramOptions const & options)
    {
        SetHistogramBincount<GlobalRangeHistogram<BinCount> >::exec(a, options);
        if(a.scale_ == 0.0)
        {
            if(options.validMinMax())
                a.setMinMax(options.minimum, options.maximum);
            else
                a.setRegionAutoInit(options.local_auto_init);
        }
    }
};

/****************************************************************************/
/*                                                                          */
/*                   internal accumulator chain classes                     */
/*                                                                          */
/****************************************************************************/

    // AccumulatorEndImpl has the following functionalities:
    //  * marks end of accumulator chain by the AccumulatorEnd tag
    //  * provides empty implementation of standard accumulator functions
    //  * provides active_accumulators_ flags for run-time activation of dynamic accumulators
    //  * provides is_dirty_ flags for caching accumulators
    //  * hold the GlobalAccumulatorHandle for global accumulator lookup from region accumulators
template <unsigned LEVEL, class GlobalAccumulatorHandle>
struct AccumulatorEndImpl
{
    typedef typename GlobalAccumulatorHandle::type  GlobalAccumulatorType;
    
    typedef AccumulatorEnd     Tag;
    typedef void               value_type;
    typedef bool               result_type;
    typedef BitArray<LEVEL>    AccumulatorFlags;
    
    static const unsigned int  workInPass = 0; 
    static const int           index = -1;
    static const unsigned      level = LEVEL;
    
    AccumulatorFlags            active_accumulators_;
    mutable AccumulatorFlags    is_dirty_;
    GlobalAccumulatorHandle     globalAccumulator_;
        
    template <class GlobalAccumulator>
    void setGlobalAccumulator(GlobalAccumulator const * a)
    {
        globalAccumulator_.pointer_ = a;
    }

    static std::string name()
    {
        return "AccumulatorEnd (internal)";
    }
        
    bool operator()() const { return false; }
    bool get() const { return false; }
    
    template <unsigned, class U>
    void pass(U const &) 
    {}
    
    template <unsigned, class U>
    void pass(U const &, double) 
    {}
    
    template <class U>
    void mergeImpl(U const &) 
    {}
    
    template <class U>
    void resize(U const &) 
    {}
        
    template <class U>
    void setCoordinateOffsetImpl(U const &)
    {}
    
    void activate() 
    {}
    
    bool isActive() const 
    { 
        return false;
    }
    
    template <class Flags>
    static void activateImpl(Flags &)
    {}
    
    template <class Accu, class Flags1, class Flags2>
    static void activateImpl(Flags1 &, Flags2 &)
    {}
    
    template <class Flags>
    static bool isActiveImpl(Flags const &)
    {
        return true;
    }
    
    void applyHistogramOptions(HistogramOptions const &)
    {}
    
    static unsigned int passesRequired()
    {
        return 0;
    }
    
    static unsigned int passesRequired(AccumulatorFlags const &)
    {
        return 0;
    }

    void reset()
    {
        active_accumulators_.clear();
        is_dirty_.clear();
    }
        
    template <int which>
    void setDirtyImpl() const
    {
        is_dirty_.template set<which>();
    }
    
    template <int which>
    void setCleanImpl() const
    {
        is_dirty_.template reset<which>();
    }
    
    template <int which>
    bool isDirtyImpl() const
    {
        return is_dirty_.template test<which>();
    }
};

    // DecoratorImpl implement the functionality of Decorator below
template <class A, unsigned CurrentPass, bool allowRuntimeActivation, unsigned WorkPass=A::workInPass>
struct DecoratorImpl
{
    template <class T>
    static void exec(A & a, T const & t)
    {}

    template <class T>
    static void exec(A & a, T const & t, double weight)
    {}
};

template <class A, unsigned CurrentPass>
struct DecoratorImpl<A, CurrentPass, false, CurrentPass>
{
    template <class T>
    static void exec(A & a, T const & t)
    {
        a.update(t);
    }
    
    template <class T>
    static void exec(A & a, T const & t, double weight)
    {
        a.update(t, weight);
    }

    static typename A::result_type get(A const & a)
    {
        return a();
    }

    static void mergeImpl(A & a, A const & o)
    {
        a += o;
    }

    template <class T>
    static void resize(A & a, T const & t)
    {
        a.reshape(t);
    }
    
    static void applyHistogramOptions(A & a, HistogramOptions const & options)
    {
        ApplyHistogramOptions<typename A::Tag>::exec(a, options);
    }

    static unsigned int passesRequired()
    {
        static const unsigned int A_workInPass = A::workInPass;
        return std::max(A_workInPass, A::InternalBaseType::passesRequired());
    }
};

template <class A, unsigned CurrentPass>
struct DecoratorImpl<A, CurrentPass, true, CurrentPass>
{
    static bool isActive(A const & a)
    {
        return A::isActiveImpl(getAccumulator<AccumulatorEnd>(a).active_accumulators_);
    }
    
    template <class T>
    static void exec(A & a, T const & t)
    {
        if(isActive(a))
            a.update(t);
    }

    template <class T>
    static void exec(A & a, T const & t, double weight)
    {
        if(isActive(a))
            a.update(t, weight);
    }

    static typename A::result_type get(A const & a)
    {
        if(!isActive(a))
        {
            std::string message = std::string("get(accumulator): attempt to access inactive statistic '") +
                                              A::Tag::name() + "'.";
            vigra_precondition(false, message);
        }
        return a();
    }

    static void mergeImpl(A & a, A const & o)
    {
        if(isActive(a))
            a += o;
    }

    template <class T>
    static void resize(A & a, T const & t)
    {
        if(isActive(a))
            a.reshape(t);
    }
    
    static void applyHistogramOptions(A & a, HistogramOptions const & options)
    {
        if(isActive(a))
            ApplyHistogramOptions<typename A::Tag>::exec(a, options);
    }
    
    template <class ActiveFlags>
    static unsigned int passesRequired(ActiveFlags const & flags)
    {
        static const unsigned int A_workInPass = A::workInPass;
        return A::isActiveImpl(flags)
                   ? std::max(A_workInPass, A::InternalBaseType::passesRequired(flags))
                   : A::InternalBaseType::passesRequired(flags);
    }
};

    // Generic reshape function (expands to a no-op when T has fixed shape, and to
    // the appropriate specialized call otherwise). Shape is an instance of MultiArrayShape<N>::type.
template <class T, class Shape>
void reshapeImpl(T &, Shape const &)
{}

template <class T, class Shape, class Initial>
void reshapeImpl(T &, Shape const &, Initial const & = T())
{}

template <unsigned int N, class T, class Alloc, class Shape>
void reshapeImpl(MultiArray<N, T, Alloc> & a, Shape const & s, T const & initial = T())
{
    MultiArray<N, T, Alloc>(s, initial).swap(a);
}

template <class T, class Alloc, class Shape>
void reshapeImpl(Matrix<T, Alloc> & a, Shape const & s, T const & initial = T())
{
    Matrix<T, Alloc>(s, initial).swap(a);
}

template <class T, class U>
void copyShapeImpl(T const &, U const &)   // to be used for scalars and static arrays
{}

template <unsigned int N, class T, class Alloc, class U>
void copyShapeImpl(MultiArray<N, T, Alloc> const & from, U & to) 
{
    to.reshape(from.shape());
}

template <class T, class Alloc, class U>
void copyShapeImpl(Matrix<T, Alloc> const & from, U & to) 
{
    to.reshape(from.shape());
}

template <class T, class U>
bool hasDataImpl(T const &)   // to be used for scalars and static arrays
{
    return true;
}

template <unsigned int N, class T, class Stride>
bool hasDataImpl(MultiArrayView<N, T, Stride> const & a) 
{
    return a.hasData();
}

    // generic functions to create suitable shape objects from various input data types 
template <unsigned int N, class T, class Stride>
inline typename MultiArrayShape<N>::type
shapeOf(MultiArrayView<N, T, Stride> const & a)
{
    return a.shape();
}

template <class T, int N>
inline Shape1
shapeOf(TinyVector<T, N> const &)
{
    return Shape1(N);
}

template <class T, class NEXT>
inline CoupledHandle<T, NEXT> const &
shapeOf(CoupledHandle<T, NEXT> const & t)
{
    return t;
}

#define VIGRA_SHAPE_OF(type) \
inline Shape1 \
shapeOf(type) \
{ \
    return Shape1(1); \
}

VIGRA_SHAPE_OF(unsigned char)
VIGRA_SHAPE_OF(signed char)
VIGRA_SHAPE_OF(unsigned short)
VIGRA_SHAPE_OF(short)
VIGRA_SHAPE_OF(unsigned int)
VIGRA_SHAPE_OF(int)
VIGRA_SHAPE_OF(unsigned long)
VIGRA_SHAPE_OF(long)
VIGRA_SHAPE_OF(unsigned long long)
VIGRA_SHAPE_OF(long long)
VIGRA_SHAPE_OF(float)
VIGRA_SHAPE_OF(double)
VIGRA_SHAPE_OF(long double)

#undef VIGRA_SHAPE_OF

    // LabelDispatch is only used in AccumulatorChainArrays and has the following functionalities:
    //  * hold an accumulator chain for global statistics
    //  * hold an array of accumulator chains (one per region) for region statistics
    //  * forward data to the appropriate chains
    //  * allocate the region array with appropriate size
    //  * store and forward activation requests
    //  * compute required number of passes as maximum from global and region accumulators
template <class T, class GlobalAccumulators, class RegionAccumulators>
struct LabelDispatch
{
    typedef LabelDispatchTag Tag;
    typedef GlobalAccumulators GlobalAccumulatorChain;
    typedef RegionAccumulators RegionAccumulatorChain;
    typedef typename LookupTag<AccumulatorEnd, RegionAccumulatorChain>::type::AccumulatorFlags ActiveFlagsType;
    typedef ArrayVector<RegionAccumulatorChain> RegionAccumulatorArray;
        
    typedef LabelDispatch type;
    typedef LabelDispatch & reference;
    typedef LabelDispatch const & const_reference;
    typedef GlobalAccumulatorChain InternalBaseType;
    
    typedef T const & argument_type;
    typedef argument_type first_argument_type;
    typedef double second_argument_type;
    typedef RegionAccumulatorChain & result_type;
    
    static const int index = GlobalAccumulatorChain::index + 1;
    
    template <class IndexDefinition, class TagFound=typename IndexDefinition::Tag>
    struct CoordIndexSelector
    {
        static const int value = 0; // default: CoupledHandle holds coordinates at index 0 
    };
    
    template <class IndexDefinition>
    struct CoordIndexSelector<IndexDefinition, CoordArgTag>
    {
        static const int value = IndexDefinition::value;
    };
    
    static const int coordIndex = CoordIndexSelector<typename LookupTag<CoordArgTag, GlobalAccumulatorChain>::type>::value;
    static const int coordSize  = CoupledHandleCast<coordIndex, T>::type::value_type::static_size;
    typedef TinyVector<double, coordSize> CoordinateType;
    
    GlobalAccumulatorChain next_;
    RegionAccumulatorArray regions_;
    HistogramOptions region_histogram_options_;
    MultiArrayIndex ignore_label_;
    ActiveFlagsType active_region_accumulators_;
    CoordinateType coordinateOffset_;
    
    template <class TAG>
    struct ActivateImpl
    {
        typedef typename LookupTag<TAG, type>::type TargetAccumulator;
        
        static void activate(GlobalAccumulatorChain & globals, RegionAccumulatorArray & regions, 
                             ActiveFlagsType & flags)
        {
            TargetAccumulator::template activateImpl<LabelDispatch>(
                      flags, getAccumulator<AccumulatorEnd>(globals).active_accumulators_);
            for(unsigned int k=0; k<regions.size(); ++k)
                getAccumulator<AccumulatorEnd>(regions[k]).active_accumulators_ = flags;
        }
        
        static bool isActive(GlobalAccumulatorChain const &, ActiveFlagsType const & flags)
        {
            return TargetAccumulator::isActiveImpl(flags);
        }
    };
    
    template <class TAG>
    struct ActivateImpl<Global<TAG> >
    {
        static void activate(GlobalAccumulatorChain & globals, RegionAccumulatorArray &, ActiveFlagsType &)
        {
            LookupTag<TAG, GlobalAccumulatorChain>::type::activateImpl(getAccumulator<AccumulatorEnd>(globals).active_accumulators_);
        }
        
        static bool isActive(GlobalAccumulatorChain const & globals, ActiveFlagsType const &)
        {
            return LookupTag<TAG, GlobalAccumulatorChain>::type::isActiveImpl(getAccumulator<AccumulatorEnd>(globals).active_accumulators_);
        }
    };
    
    template <int INDEX>
    struct ActivateImpl<LabelArg<INDEX> >
    {
        static void activate(GlobalAccumulatorChain &, RegionAccumulatorArray &, ActiveFlagsType &)
        {}
        
        static bool isActive(GlobalAccumulatorChain const & globals, ActiveFlagsType const &)
        {
            return getAccumulator<LabelArg<INDEX> >(globals).isActive();
        }
    };
    
    LabelDispatch()
    : next_(),
      regions_(),
      region_histogram_options_(),
      ignore_label_(-1),
      active_region_accumulators_()
    {}
    
    LabelDispatch(LabelDispatch const & o)
    : next_(o.next_),
      regions_(o.regions_),
      region_histogram_options_(o.region_histogram_options_),
      ignore_label_(o.ignore_label_),
      active_region_accumulators_(o.active_region_accumulators_)
    {
        for(unsigned int k=0; k<regions_.size(); ++k)
        {
            getAccumulator<AccumulatorEnd>(regions_[k]).setGlobalAccumulator(&next_);
        }
    }
    
    MultiArrayIndex maxRegionLabel() const
    {
        return (MultiArrayIndex)regions_.size() - 1;
    }
    
    void setMaxRegionLabel(unsigned maxlabel)
    {
        if(maxRegionLabel() == (MultiArrayIndex)maxlabel)
            return;
        unsigned int oldSize = regions_.size();
        regions_.resize(maxlabel + 1);
        for(unsigned int k=oldSize; k<regions_.size(); ++k)
        {
            getAccumulator<AccumulatorEnd>(regions_[k]).setGlobalAccumulator(&next_);
            getAccumulator<AccumulatorEnd>(regions_[k]).active_accumulators_ = active_region_accumulators_;
            regions_[k].applyHistogramOptions(region_histogram_options_);
            regions_[k].setCoordinateOffsetImpl(coordinateOffset_);
        }
    }
    
    void ignoreLabel(MultiArrayIndex l)
    {
        ignore_label_ = l;
    }
    
    MultiArrayIndex ignoredLabel() const
    {
        return ignore_label_;
    }
    
    void applyHistogramOptions(HistogramOptions const & options)
    {
        applyHistogramOptions(options, options);
    }
    
    void applyHistogramOptions(HistogramOptions const & regionoptions, 
                               HistogramOptions const & globaloptions)
    {
        region_histogram_options_ = regionoptions;
        for(unsigned int k=0; k<regions_.size(); ++k)
        {
            regions_[k].applyHistogramOptions(region_histogram_options_);
        }
        next_.applyHistogramOptions(globaloptions);
    }
    
    void setCoordinateOffsetImpl(CoordinateType const & offset)
    {
        coordinateOffset_ = offset;
        for(unsigned int k=0; k<regions_.size(); ++k)
        {
            regions_[k].setCoordinateOffsetImpl(coordinateOffset_);
        }
        next_.setCoordinateOffsetImpl(coordinateOffset_);
    }
    
    void setCoordinateOffsetImpl(MultiArrayIndex k, CoordinateType const & offset)
    {
        vigra_precondition(0 <= k && k < (MultiArrayIndex)regions_.size(),
             "Accumulator::setCoordinateOffset(k, offset): region k does not exist.");
        regions_[k].setCoordinateOffsetImpl(offset);
    }
    
    template <class U>
    void resize(U const & t)
    {
        if(regions_.size() == 0)
        {
            typedef HandleArgSelector<U, LabelArgTag, GlobalAccumulatorChain> LabelHandle;
            typedef typename LabelHandle::value_type LabelType;
            typedef MultiArrayView<LabelHandle::size, LabelType, StridedArrayTag> LabelArray;
            LabelArray labelArray(t.shape(), LabelHandle::getHandle(t).strides(), 
                                  const_cast<LabelType *>(LabelHandle::getHandle(t).ptr()));
            
            LabelType minimum, maximum;
            labelArray.minmax(&minimum, &maximum);
            setMaxRegionLabel(maximum);
        }
        next_.resize(t);
        // FIXME: only call resize when label k actually exists?
        for(unsigned int k=0; k<regions_.size(); ++k)
            regions_[k].resize(t);
    }
    
    template <unsigned N>
    void pass(T const & t)
    {
        typedef HandleArgSelector<T, LabelArgTag, GlobalAccumulatorChain> LabelHandle;
        if(LabelHandle::getValue(t) != ignore_label_)
        {
            next_.template pass<N>(t);
            regions_[LabelHandle::getValue(t)].template pass<N>(t);
        }
    }
    
    template <unsigned N>
    void pass(T const & t, double weight)
    {
        typedef HandleArgSelector<T, LabelArgTag, GlobalAccumulatorChain> LabelHandle;
        if(LabelHandle::getValue(t) != ignore_label_)
        {
            next_.template pass<N>(t, weight);
            regions_[LabelHandle::getValue(t)].template pass<N>(t, weight);
        }
    }
    
    static unsigned int passesRequired()
    {
        return std::max(GlobalAccumulatorChain::passesRequired(), RegionAccumulatorChain::passesRequired());
    }
    
    unsigned int passesRequiredDynamic() const
    {
        return std::max(GlobalAccumulatorChain::passesRequired(getAccumulator<AccumulatorEnd>(next_).active_accumulators_), 
                        RegionAccumulatorChain::passesRequired(active_region_accumulators_));
    }
    
    void reset()
    {
        next_.reset();
        
        active_region_accumulators_.clear();
        RegionAccumulatorArray().swap(regions_);
        // FIXME: or is it better to just reset the region accumulators?
        // for(unsigned int k=0; k<regions_.size(); ++k)
            // regions_[k].reset();
    }
    
    template <class TAG>
    void activate()
    {
        ActivateImpl<TAG>::activate(next_, regions_, active_region_accumulators_);
    }
    
    void activateAll()
    {
        getAccumulator<AccumulatorEnd>(next_).active_accumulators_.set();
        active_region_accumulators_.set();
        for(unsigned int k=0; k<regions_.size(); ++k)
            getAccumulator<AccumulatorEnd>(regions_[k]).active_accumulators_.set();
    }
    
    template <class TAG>
    bool isActive() const
    {
        return ActivateImpl<TAG>::isActive(next_, active_region_accumulators_);
    }
    
    void mergeImpl(LabelDispatch const & o)
    {
        for(unsigned int k=0; k<regions_.size(); ++k)
            regions_[k].mergeImpl(o.regions_[k]);
        next_.mergeImpl(o.next_);
    }
    
    void mergeImpl(unsigned i, unsigned j)
    {
        regions_[i].mergeImpl(regions_[j]);
        regions_[j].reset();
        getAccumulator<AccumulatorEnd>(regions_[j]).active_accumulators_ = active_region_accumulators_;
    }
    
    template <class ArrayLike>
    void mergeImpl(LabelDispatch const & o, ArrayLike const & labelMapping)
    {
        MultiArrayIndex newMaxLabel = std::max<MultiArrayIndex>(maxRegionLabel(), *argMax(labelMapping.begin(), labelMapping.end()));
        setMaxRegionLabel(newMaxLabel);
        for(unsigned int k=0; k<labelMapping.size(); ++k)
            regions_[labelMapping[k]].mergeImpl(o.regions_[k]);
        next_.mergeImpl(o.next_);
    }
};

template <class TargetTag, class TagList>
struct FindNextTag;

template <class TargetTag, class HEAD, class TAIL>
struct FindNextTag<TargetTag, TypeList<HEAD, TAIL> >
{
    typedef typename FindNextTag<TargetTag, TAIL>::type type;
};

template <class TargetTag, class TAIL>
struct FindNextTag<TargetTag, TypeList<TargetTag, TAIL> >
{
    typedef typename TAIL::Head type;
};

template <class TargetTag>
struct FindNextTag<TargetTag, TypeList<TargetTag, void> >
{
    typedef void type;
};

template <class TargetTag>
struct FindNextTag<TargetTag, void>
{
    typedef void type;
};

    // AccumulatorFactory creates the decorator hierarchy for the given TAG and configuration CONFIG
template <class TAG, class CONFIG, unsigned LEVEL=0>
struct AccumulatorFactory
{
    typedef typename FindNextTag<TAG, typename CONFIG::TagList>::type NextTag;
    typedef typename AccumulatorFactory<NextTag, CONFIG, LEVEL+1>::type NextType;
    typedef typename CONFIG::InputType InputType;
    
    template <class T>
    struct ConfigureTag
    {
        typedef TAG type;
    };
    
        // When InputType is a CoupledHandle, some tags need to be wrapped into 
        // DataFromHandle<> and/or Weighted<> modifiers. The following code does
        // this when appropriate.
    template <class T, class NEXT>
    struct ConfigureTag<CoupledHandle<T, NEXT> >
    {
        typedef typename StandardizeTag<DataFromHandle<TAG> >::type WrappedTag;
        typedef typename IfBool<(!HasModifierPriority<WrappedTag, WeightingPriority>::value && ShouldBeWeighted<WrappedTag>::value),
                                 Weighted<WrappedTag>, WrappedTag>::type type;
    };
    
    typedef typename ConfigureTag<InputType>::type UseTag;
    
        // base class of the decorator hierarchy: default (possibly empty) 
        // implementations of all members
    struct AccumulatorBase
    {
        typedef AccumulatorBase              ThisType;
        typedef TAG                          Tag;
        typedef NextType                     InternalBaseType;
        typedef InputType                    input_type;
        typedef input_type const &           argument_type;
        typedef argument_type                first_argument_type;
        typedef double                       second_argument_type;
        typedef void                         result_type;
        
        static const unsigned int            workInPass = 1;
        static const int                     index = InternalBaseType::index + 1;
        
        InternalBaseType next_;
        
        static std::string name()
        {
            return TAG::name();
        }
        
        template <class ActiveFlags>
        static void activateImpl(ActiveFlags & flags)
        {
            flags.template set<index>();
            typedef typename StandardizeDependencies<Tag>::type StdDeps;
            acc_detail::ActivateDependencies<StdDeps>::template exec<ThisType>(flags);
        }
        
        template <class Accu, class ActiveFlags, class GlobalFlags>
        static void activateImpl(ActiveFlags & flags, GlobalFlags & gflags)
        {
            flags.template set<index>();
            typedef typename StandardizeDependencies<Tag>::type StdDeps;
            acc_detail::ActivateDependencies<StdDeps>::template exec<Accu>(flags, gflags);
        }
        
        template <class ActiveFlags>
        static bool isActiveImpl(ActiveFlags & flags)
        {
            return flags.template test<index>();
        }
        
        void setDirty() const
        {
            next_.template setDirtyImpl<index>();
        }
        
        template <int INDEX>
        void setDirtyImpl() const
        {
            next_.template setDirtyImpl<INDEX>();
        }
        
        void setClean() const
        {
            next_.template setCleanImpl<index>();
        }
        
        template <int INDEX>
        void setCleanImpl() const
        {
            next_.template setCleanImpl<INDEX>();
        }
        
        bool isDirty() const
        {
            return next_.template isDirtyImpl<index>();
        }
        
        template <int INDEX>
        bool isDirtyImpl() const
        {
            return next_.template isDirtyImpl<INDEX>();
        }
        
        void reset()
        {}
        
        template <class Shape>
        void setCoordinateOffset(Shape const &)
        {}
        
        template <class Shape>
        void reshape(Shape const &)
        {}
        
        void operator+=(AccumulatorBase const &)
        {}
        
        template <class U>
        void update(U const &)
        {}
        
        template <class U>
        void update(U const &, double)
        {}
        
        template <class TargetTag>
        typename LookupDependency<TargetTag, ThisType>::result_type
        call_getDependency() const
        {
            return getDependency<TargetTag>(*this);
        }
    };

        // The middle class(es) of the decorator hierarchy implement the actual feature computation.
    typedef typename UseTag::template Impl<InputType, AccumulatorBase> AccumulatorImpl;
    
        // outer class of the decorator hierarchy. It has the following functionalities
        //  * ensure that only active accumulators are called in a dynamic accumulator chain
        //  * ensure that each accumulator is only called in its desired pass as defined in A::workInPass
        //  * determine how many passes through the data are required
    struct Accumulator
    : public AccumulatorImpl
    {
        typedef Accumulator type;
        typedef Accumulator & reference;
        typedef Accumulator const & const_reference;
        typedef AccumulatorImpl A;
        
        static const unsigned int workInPass = A::workInPass;
        static const bool allowRuntimeActivation = CONFIG::allowRuntimeActivation;
        
        template <class T>
        void resize(T const & t)
        {
            this->next_.resize(t);
            DecoratorImpl<Accumulator, workInPass, allowRuntimeActivation>::resize(*this, t);
        }
        
        void reset()
        {
            this->next_.reset();
            A::reset();
        }
        
        typename A::result_type get() const
        {
            return DecoratorImpl<A, workInPass, allowRuntimeActivation>::get(*this);
        }
        
        template <unsigned N, class T>
        void pass(T const & t)
        {
            this->next_.template pass<N>(t);
            DecoratorImpl<Accumulator, N, allowRuntimeActivation>::exec(*this, t);
        }
        
        template <unsigned N, class T>
        void pass(T const & t, double weight)
        {
            this->next_.template pass<N>(t, weight);
            DecoratorImpl<Accumulator, N, allowRuntimeActivation>::exec(*this, t, weight);
        }
        
        void mergeImpl(Accumulator const & o)
        {
            DecoratorImpl<Accumulator, Accumulator::workInPass, allowRuntimeActivation>::mergeImpl(*this, o);
            this->next_.mergeImpl(o.next_);
        }
        
        void applyHistogramOptions(HistogramOptions const & options)
        {
            DecoratorImpl<Accumulator, workInPass, allowRuntimeActivation>::applyHistogramOptions(*this, options);
            this->next_.applyHistogramOptions(options);
        }
        
        template <class SHAPE>
        void setCoordinateOffsetImpl(SHAPE const & offset)
        {
            this->setCoordinateOffset(offset);
            this->next_.setCoordinateOffsetImpl(offset);
        }
        
        static unsigned int passesRequired()
        {
            return DecoratorImpl<Accumulator, workInPass, allowRuntimeActivation>::passesRequired();
        }
        
        template <class ActiveFlags>
        static unsigned int passesRequired(ActiveFlags const & flags)
        {
            return DecoratorImpl<Accumulator, workInPass, allowRuntimeActivation>::passesRequired(flags);
        }
    };

    typedef Accumulator type;
};

template <class CONFIG, unsigned LEVEL>
struct AccumulatorFactory<void, CONFIG, LEVEL>
{
    typedef AccumulatorEndImpl<LEVEL, typename CONFIG::GlobalAccumulatorHandle> type;
};

struct InvalidGlobalAccumulatorHandle
{
    typedef Error__Global_statistics_are_only_defined_for_AccumulatorChainArray type;
    
    InvalidGlobalAccumulatorHandle()
    : pointer_(0)
    {}
    
    type const * pointer_;
};

    // helper classes to create an accumulator chain from a TypeList
    // if dynamic=true,  a dynamic accumulator will be created
    // if dynamic=false, a plain accumulator will be created
template <class T, class Selected, bool dynamic=false, class GlobalHandle=InvalidGlobalAccumulatorHandle>
struct ConfigureAccumulatorChain
#ifndef DOXYGEN
: public ConfigureAccumulatorChain<T, typename AddDependencies<typename Selected::type>::type, dynamic>
#endif
{};

template <class T, class HEAD, class TAIL, bool dynamic, class GlobalHandle>
struct ConfigureAccumulatorChain<T, TypeList<HEAD, TAIL>, dynamic, GlobalHandle>
{
    typedef TypeList<HEAD, TAIL> TagList;
    typedef T InputType;
    static const bool allowRuntimeActivation = dynamic;
    typedef GlobalHandle GlobalAccumulatorHandle;
 
    typedef typename AccumulatorFactory<HEAD, ConfigureAccumulatorChain>::type type;
};

template <class T, class Selected, bool dynamic=false>
struct ConfigureAccumulatorChainArray
#ifndef DOXYGEN
: public ConfigureAccumulatorChainArray<T, typename AddDependencies<typename Selected::type>::type, dynamic>
#endif
{};

template <class T, class HEAD, class TAIL, bool dynamic>
struct ConfigureAccumulatorChainArray<T, TypeList<HEAD, TAIL>, dynamic>
{
    typedef TypeList<HEAD, TAIL> TagList;
    typedef SeparateGlobalAndRegionTags<TagList> TagSeparator;
    typedef typename TagSeparator::GlobalTags GlobalTags;
    typedef typename TagSeparator::RegionTags RegionTags;
    typedef typename ConfigureAccumulatorChain<T, GlobalTags, dynamic>::type GlobalAccumulatorChain;

    struct GlobalAccumulatorHandle
    {
        typedef GlobalAccumulatorChain type;
        
        GlobalAccumulatorHandle()
        : pointer_(0)
        {}
        
        type const * pointer_;
    };
    
    typedef typename ConfigureAccumulatorChain<T, RegionTags, dynamic, GlobalAccumulatorHandle>::type RegionAccumulatorChain;
    
    typedef LabelDispatch<T, GlobalAccumulatorChain, RegionAccumulatorChain> type;
};

} // namespace acc_detail 

/****************************************************************************/
/*                                                                          */
/*                            accumulator chain                             */
/*                                                                          */
/****************************************************************************/

// Implement the high-level interface of an accumulator chain
template <class T, class NEXT>
class AccumulatorChainImpl
{
  public:
    typedef NEXT                                             InternalBaseType;
    typedef AccumulatorBegin                                 Tag;
    typedef typename InternalBaseType::argument_type         argument_type;
    typedef typename InternalBaseType::first_argument_type   first_argument_type;
    typedef typename InternalBaseType::second_argument_type  second_argument_type;
    typedef void                                             value_type;
    typedef typename InternalBaseType::result_type           result_type;
    
    static const int staticSize = InternalBaseType::index;

    InternalBaseType next_;

    /** \brief Current pass of the accumulator chain.
    */
    unsigned int current_pass_;
    
    AccumulatorChainImpl()
    : current_pass_(0)
    {}

    /** Set options for all histograms in the accumulator chain. See histogram accumulators for possible options. The function is ignored if there is no histogram in the accumulator chain.
    */
    void setHistogramOptions(HistogramOptions const & options)
    {
        next_.applyHistogramOptions(options);
    }
    

    /** Set regional and global options for all histograms in the accumulator chain.
    */
    void setHistogramOptions(HistogramOptions const & regionoptions, HistogramOptions const & globaloptions)
    {
        next_.applyHistogramOptions(regionoptions, globaloptions);
    }
    
    /** Set an offset for <tt>Coord<...></tt> statistics.
    
        If the offset is non-zero, coordinate statistics such as <tt>RegionCenter</tt> are computed
        in the global coordinate system defined by the \a offset. Without an offset, these statistics
        are computed in the local coordinate system of the current region of interest.
    */    
    template <class SHAPE>
    void setCoordinateOffset(SHAPE const & offset)
    {
        next_.setCoordinateOffsetImpl(offset);
    }
    
    /** Reset current_pass_ of the accumulator chain to 'reset_to_pass'.
    */
    void reset(unsigned int reset_to_pass = 0)
    {
        current_pass_ = reset_to_pass;
        if(reset_to_pass == 0)
            next_.reset();
    }
    
    template <unsigned N>
    void update(T const & t)
    {
        if(current_pass_ == N)
        {
            next_.template pass<N>(t);
        }
        else if(current_pass_ < N)
        {
            current_pass_ = N;
            if(N == 1)
                next_.resize(acc_detail::shapeOf(t));
            next_.template pass<N>(t);
        }
        else
        {
            std::string message("AccumulatorChain::update(): cannot return to pass ");
            message << N << " after working on pass " << current_pass_ << ".";
            vigra_precondition(false, message);
        }
    }
    
    template <unsigned N>
    void update(T const & t, double weight)
    {
        if(current_pass_ == N)
        {
            next_.template pass<N>(t, weight);
        }
        else if(current_pass_ < N)
        {
            current_pass_ = N;
            if(N == 1)
                next_.resize(acc_detail::shapeOf(t));
            next_.template pass<N>(t, weight);
        }
        else
        {
            std::string message("AccumulatorChain::update(): cannot return to pass ");
            message << N << " after working on pass " << current_pass_ << ".";
            vigra_precondition(false, message);
       }
    }
    
    /** Equivalent to merge(o) .
    */
    void operator+=(AccumulatorChainImpl const & o)
    {
        merge(o);
    }
    
    /** Merge the accumulator chain with accumulator chain 'o'. This only works if all selected statistics in the accumulator chain support the '+=' operator. See the documentations of the particular statistics for support information.
    */
    void merge(AccumulatorChainImpl const & o)
    {
        next_.mergeImpl(o.next_);
    }

    result_type operator()() const
    {
        return next_.get();
    }

    void operator()(T const & t)
    {
        update<1>(t);
    }
    
    void operator()(T const & t, double weight)
    {
        update<1>(t, weight);
    }

    void updatePass2(T const & t)
    {
        update<2>(t);
    }
    
    void updatePass2(T const & t, double weight)
    {
        update<2>(t, weight);
    }

    /** Upate all accumulators in the accumulator chain that work in pass N with data t. Requirement: 0 < N < 6 and N >= current_pass_ . If N < current_pass_ call reset() first.  
    */
    void updatePassN(T const & t, unsigned int N)
    {
        switch (N)
        {
            case 1: update<1>(t); break;
            case 2: update<2>(t); break;
            case 3: update<3>(t); break;
            case 4: update<4>(t); break;
            case 5: update<5>(t); break;
            default:
                vigra_precondition(false,
                     "AccumulatorChain::updatePassN(): 0 < N < 6 required.");
        }
    }
    
    /** Upate all accumulators in the accumulator chain that work in pass N with data t and weight. Requirement: 0 < N < 6 and N >= current_pass_ . If N < current_pass_ call reset() first. 
    */
    void updatePassN(T const & t, double weight, unsigned int N)
    {
        switch (N)
        {
            case 1: update<1>(t, weight); break;
            case 2: update<2>(t, weight); break;
            case 3: update<3>(t, weight); break;
            case 4: update<4>(t, weight); break;
            case 5: update<5>(t, weight); break;
            default:
                vigra_precondition(false,
                     "AccumulatorChain::updatePassN(): 0 < N < 6 required.");
        }
    }
  
    /** Return the number of passes required to compute all statistics in the accumulator chain.
    */
    unsigned int passesRequired() const
    {
        return InternalBaseType::passesRequired();
    }
};



   // Create an accumulator chain containing the Selected statistics and their dependencies.

/** \brief Create an accumulator chain containing the selected statistics and their dependencies.

    AccumulatorChain is used to compute global statistics which have to be selected at compile time. 

    The template parameters are as follows:
    - T: The input type
        - either element type of the data(e.g. double, int, RGBValue, ...)
        - or type of CoupledHandle (for simultaneous access to coordinates and/or weights)
    - Selected: statistics to be computed and index specifier for the CoupledHandle, wrapped with Select
    
    Usage:
    \code
    typedef double DataType;
    AccumulatorChain<DataType, Select<Variance, Mean, Minimum, ...> > accumulator;
    \endcode

    Usage, using CoupledHandle:
    \code
    const int dim = 3; //dimension of MultiArray
    typedef double DataType;
    typedef double WeightType;
    typedef vigra::CoupledIteratorType<dim, DataType, WeightType>::HandleType Handle;
    AccumulatorChain<Handle, Select<DataArg<1>, WeightArg<2>, Mean,...> > a;
    \endcode

    See \ref FeatureAccumulators for more information and examples of use.
 */
template <class T, class Selected, bool dynamic=false>
class AccumulatorChain
#ifndef DOXYGEN // hide AccumulatorChainImpl from documentation
: public AccumulatorChainImpl<T, typename acc_detail::ConfigureAccumulatorChain<T, Selected, dynamic>::type>
#endif
{
  public:
  // \brief TypeList of Tags in the accumulator chain (?).
    typedef typename acc_detail::ConfigureAccumulatorChain<T, Selected, dynamic>::TagList AccumulatorTags;
  
    /** Before having seen data (current_pass_==0), the shape of the data can be changed... (?)
    */
    template <class U, int N>
    void reshape(TinyVector<U, N> const & s)
    {
        vigra_precondition(this->current_pass_ == 0,
             "AccumulatorChain::reshape(): cannot reshape after seeing data. Call AccumulatorChain::reset() first.");
        this->next_.resize(s);
        this->current_pass_ = 1;
    }
     
    /** Return the names of all tags in the accumulator chain (selected statistics and their dependencies).
    */
    static ArrayVector<std::string> const & tagNames()
    {
        static ArrayVector<std::string> * n = VIGRA_SAFE_STATIC(n, new ArrayVector<std::string>(collectTagNames()));
        return *n;
    }


#ifdef DOXYGEN // hide AccumulatorChainImpl from documentation
  
  /** Set options for all histograms in the accumulator chain. See histogram accumulators for possible options. The function is ignored if there is no histogram in the accumulator chain.
   */
  void setHistogramOptions(HistogramOptions const & options);
    
  /** Set an offset for <tt>Coord<...></tt> statistics.
  
      If the offset is non-zero, coordinate statistics such as <tt>RegionCenter</tt> are computed
      in the global coordinate system defined by the \a offset. Without an offset, these statistics
      are computed in the local coordinate system of the current region of interest.
  */    
  template <class SHAPE>
  void setCoordinateOffset(SHAPE const & offset);
    
  /** Reset current_pass_ of the accumulator chain to 'reset_to_pass'. */
  void reset(unsigned int reset_to_pass = 0);

  /** Equivalent to merge(o) . */
  void operator+=(AccumulatorChainImpl const & o);
  
  /** Merge the accumulator chain with accumulator chain 'o'. This only works if all selected statistics in the accumulator chain support the '+=' operator. See the documentations of the particular statistics for support information.
   */
  void merge(AccumulatorChainImpl const & o);
  
  /** Upate all accumulators in the accumulator chain that work in pass N with data t. Requirement: 0 < N < 6 and N >= current_pass_ . If N < current_pass_ call reset first.  
   */
  void updatePassN(T const & t, unsigned int N);
  
  /** Upate all accumulators in the accumulator chain that work in pass N with data t and weight. Requirement: 0 < N < 6 and N >= current_pass_ . If N < current_pass_ call reset first. 
   */
  void updatePassN(T const & t, double weight, unsigned int N);
  
  /** Return the number of passes required to compute all statistics in the accumulator chain.
   */
  unsigned int passesRequired() const;
  
#endif   
 
  private:
    static ArrayVector<std::string> collectTagNames()
    {
        ArrayVector<std::string> n;
        acc_detail::CollectAccumulatorNames<AccumulatorTags>::exec(n);
        std::sort(n.begin(), n.end());
        return n;
    }
}; 

template <unsigned int N, class T1, class T2, class T3, class T4, class T5, class Selected, bool dynamic>
class AccumulatorChain<CoupledArrays<N, T1, T2, T3, T4, T5>, Selected, dynamic>
: public AccumulatorChain<typename CoupledArrays<N, T1, T2, T3, T4, T5>::HandleType, Selected, dynamic>
{};


    // Create a dynamic accumulator chain containing the Selected statistics and their dependencies.
    // Statistics will only be computed if activate<Tag>() is called at runtime.
/** \brief Create a dynamic accumulator chain containing the selected statistics and their dependencies.

    DynamicAccumulatorChain is used to compute global statistics with run-time activation. A set of statistics is selected at run-time and from this set statistics can be activated at run-time by calling activate<stat>() or activate(std::string stat).

    The template parameters are as follows:
    - T: The input type
        - either element type of the data(e.g. double, int, RGBValue, ...)
        - or type of CoupledHandle (for access to coordinates and/or weights)
    - Selected: statistics to be computed and index specifier for the CoupledHandle, wrapped with Select
    
    Usage:
    \code
    typedef double DataType;
    DynamicAccumulatorChain<DataType, Select<Variance, Mean, Minimum, ...> > accumulator;
    \endcode

    Usage, using CoupledHandle:
    \code
    const int dim = 3; //dimension of MultiArray
    typedef double DataType;
    typedef double WeightType;
    typedef vigra::CoupledIteratorType<dim, DataType, WeightType>::HandleType Handle;
    DynamicAccumulatorChain<Handle, Select<DataArg<1>, WeightArg<2>, Mean,...> > a;
    \endcode

    See \ref FeatureAccumulators for more information and examples of use.
 */
template <class T, class Selected>
class DynamicAccumulatorChain
: public AccumulatorChain<T, Selected, true>
{
  public:
    typedef typename AccumulatorChain<T, Selected, true>::InternalBaseType InternalBaseType;
    typedef typename DynamicAccumulatorChain::AccumulatorTags AccumulatorTags;
       
    /** Activate statistic 'tag'. Alias names are not recognized. If the statistic is not in the accumulator chain a PreconditionViolation is thrown.
    */
    void activate(std::string tag)
    {
        vigra_precondition(activateImpl(tag),
            std::string("DynamicAccumulatorChain::activate(): Tag '") + tag + "' not found.");
    }
    
    /** %activate\<TAG\>() activates statistic 'TAG'. If the statistic is not in the accumulator chain it is ignored. (?)
    */
    template <class TAG>
    void activate()
    {
        LookupTag<TAG, DynamicAccumulatorChain>::type::activateImpl(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }
    
    /** Activate all statistics in the accumulator chain.
    */
    void activateAll()
    {
        getAccumulator<AccumulatorEnd>(*this).active_accumulators_.set();
    }
    /** Return true if the statistic 'tag' is active, i.e. activate(std::string tag) or activate<TAG>() has been called. If the statistic is not in the accumulator chain a PreconditionViolation is thrown. (Note that alias names are not recognized.)
    */
    bool isActive(std::string tag) const
    {
        acc_detail::TagIsActive_Visitor v;
        vigra_precondition(isActiveImpl(tag, v),
            std::string("DynamicAccumulatorChain::isActive(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    /** %isActive\<TAG\>() returns true if statistic 'TAG' is active, i.e. activate(std::string tag) or activate<TAG>() has been called. If the statistic is not in the accumulator chain, true is returned. (?)
    */
    template <class TAG>
    bool isActive() const
    {
        return LookupTag<TAG, DynamicAccumulatorChain>::type::isActiveImpl(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }

    /** Return names of all statistics in the accumulator chain that are active.
    */
    ArrayVector<std::string> activeNames() const
    {
        ArrayVector<std::string> res;
        for(unsigned k=0; k<DynamicAccumulatorChain::tagNames().size(); ++k)
            if(isActive(DynamicAccumulatorChain::tagNames()[k]))
                res.push_back(DynamicAccumulatorChain::tagNames()[k]);
        return res;
    }
    
    /** Return number of passes required to compute the active statistics in the accumulator chain.
    */
    unsigned int passesRequired() const
    {
        return InternalBaseType::passesRequired(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }
    
  protected:
  
    bool activateImpl(std::string tag)
    {
        return acc_detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, 
                                         normalizeString(tag), acc_detail::ActivateTag_Visitor());
    }
    
    bool isActiveImpl(std::string tag, acc_detail::TagIsActive_Visitor & v) const
    {
        return acc_detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, normalizeString(tag), v);
    }
};

template <unsigned int N, class T1, class T2, class T3, class T4, class T5, class Selected>
class DynamicAccumulatorChain<CoupledArrays<N, T1, T2, T3, T4, T5>, Selected>
: public DynamicAccumulatorChain<typename CoupledArrays<N, T1, T2, T3, T4, T5>::HandleType, Selected>
{};




template<unsigned int N, class T, class SELECT>
class StandAloneAccumulatorChain : public
AccumulatorChain<
    typename CoupledHandleType<N, T>::type,
    SELECT
>
{
public:
    typedef typename CoupledHandleType<N, T>::type HandleType;
    typedef typename HandleType::base_type CoordHandle;
    typedef typename MultiArrayShape<N>::type CoordType;

    typedef SELECT SelectType;
    typedef AccumulatorChain<HandleType, SelectType>  BaseType;


    StandAloneAccumulatorChain()
    :   BaseType(),
        val_(0.0),
        handle_(&val_, CoordType(), CoordHandle(CoordType()))
    {
    }

    
    

    void updatePassN(const T & val, const CoordType & coord, unsigned int p){
        //val_ = val;
        handle_.updatePtrAdresse(val_);
        handle_. template get<1>() = val;
        const CoordType & constP = handle_. template get<0>();
        CoordType & nonConstP = * const_cast<CoordType*>(&constP);
        
        std::copy(coord.begin(), coord.end(), nonConstP.begin());
        BaseType::updatePassN(handle_, p);
    }
private:
    T val_;
    HandleType handle_;
    //AcculmatorChainType accChain_;

};

template<unsigned int N, class SELECT>
class StandAloneDataFreeAccumulatorChain : public
AccumulatorChain<
    typename CoupledHandleType<N>::type,
    SELECT
>
{
public:
    typedef typename CoupledHandleType<N>::type HandleType;
    typedef typename MultiArrayShape<N>::type CoordType;

    typedef SELECT SelectType;
    typedef AccumulatorChain<HandleType, SelectType>  BaseType;


    StandAloneDataFreeAccumulatorChain()
    :   BaseType(),
        handle_(CoordType())
    {
        
    }

    template<class IGNORED_DATA>
    void updatePassN(
        const IGNORED_DATA & ignoreData,
        const CoordType & coord, unsigned int p){
        this->updatePassN(coord,p);
    }
    

    void updatePassN(const CoordType & coord, unsigned int p){
        const CoordType & constP = handle_. template get<0>();
        CoordType & nonConstP = * const_cast<CoordType*>(&constP);
        std::copy(coord.begin(), coord.end(), nonConstP.begin());
        BaseType::updatePassN(handle_, p);
    }
private:
    HandleType handle_;
    //AcculmatorChainType accChain_;

};





/** \brief Create an array of accumulator chains containing the selected per-region and global statistics and their dependencies.

    AccumulatorChainArray is used to compute per-region statistics (as well as global statistics). The statistics are selected at compile-time. An array of accumulator chains (one per region) for region statistics is created and one accumulator chain for global statistics. The region labels always start at 0. Use the Global modifier to compute global statistics (by default per-region statistics are computed). 

    The template parameters are as follows:
    - T: The input type, type of CoupledHandle (for access to coordinates, labels and weights)
    - Selected: statistics to be computed and index specifier for the CoupledHandle, wrapped with Select

    Usage:
    \code
    const int dim = 3; //dimension of MultiArray
    typedef double DataType;
    typedef double WeightType;
    typedef unsigned int LabelType;
    typedef vigra::CoupledIteratorType<dim, DataType, WeightType, LabelType>::HandleType Handle;
    AccumulatorChainArray<Handle, Select<DataArg<1>, WeightArg<2>, LabelArg<3>, Mean, Variance, ...> > a;
    \endcode

    See \ref FeatureAccumulators for more information and examples of use.
*/
template <class T, class Selected, bool dynamic=false>
class AccumulatorChainArray
#ifndef DOXYGEN //hide AccumulatorChainImpl vom documentation
: public AccumulatorChainImpl<T, typename acc_detail::ConfigureAccumulatorChainArray<T, Selected, dynamic>::type>
#endif
{
  public:
    typedef AccumulatorChainImpl<T, typename acc_detail::ConfigureAccumulatorChainArray<T, Selected, dynamic>::type> base_type;
    typedef typename acc_detail::ConfigureAccumulatorChainArray<T, Selected, dynamic> Creator;
    typedef typename Creator::TagList AccumulatorTags;
    typedef typename Creator::GlobalTags GlobalTags;
    typedef typename Creator::RegionTags RegionTags;
    
    /** Statistics will not be computed for label l. Note that only one label can be ignored.
    */
    void ignoreLabel(MultiArrayIndex l)
    {
        this->next_.ignoreLabel(l);
    }
    
    /** Ask for a label to be ignored. Default: -1 (meaning that no label is ignored).
    */
    MultiArrayIndex ignoredLabel() const
    {
        return this->next_.ignoredLabel();
    }
    
    /** Set the maximum region label (e.g. for merging two accumulator chains).
    */
    void setMaxRegionLabel(unsigned label)
    {
        this->next_.setMaxRegionLabel(label);
    }
    
    /** Maximum region label. (equal to regionCount() - 1)
    */
    MultiArrayIndex maxRegionLabel() const
    {
        return this->next_.maxRegionLabel();
    }
    
    /** Number of Regions. (equal to maxRegionLabel() + 1)
    */
    unsigned int regionCount() const
    {
        return this->next_.regions_.size();
    }
    
    /** Equivalent to <tt>merge(o)</tt>.
    */
    void operator+=(AccumulatorChainArray const & o)
    {
        merge(o);
    }
    
    /** Merge region i with region j. 
    */
    void merge(unsigned i, unsigned j)
    {
        vigra_precondition(i <= maxRegionLabel() && j <= maxRegionLabel(),
            "AccumulatorChainArray::merge(): region labels out of range.");
        this->next_.mergeImpl(i, j);
    }
    
    /** Merge with accumulator chain o. maxRegionLabel() of the two accumulators must be equal.
    */
    void merge(AccumulatorChainArray const & o)
    {
        if(maxRegionLabel() == -1)
            setMaxRegionLabel(o.maxRegionLabel());
        vigra_precondition(maxRegionLabel() == o.maxRegionLabel(),
            "AccumulatorChainArray::merge(): maxRegionLabel must be equal.");
        this->next_.mergeImpl(o.next_);
    }

    /** Merge with accumulator chain o using a mapping between labels of the two accumulators. Label l of accumulator chain o is mapped to labelMapping[l]. Hence, all elements of labelMapping must be <= maxRegionLabel() and size of labelMapping must match o.regionCount().
    */
    template <class ArrayLike>
    void merge(AccumulatorChainArray const & o, ArrayLike const & labelMapping)
    {
        vigra_precondition(labelMapping.size() == o.regionCount(),
            "AccumulatorChainArray::merge(): labelMapping.size() must match regionCount() of RHS.");
        this->next_.mergeImpl(o.next_, labelMapping);
    }

    /** Return names of all tags in the accumulator chain (selected statistics and their dependencies).
    */
    static ArrayVector<std::string> const & tagNames()
    {
        static const ArrayVector<std::string> n = collectTagNames();
        return n;
    }

    using base_type::setCoordinateOffset;
    
    /** Set an offset for <tt>Coord<...></tt> statistics for region \a k.
    
        If the offset is non-zero, coordinate statistics such as <tt>RegionCenter</tt> are computed
        in the global coordinate system defined by the \a offset. Without an offset, these statistics
        are computed in the local coordinate system of the current region of interest.
    */    
    template <class SHAPE>
    void setCoordinateOffset(MultiArrayIndex k, SHAPE const & offset)
    {
        this->next_.setCoordinateOffsetImpl(k, offset);
    }

#ifdef DOXYGEN // hide AccumulatorChainImpl from documentation

  /** \copydoc vigra::acc::AccumulatorChain::setHistogramOptions(HistogramOptions const &) */
  void setHistogramOptions(HistogramOptions const & options);

  /** Set regional and global options for all histograms in the accumulator chain.
   */
  void setHistogramOptions(HistogramOptions const & regionoptions, HistogramOptions const & globaloptions);
    
  /** \copydoc vigra::acc::AccumulatorChain::setCoordinateOffset(SHAPE const &)
  */    
  template <class SHAPE>
  void setCoordinateOffset(SHAPE const & offset)
  
  /** \copydoc vigra::acc::AccumulatorChain::reset() */
  void reset(unsigned int reset_to_pass = 0);

  /** \copydoc vigra::acc::AccumulatorChain::operator+=() */
  void operator+=(AccumulatorChainImpl const & o);
    
  /** \copydoc vigra::acc::AccumulatorChain::updatePassN(T const &,unsigned int) */
  void updatePassN(T const & t, unsigned int N);
  
  /** \copydoc vigra::acc::AccumulatorChain::updatePassN(T const &,double,unsigned int) */
  void updatePassN(T const & t, double weight, unsigned int N);
  
#endif
    
  private:
    static ArrayVector<std::string> collectTagNames()
    {
        ArrayVector<std::string> n;
        acc_detail::CollectAccumulatorNames<AccumulatorTags>::exec(n);
        std::sort(n.begin(), n.end());
        return n;
    }
};

template <unsigned int N, class T1, class T2, class T3, class T4, class T5, class Selected, bool dynamic>
class AccumulatorChainArray<CoupledArrays<N, T1, T2, T3, T4, T5>, Selected, dynamic>
: public AccumulatorChainArray<typename CoupledArrays<N, T1, T2, T3, T4, T5>::HandleType, Selected, dynamic>
{};

/** \brief Create an array of dynamic accumulator chains containing the selected per-region and global statistics and their dependencies.


    DynamicAccumulatorChainArray is used to compute per-region statistics (as well as global statistics) with run-time activation. A set of statistics is selected at run-time and from this set statistics can be activated at run-time by calling activate<stat>() or activate(std::string stat).

     The template parameters are as follows:
    - T: The input type, type of CoupledHandle (for access to coordinates, labels and weights)
    - Selected: statistics to be computed and index specifier for the CoupledHandle, wrapped with Select

    Usage:
    \code
    const int dim = 3; //dimension of MultiArray
    typedef double DataType;
    typedef double WeightType;
    typedef unsigned int LabelType;
    typedef vigra::CoupledIteratorType<dim, DataType, WeightType, LabelType>::HandleType Handle;
    DynamicAccumulatorChainArray<Handle, Select<DataArg<1>, WeightArg<2>, LabelArg<3>, Mean, Variance, ...> > a;
    \endcode

    See \ref FeatureAccumulators for more information and examples of use.
*/
template <class T, class Selected>
class DynamicAccumulatorChainArray
: public AccumulatorChainArray<T, Selected, true>
{
  public:
    typedef typename DynamicAccumulatorChainArray::AccumulatorTags AccumulatorTags;

    /** \copydoc DynamicAccumulatorChain::activate(std::string tag) */
    void activate(std::string tag)
    {
        vigra_precondition(activateImpl(tag),
            std::string("DynamicAccumulatorChainArray::activate(): Tag '") + tag + "' not found.");
    }
    
    /** \copydoc DynamicAccumulatorChain::activate() */
    template <class TAG>
    void activate()
    {
        this->next_.template activate<TAG>();
    }
    
    /** \copydoc DynamicAccumulatorChain::activateAll() */
    void activateAll()
    {
        this->next_.activateAll();
    }
    
    /** Return true if the statistic 'tag' is active, i.e. activate(std::string tag) or activate<TAG>() has been called. If the statistic is not in the accumulator chain a PreconditionViolation is thrown. (Note that alias names are not recognized.)
     */
    bool isActive(std::string tag) const
    {
        acc_detail::TagIsActive_Visitor v;
        vigra_precondition(isActiveImpl(tag, v),
            std::string("DynamicAccumulatorChainArray::isActive(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    /** %isActive\<TAG\>() returns true if statistic 'TAG' is active, i.e. activate(std::string tag) or activate<TAG>() has been called. If the statistic is not in the accumulator chain, true is returned. (?)
     */
    template <class TAG>
    bool isActive() const
    {
        return this->next_.template isActive<TAG>();
    }
    
    /** \copydoc DynamicAccumulatorChain::activeNames() */
    ArrayVector<std::string> activeNames() const
    {
        ArrayVector<std::string> res;
        for(unsigned k=0; k<DynamicAccumulatorChainArray::tagNames().size(); ++k)
            if(isActive(DynamicAccumulatorChainArray::tagNames()[k]))
                res.push_back(DynamicAccumulatorChainArray::tagNames()[k]);
        return res;
    }
    
    /** \copydoc DynamicAccumulatorChain::passesRequired() */
    unsigned int passesRequired() const
    {
        return this->next_.passesRequiredDynamic();
    }

  protected:
  
    bool activateImpl(std::string tag)
    {
        return acc_detail::ApplyVisitorToTag<AccumulatorTags>::exec(this->next_, 
                                         normalizeString(tag), acc_detail::ActivateTag_Visitor());
    }
    
    bool isActiveImpl(std::string tag, acc_detail::TagIsActive_Visitor & v) const
    {
        return acc_detail::ApplyVisitorToTag<AccumulatorTags>::exec(this->next_, normalizeString(tag), v);
    }
};

template <unsigned int N, class T1, class T2, class T3, class T4, class T5, class Selected>
class DynamicAccumulatorChainArray<CoupledArrays<N, T1, T2, T3, T4, T5>, Selected>
: public DynamicAccumulatorChainArray<typename CoupledArrays<N, T1, T2, T3, T4, T5>::HandleType, Selected>
{};

/****************************************************************************/
/*                                                                          */
/*                        generic access functions                          */
/*                                                                          */
/****************************************************************************/

template <class TAG>
struct Error__Attempt_to_access_inactive_statistic;

namespace acc_detail {

    // accumulator lookup rules: find the accumulator that implements TAG
    
    // When A does not implement TAG, continue search in A::InternalBaseType.
template <class TAG, class A, class FromTag=typename A::Tag>
struct LookupTagImpl
#ifndef DOXYGEN
: public LookupTagImpl<TAG, typename A::InternalBaseType>
#endif
{};

    // 'const A' is treated like A, except that the reference member is now const.
template <class TAG, class A, class FromTag>
struct LookupTagImpl<TAG, A const, FromTag>
: public LookupTagImpl<TAG, A>
{
    typedef typename LookupTagImpl<TAG, A>::type const & reference;
    typedef typename LookupTagImpl<TAG, A>::type const * pointer;
};

    // When A implements TAG, report its type and associated information.
template <class TAG, class A>
struct LookupTagImpl<TAG, A, TAG>
{
    typedef TAG Tag;
    typedef A type;
    typedef A & reference;
    typedef A * pointer;
    typedef typename A::value_type value_type;
    typedef typename A::result_type result_type;
};

    // Again, 'const A' is treated like A, except that the reference member is now const.
template <class TAG, class A>
struct LookupTagImpl<TAG, A const, TAG>
: public LookupTagImpl<TAG, A, TAG>
{
    typedef typename LookupTagImpl<TAG, A, TAG>::type const & reference;
    typedef typename LookupTagImpl<TAG, A, TAG>::type const * pointer;
};

    // Recursion termination: when we end up in AccumulatorEnd without finding a 
    // suitable A, we stop and report an error
template <class TAG, class A>
struct LookupTagImpl<TAG, A, AccumulatorEnd>
{
    typedef TAG Tag;
    typedef A type;
    typedef A & reference;
    typedef A * pointer;
    typedef Error__Attempt_to_access_inactive_statistic<TAG> value_type;
    typedef Error__Attempt_to_access_inactive_statistic<TAG> result_type;
};

    // ... except when we are actually looking for AccumulatorEnd
template <class A>
struct LookupTagImpl<AccumulatorEnd, A, AccumulatorEnd>
{
    typedef AccumulatorEnd Tag;
    typedef A type;
    typedef A & reference;
    typedef A * pointer;
    typedef void value_type;
    typedef void result_type;
};

    // ... or we are looking for a global statistic, in which case
    // we continue the serach via A::GlobalAccumulatorType, but remember that 
    // we are actually looking for a global tag. 
template <class TAG, class A>
struct LookupTagImpl<Global<TAG>, A, AccumulatorEnd>
: public LookupTagImpl<TAG, typename A::GlobalAccumulatorType>
{
    typedef Global<TAG> Tag;
};

    // When we encounter the LabelDispatch accumulator, we continue the
    // search via LabelDispatch::RegionAccumulatorChain by default
template <class TAG, class A>
struct LookupTagImpl<TAG, A, LabelDispatchTag>
: public LookupTagImpl<TAG, typename A::RegionAccumulatorChain>
{};

    // ... except when we are looking for a global statistic, in which case
    // we continue via LabelDispatch::GlobalAccumulatorChain, but remember that 
    // we are actually looking for a global tag.
template <class TAG, class A>
struct LookupTagImpl<Global<TAG>, A, LabelDispatchTag>
: public LookupTagImpl<TAG, typename A::GlobalAccumulatorChain>
{
    typedef Global<TAG> Tag;
};

    // ... or we are looking for the LabelDispatch accumulator itself
template <class A>
struct LookupTagImpl<LabelDispatchTag, A, LabelDispatchTag>
{
    typedef LabelDispatchTag Tag;
    typedef A type;
    typedef A & reference;
    typedef A * pointer;
    typedef void value_type;
    typedef void result_type;
};

} // namespace acc_detail

    // Lookup the accumulator in the chain A that implements the given TAG.
template <class Tag, class A>
struct LookupTag
: public acc_detail::LookupTagImpl<typename StandardizeTag<Tag>::type, A>
{};

    // Lookup the dependency TAG of the accumulator A.
    // This template ensures that dependencies are used with matching modifiers.
    // Specifically, if you search for Count as a dependency of Weighted<Mean>, the search
    // actually returns Weighted<Count>, wheras Count will be returned for plain Mean.
template <class Tag, class A, class TargetTag>
struct LookupDependency
: public acc_detail::LookupTagImpl<
       typename TransferModifiers<TargetTag, typename StandardizeTag<Tag>::type>::type, A>
{};
 

namespace acc_detail {

    // CastImpl applies the same rules as LookupTagImpl, but returns a reference to an 
    // accumulator instance rather than an accumulator type
template <class Tag, class FromTag, class reference>
struct CastImpl
{
    template <class A>
    static reference exec(A & a)
    {
        return CastImpl<Tag, typename A::InternalBaseType::Tag, reference>::exec(a.next_);
    }
    
    template <class A>
    static reference exec(A & a, MultiArrayIndex label)
    {
        return CastImpl<Tag, typename A::InternalBaseType::Tag, reference>::exec(a.next_, label);
    }
};

template <class Tag, class reference>
struct CastImpl<Tag, Tag, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return const_cast<reference>(a);
    }
    
    template <class A>
    static reference exec(A & a, MultiArrayIndex)
    {
        vigra_precondition(false, 
            "getAccumulator(): region accumulators can only be queried for AccumulatorChainArray.");
        return a;
    }
};

template <class Tag, class reference>
struct CastImpl<Tag, AccumulatorEnd, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return a;
    }
    
    template <class A>
    static reference exec(A & a, MultiArrayIndex)
    {
        return a;
    }
};

template <class Tag, class reference>
struct CastImpl<Global<Tag>, AccumulatorEnd, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return CastImpl<Tag, typename A::GlobalAccumulatorType::Tag, reference>::exec(*a.globalAccumulator_.pointer_);
    }
};

template <class reference>
struct CastImpl<AccumulatorEnd, AccumulatorEnd, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return a;
    }
    
    template <class A>
    static reference exec(A & a, MultiArrayIndex)
    {
        return a;
    }
};

template <class Tag, class reference>
struct CastImpl<Tag, LabelDispatchTag, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        vigra_precondition(false, 
            "getAccumulator(): a region label is required when a region accumulator is queried.");
        return CastImpl<Tag, typename A::RegionAccumulatorChain::Tag, reference>::exec(a.regions_[0]);
    }
    
    template <class A>
    static reference exec(A & a, MultiArrayIndex label)
    {
        return CastImpl<Tag, typename A::RegionAccumulatorChain::Tag, reference>::exec(a.regions_[label]);
    }
};

template <class Tag, class reference>
struct CastImpl<Global<Tag>, LabelDispatchTag, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return CastImpl<Tag, typename A::GlobalAccumulatorChain::Tag, reference>::exec(a.next_);
    }
};

template <class reference>
struct CastImpl<LabelDispatchTag, LabelDispatchTag, reference>
{
    template <class A>
    static reference exec(A & a)
    {
        return a;
    }
};

} // namespace acc_detail

    // Get a reference to the accumulator TAG in the accumulator chain A
/** Get a reference to the accumulator 'TAG' in the accumulator chain 'a'. This can be useful for example to update a certain accumulator with data, set individual options or get information about a certain accumulator.\n
Example of use (set options):
\code
    vigra::MultiArray<2, double> data(...);   
    typedef UserRangeHistogram<40> SomeHistogram;   //binCount set at compile time
    typedef UserRangeHistogram<0> SomeHistogram2; // binCount must be set at run-time
    AccumulatorChain<DataType, Select<SomeHistogram, SomeHistogram2> > a;
    
    getAccumulator<SomeHistogram>(a).setMinMax(0.1, 0.9);
    getAccumulator<SomeHistogram2>(a).setMinMax(0.0, 1.0);

    extractFeatures(data.begin(), data.end(), a);
\endcode

Example of use (get information):
\code
  vigra::MultiArray<2, double> data(...));
  AccumulatorChain<double, Select<Mean, Skewness> > a;

  std::cout << "passes required for all statistics: " << a.passesRequired() << std::endl; //skewness needs two passes
  std::cout << "passes required by Mean: " << getAccumulator<Mean>(a).passesRequired() << std::endl;
\endcode
See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
template <class TAG, class A>
inline typename LookupTag<TAG, A>::reference
getAccumulator(A & a)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::reference reference;
    return acc_detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a);
}

    // Get a reference to the accumulator TAG for region 'label' in the accumulator chain A
/** Get a reference to the accumulator 'TAG' for region 'label' in the accumulator chain 'a'.
*/
template <class TAG, class A>
inline typename LookupTag<TAG, A>::reference
getAccumulator(A & a, MultiArrayIndex label)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::reference reference;
    return acc_detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a, label);
}

    // get the result of the accumulator specified by TAG
/** Get the result of the accumulator 'TAG' in the accumulator chain 'a'.\n
Example of use:
\code
    vigra::MultiArray<2, double> data(...);
    AccumulatorChain<DataType, Select<Variance, Mean, StdDev> > a;
    extractFeatures(data.begin(), data.end(), a); 
    double mean = get<Mean>(a);
\endcode
See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
template <class TAG, class A>
inline typename LookupTag<TAG, A>::result_type
get(A const & a)
{
    return getAccumulator<TAG>(a).get();
}

    // get the result of the accumulator TAG for region 'label'
/** Get the result of the accumulator 'TAG' for region 'label' in the accumulator chain 'a'.\n
Example of use:
\code
    vigra::MultiArray<2, double> data(...);
    vigra::MultiArray<2, int> labels(...);
    typedef vigra::CoupledIteratorType<2, double, int>::type Iterator;
    typedef Iterator::value_type Handle;

    AccumulatorChainArray<Handle, 
        Select<DataArg<1>, LabelArg<2>, Mean, Variance> > a;

    Iterator start = createCoupledIterator(data, labels);
    Iterator end = start.getEndIterator();
    extractFeatures(start,end,a);

    double mean_of_region_1 = get<Mean>(a,1);
    double mean_of_background = get<Mean>(a,0);
\endcode
See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
template <class TAG, class A>
inline typename LookupTag<TAG, A>::result_type
get(A const & a, MultiArrayIndex label)
{
    return getAccumulator<TAG>(a, label).get();
}

    // Get the result of the accumulator specified by TAG without checking if the accumulator is active.
    // This must be used within an accumulator implementation to access dependencies because
    // it applies the approprate modifiers to the given TAG. It must not be used in other situations.
    // FIXME: is there a shorter name?
template <class TAG, class A>
inline typename LookupDependency<TAG, A>::result_type
getDependency(A const & a)
{
    typedef typename LookupDependency<TAG, A>::Tag StandardizedTag;
    typedef typename LookupDependency<TAG, A>::reference reference;
    return acc_detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a)();
}

    // activate the dynamic accumulator specified by Tag
/** Activate the dynamic accumulator 'Tag' in the dynamic accumulator chain 'a'. Same as a.activate<Tag>() (see DynamicAccumulatorChain::activate<Tag>() or DynamicAccumulatorChainArray::activate<Tag>()). For run-time activation use DynamicAccumulatorChain::activate(std::string tag) or DynamicAccumulatorChainArray::activate(std::string tag) instead.\n
See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
template <class Tag, class A>
inline void
activate(A & a)
{
    a.template activate<Tag>();
}

    // check if the dynamic accumulator specified by Tag is active
/** Check if the dynamic accumulator 'Tag' in the accumulator chain 'a' is active. Same as a.isActive<Tag>() (see DynamicAccumulatorChain::isActive<Tag>() or DynamicAccumulatorChainArray::isActive<Tag>()). At run-time, use DynamicAccumulatorChain::isActive(std::string tag) const or DynamicAccumulatorChainArray::isActive(std::string tag) const instead.\n
See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
template <class Tag, class A>
inline bool
isActive(A const & a)
{
    return a.template isActive<Tag>();
}

/****************************************************************************/
/*                                                                          */
/*                               generic loops                              */
/*                                                                          */
/****************************************************************************/

/** Generic loop to collect statistics from one or several arrays.

This function automatically performs as many passes over the data as necessary for the selected statistics. The basic version of <tt>extractFeatures()</tt> takes an iterator pair and a reference to an accumulator chain:
\code
namespace vigra { namespace acc {

    template <class ITERATOR, class ACCUMULATOR>
    void extractFeatures(ITERATOR start, ITERATOR end, ACCUMULATOR & a);
}}
\endcode
The <tt>ITERATOR</tt> can be any STL-conforming <i>forward iterator</i> (including raw pointers and \ref vigra::CoupledScanOrderIterator). The <tt>ACCUMULATOR</tt> must be instantiated with the <tt>ITERATOR</tt>'s <tt>value_type</tt> as its first template argument. For example, to use a raw pointer you write:
\code
    AccumulatorChain<double, Select<Mean, Variance> > a;

    double * start = ...,
           * end   = ...;
    extractFeatures(start, end, a);
\endcode
Similarly, you can use MultiArray's scan-order iterator:
\code    
    AccumulatorChain<TinyVector<float, 2>, Select<Mean, Variance> > a;

    MultiArray<3, TinyVector<float, 2> > data(...);
    extractFeatures(data.begin(), data.end(), a);
\endcode
An alternative syntax is used when you want to compute weighted or region statistics (or both). Then it is necessary to iterate over several arrays simultaneously. This fact is best conveyed to the accumulator via the helper class \ref vigra::CoupledArrays that is used as the accumulator's first template argument and holds the dimension and value types of the arrays involved. To actually compute the features, you then pass appropriate arrays to the <tt>extractfeatures()</tt> function directly. For example, region statistics can be obtained like this:
\code
    MultiArray<3, double> data(...);
    MultiArray<3, int> labels(...);

    AccumulatorChainArray<CoupledArrays<3, double, int>,
                          Select<DataArg<1>, LabelArg<2>, // where to look for data and region labels
                                 Mean, Variance> >        // what statistics to compute
        a;

    extractFeatures(data, labels, a);
\endcode
This form of <tt>extractFeatures()</tt> is supported for up to five arrays (although at most three are currently making sense in practice):
\code
namespace vigra { namespace acc {

    template <unsigned int N, class T1, class S1,
              class ACCUMULATOR>
    void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                         ACCUMULATOR & a);
                         
    ...

    template <unsigned int N, class T1, class S1,
                              class T2, class S2,
                              class T3, class S3,
                              class T4, class S4,
                              class T5, class S5,
              class ACCUMULATOR>
    void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                         MultiArrayView<N, T2, S2> const & a2, 
                         MultiArrayView<N, T3, S3> const & a3, 
                         MultiArrayView<N, T4, S4> const & a4, 
                         MultiArrayView<N, T5, S5> const & a5, 
                         ACCUMULATOR & a);
}}
\endcode
Of course, the number and types of the arrays specified in <tt>CoupledArrays</tt> must conform to the number and types of the arrays passed to <tt>extractFeatures()</tt>.

See \ref FeatureAccumulators for more information about feature computation via accumulators.
*/
doxygen_overloaded_function(template <...> void extractFeatures)


template <class ITERATOR, class ACCUMULATOR>
void extractFeatures(ITERATOR start, ITERATOR end, ACCUMULATOR & a)
{
    for(unsigned int k=1; k <= a.passesRequired(); ++k)
        for(ITERATOR i=start; i < end; ++i)
            a.updatePassN(*i, k);
}

template <unsigned int N, class T1, class S1,
          class ACCUMULATOR>
void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                     ACCUMULATOR & a)
{
    typedef typename CoupledIteratorType<N, T1>::type Iterator;
    Iterator start = createCoupledIterator(a1),
             end   = start.getEndIterator();
    extractFeatures(start, end, a);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class ACCUMULATOR>
void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                     MultiArrayView<N, T2, S2> const & a2, 
                     ACCUMULATOR & a)
{
    typedef typename CoupledIteratorType<N, T1, T2>::type Iterator;
    Iterator start = createCoupledIterator(a1, a2),
             end   = start.getEndIterator();
    extractFeatures(start, end, a);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3,
          class ACCUMULATOR>
void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                     MultiArrayView<N, T2, S2> const & a2, 
                     MultiArrayView<N, T3, S3> const & a3, 
                     ACCUMULATOR & a)
{
    typedef typename CoupledIteratorType<N, T1, T2, T3>::type Iterator;
    Iterator start = createCoupledIterator(a1, a2, a3),
             end   = start.getEndIterator();
    extractFeatures(start, end, a);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3,
                          class T4, class S4,
          class ACCUMULATOR>
void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                     MultiArrayView<N, T2, S2> const & a2, 
                     MultiArrayView<N, T3, S3> const & a3, 
                     MultiArrayView<N, T4, S4> const & a4, 
                     ACCUMULATOR & a)
{
    typedef typename CoupledIteratorType<N, T1, T2, T3, T4>::type Iterator;
    Iterator start = createCoupledIterator(a1, a2, a3, a4),
             end   = start.getEndIterator();
    extractFeatures(start, end, a);
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3,
                          class T4, class S4,
                          class T5, class S5,
          class ACCUMULATOR>
void extractFeatures(MultiArrayView<N, T1, S1> const & a1, 
                     MultiArrayView<N, T2, S2> const & a2, 
                     MultiArrayView<N, T3, S3> const & a3, 
                     MultiArrayView<N, T4, S4> const & a4, 
                     MultiArrayView<N, T5, S5> const & a5, 
                     ACCUMULATOR & a)
{
    typedef typename CoupledIteratorType<N, T1, T2, T3, T4, T5>::type Iterator;
    Iterator start = createCoupledIterator(a1, a2, a3, a4, a5),
             end   = start.getEndIterator();
    extractFeatures(start, end, a);
}

/****************************************************************************/
/*                                                                          */
/*                          AccumulatorResultTraits                         */
/*                                                                          */
/****************************************************************************/

template <class T>
struct AccumulatorResultTraits
{
    typedef T                                       type;
    typedef T                                       element_type;
    typedef double                                  element_promote_type;
    typedef T                                       MinmaxType;
    typedef element_promote_type                    SumType;
    typedef element_promote_type                    FlatCovarianceType;
    typedef element_promote_type                    CovarianceType;
};

template <class T, int N>
struct AccumulatorResultTraits<TinyVector<T, N> >
{
    typedef TinyVector<T, N>                             type;
    typedef T                                            element_type;
    typedef double                                       element_promote_type;
    typedef TinyVector<T, N>                             MinmaxType;
    typedef TinyVector<element_promote_type, N>          SumType;
    typedef TinyVector<element_promote_type, N*(N+1)/2>  FlatCovarianceType;
    typedef Matrix<element_promote_type>                 CovarianceType;
};

// (?) beign change
template <class T, unsigned int RED_IDX, unsigned int GREEN_IDX, unsigned int BLUE_IDX>
struct AccumulatorResultTraits<RGBValue<T, RED_IDX, GREEN_IDX, BLUE_IDX> >
{
    typedef RGBValue<T>                                  type;
    typedef T                                            element_type;
    typedef double                                       element_promote_type;
    typedef RGBValue<T>                                  MinmaxType;
    typedef RGBValue<element_promote_type>               SumType;
    typedef TinyVector<element_promote_type, 3*(3+1)/2>  FlatCovarianceType;
    typedef Matrix<element_promote_type>                 CovarianceType;
};
// end change


template <unsigned int N, class T, class Stride>
struct AccumulatorResultTraits<MultiArrayView<N, T, Stride> >
{
    typedef MultiArrayView<N, T, Stride>            type;
    typedef T                                       element_type;
    typedef double                                  element_promote_type;
    typedef MultiArray<N, T>                        MinmaxType;
    typedef MultiArray<N, element_promote_type>     SumType;
    typedef MultiArray<1, element_promote_type>     FlatCovarianceType;
    typedef Matrix<element_promote_type>            CovarianceType;
};

template <unsigned int N, class T, class Alloc>
struct AccumulatorResultTraits<MultiArray<N, T, Alloc> >
{
    typedef MultiArrayView<N, T, Alloc>             type;
    typedef T                                       element_type;
    typedef double                                  element_promote_type;
    typedef MultiArray<N, T>                        MinmaxType;
    typedef MultiArray<N, element_promote_type>     SumType;
    typedef MultiArray<1, element_promote_type>     FlatCovarianceType;
    typedef Matrix<element_promote_type>            CovarianceType;
};

/****************************************************************************/
/*                                                                          */
/*                           modifier implementations                       */
/*                                                                          */
/****************************************************************************/

/** \brief Modifier. Compute statistic globally rather than per region. 

This modifier only works when labels are given (with (Dynamic)AccumulatorChainArray), in which case statistics are computed per-region by default.
*/
template <class TAG>
class Global
{
  public:
    typedef typename StandardizeTag<TAG>::type  TargetTag;
    typedef typename TargetTag::Dependencies    Dependencies;
    
    static std::string name() 
    { 
        return std::string("Global<") + TargetTag::name() + " >";
        // static const std::string n = std::string("Global<") + TargetTag::name() + " >";
        // return n;
    }
};

/** \brief Specifies index of data in CoupledHandle. 

    If AccumulatorChain is used with CoupledIterator, DataArg<INDEX> tells the accumulator chain which index of the Handle contains the data. (Coordinates are always index 0)
*/
template <int INDEX>
class DataArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("DataArg<") + asString(INDEX) + "> (internal)";
        // static const std::string n = std::string("DataArg<") + asString(INDEX) + "> (internal)";
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef DataArgTag Tag;
        typedef void value_type;
        typedef void result_type;

        static const int value = INDEX;
        static const unsigned int workInPass = 0;
    };
};

namespace acc_detail {

template <class T, int DEFAULT, class TAG, class IndexDefinition, 
          class TagFound=typename IndexDefinition::Tag>
struct HandleArgSelectorImpl
{
    static const int value = DEFAULT;
    typedef typename CoupledHandleCast<value, T>::type type;
    typedef typename CoupledHandleCast<value, T>::value_type value_type;
    static const int size = type::dimensions;
    
    template <class U, class NEXT>
    static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type const & 
    getHandle(CoupledHandle<U, NEXT> const & t)
    {
        return vigra::cast<value>(t);
    }
    
    template <class U, class NEXT>
    static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference 
    getValue(CoupledHandle<U, NEXT> const & t)
    {
        return vigra::get<value>(t);
    }
};

template <class T, int DEFAULT, class TAG, class IndexDefinition>
struct HandleArgSelectorImpl<T, DEFAULT, TAG, IndexDefinition, TAG>
{
    static const int value = IndexDefinition::value;
    typedef typename CoupledHandleCast<value, T>::type type;
    typedef typename CoupledHandleCast<value, T>::value_type value_type;
    static const int size = type::dimensions;
    
    template <class U, class NEXT>
    static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type const & 
    getHandle(CoupledHandle<U, NEXT> const & t)
    {
        return vigra::cast<value>(t);
    }
    
    template <class U, class NEXT>
    static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference 
    getValue(CoupledHandle<U, NEXT> const & t)
    {
        return vigra::get<value>(t);
    }
};

} // namespace acc_detail

template <class T, class CHAIN>
struct HandleArgSelector<T, LabelArgTag, CHAIN>
: public acc_detail::HandleArgSelectorImpl<T, 2, LabelArgTag,
                                           typename LookupTag<LabelArgTag, CHAIN>::type>
{};

template <class T, class CHAIN>
struct HandleArgSelector<T, DataArgTag, CHAIN>
: public acc_detail::HandleArgSelectorImpl<T, 1, DataArgTag,
                                           typename LookupTag<DataArgTag, CHAIN>::type>
{};

template <class T, class CHAIN>
struct HandleArgSelector<T, CoordArgTag, CHAIN>
: public acc_detail::HandleArgSelectorImpl<T, 0, CoordArgTag,
                                           typename LookupTag<CoordArgTag, CHAIN>::type>
{
    typedef acc_detail::HandleArgSelectorImpl<T, 0, CoordArgTag,
                         typename LookupTag<CoordArgTag, CHAIN>::type> base_type;
    typedef TinyVector<double, base_type::size> value_type;
};

// Tags are automatically wrapped with DataFromHandle if CoupledHandle used
template <class TAG>
class DataFromHandle
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef typename TargetTag::Dependencies Dependencies;
    
    static std::string name() 
    { 
        return std::string("DataFromHandle<") + TargetTag::name() + " > (internal)";
        // static const std::string n = std::string("DataFromHandle<") + TargetTag::name() + " > (internal)";
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename HandleArgSelector<T, DataArgTag, BASE>::value_type, BASE>
    {
        typedef HandleArgSelector<T, DataArgTag, BASE>   DataHandle;
        typedef typename DataHandle::value_type          input_type;
        typedef input_type const &                       argument_type;
        typedef argument_type                            first_argument_type;
        
        typedef typename TargetTag::template Impl<input_type, BASE> ImplType;
        
        using ImplType::reshape;
        
        template <class U, class NEXT>
        void reshape(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::reshape(acc_detail::shapeOf(DataHandle::getValue(t)));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::update(DataHandle::getValue(t));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            ImplType::update(DataHandle::getValue(t), weight);
        }
    };
};

/** \brief Modifier. Compute statistic from pixel coordinates rather than from pixel values. 

    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
template <class TAG>
class Coord
{
  public:
    typedef typename StandardizeTag<TAG>::type   TargetTag;
    typedef typename TargetTag::Dependencies     Dependencies;
    
    static std::string name() 
    { 
        return std::string("Coord<") + TargetTag::name() + " >";
        // static const std::string n = std::string("Coord<") + TargetTag::name() + " >";
        // return n;
    }
        
    template <class T, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename HandleArgSelector<T, CoordArgTag, BASE>::value_type, BASE>
    {
        typedef HandleArgSelector<T, CoordArgTag, BASE>   CoordHandle;
        typedef typename CoordHandle::value_type          input_type;
        typedef input_type const &                        argument_type;
        typedef argument_type                             first_argument_type;
        
        typedef typename TargetTag::template Impl<input_type, BASE> ImplType;
        
        input_type offset_;
        
        Impl()
        : offset_()
        {}
        
        void setCoordinateOffset(input_type const & offset)
        {
            offset_ = offset;
        }
        
        using ImplType::reshape;
        
        template <class U, class NEXT>
        void reshape(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::reshape(acc_detail::shapeOf(CoordHandle::getValue(t)));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::update(CoordHandle::getValue(t)+offset_);
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            ImplType::update(CoordHandle::getValue(t)+offset_, weight);
        }
    };
};

/** \brief Specifies index of data in CoupledHandle. 

    If AccumulatorChain is used with CoupledIterator, WeightArg<INDEX> tells the accumulator chain which index of the Handle contains the weights. (Note that coordinates are always index 0.)
*/
template <int INDEX>
class WeightArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("WeightArg<") + asString(INDEX) + "> (internal)";
        // static const std::string n = std::string("WeightArg<") + asString(INDEX) + "> (internal)";
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef WeightArgTag Tag;
        typedef void value_type;
        typedef void result_type;

        static const int value = INDEX;
        static const unsigned int workInPass = 0;
    };
};

/** \brief Compute weighted version of the statistic.
*/
template <class TAG>
class Weighted
{
  public:
    typedef typename StandardizeTag<TAG>::type   TargetTag;
    typedef typename TargetTag::Dependencies     Dependencies;
    
    static std::string name() 
    { 
        return std::string("Weighted<") + TargetTag::name() + " >";
        // static const std::string n = std::string("Weighted<") + TargetTag::name() + " >";
        // return n;
    }
    
    template <class IndexDefinition, class TagFound=typename IndexDefinition::Tag>
    struct WeightIndexSelector
    {
        template <class U, class NEXT>
        static double exec(CoupledHandle<U, NEXT> const & t)
        {
            return (double)*t; // default: CoupledHandle holds weights at the last (outermost) index 
        }
    };
    
    template <class IndexDefinition>
    struct WeightIndexSelector<IndexDefinition, WeightArgTag>
    {
        template <class U, class NEXT>
        static double exec(CoupledHandle<U, NEXT> const & t)
        {
            return (double)get<IndexDefinition::value>(t);
        }
    };
    
    template <class T, class BASE>
    struct Impl
    : public TargetTag::template Impl<T, BASE>
    {
        typedef typename TargetTag::template Impl<T, BASE> ImplType;
        
        typedef typename LookupTag<WeightArgTag, BASE>::type FindWeightIndex;
                
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::update(t, WeightIndexSelector<FindWeightIndex>::exec(t));
        }
    };
};

// Centralize by subtracting the mean and cache the result
class Centralize
{
  public:
    typedef Select<Mean> Dependencies;
    
    static std::string name() 
    { 
         return "Centralize (internal)";
        // static const std::string n("Centralize (internal)");
        // return n;
    }
   
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename AccumulatorResultTraits<U>::element_promote_type element_type;
        typedef typename AccumulatorResultTraits<U>::SumType              value_type;
        typedef value_type const &                                  result_type;

        mutable value_type value_;
        
        Impl()
        : value_()  // call default constructor explicitly to ensure zero initialization
        {}
        
        void reset()
        {
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        void update(U const & t) const
        {
            using namespace vigra::multi_math;
            value_ = t - getDependency<Mean>(*this);
        }
        
        void update(U const & t, double) const
        {
            update(t);
        }
        
        result_type operator()(U const & t) const
        {
            update(t);
            return value_;
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

/** \brief Modifier. Substract mean before computing statistic. 

Works in pass 2, %operator+=() not supported (merging not supported).
*/
template <class TAG>
class Central
{
  public:
    typedef typename StandardizeTag<TAG>::type                    TargetTag;
    typedef Select<Centralize, typename TargetTag::Dependencies>  Dependencies;
    
    static std::string name() 
    { 
        return std::string("Central<") + TargetTag::name() + " >";
        // static const std::string n = std::string("Central<") + TargetTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE>
    {
        typedef typename TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE> ImplType;
        
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "Central<...>::operator+=(): not supported.");
        }
    
        template <class T>
        void update(T const & t)
        {
            ImplType::update(getDependency<Centralize>(*this));
        }
        
        template <class T>
        void update(T const & t, double weight)
        {
            ImplType::update(getDependency<Centralize>(*this), weight);
        }
    };
};

    // alternative implementation without caching 
    //
// template <class TAG>
// class Central
// {
  // public:
    // typedef typename StandardizeTag<TAG>::type TargetTag;
    // typedef TypeList<Mean, typename TransferModifiers<Central<TargetTag>, typename TargetTag::Dependencies::type>::type> Dependencies;
    
    // template <class U, class BASE>
    // struct Impl
    // : public TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE>
    // {
        // typedef typename TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE> ImplType;
        
        // static const unsigned int workInPass = 2;
        
        // void operator+=(Impl const & o)
        // {
            // vigra_precondition(false,
                // "Central<...>::operator+=(): not supported.");
        // }
    
        // template <class T>
        // void update(T const & t)
        // {
            // ImplType::update(t - getDependency<Mean>(*this));
        // }
        
        // template <class T>
        // void update(T const & t, double weight)
        // {
            // ImplType::update(t - getDependency<Mean>(*this), weight);
        // }
    // };
// };


class PrincipalProjection
{
  public:
    typedef Select<Centralize, Principal<CoordinateSystem> > Dependencies;
    
    static std::string name() 
    { 
        return "PrincipalProjection (internal)";
        // static const std::string n("PrincipalProjection (internal)");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename AccumulatorResultTraits<U>::element_promote_type element_type;
        typedef typename AccumulatorResultTraits<U>::SumType              value_type;
        typedef value_type const &                                  result_type;

        mutable value_type value_;
        
        Impl()
        : value_()  // call default constructor explicitly to ensure zero initialization
        {}
        
        void reset()
        {
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        void update(U const & t) const
        {
            for(unsigned int k=0; k<t.size(); ++k)
            {
                value_[k] = getDependency<Principal<CoordinateSystem> >(*this)(0, k)*getDependency<Centralize>(*this)[0];
                for(unsigned int d=1; d<t.size(); ++d)
                    value_[k] += getDependency<Principal<CoordinateSystem> >(*this)(d, k)*getDependency<Centralize>(*this)[d];
            }
        }
        
        void update(U const & t, double) const
        {
            update(t);
        }
        
        result_type operator()(U const & t) const
        {
            getAccumulator<Centralize>(*this).update(t);
            update(t);
            return value_;
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

/** \brief Modifier. Project onto PCA eigenvectors.

    Works in pass 2, %operator+=() not supported (merging not supported).
*/
template <class TAG>
class Principal
{
  public:
    typedef typename StandardizeTag<TAG>::type                             TargetTag;
    typedef Select<PrincipalProjection, typename TargetTag::Dependencies>  Dependencies;
    
    static std::string name() 
    { 
        return std::string("Principal<") + TargetTag::name() + " >";
        // static const std::string n = std::string("Principal<") + TargetTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE>
    {
        typedef typename TargetTag::template Impl<typename AccumulatorResultTraits<U>::SumType, BASE> ImplType;
        
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "Principal<...>::operator+=(): not supported.");
        }
    
        template <class T>
        void update(T const & t)
        {
            ImplType::update(getDependency<PrincipalProjection>(*this));
        }
        
        template <class T>
        void update(T const & t, double weight)
        {
            ImplType::update(getDependency<PrincipalProjection>(*this), weight);
        }
    };
};

/*
important notes on modifiers:
 * upon accumulator creation, modifiers are reordered so that data preparation is innermost, 
   and data access is outermost, e.g.:
        Coord<DivideByCount<Principal<PowerSum<2> > > >
 * modifiers are automatically transfered to dependencies as appropriate
 * modifiers for lookup (getAccumulator and get) of dependent accumulators are automatically adjusted
 * modifiers must adjust workInPass for the contained accumulator as appropriate
 * we may implement convenience versions of Select that apply a modifier to all 
   contained tags at once
 * weighted accumulators have their own Count object when used together
   with unweighted ones (this is as yet untested - FIXME)
 * certain accumulators must remain unchanged when wrapped in certain modifiers: 
    * Count: always except for Weighted<Count> and CoordWeighted<Count>
    * Sum: data preparation modifiers
    * FlatScatterMatrixImpl, CovarianceEigensystemImpl: Principal and Whitened
 * will it be useful to implement initPass<N>() or finalizePass<N>() ?
*/

/****************************************************************************/
/*                                                                          */
/*                        the actual accumulators                           */
/*                                                                          */
/****************************************************************************/

/** \brief Basic statistic. Identity matrix of appropriate size.
*/
class CoordinateSystem
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "CoordinateSystem";
        // static const std::string n("CoordinateSystem");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef double              element_type;
        typedef Matrix<double>      value_type;
        typedef value_type const &  result_type;

        value_type value_;
        
        Impl()
        : value_()  // call default constructor explicitly to ensure zero initialization
        {}
        
        void reset()
        {
            value_ = element_type();
        }

        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

template <class BASE, class T, 
          class ElementType=typename AccumulatorResultTraits<T>::element_promote_type, 
          class SumType=typename AccumulatorResultTraits<T>::SumType>
struct SumBaseImpl
: public BASE
{
    typedef ElementType         element_type;
    typedef SumType             value_type;
    typedef value_type const &  result_type;

    value_type value_;
    
    SumBaseImpl()
    : value_()  // call default constructor explicitly to ensure zero initialization
    {}
    
    void reset()
    {
        value_ = element_type();
    }

    template <class Shape>
    void reshape(Shape const & s)
    {
        acc_detail::reshapeImpl(value_, s);
    }
    
    void operator+=(SumBaseImpl const & o)
    {
        value_ += o.value_;
    }

    result_type operator()() const
    {
        return value_;
    }
};

// Count
template <>
class PowerSum<0>
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "PowerSum<0>";
        // static const std::string n("PowerSum<0>");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, T, double, double>
    {
        void update(T const & t)
        {
            ++this->value_;
        }
        
        void update(T const & t, double weight)
        {
            this->value_ += weight;
        }
    };
};

// Sum
template <>
class PowerSum<1>
{
  public:
    typedef Select<> Dependencies;
     
    static std::string name() 
    { 
        return "PowerSum<1>";
        // static const std::string n("PowerSum<1>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        void update(U const & t)
        {
            this->value_ += t;
        }
        
        void update(U const & t, double weight)
        {
            using namespace multi_math;

            this->value_ += weight*t;
        }
    };
};

/** \brief Basic statistic. PowerSum<N> =@f$ \sum_i x_i^N @f$

    Works in pass 1, %operator+=() supported (merging supported).
*/
template <unsigned N>
class PowerSum
{
  public:
    typedef Select<> Dependencies;
     
    static std::string name() 
    { 
        return std::string("PowerSum<") + asString(N) + ">";
        // static const std::string n = std::string("PowerSum<") + asString(N) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        void update(U const & t)
        {
            using namespace vigra::multi_math;            
            this->value_ += pow(t, (int)N);
        }
        
        void update(U const & t, double weight)
        {
            using namespace vigra::multi_math;            
            this->value_ += weight*pow(t, (int)N);
        }
    };
};

template <>
class AbsPowerSum<1>
{
  public:
    typedef Select<> Dependencies;
     
    static std::string name() 
    { 
        return "AbsPowerSum<1>";
        // static const std::string n("AbsPowerSum<1>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        void update(U const & t)
        {
            using namespace vigra::multi_math;            
            this->value_ += abs(t);
        }
        
        void update(U const & t, double weight)
        {
            using namespace vigra::multi_math;            
            this->value_ += weight*abs(t);
        }
    };
};

/** \brief Basic statistic. AbsPowerSum<N> =@f$ \sum_i |x_i|^N @f$

    Works in pass 1, %operator+=() supported (merging supported).
*/
template <unsigned N>
class AbsPowerSum
{
  public:
    typedef Select<> Dependencies;
     
    static std::string name() 
    { 
        return std::string("AbsPowerSum<") + asString(N) + ">";
        // static const std::string n = std::string("AbsPowerSum<") + asString(N) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        void update(U const & t)
        {
            using namespace vigra::multi_math;            
            this->value_ += pow(abs(t), (int)N);
        }
        
        void update(U const & t, double weight)
        {
            using namespace vigra::multi_math;            
            this->value_ += weight*pow(abs(t), (int)N);
        }
    };
};

template <class BASE, class VALUE_TYPE, class U>
struct CachedResultBase
: public BASE
{
    typedef typename AccumulatorResultTraits<U>::element_type  element_type;
    typedef VALUE_TYPE                                         value_type;
    typedef value_type const &                                 result_type;

    mutable value_type value_;
    
    CachedResultBase()
    : value_()  // call default constructor explicitly to ensure zero initialization
    {}
    
    void reset()
    {
        value_ = element_type();
        this->setClean();
    }

    template <class Shape>
    void reshape(Shape const & s)
    {
        acc_detail::reshapeImpl(value_, s);
    }

    void operator+=(CachedResultBase const &)
    {
        this->setDirty();
    }

    void update(U const &)
    {
        this->setDirty();
    }
    
    void update(U const &, double)
    {
         this->setDirty();
    }
};

// cached Mean and Variance
/** \brief Modifier. Divide statistic by Count:  DivideByCount<TAG> = TAG / Count .
*/
template <class TAG>
class DivideByCount
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef Select<TargetTag, Count> Dependencies;
  
    static std::string name() 
    { 
        return std::string("DivideByCount<") + TargetTag::name() + " >";
        // static const std::string n = std::string("DivideByCount<") + TargetTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public CachedResultBase<BASE, typename LookupDependency<TargetTag, BASE>::value_type, U> 
    {
        typedef typename CachedResultBase<BASE, typename LookupDependency<TargetTag, BASE>::value_type, U>::result_type result_type;
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                using namespace multi_math;
                this->value_ = getDependency<TargetTag>(*this) / getDependency<Count>(*this);
                this->setClean();
            }
            return this->value_;
        }
    };
};

// UnbiasedVariance
/** \brief Modifier. Divide statistics by Count-1:  DivideUnbiased<TAG> = TAG / (Count-1)
*/
template <class TAG>
class DivideUnbiased
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef Select<TargetTag, Count> Dependencies;
      
    static std::string name() 
    { 
        return std::string("DivideUnbiased<") + TargetTag::name() + " >";
        // static const std::string n = std::string("DivideUnbiased<") + TargetTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<TargetTag, BASE>::value_type  value_type;
        typedef value_type                                       result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return getDependency<TargetTag>(*this) / (getDependency<Count>(*this) - 1.0);
        }
    };
};

// RootMeanSquares and StdDev
/** \brief Modifier. RootDivideByCount<TAG> = sqrt( TAG/Count )
*/
template <class TAG>
class RootDivideByCount
{
  public:
    typedef typename StandardizeTag<DivideByCount<TAG> >::type TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    static std::string name() 
    { 
        typedef typename StandardizeTag<TAG>::type InnerTag;
        return std::string("RootDivideByCount<") + InnerTag::name() + " >";
        // static const std::string n = std::string("RootDivideByCount<") + InnerTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<TargetTag, BASE>::value_type  value_type;
        typedef value_type                                       result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return sqrt(getDependency<TargetTag>(*this));
        }
    };
};

// UnbiasedStdDev
/** \brief Modifier. RootDivideUnbiased<TAG> = sqrt( TAG / (Count-1) )
*/
template <class TAG>
class RootDivideUnbiased
{
  public:
    typedef typename StandardizeTag<DivideUnbiased<TAG> >::type TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    static std::string name() 
    { 
        typedef typename StandardizeTag<TAG>::type InnerTag;
        return std::string("RootDivideUnbiased<") + InnerTag::name() + " >";
        // static const std::string n = std::string("RootDivideUnbiased<") + InnerTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<TargetTag, BASE>::value_type  value_type;
        typedef value_type                                       result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return sqrt(getDependency<TargetTag>(*this));
        }
    };
};

/** \brief Spezialization: works in pass 1, %operator+=() supported (merging supported).
*/
template <>
class Central<PowerSum<2> >
{
  public:
    typedef Select<Mean, Count> Dependencies;
     
    static std::string name() 
    { 
        return "Central<PowerSum<2> >";
        // static const std::string n("Central<PowerSum<2> >");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        void operator+=(Impl const & o)
        {
            using namespace vigra::multi_math;
            double n1 = getDependency<Count>(*this), n2 = getDependency<Count>(o);
            if(n1 == 0.0)
            {
                this->value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                this->value_ += o.value_ + n1 * n2 / (n1 + n2) * sq(getDependency<Mean>(*this) - getDependency<Mean>(o));
            }
        }
    
        void update(U const & t)
        {
            double n = getDependency<Count>(*this);
            if(n > 1.0)
            {
                using namespace vigra::multi_math;
                this->value_ += n / (n - 1.0) * sq(getDependency<Mean>(*this) - t);
            }
        }
        
        void update(U const & t, double weight)
        {
            double n = getDependency<Count>(*this);
            if(n > weight)
            {
                using namespace vigra::multi_math;
                this->value_ += n / (n - weight) * sq(getDependency<Mean>(*this) - t);
            }
        }
    };
};

/** \brief Specialization: works in pass 2, %operator+=() supported (merging supported).
*/
template <>
class Central<PowerSum<3> >
{
  public:
    typedef Select<Centralize, Count, Mean, Central<PowerSum<2> > > Dependencies;
     
    static std::string name() 
    { 
        return "Central<PowerSum<3> >";
        // static const std::string n("Central<PowerSum<3> >");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        typedef typename SumBaseImpl<BASE, U>::value_type value_type;

        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            typedef Central<PowerSum<2> > Sum2Tag;
            
            using namespace vigra::multi_math;
            double n1 = getDependency<Count>(*this), n2 = getDependency<Count>(o);
            if(n1 == 0.0)
            {
                this->value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                double n = n1 + n2;
                double weight = n1 * n2 * (n1 - n2) / sq(n);
                value_type delta = getDependency<Mean>(o) - getDependency<Mean>(*this);
                this->value_ += o.value_ + weight * pow(delta, 3) +
                               3.0 / n * delta * (n1 * getDependency<Sum2Tag>(o) - n2 * getDependency<Sum2Tag>(*this));
            }
        }
    
        void update(U const & t)
        {
            using namespace vigra::multi_math;            
            this->value_ += pow(getDependency<Centralize>(*this), 3);
        }
        
        void update(U const & t, double weight)
        {
            using namespace vigra::multi_math;            
            this->value_ += weight*pow(getDependency<Centralize>(*this), 3);
        }
    };
};
/** \brief Specialization: works in pass 2, %operator+=() supported (merging supported).
*/
template <>
class Central<PowerSum<4> >
{
  public:
    typedef Select<Centralize, Central<PowerSum<3> > > Dependencies;
     
    static std::string name() 
    { 
        return "Central<PowerSum<4> >";
        // static const std::string n("Central<PowerSum<4> >");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public SumBaseImpl<BASE, U>
    {
        typedef typename SumBaseImpl<BASE, U>::value_type value_type;

        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            typedef Central<PowerSum<2> > Sum2Tag;
            typedef Central<PowerSum<3> > Sum3Tag;

            using namespace vigra::multi_math;
            double n1 = getDependency<Count>(*this), n2 = getDependency<Count>(o);
            if(n1 == 0.0)
            {
                this->value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                double n = n1 + n2;
                double n1_2 = sq(n1);
                double n2_2 = sq(n2);
                double n_2 = sq(n);
                double weight = n1 * n2 * (n1_2 - n1*n2 + n2_2) / n_2 / n;
                value_type delta = getDependency<Mean>(o) - getDependency<Mean>(*this);
                this->value_ += o.value_ + weight * pow(delta, 4) +
                              6.0 / n_2 * sq(delta) * (n1_2 * getDependency<Sum2Tag>(o) + n2_2 * getDependency<Sum2Tag>(*this)) +
                              4.0 / n * delta * (n1 * getDependency<Sum3Tag>(o) - n2 * getDependency<Sum3Tag>(*this));
            }
        }
    
        void update(U const & t)
        {
            using namespace vigra::multi_math;            
            this->value_ += pow(getDependency<Centralize>(*this), 4);
        }
        
        void update(U const & t, double weight)
        {
            using namespace vigra::multi_math;            
            this->value_ += weight*pow(getDependency<Centralize>(*this), 4);
        }
    };
};

/** \brief Basic statistic. Skewness. 

    %Skewness =@f$ \frac{ \frac{1}{n}\sum_i (x_i-\hat{x})^3 }{ (\frac{1}{n}\sum_i (x_i-\hat{x})^2)^{3/2} } @f$ .
    Works in pass 2, %operator+=() supported (merging supported).
*/
class Skewness
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<3> > > Dependencies;
    
    static std::string name() 
    { 
        return "Skewness";
        // static const std::string n("Skewness");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<3> >, BASE>::value_type   value_type;
        typedef value_type                                                    result_type;

        result_type operator()() const
        {
            typedef Central<PowerSum<3> > Sum3;
            typedef Central<PowerSum<2> > Sum2;
        
            using namespace multi_math;
            return sqrt(getDependency<Count>(*this)) * getDependency<Sum3>(*this) / pow(getDependency<Sum2>(*this), 1.5);
        }
    };
};

/** \brief Basic statistic. Unbiased Skewness.

    Works in pass 2, %operator+=() supported (merging supported).
*/
class UnbiasedSkewness
{
  public:
    typedef Select<Skewness> Dependencies;
    
    static std::string name() 
    { 
        return "UnbiasedSkewness";
        // static const std::string n("UnbiasedSkewness");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<3> >, BASE>::value_type   value_type;
        typedef value_type                                                    result_type;

        result_type operator()() const
        {
            using namespace multi_math;
            double n = getDependency<Count>(*this);
            return sqrt(n*(n-1.0)) / (n - 2.0) * getDependency<Skewness>(*this);
        }
    };
};

/** \brief Basic statistic. Kurtosis. 

    %Kurtosis = @f$ \frac{ \frac{1}{n}\sum_i (x_i-\bar{x})^4 }{
    (\frac{1}{n} \sum_i(x_i-\bar{x})^2)^2 } - 3 @f$ . 
    Works in pass 2, %operator+=() supported (merging supported).
*/
class Kurtosis
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<4> > > Dependencies;
    
    static std::string name() 
    { 
        return "Kurtosis";
        // static const std::string n("Kurtosis");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<4> >, BASE>::value_type  value_type;
        typedef value_type                                                          result_type;

        result_type operator()() const
        {
            typedef Central<PowerSum<4> > Sum4;
            typedef Central<PowerSum<2> > Sum2;
        
            using namespace multi_math;
            return getDependency<Count>(*this) * getDependency<Sum4>(*this) / sq(getDependency<Sum2>(*this)) - 3.0;
        }
    };
};

/** \brief Basic statistic. Unbiased Kurtosis.

    Works in pass 2, %operator+=() supported (merging supported).
*/
class UnbiasedKurtosis
{
  public:
    typedef Select<Kurtosis> Dependencies;
    
    static std::string name() 
    { 
        return "UnbiasedKurtosis";
        // static const std::string n("UnbiasedKurtosis");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<4> >, BASE>::value_type value_type;
        typedef value_type                                                         result_type;

        result_type operator()() const
        {
            using namespace multi_math;
            double n = getDependency<Count>(*this);
            return (n-1.0)/((n-2.0)*(n-3.0))*((n+1.0)*getDependency<Kurtosis>(*this) + value_type(6.0));
        }
    };
};

namespace acc_detail {

template <class Scatter, class Sum>
void updateFlatScatterMatrix(Scatter & sc, Sum const & s, double w)
{
    int size = s.size();
    for(MultiArrayIndex j=0, k=0; j<size; ++j)
        for(MultiArrayIndex i=j; i<size; ++i, ++k)
            sc[k] += w*s[i]*s[j];
}

template <class Sum>
void updateFlatScatterMatrix(double & sc, Sum const & s, double w)
{
    sc += w*s*s;
}

template <class Cov, class Scatter>
void flatScatterMatrixToScatterMatrix(Cov & cov, Scatter const & sc)
{
    int size = cov.shape(0), k=0;
    for(MultiArrayIndex j=0; j<size; ++j)
    {
        cov(j,j) = sc[k++];
        for(MultiArrayIndex i=j+1; i<size; ++i)
        {
            cov(i,j) = sc[k++];
            cov(j,i) = cov(i,j);
        }
    }
}

template <class Scatter>
void flatScatterMatrixToScatterMatrix(double & cov, Scatter const & sc)
{
    cov = sc;
}

template <class Cov, class Scatter>
void flatScatterMatrixToCovariance(Cov & cov, Scatter const & sc, double n)
{
    int size = cov.shape(0), k=0;
    for(MultiArrayIndex j=0; j<size; ++j)
    {
        cov(j,j) = sc[k++] / n;
        for(MultiArrayIndex i=j+1; i<size; ++i)
        {
            cov(i,j) = sc[k++] / n;
            cov(j,i) = cov(i,j);
        }
    }
}

template <class Scatter>
void flatScatterMatrixToCovariance(double & cov, Scatter const & sc, double n)
{
    cov = sc / n;
}

} // namespace acc_detail

// we only store the flattened upper triangular part of the scatter matrix
/** \brief Basic statistic. Flattened uppter-triangular part of scatter matrix.

    Works in pass 1, %operator+=() supported (merging supported).
*/
class FlatScatterMatrix
{
  public:
    typedef Select<Mean, Count> Dependencies;
    
    static std::string name() 
    { 
        return "FlatScatterMatrix";
        // static const std::string n("FlatScatterMatrix");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<U>::FlatCovarianceType    value_type;
        typedef value_type const &                                   result_type;
       
        typedef typename AccumulatorResultTraits<U>::SumType        SumType;

        value_type value_;
        SumType     diff_;
        
        Impl()
        : value_(),  // call default constructor explicitly to ensure zero initialization
          diff_()
        {}
        
        void reset()
        {
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            acc_detail::reshapeImpl(value_, Shape1(size*(size+1)/2));
            acc_detail::reshapeImpl(diff_, s);
        }
        
        void operator+=(Impl const & o)
        {
            double n1 = getDependency<Count>(*this), n2 = getDependency<Count>(o);
            if(n1 == 0.0)
            {
                value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                using namespace vigra::multi_math;
                diff_ = getDependency<Mean>(*this) - getDependency<Mean>(o);
                acc_detail::updateFlatScatterMatrix(value_, diff_, n1 * n2 / (n1 + n2));
                value_ += o.value_;
            }
        }
    
        void update(U const & t)
        {
            compute(t);
        }
        
        void update(U const & t, double weight)
        {
            compute(t, weight);
        }
        
        result_type operator()() const
        {
            return value_;
        }
        
      private:
        void compute(U const & t, double weight = 1.0)
        {
            double n = getDependency<Count>(*this);
            if(n > weight)
            {
                using namespace vigra::multi_math;
                diff_ = getDependency<Mean>(*this) - t;
                acc_detail::updateFlatScatterMatrix(value_, diff_, n * weight / (n - weight));
            }
        }
    };
};

// Covariance
template <>
class DivideByCount<FlatScatterMatrix>
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    static std::string name() 
    { 
        return "DivideByCount<FlatScatterMatrix>";
        // static const std::string n("DivideByCount<FlatScatterMatrix>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public CachedResultBase<BASE, typename AccumulatorResultTraits<U>::CovarianceType, U>
    {
        typedef CachedResultBase<BASE, typename AccumulatorResultTraits<U>::CovarianceType, U> BaseType;      
        typedef typename BaseType::result_type result_type;
        
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            acc_detail::reshapeImpl(this->value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                acc_detail::flatScatterMatrixToCovariance(this->value_, getDependency<FlatScatterMatrix>(*this), getDependency<Count>(*this));
                this->setClean();
            }
            return this->value_;
        }
    };
};

// UnbiasedCovariance
template <>
class DivideUnbiased<FlatScatterMatrix>
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    static std::string name() 
    { 
        return "DivideUnbiased<FlatScatterMatrix>";
        // static const std::string n("DivideUnbiased<FlatScatterMatrix>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public CachedResultBase<BASE, typename AccumulatorResultTraits<U>::CovarianceType, U>
    {
        typedef CachedResultBase<BASE, typename AccumulatorResultTraits<U>::CovarianceType, U> BaseType;      
        typedef typename BaseType::result_type result_type;
        
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            acc_detail::reshapeImpl(this->value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                acc_detail::flatScatterMatrixToCovariance(this->value_, getDependency<FlatScatterMatrix>(*this), getDependency<Count>(*this) - 1.0);
                this->setClean();
            }
            return this->value_;
        }
    };
};

/** Basic statistic. ScatterMatrixEigensystem (?)
*/
class ScatterMatrixEigensystem
{
  public:
    typedef Select<FlatScatterMatrix> Dependencies;
    
    static std::string name() 
    { 
        return "ScatterMatrixEigensystem";
        // static const std::string n("ScatterMatrixEigensystem");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<U>::SumType               EigenvalueType;
        typedef typename AccumulatorResultTraits<U>::CovarianceType        EigenvectorType;
        typedef std::pair<EigenvalueType, EigenvectorType>                 value_type;
        typedef value_type const &                                         result_type;

        mutable value_type value_;
        
        Impl()
        : value_()
        {}
        
        void operator+=(Impl const & o)
        {
            if(!acc_detail::hasDataImpl(value_.second))
            {
                acc_detail::copyShapeImpl(o.value_.first, value_.first);
                acc_detail::copyShapeImpl(o.value_.second, value_.second);
            }
            this->setDirty();
        }

        void update(U const &)
        {
            this->setDirty();
        }
        
        void update(U const &, double)
        {
             this->setDirty();
        }

        void reset()
        {
            value_.first = element_type();
            value_.second = element_type();
            this->setClean();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            acc_detail::reshapeImpl(value_.first, Shape1(size));
            acc_detail::reshapeImpl(value_.second, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                compute(getDependency<FlatScatterMatrix>(*this), value_.first, value_.second);
                this->setClean();
            }
            return value_;
        }
        
      private:
        template <class Flat, class EW, class EV>
        static void compute(Flat const & flatScatter, EW & ew, EV & ev)
        {
            EigenvectorType scatter(ev.shape());
            acc_detail::flatScatterMatrixToScatterMatrix(scatter, flatScatter);
            // create a view because EW could be a TinyVector
            MultiArrayView<2, element_type> ewview(Shape2(ev.shape(0), 1), &ew[0]);
            symmetricEigensystem(scatter, ewview, ev);
        }
        
        static void compute(double v, double & ew, double & ev)
        {
            ew = v;
            ev = 1.0;
        }
    };
};

// CovarianceEigensystem
template <>
class DivideByCount<ScatterMatrixEigensystem>
{
  public:
    typedef Select<ScatterMatrixEigensystem, Count> Dependencies;
    
    static std::string name() 
    { 
        return "DivideByCount<ScatterMatrixEigensystem>";
        // static const std::string n("DivideByCount<ScatterMatrixEigensystem>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<ScatterMatrixEigensystem, BASE>::type  SMImpl;
        typedef typename SMImpl::element_type                             element_type;
        typedef typename SMImpl::EigenvalueType                           EigenvalueType;
        typedef typename SMImpl::EigenvectorType                          EigenvectorType;
        typedef std::pair<EigenvalueType, EigenvectorType const &>        value_type;
        typedef value_type const &                                        result_type;

        mutable value_type value_;
        
        Impl()
        : value_(EigenvalueType(), BASE::template call_getDependency<ScatterMatrixEigensystem>().second)
        {}
        
        void operator+=(Impl const &)
        {
            this->setDirty();
        }

        void update(U const &)
        {
            this->setDirty();
        }
        
        void update(U const &, double)
        {
             this->setDirty();
        }

        void reset()
        {
            value_.first = element_type();
            this->setClean();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            acc_detail::reshapeImpl(value_.first, Shape2(size,1));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                value_.first = getDependency<ScatterMatrixEigensystem>(*this).first / getDependency<Count>(*this);
                this->setClean();
            }
            return value_;
        }
    };
};

// alternative implementation of CovarianceEigensystem - solve eigensystem directly
//
// template <>
// class DivideByCount<ScatterMatrixEigensystem>
// {
  // public:
    // typedef Select<Covariance> Dependencies;
    
    // template <class U, class BASE>
    // struct Impl
    // : public BASE
    // {
        // typedef typename AccumulatorResultTraits<U>::element_promote_type  element_type;
        // typedef typename AccumulatorResultTraits<U>::SumType               EigenvalueType;
        // typedef typename AccumulatorResultTraits<U>::CovarianceType        EigenvectorType;
        // typedef std::pair<EigenvalueType, EigenvectorType>                 value_type;
        // typedef value_type const &                                         result_type;

        // mutable value_type value_;
        
        // Impl()
        // : value_()
        // {}
        
        // void operator+=(Impl const &)
        // {
            // this->setDirty();
        // }

        // void update(U const &)
        // {
            // this->setDirty();
        // }
        
        // void update(U const &, double)
        // {
             // this->setDirty();
        // }

        // void reset()
        // {
            // value_.first = element_type();
            // value_.second = element_type();
            // this->setClean();
        // }
    
        // template <class Shape>
        // void reshape(Shape const & s)
        // {
            // int size = prod(s);
            // acc_detail::reshapeImpl(value_.first, Shape2(size,1));
            // acc_detail::reshapeImpl(value_.second, Shape2(size,size));
        // }
        
        // result_type operator()() const
        // {
            // if(this->isDirty())
            // {
                // compute(getDependency<Covariance>(*this), value_.first, value_.second);
                // this->setClean();
            // }
            // return value_;
        // }
        
      // private:
        // template <class Cov, class EW, class EV>
        // static void compute(Cov const & cov, EW & ew, EV & ev)
        // {
            // // create a view because EW could be a TinyVector
            // MultiArrayView<2, element_type> ewview(Shape2(cov.shape(0), 1), &ew[0]);
            // symmetricEigensystem(cov, ewview, ev);
        // }
        
        // static void compute(double cov, double & ew, double & ev)
        // {
            // ew = cov;
            // ev = 1.0;
        // }
    // };
// };

// covariance eigenvalues
/** \brief Specialization (covariance eigenvalues): works in pass 1, %operator+=() supported (merging).
*/
template <>
class Principal<PowerSum<2> >
{
  public:
    typedef Select<ScatterMatrixEigensystem> Dependencies;
     
    static std::string name() 
    { 
        return "Principal<PowerSum<2> >";
        // static const std::string n("Principal<PowerSum<2> >");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<ScatterMatrixEigensystem, BASE>::type::EigenvalueType value_type;
        typedef value_type const &                                                       result_type;
        
        result_type operator()() const
        {
            return getDependency<ScatterMatrixEigensystem>(*this).first;
        }
    };
};


// Principal<CoordinateSystem> == covariance eigenvectors
/** \brief Specialization (covariance eigenvectors): works in pass 1, %operator+=() supported (merging).
*/
template <>
class Principal<CoordinateSystem>
{
  public:
    typedef Select<ScatterMatrixEigensystem> Dependencies;
     
    static std::string name() 
    { 
        return "Principal<CoordinateSystem>";
        // static const std::string n("Principal<CoordinateSystem>");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupDependency<ScatterMatrixEigensystem, BASE>::type::EigenvectorType value_type;
        typedef value_type const &                                                        result_type;
        
        result_type operator()() const
        {
            return getDependency<ScatterMatrixEigensystem>(*this).second;
        }
    };
};

/** \brief Basic statistic. %Minimum value.

    Works in pass 1, %operator+=() supported (merging supported).
*/
class Minimum
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "Minimum";
        // static const std::string n("Minimum");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_type element_type;
        typedef typename AccumulatorResultTraits<U>::MinmaxType   value_type;
        typedef value_type const &                                result_type;

        value_type value_;
        
        Impl()
        {
            value_ = NumericTraits<element_type>::max();
        }
        
        void reset()
        {
            value_ = NumericTraits<element_type>::max();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s, NumericTraits<element_type>::max());
        }
        
        void operator+=(Impl const & o)
        {
            updateImpl(o.value_); // necessary because std::min causes ambiguous overload
        }
    
        void update(U const & t)
        {
            updateImpl(t);
        }
        
        void update(U const & t, double)
        {
            updateImpl(t);
        }
        
        result_type operator()() const
        {
            return value_;
        }
        
      private:
        template <class T>
        void updateImpl(T const & o)
        {
            using namespace multi_math;
            value_ = min(value_, o);
        }
        
        template <class T, class Alloc>
        void updateImpl(MultiArray<1, T, Alloc> const & o)
        {
            value_ = multi_math::min(value_, o);
        }
    };
};

/** \brief Basic statistic. %Maximum value.

    Works in pass 1, %operator+=() supported (merging supported).
*/
class Maximum
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "Maximum";
        // static const std::string n("Maximum");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_type element_type;
        typedef typename AccumulatorResultTraits<U>::MinmaxType   value_type;
        typedef value_type const &                                result_type;

        value_type value_;
        
        Impl()
        {
            value_ = NumericTraits<element_type>::min();
        }
        
        void reset()
        {
            value_ = NumericTraits<element_type>::min();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s, NumericTraits<element_type>::min());
        }
        
        void operator+=(Impl const & o)
        {
            updateImpl(o.value_); // necessary because std::max causes ambiguous overload
        }
    
        void update(U const & t)
        {
            updateImpl(t);
        }
        
        void update(U const & t, double)
        {
            updateImpl(t);
        }
        
        result_type operator()() const
        {
            return value_;
        }
        
      private:
        template <class T>
        void updateImpl(T const & o)
        {
            using namespace multi_math;
            value_ = max(value_, o);
        }
        
        template <class T, class Alloc>
        void updateImpl(MultiArray<1, T, Alloc> const & o)
        {
            value_ = multi_math::max(value_, o);
        }
    };
};

/** \brief Basic statistic. First data value seen of the object. 

    Usually used as <tt>Coord<FirstSeen></tt> (alias <tt>RegionAnchor</tt>) 
    which provides a well-defined anchor point for the region.
*/
class FirstSeen
{
  public:
    typedef Select<Count> Dependencies;
    
    static std::string name() 
    { 
        return "FirstSeen";
        // static const std::string n("FirstSeen");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_type element_type;
        typedef typename AccumulatorResultTraits<U>::MinmaxType   value_type;
        typedef value_type const &                                result_type;

        value_type value_;
        
        Impl()
        : value_()
        {}
        
        void reset()
        {
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        void operator+=(Impl const & o)
        {
            // FIXME: only works for Coord<FirstSeen>
            if(reverse(o.value_) < reverse(value_))
                value_ = o.value_;
        }
    
        void update(U const & t)
        {
            if(getDependency<Count>(*this) == 1)
                value_ = t;
        }
        
        void update(U const & t, double weight)
        {
            update(t);
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

/** \brief Return both the minimum and maximum in <tt>std::pair</tt>. 

    Usually used as <tt>Coord<Range></tt> (alias <tt>BoundingBox</tt>).
    Note that <tt>Range</tt> returns a closed interval, i.e. the upper
    limit is part of the range.
*/
class Range
{
  public:
    typedef Select<Minimum, Maximum> Dependencies;
    
    static std::string name() 
    { 
        return "Range";
        // static const std::string n("Range");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::MinmaxType   minmax_type;
        typedef std::pair<minmax_type, minmax_type>               value_type;
        typedef value_type                                        result_type;

        result_type operator()() const
        {
            return value_type(getDependency<Minimum>(*this), getDependency<Maximum>(*this));
        }
    };
};

/** \brief Basic statistic. Data value where weight assumes its minimal value. 

    Weights must be given. Coord<ArgMinWeight> gives coordinate where weight assumes its minimal value. Works in pass 1, %operator+=() supported (merging supported).
*/
class ArgMinWeight
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "ArgMinWeight";
        // static const std::string n("ArgMinWeight");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_type element_type;
        typedef typename AccumulatorResultTraits<U>::MinmaxType   value_type;
        typedef value_type const &                                result_type;

        double min_weight_;
        value_type value_;
        
        Impl()
        : min_weight_(NumericTraits<double>::max()),
          value_()
        {}
        
        void reset()
        {
            min_weight_ = NumericTraits<double>::max();
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        void operator+=(Impl const & o)
        {
            using namespace multi_math;
            if(o.min_weight_ < min_weight_)
            {
                min_weight_ = o.min_weight_;
                value_ = o.value_;
            }
        }
    
        void update(U const & t)
        {
            vigra_precondition(false, "ArgMinWeight::update() needs weights.");
        }
        
        void update(U const & t, double weight)
        {
            if(weight < min_weight_)
            {
                min_weight_ = weight;
                value_ = t;
            }
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

/** \brief Basic statistic. Data where weight assumes its maximal value. 

    Weights must be given. Coord<ArgMinWeight> gives coordinate where weight assumes its maximal value. Works in pass 1, %operator+=() supported (merging supported).
*/
class ArgMaxWeight
{
  public:
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return "ArgMaxWeight";
        // static const std::string n("ArgMaxWeight");
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<U>::element_type element_type;
        typedef typename AccumulatorResultTraits<U>::MinmaxType   value_type;
        typedef value_type const &                                result_type;

        double max_weight_;
        value_type value_;
        
        Impl()
        : max_weight_(NumericTraits<double>::min()),
          value_()
        {}
        
        void reset()
        {
            max_weight_ = NumericTraits<double>::min();
            value_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            acc_detail::reshapeImpl(value_, s);
        }
        
        void operator+=(Impl const & o)
        {
            using namespace multi_math;
            if(o.max_weight_ > max_weight_)
            {
                max_weight_ = o.max_weight_;
                value_ = o.value_;
            }
        }
    
        void update(U const & t)
        {
            vigra_precondition(false, "ArgMaxWeight::update() needs weights.");
        }
        
        void update(U const & t, double weight)
        {
            if(weight > max_weight_)
            {
                max_weight_ = weight;
                value_ = t;
            }
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};


template <class BASE, int BinCount>
class HistogramBase
: public BASE
{
  public:
  
    typedef double                        element_type;
    typedef TinyVector<double, BinCount>  value_type;
    typedef value_type const &            result_type;
    
    value_type value_;
    double left_outliers, right_outliers;
    
    HistogramBase()
    : value_(),
      left_outliers(), 
      right_outliers()
    {}
    
    void reset()
    {
        value_ = element_type();
        left_outliers = 0.0;
        right_outliers = 0.0;
    }

    void operator+=(HistogramBase const & o)
    {
        value_ += o.value_;
        left_outliers += o.left_outliers;
        right_outliers += o.right_outliers;
    }
        
    result_type operator()() const
    {
        return value_;
    }
};

template <class BASE>
class HistogramBase<BASE, 0>
: public BASE
{
  public:
  
    typedef double                        element_type;
    typedef MultiArray<1, double>         value_type;
    typedef value_type const &            result_type;
    
    value_type value_;
    double left_outliers, right_outliers;
    
    HistogramBase()
    : value_(),
      left_outliers(), 
      right_outliers()
    {}
    
    void reset()
    {
        value_ = element_type();
        left_outliers = 0.0;
        right_outliers = 0.0;
    }
    
    void operator+=(HistogramBase const & o)
    {
        if(value_.size() == 0)
        {
            value_ = o.value_;
        }
        else if(o.value_.size() > 0)
        {
            vigra_precondition(value_.size() == o.value_.size(),
                "HistogramBase::operator+=(): bin counts must be equal.");
            value_ += o.value_;
        }
        left_outliers += o.left_outliers;
        right_outliers += o.right_outliers;
    }
        
    void setBinCount(int binCount)
    {
        vigra_precondition(binCount > 0,
            "HistogramBase:.setBinCount(): binCount > 0 required.");
        value_type(Shape1(binCount)).swap(value_);
    }

    result_type operator()() const
    {
        return value_;
    }
};

template <class BASE, int BinCount, class U=typename BASE::input_type>
class RangeHistogramBase
: public HistogramBase<BASE, BinCount>
{
  public:
    double scale_, offset_, inverse_scale_;
    
    RangeHistogramBase()
    : scale_(),
      offset_(), 
      inverse_scale_()
    {}
    
    void reset()
    {
        scale_ = 0.0;
        offset_ = 0.0;
        inverse_scale_ = 0.0;
        HistogramBase<BASE, BinCount>::reset();
    }

    void operator+=(RangeHistogramBase const & o)
    {
        vigra_precondition(scale_ == 0.0 || o.scale_ == 0.0 || (scale_ == o.scale_ && offset_ == o.offset_),
            "RangeHistogramBase::operator+=(): cannot merge histograms with different data mapping.");
        
        HistogramBase<BASE, BinCount>::operator+=(o);
        if(scale_ == 0.0)
        {
            scale_ = o.scale_;
            offset_ = o.offset_;
            inverse_scale_ = o.inverse_scale_;
        }
    }

    void update(U const & t)
    {
        update(t, 1.0);
    }
    
    void update(U const & t, double weight)
    {
        double m = mapItem(t);
        int index =  (m == (double)this->value_.size())
                       ? (int)m - 1
                       : (int)m;
        if(index < 0)
            this->left_outliers += weight;
        else if(index >= (int)this->value_.size())
            this->right_outliers += weight;
        else
            this->value_[index] += weight;
    }
    
    void setMinMax(double mi, double ma)
    {
        vigra_precondition(this->value_.size() > 0,
            "RangeHistogramBase::setMinMax(...): setBinCount(...) has not been called.");
        vigra_precondition(mi <= ma,
            "RangeHistogramBase::setMinMax(...): min <= max required.");
        if(mi == ma)
            ma += this->value_.size() * NumericTraits<double>::epsilon();
        offset_ = mi;
        scale_ = (double)this->value_.size() / (ma - mi);
        inverse_scale_ = 1.0 / scale_;
    }
    
    double mapItem(double t) const
    {
        return scale_ * (t - offset_);
    }
    
    double mapItemInverse(double t) const
    {
        return inverse_scale_ * t + offset_;
    }
    
    template <class ArrayLike>
    void computeStandardQuantiles(double minimum, double maximum, double count, 
                                  ArrayLike const & desiredQuantiles, ArrayLike & res) const
    {
        if(count == 0.0) {
            return;
        }
        
        ArrayVector<double> keypoints, cumhist;
        double mappedMinimum = mapItem(minimum);
        double mappedMaximum = mapItem(maximum);
        
        keypoints.push_back(mappedMinimum);
        cumhist.push_back(0.0);
        
        if(this->left_outliers > 0.0)
        {
            keypoints.push_back(0.0);
            cumhist.push_back(this->left_outliers);
        }
        
        int size = (int)this->value_.size();
        double cumulative = this->left_outliers;
        for(int k=0; k<size; ++k)
        {
            if(this->value_[k] > 0.0)
            {
                if(keypoints.back() <= k)
                {
                    keypoints.push_back(k);
                    cumhist.push_back(cumulative);
                }
                cumulative += this->value_[k];
                keypoints.push_back(k+1);
                cumhist.push_back(cumulative);
            }
        }
        
        if(this->right_outliers > 0.0)
        {
            if(keypoints.back() != size)
            {
                keypoints.push_back(size);
                cumhist.push_back(cumulative);
            }
            keypoints.push_back(mappedMaximum);
            cumhist.push_back(count);
        }
        else
        {
            keypoints.back() = mappedMaximum;
            cumhist.back() = count;
        }
        
        int quantile = 0, end = (int)desiredQuantiles.size();
        
        if(desiredQuantiles[0] == 0.0)
        {
            res[0] = minimum;
            ++quantile;
        }
        if(desiredQuantiles[end-1] == 1.0)
        {
            res[end-1] = maximum;
            --end;
        }
        
        int point = 0;
        double qcount = count * desiredQuantiles[quantile];
        while(quantile < end)
        {
            if(cumhist[point] < qcount && cumhist[point+1] >= qcount)
            {
                double t = (qcount - cumhist[point]) / (cumhist[point+1] - cumhist[point]) * (keypoints[point+1] - keypoints[point]);
                res[quantile] = mapItemInverse(t + keypoints[point]);
                ++quantile;
                qcount = count * desiredQuantiles[quantile];
            }
            else
            {
                ++point;
            }
        }
    }
};

/** \brief Histogram where data values are equal to bin indices.

    - If BinCount != 0, the return type of the accumulator is TinyVector<double, BinCount> .
    - If BinCount == 0, the return type of the accumulator is MultiArray<1, double> . BinCount can be set by calling getAccumulator<IntegerHistogram<0> >(acc_chain).setBinCount(bincount).  
    - Outliers can be accessed via getAccumulator<IntegerHistogram<Bincount>>(a).left_outliers and getAccumulator<...>(acc_chain).right_outliers.
    - Note that histogram options (for all histograms in the accumulator chain) can also be set by passing an instance of HistogramOptions to the accumulator chain via acc_chain.setHistogramOptions().
    Works in pass 1, %operator+=() supported (merging supported).
*/
template <int BinCount>
class IntegerHistogram
{
  public:
    
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("IntegerHistogram<") + asString(BinCount) + ">";
        // static const std::string n = std::string("IntegerHistogram<") + asString(BinCount) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public HistogramBase<BASE, BinCount>
    {
        void update(int index)
        {
            if(index < 0)
                ++this->left_outliers;
            else if(index >= (int)this->value_.size())
                ++this->right_outliers;
            else
                ++this->value_[index];
        }
        
        void update(int index, double weight)
        {
            // cannot compute quantile from weighted integer histograms,
            // so force people to use UserRangeHistogram or AutoRangeHistogram
            vigra_precondition(false, "IntegerHistogram::update(): weighted histograms not supported, use another histogram type.");
        }
    
        template <class ArrayLike>
        void computeStandardQuantiles(double minimum, double maximum, double count, 
                                      ArrayLike const & desiredQuantiles, ArrayLike & res) const
        {
            int quantile = 0, end = (int)desiredQuantiles.size();
            
            if(desiredQuantiles[0] == 0.0)
            {
                res[0] = minimum;
                ++quantile;
            }
            if(desiredQuantiles[end-1] == 1.0)
            {
                res[end-1] = maximum;
                --end;
            }
            
            count -= 1.0;
            int currentBin = 0, size = (int)this->value_.size();
            double cumulative1 = this->left_outliers,
                   cumulative2 = this->value_[currentBin] + cumulative1;
            
            // add a to the quantiles to account for the fact that counting
            // corresponds to 1-based indexing (one element == index 1)
            double qcount = desiredQuantiles[quantile]*count + 1.0;
            
            while(quantile < end)
            {
                if(cumulative2 == qcount)
                {
                    res[quantile] = currentBin;
                    ++quantile;
                    qcount = desiredQuantiles[quantile]*count + 1.0;
                }
                else if(cumulative2 > qcount)
                {
                    if(cumulative1 > qcount) // in left_outlier bin
                    {
                        res[quantile] = minimum;
                    }
                    if(cumulative1 + 1.0 > qcount) // between bins
                    {
                        res[quantile] = currentBin - 1 + qcount - std::floor(qcount);
                    }
                    else // standard case
                    {
                        res[quantile] = currentBin;
                    }
                    ++quantile;
                    qcount = desiredQuantiles[quantile]*count + 1.0;
                }
                else if(currentBin == size-1) // in right outlier bin
                {
                    res[quantile] = maximum;
                    ++quantile;
                    qcount = desiredQuantiles[quantile]*count + 1.0;
                }
                else
                {
                    ++currentBin;
                    cumulative1 = cumulative2;
                    cumulative2 += this->value_[currentBin];
                }
            }
        }
    };
};

/** \brief Histogram where user provides bounds for linear range mapping from values to indices.

    - If BinCount != 0, the return type of the accumulator is TinyVector<double, BinCount> .
    - If BinCount == 0, the return type of the accumulator is MultiArray<1, double> . BinCount can be set by calling getAccumulator<UserRangeHistogram<0> >(acc_chain).setBinCount(bincount).
    - Bounds for the mapping (min/max) must be set before seeing data by calling getAccumulator<UserRangeHistogram<BinCount> >.setMinMax(min, max).
    - Options can also be passed to the accumulator chain via an instance of HistogramOptions .
    - Works in pass 1, %operator+=() is supported (merging) if both histograms have the same data mapping.
    - Outliers can be accessed via getAccumulator<...>(a).left_outliers and getAccumulator<...>(a).right_outliers.
    - Note that histogram options (for all histograms in the accumulator chain) can also be set by passing an instance of HistogramOptions to the accumulator chain via acc_chain.setHistogramOptions().
*/
template <int BinCount>
class UserRangeHistogram
{
  public:
    
    typedef Select<> Dependencies;
    
    static std::string name() 
    { 
        return std::string("UserRangeHistogram<") + asString(BinCount) + ">";
        // static const std::string n = std::string("UserRangeHistogram<") + asString(BinCount) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public RangeHistogramBase<BASE, BinCount, U>
    {
        void update(U const & t)
        {
            update(t, 1.0);
        }
        
        void update(U const & t, double weight)
        {
            vigra_precondition(this->scale_ != 0.0,
                "UserRangeHistogram::update(): setMinMax(...) has not been called.");
                
            RangeHistogramBase<BASE, BinCount, U>::update(t, weight);
        }
    };
};

/** \brief Histogram where range mapping bounds are defined by minimum and maximum of data.

    - If BinCount != 0, the return type of the accumulator is TinyVector<double, BinCount> .
    - If BinCount == 0, the return type of the accumulator is MultiArray<1, double> . BinCount can be set by calling getAccumulator<AutoRangeHistogram>(acc_chain).setBinCount(bincount).
    - Becomes a UserRangeHistogram if min/max is set.
    - Works in pass 2, %operator+=() is supported (merging) if both histograms have the same data mapping.
    - Outliers can be accessed via getAccumulator<...>(acc_chain).left_outliers and getAccumulator<...>(acc_chain).right_outliers .
    - Note that histogram options (for all histograms in the accumulator chain) can also be set by passing an instance of HistogramOptions to the accumulator chain via acc_chain.setHistogramOptions().
*/
template <int BinCount>
class AutoRangeHistogram
{
  public:
    
    typedef Select<Minimum, Maximum> Dependencies;
    
    static std::string name() 
    { 
        return std::string("AutoRangeHistogram<") + asString(BinCount) + ">";
        // static const std::string n = std::string("AutoRangeHistogram<") + asString(BinCount) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public RangeHistogramBase<BASE, BinCount, U>
    {
        static const unsigned int workInPass = LookupDependency<Minimum, BASE>::type::workInPass + 1;
        
        void update(U const & t)
        {
            update(t, 1.0);
        }
        
        void update(U const & t, double weight)
        {
            if(this->scale_ == 0.0)
                this->setMinMax(getDependency<Minimum>(*this), getDependency<Maximum>(*this));
                
            RangeHistogramBase<BASE, BinCount, U>::update(t, weight);
        }
    };
};

/** \brief Like AutoRangeHistogram, but use global min/max rather than region min/max.

    - If BinCount != 0, the return type of the accumulator is TinyVector<double, BinCount> .
    - If BinCount == 0, the return type of the accumulator is MultiArray<1, double> . BinCount can be set by calling getAccumulator<GlobalRangeHistogram<0>>(acc_chain).setBinCount(bincount).
    - Becomes a UserRangeHistogram if min/max is set.
    - Works in pass 2, %operator+=() is supported (merging) if both histograms have the same data mapping.
    - Outliers can be accessed via getAccumulator<GlobalRangeHistogram<Bincount>>(acc_chain).left_outliers and getAccumulator<...>(acc_chain).right_outliers .
    - Histogram options (for all histograms in the accumulator chain) can also be set by passing an instance of HistogramOptions to the accumulator chain via acc_chain.setHistogramOptions().
*/
template <int BinCount>
class GlobalRangeHistogram
{
  public:
    
    typedef Select<Global<Minimum>, Global<Maximum>, Minimum, Maximum> Dependencies;
    
    static std::string name() 
    { 
        return std::string("GlobalRangeHistogram<") + asString(BinCount) + ">";
        // static const std::string n = std::string("GlobalRangeHistogram<") + asString(BinCount) + ">";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public RangeHistogramBase<BASE, BinCount, U>
    {
        static const unsigned int workInPass = LookupDependency<Minimum, BASE>::type::workInPass + 1;
        
        bool useLocalMinimax_;
        
        Impl()
        : useLocalMinimax_(false)
        {}
        
        void setRegionAutoInit(bool locally)
        {
            this->scale_ = 0.0;
            useLocalMinimax_ = locally;
        }
        
        void update(U const & t)
        {
            update(t, 1.0);
        }
        
        void update(U const & t, double weight)
        {
            if(this->scale_ == 0.0)
            {
                if(useLocalMinimax_)
                    this->setMinMax(getDependency<Minimum>(*this), getDependency<Maximum>(*this));
                else
                    this->setMinMax(getDependency<Global<Minimum> >(*this), getDependency<Global<Maximum> >(*this));
            }
            
            RangeHistogramBase<BASE, BinCount, U>::update(t, weight);
        }
    };
};

/** \brief Compute (0%, 10%, 25%, 50%, 75%, 90%, 100%) quantiles from given histogram.

    Return type is TinyVector<double, 7> . 
*/
template <class HistogramAccumulator> 
class StandardQuantiles
{
  public:
    
    typedef typename StandardizeTag<HistogramAccumulator>::type HistogramTag;
    typedef Select<HistogramTag, Minimum, Maximum, Count> Dependencies;

    static std::string name() 
    { 
        return std::string("StandardQuantiles<") + HistogramTag::name() + " >";
        // static const std::string n = std::string("StandardQuantiles<") + HistogramTag::name() + " >";
        // return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public CachedResultBase<BASE, TinyVector<double, 7>, U>
    {
        typedef typename CachedResultBase<BASE, TinyVector<double, 7>, U>::result_type result_type;
        typedef typename CachedResultBase<BASE, TinyVector<double, 7>, U>::value_type  value_type;
        
        static const unsigned int workInPass = LookupDependency<HistogramTag, BASE>::type::workInPass;
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                double desiredQuantiles[] = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
                getAccumulator<HistogramTag>(*this).computeStandardQuantiles(getDependency<Minimum>(*this), getDependency<Maximum>(*this), 
                                                                             getDependency<Count>(*this), value_type(desiredQuantiles), 
                                                                             this->value_);
                this->setClean();
            }
            return this->value_;
        }
    };
};

template <int N>
struct feature_RegionContour_can_only_be_computed_for_2D_arrays
: vigra::staticAssert::AssertBool<N==2>
{};

/** \brief Compute the contour of a 2D region. 

    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
class RegionContour
{
  public:
    typedef Select<Count> Dependencies;
    
    static std::string name() 
    { 
        return std::string("RegionContour");
        // static const std::string n = std::string("RegionContour");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef HandleArgSelector<T, LabelArgTag, BASE>               LabelHandle;
        typedef TinyVector<double, 2>                                 point_type;
        typedef Polygon<point_type>                                   value_type;
        typedef value_type const &                                    result_type;
        
        point_type offset_;
        value_type contour_;
        
        Impl()
        : offset_()
        , contour_()
        {}
        
        void setCoordinateOffset(point_type const & offset)
        {
            offset_ = offset;
        }
                
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            VIGRA_STATIC_ASSERT((feature_RegionContour_can_only_be_computed_for_2D_arrays<
                                 CoupledHandle<U, NEXT>::dimensions>));
            if(getDependency<Count>(*this) == 1)
            {
                contour_.clear();
                extractContour(LabelHandle::getHandle(t).arrayView(), t.point(), contour_);
                contour_ += offset_;
            }
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            update(t);
        }
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "RegionContour::operator+=(): RegionContour cannot be merged.");
        }
        
        result_type operator()() const
        {
            return contour_;
        }
    };
};


/** \brief Compute the perimeter of a 2D region. 

    This is the length of the polygon returned by RegionContour.

    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
class RegionPerimeter
{
  public:
    typedef Select<RegionContour> Dependencies;
    
    static std::string name() 
    { 
        return std::string("RegionPerimeter");
        // static const std::string n = std::string("RegionPerimeter");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef double       value_type;
        typedef value_type   result_type;
        
        result_type operator()() const
        {
            return getDependency<RegionContour>(*this).length();
        }
    };
};

/** \brief Compute the circularity of a 2D region. 

    The is the ratio between the perimeter of a circle with the same area as the 
    present region and the perimeter of the region, i.e. \f[c = \frac{2 \sqrt{\pi a}}{p} \f], where a and p are the area and length of the polygon returned by RegionContour.
    
    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
class RegionCircularity
{
  public:
    typedef Select<Count, RegionContour> Dependencies;
    
    static std::string name() 
    { 
        return std::string("RegionCircularity");
        // static const std::string n = std::string("RegionCircularity");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef double       value_type;
        typedef value_type   result_type;
        
        result_type operator()() const
        {
            return 2.0*sqrt(M_PI*getDependency<RegionContour>(*this).area()) / getDependency<RegionContour>(*this).length();
        }
    };
};

/** \brief Compute the eccentricity of a 2D region in terms of its prinipal radii. 

    Formula: \f[ e = \sqrt{1 - m^2 / M^2 } \f], where m and M are the minor and major principal radius.

    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
class RegionEccentricity
{
  public:
    typedef Select<RegionRadii> Dependencies;
    
    static std::string name() 
    { 
        return std::string("RegionEccentricity");
        // static const std::string n = std::string("RegionEccentricity");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef double       value_type;
        typedef value_type   result_type;
        
        result_type operator()() const
        {
            double M = getDependency<RegionRadii>(*this).front(),
                   m = getDependency<RegionRadii>(*this).back();
            return sqrt(1.0 - sq(m/M));
        }
    };
};

template <int N>
struct feature_ConvexHull_can_only_be_computed_for_2D_arrays
: vigra::staticAssert::AssertBool<N==2>
{};

/** \brief Compute the contour of a 2D region. 

    AccumulatorChain must be used with CoupledIterator in order to have access to pixel coordinates.
 */
class ConvexHull
{
  public:
    typedef Select<BoundingBox, RegionContour, RegionCenter> Dependencies;
    
    static std::string name() 
    { 
        return std::string("ConvexHull");
        // static const std::string n = std::string("ConvexHull");
        // return n;
    }
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int            workInPass = 2;
        
        typedef HandleArgSelector<T, LabelArgTag, BASE>               LabelHandle;
        typedef TinyVector<double, 2>                                 point_type;
        typedef Polygon<point_type>                                   polygon_type;
        typedef Impl                                                  value_type;
        typedef value_type const &                                    result_type;
        
        polygon_type convex_hull_;
        point_type input_center_, convex_hull_center_, defect_center_;
        double convexity_, rugosity_, mean_defect_displacement_,
               defect_area_mean_, defect_area_variance_, defect_area_skewness_, defect_area_kurtosis_;
        int convexity_defect_count_;
        ArrayVector<MultiArrayIndex> convexity_defect_area_;
        bool features_computed_;
        
        Impl()
        : convex_hull_()
        , input_center_()
        , convex_hull_center_()
        , defect_center_()
        , convexity_()
        , rugosity_()
        , mean_defect_displacement_()
        , defect_area_mean_()
        , defect_area_variance_()
        , defect_area_skewness_()
        , defect_area_kurtosis_()
        , convexity_defect_count_()
        , convexity_defect_area_()
        , features_computed_(false)
        {}
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            VIGRA_STATIC_ASSERT((feature_ConvexHull_can_only_be_computed_for_2D_arrays<
                                  CoupledHandle<U, NEXT>::dimensions>));
            if(!features_computed_)
            {
                using namespace functor;
                Shape2 start = getDependency<Coord<Minimum> >(*this),
                       stop  = getDependency<Coord<Maximum> >(*this) + Shape2(1);
                point_type offset(start);
                input_center_ = getDependency<RegionCenter>(*this);
                MultiArrayIndex label = LabelHandle::getValue(t);
                
                convex_hull_.clear();
                convexHull(getDependency<RegionContour>(*this), convex_hull_);
                convex_hull_center_ = centroid(convex_hull_);
                
                convexity_ = getDependency<RegionContour>(*this).area() / convex_hull_.area();
                rugosity_ = getDependency<RegionContour>(*this).length() / convex_hull_.length();
                
                MultiArray<2, UInt8> convex_hull_difference(stop-start);
                fillPolygon(convex_hull_ - offset, convex_hull_difference, 1);
                combineTwoMultiArrays(convex_hull_difference, 
                                      LabelHandle::getHandle(t).arrayView().subarray(start, stop),
                                      convex_hull_difference,
                                      ifThenElse(Arg2() == Param(label), Param(0), Arg1()));
                                      
                MultiArray<2, UInt32> convexity_defects(stop-start);
                convexity_defect_count_ = 
                   labelImageWithBackground(convex_hull_difference, convexity_defects, false, 0);
                
                if (convexity_defect_count_ != 0)
                {
                    AccumulatorChainArray<CoupledArrays<2, UInt32>,
                                          Select<LabelArg<1>, Count, RegionCenter> > convexity_defects_stats;
                    convexity_defects_stats.ignoreLabel(0);
                    extractFeatures(convexity_defects, convexity_defects_stats);
                    
                    double total_defect_area = 0.0;
                    mean_defect_displacement_ = 0.0;
                    defect_center_ = point_type();
                    for (int k = 1; k <= convexity_defect_count_; ++k)
                    {
                        double area = get<Count>(convexity_defects_stats, k);
                        point_type center = get<RegionCenter>(convexity_defects_stats, k) + offset;
                        
                        convexity_defect_area_.push_back(area);
                        total_defect_area += area;
                        defect_center_ += area*center;
                        mean_defect_displacement_ += area*norm(input_center_ - center);
                    }
                    sort(convexity_defect_area_.begin(), convexity_defect_area_.end(),
                         std::greater<MultiArrayIndex>());
                    mean_defect_displacement_ /= total_defect_area;
                    defect_center_ /= total_defect_area;
                    
                    AccumulatorChain<MultiArrayIndex,
                                     Select<Mean, UnbiasedVariance, UnbiasedSkewness, UnbiasedKurtosis> > defect_area_stats;
                    extractFeatures(convexity_defect_area_.begin(),
                                    convexity_defect_area_.end(), defect_area_stats);
                    
                    defect_area_mean_ = convexity_defect_count_ > 0
                        ? get<Mean>(defect_area_stats)
                        : 0.0;
                    defect_area_variance_ = convexity_defect_count_ > 1
                        ? get<UnbiasedVariance>(defect_area_stats)
                        : 0.0;
                    defect_area_skewness_ = convexity_defect_count_ > 2
                        ? get<UnbiasedSkewness>(defect_area_stats)
                        : 0.0;
                    defect_area_kurtosis_ = convexity_defect_count_ > 3
                        ? get<UnbiasedKurtosis>(defect_area_stats)
                        : 0.0;
                }
                /**********************************************/
                features_computed_ = true;
            }
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            update(t);
        }
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "ConvexHull::operator+=(): ConvexHull features cannot be merged.");
        }
        
        result_type operator()() const
        {
            return *this;
        }
        
        /*
         * Returns the convex hull polygon.
         */
        polygon_type const & hull() const
        {
            return convex_hull_;
        }

        /*
         * Returns the area enclosed by the input polygon.
         */
        double inputArea() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return getDependency<RegionContour>(*this).area();
        }
        
        /*
         * Returns the area enclosed by the convex hull polygon.
         */
        double hullArea() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return convex_hull_.area();
        }

        /*
         * Returns the perimeter of the input polygon.
         */
        double inputPerimeter() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return getDependency<RegionContour>(*this).length();
        }

        /*
         * Returns the perimeter of the convex hull polygon.
         */
        double hullPerimeter() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return convex_hull_.length();
        }
        
        /*
         * Center of the original region.
         */
        point_type const & inputCenter() const
        {
            return input_center_;
        }
        
        /*
         * Center of the region enclosed by the convex hull.
         */
        point_type const & hullCenter() const
        {
            return convex_hull_center_;
        }
        
        /*
         * Center of difference between the convex hull and the original region.
         */
        point_type const & convexityDefectCenter() const
        {
            return defect_center_;
        }

        /*
         * Returns the ratio between the input area and the convex hull area.
         * This is always <tt><= 1</tt>, and the smaller the value is, 
         * the less convex is the input polygon.
         */
        double convexity() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return convexity_;
        }

        /*
         * Returns the ratio between the input perimeter and the convex perimeter.
         * This is always <tt>>= 1</tt>, and the higher the value is, the less 
         * convex is the input polygon.
         */
        double rugosity() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return rugosity_;
        }

        /*
         * Returns the number of convexity defects (i.e. number of connected components
         * of the difference between convex hull and input region).
         */
        int convexityDefectCount() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return convexity_defect_count_;
        }

        /*
         * Returns the mean area of the convexity defects.
         */
        double convexityDefectAreaMean() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return defect_area_mean_;
        }

        /*
         * Returns the variance of the convexity defect areas.
         */
        double convexityDefectAreaVariance() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return defect_area_variance_;
        }

        /*
         * Returns the skewness of the convexity defect areas.
         */
        double convexityDefectAreaSkewness() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return defect_area_skewness_;
        }

        /*
         * Returns the kurtosis of the convexity defect areas.
         */
        double convexityDefectAreaKurtosis() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return defect_area_kurtosis_;
        }

        /*
         * Returns the mean distance between the defect areas and the center of 
         * the input region, weighted by the area of each defect region.
         */
        double meanDefectDisplacement() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return mean_defect_displacement_;
        }

        /*
         * Returns the areas of the convexity defect regions (ordered descending).
         */
        ArrayVector<MultiArrayIndex> const & defectAreaList() const
        {
            vigra_precondition(features_computed_,
                    "ConvexHull: features must be calculated first.");
            return convexity_defect_area_;
        }
    };
};


}} // namespace vigra::acc

#endif // VIGRA_ACCUMULATOR_HXX
