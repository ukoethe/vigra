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

#ifndef VIGRA_ACCUMULATOR_GRAMMAR_HXX
#define VIGRA_ACCUMULATOR_GRAMMAR_HXX

#ifdef _MSC_VER
#pragma warning (disable: 4503)
#endif

#include "config.hxx"
#include "metaprogramming.hxx"

namespace vigra {

namespace acc {

/**************************************************************************/
/*                                                                        */
/*                      irreducible basic accumulators                    */
/*                                                                        */
/**************************************************************************/
    
   
class CoordinateSystem;                        // returns an identity matrix of appropriate size

template <unsigned N> class PowerSum;          // sum over powers of values
template <unsigned N> class AbsPowerSum;       // sum over powers of absolute values

class Skewness;                                // skewness
class UnbiasedSkewness;                        // unbiased estimator for skewness
class Kurtosis;                                // excess kurtosis
class UnbiasedKurtosis;                        // unbiased estimator for excess kurtosis
class FlatScatterMatrix;                       // flattened upper-triangular part of the scatter matrix
class ScatterMatrixEigensystem;                // eigenvalues and eigenvectors of the scatter matrix

    // in all histogram classes: set bin count at runtime if BinCount == 0
template <int BinCount> class IntegerHistogram;      // use data values directly as bin indices
template <int BinCount> class UserRangeHistogram;    // set min/max explicitly at runtime
template <int BinCount> class AutoRangeHistogram;    // get min/max from accumulators
template <int BinCount> class GlobalRangeHistogram;  // like AutoRangeHistogram, but use global min/max rather than region min/max

class FirstSeen;                               // remember the first value seen
class Minimum;                                 // minimum
class Maximum;                                 // maximum
class Range;                                   // minimum and maximum as a <tt>std::pair</tt>
template <class Hist> class StandardQuantiles; // compute (min, 10%, 25%, 50%, 75%, 90%, max) quantiles from 
                                               // min/max accumulators and given histogram

class ArgMinWeight;                            // store the value (or coordinate) where weight was minimal
class ArgMaxWeight;                            // store the value (or coordinate) where weight was maximal

    // FIXME: not yet implemented
template <unsigned NDim> class MultiHistogram; // multi-dimensional histogram
                                               // (always specify number of bins at runtime)

class Centralize;                              // cache centralized values
class PrincipalProjection;                     // cache values after principal projection
    // FIXME: not yet implemented
class Whiten;                                  // cache values after whitening
class RangeMapping;                            // map value from [min, max] to another range and cache result (e.g. for histogram creation)

template <int INDEX>  class DataArg;           // specifiy the index of the data member in a CoupledHandle
template <int INDEX>  class WeightArg;         // specifiy the index of the weight member in a CoupledHandle
template <int INDEX>  class LabelArg;          // specifiy the index of the label member in a CoupledHandle
template <int INDEX>  class CoordArg;          // specifiy the index of the coord member in a CoupledHandle

class RegionContour;                           // compute the contour of a 2D region
class RegionPerimeter;                         // compute the perimeter of a 2D region
class RegionCircularity;                       // compare perimeter of a 2D region with a circle of same area
class RegionEccentricity;                      // ecentricity of a 2D region from major and minor axis

/* 
Quantiles other than minimum and maximum require more thought:
--------------------------------------------------------------
 * Exact quantiles can be found in time O(n) using recursive partitioning (quickselect),
   but this requires the data (or an auxiliary index array) to be re-arranged.
 * Exact quantiles can be found in time O(k*n) using recursive histogram refinement,
   were k = O(log(n)/log(BinCount)) is the expected number of required passes over the 
   data, provided that bins are filled approximately evenly. If data can be re-arranged,
   such that we remember the items corresponding to each bin, running time reduces to O(n) 
   (this is Tibshirani's 'binmedian' algorithm). For the median, Tibshirani proves that
   the initial histogram only needs to cover the interval [Mean-StdDev, Mean+StdDev].
 * Both strategies can be combined: perform k passes to reduce bin size to about 
   n/(BinCount)^k, and then perform quickselect on an auxiliary array of size 
   O(n/(BinCount)^k) which has been filled during the final pass.
 * Good approximate results can be obtained by early stopping of histogram refinement
   (Tibshirani's 'binapprox' algorithm). A 2-pass algorithm for the median achieves
   accuracy of StdDev/BinCount: Mean and StdDev are computed during pass 1, 
   and a histogram over [Mean-StdDev, Mean+StdDev] during pass 2.
 * A 1-pass approximation method is described in Chen et al. However, it assumes that
   samples arrive in random order which is usually not true in image data.
*/

/**************************************************************************/
/*                                                                        */
/*                  modifiers for composite accumulators                  */
/*                                                                        */
/**************************************************************************/

   // data normalization w.r.t. number of samples
template <class A>    class DivideByCount;       //  A / count
template <class A>    class RootDivideByCount;   //  sqrt(A / count)
template <class A>    class DivideUnbiased;      //  A / (count - 1)  
template <class A>    class RootDivideUnbiased;  //  sqrt(A / (count - 1))

    // data access
template <class A> class Coord;           // use pixel coordinate instead of pixel value (index 0 of CoupledHandle)
template <class A> class Weighted;        // use (value, weight) pairs (index 1 and 2 of CoupledHandle)
template <class A> class CoordWeighted;   // use (coord, weight) pairs(index 0 and end of CoupledHandle)
template <class A> class DataFromHandle;  // extract data from index 1 of a CoupledHandle

    // data preparation
template <class A> class Central;    // subtract mean
template <class A> class Principal;  // subtract mean and rotate to principal coordinates

    // FIXME: not implemented yet
template <class A> class Whitened;   // transform to principal coordinates and scale to unit variance

template <class A> class Global;     // compute statistic A globally rather than per region

/**************************************************************************/
/*                                                                        */
/*                   alias names for important features                   */
/*                                                                        */
/**************************************************************************/

/** \brief Alias. Count. */
typedef PowerSum<0>                                 Count;
/** \brief Alias. Sum. */
typedef PowerSum<1>                                 Sum;
/** \brief Alias. Sum of squares. */
typedef PowerSum<2>                                 SumOfSquares;

/** \brief Alias. Mean. */
typedef DivideByCount<Sum>                          Mean;
/** \brief Alias. Root mean square. */
typedef RootDivideByCount<SumOfSquares>             RootMeanSquares;

// desired pseudocode (unfortunately not legal in C++)
//
//     template <unsigned N> 
//     typedef DivideByCount<PowerSum<N> >          Moment;
//
// actual definition (desired behavior is realised by rules below)
//
/** \brief Alias. Moment<N>. */
template <unsigned N> class                         Moment;  
/** \brief Alias. CentralMoment<N>. */
template <unsigned N> class                         CentralMoment;  

/** \brief Alias. Sum of squared differences. */
typedef Central<PowerSum<2> >                       SumOfSquaredDifferences;
/** \brief Alias. Sum of squared differences. */
typedef SumOfSquaredDifferences                     SSD;

/** \brief Alias. Variance. */
typedef DivideByCount<Central<PowerSum<2> > >       Variance;
/** \brief Alias. Standard deviation. */
typedef RootDivideByCount<Central<PowerSum<2> > >   StdDev;
/** \brief Alias. Unbiased variance. */
typedef DivideUnbiased<Central<PowerSum<2> > >      UnbiasedVariance;
/** \brief Alias. Unbiased standard deviation. */
typedef RootDivideUnbiased<Central<PowerSum<2> > >  UnbiasedStdDev;

/** \brief Alias. Covariance. */
typedef DivideByCount<FlatScatterMatrix>            Covariance;
/** \brief Alias. Unbiased covariance. */
typedef DivideUnbiased<FlatScatterMatrix>           UnbiasedCovariance;
/** \brief Alias. Covariance eigensystem. */
typedef DivideByCount<ScatterMatrixEigensystem>     CovarianceEigensystem;

/** \brief Alias. Absolute sum. */
typedef AbsPowerSum<1>                              AbsSum;
/** \brief Alias. Sum of absolute differences. */
typedef Central<AbsSum>                             SumOfAbsDifferences;
/** \brief Alias. Mean absolute deviation. */
typedef DivideByCount<SumOfAbsDifferences>          MeanAbsoluteDeviation;

/** \brief Alias. Rectangle enclosing the region, as a <tt>std::pair</tt> of coordinates. */
typedef Coord<Range>                                BoundingBox;
/** \brief Alias. Anchor point (first point of the region seen by scan-order traversal. */
typedef Coord<FirstSeen>                            RegionAnchor;

/** \brief Alias. Region center. */
typedef Coord<Mean>                                 RegionCenter;
/** \brief Alias. Region radii. */
typedef Coord<Principal<StdDev> >                   RegionRadii;
/** \brief Alias. Region axes. */
typedef Coord<Principal<CoordinateSystem> >         RegionAxes;

/** \brief Alias. Center of mass. */
typedef Weighted<RegionCenter>                      CenterOfMass;
/** \brief Alias. Moments of inertia. */
typedef Weighted<Coord<Principal<Variance> > >      MomentsOfInertia;
/** \brief Alias. Axes of inertia. */
typedef Weighted<RegionAxes>                        AxesOfInertia;

/**************************************************************************/
/*                                                                        */
/*                        Tag standardization rules                       */
/*                                                                        */
/**************************************************************************/

namespace acc_detail {

template <class A>
struct ModifierRule
{
    typedef A type;
};

} // namespace acc_detail

template <class A>
struct Error___Tag_modifiers_of_same_kind_must_not_be_combined;

    // apply rules as long as the Tag type changes ...
template <class A, class S=typename acc_detail::ModifierRule<A>::type>
struct StandardizeTag
{
    typedef typename StandardizeTag<S>::type type;
};

    // ... and stop otherwise ...
template <class A>
struct StandardizeTag<A, A>
{
    typedef A type;
};

    // ... or fail when the tag spec was non-conforming 
template <class A, class B>
struct StandardizeTag<A, Error___Tag_modifiers_of_same_kind_must_not_be_combined<B> >
    : public Error___Tag_modifiers_of_same_kind_must_not_be_combined<B>
{};

namespace acc_detail {

    // Assign priorities to modifiers to determine their standard order (by ascending priority).
    // SubstitutionMask determines which modifiers must be automatically transferred to dependencies.
enum { MinPriority = 1, 
       AccumulatorPriority = 32,
       PrepareDataPriority = 16,
       NormalizePriority = 8,
       AccessDataPriority = 4,
       WeightingPriority = 2,
       GlobalPriority = 1,
       MaxPriority = 32,
       SubstitutionMask = PrepareDataPriority | AccessDataPriority | WeightingPriority | GlobalPriority };

template <class A>
struct ModifierPriority
{
    static const int value = AccumulatorPriority;
};

#define VIGRA_MODIFIER_PRIORITY(MODIFIER, VALUE) \
template <class A> \
struct ModifierPriority<MODIFIER<A> > \
{ \
    static const int value = VALUE; \
};

VIGRA_MODIFIER_PRIORITY(Global, GlobalPriority)

VIGRA_MODIFIER_PRIORITY(Weighted, WeightingPriority)

VIGRA_MODIFIER_PRIORITY(Coord, AccessDataPriority)
VIGRA_MODIFIER_PRIORITY(DataFromHandle, AccessDataPriority)

VIGRA_MODIFIER_PRIORITY(DivideByCount, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(RootDivideByCount, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(DivideUnbiased, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(RootDivideUnbiased, NormalizePriority)

VIGRA_MODIFIER_PRIORITY(Central, PrepareDataPriority)
VIGRA_MODIFIER_PRIORITY(Principal, PrepareDataPriority)
VIGRA_MODIFIER_PRIORITY(Whitened, PrepareDataPriority)

    // explicitly set priority for base accumulators that look like modifiers
VIGRA_MODIFIER_PRIORITY(StandardQuantiles, AccumulatorPriority)

#undef VIGRA_MODIFIER_PRIORITY

    // check if the tag A contains a modifier with TARGET_PRIORITY
template <class A, int TARGET_PRIORITY, int PRIORITY=ModifierPriority<A>::value>
struct HasModifierPriority
{
    typedef VigraFalseType type;
    static const bool value = false;
};

template <class A, int TARGET_PRIORITY>
struct HasModifierPriority<A, TARGET_PRIORITY, TARGET_PRIORITY>
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class A, template <class> class B, int TARGET_PRIORITY, int PRIORITY>
struct HasModifierPriority<B<A>, TARGET_PRIORITY, PRIORITY>
: public HasModifierPriority<A, TARGET_PRIORITY>
{};

template <class A, template <class> class B, int TARGET_PRIORITY>
struct HasModifierPriority<B<A>, TARGET_PRIORITY, TARGET_PRIORITY>
{
    typedef VigraTrueType type;
    static const bool value = true;
};

    // three-way compare 
template <class A, class B>
struct ModifierCompare
{
    static const int p1 = ModifierPriority<A>::value;
    static const int p2 = ModifierPriority<B>::value;
    static const int value = p1 < p2 
                                ? -1
                                : p2 < p1
                                     ? 1
                                     : 0;
};

template <class A>
struct ModifierCompareToInner;

template <class A, template <class> class B>
struct ModifierCompareToInner<B<A> >
: public ModifierCompare<B<A>, A>
{};

    // sort modifiers by ascending priority
template <class A, int compare=ModifierCompareToInner<A>::value>
struct ModifierOrder;

    // do nothing if the order is correct (compare == -1)
template <class A>
struct ModifierOrder<A, -1>
{
    typedef A type;
};

    // fail if there are two modifiers with the same priority (compare == 0)
template <class A, template <class> class B, template <class> class C>
struct ModifierOrder<C<B<A> >, 0>
{
    typedef Error___Tag_modifiers_of_same_kind_must_not_be_combined<C<B<A> > > type;
};

    // sort if the order is reversed (compare == 1)
template <class A, template <class> class B, template <class> class C>
struct ModifierOrder<C<B<A> >, 1>
{
    typedef B<C<A> > type;
};

#define VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(OUTER, INNER, RESULT) \
template <class A> \
struct ModifierOrder<OUTER<INNER<A > >, 0> \
{ \
    typedef RESULT<A > type; \
};

    // drop duplicates
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Central, Central, Central)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Principal, Principal, Principal)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Whitened, Whitened)

    // the strongest data preparation modifier takes precendence
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Principal, Central, Principal)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Central, Whitened)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Principal, Whitened)

    // Coord takes precendence over DataFromHandle
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(DataFromHandle, Coord, Coord)

#undef VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS

    // drop duplicates
template <class A, template <class> class B>
struct ModifierRule<B<B<A> > >
{
    typedef B<A> type;
};

template <class A, int PRIORITY=ModifierPriority<A>::value>
struct RecurseModifier;

template <class A, template <class> class B, int PRIORITY>
struct RecurseModifier<B<A>, PRIORITY>
{
    typedef typename ModifierOrder<B<typename StandardizeTag<A>::type> >::type type;
};

template <class A, template <class> class B>
struct RecurseModifier<B<A>, AccumulatorPriority>
{
    typedef B<A> type;
};

    // recurse down the modifier chain, but only of B is actually a modifier,
    // and not a templated base accumulator (i.e. do not recurse if B's
    // priority is 'AccumulatorPriority')
template <class A, template <class> class B>
struct ModifierRule<B<A> >
: public RecurseModifier<B<A> >
{};

    // reduce the SOURCE modifier to the TARGET modifier,
    // using the given TEMPLATE arguments
    // (this is a work-around for the lack of templated typedef in C++)
#define VIGRA_REDUCE_MODFIER(TEMPLATE, SOURCE, TARGET) \
template <TEMPLATE > \
struct ModifierRule<SOURCE > \
{ \
    typedef TARGET type; \
};

#define VIGRA_VOID

    // centralizing doesn't change the CoordinateSystem
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<CoordinateSystem>, CoordinateSystem)
    // whitened CoordinateSystem are the same as principal CoordinateSystem
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Whitened<CoordinateSystem>, Principal<CoordinateSystem>)

    // counting modified data is the same as counting data, except for weighted data and global counting
VIGRA_REDUCE_MODFIER(template <class> class A, A<Count>, Count)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Weighted<Count>, Weighted<Count>)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, CoordWeighted<Count>, Weighted<Count>)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Global<Count>, Global<Count>)

    // reduce aliases that typedef can't handle
VIGRA_REDUCE_MODFIER(unsigned N, Moment<N>, DivideByCount<PowerSum<N> >)
VIGRA_REDUCE_MODFIER(unsigned N, CentralMoment<N>, DivideByCount<Central<PowerSum<N> > >)
VIGRA_REDUCE_MODFIER(class A, CoordWeighted<A>, Weighted<Coord<A> >)

    // reduce statistics that are inherently centered
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<Centralize>, Centralize)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<Skewness>, Skewness)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<Kurtosis>, Kurtosis)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<FlatScatterMatrix>, FlatScatterMatrix)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Central<ScatterMatrixEigensystem>, ScatterMatrixEigensystem)

VIGRA_REDUCE_MODFIER(VIGRA_VOID, Principal<Centralize>, PrincipalProjection)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Whitened<Centralize>, Whiten)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Principal<PrincipalProjection>, PrincipalProjection)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Whitened<PrincipalProjection>, Whiten)
VIGRA_REDUCE_MODFIER(VIGRA_VOID, Whitened<Whiten>, Whiten)

    // ignore all modifiers of RegionContour and related features
VIGRA_REDUCE_MODFIER(template <class> class A, A<RegionContour>, RegionContour)
VIGRA_REDUCE_MODFIER(template <class> class A, A<RegionPerimeter>, RegionPerimeter)
VIGRA_REDUCE_MODFIER(template <class> class A, A<RegionCircularity>, RegionCircularity)
VIGRA_REDUCE_MODFIER(template <class> class A, A<RegionEccentricity>, RegionEccentricity)

    // reduce even absolute powers to plain powers
template <unsigned N>
struct ModifierRule<AbsPowerSum<N> >
{
    typedef typename IfBool<(N % 2 == 0), PowerSum<N>, AbsPowerSum<N> >::type type;
};

#undef VIGRA_VOID
#undef VIGRA_REDUCE_MODFIER

template <class A>
struct ShouldBeWeighted
{
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct ShouldBeWeighted<ArgMinWeight>
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <>
struct ShouldBeWeighted<ArgMaxWeight>
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class A, template <class> class B>
struct ShouldBeWeighted<B<A> >
: public ShouldBeWeighted<A>
{};

} // namespace acc_detail

template <class A>
struct IsCoordinateFeature
{
    typedef VigraFalseType type;
    static const bool value = false;
};

template <class A, template <class> class B>
struct IsCoordinateFeature<B<A> >
{
    typedef typename IsCoordinateFeature<A>::type type;
    static const bool value = IsCoordinateFeature<A>::value;
};

template <class A>
struct IsCoordinateFeature<Coord<A> >
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class A>
struct IsPrincipalFeature
{
    typedef VigraFalseType type;
    static const bool value = false;
};

template <class A, template <class> class B>
struct IsPrincipalFeature<B<A> >
{
    typedef typename IsPrincipalFeature<A>::type type;
    static const bool value = IsPrincipalFeature<A>::value;
};

template <class A>
struct IsPrincipalFeature<Principal<A> >
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class A>
struct IsPrincipalFeature<Whitened<A> >
{
    typedef VigraTrueType type;
    static const bool value = true;
};

/**************************************************************************/
/*                                                                        */
/*                           Tag transfer rules                           */
/*                                                                        */
/**************************************************************************/

namespace acc_detail {

template <class A>
struct DefaultModifier;

template <class A>
struct ModifierPriority<DefaultModifier<A> >
{
    static const int value = ModifierPriority<A>::value << 1;
};

template <class A, int TargetPriority, int Priority=ModifierPriority<A>::value>
struct InsertDefaultModifier
{
    typedef DefaultModifier<typename InsertDefaultModifier<A, (TargetPriority << 1)>::type> type;
};

template <class A, int TargetPriority>
struct InsertDefaultModifier<A, TargetPriority, TargetPriority>
{
    typedef A type;
};

template <class A, int TargetPriority, int Priority=ModifierPriority<A>::value>
struct TagLongForm;

template <class A, int TargetPriority>
struct TagLongForm<A, TargetPriority, MaxPriority>
{
    typedef typename InsertDefaultModifier<A, TargetPriority>::type type;
};

template <class A, template <class> class B, int TargetPriority>
struct TagLongForm<B<A>, TargetPriority, MaxPriority>
{
    typedef typename InsertDefaultModifier<B<A>, TargetPriority>::type type;
};

template <class A, template <class> class B, int TargetPriority, int Priority>
struct TagLongForm<B<A>, TargetPriority, Priority>
{
    typedef typename TagLongForm<A, (Priority << 1)>::type Inner;
    typedef typename InsertDefaultModifier<B<Inner>, TargetPriority>::type type;
};

template <class A, template <class> class B, int TargetPriority>
struct TagLongForm<B<A>, TargetPriority, TargetPriority>
{
    typedef typename TagLongForm<A, (TargetPriority << 1)>::type Inner;
    typedef B<Inner> type;
};

template <class A>
struct LongModifierRule
{
    typedef A type;
};

    // apply rules as long as the Tag type changes ...
template <class A, class S=typename LongModifierRule<A>::type>
struct StandardizeTagLongForm
{
    typedef typename StandardizeTagLongForm<S>::type type;
};

    // ... and stop otherwise ...
template <class A>
struct StandardizeTagLongForm<A, A>
{
    typedef A type;
};

template <class A, template <class> class B>
struct LongModifierRule<B<A> >
{
    typedef B<typename LongModifierRule<A>::type> type;
};

template <class A>
struct LongModifierRule<DefaultModifier<A> >
{
    typedef A type;
};

#define VIGRA_DROP_DATA_PREPARATION_MODIFIERS(SOURCE, TARGET) \
template <> \
struct LongModifierRule<SOURCE > \
{ \
    typedef TARGET type; \
};

VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Central<Sum>, Sum)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Principal<Sum>, Sum)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Whitened<Sum>, Sum)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Principal<FlatScatterMatrix>, FlatScatterMatrix)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Whitened<FlatScatterMatrix>, FlatScatterMatrix)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Principal<ScatterMatrixEigensystem>, ScatterMatrixEigensystem)
VIGRA_DROP_DATA_PREPARATION_MODIFIERS(Whitened<ScatterMatrixEigensystem>, ScatterMatrixEigensystem)

#undef VIGRA_DROP_DATA_PREPARATION_MODIFIERS

template <class A>
struct CheckSubstitutionFlag
{
    static const bool value = (ModifierPriority<A>::value & SubstitutionMask) != 0;
};

template <class A, class B, 
          bool substitute=CheckSubstitutionFlag<A>::value>
struct SubstituteModifiers;

template <class A, class B>
struct SubstituteModifiers<A, B, false>
{
    typedef B type;
};

template <class A, template <class> class AA, class B, template <class> class BB>
struct SubstituteModifiers<AA<A>, BB<B>, true>
{
    typedef AA<typename SubstituteModifiers<A, B>::type> type;
};

template <class A, class B, template <class> class BB>
struct SubstituteModifiers<DefaultModifier<A>, BB<B>, true>
{
    typedef BB<typename SubstituteModifiers<A, B>::type> type;
};

template <class A, template <class> class AA, class B, template <class> class BB>
struct SubstituteModifiers<AA<A>, BB<B>, false>
{
    typedef BB<typename SubstituteModifiers<A, B>::type> type;
};

} // namespace acc_detail

template <class A, class B>
struct TransferModifiers
{
    typedef typename StandardizeTag<A>::type StdA;
    typedef typename StandardizeTag<B>::type StdB;
    typedef typename acc_detail::TagLongForm<StdA, acc_detail::MinPriority>::type AA;
    typedef typename acc_detail::TagLongForm<StdB, acc_detail::MinPriority>::type BB;
    typedef typename acc_detail::SubstituteModifiers<AA, BB>::type AB;
    typedef typename acc_detail::StandardizeTagLongForm<AB>::type StdAB;
    typedef typename StandardizeTag<StdAB>::type type;
};

template <class A, class HEAD, class TAIL>
struct TransferModifiers<A, TypeList<HEAD, TAIL> >
{
    typedef TypeList<typename TransferModifiers<A, HEAD>::type,
                     typename TransferModifiers<A, TAIL>::type> type;
};

template <class A>
struct TransferModifiers<A, void>
{
    typedef void type;
};

template <class TargetTag, class A=typename TargetTag::Dependencies>
struct StandardizeDependencies
#ifndef DOXYGEN
: public StandardizeDependencies<TargetTag, typename A::type>
#endif
{};

template <class TargetTag, class HEAD, class TAIL>
struct StandardizeDependencies<TargetTag, TypeList<HEAD, TAIL> >
{
    typedef typename StandardizeTag<TargetTag>::type Target;
    typedef typename TransferModifiers<Target, TypeList<HEAD, TAIL> >::type type;
};

template <class TargetTag>
struct StandardizeDependencies<TargetTag, void>
{
    typedef void type;
};

}} // namespace vigra::acc

#endif // VIGRA_ACCUMULATOR_GRAMMAR_HXX
