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

namespace acc1 {

/**************************************************************************/
/*                                                                        */
/*                      irreducible basic accumulators                    */
/*                                                                        */
/**************************************************************************/
    
   
class Axes;                                    // returns an identity matrix of appropriate size

template <unsigned N> class PowerSum;          // sum over powers of values
template <unsigned N> class AbsPowerSum;       // sum over powers of absolute values

class SkewnessImpl;                            // skewness
class KurtosisImpl;                            // kurtosis
class FlatScatterMatrixImpl;                   // flattened upper-triangular part of the scatter matrix
class CovarianceEigensystemImpl;               // eigenvalues and eigenvectors of the covariance matrix

template <unsigned Percent> class  Quantile;   // quantiles (including minimum and maximum)
template <unsigned BinCount> class Histogram;  // histogram with fixed number of bins
template <>  class Histogram<0>;               // histogram where number of bins is specified at runtime
template <unsigned Percent, class Hist> 
class HistogramQuantile;                       // quantiles computed from Histogram
template <unsigned NDim> class MultiHistogram; // multi-dimensional histogram
                                               // (always specify number of bins at runtime)

class AccumulatorArray;                        // one accumulator for each region (may be better implemented by other means)

class CachePreparedData;                       // cache the prepared (centralized, whitened etc.) version of the last data value 
class RangeMapping;                            // map value from [min, max] to another range and cache result (e.g. for histogram creation)


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
template <class A> class Coord;          // use pixel coordinate instead of pixel value
template <class A> class Weighted;       // use (value, weight) pairs
template <class A> class CoordWeighted;  // use (coord, weight) pairs

    // data preparation
template <class A> class Central;    // subtract mean
template <class A> class Principal;  // subtract mean and rotate to principle coordinates
template <class A> class Whitened;   // transform to principal coordinates and scale to unit variance

/**************************************************************************/
/*                                                                        */
/*                              synonyms                                  */
/*                                                                        */
/**************************************************************************/

typedef PowerSum<0>                                 Count;
typedef PowerSum<1>                                 Sum;
typedef PowerSum<2>                                 SumOfSquares;

typedef DivideByCount<Sum>                          Mean;
typedef RootDivideByCount<SumOfSquares>             RootMeanSquares;

// desired pseudocode (unfortunately not legal in C++)
//
//     template <unsigned N> 
//     typedef DivideByCount<PowerSum<N> >          Moment;
//
// actual definition (desired behavior is realised by rules below)
//
template <unsigned N> class                         Moment;  
template <unsigned N> class                         CentralMoment;  

typedef Central<PowerSum<2> >                       SumOfSquaredDifferences;
typedef SumOfSquaredDifferences                     SSD;

typedef DivideByCount<Central<PowerSum<2> > >       Variance;
typedef RootDivideByCount<Central<PowerSum<2> > >   StdDev;
typedef DivideUnbiased<Central<PowerSum<2> > >      UnbiasedVariance;
typedef RootDivideUnbiased<Central<PowerSum<2> > >  UnbiasedStdDev;

typedef Central<SkewnessImpl>                       Skewness;
typedef Central<KurtosisImpl>                       Kurtosis;

typedef Central<FlatScatterMatrixImpl>              FlatScatterMatrix;
typedef DivideByCount<FlatScatterMatrix>            Covariance;
typedef DivideUnbiased<FlatScatterMatrix>           UnbiasedCovariance;
typedef Central<CovarianceEigensystemImpl>          CovarianceEigensystem;

typedef AbsPowerSum<1>                              AbsSum;
typedef Central<AbsSum>                             SumOfAbsDifferences;
typedef DivideByCount<SumOfAbsDifferences>          MeanAbsoluteDeviation;

typedef Quantile<0>                                 Minimum;
typedef Quantile<100>                               Maximum;

typedef Central<CachePreparedData>                  Centralize;
typedef Principal<CachePreparedData>                PrincipalProjection;
typedef Whitened<CachePreparedData>                 Whiten;

/**************************************************************************/
/*                                                                        */
/*                        Tag standardization rules                       */
/*                                                                        */
/**************************************************************************/

namespace detail {

template <class A>
struct ModifierRule
{
    typedef A type;
};

} // namespace detail

template <class A>
struct Error___Tag_modifiers_of_same_kind_must_not_be_combined;

    // apply rules as long as the Tag type changes ...
template <class A, class S=typename detail::ModifierRule<A>::type>
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

    // fail when the tag spec was non-conforming 
template <class A, class B>
struct StandardizeTag<A, Error___Tag_modifiers_of_same_kind_must_not_be_combined<B> >
    : public Error___Tag_modifiers_of_same_kind_must_not_be_combined<B>
{};


    // Select is a synonym for MakeTypeList
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
    typename StandardizeTag<T19>::type, typename StandardizeTag<T20>::type >
{};

namespace detail {

enum { MinPriority = 1, 
       AccumulatorPriority = 8,
       PrepareDataPriority = 4,
       NormalizePriority = 2,
       AccessDataPriority = 1,
       MaxPriority = 8,
       SubstitutionMask = PrepareDataPriority | AccessDataPriority };

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

VIGRA_MODIFIER_PRIORITY(Coord, AccessDataPriority)
VIGRA_MODIFIER_PRIORITY(Weighted, AccessDataPriority)
VIGRA_MODIFIER_PRIORITY(CoordWeighted, AccessDataPriority)

VIGRA_MODIFIER_PRIORITY(DivideByCount, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(RootDivideByCount, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(DivideUnbiased, NormalizePriority)
VIGRA_MODIFIER_PRIORITY(RootDivideUnbiased, NormalizePriority)

VIGRA_MODIFIER_PRIORITY(Central, PrepareDataPriority)
VIGRA_MODIFIER_PRIORITY(Principal, PrepareDataPriority)
VIGRA_MODIFIER_PRIORITY(Whitened, PrepareDataPriority)

#undef VIGRA_MODIFIER_PRIORITY

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

    // modifers are sorted by ascending priority
template <class A, int compare=ModifierCompareToInner<A>::value>
struct ModifierOrder;

    // do nothing of the order is correct (compare == -1)
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

#define VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(OUTER, INNER) \
template <class A> \
struct ModifierOrder<OUTER<INNER<A > >, 0> \
{ \
    typedef OUTER<A > type; \
};

    // the strongest data preparation modifier takes precendence
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Principal, Central)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Central)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Principal)

    // drop one modifier when duplication occurs
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Central, Central)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Principal, Principal)
VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS(Whitened, Whitened)

#undef VIGRA_CLEANUP_DATA_PREPARATION_MODIFIERS

    // recurse down the modifier chain
template <class A, template <class> class B>
struct ModifierRule<B<A> >
{
    typedef typename ModifierOrder<B<typename StandardizeTag<A>::type> >::type type;
};

    // counting modified data is the same as counting data ...
template <template <class> class A>
struct ModifierRule<A<Count> >
{
    typedef Count type;
};

    // ... except for weighted data
template <>
struct ModifierRule<Weighted<Count> >
{
    typedef Weighted<Count> type;
};

template <>
struct ModifierRule<CoordWeighted<Count> >
{
    typedef Weighted<Count> type;
};

    // expand the Moment<N> synonym
template <unsigned N>
struct ModifierRule<Moment<N> >
{
    typedef DivideByCount<PowerSum<N> > type;
};

    // expand the CentralMoment<N> synonym
template <unsigned N>
struct ModifierRule<CentralMoment<N> >
{
    typedef DivideByCount<Central<PowerSum<N> > > type;
};

    // reduce even absolute powers to plain powers
template <unsigned N>
struct ModifierRule<AbsPowerSum<N> >
{
    typedef typename IfBool<(N % 2 == 0), PowerSum<N>, AbsPowerSum<N> >::type type;
};

    // reductions to CoordWeighted<A>
template <class A>
struct ModifierRule<Weighted<Coord<A> > >
{
    typedef CoordWeighted<A> type;
};

template <class A>
struct ModifierRule<Coord<Weighted<A> > >
{
    typedef CoordWeighted<A> type;
};

} // namespace detail

/**************************************************************************/
/*                                                                        */
/*                           Tag transfer rules                           */
/*                                                                        */
/**************************************************************************/

namespace detail {

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

template <>
struct LongModifierRule<Central<Sum> >
{
    typedef Sum type;
};

template <>
struct LongModifierRule<Principal<Sum> >
{
    typedef Sum type;
};

template <>
struct LongModifierRule<Whitened<Sum> >
{
    typedef Sum type;
};

template <>
struct LongModifierRule<Principal<CovarianceEigensystemImpl> >
{
    typedef CovarianceEigensystem type;
};

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

template <class A0, template <class> class A1, class B0, template <class> class B1>
struct SubstituteModifiers<A1<A0>, B1<B0>, true>
{
    typedef A1<typename SubstituteModifiers<A0, B0>::type> type;
};

template <class A0, template <class> class A1, class B0, template <class> class B1>
struct SubstituteModifiers<A1<A0>, B1<B0>, false>
{
    typedef B1<typename SubstituteModifiers<A0, B0>::type> type;
};

} // namespace detail

template <class A, class B>
struct TransferModifiers
{
    typedef typename StandardizeTag<A>::type StdA;
    typedef typename StandardizeTag<B>::type StdB;
    typedef typename detail::TagLongForm<StdA, detail::MinPriority>::type AA;
    typedef typename detail::TagLongForm<StdB, detail::MinPriority>::type BB;
    typedef typename detail::SubstituteModifiers<AA, BB>::type AB;
    typedef typename detail::StandardizeTagLongForm<AB>::type StdAB;
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

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_GRAMMAR_HXX
