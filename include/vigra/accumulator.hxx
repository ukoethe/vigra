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
//#include "accumulator-detail.hxx"
#include "config.hxx"
#include "metaprogramming.hxx"
#include "bit_array.hxx"
#include "static_assert.hxx"
#include "mathutil.hxx"
#include "utilities.hxx"
#include "matrix.hxx"
#include "multi_math.hxx"
#include "eigensystem.hxx"
#include <algorithm>
#include <iostream>

#if 0
namespace vigra {

namespace acc1 {

struct AccumulatorBegin;
struct AccumulatorEnd;

    // a poor man's replacement for template typedefs
template <class A>
struct Synonym
{
    typedef A type;
};

template <class A, class B=typename Synonym<A>::type>
struct ExpandSynonym
{
    typedef typename ExpandSynonym<B>::type type;
};

template <class A>
struct ExpandSynonym<A, A>
{
    typedef A type;
};

template <class A>
struct AccumulatorTraits
{
    typedef A type;
    typedef void Contained;
    static const int priority = INT_MAX;
};

namespace detail  {

// FIXME: this should be refactored

template <int N1, int N2>
struct Less
{
    static const bool value = N1 < N2;
};

template <class T1, int Priority1, class T2, int Priority2>
struct Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined {};

template <class T1, int Priority, class T2>
struct Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined<T1, Priority, T2, Priority>;

template <class TAG, 
          int BOUND=INT_MIN, 
          bool SKIP=Less<AccumulatorTraits<TAG>::priority, BOUND>::value, 
          class Contained=typename ExpandSynonym<typename AccumulatorTraits<TAG>::Contained>::type>
struct LowestPriority
: public Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined<
                 typename LowestPriority<Contained, BOUND>::type, 
                 LowestPriority<Contained, BOUND>::priority,
                 TAG, 
                 AccumulatorTraits<TAG>::priority>
{
    typedef typename LowestPriority<Contained, BOUND>::type BestContained;
    static const int P1 = LowestPriority<Contained, BOUND>::priority;
    static const int P2 = AccumulatorTraits<TAG>::priority;
    static const int priority = P1 < P2 ? P1 : P2;    
    typedef typename IfBool<(P1 < P2), BestContained, TAG>::type type;
};

template <class TAG, int BOUND, class Contained>
struct LowestPriority<TAG, BOUND, true, Contained>
{
    typedef typename LowestPriority<Contained, BOUND>::type type;
    static const int priority = LowestPriority<Contained, BOUND>::priority;
};

template <class TAG, int BOUND, bool SKIP>
struct LowestPriority<TAG, BOUND, SKIP, void>
{
    typedef TAG type;
    static const int priority = AccumulatorTraits<TAG>::priority;
};

template <class TAG, int BOUND>
struct LowestPriority<TAG, BOUND, true, void>
{
    typedef TAG type;
    static const int priority = AccumulatorTraits<TAG>::priority;
};

template <class TAG, 
          int BOUND=INT_MIN, 
          int PRIORITY=LowestPriority<TAG, BOUND>::priority>
struct SortModifiers
{
    typedef AccumulatorTraits<typename LowestPriority<TAG, BOUND>::type> Traits;
    typedef typename Traits::template rebind<typename SortModifiers<TAG, PRIORITY+1>::type>::type RawType;
    typedef typename ExpandSynonym<RawType>::type type;
};

template <class TAG, int BOUND>
struct SortModifiers<TAG, BOUND, INT_MAX>
{
    typedef typename LowestPriority<TAG, BOUND>::type type;
};

template <class TAG, 
          int BOUND=INT_MIN, 
          bool SKIP=Less<AccumulatorTraits<TAG>::priority, BOUND>::value, 
          class Contained=typename AccumulatorTraits<TAG>::Contained>
struct LowestPriorityRaw
{
    typedef typename LowestPriorityRaw<Contained, BOUND>::type BestContained;
    static const int P1 = LowestPriorityRaw<Contained, BOUND>::priority;
    static const int P2 = AccumulatorTraits<TAG>::priority;
    static const int priority = P1 < P2 ? P1 : P2;    
    typedef typename IfBool<(P1 < P2), BestContained, TAG>::type type;
};

template <class TAG, int BOUND, class Contained>
struct LowestPriorityRaw<TAG, BOUND, true, Contained>
{
    typedef typename LowestPriorityRaw<Contained, BOUND>::type type;
    static const int priority = LowestPriorityRaw<Contained, BOUND>::priority;
};

template <class TAG, int BOUND, bool SKIP>
struct LowestPriorityRaw<TAG, BOUND, SKIP, void>
{
    typedef TAG type;
    static const int priority = AccumulatorTraits<TAG>::priority;
};

template <class TAG, 
          int BOUND=INT_MIN, 
          int PRIORITY=LowestPriorityRaw<TAG, BOUND>::priority>
struct SortRaw
{
    typedef AccumulatorTraits<typename LowestPriorityRaw<TAG, BOUND>::type> Traits;
    typedef typename Traits::template rebind<typename SortRaw<TAG, PRIORITY+1>::type>::type type;
};

template <class TAG, int BOUND>
struct SortRaw<TAG, BOUND, INT_MAX>
{
    typedef typename LowestPriorityRaw<TAG, BOUND>::type type;
};

} // namespace detail

template <class T>
struct StandardizeTag
{
        // since synonyms are only defined for the canonical order,
        // we must first sort without synonym expansion
    typedef typename detail::SortRaw<T>::type RawSorted;
    typedef typename ExpandSynonym<RawSorted>::type Expanded;
    typedef typename detail::SortModifiers<Expanded>::type type;
};

/****************************************************************************/
/*                                                                          */
/*                          accumulator modifiers                           */
/*                                                                          */
/****************************************************************************/

    // AccumulatorTraits::priority determines the modifier order.
    // Priorities have been choosen so that the canonical ascending order will be
    //     Coord< Normalize< Central< Plain > > >
    // This order is important for modifier transfer to work correctly.
    
    // modifier that forwards data after applying Centralize
template <class A>
class Central {};

    // last-seen-value accumulator that actually performs the centralization
class Centralize;

template <class A>
struct AccumulatorTraits<Central<A> >
{
    typedef A Contained;
    typedef Central<A> type;
    static const int priority = 20;
    
    template <class T>
    struct rebind
    {
        typedef Central<T> type;
    };
};

    // modifier that forwards data after applying PrincipalProjection
template <class A>
class Principal {};

    // last-seen-value accumulator that applies principal projection
class PrincipalProjection;

template <class A>
struct AccumulatorTraits<Principal<A> >
{
    typedef A Contained;
    typedef Principal<A> type;
    static const int priority = 20;
    
    template <class T>
    struct rebind
    {
        typedef Principal<T> type;
    };
};

    // Standard normalization method and normalization modifier
    // FIXME: figure out when caching is appropriate
struct DivideByCount;

template <class A, class NormlizationMethod=DivideByCount>
class Normalize;

template <class A, class NormlizationMethod>
struct AccumulatorTraits<Normalize<A, NormlizationMethod> >
{
    typedef A Contained;
    typedef Normalize<Contained, NormlizationMethod> type;
    static const int priority = 10;
    
    template <class T>
    struct rebind
    {
        typedef Normalize<T, NormlizationMethod> type;
    };
};

    // data access modifier that forwards the coordinate part of a compound data type
template <class A>
class Coord {};

template <class A>
struct AccumulatorTraits<Coord<A> >
{
    typedef A Contained;
    typedef Coord<Contained> type;
    static const int priority = 0;
    
    template <class T>
    struct rebind
    {
        typedef Coord<T> type;
    };
};

    // data access modifier that forwards data/weight pairs from a compound data type
template <class A>
class Weighted {};

template <class A>
struct AccumulatorTraits<Weighted<A> >
{
    typedef A Contained;
    typedef Weighted<Contained> type;
    static const int priority = 0;
    
    template <class T>
    struct rebind
    {
        typedef Weighted<T> type;
    };
};

    // data access modifier that forwards coordinate/weight pairs from a compound data type
template <class A>
class CoordWeighted {};

template <class A>
struct AccumulatorTraits<CoordWeighted<A> >
{
    typedef A Contained;
    typedef CoordWeighted<Contained>  type;
    static const int priority = 0;
    
    template <class T>
    struct rebind
    {
        typedef CoordWeighted<T> type;
    };
};

/****************************************************************************/
/*                                                                          */
/*                     accumulators to be implemented                       */
/*                     and their relationships                              */
/*                                                                          */
/****************************************************************************/

    // returns an identity matrix of appropriate size
class Axes;

    // sum over the power of values
template <unsigned N>
class PowerSum;

    // explicitly specialize PowerSum<0> to always use value_type=double, 
    // regardless of the data type
template <>
class PowerSum<0>;

typedef PowerSum<0> Count;

#define VIGRA_SIMPLE_SYNONYM(TAG, REPLACEMENT) \
template <> \
struct Synonym<TAG > \
{ \
    typedef REPLACEMENT type; \
};

    // Count ignores all modifers except for Weighted<>
VIGRA_SIMPLE_SYNONYM(Central<Count>, Count)
VIGRA_SIMPLE_SYNONYM(Principal<Count>, Count)
VIGRA_SIMPLE_SYNONYM(Coord<Count>, Count)
VIGRA_SIMPLE_SYNONYM(CoordWeighted<Count>, Weighted<Count>)

template <class NormlizationMethod>
struct Synonym<Normalize<Count, NormlizationMethod> >
{
    typedef Count type;
};

    // explicitly specialize PowerSum<1> to get rid of pow(value, 1)
template <>
class PowerSum<1>;
typedef PowerSum<1> Sum;

typedef PowerSum<2> SumOfSquares;

    // more normalization methods
struct RootDivideByCount;
struct DivideUnbiased;
struct RootDivideUnbiased;

typedef Normalize<Sum>                              Mean;
typedef Normalize<SumOfSquares, RootDivideByCount>  RootMeanSquares;

template <unsigned N>
class Moment;

template <unsigned N>
struct Synonym<Moment<N> >
{
    typedef Normalize<PowerSum<N> > type;
};

    // explicitly specialize Central<PowerSum<2> > to use incremental algorithms
    // (standard centralized accumulators use two passes)
template <>
class Central<PowerSum<2> >;

typedef Central<PowerSum<2> > SumOfSquaredDifferences;
typedef Central<PowerSum<2> > SSD;

class Variance;
class StdDev;
class UnbiasedVariance;
class UnbiasedStdDev;

VIGRA_SIMPLE_SYNONYM(Variance, Central<Moment<2> >)
VIGRA_SIMPLE_SYNONYM(Central<Variance>, Variance)

#define VIGRA_NORMALIZE_REPLACEMENT(INNER, METHOD) \
    Normalize<INNER, METHOD>

VIGRA_SIMPLE_SYNONYM(StdDev, VIGRA_NORMALIZE_REPLACEMENT(SSD, RootDivideByCount))
VIGRA_SIMPLE_SYNONYM(Central<StdDev>, VIGRA_NORMALIZE_REPLACEMENT(SSD, RootDivideByCount))

VIGRA_SIMPLE_SYNONYM(UnbiasedVariance, VIGRA_NORMALIZE_REPLACEMENT(SSD, DivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Central<UnbiasedVariance>, VIGRA_NORMALIZE_REPLACEMENT(SSD, DivideUnbiased))

VIGRA_SIMPLE_SYNONYM(UnbiasedStdDev, VIGRA_NORMALIZE_REPLACEMENT(SSD, RootDivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Central<UnbiasedStdDev>, VIGRA_NORMALIZE_REPLACEMENT(SSD, RootDivideUnbiased))

template <unsigned N>
class CentralMoment;

template <unsigned N>
struct Synonym<CentralMoment<N> >
{
    typedef Central<Moment<N> > type;
};

    // explicitly specialize Central<PowerSum<3> > and Central<PowerSum<4> > to implement 
    // merge functions (standard centralized accumulators do not support merge)
template <>
class Central<PowerSum<3> >;
template <>
class Central<PowerSum<4> >;

    // make sure that inherently centralized accumulators contain
    // the Central<> modifier
class SkewnessImpl;
typedef Central<SkewnessImpl> Skewness;

class KurtosisImpl;
typedef Central<KurtosisImpl> Kurtosis;

class FlatScatterMatrixImpl;
typedef Central<FlatScatterMatrixImpl> FlatScatterMatrix;

class CovarianceEigensystemImpl;
typedef Central<CovarianceEigensystemImpl> CovarianceEigensystem;

typedef Normalize<FlatScatterMatrix>         Covariance;
typedef Normalize<FlatScatterMatrix, DivideUnbiased> UnbiasedCovariance;

    // explicitly specialize Principal<PowerSum<2> >, Principal<Variance>, and 
    // Principal<Axes> to access  CovarianceEigensystem directly
template <> class Principal<PowerSum<2> >;
template <> class Normalize<Principal<PowerSum<2> > >;
template <> class Principal<Axes>;

VIGRA_SIMPLE_SYNONYM(Principal<Variance>, Principal<Moment<2> >)
VIGRA_SIMPLE_SYNONYM(Principal<StdDev>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, RootDivideByCount))
VIGRA_SIMPLE_SYNONYM(Principal<UnbiasedVariance>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, DivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Principal<UnbiasedStdDev>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, RootDivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Principal<Skewness>, Principal<SkewnessImpl>)
VIGRA_SIMPLE_SYNONYM(Principal<Kurtosis>, Principal<KurtosisImpl>)

    // sum of absolute values
class AbsSum;

    // explicitly specialize Central<AbsSum> and Principal<AbsSum> because 
    // they need two passes
template <>
class Central<AbsSum>;

typedef Central<AbsSum> SumOfAbsDifferences;

class MeanAbsoluteDeviation;
VIGRA_SIMPLE_SYNONYM(MeanAbsoluteDeviation, Normalize<SumOfAbsDifferences>)

template <>
class Principal<AbsSum>;

template <unsigned Percent>
class IncrementalQuantile;

typedef IncrementalQuantile<0>   Minimum;
typedef IncrementalQuantile<50>  IncrementalMedian;
typedef IncrementalQuantile<100> Maximum;

    // last-seen-value accumulator that maps [min, max] to another range (e.g. for histogram creation)
class RangeMapping;

    // works in 2-pass mode when RangeMapping is calculated from data
    // and in 1-pass mode when RangeMapping is given by the user.
template <unsigned BinCount>
class Histogram;

    // specify number of bins at runtime when BinCount=0
template <>
class Histogram<0>;

template <unsigned Percent, unsigned BinCount>
class HistogramQuantile;

template <unsigned NDim>
class MultiHistogram;

class AccumulatorArray;

#undef VIGRA_SIMPLE_SYNONYM
#undef VIGRA_NORMALIZE_REPLACEMENT

/****************************************************************************/
/*                                                                          */
/*                rules to transfer modifiers to dependencies               */
/*                                                                          */
/****************************************************************************/

    // basic rule: transfer all modifiers
template <class From, class To>
struct DontTransferModifier
{
    typedef VigraFalseType type;
};

    // chain begin rule: never transfer anything to AccumulatorBegin
template <class From>
struct DontTransferModifier<From, AccumulatorBegin>
{
    typedef VigraFalseType type;
};

    // normalization rule: never transfer normalization
template <class A, class NormlizationFunctor, class To>
struct DontTransferModifier<Normalize<A, NormlizationFunctor>, To>
{
    typedef VigraTrueType type;
};

    // specific rules
#define VIGRA_DONT_TRANSFER(FROM, TO)\
template <class A> \
struct DontTransferModifier<FROM<A>, TO> \
{\
    typedef VigraTrueType type; \
}; 

    // counting is always the same
VIGRA_DONT_TRANSFER(Coord, Count)
VIGRA_DONT_TRANSFER(Central, Count)
VIGRA_DONT_TRANSFER(Principal, Count)

    // Centralized sums are identical to zero
VIGRA_DONT_TRANSFER(Central, Sum)
VIGRA_DONT_TRANSFER(Principal, Sum)

    // modifying accumulators themselves must not be modified
VIGRA_DONT_TRANSFER(Central, Centralize)
VIGRA_DONT_TRANSFER(Central, PrincipalProjection)
VIGRA_DONT_TRANSFER(Central, RangeMapping)
VIGRA_DONT_TRANSFER(Principal, Centralize)
VIGRA_DONT_TRANSFER(Principal, PrincipalProjection)
VIGRA_DONT_TRANSFER(Principal, RangeMapping)

    // covariances are inherently central
VIGRA_DONT_TRANSFER(Principal, Covariance)
VIGRA_DONT_TRANSFER(Principal, UnbiasedCovariance)
VIGRA_DONT_TRANSFER(Principal, CovarianceEigensystem)

#undef VIGRA_DONT_TRANSFER

namespace detail {

// FIXME: this should be refactored

template <class T1, class T2, 
          class Next=typename AccumulatorTraits<T1>::Contained,
          class Contained=typename AccumulatorTraits<T2>::Contained>
struct ContainsModifier
{
    typedef typename IsSameType<typename AccumulatorTraits<T1>::template rebind<Contained>::type, T2>::type Same;
    typedef typename Or<Same, typename ContainsModifier<T1, Contained>::type>::type type;
};

template <class T1, class T2, class Next>
struct ContainsModifier<T1, T2, Next, void>
{
    typedef VigraFalseType type;
};

template <class T1, class T2, class Next>
struct ContainsModifier<T1, T2, void, Next>
{
    typedef VigraFalseType type;
};

template <class T1, class T2>
struct ContainsModifier<T1, T2, void, void>
{
    typedef VigraFalseType type;
};

template <class From, class To, class Next=typename AccumulatorTraits<To>::Contained>
struct ForbiddenTransfer
{
    typedef typename Or<typename DontTransferModifier<From, To>::type,
                        typename ForbiddenTransfer<From, Next>::type>::type type;
};

template <class From, class To>
struct ForbiddenTransfer<From, To, void>
{
    typedef typename DontTransferModifier<From, To>::type type;
};

template <class From, class To, class Next=typename AccumulatorTraits<From>::Contained,
         class Dont=typename Or<typename ContainsModifier<From,To>::type,
                                typename ForbiddenTransfer<From, To>::type>::type>
struct TransferModifiersImpl
{
    typedef typename TransferModifiersImpl<Next, To>::type Inner;
    typedef typename AccumulatorTraits<From>::template rebind<Inner>::type type;
};

template <class From, class To, class Next>
struct TransferModifiersImpl<From, To, Next, VigraTrueType>
{
    typedef typename TransferModifiersImpl<Next, To>::type Inner;
    typedef Inner type;
};

template <class From, class To, class Contains>
struct TransferModifiersImpl<From, To, void, Contains>
{
    typedef To type;
};

template<class A, class Next=typename AccumulatorTraits<A>::Contained >
struct RemoveCentralPrincipalPairs
{
    typedef typename RemoveCentralPrincipalPairs<Next>::type Cleaned;
    typedef typename AccumulatorTraits<A>::template rebind<Cleaned>::type type;
};

template<class A>
struct RemoveCentralPrincipalPairs<Principal<Central<A> >, Central<A> >
{
    typedef Principal<A> type;
};

template<class A>
struct RemoveCentralPrincipalPairs<Central<Principal<A> >, Principal<A> >
{
    typedef Central<A> type;
};

template<class A>
struct RemoveCentralPrincipalPairs<A, void >
{
    typedef A type;
};

} // namespace detail

template <class From, class To>
struct TransferModifiers
{
    typedef typename detail::TransferModifiersImpl<From, To>::type Unsorted;
    typedef typename detail::SortRaw<Unsorted>::type Sorted;
    typedef typename detail::RemoveCentralPrincipalPairs<Sorted>::type Cleaned;
    typedef typename StandardizeTag<Cleaned>::type type;
};

template <class From, class Head, class Tail>
struct TransferModifiers<From, TypeList<Head, Tail> >
{
    typedef TypeList<typename TransferModifiers<From, Head>::type,
                     typename TransferModifiers<From, Tail>::type> type;
};

template <class From>
struct TransferModifiers<From, void>
{
    typedef void type;
};

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

/****************************************************************************/
/*                                                                          */
/*                             infrastructure                               */
/*                                                                          */
/****************************************************************************/

struct AccumulatorBegin;

struct AccumulatorEnd 
{
    typedef AccumulatorEnd Tag;
    typedef void value_type;
    typedef bool result_type;
    static const unsigned int workInPass = 0; 
    static const int index = -1;
    
    bool operator()() const { return false; }
    
    template <unsigned, class T>
    void pass(T const &) {}
    template <unsigned, class T>
    void pass(T const &, double) {}
    
    template <class T>
    void merge(T const &) {}
    
    template <class T>
    void resize(T const &) {}
    
    void reset() {}
    
    void activate() {}
    bool isActive() const { return false; }
    
    unsigned int passesRequired() const
    {
        return 0;
    }
};

template <class T, class TAG, class NEXT=AccumulatorEnd>
struct AccumulatorBase;

namespace detail {

template <bool dynamic, unsigned LEVEL>
struct AccumulatorFlags
: public AccumulatorFlags<false, LEVEL>
{
    BitArray<LEVEL> active_accumulators_;
   
    template <int which>
    void activateImpl()
    {
        active_accumulators_.set<which>();
    }
    
    template <int which>
    bool isActiveImpl() const
    {
        return active_accumulators_.test<which>();
    }
    
    void reset()
    {
        active_accumulators_.clear();
    }
};

template <unsigned LEVEL>
struct AccumulatorFlags<false, LEVEL>
: public AccumulatorEnd
{
    typedef AccumulatorEnd BaseType;
    static const unsigned level = LEVEL;
    
    mutable BitArray<LEVEL> is_dirty_;
    
    template <int which>
    void setDirtyImpl() const
    {
        is_dirty_.set<which>();
    }
    
    template <int which>
    void setCleanImpl() const
    {
        is_dirty_.reset<which>();
    }
    
    template <int which>
    bool isDirtyImpl() const
    {
        return is_dirty_.test<which>();
    }
};

    // Insert the dependencies of the selected functors into the TypeList and sort
    // the list such that dependencies come after the functors using them. Make sure 
    // that each functor is contained only once.
template <class T>
struct AddDependencies
{
        // transform Selected<...> into TypeList<...>
    typedef typename AddDependencies<typename T::type>::type type;
};

template <class HEAD, class TAIL>
struct AddDependencies<TypeList<HEAD, TAIL> >
{
    typedef typename AddDependencies<TAIL>::type                                   TailWithDependencies;
    typedef typename AddDependencies<typename HEAD::Dependencies>::type            HeadDependencies;
    typedef TypeList<HEAD, HeadDependencies>                                       HeadWithDependencies;
    typedef typename PushUnique<HeadWithDependencies, TailWithDependencies>::type  type;
};

template <>
struct AddDependencies<void>
{
    typedef void type;
};

    // Helper class to activate dependencies at runtime (i.e. when activate<Tag>(accu) is called,
    // activate() must also be called for Tag's dependencies).
template <class Dependencies>
struct ActivateDependencies
{
    template <class Accumulator>
    static void exec(Accumulator & a)
    {
        activate<typename Dependencies::Head>(a);
        ActivateDependencies<typename Dependencies::Tail>::exec(a);
    }
};

template <>
struct ActivateDependencies<void>
{
    template <class Accumulator>
    static void exec(Accumulator & a)
    {}
};

template <class A, unsigned CurrentPass, bool Dynamic, unsigned WorkPass=A::workInPass>
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

    static void merge(A & a, A const & o)
    {
        a += o;
    }

    template <class Shape>
    static void resize(A & a, Shape const & o)
    {
        a.reshape(o);
    }
    
    static unsigned int passesRequired(A const & a)
    {
        return std::max(A::workInPass, a.next_.passesRequired());
    }
};

template <class A, unsigned CurrentPass, bool Dynamic>
struct DecoratorImpl<A, CurrentPass, Dynamic, CurrentPass>
{
    template <class T>
    static void exec(A & a, T const & t)
    {
        if(a.isActive())
            a.update(t);
    }

    template <class T>
    static void exec(A & a, T const & t, double weight)
    {
        if(a.isActive())
            a.update(t, weight);
    }

    static typename A::result_type get(A const & a)
    {
        vigra_precondition(a.isActive(),
            std::string("get(accumulator): attempt to access inactive statistic '")
                                 << typeid(A::Tag).name() << "'.");
        return a();
    }

    static void merge(A & a, A const & o)
    {
        if(a.isActive())
            a += o;
    }

    template <class Shape>
    static void resize(A & a, Shape const & o)
    {
        if(a.isActive())
            a.reshape(o);
    }
    
    static unsigned int passesRequired(A const & a)
    {
        return a.isActive()
                   ? std::max(A::workInPass, a.next_.passesRequired())
                   : a.next_.passesRequired();
    }
};

template <class T, class A, bool Dynamic, unsigned level>
struct Decorator
: public A
{
    typedef Decorator type;
    typedef Decorator & reference;
    typedef Decorator const & const_reference;
    
    template <class Shape>
    void resize(Shape const & s)
    {
        this->next_.resize(s);
        DecoratorImpl<Decorator, Decorator::workInPass, Dynamic>::resize(*this, s);
    }
    
    void reset()
    {
        this->next_.reset();
        A::reset();
    }
    
    typename A::result_type operator()() const
    {
        return DecoratorImpl<A, A::workInPass, Dynamic>::get(*this);
    }
    
    template <unsigned N>
    void pass(T const & t)
    {
        this->next_.pass<N>(t);
        DecoratorImpl<Decorator, N, Dynamic>::exec(*this, t);
    }
    
    template <unsigned N>
    void pass(T const & t, double weight)
    {
        this->next_.pass<N>(t, weight);
        DecoratorImpl<Decorator, N, Dynamic>::exec(*this, t, weight);
    }
    
    void merge(Decorator const & o)
    {
        DecoratorImpl<Decorator, Decorator::workInPass, Dynamic>::merge(*this, o);
        this->next_.merge(o.next_);
    }
    
    unsigned int passesRequired() const
    {
        return DecoratorImpl<Decorator, Decorator::workInPass, Dynamic>::passesRequired(*this);
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

template <unsigned int N, class T, class Alloc>
void reshapeImpl(MultiArray<N, T, Alloc> & a, typename MultiArrayShape<N>::type const & s, T const & initial = T())
{
    MultiArray<N, T, Alloc>(s, initial).swap(a);
}

template <class T, class Alloc>
void reshapeImpl(Matrix<T, Alloc> & a, Shape2 const & s, T const & initial = T())
{
    Matrix<T, Alloc>(s, initial).swap(a);
}

    // check if the accumulator chain A contains result_types that require runtime reshape
template <class A, class ResultType = typename A::value_type>
struct NeedsReshape
{
    typedef typename NeedsReshape<typename A::BaseType>::type type;
};

template <class A, unsigned int N, class T, class Alloc>
struct NeedsReshape<A, MultiArray<N, T, Alloc> >
{
    typedef VigraTrueType type;
};

template <class A, class T>
struct NeedsReshape<A, Matrix<T> >
{
    typedef VigraTrueType type;
};

template <>
struct NeedsReshape<AccumulatorEnd, AccumulatorEnd::value_type>
{
    typedef VigraFalseType type;
};

    // This functor is inserted on top of an accumulator chain to call reshape() 
    // when the first data item arrives. This is necessary if the shape of the result 
    // depends on the shape of the input and cannot be determined at compile time. 
    // The above NeedsReshape traits specify the types where this applies. If the chain
    // doesn't contain such types, ReshapeImpl will expand into a do-nothing version.
template <class NeedsReshape>
struct ReshapeImpl
{
    bool done_;
    
    ReshapeImpl()
    : done_(false)
    {}
    
    void reset()
    {
        done_ = false;
    }
    
    template <class A, class Shape>
    void operator()(A & a, Shape const & s)
    {
        if(!done_)
        {
            a.resize(s);
            done_ = true;
        }
    }
    
    template <class A, class T>
    void operator()(A & a, T const & t, MetaInt<1>)
    {
        if(!done_)
        {
            a.resize(shape(t));
            done_ = true;
        }
    }
    
    template <class A, class T, unsigned N>
    void operator()(A &, T const &, MetaInt<N>)
    {}
        
    template <unsigned int N, class U, class Stride>
    typename MultiArrayShape<N>::type
    shape(MultiArrayView<N, U, Stride> const & a)
    {
        return a.shape();
    }
    
    template <class U, int N>
    Shape1
    shape(TinyVector<U, N> const &)
    {
        return Shape1(N);
    }
};

template <>
struct ReshapeImpl<VigraFalseType>
{
    void reset()
    {}

    template <class A, class Shape>
    void operator()(A &, Shape const &)
    {}

    template <class A, class T, unsigned N>
    void operator()(A &, T const &, MetaInt<N>)
    {}
};

    // helper classes to create an accumulator chain from a TypeList
    // if dynamic=true,  a dynamic accumulator will be created
    // if dynamic=false, a plain accumulator will be created
template <class T, class Accumulators, bool dynamic=false, unsigned level = 0>
struct Compose
{
    typedef typename Accumulators::Head Tag; 
    typedef typename Compose<T, typename Accumulators::Tail, dynamic, level+1>::type BaseType;
    typedef Decorator<T, typename Tag::template Impl<T, AccumulatorBase<T, Tag, BaseType> >, dynamic, level>  type;
};

template <class T, bool dynamic, unsigned level>
struct Compose<T, void, dynamic, level> 
{ 
    typedef AccumulatorFlags<dynamic, level> type;
};

} // namespace detail 

/****************************************************************************/
/*                                                                          */
/*                            accumulator chain                             */
/*                                                                          */
/****************************************************************************/

    // Create an accumulator chain containing the Selected statistics and their dependencies.
template <class T, class Selected, bool dynamic = false>
struct Accumulator
: public detail::Compose<T, typename detail::AddDependencies<typename Selected::type>::type, dynamic>::type
{
    typedef typename detail::AddDependencies<typename Selected::type>::type AccumulatorTags;
    typedef typename detail::Compose<T, AccumulatorTags, dynamic>::type BaseType;
    
    typedef AccumulatorBegin                         Tag;
    typedef typename BaseType::argument_type         argument_type;
    typedef typename BaseType::first_argument_type   first_argument_type;
    typedef typename BaseType::second_argument_type  second_argument_type;
    typedef typename BaseType::result_type           result_type;

    BaseType next_;
    detail::ReshapeImpl<typename detail::NeedsReshape<BaseType>::type> reshape_;
    
    void reset()
    {
        reshape_.reset();
        next_.reset();
    }
    
    template <class Shape>
    void reshape(Shape const & s)
    {
        reshape_.reset();
        reshape_(next_, s);
    }

    template <unsigned N>
    void update(T const & t)
    {
        reshape_(next_, t, MetaInt<N>());
        next_.pass<N>(t);
    }
    
    template <unsigned N>
    void update(T const & t, double weight)
    {
        reshape_(*this, t, MetaInt<N>());
        next_.pass<N>(t, weight);
    }
    
    void operator+=(Accumulator const & o)
    {
        merge(o);
    }
    
    void merge(Accumulator const & o)
    {
        next_.merge(o.next_);
    }

    result_type operator()() const
    {
        return next_();
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
                     "Accumulator::updatePassN(): 0 < N < 6 required.");
        }
    }
    
	void updatePassN(T const & tt, double weight, unsigned int N)
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
                     "Accumulator::updatePassN(): 0 < N < 6 required.");
        }
    }
    
    unsigned int passesRequired() const
    {
        return next_.passesRequired();
    }
};

    // Create a dynamic accumulator chain containing the Selected statistics and their dependencies.
    // Statistics will only be computed if activate<Tag>() is called at runtime.
template <class T, class Selected>
struct DynamicAccumulator
: public Accumulator<T, Selected, true>
{};

/****************************************************************************/
/*                                                                          */
/*                        generic access functions                          */
/*                                                                          */
/****************************************************************************/

namespace detail {

template <class TAG, class A, class FromTag=typename A::Tag>
struct LookupTagImpl
: public LookupTagImpl<TAG, typename A::BaseType>
{};

template <class TAG, class A, class FromTag>
struct LookupTagImpl<TAG, A const, FromTag>
: public LookupTagImpl<TAG, A>
{
    typedef typename LookupTagImpl<TAG, A>::type const & reference;
};

template <class TAG, class A>
struct LookupTagImpl<TAG, A, TAG>
{
    typedef TAG Tag;
    typedef A type;
    typedef A & reference;
    typedef typename A::value_type value_type;
    typedef typename A::result_type result_type;
};

template <class TAG, class A>
struct LookupTagImpl<TAG, A, AccumulatorEnd>
{
    typedef TAG Tag;
    typedef A type;
    typedef A & reference;
    typedef typename A::value_type value_type;
    typedef typename A::result_type result_type;
};

} // namespace detail

template <class Tag, class A>
struct LookupTag
: public detail::LookupTagImpl<typename TransferModifiers<typename A::Tag, 
                                                         typename StandardizeTag<Tag>::type>::type, A>
{}; 

namespace detail {

template <class Tag, class FromTag, class reference>
struct CastImpl
{
    template <class A>
    static reference exec(A & a)
    {
        return CastImpl<Tag, typename A::BaseType::Tag, reference>::exec(a.next_);
    }
};

template <class Tag, class reference>
struct CastImpl<Tag, Tag, reference>
{
    template <class A>
    static reference exec(A & a)
    {
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
};

template <class Tag, class FromTag, class result_type>
struct GetImpl
{
    template <class A>
    static result_type exec(A const & a)
    {
        return GetImpl<Tag, typename A::BaseType::Tag, result_type>::exec(a.next_);
    }
};

template <class Tag, class result_type>
struct GetImpl<Tag, Tag, result_type>
{
    template <class A>
    static result_type exec(A const & a)
    {
        return a();
    }
};

template <class Tag, class result_type>
struct GetImpl<Tag, AccumulatorEnd, result_type>
{
    template <class A>
    static bool exec(A const & a)
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
        return false;
    }
};

} // namespace detail

    // cast an accumulator chain to the type specified by Tag
template <class TAG, class A>
typename LookupTag<TAG, A>::reference
cast(A & a)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::reference reference;
    return detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a);
}

    // get the result of the accumulator specified by Tag
template <class TAG, class A>
typename LookupTag<TAG, A>::result_type
get(A const & a)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::result_type result_type;
    return detail::GetImpl<StandardizedTag, typename A::Tag, result_type>::exec(a);
}

    // activate the dynamic accumulator specified by Tag
template <class Tag, class A>
void
activate(A & a)
{
    cast<Tag>(a).activate();
}

    // activate the dynamic accumulator specified by Tag
template <class Tag, class A>
bool
isActive(A const & a)
{
    return cast<Tag>(a).isActive();
}

/****************************************************************************/
/*                                                                          */
/*                        the actual accumulators                           */
/*                                                                          */
/****************************************************************************/

    // this class closes an accumulator chain
template <class T, class TAG, class NEXT>
struct AccumulatorBase 
{
	typedef TAG            Tag;
        // FIXME: this will become more sophisticated to support modifiers
    typedef typename TransferModifiers<TAG, typename Tag::Dependencies::type>::type Dependencies;
	typedef NEXT           BaseType;
    typedef T const &      argument_type;
    typedef argument_type  first_argument_type;
    typedef double         second_argument_type;
    
    static const unsigned int workInPass = 1;
    static const int index = NEXT::index + 1;
    
    NEXT next_;
    
    void activate()
    {
        next_.activateImpl<index>();
        detail::ActivateDependencies<Dependencies>::exec(*this);
    }
    
    template <int INDEX>
    void activateImpl()
    {
        next_.activateImpl<INDEX>();
    }
    
    bool isActive() const
    {
        return next_.isActiveImpl<index>();
    }
    
    template <int INDEX>
    bool isActiveImpl() const
    {
        return next_.isActiveImpl<INDEX>();
    }
    
    void setDirty() const
    {
        next_.setDirtyImpl<index>();
    }
    
    template <int INDEX>
    void setDirtyImpl() const
    {
        next_.setDirtyImpl<INDEX>();
    }
    
    void setClean() const
    {
        next_.setCleanImpl<index>();
    }
    
    template <int INDEX>
    void setCleanImpl() const
    {
        next_.setCleanImpl<INDEX>();
    }
    
    bool isDirty() const
    {
        return next_.isDirtyImpl<index>();
    }
    
    template <int INDEX>
    bool isDirtyImpl() const
    {
        return next_.isDirtyImpl<INDEX>();
    }
    
    void reset()
    {}
    
    template <class Shape>
    void reshape(Shape const &)
    {}
    
	void operator+=(AccumulatorBase const &)
    {}
    
	void update(argument_type)
    {}
    
    void update(first_argument_type, second_argument_type)
    {}
};

template <class T>
struct AccumulatorResultTraits
{
    typedef T                                       type;
    typedef T                                       element_type;
    typedef typename NumericTraits<T>::RealPromote  element_promote_type;
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
    typedef typename NumericTraits<T>::RealPromote       element_promote_type;
    typedef TinyVector<T, N>                             MinmaxType;
    typedef TinyVector<element_promote_type, N>          SumType;
    typedef TinyVector<element_promote_type, N*(N+1)/2>  FlatCovarianceType;
    typedef Matrix<element_promote_type>                 CovarianceType;
};

template <unsigned int N, class T, class Stride>
struct AccumulatorResultTraits<MultiArrayView<N, T, Stride> >
{
    typedef MultiArrayView<N, T, Stride>            type;
    typedef T                                       element_type;
    typedef typename NumericTraits<T>::RealPromote  element_promote_type;
    typedef MultiArray<N, T>                        MinmaxType;
    typedef MultiArray<N, element_promote_type>     SumType;
    typedef MultiArray<1, element_promote_type>     FlatCovarianceType;
    typedef Matrix<element_promote_type>            CovarianceType;
};

/*
important notes on modifiers:
 * upon accumulator creation, reorder modifiers so that data access is innermost, 
   and normalization is outermost, e.g.:
        Normalize<Principal<Coord<Skewness> > >
   this mcould be done using priority traits
 * automatically transfer modifiers to dependencies as appropriate
 * automatically adjust modifiers for lookup (cast and get) of dependent accumulators
 * modifiers must adjust workInPass for the contained accumulator as appropriate
 * we may implement convenience versions of Select that apply a modifier to all 
   contained tags at once
 * make sure that weighted accumulators have their own Count object when used together
   with unweighted ones (or disallow combination of weighted and unweighted processing)
 * certain accumulators must remain unchanged when wrapped in certain modifiers: 
    * Count: always
    * Centralize, PricipleProjection sometimes
 * will it be useful to implement initPass<N>() or finalizePass<N>() ?
*/

template <class T, class BASE, 
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
        detail::reshapeImpl(value_, s);
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
    
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE, double, double>
    {
        void update(T const & t)
        {
            ++value_;
        }
        
        void update(T const & t, double weight)
        {
            value_ += weight;
        }
    };
};

// Sum
template <>
class PowerSum<1>
{
  public:
    typedef Select<> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        void update(T const & t)
        {
            value_ += t;
        }
        
        void update(T const & t, double weight)
        {
            value_ += weight*t;
        }
    };
};

template <unsigned N>
class PowerSum<N>
{
  public:
    typedef Select<> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += pow(t, (int)N);
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*pow(t, (int)N);
        }
    };
};

class AbsSum
{
  public:
    typedef Select<> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += abs(t);
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*abs(t);
        }
    };
};

template <>
class IncrementalQuantile<0>
{
  public:
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_type element_type;
        typedef typename AccumulatorResultTraits<T>::MinmaxType   value_type;
        typedef value_type const &                          result_type;

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
            detail::reshapeImpl(value_, s, NumericTraits<element_type>::max());
        }
        
       void operator+=(Impl const & o)
        {
            using namespace multi_math;
            value_ = min(value_, o.value_);
        }
    
        void update(T const & t)
        {
            using namespace multi_math;
            value_ = min(value_, t);
        }
        
        void update(T const & t, double weight)
        {
            vigra_precondition(false, "Minimum accumulator does not support weights.");
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

template <>
class IncrementalQuantile<100>
{
  public:
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_type element_type;
        typedef typename AccumulatorResultTraits<T>::MinmaxType   value_type;
        typedef value_type const &                          result_type;

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
            detail::reshapeImpl(value_, s, NumericTraits<element_type>::min());
        }
        
        void operator+=(Impl const & o)
        {
            using namespace multi_math;
            value_ = max(value_, o.value_);
        }
    
        void update(T const & t)
        {
            using namespace multi_math;
            value_ = max(value_, t);
        }
        
        void update(T const & t, double weight)
        {
            vigra_precondition(false, "Maximum accumulator does not support weights.");
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

template <class T, class BASE, 
          class VALUE_TYPE = typename AccumulatorResultTraits<T>::SumType>
struct CachedResultBase
: public BASE
{
    typedef typename AccumulatorResultTraits<T>::element_type  element_type;
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
        detail::reshapeImpl(value_, s);
    }

    void operator+=(CachedResultBase const &)
    {
        this->setDirty();
    }

    void update(T const &)
    {
        this->setDirty();
    }
    
    void update(T const &, double)
    {
         this->setDirty();
    }
};

// cached Mean and Variance
template <class TAG>
class Normalize<TAG, DivideByCount>
{
  public:
    typedef Select<TAG, Count> Dependencies;
  
    template <class T, class BASE>
    struct Impl
    : public CachedResultBase<T, BASE> 
    {
        result_type operator()() const
        {
            if(this->isDirty())
            {
                using namespace multi_math;
                value_ = get<TAG>(*this) / get<Count>(*this);
                this->setClean();
            }
            return value_;
        }
    };
};

// UnbiasedVariance
template <class TAG>
class Normalize<TAG, DivideUnbiased>
{
  public:
    typedef Select<TAG, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<TAG, BASE>::value_type  value_type;
        typedef value_type                                 result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return get<TAG>(*this) / (get<Count>(*this) - 1.0);
        }
    };
};

// RootMeanSquares and StdDev
template <class TAG>
class Normalize<TAG, RootDivideByCount>
{
  public:
    typedef Normalize<TAG, DivideByCount> TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<TargetTag, BASE>::value_type  value_type;
        typedef value_type                                       result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return sqrt(get<TargetTag>(*this));
        }
    };
};

// UnbiasedStdDev
template <class TAG>
class Normalize<TAG, RootDivideUnbiased>
{
  public:
    typedef Normalize<TAG, DivideUnbiased> TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<TargetTag, BASE>::value_type  value_type;
        typedef value_type                                       result_type;
        
        result_type operator()() const
        {
            using namespace multi_math;
            return sqrt(get<TargetTag>(*this));
        }
    };
};

class Centralize
{
  public:
    typedef Select<Mean> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename AccumulatorResultTraits<T>::element_promote_type element_type;
        typedef typename AccumulatorResultTraits<T>::SumType              value_type;
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
            detail::reshapeImpl(value_, s);
        }
        
        void update(T const & t)
        {
            using namespace vigra::multi_math;
            value_ = t - get<Mean>(*this);
        }
        
        void update(T const & t, double)
        {
            using namespace vigra::multi_math;
            value_ = t - get<Mean>(*this);
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

template <>
class Central<PowerSum<2> >
{
  public:
    typedef Select<Mean, Count> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        void operator+=(Impl const & o)
        {
            using namespace vigra::multi_math;
            double n1 = get<Count>(*this), n2 = get<Count>(o);
            if(n1 == 0.0)
            {
                this->value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                this->value_ += o.value_ + n1 * n2 / (n1 + n2) * sq(get<Mean>(*this) - get<Mean>(o));
            }
        }
    
        void update(T const & t)
        {
            double n = get<Count>(*this);
            if(n > 1.0)
            {
                using namespace vigra::multi_math;
                value_ += n / (n - 1.0) * sq(get<Mean>(*this) - t);
            }
        }
        
        void update(T const & t, double weight)
        {
            double n = get<Count>(*this);
            if(n > weight)
            {
                using namespace vigra::multi_math;
                value_ += n / (n - weight) * sq(get<Mean>(*this) - t);
            }
        }
    };
};

    // FIXME: automatic forwarding of dependencies to Central<Moment<3> > doesn't work
template <>
class Central<PowerSum<3> >
{
  public:
    typedef Select<Centralize, Count, Mean, Central<PowerSum<2> > > Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            typedef Central<PowerSum<2> > Sum2Tag;
            
            using namespace vigra::multi_math;
            double n1 = get<Count>(*this), n2 = get<Count>(o);
            if(n1 == 0.0)
            {
                this->value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                double n = n1 + n2;
                double weight = n1 * n2 * (n1 - n2) / sq(n);
                value_type delta = get<Mean>(o) - get<Mean>(*this);
                this->value_ += o.value_ + weight * pow(delta, 3) +
                               3.0 / n * delta * (n1 * get<Sum2Tag>(o) - n2 * get<Sum2Tag>(*this));
            }
        }
    
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += pow(get<Centralize>(*this), 3);
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*pow(get<Centralize>(*this), 3);
        }
    };
};

template <>
class Central<PowerSum<4> >
{
  public:
    typedef Select<Centralize, Central<PowerSum<3> > > Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            typedef Central<PowerSum<2> > Sum2Tag;
            typedef Central<PowerSum<3> > Sum3Tag;

            using namespace vigra::multi_math;
            double n1 = get<Count>(*this), n2 = get<Count>(o);
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
                value_type delta = get<Mean>(o) - get<Mean>(*this);
                this->value_ += o.value_ + weight * pow(delta, 4) +
                              6.0 / n_2 * sq(delta) * (n1_2 * get<Sum2Tag>(o) + n2_2 * get<Sum2Tag>(*this)) +
                              4.0 / n * delta * (n1 * get<Sum3Tag>(o) - n2 * get<Sum3Tag>(*this));
            }
        }
    
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += pow(get<Centralize>(*this), 4);
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*pow(get<Centralize>(*this), 4);
        }
    };
};

template <unsigned int N>
class Central<PowerSum<N> >
{
  public:
    typedef Select<Centralize> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "Central<PowerSum<N> >::operator+=(): not implemented for N > 4.");
        }
    
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += pow(get<Centralize>(*this), (int)N);
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*pow(get<Centralize>(*this), (int)N);
        }
    };
};

template <>
class Central<AbsSum>
{
  public:
    typedef Select<Centralize> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "Central<AbsSum>::operator+=(): not supported.");
        }
    
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += abs(get<Centralize>(*this));
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*abs(get<Centralize>(*this));
        }
    };
};

template <>
class Central<SkewnessImpl>
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<3> > > Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupTag<Central<PowerSum<3> >, BASE>::value_type   value_type;
        typedef value_type                                                    result_type;

        result_type operator()() const
        {
            typedef Central<PowerSum<3> > Sum3;
            typedef Central<PowerSum<2> > Sum2;
        
			using namespace multi_math;
            return sqrt(get<Count>(*this)) * get<Sum3>(*this) / pow(get<Sum2>(*this), 1.5);
        }
    };
};

template <>
class Central<KurtosisImpl>
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<4> > > Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupTag<Central<PowerSum<4> >, BASE>::value_type value_type;
        typedef value_type                                                  result_type;

        result_type operator()() const
        {
            typedef Central<PowerSum<4> > Sum4;
            typedef Central<PowerSum<2> > Sum2;
        
			using namespace multi_math;
            return get<Count>(*this) * get<Sum4>(*this) / sq(get<Sum2>(*this));
        }
    };
};

namespace detail {

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

} // namespace detail

// we only store the flattened upper triangular part of the scatter matrix
template <>
class Central<FlatScatterMatrixImpl>
{
  public:
    typedef Select<Mean, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<T>::FlatCovarianceType    value_type;
        typedef value_type const &                                   result_type;
       
        typedef typename AccumulatorResultTraits<T>::SumType        SumType;

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
            detail::reshapeImpl(value_, Shape1(size*(size+1)/2));
            detail::reshapeImpl(diff_, s);
        }
        
        void operator+=(Impl const & o)
        {
            double n1 = get<Count>(*this), n2 = get<Count>(o);
            if(n1 == 0.0)
            {
                value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                using namespace vigra::multi_math;
                diff_ = get<Mean>(*this) - get<Mean>(o);
                detail::updateFlatScatterMatrix(value_, diff_, n1 * n2 / (n1 + n2));
                value_ += o.value_;
            }
        }
    
        void update(T const & t)
        {
            compute(t);
        }
        
        void update(T const & t, double weight)
        {
            compute(t, weight);
        }
        
        result_type operator()() const
        {
            return value_;
        }
        
      private:
        void compute(T const & t, double weight = 1.0)
        {
            double n = get<Count>(*this);
            if(n > weight)
            {
                using namespace vigra::multi_math;
                diff_ = get<Mean>(*this) - t;
                detail::updateFlatScatterMatrix(value_, diff_, n * weight / (n - weight));
            }
        }
    };
};

// Covariance
template <>
class Normalize<FlatScatterMatrix, DivideByCount>
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public CachedResultBase<T, BASE,
                              typename AccumulatorResultTraits<T>::CovarianceType>
    {
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            detail::reshapeImpl(value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                detail::flatScatterMatrixToCovariance(value_, get<FlatScatterMatrix>(*this), get<Count>(*this));
                this->setClean();
            }
            return value_;
        }
    };
};

// UnbiasedCovariance
template <>
class Normalize<FlatScatterMatrix, DivideUnbiased>
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public CachedResultBase<T, BASE,
                              typename AccumulatorResultTraits<T>::CovarianceType>
    {
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            detail::reshapeImpl(value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                detail::flatScatterMatrixToCovariance(value_, get<FlatScatterMatrix>(*this), get<Count>(*this) - 1.0);
                this->setClean();
            }
            return value_;
        }
    };
};

// CovarianceEigensystem
template <>
class Central<CovarianceEigensystemImpl>
{
  public:
    typedef Select<Covariance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<T>::SumType               EigenvalueType;
        typedef typename AccumulatorResultTraits<T>::CovarianceType        EigenvectorType;
        typedef std::pair<EigenvalueType, EigenvectorType>                 value_type;
        typedef value_type const &                                         result_type;

        mutable value_type value_;
        
        Impl()
        : value_()
        {}
        
        void operator+=(Impl const &)
        {
            this->setDirty();
        }

        void update(T const &)
        {
            this->setDirty();
        }
        
        void update(T const &, double)
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
            detail::reshapeImpl(value_.first, Shape2(size,1));
            detail::reshapeImpl(value_.second, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                compute(get<Covariance>(*this), value_.first, value_.second);
                this->setClean();
            }
            return value_;
        }
        
      private:
        template <class Cov, class EW, class EV>
        static void compute(Cov const & cov, EW & ew, EV & ev)
        {
            // create a view because EW could be a TinyVector
            MultiArrayView<2, element_type> ewview(Shape2(cov.shape(0), 1), &ew[0]);
            symmetricEigensystem(cov, ewview, ev);
        }
        
        static void compute(double cov, double & ew, double & ev)
        {
            ew = cov;
            ev = 1.0;
        }
    };
};

// covariance eigenvalues
template <>
class Principal<PowerSum<2> >
{
  public:
    typedef Select<CovarianceEigensystem, Count> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<CovarianceEigensystem, BASE>::type::EigenvalueType value_type;
        typedef value_type result_type;
        
        result_type operator()() const
        {
            using namespace vigra::multi_math;            
            return get<Count>(*this)*get<CovarianceEigensystem>(*this).first;
        }
    };
};

// Principal<Variance> == covariance eigenvalues
template <>
class Normalize<Principal<PowerSum<2> > >
{
  public:
    typedef Select<CovarianceEigensystem> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<CovarianceEigensystem, BASE>::type::EigenvalueType value_type;
        typedef value_type const & result_type;
        
        result_type operator()() const
        {
            using namespace vigra::multi_math;            
            return get<CovarianceEigensystem>(*this).first;
        }
    };
};

template <>
class Principal<AbsSum>
{
  public:
    typedef Select<PrincipalProjection> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public SumBaseImpl<T, BASE>
    {
        static const unsigned int workInPass = 2;
        
        void operator+=(Impl const & o)
        {
            vigra_precondition(false,
                "Principal<AbsSum>::operator+=(): not supported.");
        }
    
        void update(T const & t)
        {
            using namespace vigra::multi_math;            
            value_ += abs(get<PrincipalProjection>(*this));
        }
        
        void update(T const & t, double weight)
        {
            using namespace vigra::multi_math;            
            value_ += weight*abs(get<PrincipalProjection>(*this));
        }
    };
};

}} // namespace vigra::acc1

#endif

#endif // VIGRA_ACCUMULATOR_HXX
