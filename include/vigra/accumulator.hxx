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

namespace vigra {

namespace acc1 {

#if 1
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

template <class A, bool CacheResult=true, class NormlizationMethod=DivideByCount>
class Normalize {};

template <class A, bool CacheResult, class NormlizationMethod>
struct AccumulatorTraits<Normalize<A, CacheResult, NormlizationMethod> >
{
    typedef A Contained;
    typedef Normalize<Contained, CacheResult, NormlizationMethod> type;
    static const int priority = 10;
    
    template <class T>
    struct rebind
    {
        typedef Normalize<T, CacheResult, NormlizationMethod> type;
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

template <bool CacheResult, class NormlizationMethod>
struct Synonym<Normalize<Count, CacheResult, NormlizationMethod> >
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
struct NormalizeCovariance;
struct NormalizeCovarianceUnbiased;

typedef Normalize<Sum, true>                               Mean;
typedef Normalize<SumOfSquares, false, RootDivideByCount>  RootMeanSquares;

template <unsigned N>
class Moment;

template <unsigned N>
struct Synonym<Moment<N> >
{
    typedef Normalize<PowerSum<N>, true> type;
};

    // explicitly specialize Central<PowerSum<2> > to use incremental algorithms
    // (standard centralize accumulators use two passes)
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

#define VIGRA_NORMALIZE_REPLACEMENT(INNER, CACHE, METHOD) \
    Normalize<INNER, CACHE, METHOD>

VIGRA_SIMPLE_SYNONYM(StdDev, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, RootDivideByCount))
VIGRA_SIMPLE_SYNONYM(Central<StdDev>, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, RootDivideByCount))

VIGRA_SIMPLE_SYNONYM(UnbiasedVariance, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, DivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Central<UnbiasedVariance>, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, DivideUnbiased))

VIGRA_SIMPLE_SYNONYM(UnbiasedStdDev, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, RootDivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Central<UnbiasedStdDev>, VIGRA_NORMALIZE_REPLACEMENT(SSD, false, RootDivideUnbiased))

// template <>
// struct Synonym<Variance>
// {
    // typedef Central<Moment<2> > type;
// };

// template <>
// struct Synonym<StdDev>
// {
    // typedef Normalize<SSD, false, RootDivideByCount> type;
// };

// template <>
// struct Synonym<UnbiasedVariance>
// {
    // typedef Normalize<SSD, false, DivideUnbiased> type;
// };

// template <>
// struct Synonym<UnbiasedStdDev>
// {
    // typedef Normalize<SSD, false, RootDivideUnbiased> type;
// };

    // explicitly specialize Central<PowerSum<3> > and Central<PowerSum<4> > to implement 
    // merge function (standard centralize accumulators do not support merge)
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

typedef Normalize<FlatScatterMatrix, true, NormalizeCovariance>         Covariance;
typedef Normalize<FlatScatterMatrix, true, NormalizeCovarianceUnbiased> UnbiasedCovariance;

    // explicitly specialize to access CovarianceEigensystem directly
template <> class Principal<PowerSum<2> >;
template <> class Principal<Axes>;

VIGRA_SIMPLE_SYNONYM(Principal<Variance>, Principal<Moment<2> >)
VIGRA_SIMPLE_SYNONYM(Principal<StdDev>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, false, RootDivideByCount))
VIGRA_SIMPLE_SYNONYM(Principal<UnbiasedVariance>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, false, DivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Principal<UnbiasedStdDev>, VIGRA_NORMALIZE_REPLACEMENT(Principal<PowerSum<2> >, false, RootDivideUnbiased))
VIGRA_SIMPLE_SYNONYM(Principal<Skewness>, Principal<SkewnessImpl>)
VIGRA_SIMPLE_SYNONYM(Principal<Kurtosis>, Principal<Kurtosis>)

// template <>
// struct Synonym<Principal<Variance> >
// {
    // typedef Principal<Moment<2> > type;
// };

// template <>
// struct Synonym<Principal<StdDev> >
// {
    // typedef Normalize<Principal<PowerSum<2> >, false, RootDivideByCount> type;
// };

// template <>
// struct Synonym<Principal<UnbiasedVariance> >
// {
    // typedef Normalize<Principal<PowerSum<2> >, false, DivideUnbiased> type;
// };

// template <>
// struct Synonym<Principal<UnbiasedStdDev> >
// {
    // typedef Normalize<Principal<PowerSum<2> >, false, RootDivideUnbiased> type;
// };

// template <>
// struct Synonym<Principal<Skewness> >
// {
    // typedef Principal<SkewnessImpl> type;
// };

// template <>
// struct Synonym<Principal<Kurtosis> >
// {
    // typedef Principal<KurtosisImpl> type;
// };

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

    // normalization rule: never transfer normalization
template <class A, bool CacheResult, class NormlizationFunctor, class To>
struct DontTransferModifier<Normalize<A, CacheResult, NormlizationFunctor>, To>
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

    // modifying accumulators themselves must not be modified
VIGRA_DONT_TRANSFER(Central, Centralize)
VIGRA_DONT_TRANSFER(Central, PrincipalProjection)
VIGRA_DONT_TRANSFER(Central, RangeMapping)
VIGRA_DONT_TRANSFER(Principal, Centralize)
VIGRA_DONT_TRANSFER(Principal, PrincipalProjection)
VIGRA_DONT_TRANSFER(Principal, RangeMapping)

    // // we cannot make something principal that is already central and vice versa
// template <class A, class To>
// struct DontTransferModifier<Central<A>, To>
// {
    // static const int priority = AccumulatorTraits<Central<A> >::priority;
    // typedef typename ContainsPriority<priority, To>::type type;
// };

// template <class A, class To>
// struct DontTransferModifier<Principal<A>, To>
// {
    // static const int priority = AccumulatorTraits<Principal<A> >::priority;
    // typedef typename ContainsPriority<priority, To>::type type;
// };

#undef VIGRA_DONT_TRANSFER


#endif

/****************************************************************************/
/*                                                                          */
/*                             infrastructure                               */
/*                                                                          */
/****************************************************************************/

    // this class closes an accumulator chain
template <class T, class TAG, class NEXT>
struct AccumulatorBase 
{
	typedef TAG            Tag;
        // FIXME: this will become more sophisticated to support modifiers
    typedef typename Tag::Dependencies::type Dependencies;
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
    
	void merge(AccumulatorBase const &)
    {}
    
	void update(argument_type)
    {}
    
    void update(first_argument_type, second_argument_type)
    {}
};

/****************************************************************************/
/*                                                                          */
/*                        the actual accumulators                           */
/*                                                                          */
/****************************************************************************/

#if 0

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

template <unsigned N>
class PowerSum;

typedef PowerSum<0> Count;
typedef PowerSum<1> Sum;
typedef PowerSum<2> SumOfSquares;

    // normalization functors
struct DivideByCount;
struct RootDivideByCount;
struct DivideUnbiased;
struct RootDivideUnbiased;
struct NormalizeCovariance;
struct NormalizeCovarianceUnbiased;
    // Normalization modifier
    // FIXME: figure out when caching is appropriate
template <class A, bool CacheResult, class NormlizationFunctor=DivideByCount>
class Normalize;

typedef Normalize<Sum, true> Mean;
typedef Normalize<SumOfSquares, false, RootDivideByCount> RootMeanSquares;

    // FIXME: this may pose problems to modifier reordering
template <unsigned N>
class Moment
: Normalize<PowerSum<N>, true> 
{};

class SumOfSquaredDifferences;
typedef SumOfSquaredDifferences SSD;
typedef Normalize<SSD, false> Variance;
typedef Normalize<SSD, false, RootDivideByCount> StdDev;
typedef Normalize<SSD, false, DivideUnbiased> UnbiasedVariance;
typedef Normalize<SSD, false, RootDivideUnbiased> UnbiasedStdDev;

class FlatScatterMatrix;
typedef Normalize<FlatScatterMatrix, true, NormalizeCovariance> Covariance;
typedef Normalize<FlatScatterMatrix, true, NormalizeCovarianceUnbiased> UnbiasedCovariance;
class CovarianceEigensystem;
class PrincipalAxes;

    // last-seen-value accumulator that applies centralization
class Centralize;

    // modifier that forwards data after applying Centralize
template <class A>
class Central;

    // explicitly specialize to implement merge functions
    // (standard centralize accumulators will not support merge)
template <unsigned N>
class Central<PowerSum<N> >;

    // FIXME: this may pose problems to modifier reordering
template <unsigned N>
class Central<Moment<N> >
: Normalize<Central<PowerSum<N> >, true>
{};

class Skewness;
class Kurtosis;

    // last-seen-value accumulator that applies principal projection
class PrincipalProjection;

    // modifier that forwards data after applying PrincipalProjection
template <class A>
class Principal;

    // explicitly specialize to use 1-pass CovarianceEigensystem
template<>
class Principal<Variance>;
template<>
class Principal<StdDev>;
template<>
class Principal<UnbiasedVariance>;
template<>
class Principal<UnbiasedStdDev>;

template <unsigned N>
class Principal<Moment<N> >
: Normalize<Principal<PowerSum<N> >, true>
{};

    // data access modifier that forwards the coordinate part of a compound data type
template <class A>
class Coord;

    // data access modifier that forwards data/weight pairs from a compound data type
template <class A>
class Weighted;

    // data access modifier that forwards coordinate/weight pairs from a compound data type
template <class A>
class CoordWeighted;

template <unsigned Percent>
class IncrementalQuantile;

typedef IncrementalQuantile<0>   Minimum;
typedef IncrementalQuantile<50>  IncrementalMedian;
typedef IncrementalQuantile<100> Maximum;

    // last-seen-value accumulator that maps [min, max] to another range (e.g. for histogram creation)
template <class A>
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
typedef HistogramQuantile<50> HistogramMedian;

template <unsigned NDim>
class MultiHistogram;

class AccumulatorArray;

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

#endif

#if 0
class Count
{
  public:
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
            // All accumulators must contain appropriate versions of these typedefs.
        typedef double value_type;
        typedef double result_type;
        
        value_type value_;
        
            // All accumulators must implement a default constructor, reset(), merge() 
            // and operator()(). The default constructor must initialize the data.
            // If the shape of value_ can only be determined at runtime, reshape() 
            // must be implemented as well, resizing and initializing the data 
            // appropriately (see class Minimum below).
        Impl()
        : value_(0.0)
        {}
        
        void reset()
        {
            value_ = 0.0;
        }
        
        void merge(Impl const & o)
        {
            value_ += o.value_;
        }
        
        result_type operator()() const
        {
            return value_;
        }

            // All accumulators must implement either updatePass1() or updatePass2().
            // In the latter case, passesRequired() must be implemented as well, returning 2.
            // It is also possible to implement both update functions and decide at runtime
            // which one to use (e.g. depending on some initialization function). The
            // result of passesRequired() should then reflect this choice.
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

class Minimum
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
        
       void merge(Impl const & o)
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

class Maximum
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
        
        void merge(Impl const & o)
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

class Sum
{
  public:
    typedef Select<> Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type element_type;
        typedef typename AccumulatorResultTraits<T>::SumType              value_type;
        typedef value_type const &                                  result_type;

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
            detail::reshapeImpl(value_, s);
        }
        
        void merge(Impl const & o)
        {
            value_ += o.value_;
        }
    
        void update(T const & t)
        {
            value_ += t;
        }
        
        void update(T const & t, double weight)
        {
            value_ += weight*t;
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

class Mean
{
  public:
    typedef Select<Sum, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
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
            this->setClean();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            detail::reshapeImpl(value_, s);
        }

        void merge(Impl const &)
        {
            this->setDirty();
        }
    
        void update(T const &)
        {
            this->setDirty();
        }
        
        void update(T const & t, double)
        {
             this->setDirty();
       }
        
        result_type operator()() const
        {
			if(this->isDirty())
            {
                using namespace multi_math;
                value_ = get<Sum>(*this) / get<Count>(*this);
                this->setClean();
            }
            return value_;
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

namespace detail
{

template <unsigned int N>
struct CentralMomentsHelper
{
    template <class Left, class Right>
    static void merge(Left &, Right const &)
    {
        vigra_precondition(false,
            "CentralMoment<N>::merge(): not implemented for N > 4.");
    }
};

template <>
struct CentralMomentsHelper<2u>
{
    template <class Accu>
    static void merge(Accu & l, Accu const & r)
    {
        using namespace vigra::multi_math;
        double n1 = get<Count>(l), n2 = get<Count>(r);
        if(n1 == 0.0)
        {
            l.value_ = r.value_;
        }
        else if(n2 != 0.0)
        {
            l.value_ += r.value_ + n1 * n2 / (n1 + n2) * sq(get<Sum>(l) / n1 - get<Sum>(r) / n2);
        }
    }
};

template <>
struct CentralMomentsHelper<3u>
{
    template <class Accu>
    static void merge(Accu & l, Accu const & r);
};

template <class Accu>
void CentralMomentsHelper<3u>::merge(Accu & l, Accu const & r)
{
    using namespace vigra::multi_math;
    typedef typename LookupTag<Sum, Accu>::value_type SumType;
    double n1 = get<Count>(l), n2 = get<Count>(r);
    if(n1 == 0.0)
    {
        l.value_ = r.value_;
    }
    else if(n2 != 0.0)
    {
        double n = n1 + n2;
        double weight = n1 * n2 * (n1 - n2) / sq(n);
        SumType delta = get<Sum>(r) / n2 - get<Sum>(l) / n1;
        l.value_ += r.value_ + weight * pow(delta, 3) +
                       3.0 / n * delta * (  n1 * cast<CentralMoment<2> >(r).value_
                                          - n2 * cast<CentralMoment<2> >(l).value_);
    }
}

template <>
struct CentralMomentsHelper<4u>
{
    template <class Accu>
    static void merge(Accu & l, Accu const & r);
};

template <class Accu>
void CentralMomentsHelper<4u>::merge(Accu & l, Accu const & r)
{
    using namespace vigra::multi_math;
    typedef typename LookupTag<Sum, Accu>::value_type SumType;
    double n1 = get<Count>(l), n2 = get<Count>(r);
    if(n1 == 0.0)
    {
        l.value_ = r.value_;
    }
    else if(n2 != 0.0)
    {
        double n = n1 + n2;
        double n1_2 = sq(n1);
        double n2_2 = sq(n2);
        double n_2 = sq(n);
        double weight = n1 * n2 * (n1_2 - n1*n2 + n2_2) / n_2 / n;
        SumType delta = get<Sum>(r) / n2 - get<Sum>(l) / n1;
        l.value_ += r.value_ + weight * pow(delta, 4) +
                      6.0 / n_2 * sq(delta) * (  n1_2 * cast<CentralMoment<2> >(r).value_
                                               + n2_2 * cast<CentralMoment<2> >(l).value_ ) +
                      4.0 / n * delta * (  n1 * cast<CentralMoment<3> >(r).value_
                                         - n2 * cast<CentralMoment<3> >(l).value_);
    }
}

} // namsspace detail

template <unsigned int N>
class CentralMoment
{
  public:
    typedef Select<Centralize, typename IfBool<(N > 2), CentralMoment<N-1>, void>::type > Dependencies;
     
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename AccumulatorResultTraits<T>::element_promote_type element_type;
        typedef typename AccumulatorResultTraits<T>::SumType              value_type;
        typedef value_type                                          result_type;

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
            detail::reshapeImpl(value_, s);
        }
        
        void merge(Impl const & o)
        {
            detail::CentralMomentsHelper<N>::merge(*this, o);
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
        
        result_type operator()() const
        {
            using namespace vigra::multi_math;
            return value_ / get<Count>(*this);
        }
    };
};

class Skewness
{
  public:
    typedef Select<CentralMoment<3> > Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupTag<CentralMoment<3>, BASE>::value_type value_type;
        typedef value_type                                             result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<Count>(*this)) * 
                     cast<CentralMoment<3> >(*this).value_ / pow(cast<CentralMoment<2> >(*this).value_, 1.5);
        }
    };
};

class Kurtosis
{
  public:
    typedef Select<CentralMoment<4> > Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupTag<CentralMoment<4>, BASE>::value_type value_type;
        typedef value_type                                             result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return get<Count>(*this) * 
                     cast<CentralMoment<4> >(*this).value_ / sq(cast<CentralMoment<2> >(*this).value_);
        }
    };
};

class SumSquaredDifferences
{
  public:
    typedef Select<Mean, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<T>::SumType               value_type;
        typedef value_type const &                                   result_type;
       
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
            detail::reshapeImpl(value_, s);
        }
        
        void merge(Impl const & o)
        {
            detail::CentralMomentsHelper<2>::merge(*this, o);
        }
    
        void update(T const & t)
        {
            double n = get<Count>(*this);
            if(n > 1.0)
            {
                using namespace vigra::multi_math;
                value_ += n / (n - 1.0) * sq(get<Sum>(*this) / n - t);
            }
        }
        
        void update(T const & t, double weight)
        {
            double n = get<Count>(*this);
            if(n > weight)
            {
                using namespace vigra::multi_math;
                value_ += n / (n - weight) * sq(get<Sum>(*this) / n - t);
            }
        }
        
        result_type operator()() const
        {
            return value_;
        }
    };
};

typedef SumSquaredDifferences SSD;

class Variance
{
  public:
    typedef Select<SumSquaredDifferences, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<SumSquaredDifferences, BASE>::value_type value_type;
        typedef value_type                                                  result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return get<SumSquaredDifferences>(*this) / get<Count>(*this);
        }
    };
};

class StdDev
{
  public:
    typedef Select<Variance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<Variance, BASE>::value_type value_type;
        typedef value_type                                     result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<SumSquaredDifferences>(*this) / get<Count>(*this));
        }
    };
};

class UnbiasedVariance
{
  public:
    typedef Select<SumSquaredDifferences, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<SumSquaredDifferences, BASE>::value_type value_type;
        typedef value_type                                                  result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return get<SumSquaredDifferences>(*this) / (get<Count>(*this) - 1.0);
        }
    };
};

class UnbiasedStdDev
{
  public:
    typedef Select<UnbiasedVariance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename LookupTag<Variance, BASE>::value_type value_type;
        typedef value_type                                     result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<SumSquaredDifferences>(*this) / (get<Count>(*this) - 1.0));
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
class FlatScatterMatrix
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
        
        void merge(Impl const & o)
        {
            double n1 = get<Count>(*this), n2 = get<Count>(o);
            if(n1 == 0.0)
            {
                value_ = o.value_;
            }
            else if(n2 != 0.0)
            {
                using namespace vigra::multi_math;
                diff_ = get<Sum>(*this) / n1 - get<Sum>(o) / n2;
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
                diff_ = get<Sum>(*this) / n - t;
                detail::updateFlatScatterMatrix(value_, diff_, n * weight / (n - weight));
            }
        }
    };
};

class Covariance
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<T>::CovarianceType        value_type;
        typedef value_type const &                                   result_type;

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
            int size = prod(s);
            detail::reshapeImpl(value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            detail::flatScatterMatrixToCovariance(value_, cast<FlatScatterMatrix>(*this).value_, get<Count>(*this));
            return value_;
        }
    };
};

class CovarianceEigensystem
{
  public:
    typedef Select<Covariance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type        element_type;
        typedef typename AccumulatorResultTraits<T>::SumType                     EigenvalueType;
        typedef typename AccumulatorResultTraits<T>::CovarianceType              EigenvectorType;
        typedef EigenvectorType                                            value_type;
        typedef std::pair<EigenvalueType const &, EigenvectorType const &> result_type;

        mutable EigenvalueType eigenvalues_;
        mutable EigenvectorType eigenvectors_;
        
        Impl()
        : eigenvalues_(),  // call default constructor explicitly to ensure zero initialization
          eigenvectors_()
        {}
        
        void reset()
        {
            eigenvalues_ = element_type();
            eigenvectors_ = element_type();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            detail::reshapeImpl(eigenvalues_, Shape2(size,1));
            detail::reshapeImpl(eigenvectors_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            compute(get<Covariance>(*this), eigenvalues_, eigenvectors_);
            return result_type(eigenvalues_, eigenvectors_);
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

class UnbiasedCovariance
{
  public:
    typedef Select<FlatScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename AccumulatorResultTraits<T>::element_promote_type  element_type;
        typedef typename AccumulatorResultTraits<T>::CovarianceType        value_type;
        typedef value_type const &                                   result_type;

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
            int size = prod(s);
            detail::reshapeImpl(value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            detail::flatScatterMatrixToCovariance(value_, cast<FlatScatterMatrix>(*this).value_, get<Count>(*this)-1.0);
            return value_;
        }
    };
};

#endif

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
