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
#include "multi_iterator_coupled.hxx"
#include "matrix.hxx"
#include "multi_math.hxx"
#include "eigensystem.hxx"
#include "histogram.hxx"
#include <algorithm>
#include <iostream>

namespace vigra {

namespace acc1 {

/****************************************************************************/
/*                                                                          */
/*                             infrastructure                               */
/*                                                                          */
/****************************************************************************/

    // Select is a wrapper for MakeTypeList that additionally performs tag standardization
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
    
    static std::string const & name() 
    { 
        static const std::string n("AccumulatorBegin (internal)");
        return n;
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

struct Error__Global_statistics_are_only_defined_for_AccumulatorChainArray;

template <int INDEX>
class LabelArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("LabelArg<") + asString(INDEX) + "> (internal)";
        return n;
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
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("CoordArg<") + asString(INDEX) + "> (internal)";
        return n;
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

namespace detail {

/****************************************************************************/
/*                                                                          */
/*                   internal tag handling meta-functions                   */
/*                                                                          */
/****************************************************************************/

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
    typedef typename Inner::RegionTags                  RegionTags;
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
        static const std::string name = normalizeString(HEAD::name());
        if(name == tag)
        {
            v.exec<HEAD>(a);
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
    void merge(U const &) 
    {}
    
    template <class U>
    void resize(U const &) 
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

    static void merge(A & a, A const & o)
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
        return std::max(A::workInPass, A::InternalBaseType::passesRequired());
    }
};

template <class A, unsigned CurrentPass, bool allowRuntimeActivation>
struct DecoratorImpl<A, CurrentPass, allowRuntimeActivation, CurrentPass>
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
        static const std::string message = std::string("get(accumulator): attempt to access inactive statistic '") +
                                                                                   typeid(typename A::Tag).name() + "'.";
        vigra_precondition(isActive(a), message);
        return a();
    }

    static void merge(A & a, A const & o)
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
        return A::isActiveImpl(flags)
                   ? std::max(A::workInPass, A::InternalBaseType::passesRequired(flags))
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
    
    GlobalAccumulatorChain next_;
    RegionAccumulatorArray regions_;
    HistogramOptions region_histogram_options_;
    MultiArrayIndex ignore_label_;
    ActiveFlagsType active_region_accumulators_;
    
    template <class IndexDefinition, class TagFound=typename IndexDefinition::Tag>
    struct LabelIndexSelector
    {
        static const int value = 2; // default: CoupledHandle holds labels at index 2

        template <class U, class NEXT>
        static MultiArrayIndex exec(CoupledHandle<U, NEXT> const & t)
        {
            return (MultiArrayIndex)get<value>(t); 
        }
    };
    
    template <class IndexDefinition>
    struct LabelIndexSelector<IndexDefinition, LabelArgTag>
    {
        static const int value = IndexDefinition::value;

        template <class U, class NEXT>
        static MultiArrayIndex exec(CoupledHandle<U, NEXT> const & t)
        {
            return (MultiArrayIndex)get<value>(t);
        }
    };
    
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
    
    typedef typename LookupTag<LabelArgTag, GlobalAccumulatorChain>::type FindLabelIndex;
    
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
        }
    }
    
    void ignoreLabel(MultiArrayIndex l)
    {
        ignore_label_ = l;
    }
    
    void applyHistogramOptions(HistogramOptions const & options)
    {
        applyHistogramOptions(options, options);
    }
    
    void applyHistogramOptions(HistogramOptions const & regionoptions, HistogramOptions const & globaloptions)
    {
        region_histogram_options_ = regionoptions;
        next_.applyHistogramOptions(globaloptions);
    }
    
    template <class U>
    void resize(U const & t)
    {
        if(regions_.size() == 0)
        {
            static const int labelIndex = LabelIndexSelector<FindLabelIndex>::value;
            typedef typename CoupledHandleCast<labelIndex, T>::type LabelHandle;
            typedef typename LabelHandle::value_type LabelType;
            typedef MultiArrayView<LabelHandle::dimensions, LabelType, StridedArrayTag> LabelArray;
            LabelArray labelArray(t.shape(), cast<labelIndex>(t).strides(), const_cast<LabelType *>(cast<labelIndex>(t).ptr()));
            
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
        if(LabelIndexSelector<FindLabelIndex>::exec(t) != ignore_label_)
        {
            next_.pass<N>(t);
            regions_[LabelIndexSelector<FindLabelIndex>::exec(t)].pass<N>(t);
        }
    }
    
    template <unsigned N>
    void pass(T const & t, double weight)
    {
        if(LabelIndexSelector<FindLabelIndex>::exec(t) != ignore_label_)
        {
            next_.pass<N>(t, weight);
            regions_[LabelIndexSelector<FindLabelIndex>::exec(t)].pass<N>(t, weight);
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
    
    void merge(LabelDispatch const & o)
    {
        for(unsigned int k=0; k<regions_.size(); ++k)
            regions_[k].merge(o.regions_[k]);
        next_.merge(o.next_);
    }
    
    void merge(unsigned i, unsigned j)
    {
        regions_[i].merge(regions_[j]);
        regions_[j].reset();
        getAccumulator<AccumulatorEnd>(regions_[j]).active_accumulators_ = active_region_accumulators_;
    }
    
    template <class ArrayLike>
    void merge(LabelDispatch const & o, ArrayLike const & labelMapping)
    {
        MultiArrayIndex newMaxLabel = std::max<MultiArrayIndex>(maxRegionLabel(), *argMax(labelMapping.begin(), labelMapping.end()));
        setMaxRegionLabel(newMaxLabel);
        for(unsigned int k=0; k<labelMapping.size(); ++k)
            regions_[labelMapping[k]].merge(o.regions_[k]);
        next_.merge(o.next_);
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
            detail::ActivateDependencies<StdDeps>::template exec<ThisType>(flags);
        }
        
        template <class Accu, class ActiveFlags, class GlobalFlags>
        static void activateImpl(ActiveFlags & flags, GlobalFlags & gflags)
        {
            flags.template set<index>();
            typedef typename StandardizeDependencies<Tag>::type StdDeps;
            detail::ActivateDependencies<StdDeps>::template exec<Accu>(flags, gflags);
        }
        
        template <class ActiveFlags>
        static bool isActiveImpl(ActiveFlags & flags)
        {
            return flags.template test<index>();
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
            this->next_.pass<N>(t);
            DecoratorImpl<Accumulator, N, allowRuntimeActivation>::exec(*this, t);
        }
        
        template <unsigned N, class T>
        void pass(T const & t, double weight)
        {
            this->next_.pass<N>(t, weight);
            DecoratorImpl<Accumulator, N, allowRuntimeActivation>::exec(*this, t, weight);
        }
        
        void merge(Accumulator const & o)
        {
            DecoratorImpl<Accumulator, Accumulator::workInPass, allowRuntimeActivation>::merge(*this, o);
            this->next_.merge(o.next_);
        }
        
        void applyHistogramOptions(HistogramOptions const & options)
        {
            DecoratorImpl<Accumulator, workInPass, allowRuntimeActivation>::applyHistogramOptions(*this, options);
            this->next_.applyHistogramOptions(options);
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
: public ConfigureAccumulatorChain<T, typename AddDependencies<typename Selected::type>::type, dynamic>
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
: public ConfigureAccumulatorChainArray<T, typename AddDependencies<typename Selected::type>::type, dynamic>
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

} // namespace detail 

/****************************************************************************/
/*                                                                          */
/*                            accumulator chain                             */
/*                                                                          */
/****************************************************************************/

    // Implement the high-level interface of an accumulator chain
template <class T, class NEXT>
struct AccumulatorChainImpl
{
    typedef NEXT                                             InternalBaseType;
    typedef AccumulatorBegin                                 Tag;
    typedef typename InternalBaseType::argument_type         argument_type;
    typedef typename InternalBaseType::first_argument_type   first_argument_type;
    typedef typename InternalBaseType::second_argument_type  second_argument_type;
    typedef void                                             value_type;
    typedef typename InternalBaseType::result_type           result_type;
    
    static const int staticSize = InternalBaseType::index;

    InternalBaseType next_;
    unsigned int current_pass_;
    
    AccumulatorChainImpl()
    : current_pass_(0)
    {}
    
    void setHistogramOptions(HistogramOptions const & options)
    {
        next_.applyHistogramOptions(options);
    }
    
    void setHistogramOptions(HistogramOptions const & regionoptions, HistogramOptions const & globaloptions)
    {
        next_.applyHistogramOptions(regionoptions, globaloptions);
    }
    
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
                next_.resize(detail::shapeOf(t));
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
                next_.resize(detail::shapeOf(t));
            next_.template pass<N>(t, weight);
        }
        else
        {
            std::string message("AccumulatorChain::update(): cannot return to pass ");
            message << N << " after working on pass " << current_pass_ << ".";
            vigra_precondition(false, message);
       }
    }
    
    void operator+=(AccumulatorChainImpl const & o)
    {
        merge(o);
    }
    
    void merge(AccumulatorChainImpl const & o)
    {
        next_.merge(o.next_);
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
    
    unsigned int passesRequired() const
    {
        return InternalBaseType::passesRequired();
    }
};

   // Create an accumulator chain containing the Selected statistics and their dependencies.
template <class T, class Selected, bool dynamic=false>
struct AccumulatorChain
: public AccumulatorChainImpl<T, typename detail::ConfigureAccumulatorChain<T, Selected, dynamic>::type>
{
    typedef typename detail::ConfigureAccumulatorChain<T, Selected, dynamic>::TagList AccumulatorTags;
    
    template <class U, int N>
    void reshape(TinyVector<U, N> const & s)
    {
        vigra_precondition(this->current_pass_ == 0,
             "AccumulatorChain::reshape(): cannot reshape after seeing data. Call AccumulatorChain::reset() first.");
        this->next_.resize(s);
        this->current_pass_ = 1;
    }
     
    static ArrayVector<std::string> const & tagNames()
    {
        static const ArrayVector<std::string> n = collectTagNames();
        return n;
    }
    
  private:
    static ArrayVector<std::string> collectTagNames()
    {
        ArrayVector<std::string> n;
        detail::CollectAccumulatorNames<AccumulatorTags>::exec(n);
        std::sort(n.begin(), n.end());
        return n;
    }
};   

    // Create a dynamic accumulator chain containing the Selected statistics and their dependencies.
    // Statistics will only be computed if activate<Tag>() is called at runtime.
template <class T, class Selected>
struct DynamicAccumulatorChain
: public AccumulatorChain<T, Selected, true>
{
    typedef typename AccumulatorChain<T, Selected, true>::InternalBaseType InternalBaseType;
    typedef typename DynamicAccumulatorChain::AccumulatorTags AccumulatorTags;
        
    void activate(std::string tag)
    {
        vigra_precondition(activateImpl(tag),
            std::string("DynamicAccumulatorChain::activate(): Tag '") + tag + "' not found.");
    }
    
    template <class TAG>
    void activate()
    {
        LookupTag<TAG, DynamicAccumulatorChain>::type::activateImpl(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }
    
    void activateAll()
    {
        getAccumulator<AccumulatorEnd>(*this).active_accumulators_.set();
    }
    
    bool isActive(std::string tag) const
    {
        detail::TagIsActive_Visitor v;
        vigra_precondition(isActiveImpl(tag, v),
            std::string("DynamicAccumulatorChain::isActive(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    template <class TAG>
    bool isActive() const
    {
        return LookupTag<TAG, DynamicAccumulatorChain>::type::isActiveImpl(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }
    
    ArrayVector<std::string> activeNames() const
    {
        ArrayVector<std::string> res;
        for(unsigned k=0; k<DynamicAccumulatorChain::tagNames().size(); ++k)
            if(isActive(DynamicAccumulatorChain::tagNames()[k]))
                res.push_back(DynamicAccumulatorChain::tagNames()[k]);
        return res;
    }
    
    unsigned int passesRequired() const
    {
        return InternalBaseType::passesRequired(getAccumulator<AccumulatorEnd>(*this).active_accumulators_);
    }
    
  protected:
  
    bool activateImpl(std::string tag)
    {
        return detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, 
                                         normalizeString(tag), detail::ActivateTag_Visitor());
    }
    
    bool isActiveImpl(std::string tag, detail::TagIsActive_Visitor & v) const
    {
        return detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, normalizeString(tag), v);
    }
};

template <class T, class Selected, bool dynamic=false>
struct AccumulatorChainArray
: public AccumulatorChainImpl<T, typename detail::ConfigureAccumulatorChainArray<T, Selected, dynamic>::type>
{
    typedef typename detail::ConfigureAccumulatorChainArray<T, Selected, dynamic> Creator;
    typedef typename Creator::TagList AccumulatorTags;
    typedef typename Creator::GlobalTags GlobalTags;
    typedef typename Creator::RegionTags RegionTags;
    
    void ignoreLabel(MultiArrayIndex l)
    {
        this->next_.ignoreLabel(l);
    }
    
    void setMaxRegionLabel(unsigned label)
    {
        this->next_.setMaxRegionLabel(label);
    }
    
    MultiArrayIndex maxRegionLabel() const
    {
        return this->next_.maxRegionLabel();
    }
     
    unsigned int regionCount() const
    {
        return this->next_.regions_.size();
    }
     
    void merge(unsigned i, unsigned j)
    {
        vigra_precondition(i <= maxRegionLabel() && j <= maxRegionLabel(),
            "AccumulatorChainArray::merge(): region labels out of range.");
        this->next_.merge(i, j);
    }
    
    void merge(AccumulatorChainArray const & o)
    {
        vigra_precondition(maxRegionLabel() == o.maxRegionLabel(),
            "AccumulatorChainArray::merge(): maxRegionLabel must be equal.");
        this->next_.merge(o.next_);
    }

    template <class ArrayLike>
    void merge(AccumulatorChainArray const & o, ArrayLike const & labelMapping)
    {
        vigra_precondition(labelMapping.size() == o.regionCount(),
            "AccumulatorChainArray::merge(): labelMapping.size() must match regionCount() of RHS.");
        this->next_.merge(o.next_, labelMapping);
    }
    
    static ArrayVector<std::string> const & tagNames()
    {
        static const ArrayVector<std::string> n = collectTagNames();
        return n;
    }
    
  private:
    static ArrayVector<std::string> collectTagNames()
    {
        ArrayVector<std::string> n;
        detail::CollectAccumulatorNames<AccumulatorTags>::exec(n);
        std::sort(n.begin(), n.end());
        return n;
    }
};   

template <class T, class Selected>
struct DynamicAccumulatorChainArray
: public AccumulatorChainArray<T, Selected, true>
{
    typedef typename DynamicAccumulatorChainArray::AccumulatorTags AccumulatorTags;
        
    void activate(std::string tag)
    {
        vigra_precondition(activateImpl(tag),
            std::string("DynamicAccumulatorChainArray::activate(): Tag '") + tag + "' not found.");
    }

    template <class TAG>
    void activate()
    {
        this->next_.activate<TAG>();
    }
    
    void activateAll()
    {
        this->next_.activateAll();
    }
    
    bool isActive(std::string tag) const
    {
        detail::TagIsActive_Visitor v;
        vigra_precondition(isActiveImpl(tag, v),
            std::string("DynamicAccumulatorChainArray::isActive(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    template <class TAG>
    bool isActive() const
    {
        return this->next_.isActive<TAG>();
    }
    
    ArrayVector<std::string> activeNames() const
    {
        ArrayVector<std::string> res;
        for(unsigned k=0; k<DynamicAccumulatorChainArray::tagNames().size(); ++k)
            if(isActive(DynamicAccumulatorChainArray::tagNames()[k]))
                res.push_back(DynamicAccumulatorChainArray::tagNames()[k]);
        return res;
    }
    
    unsigned int passesRequired() const
    {
        return this->next_.passesRequiredDynamic();
    }
    
  protected:
  
    bool activateImpl(std::string tag)
    {
        return detail::ApplyVisitorToTag<AccumulatorTags>::exec(this->next_, 
                                         normalizeString(tag), detail::ActivateTag_Visitor());
    }
    
    bool isActiveImpl(std::string tag, detail::TagIsActive_Visitor & v) const
    {
        return detail::ApplyVisitorToTag<AccumulatorTags>::exec(this->next_, normalizeString(tag), v);
    }
};

/****************************************************************************/
/*                                                                          */
/*                        generic access functions                          */
/*                                                                          */
/****************************************************************************/

template <class TAG>
struct Error__Attempt_to_access_inactive_statistic;

namespace detail {

    // accumulator lookup rules: find the accumulator that implements TAG
    
    // When A does not implement TAG, continue search in A::InternalBaseType.
template <class TAG, class A, class FromTag=typename A::Tag>
struct LookupTagImpl
: public LookupTagImpl<TAG, typename A::InternalBaseType>
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

} // namespace detail

    // Lookup the accumulator in the chain A that implements the given TAG.
template <class Tag, class A>
struct LookupTag
: public detail::LookupTagImpl<typename StandardizeTag<Tag>::type, A>
{};

    // Lookup the dependency TAG of the accumulator A.
    // This template ensures that dependencies are used with matching modifiers.
    // Specifically, if you search for Count as a dependency of Weighted<Mean>, the search
    // actually returns Weighted<Count>, wheras Count will be returned for plain Mean.
template <class Tag, class A, class TargetTag>
struct LookupDependency
: public detail::LookupTagImpl<
       typename TransferModifiers<TargetTag, typename StandardizeTag<Tag>::type>::type, A>
{};
 

namespace detail {

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

} // namespace detail

    // Get a reference to the accumulator TAG in the accumulator chain A
template <class TAG, class A>
inline typename LookupTag<TAG, A>::reference
getAccumulator(A & a)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::reference reference;
    return detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a);
}

    // Get a reference to the accumulator TAG for region 'label' in the accumulator chain A
template <class TAG, class A>
inline typename LookupTag<TAG, A>::reference
getAccumulator(A & a, MultiArrayIndex label)
{
    typedef typename LookupTag<TAG, A>::Tag StandardizedTag;
    typedef typename LookupTag<TAG, A>::reference reference;
    return detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a, label);
}

    // get the result of the accumulator specified by TAG
template <class TAG, class A>
inline typename LookupTag<TAG, A>::result_type
get(A const & a)
{
    return getAccumulator<TAG>(a).get();
}

    // get the result of the accumulator TAG for region 'label'
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
    return detail::CastImpl<StandardizedTag, typename A::Tag, reference>::exec(a)();
}

    // activate the dynamic accumulator specified by Tag
template <class Tag, class A>
inline void
activate(A & a)
{
    a.activate<Tag>();
}

    // check if the dynamic accumulator specified by Tag is active
template <class Tag, class A>
inline bool
isActive(A const & a)
{
    return a.isActive<Tag>();
}

/****************************************************************************/
/*                                                                          */
/*                               generic loops                              */
/*                                                                          */
/****************************************************************************/

template <class ITERATOR, class ACCUMULATOR>
void collectStatistics(ITERATOR start, ITERATOR end, ACCUMULATOR & a)
{
    for(unsigned int k=1; k <= a.passesRequired(); ++k)
        for(ITERATOR i=start; i < end; ++i)
            a.updatePassN(*i, k);
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

template <class TAG>
class Global
{
  public:
    typedef typename StandardizeTag<TAG>::type  TargetTag;
    typedef typename TargetTag::Dependencies    Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("Global<") + TargetTag::name() + " >";
        return n;
    }
};

template <int INDEX>
class DataArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("DataArg<") + asString(INDEX) + "> (internal)";
        return n;
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

template <class TAG>
class DataFromHandle
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef typename TargetTag::Dependencies Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("DataFromHandle<") + TargetTag::name() + " > (internal)";
        return n;
    }
    
    template <class IndexDefinition, class TagFound=typename IndexDefinition::Tag>
    struct DataIndexSelector
    {
        static const int value = 1; // default: CoupledHandle holds data at index 1 
        
        template <class U, class NEXT>
        static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference 
        exec(CoupledHandle<U, NEXT> const & t)
        {
            return get<value>(t);
        }
    };
    
    template <class IndexDefinition>
    struct DataIndexSelector<IndexDefinition, DataArgTag>
    {
        static const int value = IndexDefinition::value;
        
        template <class U, class NEXT>
        static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference
        exec(CoupledHandle<U, NEXT> const & t)
        {
            return get<value>(t);
        }
    };
    
    template <class T, class BASE>
    struct SelectInputType
    {
        typedef typename LookupTag<DataArgTag, BASE>::type FindDataIndex;
        typedef DataIndexSelector<FindDataIndex> DataIndex;
        typedef typename CoupledHandleCast<DataIndex::value, T>::type::value_type type;
    };
    
    template <class T, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename SelectInputType<T, BASE>::type, BASE>
    {
        typedef SelectInputType<T, BASE>                InputTypeSelector;
        typedef typename InputTypeSelector::DataIndex   DataIndex;
        typedef typename InputTypeSelector::type        input_type;
        typedef input_type const &                      argument_type;
        typedef argument_type                           first_argument_type;
        
        typedef typename TargetTag::template Impl<input_type, BASE> ImplType;
        
        using ImplType::reshape;
        
        template <class U, class NEXT>
        void reshape(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::reshape(detail::shapeOf(DataIndex::exec(t)));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::update(DataIndex::exec(t));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            ImplType::update(DataIndex::exec(t), weight);
        }
    };
};

template <class TAG>
class Coord
{
  public:
    typedef typename StandardizeTag<TAG>::type   TargetTag;
    typedef typename TargetTag::Dependencies     Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("Coord<") + TargetTag::name() + " >";
        return n;
    }
    
    template <class IndexDefinition, class TagFound=typename IndexDefinition::Tag>
    struct CoordIndexSelector
    {
        static const int value = 0; // default: CoupledHandle holds coordinates at index 0 
        
        template <class U, class NEXT>
        static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference 
        exec(CoupledHandle<U, NEXT> const & t)
        {
            return get<value>(t);
        }
    };
    
    template <class IndexDefinition>
    struct CoordIndexSelector<IndexDefinition, CoordArgTag>
    {
        static const int value = IndexDefinition::value;
        
        template <class U, class NEXT>
        static typename CoupledHandleCast<value, CoupledHandle<U, NEXT> >::type::const_reference
        exec(CoupledHandle<U, NEXT> const & t)
        {
            return get<value>(t);
        }
    };
     
    template <class T, class BASE>
    struct SelectInputType
    {
        typedef typename LookupTag<CoordArgTag, BASE>::type FindDataIndex;
        typedef CoordIndexSelector<FindDataIndex> CoordIndex;
        typedef typename CoupledHandleCast<CoordIndex::value, T>::type::value_type type;
    };
    
    template <class T, class BASE>
    struct Impl
    : public TargetTag::template Impl<typename SelectInputType<T, BASE>::type, BASE>
    {
        typedef SelectInputType<T, BASE>                InputTypeSelector;
        typedef typename InputTypeSelector::CoordIndex  CoordIndex;
        typedef typename InputTypeSelector::type        input_type;
        typedef input_type const &                      argument_type;
        typedef argument_type                           first_argument_type;
        
        typedef typename TargetTag::template Impl<input_type, BASE> ImplType;
        
        using ImplType::reshape;
        
        template <class U, class NEXT>
        void reshape(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::reshape(detail::shapeOf(CoordIndex::exec(t)));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t)
        {
            ImplType::update(CoordIndex::exec(t));
        }
        
        template <class U, class NEXT>
        void update(CoupledHandle<U, NEXT> const & t, double weight)
        {
            ImplType::update(CoordIndex::exec(t), weight);
        }
    };
};

template <int INDEX>
class WeightArg
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("WeightArg<") + asString(INDEX) + "> (internal)";
        return n;
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

template <class TAG>
class Weighted
{
  public:
    typedef typename StandardizeTag<TAG>::type   TargetTag;
    typedef typename TargetTag::Dependencies     Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("Weighted<") + TargetTag::name() + " >";
        return n;
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
    
    static std::string const & name() 
    { 
        static const std::string n("Centralize (internal)");
        return n;
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
            detail::reshapeImpl(value_, s);
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

template <class TAG>
class Central
{
  public:
    typedef typename StandardizeTag<TAG>::type                    TargetTag;
    typedef Select<Centralize, typename TargetTag::Dependencies>  Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("Central<") + TargetTag::name() + " >";
        return n;
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

// Compute principal projection and cache the result
class PrincipalProjection
{
  public:
    typedef Select<Centralize, Principal<CoordinateSystem> > Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("PrincipalProjection (internal)");
        return n;
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
            detail::reshapeImpl(value_, s);
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

template <class TAG>
class Principal
{
  public:
    typedef typename StandardizeTag<TAG>::type                             TargetTag;
    typedef Select<PrincipalProjection, typename TargetTag::Dependencies>  Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("Principal<") + TargetTag::name() + " >";
        return n;
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

class CoordinateSystem
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("CoordinateSystem");
        return n;
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
            detail::reshapeImpl(value_, s);
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
    
    static std::string const & name() 
    { 
        static const std::string n("PowerSum<0>");
        return n;
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
     
    static std::string const & name() 
    { 
        static const std::string n("PowerSum<1>");
        return n;
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
            this->value_ += weight*t;
        }
    };
};

template <unsigned N>
class PowerSum
{
  public:
    typedef Select<> Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n = std::string("PowerSum<") + asString(N) + ">";
        return n;
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
     
    static std::string const & name() 
    { 
        static const std::string n("AbsPowerSum<1>");
        return n;
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

template <unsigned N>
class AbsPowerSum
{
  public:
    typedef Select<> Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n = std::string("AbsPowerSum<") + asString(N) + ">";
        return n;
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
        detail::reshapeImpl(value_, s);
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
template <class TAG>
class DivideByCount
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef Select<TargetTag, Count> Dependencies;
  
    static std::string const & name() 
    { 
        static const std::string n = std::string("DivideByCount<") + TargetTag::name() + " >";
        return n;
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
template <class TAG>
class DivideUnbiased
{
  public:
    typedef typename StandardizeTag<TAG>::type TargetTag;
    typedef Select<TargetTag, Count> Dependencies;
      
    static std::string const & name() 
    { 
        static const std::string n = std::string("DivideUnbiased<") + TargetTag::name() + " >";
        return n;
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
template <class TAG>
class RootDivideByCount
{
  public:
    typedef typename StandardizeTag<DivideByCount<TAG> >::type TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    static std::string const & name() 
    { 
        typedef typename StandardizeTag<TAG>::type InnerTag;
        static const std::string n = std::string("RootDivideByCount<") + InnerTag::name() + " >";
        return n;
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
template <class TAG>
class RootDivideUnbiased
{
  public:
    typedef typename StandardizeTag<DivideUnbiased<TAG> >::type TargetTag;
    typedef Select<TargetTag> Dependencies;
    
    static std::string const & name() 
    { 
        typedef typename StandardizeTag<TAG>::type InnerTag;
        static const std::string n = std::string("RootDivideUnbiased<") + InnerTag::name() + " >";
        return n;
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

template <>
class Central<PowerSum<2> >
{
  public:
    typedef Select<Mean, Count> Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n("Central<PowerSum<2> >");
        return n;
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

template <>
class Central<PowerSum<3> >
{
  public:
    typedef Select<Centralize, Count, Mean, Central<PowerSum<2> > > Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n("Central<PowerSum<3> >");
        return n;
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

template <>
class Central<PowerSum<4> >
{
  public:
    typedef Select<Centralize, Central<PowerSum<3> > > Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n("Central<PowerSum<4> >");
        return n;
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

class Skewness
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<3> > > Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("Skewness");
        return n;
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

class UnbiasedSkewness
{
  public:
    typedef Select<Skewness> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("UnbiasedSkewness");
        return n;
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

class Kurtosis
{
  public:
    typedef Select<Central<PowerSum<2> >, Central<PowerSum<4> > > Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("Kurtosis");
        return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<4> >, BASE>::value_type value_type;
        typedef value_type                                                  result_type;

        result_type operator()() const
        {
            typedef Central<PowerSum<4> > Sum4;
            typedef Central<PowerSum<2> > Sum2;
        
			using namespace multi_math;
            return getDependency<Count>(*this) * getDependency<Sum4>(*this) / sq(getDependency<Sum2>(*this)) - value_type(3.0);
        }
    };
};

class UnbiasedKurtosis
{
  public:
    typedef Select<Kurtosis> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("UnbiasedKurtosis");
        return n;
    }
    
    template <class U, class BASE>
    struct Impl
    : public BASE
    {
        static const unsigned int workInPass = 2;
        
        typedef typename LookupDependency<Central<PowerSum<4> >, BASE>::value_type value_type;
        typedef value_type                                                  result_type;

        result_type operator()() const
        {
			using namespace multi_math;
            double n = getDependency<Count>(*this);
            return (n-1.0)/((n-2.0)*(n-3.0))*((n+1.0)*getDependency<Kurtosis>(*this) + value_type(6.0));
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

} // namespace detail

// we only store the flattened upper triangular part of the scatter matrix
class FlatScatterMatrix
{
  public:
    typedef Select<Mean, Count> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("FlatScatterMatrix");
        return n;
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
            detail::reshapeImpl(value_, Shape1(size*(size+1)/2));
            detail::reshapeImpl(diff_, s);
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
                detail::updateFlatScatterMatrix(value_, diff_, n1 * n2 / (n1 + n2));
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
                detail::updateFlatScatterMatrix(value_, diff_, n * weight / (n - weight));
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
    
    static std::string const & name() 
    { 
        static const std::string n("DivideByCount<FlatScatterMatrix>");
        return n;
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
            detail::reshapeImpl(this->value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                detail::flatScatterMatrixToCovariance(this->value_, getDependency<FlatScatterMatrix>(*this), getDependency<Count>(*this));
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
    
    static std::string const & name() 
    { 
        static const std::string n("DivideUnbiased<FlatScatterMatrix>");
        return n;
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
            detail::reshapeImpl(this->value_, Shape2(size,size));
        }
        
        result_type operator()() const
        {
            if(this->isDirty())
            {
                detail::flatScatterMatrixToCovariance(this->value_, getDependency<FlatScatterMatrix>(*this), getDependency<Count>(*this) - 1.0);
                this->setClean();
            }
            return this->value_;
        }
    };
};

class ScatterMatrixEigensystem
{
  public:
    typedef Select<FlatScatterMatrix> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("ScatterMatrixEigensystem");
        return n;
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
            value_.second = element_type();
            this->setClean();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            int size = prod(s);
            detail::reshapeImpl(value_.first, Shape1(size));
            detail::reshapeImpl(value_.second, Shape2(size,size));
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
            detail::flatScatterMatrixToScatterMatrix(scatter, flatScatter);
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
    
    static std::string const & name() 
    { 
        static const std::string n("DivideByCount<ScatterMatrixEigensystem>");
        return n;
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
            detail::reshapeImpl(value_.first, Shape2(size,1));
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
            // detail::reshapeImpl(value_.first, Shape2(size,1));
            // detail::reshapeImpl(value_.second, Shape2(size,size));
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
template <>
class Principal<PowerSum<2> >
{
  public:
    typedef Select<ScatterMatrixEigensystem> Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n("Principal<PowerSum<2> >");
        return n;
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
template <>
class Principal<CoordinateSystem>
{
  public:
    typedef Select<ScatterMatrixEigensystem> Dependencies;
     
    static std::string const & name() 
    { 
        static const std::string n("Principal<CoordinateSystem>");
        return n;
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

class Minimum
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("Minimum");
        return n;
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
            detail::reshapeImpl(value_, s, NumericTraits<element_type>::max());
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

class Maximum
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("Maximum");
        return n;
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
            detail::reshapeImpl(value_, s, NumericTraits<element_type>::min());
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

class ArgMinWeight
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("ArgMinWeight");
        return n;
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
            detail::reshapeImpl(value_, s);
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

class ArgMaxWeight
{
  public:
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n("ArgMaxWeight");
        return n;
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
            detail::reshapeImpl(value_, s);
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
        value_ += o.value_;
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
        vigra_precondition(mi < ma,
            "RangeHistogramBase::setMinMax(...): min < max required.");
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

template <int BinCount>
class IntegerHistogram
{
  public:
    
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("IntegerHistogram<") + asString(BinCount) + ">";
        return n;
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

template <int BinCount>
class UserRangeHistogram
{
  public:
    
    typedef Select<> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("UserRangeHistogram<") + asString(BinCount) + ">";
        return n;
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

template <int BinCount>
class AutoRangeHistogram
{
  public:
    
    typedef Select<Minimum, Maximum> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("AutoRangeHistogram<") + asString(BinCount) + ">";
        return n;
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

template <int BinCount>
class GlobalRangeHistogram
{
  public:
    
    typedef Select<Global<Minimum>, Global<Maximum>, Minimum, Maximum> Dependencies;
    
    static std::string const & name() 
    { 
        static const std::string n = std::string("GlobalRangeHistogram<") + asString(BinCount) + ">";
        return n;
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

template <class HistogramAccumulator> 
class StandardQuantiles
{
  public:
    
    typedef typename StandardizeTag<HistogramAccumulator>::type HistogramTag;
    typedef Select<HistogramTag, Minimum, Maximum, Count> Dependencies;

    static std::string const & name() 
    { 
        static const std::string n = std::string("StandardQuantiles<") + HistogramTag::name() + " >";
        return n;
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
                static const double desiredQuantiles[] = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
                getAccumulator<HistogramTag>(*this).computeStandardQuantiles(getDependency<Minimum>(*this), getDependency<Maximum>(*this), 
                                                                             getDependency<Count>(*this), value_type(desiredQuantiles), 
                                                                             this->value_);
                this->setClean();
            }
            return this->value_;
        }
    };
};

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
