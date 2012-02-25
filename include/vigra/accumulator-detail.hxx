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

#ifndef VIGRA_ACCUMULATOR_DETAIL_HXX
#define VIGRA_ACCUMULATOR_DETAIL_HXX

#include "config.hxx"
#include "metaprogramming.hxx"
#include "matrix.hxx"
#include "bit_array.hxx"
#include "static_assert.hxx"
#include <algorithm>
#include <iostream>

namespace vigra {

namespace acc1 {

struct AccumulatorEnd 
{
    typedef AccumulatorEnd Tag;
    typedef void value_type;
    typedef void result_type;
    static const unsigned int workInPass = 0; 
    static const int index = -1;
    
    template <class T>
    void operator()(T const &) {}
    template <class T>
    void operator()(T const &, double) {}
    void operator()() const {}
    
    template <unsigned, class T>
    void pass(T const &) {}
    template <unsigned, class T>
    void pass(T const &, double) {}
    
    template <class T>
    void operator+=(T const &) {}
    
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

    // Select is a synonym for MakeTypeList
template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct Select
: public MakeTypeList<T01, T02, T03, T04, T05, T06, T07, T08, T09, T10,
                      T11, T12, T13, T14, T15, T16, T17, T18, T19, T20>
{};

template <class A>
struct ModifierTraits;

template <class From, class To>
struct DontTransferModifier;

namespace detail {

template <int N1, int N2>
struct Less
{
    static const bool value = N1 < N2;
};

template <class T1, int Priority1, class T2, int Priority2>
struct Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined {};

template <class T1, int Priority, class T2>
struct Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined<T1, Priority, T2, Priority>;

template <class TAG, int BOUND=INT_MIN, bool SKIP=Less<ModifierTraits<TAG>::priority, BOUND>::value, 
          class Contained=typename ModifierTraits<TAG>::ContainedTag>
struct LowestPriority
: public Accumulator___Tag_modifiers_with_same_priority_may_not_be_combined<
                 typename LowestPriority<Contained, BOUND>::type, LowestPriority<Contained, BOUND>::priority,
                 typename ModifierTraits<TAG>::type, ModifierTraits<TAG>::priority>
{
    typedef typename LowestPriority<Contained, BOUND>::type BestContained;
    static const int P1 = LowestPriority<Contained, BOUND>::priority;
    static const int P2 = ModifierTraits<TAG>::priority;
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
    static const int priority = ModifierTraits<TAG>::priority;
};

template <class TAG, int BOUND=INT_MIN, int PRIORITY=LowestPriority<TAG, BOUND>::priority>
struct SortModifiers
{
    typedef ModifierTraits<typename LowestPriority<TAG, BOUND>::type> Traits;
    typedef typename Traits::template rebind<typename SortModifiers<TAG, PRIORITY+1>::type>::type type;
};

template <class TAG, int BOUND>
struct SortModifiers<TAG, BOUND, INT_MAX>
{
    typedef typename LowestPriority<TAG, BOUND>::type type;
};

template <class T1, class T2, 
          class Next=typename ModifierTraits<T1>::ContainedTag,
          class Contained=typename ModifierTraits<T2>::ContainedTag>
struct ContainsModifier
{
    typedef typename IsSameType<typename ModifierTraits<T1>::template rebind<Contained>::type, T2>::type Same;
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

template <class From, class To, class Next=typename ModifierTraits<To>::ContainedTag>
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

template <class From, class To, class Next=typename ModifierTraits<From>::ContainedTag,
         class Dont=typename Or<typename ContainsModifier<From,To>::type,
                                typename ForbiddenTransfer<From, To>::type>::type>
struct TransferModifiersImpl
{
    typedef typename TransferModifiersImpl<Next, To>::type Inner;
    typedef typename ModifierTraits<From>::template rebind<Inner>::type type;
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

template <class From, class To>
struct TransferModifiers
{
    typedef typename TransferModifiersImpl<From, To>::type Unsorted;
    typedef typename SortModifiers<Unsorted>::type type;
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

    // Insert the dependencies of the selected functors into the TypeList and sort
    // the list such that dependencies come after the functors using them. Make sure 
    // that each functor is contained only once.
template <class T>
struct AddDependencies
{
    typedef typename AddDependencies<typename T::type>::type type;
};

template <class HEAD, class TAIL>
struct AddDependencies<TypeList<HEAD, TAIL> >
{
    typedef typename AddDependencies<TAIL>::type                                   TailWithDependencies;
    typedef typename AddDependencies<typename HEAD::Dependencies>::type            HeadDependencies;
    typedef TypeList<HEAD, HeadDependencies>                                       HeadWithDependencies;
    typedef typename PushUnique<HeadWithDependencies, TailWithDependencies>::type  type;
    typedef typename type::Head Head;
    typedef typename type::Tail Tail;
};

template <>
struct AddDependencies<void>
{
    typedef void type;
    typedef void Head;
    typedef void Tail;
};

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
        a.merge(o);
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
            a.merge(o);
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
    
    void operator+=(Decorator const & o)
    {
        DecoratorImpl<Decorator, Decorator::workInPass, Dynamic>::merge(*this, o);
        this->next_ += o.next_;
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
    void operator()(A & a, T const & t, MetaInt<N>)
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
    typedef typename detail::Compose<T, AccumulatorTags, dynamic>::type Accumulators;

    detail::ReshapeImpl<typename detail::NeedsReshape<Accumulators>::type> reshape_;
    
    void reset()
    {
        reshape_.reset();
        Accumulators::reset();
    }

    template <unsigned N>
    void update(T const & t)
    {
        reshape_(*this, t, MetaInt<N>());
        Accumulators::pass<N>(t);
    }
    
    template <unsigned N>
    void update(T const & t, double weight)
    {
        reshape_(*this, t, MetaInt<N>());
        Accumulators::pass<N>(t, weight);
    }

    using Accumulators::operator();
	
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

template <class Tag, class A, class FromTag=typename A::Tag>
struct LookupTag
: public LookupTag<Tag, typename A::BaseType>
{};

template <class Tag, class A, class FromTag>
struct LookupTag<Tag, A const, FromTag>
: public LookupTag<Tag, A>
{
    typedef typename LookupTag<Tag, A>::type const & reference;
};

template <class Tag, class A>
struct LookupTag<Tag, A, Tag>
{
    typedef A type;
    typedef A & reference;
    typedef typename A::value_type value_type;
    typedef typename A::result_type result_type;
};

template <class Tag, class A>
struct LookupTag<Tag, A const, Tag>
: public LookupTag<Tag, A, Tag>
{
    typedef A const & reference;
};

template <class Tag, class A>
struct LookupTag<Tag, A, AccumulatorEnd>
{
    typedef A type;
    typedef A & reference;
    typedef void value_type;
    typedef void result_type;
};

namespace detail {

template <class Tag, class FromTag>
struct CastImpl
{
    template <class A>
    static typename LookupTag<Tag, A>::reference
    cast(A & a)
    {
        return CastImpl<Tag, typename A::BaseType::Tag>::cast(a.next_);
    }
    
    template <class A>
    static typename LookupTag<Tag, A>::result_type
    get(A const & a)
    {
        return CastImpl<Tag, typename A::BaseType::Tag>::get(a.next_);
    }
};

template <class Tag>
struct CastImpl<Tag, Tag>
{
    template <class A>
    static typename LookupTag<Tag, A>::reference
    cast(A & a)
    {
        return a;
    }
    
    template <class A>
    static typename LookupTag<Tag, A>::result_type
    get(A const & a)
    {
        return a();
    }
};

template <class Tag>
struct CastImpl<Tag, AccumulatorEnd>
{
    template <class A>
    static typename LookupTag<Tag, A>::reference
    cast(A & a)
    {
        return a;
    }
    
    template <class A>
    static void
    get(A const & a)
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
    }
};

} // namespace detail

    // cast an accumulator chain to the type specified by Tag
template <class Tag, class A>
typename LookupTag<Tag, A>::reference
cast(A & a)
{
    return detail::CastImpl<Tag, typename A::Tag>::cast(a);
}

    // get the result of the accumulator specified by Tag
template <class Tag, class A>
typename LookupTag<Tag, A>::result_type
get(A const & a)
{
    return detail::CastImpl<Tag, typename A::Tag>::get(a);
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

namespace detail {

#if 0
    // this class closes an accumulator chain
struct AccumulatorBase 
{
	typedef AccumulatorBase Tag;
    typedef void            value_type;
        // 'result_type' can be 'value_type' or 'value_type const &'
        // depending on the return value of operator()()
    typedef void            result_type;
    
    void operator()() const
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access inactive statistic."));
    }
    
    void activate()
    {}
    
    void reset()
    {}
    
    unsigned int passCount() const
    {
        return 1;
    }

    unsigned int passesRequired() const
    {
        return 1;
    }

    template <class Shape>
    void reshape(Shape const &)
    {}

    template <class Shape>
    void resize(Shape const &)
    {}
};

template <class T, unsigned LEVEL=0>
struct TypedAccumulatorBase
: public AccumulatorBase
{
	typedef AccumulatorBase Tag;
    typedef AccumulatorBase BaseType;
    
    static const UInt64 level = LEVEL;

    typedef T const &       argument_type;
    typedef argument_type   first_argument_type;
    typedef double          second_argument_type;
    
        // used for dynamic accumulators to keep track of the active ones.
    BitArray<LEVEL> active_accumulators_;

    void reset()
    {
        active_accumulators_.clear();
    }
    
	void operator+=(TypedAccumulatorBase const &)
    {}
    
	void pass1(argument_type)
    {}
    
    void pass1(first_argument_type, second_argument_type)
    {}
    
    void pass2(argument_type)
    {}
    
    void pass2(first_argument_type, second_argument_type)
    {}
    
	void updatePass1(argument_type)
    {}
    
    void updatePass1(first_argument_type, second_argument_type)
    {}
    
    void updatePass2(argument_type)
    {}
    
    void updatePass2(first_argument_type, second_argument_type)
    {}
};

// class Count;

// template <class TAG = void>
// class Coord;

// template <class Accumulators>
// struct CoordDependencies
// {
    // typedef typename Push<typename Coord<typename Accumulators::Head>::type, 
                          // typename CoordDependencies<typename Accumulators::Tail>::type>::type type;
// };

// template <>
// struct CoordDependencies<void>
// {
    // typedef void type;
// };

// template <class TAG>
// class Coord
// {
  // public:
  
    // typedef typename CoordDependencies<typename TAG::Dependencies::type>::type Dependencies;
    // typedef Coord<TAG> type;
    
    // template<class U, class Base>
    // struct Impl
    // : public TAG::template Impl<typename U::CoordType, Base>
    // {
        // typedef Coord<TAG> Tag;
        // typedef Base BaseType;
        // typedef typename TAG::template Impl<typename U::CoordType, Base> WrappedType;
    // };
// };

// template <>
// class Coord<Count>
// {
  // public:
    // typedef Count type;
// };

// template <>
// class Coord<void>
// {
  // public:
    // typedef void type;
// };

// struct Coupled
// {
    // typedef Shape2 CoordType;
    // typedef double DataType;
// };

// template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          // class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          // class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          // class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
// class UseCoord
// : public MakeTypeList<typename Coord<T01>::type, typename Coord<T02>::type, typename Coord<T03>::type, 
                      // typename Coord<T04>::type, typename Coord<T05>::type, typename Coord<T06>::type, 
                      // typename Coord<T07>::type, typename Coord<T08>::type, typename Coord<T09>::type, 
                      // typename Coord<T10>::type, typename Coord<T11>::type, typename Coord<T12>::type, 
                      // typename Coord<T13>::type, typename Coord<T14>::type, typename Coord<T15>::type, 
                      // typename Coord<T16>::type, typename Coord<T17>::type, typename Coord<T18>::type, 
                      // typename Coord<T19>::type, typename Coord<T20>::type >
// {};

// template <class T>
// struct Expand
// {
    // typedef T type;
// };

// template <class T01, class T02, class T03, class T04, class T05,
          // class T06, class T07, class T08, class T09, class T10,
          // class T11, class T12, class T13, class T14, class T15,
          // class T16, class T17, class T18, class T19, class T20>
// struct Expand<UseCoord<T01, T02, T03, T04, T05, T06, T07, T08, T09, T10,
                  // T11, T12, T13, T14, T15, T16, T17, T18, T19, T20> >
// : public UseCoord<T01, T02, T03, T04, T05, T06, T07, T08, T09, T10,
                  // T11, T12, T13, T14, T15, T16, T17, T18, T19, T20>
// {};

    // // Select is a synonym for MakeTypeList
// template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          // class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          // class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          // class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
// struct Select
// : public MakeTypeList<typename Expand<T01>::type, typename Expand<T02>::type, typename Expand<T03>::type, 
                      // typename Expand<T04>::type, typename Expand<T05>::type, typename Expand<T06>::type, 
                      // typename Expand<T07>::type, typename Expand<T08>::type, typename Expand<T09>::type, 
                      // typename Expand<T10>::type, typename Expand<T11>::type, typename Expand<T12>::type, 
                      // typename Expand<T13>::type, typename Expand<T14>::type, typename Expand<T15>::type, 
                      // typename Expand<T16>::type, typename Expand<T17>::type, typename Expand<T18>::type, 
                      // typename Expand<T19>::type, typename Expand<T20>::type >
// {};

template <class TAG>
struct Dynamic;

template <class T, class Wrapped>
struct ReshapeConfig
{
    template <class A, class Shape>
    static void reshape(A & a, Shape const & s)
    {}
};

template <unsigned int N, class T, class Alloc, class Wrapped>
struct ReshapeConfig<MultiArray<N, T, Alloc>, Wrapped>
{
    template <class A, class Shape>
    static void reshape(A & a, Shape const & s)
    {
        ((Wrapped &)a).reshape(s);
    }
};

template <class T, class Alloc, class Wrapped>
struct ReshapeConfig<Matrix<T, Alloc>, Wrapped>
{
    template <class A, class Shape>
    static void reshape(A & a, Shape const & s)
    {
        ((Wrapped &)a).reshape(s);
    }
};

template <class T, class TAG, class BaseBase, unsigned LEVEL>
struct Config
: public ReshapeConfig<typename TAG::template Impl<T, BaseBase>::value_type, typename TAG::template Impl<T, BaseBase> >
{
    typedef TAG Tag;
    typedef typename TAG::template Impl<T, BaseBase> Wrapped;
    typedef typename BaseBase BaseType;
    typedef T const & argument_type;
    typedef typename Wrapped::result_type result_type;
    static const unsigned level = LEVEL;
    
    template <class A>
    static result_type get(A const & a)
    {
        return ((Wrapped &)a)();
    }
    
    template <class A>
    static void reset(A & a)
    {
        ((Wrapped &)a).reset();
    }
   
    template <class A>
    static unsigned int passCount(A const & a)
    {
        return ((Wrapped const &)a).passesRequired();
    }
        
    template <class A>
    static void merge(A & a, A const & o)
    {
        ((Wrapped &)a).merge(o);
    }
    
    template <class A>
	static void pass1(A & a, argument_type t)
    {
        ((Wrapped &)a).updatePass1(t);
    }
    
    template <class A>
    static void pass1(A & a, argument_type t, double weight)
    {
        ((Wrapped &)a).updatePass1(t, weight);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t)
    {
        ((Wrapped &)a).updatePass2(t);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t, double weight)
    {
        ((Wrapped &)a).updatePass2(t, weight);
    }
};

template <class CONFIG>
struct RecurseBefore
: public CONFIG
{
    typedef CONFIG Base;
    typedef typename CONFIG::BaseType BaseBase;
    typedef typename CONFIG::argument_type argument_type;
    
    template <class A>
    static void activate(A & a)
    {
        Base::activate(a);
        ((BaseBase &)a).activate();
    }
    
    template <class A>
    static void reset(A & a)
    {
        Base::reset(a);
        ((BaseBase &)a).reset();
    }
   
    template <class A>
    static unsigned int passCount(A const & a)
    {
        unsigned int p1 = Base::passCount(a);
        unsigned int p2 = ((BaseBase &)a).passCount();
        return std::max(p1, p2);
    }
        
    template <class A, class Shape>
    static void reshape(A & a, Shape const & s)
    {
        Base::reshape(a, s);
        ((BaseBase &)a).resize(s);
    }
   
    template <class A>
    static void merge(A & a, A const & o)
    {
        Base::merge(a, o);
        ((BaseBase &)a).operator+=(o);
    }
    
    template <class A>
	static void pass1(A & a, argument_type t)
    {
        ((BaseBase &)a).pass1(t);
        Base::pass1(a, t);
    }
    
    template <class A>
    static void pass1(A & a, argument_type t, double weight)
    {
        ((BaseBase &)a).pass1(t, weight);
        Base::pass1(a, t, weight);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t)
    {
        ((BaseBase &)a).pass2(t);
        Base::pass2(a, t);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t, double weight)
    {
        ((BaseBase &)a).pass2(t, weight);
        Base::pass2(a, t, weight);
    }
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

template <class T, class TAG, class BaseBase, unsigned LEVEL>
struct Config<T, Dynamic<TAG>, BaseBase, LEVEL>
: public ReshapeConfig<typename TAG::template Impl<T, BaseBase>::value_type, typename TAG::template Impl<T, BaseBase> >
{
    typedef TAG Tag;
    typedef typename TAG::template Impl<T, BaseBase> Wrapped;
    typedef typename BaseBase BaseType;
    typedef T const & argument_type;
    typedef typename Wrapped::result_type result_type;
    static const unsigned level = LEVEL;
    
    template <class A>
    static result_type get(A const & a)
    {
        vigra_precondition(a.active_accumulators_.test<LEVEL>(),
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
        return ((Wrapped &)a)();
    }
    
    template <class A>
    static void activate(A & a)
    {
        a.active_accumulators_.set<LEVEL>();
        ActivateDependencies<typename Tag::Dependencies::type>::exec(a);
    }
    
    template <class A>
    static void reset(A & a)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).reset();
    }
   
    template <class A>
    static unsigned int passCount(A const & a)
    {
        if(a.active_accumulators_.test<LEVEL>())
            return ((Wrapped const &)a).passesRequired();
        else
            return 0;
    }
        
    template <class A, class Shape>
    static void reshape(A & a, Shape const & s)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ReshapeConfig<typename TAG::template Impl<T, BaseBase>::value_type>::reshape(a, s);
    }
   
    template <class A>
    static void merge(A & a, A const & o)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).merge(o);
    }
    
    template <class A>
	static void pass1(A & a, argument_type t)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).updatePass1(t);
    }
    
    template <class A>
    static void pass1(A & a, argument_type t, double weight)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).updatePass1(t, weight);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).updatePass2(t);
    }
    
    template <class A>
    static void pass2(A & a, argument_type t, double weight)
    {
        if(a.active_accumulators_.test<LEVEL>())
            ((Wrapped &)a).updatePass2(t, weight);
    }
};

template <class CONFIG>
struct Decorator
: public CONFIG::Wrapped
{
    typedef typename CONFIG::Tag Tag;
    typedef typename CONFIG::BaseType BaseType;
    typedef typename CONFIG::argument_type argument_type;
    typedef typename CONFIG::result_type result_type;
    
    void activate()
    {
        CONFIG::activate(*this);
    }
    
    void reset()
    {
        CONFIG::reset(*this);
    }
   
    unsigned int passCount() const
    {
        return CONFIG::passCount(*this);
    }
        
    template <class Shape>
    void resize(Shape const & s)
    {
        CONFIG::reshape(*this, s);
    }
    
    result_type operator()() const
    {
        return CONFIG::get(*this);
    }
   
    void operator+=(Decorator const & o)
    {
        CONFIG::merge(*this, o);
    }
    
	void pass1(argument_type t)
    {
        CONFIG::pass1(*this, t);
    }
    
    void pass1(argument_type t, double weight)
    {
        CONFIG::pass1(*this, t, weight);
    }
    
    void pass2(argument_type t)
    {
        CONFIG::pass2(*this, t);
    }
    
    void pass2(argument_type t, double weight)
    {
        CONFIG::pass2(*this, t, weight);
    }
    
        // make sure that children do not inherit active versions of these functions
    void updatePass1(argument_type) {}
    void updatePass1(argument_type, double) {}
    void updatePass2(argument_type) {}
    void updatePass2(argument_type, double) {}
    void merge(Decorator const &) {}
    template <class Shape> void reshape(Shape const & s) {}
};
#if 0
template <class T, class TAG, class BASE, unsigned LEVEL>
struct SelectConfig
{
    typedef RecurseBefore<Config<T, TAG, BASE, LEVEL> > type;
};

template <class T, class TAG, class BASE, unsigned LEVEL>
struct SelectDecorator<T, Dynamic<TAG>, BASE, LEVEL>
{
    typedef DynamicAccumulatorDecorator<T, TAG, BASE, LEVEL> Decorator;
    typedef RecurseBefore<T, Decorator, BASE, LEVEL> type;
};
#endif

    // helper classes to create an accumulator chain from a TypeList
template <class T, class Accumulators, unsigned level=0>
struct Compose
{
    typedef typename Accumulators::Head Tag; 
    typedef typename Compose<T, typename Accumulators::Tail, level+1>::type BaseType;
    typedef Decorator<RecurseBefore<Config<T, Tag, BaseType, level> > > type;
};

template <class T, unsigned level> 
struct Compose<T, void, level> 
{ 
    typedef TypedAccumulatorBase<T, level> type; 
};

    // helper classes to create a dynamic accumulator chain from a TypeList
template <class T, class Accumulators, unsigned level=0>
struct DynamicCompose
{
    typedef typename Accumulators::Head Tag; 
    typedef typename DynamicCompose<T, typename Accumulators::Tail, level+1>::type BaseType;
    typedef Decorator<RecurseBefore<Config<T, Dynamic<Tag>, BaseType, level> > > type;
};

template <class T, unsigned level> 
struct DynamicCompose<T, void, level> 
{ 
    typedef TypedAccumulatorBase<T, level> type; 
};

#endif

} // namespace detail

#if 0
    // create an accumulator chain containing the Selected statistics and their dependencies
template <class T, class Selected>
struct Accumulator
: public detail::ReshapeImpl<T, 
           typename detail::Compose<T, typename detail::AddDependencies<typename Selected::type>::type>::type>
{
    typedef typename detail::AddDependencies<typename Selected::type>::type Accumulators;
    typedef typename detail::Compose<T, Accumulators>::type BaseType;
    typedef VigraFalseType Tag;

	void updatePass2(T const & t)
    {
        BaseType::pass2(t);
    }
    
    void updatePass2(T const & t, double weight)
    {
        BaseType::pass2(t, weight);
    }
    
    unsigned int passesRequired() const
    {
        return BaseType::passCount();
    }
};

    // create a dynamic accumulator chain containing the Selected statistics and their dependencies.
    // Statistics will only be computed if activate<Tag>() is called at runtime.
template <class T, class Selected>
struct DynamicAccumulator
: public detail::ReshapeImpl<T, 
           typename detail::DynamicCompose<T, typename detail::AddDependencies<typename Selected::type>::type>::type>
{
    typedef typename detail::AddDependencies<typename Selected::type>::type Accumulators;
    typedef typename detail::DynamicCompose<T, Accumulators>::type BaseType;
    typedef VigraFalseType Tag;

	void updatePass2(T const & t)
    {
        BaseType::pass2(t);
    }
    
    void updatePass2(T const & t, double weight)
    {
        BaseType::pass2(t, weight);
    }
    
    unsigned int passesRequired() const
    {
        return BaseType::passCount();
    }
};
#endif

#if 0
    // find the type in an accumulator chain holding the given Tag
template <class Tag, class Accumulator>
struct LookupTag
{
    typedef typename IsSameType<Tag, typename Accumulator::Tag>::type Found;
    typedef LookupTag<Tag, typename Accumulator::BaseType> Base;
    typedef typename If<Found, Accumulator, typename Base::type>::type type;
    typedef typename If<Found, typename Accumulator::value_type, 
                               typename Base::value_type>::type value_type;
    typedef typename If<Found, typename Accumulator::result_type, 
                               typename Base::result_type>::type result_type;
};

template <class Tag>
struct LookupTag<Tag, AccumulatorEnd>
{
    typedef AccumulatorEnd type;
    typedef void value_type;
    typedef void result_type;
};
#endif

#if 0
Accumulator<StridedPairPointer<...>, Bind<Data, Select<...> >,    // data statistics
                                     Bind<Coord, Select<...> >,    // coordinate statstics
                                     Bind<Coord, Data, Select<...> >  // weighted coordinate statistics
CoupledAccumulator<StridedTriplePointer<...>, Bind<1, Select<...> >,    // data statistics
                                            Bind<0, Select<...> >,    // coordinate statstics
                                            Bind<0, 1, Select<...> >  // weighted coordinate statistics
                                            Bind<1, 2, Select<...> >  // weighted data statistics
#endif

#if 0
template <unsigned N, class T>
struct CoupledPointer
{
    typedef typename MultiArrayShape<N>::type Coord;
    typedef T value_type;
    
    Coord coord_;
    value_type value_;
    
    CoupledPointer(Coord const & c, value_type v)
    : coord_(c),
      value_(v)
    {}
};

template <class Handle, unsigned index>
struct CoupledHandleTraits
{
    typedef void type;
};

template <unsigned N, class T>
struct CoupledHandleTraits<CoupledPointer<N, T>, 0>
{
    typedef typename MultiArrayShape<N>::type type;
    static type exec(CoupledPointer<N, T> const & c)
    {
        return c.coord_;
    }
};

template <unsigned N, class T>
struct CoupledHandleTraits<CoupledPointer<N, T>, 1>
{
    typedef T type;
    static type exec(CoupledPointer<N, T> const & c)
    {
        return c.value_;
    }
};

template <unsigned INDEX, unsigned N, class T>
typename CoupledHandleTraits<CoupledPointer<N, T>, INDEX>::type
get(CoupledPointer<N, T> const & c)
{
    return CoupledHandleTraits<CoupledPointer<N, T>, INDEX>::exec(c);
}

template <unsigned INDEX, class Selected>
struct Simple
{
    static const unsigned index = INDEX;
    typedef VigraFalseType IsWeighted;
    typedef Selected Accumulators;
};

template <unsigned INDEX, unsigned WEIGHT_INDEX, class Selected>
struct Weighted
{
    static const unsigned index = INDEX;
    static const unsigned weight_index = WEIGHT_INDEX;
    typedef VigraTrueType IsWeighted;
    typedef Selected Accumulators;
};

template <class A, class IsWeighted=typename A::IsWeighted>
struct CallSelector
{
    template <class Accumulator, class CoupledHandle>
    static void pass1(Accumulator & a, CoupledHandle const & h)
    {
        a(get<A::index>(h));
    }
    template <class Accumulator, class CoupledHandle>
    static void pass2(Accumulator & a, CoupledHandle const & h)
    {
        a.updatePass2(get<A::index>(h));
    }
};

template <class A>
struct CallSelector<A, VigraTrueType>
{
    template <class Accumulator, class CoupledHandle>
    static void pass1(Accumulator & a, CoupledHandle const & h)
    {
        a(get<A::index>(h), get<A::weight_index>(h));
    }
    template <class Accumulator, class CoupledHandle>
    static void pass2(Accumulator & a, CoupledHandle const & h)
    {
        a.updatePass2(get<A::index>(h), get<A::weight_index>(h));
    }
};

template <class CoupledHandle, class Selected, int LEVEL>
struct ComposeCoupledAccumulator;

template <class CoupledHandle, class Head, class Tail, int LEVEL>
struct ComposeCoupledAccumulator<CoupledHandle, TypeList<Head, Tail>, LEVEL>
: public ComposeCoupledAccumulator<CoupledHandle, Tail, LEVEL+1>
{
    typedef ComposeCoupledAccumulator<CoupledHandle, Tail, LEVEL+1> BaseType;
    
    static const unsigned index = Head::index;
    static const int level = LEVEL;
    typedef typename CoupledHandleTraits<CoupledHandle, index>::type DataType;
    typedef Accumulator<DataType, typename Head::Accumulators> AccumulatorType;
    
    AccumulatorType value_;
    
    unsigned int passesRequired() const
    {
        return std::max(BaseType::passesRequired(), value_.passesRequired());
    }
    
    void operator+=(ComposeCoupledAccumulator const & o)
    {
        value_ += o.value_;
        BaseType::operator+=(o);
    }
    
    void operator()(CoupledHandle const & h)
    {
        CallSelector<Head>::pass1(value_, h);
        BaseType::operator()(h);
    }
    
    void updatePass2(CoupledHandle const & h)
    {
        CallSelector<Head>::pass2(value_, h);
        BaseType::updatePass2(h);
    }
    
    void reset()
    {
        value_.reset();
        BaseType::reset();
    }
};

template <class CoupledHandle, int LEVEL>
struct ComposeCoupledAccumulator<CoupledHandle, void, LEVEL>
{
    unsigned int passesRequired() const
    {
        return 0;
    }
    
    void operator+=(ComposeCoupledAccumulator const &)
    {}
    
    void operator()(CoupledHandle const &)
    {}
    
    void updatePass2(CoupledHandle const &)
    {}
    
    void reset()
    {}
};

template <class CoupledHandle, class A1=void, class A2=void, class A3=void, class A4=void, class A5=void>
struct CoupledAccumulator
: public ComposeCoupledAccumulator<CoupledHandle, typename MakeTypeList<A1, A2, A3, A4, A5>::type, 0>
{
    static const int level = -1;
    typedef ComposeCoupledAccumulator<CoupledHandle, typename MakeTypeList<A1, A2, A3, A4, A5>::type, 0> BaseType;
};

// template<int LEVEL, class Tag, class CoupledHandle, class Selected>
// typename LookupTag<Tag, typename ComposeCoupledAccumulator<CoupledHandle, Selected, LEVEL>::AccumulatorType>::type const &
// cast(ComposeCoupledAccumulator<CoupledHandle, Selected, LEVEL> const & a)
// {
    // return cast<Tag>(a.value_);
// }

template <int LEVEL, class Accumulator>
struct LookupLevel
{
    typedef typename IfBool<(LEVEL == Accumulator::level), 
               Accumulator, typename LookupLevel<LEVEL, typename Accumulator::BaseType>::type>::type type;
};

template <int LEVEL, class T, int ACTUAL_LEVEL>
struct LookupLevel<LEVEL, ComposeCoupledAccumulator<T,void,ACTUAL_LEVEL> >
{
    typedef void type;
};

template<int LEVEL, class Tag, class CoupledHandle, class A1, class A2, class A3, class A4, class A5>
typename LookupTag<Tag, typename LookupLevel<LEVEL, 
         CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> >::type::AccumulatorType>::type const &
cast(CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> const & a)
{
    typedef typename LookupLevel<LEVEL, CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> >::type LevelType;
    return cast<Tag>(((LevelType const &)a).value_);
}

template<int LEVEL, class Tag, class CoupledHandle, class A1, class A2, class A3, class A4, class A5>
typename LookupTag<Tag, typename LookupLevel<LEVEL, 
         CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> >::type::AccumulatorType>::result_type
get(CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> const & a)
{
    typedef typename LookupLevel<LEVEL, CoupledAccumulator<CoupledHandle, A1, A2, A3, A4, A5> >::type LevelType;
    return get<Tag>(((LevelType const &)a).value_);
}

#endif

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_DETAIL_HXX
