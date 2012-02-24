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
#include <algorithm>
#include <iostream>

namespace vigra {

namespace acc1 {

struct None 
{
    typedef None Tag;
    typedef void value_type;
    static const unsigned int workInPass = 0; 
};

template <class T, class TAG, class NEXT=None>
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

template <class Tag, class A, class FromTag=typename A::Tag>
struct LookupTag
: public LookupTag<Tag, typename A::BaseType>
{};

template <class Tag, class A, class FromTag>
struct LookupTag<Tag, A const, FromTag>
: public LookupTag<Tag, typename A::BaseType const>
{};

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
{
    typedef A const type;
    typedef A const & reference;
    typedef typename A::value_type value_type;
    typedef typename A::result_type result_type;
};

template <class Tag, class FromTag>
struct LookupTag<Tag, None, FromTag>
{
    typedef None type;
    typedef None & reference;
    typedef void value_type;
    typedef void result_type;
};

template <class Tag, class FromTag>
struct LookupTag<Tag, None const, FromTag>
{
    typedef None const type;
    typedef None const & reference;
    typedef void value_type;
    typedef void result_type;
};

namespace detail {

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

template <class A, unsigned CurrentPass, unsigned WorkPass=A::workInPass>
struct UpdateImpl
{
    template <class T>
    static void exec(A & a, T const & t)
    {
        UpdateImpl<typename A::BaseType, CurrentPass>::exec(a.next_, t);
    }

    template <class T>
    static void exec(A & a, T const & t, double weight)
    {
        UpdateImpl<typename A::BaseType>::exec(a.next_, t, weight);
    }
};

template <class A, unsigned CurrentPass>
struct UpdateImpl<A, CurrentPass, CurrentPass>
{
    template <class T>
    static void exec(A & a, T const & t)
    {
        UpdateImpl<typename A::BaseType, CurrentPass>::exec(a.next_, t);
        a.update(t);
    }

    template <class T>
    static void exec(A & a, T const & t, double weight)
    {
        UpdateImpl<typename A::BaseType, CurrentPass>::exec(a.next_, t, weight);
        a.update(t, weight);
    }
};

template <unsigned CurrentPass, unsigned WorkPass>
struct UpdateImpl<None, CurrentPass, WorkPass>
{
    template <class T>
    static void exec(None &, T const &)
    {}

    template <class T>
    static void exec(None &, T const &, double)
    {}
};

template <class A>
struct MergeImpl
{
    static void exec(A & a, A const & o)
    {
        a.merge(o);
        MergeImpl<typename A::BaseType>::exec(a.next_, o.next_);
    }
    static void exec(None & a, None const & o)
    {}
};

template <>
struct MergeImpl<None>
{
    static void exec(None & a, None const & o)
    {}
};

template <class A>
struct ReshapeImpl
{
    template <class Shape>
    static void exec(A & a, Shape const & s)
    {
        ReshapeImpl<typename A::BaseType>::exec(a.next_, s);
        a.reshape(s);
    }
};

template <>
struct ReshapeImpl<None>
{
    template <class Shape>
    static void exec(None & a, Shape const & s)
    {}
};

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
    
    static None const & cast(None const & a)
    {
        return a;
    }
    static void get(None const & a)
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
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
struct NeedsReshape<None, None::value_type>
{
    typedef VigraFalseType type;
};

    // This accumulator is inserted on top of an accumulator chain to call reshape() 
    // when the first data item arrives. This is necessary if the shape of the result 
    // depends on the shape of the input and cannot be determined at compile time. 
    // The above NeedsReshape traits specify the types where this applies. If the chain
    // doesn't contain such types, ReshapeHelper will expand into a do-nothing version.
template <class T, class BASE, class NeedsReshape=typename NeedsReshape<BASE>::type>
struct ReshapeHelper
: public BASE
{
    // typedef VigraFalseType Tag;
    // typedef BASE           BaseType;
    
    bool needs_reshape_;
    
    ReshapeHelper()
    : needs_reshape_(true)
    {}
    
    using BASE::operator();

	void operator()(T const & t)
    {
        if(needs_reshape_)
        {
            detail::ReshapeImpl<BASE>::exec(*this, shape(t));
            needs_reshape_ = false;
        }
        detail::UpdateImpl<BASE, 1>::exec(*this, t);
    }
    
    void operator()(T const & t, double weight)
    {
        if(needs_reshape_)
        {
            detail::ReshapeImpl<BASE>::exec(*this, shape(t));
            needs_reshape_ = false;
        }
        detail::UpdateImpl<BASE, 1>::exec(*this, t, weight);
    }
    
    template <unsigned int N, class U, class Stride>
    typename MultiArrayShape<N>::type
    shape(MultiArrayView<N, U, Stride> const & a)
    {
        return a.shape();
    }
    
    template <class U, int N>
    Shape1
    shape(TinyVector<U, N> const & a)
    {
        return Shape1(N);
    }
};

template <class T, class BASE>
struct ReshapeHelper<T, BASE, VigraFalseType>
: public BASE
{
    // typedef VigraFalseType Tag;
    // typedef BASE           BaseType;

    using BASE::operator();

	void operator()(T const & t)
    {
        detail::UpdateImpl<BASE, 1>::exec(*this, t);
    }
    
    void operator()(T const & t, double weight)
    {
        detail::UpdateImpl<BASE, 1>::exec(*this, t, weight);
    }
};

    // helper classes to create an accumulator chain from a TypeList
template <class T, class Accumulators, unsigned level=0>
struct Compose
{
    typedef typename Accumulators::Head Tag; 
    typedef Compose<T, typename Accumulators::Tail, level+1> ComposeType;
    typedef typename ComposeType::type BaseType;
    typedef typename Tag::template Impl<T, AccumulatorBase<T, Tag, BaseType> >  type;
    static const unsigned passes = ComposeType::passes < type::workInPass ? type::workInPass : ComposeType::passes;
};

template <class T, unsigned level> 
struct Compose<T, void, level> 
{ 
    typedef None type;
    static const unsigned passes = 0;
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

    // create an accumulator chain containing the Selected statistics and their dependencies
template <class T, class Selected>
struct Accumulator
: public detail::ReshapeHelper<T, typename detail::Compose<T, typename detail::AddDependencies<typename Selected::type>::type>::type>
{
    typedef typename detail::AddDependencies<typename Selected::type>::type AccumulatorTags;
    typedef detail::Compose<T, AccumulatorTags> ComposeAccumulators;
    typedef typename ComposeAccumulators::type Accumulators;
//    typedef VigraFalseType Tag;

    void operator+=(Accumulator const & t)
    {
        detail::MergeImpl<Accumulators>::exec(*this, t);
    }
    
    void updatePass2(T const & t)
    {
        detail::UpdateImpl<Accumulators, 2>::exec(*this, t);
    }
    
    void updatePass2(T const & t, double weight)
    {
        detail::UpdateImpl<Accumulators, 2>::exec(*this, t, weight);
    }
    
    unsigned int passesRequired() const
    {
        return ComposeAccumulators::passes;
    }
};

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
struct Wrapper
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
   
    void operator+=(Wrapper const & o)
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
    void merge(Wrapper const &) {}
    template <class Shape> void reshape(Shape const & s) {}
};
#if 0
template <class T, class TAG, class BASE, unsigned LEVEL>
struct SelectConfig
{
    typedef RecurseBefore<Config<T, TAG, BASE, LEVEL> > type;
};

template <class T, class TAG, class BASE, unsigned LEVEL>
struct SelectWrapper<T, Dynamic<TAG>, BASE, LEVEL>
{
    typedef DynamicAccumulatorWrapper<T, TAG, BASE, LEVEL> Wrapper;
    typedef RecurseBefore<T, Wrapper, BASE, LEVEL> type;
};
#endif

    // helper classes to create an accumulator chain from a TypeList
template <class T, class Accumulators, unsigned level=0>
struct Compose
{
    typedef typename Accumulators::Head Tag; 
    typedef typename Compose<T, typename Accumulators::Tail, level+1>::type BaseType;
    typedef Wrapper<RecurseBefore<Config<T, Tag, BaseType, level> > > type;
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
    typedef Wrapper<RecurseBefore<Config<T, Dynamic<Tag>, BaseType, level> > > type;
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
: public detail::ReshapeHelper<T, 
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
: public detail::ReshapeHelper<T, 
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
struct LookupTag<Tag, None>
{
    typedef None type;
    typedef void value_type;
    typedef void result_type;
};
#endif

    // activate the dynamic accumulator specified by Tag
template <class Tag, class Accumulator>
void activate(Accumulator & a)
{
    cast<Tag>(a).activate();
}

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
