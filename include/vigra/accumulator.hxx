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

#include "config.hxx"
#include "metaprogramming.hxx"
#include "mathutil.hxx"
#include <algorithm>
#include <iostream>

namespace vigra {

namespace acc1 {

template<class HEAD, class TAIL=void>
struct AccumulatorList
{
	typedef AccumulatorList<HEAD, TAIL> type;
    typedef HEAD Head;
    typedef TAIL Tail;
};

template <class List, class T>
struct Contains
{
    typedef typename IsSameType<typename List::Head, T>::type Found;
    typedef typename If<Found, VigraTrueType, typename Contains<typename List::Tail, T>::type>::type type;
};

template <class T>
struct Contains<void, T>
{
    typedef VigraFalseType type;
};

template <class A, class Tail=void>
struct Push
{
    typedef AccumulatorList<A, typename Tail::type> type;
};

template <class A>
struct Push<A, void>
{
    typedef AccumulatorList<A, void> type;
};

template <>
struct Push<void, void>
{
    typedef void type;
};

template <class List1, class List2>
struct Merge;

template <class Head, class Tail, class List>
struct Merge<AccumulatorList<Head, Tail>, List>
{
    typedef typename Merge<Tail, List>::type Rest;
    typedef typename Contains<Rest, Head>::type HeadAlreadyInList;
    typedef typename If<HeadAlreadyInList, Rest, AccumulatorList<Head, Rest> >::type type;
};

template <class List>
struct Merge<void, List>
{
    typedef List type;
};

template <class List>
struct PushDependencies;

template <class HEAD, class TAIL>
struct PushDependencies<AccumulatorList<HEAD, TAIL> >
{
    typedef typename PushDependencies<TAIL>::type TailWithDependencies;
    typedef typename HEAD::Dependencies Dependencies; // must be of type Select<...>
    typedef typename PushDependencies<typename Dependencies::type>::type HeadDependencies;
    typedef AccumulatorList<HEAD, HeadDependencies> HeadWithDependencies;
    typedef typename Merge<HeadWithDependencies, TailWithDependencies>::type type;
    typedef typename type::Head Head;
    typedef typename type::Tail Tail;
};

template <>
struct PushDependencies<void>
{
    typedef void type;
    typedef void Head;
    typedef void Tail;
};

struct AccumulatorBase 
{
	typedef AccumulatorBase Tag;
    typedef void result_type;
    typedef double second_argument_type;

    template <class T>
	void operator()(T const &)
    {}
    
    template <class T>
    void operator()(T const &, double)
    {}
    
    void operator()() const
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access statistic '") << typeid(Tag).name() << "'.");
    }
    
    void activate()
    {}
};

template <class Tag, class Accumulator>
struct LookupTag
{
    typedef typename IsSameType<Tag, typename Accumulator::Tag>::type Found;
    typedef typename If<Found, Accumulator, typename LookupTag<Tag, typename Accumulator::BaseType>::type>::type type;
    typedef typename If<Found, typename Accumulator::result_type, 
                               typename LookupTag<Tag, typename Accumulator::BaseType>::result_type>::type result_type;
};

template <class Tag>
struct LookupTag<Tag, AccumulatorBase>
{
    typedef AccumulatorBase type;
    typedef void result_type;
};

template <class Tag, class Accumulator>
typename LookupTag<Tag, Accumulator>::type &
cast(Accumulator & a)
{
    return static_cast<typename LookupTag<Tag, Accumulator>::type &>(a);
}

template <class Tag, class Accumulator>
typename LookupTag<Tag, Accumulator>::type const &
cast(Accumulator const & a)
{
    return static_cast<typename LookupTag<Tag, Accumulator>::type const &>(a);
}

template <class Tag, class Accumulator>
typename LookupTag<Tag, Accumulator>::result_type
getImpl(Accumulator const & a)
{
    return a();
}

template <class Tag>
void getImpl(AccumulatorBase const & a)
{
    vigra_precondition(false,
        std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
}

template <class Tag, class Accumulator>
typename LookupTag<Tag, Accumulator>::result_type
get(Accumulator const & a)
{
    return getImpl<Tag>(cast<Tag>(a));
}

template <class Tag, class Accumulator>
void activate(Accumulator & a)
{
    cast<Tag>(a).activate();
}

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

template <class T, class TAG, class BaseBase>
struct DynamicAccumulatorWrapper
: public TAG::template Impl<T, BaseBase>
{
    typedef TAG Tag;
    typedef typename TAG::template Impl<T, BaseBase> Base;
    
    typedef T argument_type;
    typedef T first_argument_type;
    typedef typename Base::result_type result_type;

    bool is_active_;
    
    DynamicAccumulatorWrapper()
    : is_active_(false)
    {}
    
    void activate()
    {
        is_active_ = true;
        ActivateDependencies<typename Tag::Dependencies::type>::exec(*this);
    }
    
	result_type operator()() const
    {
        vigra_precondition(is_active_,
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
        return Base::operator()();
    }
    
    template <class T>
	void operator()(T const & t)
    {
        if(is_active_)
        {
            Base::operator()(t);
        }
        else
        {
            BaseBase::operator()(t);
        }
    }
    
    template <class T>
    void operator()(T const & t, double weight)
    {
        if(is_active_)
        {
            Base::operator()(t, weight);
        }
        else
        {
            BaseBase::operator()(t, weight);
        }
    }
};

template <class T, class Tag, class Accumulators>
struct Compose
{
    typedef typename Compose<T, typename Accumulators::Head, typename Accumulators::Tail>::type BaseType;
    typedef typename Tag::template Impl<T, BaseType> type;
};

template <class T, class Tag> 
struct Compose<T, Tag, void> 
{ 
    typedef AccumulatorBase BaseType;
    typedef typename Tag::template Impl<T, AccumulatorBase> type; 
};

template <class T, class Tag, class Accumulators>
struct DynamicCompose
{
    typedef typename DynamicCompose<T, typename Accumulators::Head, typename Accumulators::Tail>::type BaseType;
    typedef DynamicAccumulatorWrapper<T, Tag, BaseType> type;
};

template <class T, class Tag> 
struct DynamicCompose<T, Tag, void> 
{ 
    typedef AccumulatorBase BaseType;
    typedef DynamicAccumulatorWrapper<T, Tag, AccumulatorBase> type; 
};

template <class T, class Selected>
struct Accumulator
: public Compose<T, typename PushDependencies<typename Selected::type>::Head,
                    typename PushDependencies<typename Selected::type>::Tail>::type
{
    typedef typename PushDependencies<typename Selected::type>::type Accumulators;
    typedef typename Compose<T, typename Accumulators::Head, typename Accumulators::Tail>::type BaseType;
    typedef VigraFalseType Tag;
        
    using BaseType::operator();
};

template <class T, class Selected>
struct DynamicAccumulator
: public DynamicCompose<T, typename PushDependencies<typename Selected::type>::Head,
                           typename PushDependencies<typename Selected::type>::Tail>::type
{
    typedef typename PushDependencies<typename Selected::type>::type Accumulators;
    typedef typename DynamicCompose<T, typename Accumulators::Head, typename Accumulators::Tail>::type BaseType;
    typedef VigraFalseType Tag;
        
    using BaseType::operator();
};

template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct Select
{
    typedef typename Push<T19, T20>::type L19;
    typedef typename Push<T18, L19>::type L18;
    typedef typename Push<T17, L18>::type L17;
    typedef typename Push<T16, L17>::type L16;
    typedef typename Push<T15, L16>::type L15;
    typedef typename Push<T14, L15>::type L14;
    typedef typename Push<T13, L14>::type L13;
    typedef typename Push<T12, L13>::type L12;
    typedef typename Push<T11, L13>::type L11;
    typedef typename Push<T10, L11>::type L10;
    typedef typename Push<T09, L10>::type L09;
    typedef typename Push<T08, L09>::type L08;
    typedef typename Push<T07, L08>::type L07;
    typedef typename Push<T06, L07>::type L06;
    typedef typename Push<T05, L06>::type L05;
    typedef typename Push<T04, L05>::type L04;
    typedef typename Push<T03, L04>::type L03;
    typedef typename Push<T02, L03>::type L02;
    typedef typename Push<T01, L02>::type L01;
    typedef L01 type;
};

/****************************************************************************/
/*                                                                          */
/*                        the actual accumulators                           */
/*                                                                          */
/****************************************************************************/

struct Count
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Count Tag;
        typedef BASE BaseType;

        typedef T argument_type;
        typedef T first_argument_type;
        typedef double result_type;
        
        result_type count_;
        
        Impl()
        : count_(0.0)
        {}
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            ++count_;
        }
        
        void operator()(T const & t, double weight)
        {
            BaseType::operator()(t, weight);
            count_ += weight;
        }
        
        result_type operator()() const
        {
            return count_;
        }
    };
};

struct Minimum
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Minimum Tag;
        typedef BASE BaseType;

        typedef T argument_type;
        typedef T first_argument_type;
        typedef T result_type;

        result_type min_;
        
        Impl()
        : min_(NumericTraits<T>::max())
        {}
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            min_ = std::min(min_, t);
        }
        
        void operator()(T const & t, double weight)
        {
            vigra_precondition(false, "Minimum accumulator does not support weights.");
        }
        
        result_type operator()() const
        {
            return min_;
        }
    };
};

struct Sum
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Sum Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename NumericTraits<T>::RealPromote result_type;

        result_type sum_;
        
        Impl()
        : sum_(result_type())
        {}
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            sum_ += t;
        }
        
        void operator()(T const & t, double weight)
        {
            BaseType::operator()(t, weight);
            sum_ += weight*t;
        }
        
        result_type operator()() const
        {
            return sum_;
        }
    };
};

struct Mean
{
    typedef Select<Sum, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Mean Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename NumericTraits<T>::RealPromote result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
            return get<Sum>(*this) / get<Count>(*this);
        }
    };
};

struct SumSquaredDifferences
{
    typedef Select<Mean, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef SumSquaredDifferences Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename NumericTraits<T>::RealPromote result_type;

        result_type sumOfSquaredDifferences_;
        
        Impl()
        : sumOfSquaredDifferences_(result_type())
        {}
        
        void operator()(T const & t)
        {
            if(get<Count>(*this) != 0.0)
            {
                result_type old_mean = get<Mean>(*this);
                BaseType::operator()(t);
                result_type t1 = t - old_mean;
                result_type t2 = t1 / get<Count>(*this);
                sumOfSquaredDifferences_ += (get<Count>(*this)-1.0)*t1*t2;
            }
            else
            {
                BaseType::operator()(t);
            }
        }
        
        void operator()(T const & t, double weight)
        {
            if(get<Count>(*this) != 0.0)
            {
                result_type old_mean = get<Mean>(*this);
                BaseType::operator()(t, weight);
                result_type t1 = t - old_mean;
                sumOfSquaredDifferences_ +=
                    (t1 * t1 * weight / get<Count>(*this)) * (get<Count>(*this) - weight );
            }
            else
            {
                BaseType::operator()(t, weight);
            }
        }
        
        result_type operator()() const
        {
            return sumOfSquaredDifferences_;
        }
    };
};

struct Variance
{
    typedef Select<SumSquaredDifferences, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Variance Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename NumericTraits<T>::RealPromote result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
            return get<SumSquaredDifferences>(*this) / get<Count>(*this);
        }
    };
};

struct StdDev
{
    typedef Select<Variance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef StdDev Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename NumericTraits<T>::RealPromote result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
            return sqrt(get<Variance>(*this));
        }
    };
};

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
