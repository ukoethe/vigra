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
struct Cons
{
    typedef AccumulatorList<A, typename Tail::type> type;
};

template <class A>
struct Cons<A, void>
{
    typedef AccumulatorList<A, void> type;
};

template <>
struct Cons<void, void>
{
    typedef void type;
};

template <class List1, class List2>
struct Merge;

template <class Head, class Tail, class List>
struct Merge<AccumulatorList<Head, Tail>, List>
{
    typedef typename Merge<Tail, List>::type Rest;
    typedef AccumulatorList<Head, Rest> type;
};

template <class List>
struct Merge<void, List>
{
    typedef List type;
};

template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct Select
{
    typedef typename Cons<T01, Cons<T02, Cons<T03, Cons<T04, Cons<T05,
                     Cons<T06, Cons<T07, Cons<T08, Cons<T09, Cons<T10,
                     Cons<T11, Cons<T12, Cons<T13, Cons<T14, Cons<T15,
                     Cons<T16, Cons<T17, Cons<T18, Cons<T19, Cons<T20
                     > > > > > > > > > > > > > > > > > > > >::type type;
};

struct Nothing;

template <class Tag, class AccuList>
struct FindTag
{
    typedef typename IsSameType<Tag, typename AccuList::Tag>::type Found;
    typedef typename If<Found, AccuList, typename FindTag<Tag, typename AccuList::BaseType>::type>::type type;
};

template <class Tag>
struct FindTag<Tag, Nothing>
{
    typedef VigraFalseType type;
};

template <class TAG, class AccuList>
double get(AccuList const & a);

template <class TAG, class AccuList>
void activate(AccuList & a);

struct Nothing 
{
	typedef Nothing Tag;

    template <class T>
	void operator()(T const &)
    {}
    
    template <class T>
	void dynamic(T const &)
    {}
    
    template <class T>
    void operator()(T const &, double)
    {}
    
    template <class T, class AccuList, class Selection>
    struct Compose
    {
        typedef Nothing type;
    };
};

#define VIGRA_COMPOSE_ACCUMULATOR(name) \
template <class T, class AccuList, class Selection> \
struct Compose \
{ \
    typedef typename Contains<Selection, name>::type SelectThis; \
    typedef typename AccuList::Head NextType; \
    typedef typename If<SelectThis, \
	          typename Merge<Dependencies::type, Selection>::type, Selection>::type NewSelection; \
    typedef typename NextType::template Compose<T, typename AccuList::Tail, NewSelection>::type BaseType; \
    typedef typename If<SelectThis, Impl<T, BaseType>, BaseType>::type type; \
}; \
 \
template <class T, class Selection> \
struct Compose<T, void, Selection> \
{ \
    typedef typename Contains<Selection, name>::type SelectThis; \
    typedef typename If<SelectThis, Impl<T, Nothing>, Nothing>::type type; \
};

struct Count
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef double result_type;
        typedef Count Tag;
        typedef BASE BaseType;

        result_type count_;
        bool is_active_;
        
        Impl()
        : count_(0.0),
          is_active_(false)
        {}
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            ++count_;
        }
        
        void dynamic(T const & t)
        {
            BaseType::dynamic(t);
            if(is_active_)
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
        
        void activate()
        {
            is_active_ = true;
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(Count)
};

struct Minimum
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef T result_type;
        typedef Minimum Tag;
        typedef BASE BaseType;

        result_type min_;
        bool is_active_;
        
        Impl()
        : min_(NumericTraits<T>::max()),
          is_active_(false)
        {}
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            min_ = std::min(min_, t);
        }
        
        void dynamic(T const & t)
        {
            BaseType::dynamic(t);
            if(is_active_)
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
        
        void activate()
        {
            is_active_ = true;
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(Minimum)
};

struct Sum
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename NumericTraits<T>::RealPromote result_type;
        typedef Sum Tag;
        typedef BASE BaseType;
        
        result_type sum_;
        bool is_active_;
        
        Impl()
        : sum_(result_type()),
          is_active_(false)
        {}
        
        void dynamic(T const & t)
        {
            BaseType::dynamic(t);
            if(is_active_)
                sum_ += t;
        }
        
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
        
        void activate()
        {
            is_active_ = true;
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(Sum)
};

struct Mean
{
    typedef Select<Sum, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename NumericTraits<T>::RealPromote result_type;
        typedef Mean Tag;
        typedef BASE BaseType;
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
        }
        
        void dynamic(T const & t)
        {
            BaseType::dynamic(t);
        }
        
        void operator()(T const & t, double weight)
        {
            BaseType::operator()(t, weight);
        }
        
        result_type operator()() const
        {
            return get<Count>(*this) > 0
                       ? get<Sum>(*this) / get<Count>(*this)
					   : result_type();
        }
        
        void activate()
        {
            activate<Dependencies::type>(*this);
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(Mean)
};

struct Variance
{
    typedef Select<Mean> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename NumericTraits<T>::RealPromote result_type;
        typedef Variance Tag;
        typedef BASE BaseType;
        
        result_type sumOfSquaredDifferences_;
        bool is_active_;
        
        Impl()
        : sumOfSquaredDifferences_(result_type()),
          is_active_(false)
        {}
        
        void operator()(T const & t)
        {
            result_type old_mean = get<Mean>(*this);
            BaseType::operator()(t);
            result_type t1 = t - old_mean;
            result_type t2 = t1 / get<Count>(*this);
            sumOfSquaredDifferences_ += (get<Count>(*this)-1.0)*t1*t2;
        }
        
        void dynamic(T const & t)
        {
            result_type old_mean = get<Mean>(*this);
            BaseType::dynamic(t);
            if(is_active_)
            {
                result_type t1 = t - old_mean;
                result_type t2 = t1 / get<Count>(*this);
                sumOfSquaredDifferences_ += (get<Count>(*this)-1.0)*t1*t2;
            }
        }
        
        void operator()(T const & t, double weight)
        {
            result_type old_mean = get<Mean>(*this);
            BaseType::operator()(t, weight);
            result_type t1 = t - old_mean;
            if(get<Count>(*this) > weight)
                sumOfSquaredDifferences_ +=
                    (t1 * t1 * weight / get<Count>(*this)) * (get<Count>(*this) - weight );
        }
        
        result_type operator()() const
        {
            // FIXME: also support unbiased variance
            return sumOfSquaredDifferences_ / get<Count>(*this);
        }
        
        void activate()
        {
            is_active_ = true;
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(Variance)
};

struct StdDev
{
    typedef Select<Variance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef typename NumericTraits<T>::RealPromote result_type;
        typedef StdDev Tag;
        typedef BASE BaseType;
        
        void operator()(T const & t)
        {
            BaseType::operator()(t);
        }
        
        void dynamic(T const & t)
        {
            BaseType::dynamic(t);
        }
        
        void operator()(T const & t, double weight)
        {
            BaseType::operator()(t, weight);
        }
        
        result_type operator()() const
        {
            return sqrt(get<Variance>(*this));
        }
        
        void activate()
        {
            activate<Dependencies::type>(*this);
        }
    };
    
    VIGRA_COMPOSE_ACCUMULATOR(StdDev)
};

typedef Cons<StdDev, 
        Cons<Variance, 
        Cons<Mean, 
        Cons<Sum, 
        Cons<Minimum, 
        Cons<Count> > > > > >::type SortedAccumulators;
        
template <class T, class Selected>
struct Accumulator
: public SortedAccumulators::Head::template 
               Compose<T, typename SortedAccumulators::Tail, typename Selected::type>::type
{
    typedef typename SortedAccumulators::Head::template 
               Compose<T, typename SortedAccumulators::Tail, typename Selected::type>::type BaseType;
    typedef Nothing Tag;
        
    void operator()(T const & t)
    {
        BaseType::operator()(t);
    }
    
    void dynamic(T const & t)
    {
        BaseType::dynamic(t);
    }
    
    void operator()(T const & t, double weight)
    {
        BaseType::operator()(t, weight);
    }
};
        
template <class Tag, class AccuType>
struct GetImpl
{
    static double exec(AccuType const & a)
    {
        return a();
    }
};

template <class Tag>
struct GetImpl<Tag, VigraFalseType>
{
    template <class AccuType>
    static double exec(AccuType const & a)
    {
        vigra_precondition(false,
		 std::string("get(accumulator): attempt to access unselected statistic '") << typeid(Tag).name() << "'.");
		return 0.0;
    }
};

// FIXME: determine correct return type
template <class Tag, class AccuList>
double get(AccuList const & a)
{
    return GetImpl<Tag, typename FindTag<Tag, AccuList>::type>::exec(a);
}

template <class Tag, class AccuList, class AccuType=typename FindTag<Tag, AccuList>::type>
struct ActivateImpl
{
    static void exec(AccuType & a)
    {
        a.activate();
    }
};

template <class Tag, class AccuList>
struct ActivateImpl<Tag, AccuList, VigraFalseType>
{
    template <class AccuType>
    static void exec(AccuType & a)
    {
        vigra_precondition(false,
		 std::string("activate(accumulator): attempt to access unselected statistic '") << typeid(Tag).name() << "'.");
    }
};

template <class Head, class Tail, class AccuList>
struct ActivateImpl<AccumulatorList<Head, Tail>, AccuList>
{
    template <class AccuType>
    static void exec(AccuType & a)
    {
        activate<Head>(a);
        ActivateImpl<Tail, AccuList>::exec(a);
    }
};

template <class AccuList>
struct ActivateImpl<void, AccuList>
{
    template <class AccuType>
    static void exec(AccuType &)
    {}
};

template <class Tag, class AccuList>
void activate(AccuList & a)
{
    ActivateImpl<Tag, AccuList>::exec(a);
}

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
