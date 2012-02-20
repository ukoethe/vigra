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
#include "utilities.hxx"
#include "matrix.hxx"
#include "multi_math.hxx"
#include <algorithm>
#include <iostream>

namespace vigra {

namespace acc1 {

// make Select a synonym for MakeTypeList
template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct Select
: public MakeTypeList<T01, T02, T03, T04, T05, T06, T07, T08, T09, T10,
                      T11, T12, T13, T14, T15, T16, T17, T18, T19, T20>
{};

template <class List>
struct PushDependencies;

template <class HEAD, class TAIL>
struct PushDependencies<TypeList<HEAD, TAIL> >
{
    typedef typename PushDependencies<TAIL>::type TailWithDependencies;
    typedef typename HEAD::Dependencies Dependencies; // must be of type Select<...>
    typedef typename PushDependencies<typename Dependencies::type>::type HeadDependencies;
    typedef TypeList<HEAD, HeadDependencies> HeadWithDependencies;
    typedef typename MergeUnique<HeadWithDependencies, TailWithDependencies>::type type;
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

template <class T>
struct AccumulatorTraits
{
    typedef T type;
    typedef T element_type;
    typedef T argument_type;
    typedef T                                      MinmaxType;
    typedef typename NumericTraits<T>::RealPromote SumType;
    typedef typename NumericTraits<T>::RealPromote CovarianceType;
    
    typedef typename Select<>::type NeedReshape;
};

template <class T, int N>
struct AccumulatorTraits<TinyVector<T, N> >
{
    typedef TinyVector<T, N>         type;
    typedef T                        element_type;
    typedef TinyVector<T, N> const & argument_type;
    typedef TinyVector<T, N>                                      MinmaxType;
    typedef TinyVector<typename NumericTraits<T>::RealPromote, N> SumType;
    typedef Matrix<typename NumericTraits<T>::RealPromote>        CovarianceType;
    
    typedef typename Select<CovarianceType>::type NeedReshape;
};

template <unsigned int N, class T, class Stride>
struct AccumulatorTraits<MultiArrayView<N, T, Stride> >
{
    typedef MultiArrayView<N, T, Stride>         type;
    typedef T                                    element_type;
    typedef MultiArrayView<N, T, Stride> const & argument_type;
    typedef MultiArray<N, T>                                      MinmaxType;
    typedef MultiArray<N, typename NumericTraits<T>::RealPromote> SumType;
    typedef Matrix<typename NumericTraits<T>::RealPromote>        CovarianceType;
    
    typedef typename Select<MinmaxType, SumType, CovarianceType>::type NeedReshape;
};

template <class NeedReshape>
struct ReshapeTraits
{
    template <class Array, class Shape, class Initial>
    static void reshape(Array & a, Shape const & s, Initial initial)
    {
        Array(s, initial).swap(a);
    }
    
    template <class Array, class Shape, class Initial>
    static void flatReshape(Array & a, Shape const & s, Initial initial)
    {
        MultiArrayIndex size = prod(s);
        Array(size, size, initial).swap(a);
    }
};

template <>
struct ReshapeTraits<VigraFalseType>
{
    template <class Array, class Shape, class Initial>
    static void reshape(Array &, Shape const &, Initial)
    {}
    
    template <class Array, class Shape, class Initial>
    static void flatReshape(Array &, Shape const &, Initial)
    {}
};

struct AccumulatorBase 
{
	typedef AccumulatorBase Tag;
    
    void operator()() const
    {
        vigra_precondition(false,
            std::string("get(accumulator): attempt to access inactive statistic."));
    }
    
    void activate()
    {}
    
    void reset()
    {}
    
    template <class Shape>
    void reshape(Shape const &)
    {}
};

template <class T>
struct TypedAccumulatorBase
: public AccumulatorBase
{
	typedef AccumulatorBase Tag;
    typedef AccumulatorBase BaseType;

    typedef AccumulatorTraits<T>                       accumulator_traits;
    typedef typename accumulator_traits::argument_type argument_type;
    typedef argument_type                              first_argument_type;
    typedef double                                     second_argument_type;
    typedef void                                       result_type;
        // 'qualified_result_type' can be 'result_type' or 'result_type const &'
        // depending on the return value of operator()()
    typedef void                                       qualified_result_type;

	void operator+=(TypedAccumulatorBase const &)
    {}
    
	void operator()(argument_type)
    {}
    
    void operator()(first_argument_type, second_argument_type)
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
    
    typedef typename Base::argument_type         argument_type;
    typedef typename Base::first_argument_type   first_argument_type;
    typedef typename Base::second_argument_type  second_argument_type;
    typedef typename Base::result_type           result_type;
    typedef typename Base::qualified_result_type qualified_result_type;

    bool is_active_;
    
    DynamicAccumulatorWrapper()
    : is_active_(false)
    {}
    
    void activate()
    {
        is_active_ = true;
        ActivateDependencies<typename Tag::Dependencies::type>::exec(*this);
    }
     
    void reset()
    {
        if(is_active_)
            Base::reset();
        else
            BaseBase::reset();
    }
   
    template <class Shape>
    void reshape(Shape const & s)
    {
        if(is_active_)
            Base::reshape(s);
        else
            BaseBase::reshape(s);
    }
   
    void operator+=(DynamicAccumulatorWrapper const & o)
    {
        if(is_active_)
            Base::operator+=(o);
        else
            BaseBase::operator+=(o);
    }
    
	qualified_result_type operator()() const
    {
        vigra_precondition(is_active_,
            std::string("get(accumulator): attempt to access inactive statistic '") << typeid(Tag).name() << "'.");
        return Base::operator()();
    }
    
	void operator()(T const & t)
    {
        if(is_active_)
            Base::operator()(t);
        else
            BaseBase::operator()(t);
    }
    
    void operator()(T const & t, double weight)
    {
        if(is_active_)
            Base::operator()(t, weight);
        else
            BaseBase::operator()(t, weight);
    }
};

template <class T, class A>
struct NeedsReshape
{
    typedef typename NeedsReshape<T, typename A::BaseType>::type BaseNeedsReshape;
    typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, 
                              typename A::result_type>::type CurrentNeedsReshape;
    typedef typename Or<BaseNeedsReshape, CurrentNeedsReshape>::type type;
};

template <class T>
struct NeedsReshape<T, AccumulatorBase>
{
    typedef VigraFalseType type;
};

template <class T, class Base, class NeedsReshape=typename NeedsReshape<T, Base>::type>
struct ReshapeAccumulator
: public Base
{
    typedef VigraFalseType Tag;
    typedef Base BaseType;
    
    bool needs_reshape_;
    
    ReshapeAccumulator()
    : needs_reshape_(true)
    {}
    
    using Base::operator();
    
	void operator()(T const & t)
    {
        if(needs_reshape_)
        {
            Base::reshape(shape(t));
            needs_reshape_ = false;
        }
        Base::operator()(t);
    }
    
    void operator()(T const & t, double weight)
    {
        if(needs_reshape_)
        {
            Base::reshape(shape(t));
            needs_reshape_ = false;
        }
        Base::operator()(t, weight);
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

template <class T, class Base>
struct ReshapeAccumulator<T, Base, VigraFalseType>
: public Base
{
    typedef VigraFalseType Tag;
    typedef Base BaseType;
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
    typedef TypedAccumulatorBase<T> BaseType;
    typedef typename Tag::template Impl<T, BaseType> type; 
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
    typedef TypedAccumulatorBase<T> BaseType;
    typedef DynamicAccumulatorWrapper<T, Tag, BaseType> type; 
};

template <class T, class Selected>
struct Accumulator
: public ReshapeAccumulator<T, 
                 typename Compose<T, typename PushDependencies<typename Selected::type>::Head,
                                     typename PushDependencies<typename Selected::type>::Tail>::type>
{
    typedef typename PushDependencies<typename Selected::type>::type Accumulators;
    typedef typename ReshapeAccumulator<T, typename Compose<T, typename Accumulators::Head, 
                                                      typename Accumulators::Tail>::type> BaseType;
    typedef VigraFalseType Tag;
};

template <class T, class Selected>
struct DynamicAccumulator
: public ReshapeAccumulator<T, typename DynamicCompose<T, typename PushDependencies<typename Selected::type>::Head,
                                                 typename PushDependencies<typename Selected::type>::Tail>::type>
{
    typedef typename PushDependencies<typename Selected::type>::type Accumulators;
    typedef typename ReshapeAccumulator<T, typename DynamicCompose<T, typename Accumulators::Head, 
                                                             typename Accumulators::Tail>::type> BaseType;
    typedef VigraFalseType Tag;
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

        typedef typename BaseType::argument_type        argument_type;
        typedef typename BaseType::first_argument_type  first_argument_type;
        typedef typename BaseType::second_argument_type second_argument_type;
        typedef double                                  result_type;
        typedef double                                  qualified_result_type;
        
        result_type count_;
        
        Impl()
        : count_(0.0)
        {}
        
        void reset()
        {
            count_ = 0.0;
            BaseType::reset();
        }
        
        void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            count_ += o.count_;
        }
    
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

        typedef typename AccumulatorTraits<T>::element_type element_type;
        typedef typename BaseType::argument_type            argument_type;
        typedef typename BaseType::first_argument_type      first_argument_type;
        typedef typename BaseType::second_argument_type     second_argument_type;
        typedef typename AccumulatorTraits<T>::MinmaxType   result_type;
        typedef result_type const &                         qualified_result_type;

        result_type min_;
        
        Impl()
        {
            min_ = NumericTraits<element_type>::max();
        }
        
        void reset()
        {
            min_ = NumericTraits<element_type>::max();
            BaseType::reset();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, result_type>::type NeedReshape;
            ReshapeTraits<NeedReshape>::reshape(min_, s, NumericTraits<element_type>::max());
            BaseType::reshape(s);
        }
        
       void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            using namespace multi_math;
            min_ = min(min_, o.min_);
        }
    
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            using namespace multi_math;
            min_ = min(min_, t);
        }
        
        void operator()(T const & t, double weight)
        {
            vigra_precondition(false, "Minimum accumulator does not support weights.");
        }
        
        qualified_result_type operator()() const
        {
            return min_;
        }
    };
};

struct Maximum
{
    typedef Select<> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Maximum Tag;
        typedef BASE BaseType;

        typedef typename AccumulatorTraits<T>::element_type element_type;
        typedef typename BaseType::argument_type            argument_type;
        typedef typename BaseType::first_argument_type      first_argument_type;
        typedef typename BaseType::second_argument_type     second_argument_type;
        typedef typename AccumulatorTraits<T>::MinmaxType   result_type;
        typedef result_type const &                         qualified_result_type;

        result_type max_;
        
        Impl()
        {
            max_ = NumericTraits<element_type>::min();
        }
        
        void reset()
        {
            max_ = NumericTraits<element_type>::min();
            BaseType::reset();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, result_type>::type NeedReshape;
            ReshapeTraits<NeedReshape>::reshape(max_, s, NumericTraits<element_type>::min());
            BaseType::reshape(s);
        }
        
        void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            using namespace multi_math;
            max_ = max(max_, o.max_);
        }
    
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            using namespace multi_math;
            max_ = max(max_, t);
        }
        
        void operator()(T const & t, double weight)
        {
            vigra_precondition(false, "Maximum accumulator does not support weights.");
        }
        
        qualified_result_type operator()() const
        {
            return max_;
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
        
        typedef typename AccumulatorTraits<T>::element_type element_type;
        typedef typename BaseType::argument_type            argument_type;
        typedef typename BaseType::first_argument_type      first_argument_type;
        typedef typename BaseType::second_argument_type     second_argument_type;
        typedef typename AccumulatorTraits<T>::SumType      result_type;
        typedef result_type const &                         qualified_result_type;

        result_type sum_;
        
        Impl()
        {
            sum_ = element_type();
        }
        
        void reset()
        {
            sum_ = element_type();
            BaseType::reset();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, result_type>::type NeedReshape;
            ReshapeTraits<NeedReshape>::reshape(sum_, s, element_type());
            BaseType::reshape(s);
        }
        
        void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            sum_ += o.sum_;
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
        
        qualified_result_type operator()() const
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
        
        typedef typename BaseType::argument_type           argument_type;
        typedef typename BaseType::first_argument_type     first_argument_type;
        typedef typename BaseType::second_argument_type    second_argument_type;
        typedef typename LookupTag<Sum, Impl>::result_type result_type;
        typedef result_type                                qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
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
        
        typedef typename AccumulatorTraits<T>::element_type   element_type;
        typedef typename BaseType::argument_type              argument_type;
        typedef typename BaseType::first_argument_type        first_argument_type;
        typedef typename BaseType::second_argument_type       second_argument_type;
        typedef typename AccumulatorTraits<T>::SumType        result_type;
        typedef result_type const &                           qualified_result_type;
       
        result_type sumOfSquaredDifferences_;
        
        Impl()
        {
            sumOfSquaredDifferences_ = element_type();
        }
        
        void reset()
        {
            sumOfSquaredDifferences_ = element_type();
            BaseType::reset();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, result_type>::type NeedReshape;
            ReshapeTraits<NeedReshape>::reshape(sumOfSquaredDifferences_, s, element_type());
            BaseType::reshape(s);
        }
        
        void operator+=(Impl const & o)
        {
            using namespace vigra::multi_math;            

            double weight = get<Count>(*this) * get<Count>(o) / (get<Count>(*this) + get<Count>(o));
            sumOfSquaredDifferences_ += o.sumOfSquaredDifferences_ + 
                                        weight * sq(get<Sum>(*this) / get<Count>(*this) - get<Sum>(o) / get<Count>(o));
            
            // this must be last because the above computation needs the old values
            BaseType::operator+=(o);
        }
    
        void operator()(T const & t)
        {
            double old_count = get<Count>(*this);
            if(old_count != 0.0)
            {
                using namespace vigra::multi_math;
                
                sumOfSquaredDifferences_ += old_count / (old_count + 1.0) * sq(get<Sum>(*this) / old_count - t);
            }
            BaseType::operator()(t);
        }
        
        void operator()(T const & t, double weight)
        {
            double old_count = get<Count>(*this);
            if(old_count != 0.0)
            {
                using namespace vigra::multi_math;
                
                sumOfSquaredDifferences_ += old_count / (old_count + weight) * weight * sq(get<Sum>(*this) / old_count - t);
            }
            BaseType::operator()(t, weight);
        }
        
        qualified_result_type operator()() const
        {
            return sumOfSquaredDifferences_;
        }
    };
};

typedef SumSquaredDifferences SSD;

struct Variance
{
    typedef Select<SumSquaredDifferences, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Variance Tag;
        typedef BASE BaseType;
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<SumSquaredDifferences, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;


        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
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
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<Variance, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<SumSquaredDifferences>(*this) / get<Count>(*this));
        }
    };
};

struct UnbiasedVariance
{
    typedef Select<SumSquaredDifferences, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef UnbiasedVariance Tag;
        typedef BASE BaseType;
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<SumSquaredDifferences, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
            return get<SumSquaredDifferences>(*this) / (get<Count>(*this) - 1.0);
        }
    };
};

struct UnbiasedStdDev
{
    typedef Select<UnbiasedVariance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef UnbiasedStdDev Tag;
        typedef BASE BaseType;
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<Variance, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<SumSquaredDifferences>(*this) / (get<Count>(*this) - 1.0));
        }
    };
};

struct ScatterMatrix
{
    typedef Select<Mean, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef ScatterMatrix Tag;
        typedef BASE BaseType;
        
        typedef typename AccumulatorTraits<T>::element_type   element_type;
        typedef typename BaseType::argument_type              argument_type;
        typedef typename BaseType::first_argument_type        first_argument_type;
        typedef typename BaseType::second_argument_type       second_argument_type;
        typedef typename AccumulatorTraits<T>::CovarianceType result_type;
        typedef result_type const &                           qualified_result_type;
       
        typedef typename AccumulatorTraits<T>::SumType        SumType;

        result_type scatter_matrix_;
        SumType     diff_;
        
        Impl()
        {
            scatter_matrix_ = element_type();
        }
        
        void reset()
        {
            scatter_matrix_ = element_type();
            BaseType::reset();
        }
    
        template <class Shape>
        void reshape(Shape const & s)
        {
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, result_type>::type NeedReshape;
            // we always compute the scatter matrix from flattened arrays
            ReshapeTraits<NeedReshape>::flatReshape(scatter_matrix_, s, element_type());
            typedef typename Contains<typename AccumulatorTraits<T>::NeedReshape, SumType>::type SumNeedReshape;
            ReshapeTraits<SumNeedReshape>::reshape(diff_, s, element_type());
            BaseType::reshape(s);
        }
        
        void operator+=(Impl const & o)
        {
            using namespace vigra::multi_math;            
            diff_ = get<Sum>(*this) / get<Count>(*this) - get<Sum>(o) / get<Count>(o);
            double weight = get<Count>(*this) * get<Count>(o) / (get<Count>(*this) + get<Count>(o));
            addWeightedOuterProduct(scatter_matrix_, diff_, weight);
            scatter_matrix_ += o.scatter_matrix_;
            
            // this must be last because the above computation needs the old values
            BaseType::operator+=(o);
        }
    
        void operator()(T const & t)
        {
            double old_count = get<Count>(*this);
            if(old_count != 0.0)
            {
                using namespace vigra::multi_math;
                diff_ = get<Sum>(*this) / old_count - t;
                double weight = old_count / (old_count + 1.0);
                addWeightedOuterProduct(scatter_matrix_, diff_, weight);
            }
            BaseType::operator()(t);
        }
        
        void operator()(T const & t, double weight)
        {
            double old_count = get<Count>(*this);
            if(old_count != 0.0)
            {
                using namespace vigra::multi_math;
                diff_ = get<Sum>(*this) / old_count - t;
                double weight = old_count / (old_count + weight) * weight;
                addWeightedOuterProduct(scatter_matrix_, diff_, weight);
            }
            BaseType::operator()(t, weight);
        }
        
        qualified_result_type operator()() const
        {
            return scatter_matrix_;
        }
        
      private:
        template <class SSD>
        static void addWeightedOuterProduct(SSD & ssd, SumType const & m, double w)
        {
            int size = m.size();
            for(MultiArrayIndex j=0; j<size; ++j)
                for(MultiArrayIndex i=0; i<size; ++i)
                    ssd(i,j) += w*m[i]*m[j];
        }
        
        static void addWeightedOuterProduct(double & ssd, SumType const & m, double w)
        {
            ssd += w*m*m;
        }
    };
};

struct Covariance
{
    typedef Select<ScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Covariance Tag;
        typedef BASE BaseType;
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<ScatterMatrix, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
            return get<ScatterMatrix>(*this) / get<Count>(*this);
        }
    };
};

struct UnbiasedCovariance
{
    typedef Select<ScatterMatrix, Count> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef UnbiasedCovariance Tag;
        typedef BASE BaseType;
        
        typedef typename BaseType::argument_type         argument_type;
        typedef typename BaseType::first_argument_type   first_argument_type;
        typedef typename BaseType::second_argument_type  second_argument_type;
        typedef typename LookupTag<ScatterMatrix, Impl>::result_type result_type;
        typedef result_type                              qualified_result_type;

        using BaseType::operator();
        
        qualified_result_type operator()() const
        {
			using namespace multi_math;
            return get<ScatterMatrix>(*this) / (get<Count>(*this) - 1.0);
        }
    };
};

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
