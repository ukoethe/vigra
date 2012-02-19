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

struct AccumulatorBase 
{
	typedef AccumulatorBase Tag;
    typedef void result_type;
    typedef double second_argument_type;

	void operator+=(AccumulatorBase const &)
    {}
    
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
    
    void reset()
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
    
    void operator+=(DynamicAccumulatorWrapper const & o)
    {
        if(is_active_)
        {
            Base::operator+=(o);
        }
        else
        {
            BaseBase::operator+=(o);
        }
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
    
    template <class U>
    struct Result
    {
        typedef U type;
        
        static void init(type & t)
        {
            t = NumericTraits<U>::max();
        }
        
        template <class V>
        static void updateMin(type & t, V const & v)
        {
            t = min(t, v);
        }
    };
    
    template <unsigned int N, class U, class Stride>
    struct Result<MultiArrayView<N, U, Stride> >
    {
        typedef MultiArray<N, U> type;
        
        static void init(type & t)
        {
            type().swap(t);
        }
        
        template <class V>
        static void updateMin(type & t, V const & v)
        {
            using namespace multi_math;
            
            if(t.size() == 0)
            {
                t.reshape(v.shape(), NumericTraits<U>::max());
            }
            t = min(t, v);
        }
    };
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Minimum Tag;
        typedef BASE BaseType;
        typedef Result<T> ResultTraits;

        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename ResultTraits::type result_type;

        result_type min_;
        
        Impl()
        {
            ResultTraits::init(min_);
        }
        
        void reset()
        {
            ResultTraits::init(min_);
            BaseType::reset();
        }
    
        void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            ResultTraits::updateMin(min_, o.min_);
        }
    
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            ResultTraits::updateMin(min_, t);
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

struct Maximum
{
    typedef Select<> Dependencies;
    
    template <class U>
    struct Result
    {
        typedef U type;
        
        static void init(type & t)
        {
            t = NumericTraits<U>::min();
        }
        
        template <class V>
        static void updateMax(type & t, V const & v)
        {
            t = max(t, v);
        }
    };
    
    template <unsigned int N, class U, class Stride>
    struct Result<MultiArrayView<N, U, Stride> >
    {
        typedef MultiArray<N, U> type;
        
        static void init(type & t)
        {
            type().swap(t);
        }
        
        template <class V>
        static void updateMax(type & t, V const & v)
        {
            using namespace multi_math;
            
            if(t.size() == 0)
            {
                t.reshape(v.shape(), NumericTraits<U>::min());
            }
            t = max(t, v);
        }
    };
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Maximum Tag;
        typedef BASE BaseType;
        typedef Result<T> ResultTraits;

        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename ResultTraits::type result_type;

        result_type max_;
        
        Impl()
        {
            ResultTraits::init(max_);
        }
        
        void reset()
        {
            ResultTraits::init(max_);
            BaseType::reset();
        }
    
        void operator+=(Impl const & o)
        {
            ResultTraits::updateMax(max_, o.max_);
        }
    
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            ResultTraits::updateMax(max_, t);
        }
        
        void operator()(T const & t, double weight)
        {
            vigra_precondition(false, "Maximum accumulator does not support weights.");
        }
        
        result_type operator()() const
        {
            return max_;
        }
    };
};

struct Sum
{
    typedef Select<> Dependencies;
    
    template <class U>
    struct Result
    {
        typedef typename NumericTraits<U>::RealPromote type;
        
        static void init(type & t)
        {
            t = type();
        }
        
        template <class V>
        static void resize(type &, V const &)
        {}
    };
    
    template <unsigned int N, class U, class Stride>
    struct Result<MultiArrayView<N, U, Stride> >
    {
        typedef typename NumericTraits<U>::RealPromote ValueType;
        typedef MultiArray<N, ValueType> type;
        
        static void init(type & t)
        {
            type().swap(t);
        }
        
        static void resize(type & t, MultiArrayView<N, U, Stride> const & u)
        {
            if(t.size() == 0)
                t.reshape(u.shape());
        }
    };
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef Sum Tag;
        typedef BASE BaseType;
        typedef Result<T> ResultTraits;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename ResultTraits::type result_type;

        result_type sum_;
        
        Impl()
        {
            ResultTraits::init(sum_);
        }
        
        void reset()
        {
            ResultTraits::init(sum_);
            BaseType::reset();
        }
    
        void operator+=(Impl const & o)
        {
            BaseType::operator+=(o);
            sum_ += o.sum_;
        }
    
        void operator()(T const & t)
        {
            BaseType::operator()(t);
            ResultTraits::resize(sum_, t);
            sum_ += t;
        }
        
        void operator()(T const & t, double weight)
        {
            BaseType::operator()(t, weight);
            ResultTraits::resize(sum_, t);
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
        typedef typename LookupTag<Sum, Impl>::result_type result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
			using namespace multi_math;
            return get<Sum>(*this) / get<Count>(*this);
        }
    };
};

struct SumSquaredDifferences
{
    typedef Select<Mean, Count> Dependencies;
    
    template <class U>
    struct Result
    {
        typedef typename NumericTraits<U>::RealPromote type;
        
        static void init(type & t)
        {
            t = type();
        }
        
        template <class V>
        static void resize(type &, V)
        {}
        
        template <class V>
        static type prod(V v)
        {
            return v*v;
        }
    };
    
    template <class U, int N>
    struct Result<TinyVector<U, N> >
    {
        typedef typename NumericTraits<U>::RealPromote ValueType;
        typedef Matrix<ValueType> type;
        
        static void init(type & t)
        {
            type(N,N).swap(t);
        }
        
        template <class V>
        static void resize(type &, V const &)
        {}
        
        template <class V>
        static linalg::TemporaryMatrix<ValueType> prod(TinyVector<V, N> const & v)
        {
            linalg::TemporaryMatrix<ValueType> res(N, N);
            for(MultiArrayIndex j=0; j<N; ++j)
                for(MultiArrayIndex i=0; i<N; ++i)
                    res(i,j) = v[i]*v[j];
            return res;
        }
    };
    
    template <unsigned int N, class U, class Stride>
    struct Result<MultiArrayView<N, U, Stride> >
    {
        typedef typename NumericTraits<U>::RealPromote ValueType;
        typedef Matrix<ValueType> type;
        
        static void init(type & t)
        {
            type().swap(t);
        }
        
        static void resize(type & t, MultiArrayView<N, U, Stride> const & u)
        {
            t.reshape(u.size(), u.size());
        }
        
        template <class V, class S>
        static linalg::TemporaryMatrix<ValueType> prod(MultiArrayView<N, V, S> const & v)
        {
            int size = v.size();
            linalg::TemporaryMatrix<ValueType> res(size, size);
            for(MultiArrayIndex j=0; j<size; ++j)
                for(MultiArrayIndex i=0; i<size; ++i)
                    res(i,j) = v[i]*v[j];
            return res;
        }
    };
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef SumSquaredDifferences Tag;
        typedef BASE BaseType;
        typedef Result<T> ResultTraits;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename ResultTraits::type result_type;

        result_type sumOfSquaredDifferences_;
        
        Impl()
        {
            ResultTraits::init(sumOfSquaredDifferences_);
        }
        
        void reset()
        {
            ResultTraits::init(sumOfSquaredDifferences_);
            BaseType::reset();
        }
    
        void operator+=(Impl const & o)
        {
            using namespace vigra::multi_math;            
            typedef typename LookupTag<Mean, Impl>::result_type MeanResult;
            MeanResult diff = get<Mean>(*this) - get<Mean>(o);
            double weight = get<Count>(*this) * get<Count>(o) / (get<Count>(*this) + get<Count>(o));
            sumOfSquaredDifferences_ += o.sumOfSquaredDifferences_ + weight * ResultTraits::prod(diff);
            
            // this must be last because the above computation needs the old values
            BaseType::operator+=(o);
        }
    
        void operator()(T const & t)
        {
			using namespace multi_math;
            if(get<Count>(*this) != 0.0)
            {
                typedef typename LookupTag<Mean, Impl>::result_type MeanResult;
                MeanResult old_mean = get<Mean>(*this);

                BaseType::operator()(t);

				old_mean -= t;
                double adjustedWeight = (get<Count>(*this) - 1.0) / get<Count>(*this);
                sumOfSquaredDifferences_ += adjustedWeight * ResultTraits::prod(old_mean);
            }
            else
            {
                ResultTraits::resize(sumOfSquaredDifferences_, t);
                BaseType::operator()(t);
            }
        }
        
        void operator()(T const & t, double weight)
        {
			using namespace multi_math;
            if(get<Count>(*this) != 0.0)
            {
                typedef typename LookupTag<Mean, Impl>::result_type MeanResult;
                MeanResult old_mean = get<Mean>(*this);

                BaseType::operator()(t, weight);

				old_mean -= t;
                double adjustedWeight = (get<Count>(*this) - weight) / get<Count>(*this) * weight;
                sumOfSquaredDifferences_ += adjustedWeight * ResultTraits::prod(old_mean);
            }
            else
            {
                ResultTraits::resize(sumOfSquaredDifferences_, t);
                BaseType::operator()(t, weight);
            }
        }
        
        result_type operator()() const
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
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename LookupTag<SumSquaredDifferences, Impl>::result_type result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
			using namespace multi_math;
            return get<SumSquaredDifferences>(*this) / get<Count>(*this);
        }
    };
};

typedef Variance Covariance;

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
        typedef typename LookupTag<Variance, Impl>::result_type result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<Variance>(*this));
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
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename LookupTag<SumSquaredDifferences, Impl>::result_type result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
			using namespace multi_math;
            return get<SumSquaredDifferences>(*this) / (get<Count>(*this) - 1.0);
        }
    };
};

typedef UnbiasedVariance UnbiasedCovariance;

struct UnbiasedStdDev
{
    typedef Select<UnbiasedVariance> Dependencies;
    
    template <class T, class BASE>
    struct Impl
    : public BASE
    {
        typedef UnbiasedStdDev Tag;
        typedef BASE BaseType;
        
        typedef T argument_type;
        typedef T first_argument_type;
        typedef typename LookupTag<Variance, Impl>::result_type result_type;

        using BaseType::operator();
        
        result_type operator()() const
        {
			using namespace multi_math;
            return sqrt(get<UnbiasedVariance>(*this));
        }
    };
};

}} // namespace vigra::acc1

#endif // VIGRA_ACCUMULATOR_HXX
