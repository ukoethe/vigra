/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2008 by Gunnar Kedenburg and Ullrich Koethe       */
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


#ifndef VIGRA_MULTI_ARRAY_HXX
#define VIGRA_MULTI_ARRAY_HXX

#include <memory>
#include <algorithm>
#include "accessor.hxx"
#include "tinyvector.hxx"
#include "rgbvalue.hxx"
#include "basicimage.hxx"
#include "imageiterator.hxx"
#include "numerictraits.hxx"
#include "multi_iterator.hxx"
#include "multi_pointoperators.hxx"
#include "metaprogramming.hxx"
#include "mathutil.hxx"
#include "algorithm.hxx"

// Bounds checking Macro used if VIGRA_CHECK_BOUNDS is defined.
#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(this->isInside(diff), "Index out of bounds")
#else
#define VIGRA_ASSERT_INSIDE(diff)
#endif

namespace vigra
{

namespace detail
{

/********************************************************/
/*                                                      */
/*                     MaybeStrided                     */
/*                                                      */
/********************************************************/

/* metatag implementing a test for marking MultiArrays that were
    indexed at the zero'th dimension as strided, and all others as
    unstrided.

    <b>\#include</b> \<vigra/multi_array.hxx\> <br/>
    Namespace: vigra::detail
*/
template <class StrideTag, unsigned int N>
struct MaybeStrided
{
    typedef StrideTag type;
};

template <class StrideTag>
struct MaybeStrided <StrideTag, 0>
{
    typedef StridedArrayTag type;
};

/********************************************************/
/*                                                      */
/*                MultiIteratorChooser                  */
/*                                                      */
/********************************************************/

/* metatag implementing a test (by pattern matching) for marking
    MultiArrays that were indexed at the zero'th dimension as strided.

    <b>\#include</b> \<vigra/multi_array.hxx\> <br/>
    Namespace: vigra::detail
*/
template <class O>
struct MultiIteratorChooser;

/********************************************************/
/*                                                      */
/*       MultiIteratorChooser <StridedArrayTag>         */
/*                                                      */
/********************************************************/

/* specialization of the MultiIteratorChooser for strided arrays.

    <b>\#include</b> \<vigra/multi_array.hxx\> <br/>
    Namespace: vigra::detail
*/
template <>
struct MultiIteratorChooser <StridedArrayTag>
{
    template <unsigned int N, class T, class REFERENCE, class POINTER>
    struct Traverser
    {
        typedef StridedMultiIterator <N, T, REFERENCE, POINTER> type;
    };
    
    template <unsigned int N, class T, class REFERENCE, class POINTER>
    struct Iterator
    {
        typedef StridedScanOrderIterator <N, T, REFERENCE, POINTER> type;
    };
    
    template <class Iter, class View>
    static Iter constructIterator(View * v)
    {
        return v->begin();
    }
};

/********************************************************/
/*                                                      */
/*      MultiIteratorChooser <UnstridedArrayTag>        */
/*                                                      */
/********************************************************/

/* specialization of the MultiIteratorChooser for unstrided arrays.

    <b>\#include</b> \<vigra/multi_array.hxx\> <br/>
    Namespace: vigra::detail
*/
template <>
struct MultiIteratorChooser <UnstridedArrayTag>
{
    template <unsigned int N, class T, class REFERENCE, class POINTER>
    struct Traverser
    {
        typedef MultiIterator <N, T, REFERENCE, POINTER> type;
    };
    
    template <unsigned int N, class T, class REFERENCE, class POINTER>
    struct Iterator
    {
        typedef POINTER type;
    };
    
    template <class Iter, class View>
    static Iter constructIterator(View * v)
    {
        return v->data();
    }
};

/********************************************************/
/*                                                      */
/*                   helper functions                   */
/*                                                      */
/********************************************************/

template <class DestIterator, class Shape, class T>
inline void
initMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<0>)
{
    DestIterator dend = d + shape[0];
    for(; d < dend; ++d)
    {
        *d = init;
    }
}

template <class DestIterator, class Shape, class T, int N>
void
initMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<N>)
{
    DestIterator dend = d + shape[N];
    for(; d < dend; ++d)
    {
        initMultiArrayData(d.begin(), shape, init, MetaInt<N-1>());
    }
}

// FIXME: the explicit overload for MultiIterator<1, UInt8, ... > works around a compiler crash in VisualStudio 2010
#define VIGRA_COPY_MULTI_ARRAY_DATA(name, op) \
template <class SrcIterator, class Shape, class DestIterator> \
inline void \
name##MultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>) \
{     \
    for(MultiArrayIndex i=0; i < shape[0]; ++i, ++s, ++d) \
    { \
        *d op detail::RequiresExplicitCast<typename DestIterator::value_type>::cast(*s); \
    } \
} \
 \
template <class Ref, class Ptr, class Shape, class DestIterator> \
inline void \
name##MultiArrayData(MultiIterator<1, UInt8, Ref, Ptr> si, Shape const & shape, DestIterator d, MetaInt<0>) \
{ \
    Ptr s = &(*si); \
    for(MultiArrayIndex i=0; i < shape[0]; ++i, ++s, ++d) \
    { \
        *d op detail::RequiresExplicitCast<typename DestIterator::value_type>::cast(*s); \
    } \
} \
\
template <class SrcIterator, class Shape, class DestIterator, int N> \
void \
name##MultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<N>) \
{ \
    for(MultiArrayIndex i=0; i < shape[N]; ++i, ++s, ++d) \
    { \
        name##MultiArrayData(s.begin(), shape, d.begin(), MetaInt<N-1>()); \
    } \
} \
\
template <class DestIterator, class Shape, class T> \
inline void \
name##ScalarMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<0>) \
{     \
    for(MultiArrayIndex i=0; i < shape[0]; ++i, ++d) \
    { \
        *d op detail::RequiresExplicitCast<typename DestIterator::value_type>::cast(init); \
    } \
} \
 \
template <class DestIterator, class Shape, class T, int N> \
void \
name##ScalarMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<N>) \
{     \
    for(MultiArrayIndex i=0; i < shape[N]; ++i, ++d) \
    { \
        name##ScalarMultiArrayData(d.begin(), shape, init, MetaInt<N-1>()); \
    } \
}

VIGRA_COPY_MULTI_ARRAY_DATA(copy, =)
VIGRA_COPY_MULTI_ARRAY_DATA(copyAdd, +=)
VIGRA_COPY_MULTI_ARRAY_DATA(copySub, -=)
VIGRA_COPY_MULTI_ARRAY_DATA(copyMul, *=)
VIGRA_COPY_MULTI_ARRAY_DATA(copyDiv, /=)

#undef VIGRA_COPY_MULTI_ARRAY_DATA

template <class SrcIterator, class Shape, class T, class ALLOC>
inline void
uninitializedCopyMultiArrayData(SrcIterator s, Shape const & shape, T * & d, ALLOC & a, MetaInt<0>)
{
    SrcIterator send = s + shape[0];
    for(; s < send; ++s, ++d)
    {
        a.construct(d, static_cast<T const &>(*s));
    }
}

// FIXME: this overload works around a compiler crash in VisualStudio 2010
template <class Ref, class Ptr, class Shape, class T, class ALLOC>
inline void
uninitializedCopyMultiArrayData(MultiIterator<1, UInt8, Ref, Ptr> si, Shape const & shape, T * & d, ALLOC & a, MetaInt<0>)
{
    Ptr s = &(*si), send = s + shape[0];
    for(; s < send; ++s, ++d)
    {
        a.construct(d, static_cast<T const &>(*s));
    }
}

template <class SrcIterator, class Shape, class T, class ALLOC, int N>
void
uninitializedCopyMultiArrayData(SrcIterator s, Shape const & shape, T * & d, ALLOC & a, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s < send; ++s)
    {
        uninitializedCopyMultiArrayData(s.begin(), shape, d, a, MetaInt<N-1>());
    }
}

template <class SrcIterator, class Shape, class T, class Functor>
inline void
reduceOverMultiArray(SrcIterator s, Shape const & shape, T & result, Functor const & f, MetaInt<0>)
{
    SrcIterator send = s + shape[0];
    for(; s < send; ++s)
    {
        f(result, *s);
    }
}

template <class SrcIterator, class Shape, class T, class Functor, int N>
void
reduceOverMultiArray(SrcIterator s, Shape const & shape, T & result, Functor const & f, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s < send; ++s)
    {
        reduceOverMultiArray(s.begin(), shape, result, f, MetaInt<N-1>());
    }
}

struct MaxNormReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        T v = norm(u);
        if(result < v)
            result = v;
    }
};

struct L1NormReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result += norm(u);
    }
};

struct SquaredL2NormReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result += squaredNorm(u);
    }
};

template <class T>
struct WeightedL2NormReduceFunctor
{
    T scale;

    WeightedL2NormReduceFunctor(T s)
    : scale(s)
    {}

    template <class U>
    void operator()(T & result, U const & u) const
    {
        result += squaredNorm(u * scale);
    }
};

struct SumReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result += u;
    }
};

struct ProdReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result *= u;
    }
};

struct MinmaxReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        if(u < result.first)
            result.first = u;
        if(result.second < u)
            result.second = u;
    }
};

struct MeanVarianceReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        ++result.first;
        typename T::second_type t1 = u - result.second;
        typename T::second_type t2 = t1 / result.first;
        result.second += t2;
        result.third += (result.first-1.0)*t1*t2;
    }
};

struct AllTrueReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result = result && (u != NumericTraits<U>::zero());
    }
};

struct AnyTrueReduceFunctor
{
    template <class T, class U>
    void operator()(T & result, U const & u) const
    {
        result = result || (u != NumericTraits<U>::zero());
    }
};

template <class SrcIterator, class Shape, class DestIterator>
inline bool
equalityOfMultiArrays(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>)
{
    SrcIterator send = s + shape[0];
    for(; s < send; ++s, ++d)
    {
        if(!(*s == *d))
            return false;
    }
    return true;
}

template <class SrcIterator, class Shape, class DestIterator, int N>
bool
equalityOfMultiArrays(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s < send; ++s, ++d)
    {
        if(!equalityOfMultiArrays(s.begin(), shape, d.begin(), MetaInt<N-1>()))
            return false;
    }
    return true;
}


template <class SrcIterator, class Shape, class DestIterator>
inline void
swapDataImpl(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>)
{
    SrcIterator send = s + shape[0];
    for(; s < send; ++s, ++d)
        std::swap(*s, *d);
}

template <class SrcIterator, class Shape, class DestIterator, int N>
void
swapDataImpl(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<N>)
{
    SrcIterator send = s + shape[N];
    for(; s < send; ++s, ++d)
        swapDataImpl(s.begin(), shape, d.begin(), MetaInt<N-1>());
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                     MultiArrayView                   */
/*                                                      */
/********************************************************/

// forward declarations

namespace multi_math {

template <class T>
struct MultiMathOperand;

namespace math_detail {

template <unsigned int N, class T, class C, class E>
void assign(MultiArrayView<N, T, C>, MultiMathOperand<E> const &);

template <unsigned int N, class T, class C, class E>
void plusAssign(MultiArrayView<N, T, C>, MultiMathOperand<E> const &);

template <unsigned int N, class T, class C, class E>
void minusAssign(MultiArrayView<N, T, C>, MultiMathOperand<E> const &);

template <unsigned int N, class T, class C, class E>
void multiplyAssign(MultiArrayView<N, T, C>, MultiMathOperand<E> const &);

template <unsigned int N, class T, class C, class E>
void divideAssign(MultiArrayView<N, T, C>, MultiMathOperand<E> const &);

template <unsigned int N, class T, class A, class E>
void assignOrResize(MultiArray<N, T, A> &, MultiMathOperand<E> const &);

template <unsigned int N, class T, class A, class E>
void plusAssignOrResize(MultiArray<N, T, A> &, MultiMathOperand<E> const &);

template <unsigned int N, class T, class A, class E>
void minusAssignOrResize(MultiArray<N, T, A> &, MultiMathOperand<E> const &);

template <unsigned int N, class T, class A, class E>
void multiplyAssignOrResize(MultiArray<N, T, A> &, MultiMathOperand<E> const &);

template <unsigned int N, class T, class A, class E>
void divideAssignOrResize(MultiArray<N, T, A> &, MultiMathOperand<E> const &);

} // namespace math_detail

} // namespace multi_math

template <class T> class FindSum;

struct UnsuitableTypeForExpandElements {};

template <class T>
struct ExpandElementResult
{
    typedef UnsuitableTypeForExpandElements type;
};

template <class T>
struct ExpandElementResult<std::complex<T> >
{
    typedef T type;
    enum { size = 2 };
};

template <class T>
class FFTWComplex;

template <class T>
struct ExpandElementResult<FFTWComplex<T> >
{
    typedef T type;
    enum { size = 2 };
};

template <class T, int SIZE>
struct ExpandElementResult<TinyVector<T, SIZE> >
{
    typedef T type;
    enum { size = SIZE };
};

template <class T, unsigned int R, unsigned int G, unsigned int B>
struct ExpandElementResult<RGBValue<T, R, G, B> >
{
    typedef T type;
    enum { size = 3 };
};

#define VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(TYPE) \
template <>  \
struct ExpandElementResult<TYPE> \
{ \
    typedef TYPE type; \
    enum { size = 1 }; \
}; \

VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(bool)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(char)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(signed char)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(signed short)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(signed int)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(signed long)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(signed long long)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(unsigned char)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(unsigned short)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(unsigned int)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(unsigned long)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(unsigned long long)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(float)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(double)
VIGRA_DEFINE_EXPAND_ELEMENT_RESULT(long double)

#undef VIGRA_DEFINE_EXPAND_ELEMENT_RESULT


/********************************************************/
/*                                                      */
/*                       NormTraits                     */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class C>
struct NormTraits<MultiArrayView<N, T, C> >
{
    typedef MultiArrayView<N, T, C>                                      Type;
    typedef typename NormTraits<T>::SquaredNormType SquaredNormType;
    typedef typename SquareRootTraits<SquaredNormType>::SquareRootResult NormType;
};

template <unsigned int N, class T, class A>
struct NormTraits<MultiArray<N, T, A> >
: public NormTraits<typename MultiArray<N, T, A>::view_type>
{
    typedef NormTraits<typename MultiArray<N, T, A>::view_type>  BaseType;
    typedef MultiArray<N, T, A>                                  Type;
    typedef typename BaseType::SquaredNormType                   SquaredNormType;
    typedef typename BaseType::NormType                          NormType;
};

/** \brief Base class for, and view to, \ref vigra::MultiArray.

This class implements the interface of both MultiArray and
MultiArrayView.  By default, MultiArrayViews are tagged as
strided (using <tt>StridedArrayTag</tt> as third template parameter). 
This means that the array elements need not be consecutive in memory,
making the view flexible to represent all kinds of subarrays and slices.
In certain cases (which have become rare due to improvements of 
optimizer and processor technology), an array may be tagged with 
<tt>UnstridedArrayTag</tt> which indicates that the first array dimension
is guaranteed to be unstrided, i.e. has consecutive elements in memory.

In addition to the member functions described here, <tt>MultiArrayView</tt>
and its subclasses support arithmetic and algebraic functions via the 
module \ref MultiMathModule.

If you want to apply an algorithm requiring an image to a
<tt>MultiArrayView</tt> of appropriate (2-dimensional) shape, you can
create a \ref vigra::BasicImageView that acts as a wrapper with the
necessary interface -- see \ref MultiArrayToImage.

The template parameter are as follows
\code
    N: the array dimension

    T: the type of the array elements

    C: a tag determining whether the array's inner dimension is strided
       or not. An array is unstrided if the array elements occupy consecutive
       memory location, strided if there is an offset in between (e.g.
       when a view is created that skips every other array element).
       The compiler can generate faster code for unstrided arrays.
       Possible values: StridedArrayTag (default), UnstridedArrayTag
\endcode

<b>\#include</b> \<vigra/multi_array.hxx\> <br/>
Namespace: vigra
*/
template <unsigned int N, class T, class StrideTag>
class MultiArrayView
{
public:

        /** the array's actual dimensionality.
            This ensures that MultiArrayView can also be used for
            scalars (that is, when <tt>N == 0</tt>). Calculated as:<br>
            \code
            actual_dimension = (N==0) ? 1 : N
            \endcode
         */
    enum ActualDimension { actual_dimension = (N==0) ? 1 : N };

        /** the array's value type
         */
    typedef T value_type;

        /** reference type (result of operator[])
         */
    typedef value_type &reference;

        /** const reference type (result of operator[] const)
         */
    typedef const value_type &const_reference;

        /** pointer type
         */
    typedef value_type *pointer;

        /** const pointer type
         */
    typedef const value_type *const_pointer;

        /** difference type (used for multi-dimensional offsets and indices)
         */
    typedef typename MultiArrayShape<actual_dimension>::type difference_type;

        /** key type (argument of index operator array[i] -- same as difference_type)
         */
    typedef difference_type key_type;

        /** size type
         */
    typedef difference_type size_type;

        /** difference and index type for a single dimension
         */
    typedef MultiArrayIndex difference_type_1;

        /** scan-order iterator (StridedScanOrderIterator) type
         */
    typedef StridedScanOrderIterator<actual_dimension, T, T &, T *> iterator;

        /** const scan-order iterator (StridedScanOrderIterator) type
         */
    typedef StridedScanOrderIterator<actual_dimension, T, T const &, T const *> const_iterator;

        /** traverser (MultiIterator) type
         */
    typedef typename vigra::detail::MultiIteratorChooser <
        StrideTag>::template Traverser <actual_dimension, T, T &, T *>::type traverser;

        /** const traverser (MultiIterator) type
         */
    typedef typename vigra::detail::MultiIteratorChooser <
        StrideTag>::template Traverser <actual_dimension, T, T const &, T const *>::type const_traverser;

        /** the view type associated with this array.
         */
    typedef MultiArrayView <N, T, StrideTag> view_type;

        /** the matrix type associated with this array.
         */
    typedef MultiArray <N, T> matrix_type;

    bool checkInnerStride(UnstridedArrayTag) const
    {
        return m_stride[0] <= 1;
    }
    
    bool checkInnerStride(StridedArrayTag) const
    {
        return true;
    }

  protected:

    typedef typename difference_type::value_type diff_zero_t;

        /** the shape of the image pointed to is stored here.
        */
    difference_type m_shape;

        /** the strides (offset of a sample to the next) for every dimension
            are stored here.
        */
    difference_type m_stride;

        /** pointer to the image.
         */
    pointer m_ptr;

    template <class CN>
    void assignImpl(const MultiArrayView <N, T, CN>& rhs);

    template <class U, class CN>
    void copyImpl(const MultiArrayView <N, U, CN>& rhs);

    template <class U, class CN>
    void swapDataImpl(MultiArrayView <N, U, CN> rhs);

    template <class CN>
    bool arraysOverlap(const MultiArrayView <N, T, CN>& rhs) const;

    template <class U, class CN>
    bool arraysOverlap(const MultiArrayView <N, U, CN>&) const
    {
        return false;
    }

public:

        /** default constructor: create an invalid view,
         * i.e. hasData() returns false and size() is zero.
         */
    MultiArrayView ()
        : m_shape (diff_zero_t(0)), m_stride (diff_zero_t(0)), m_ptr (0)
    {}

        /** construct from another array view.
            Throws a precondition error if this array has UnstridedArrayTag, but the 
            innermost dimension of \a other is strided.
         */
    template <class Stride>
    MultiArrayView (const MultiArrayView<N, T, Stride> &other)
    : m_shape (other.shape()),
      m_stride (other.stride()),
      m_ptr (other.data())
    {
        vigra_precondition(other.checkInnerStride(StrideTag()),
            "MultiArrayView<..., UnstridedArrayTag>(MultiArrayView const &): cannot create unstrided view from strided array.");
    }

        /** construct from shape and pointer
         */
    MultiArrayView (const difference_type &shape, const_pointer ptr)
    : m_shape (shape),
      m_stride (detail::defaultStride<actual_dimension>(shape)),
      m_ptr (const_cast<pointer>(ptr))
    {}

        /** Construct from shape, strides (offset of a sample to the
            next) for every dimension, and pointer.  (Note that
            strides are not given in bytes, but in offset steps of the
            respective pointer type.)
         */
    MultiArrayView (const difference_type &shape,
                    const difference_type &stride,
                    const_pointer ptr)
    : m_shape (shape),
      m_stride (stride),
      m_ptr (const_cast<pointer>(ptr))
    {
        vigra_precondition(checkInnerStride(StrideTag()),
            "MultiArrayView<..., UnstridedArrayTag>::MultiArrayView(): First dimension of given array is not unstrided.");
    }
    
        /** Construct from an old-style BasicImage.
         */
    template <class ALLOC>
    MultiArrayView (BasicImage<T, ALLOC> const & image)
    : m_shape (Shape2(image.width(), image.height())),
      m_stride (detail::defaultStride<actual_dimension>(m_shape)),
      m_ptr (const_cast<pointer>(image.data()))
    {}
    
        /** Conversion to a strided view.
         */
    operator MultiArrayView<N, T, StridedArrayTag>() const
    {
        return MultiArrayView<N, T, StridedArrayTag>(m_shape, m_stride, m_ptr);
    }

	/** Reset this <tt>MultiArrayView</tt> to an invalid state (as after default construction).
		Can e.g. be used prior to assignment to make a view object point to new data.
         */
    void reset() {
	m_shape = diff_zero_t(0);
	m_stride = diff_zero_t(0);
	m_ptr = 0;
    }


        /** Assignment. There are 3 cases:

            <ul>
            <li> When this <tt>MultiArrayView</tt> does not point to valid data
                 (e.g. after default construction), it becomes a new view of \a rhs.
            <li> Otherwise, when the shapes of the two arrays match, the contents 
                 (i.e. the elements) of \a rhs are copied.
            <li> Otherwise, a <tt>PreconditionViolation</tt> exception is thrown.
            </ul>
         */
    MultiArrayView & operator=(MultiArrayView const & rhs)
    {
        if(this != &rhs)
            assignImpl(rhs);
        return *this;
    }

    template<class Stride2>
    MultiArrayView & operator=(MultiArrayView<N, T, Stride2> const & rhs)
    {
        assignImpl(rhs);
        return *this;
    }

        /** Assignment of a differently typed MultiArrayView. It copies the elements
            of\a rhs or fails with <tt>PreconditionViolation</tt> exception when 
            the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator=(MultiArrayView<N, U, C1> const & rhs)
    {
        vigra_precondition(this->shape() == rhs.shape(),
            "MultiArrayView::operator=(): shape mismatch.");
        this->copyImpl(rhs);
        return *this;
    }

        /** Assignment of a scalar. Equivalent to MultiArrayView::init(v).
         */
    MultiArrayView & operator=(value_type const & v)
    {
        return init(v);
    }

        /** Add-assignment of a compatible MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator+=(MultiArrayView<N, U, C1> const & rhs);

        /** Subtract-assignment of a compatible MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator-=(MultiArrayView<N, U, C1> const & rhs);

        /** Multiply-assignment of a compatible MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator*=(MultiArrayView<N, U, C1> const & rhs);

        /** Divide-assignment of a compatible MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator/=(MultiArrayView<N, U, C1> const & rhs);

        /** Add-assignment of a scalar.
         */
    MultiArrayView & operator+=(T const & rhs)
    {
        detail::copyAddScalarMultiArrayData(traverser_begin(), shape(), rhs, MetaInt<actual_dimension-1>());
        return *this;
    }

        /** Subtract-assignment of a scalar.
         */
    MultiArrayView & operator-=(T const & rhs)
    {
        detail::copySubScalarMultiArrayData(traverser_begin(), shape(), rhs, MetaInt<actual_dimension-1>());
        return *this;
    }

        /** Multiply-assignment of a scalar.
         */
    MultiArrayView & operator*=(T const & rhs)
    {
        detail::copyMulScalarMultiArrayData(traverser_begin(), shape(), rhs, MetaInt<actual_dimension-1>());
        return *this;
    }

        /** Divide-assignment of a scalar.
         */
    MultiArrayView & operator/=(T const & rhs)
    {
        detail::copyDivScalarMultiArrayData(traverser_begin(), shape(), rhs, MetaInt<actual_dimension-1>());
        return *this;
    }

        /** Assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::assign(*this, rhs);
        return *this;
    }

        /** Add-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator+=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::plusAssign(*this, rhs);
        return *this;
    }

        /** Subtract-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator-=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::minusAssign(*this, rhs);
        return *this;
    }

        /** Multiply-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator*=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::multiplyAssign(*this, rhs);
        return *this;
    }

        /** Divide-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator/=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::divideAssign(*this, rhs);
        return *this;
    }

        /** array access.
         */
    reference operator[] (const difference_type &d)
    {
        VIGRA_ASSERT_INSIDE(d);
        return m_ptr [dot (d, m_stride)];
    }

        /** array access.
         */
    const_reference operator[] (const difference_type &d) const
    {
        VIGRA_ASSERT_INSIDE(d);
        return m_ptr [dot (d, m_stride)];
    }

        /** equivalent to bindInner(), when M < N.
         */
    template <int M>
    MultiArrayView <N-M, T, StridedArrayTag> operator[] (const TinyVector<MultiArrayIndex, M> &d) const
    {
        return bindInner(d);
    }

        /** Array access in scan-order sense.
            Mostly useful to support standard indexing for 1-dimensional multi-arrays,
            but works for any N. Use scanOrderIndexToCoordinate() and
            coordinateToScanOrderIndex() for conversion between indices and coordinates.
            
            <b>Note:</b> This function should not be used in the inner loop, because the
            conversion of the scan order index into a memory address is expensive
            (it must take into account that memory may not be consecutive for subarrays
            and/or strided arrays). Always prefer operator() if possible.            
         */
    reference operator[](difference_type_1 d)
    {
        VIGRA_ASSERT_INSIDE(scanOrderIndexToCoordinate(d));
        return m_ptr [detail::ScanOrderToOffset<actual_dimension>::exec(d, m_shape, m_stride)];
    }

        /** Array access in scan-order sense.
            Mostly useful to support standard indexing for 1-dimensional multi-arrays,
            but works for any N. Use scanOrderIndexToCoordinate() and
            coordinateToScanOrderIndex() for conversion between indices and coordinates.
             
            <b>Note:</b> This function should not be used in the inner loop, because the
            conversion of the scan order index into a memory address is expensive
            (it must take into account that memory may not be consecutive for subarrays
            and/or strided arrays). Always prefer operator() if possible.            
        */
    const_reference operator[](difference_type_1 d) const
    {
        VIGRA_ASSERT_INSIDE(scanOrderIndexToCoordinate(d));
        return m_ptr [detail::ScanOrderToOffset<actual_dimension>::exec(d, m_shape, m_stride)];
    }

        /** convert scan-order index to coordinate.
         */
    difference_type scanOrderIndexToCoordinate(difference_type_1 d) const
    {
        difference_type result;
        detail::ScanOrderToCoordinate<actual_dimension>::exec(d, m_shape, result);
        return result;
    }

        /** convert coordinate to scan-order index.
         */
    difference_type_1 coordinateToScanOrderIndex(const difference_type &d) const
    {
        return detail::CoordinateToScanOrder<actual_dimension>::exec(m_shape, d);
    }

        /** 1D array access. Use only if N == 1.
         */
    reference operator() (difference_type_1 x)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x));
        return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x)];
    }

        /** 2D array access. Use only if N == 2.
         */
    reference operator() (difference_type_1 x, difference_type_1 y)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y));
        return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x, y)];
    }

        /** 3D array access. Use only if N == 3.
         */
    reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y, z));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    }

        /** 4D array access. Use only if N == 4.
         */
    reference operator() (difference_type_1 x, difference_type_1 y,
                          difference_type_1 z, difference_type_1 u)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y, z, u));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    }

        /** 5D array access. Use only if N == 5.
         */
    reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z,
                          difference_type_1 u, difference_type_1 v)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y,z, u,v));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    }

        /** 1D const array access. Use only if N == 1.
         */
    const_reference operator() (difference_type_1 x) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x));
        return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x)];
    }

        /** 2D const array access. Use only if N == 2.
         */
    const_reference operator() (difference_type_1 x, difference_type_1 y) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y));
        return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x, y)];
    }

        /** 3D const array access. Use only if N == 3.
         */
    const_reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x,y,z));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    }

        /** 4D const array access. Use only if N == 4.
         */
    const_reference operator() (difference_type_1 x, difference_type_1 y,
                                difference_type_1 z, difference_type_1 u) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x,y,z,u));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    }

        /** 5D const array access. Use only if N == 5.
         */
    const_reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z,
                                difference_type_1 u, difference_type_1 v) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x,y,z,u,v));
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    }

        /** Init with a constant.
         */
    template <class U>
    MultiArrayView & init(const U & init)
    {
        if(hasData())
            detail::copyScalarMultiArrayData(traverser_begin(), shape(), init, MetaInt<actual_dimension-1>());
        return *this;
    }


        /** Copy the data of the right-hand array (array shapes must match).
         */
    void copy(const MultiArrayView & rhs)
    {
        if(this == &rhs)
            return;
        this->copyImpl(rhs);
    }

        /** Copy the data of the right-hand array (array shapes must match).
         */
    template <class U, class CN>
    void copy(const MultiArrayView <N, U, CN>& rhs)
    {
        this->copyImpl(rhs);
    }

        /** Swap the pointers, shaes and strides between two array views.
        
            This function must be used with care. Never swap a MultiArray
            (which owns data) with a MultiArrayView:
            \code
                MultiArray<2, int> a(3,2), b(3,2);
                MultiArrayView<2, int> va(a);
                
                va.swap(b);   // danger!
            \endcode
            Now, <tt>a</tt> and <tt>b</tt> refer to the same memory. This may lead
            to a crash in their destructor, and in any case leaks <tt>b</tt>'s original 
            memory. Only use swap() on copied MultiArrayViews:
            \code
                MultiArray<2, int> a(3,2), b(3,2);
                MultiArrayView<2, int> va(a), vb(b);
                
                va.swap(vb);   // OK
            \endcode
         */
    void swap(MultiArrayView & other)
    {
        if (this == &other)
            return;
        std::swap(this->m_shape,  other.m_shape);
        std::swap(this->m_stride, other.m_stride);
        std::swap(this->m_ptr,    other.m_ptr);
    }

        /** swap the data between two MultiArrayView objects.

            The shapes of the two array must match.
        */
    void swapData(MultiArrayView rhs)
    {
        if(this != &rhs)
            swapDataImpl(rhs);
    }

        /** swap the data between two MultiArrayView objects.

            The shapes of the two array must match.
        */
    template <class T2, class C2>
    void swapData(MultiArrayView <N, T2, C2> rhs)
    {
        swapDataImpl(rhs);
    }
    
        /** check whether the array is unstrided (i.e. has consecutive memory) up 
            to the given dimension.

            \a dimension can range from 0 ... N-1. If a certain dimension is unstrided, 
            all lower dimensions are also unstrided.
        */
    bool isUnstrided(unsigned int dimension = N-1) const
    {
        difference_type s = vigra::detail::defaultStride<actual_dimension>(shape());
        for(unsigned int k = 0; k <= dimension; ++k)
            if(stride(k) != s[k])
                return false;
        return true;
    }

        /** bind the M outmost dimensions to certain indices.
            this reduces the dimensionality of the image to
            max { 1, N-M }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 1D array by fixing index 1 to 12, and index 2 to 10
            MultiArrayView <1, double> array1 = array3.bindOuter(TinyVector<MultiArrayIndex, 2>(12, 10));
            \endcode
        */
    template <int M, class Index>
    MultiArrayView <N-M, T, StrideTag> bindOuter(const TinyVector <Index, M> &d) const;

        /** bind the M innermost dimensions to certain indices.
            this reduces the dimensionality of the image to
            max { 1, N-M }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 1D array by fixing index 0 to 12, and index 1 to 10
            MultiArrayView <1, double, StridedArrayTag> array1 = array3.bindInner(TinyVector<MultiArrayIndex, 2>(12, 10));
            \endcode
        */
    template <int M, class Index>
    MultiArrayView <N-M, T, StridedArrayTag> bindInner(const TinyVector <Index, M> &d) const;

        /** bind dimension M to index d.
            this reduces the dimensionality of the image to
            max { 1, N-1 }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 2D array by fixing index 1 to 12
            MultiArrayView <2, double> array2 = array3.bind<1>(12);

            // get a 2D array by fixing index 0 to 23
            MultiArrayView <2, double, StridedArrayTag> array2a = array3.bind<0>(23);
            \endcode
         */
    template <unsigned int M>
    MultiArrayView <N-1, T, typename vigra::detail::MaybeStrided<StrideTag, M>::type >
    bind (difference_type_1 d) const;

        /** bind the outmost dimension to a certain index.
            this reduces the dimensionality of the image to
            max { 1, N-1 }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 2D array by fixing the outermost index (i.e. index 2) to 12
            MultiArrayView <2, double> array2 = array3.bindOuter(12);
            \endcode
        */
    MultiArrayView <N-1, T, StrideTag> bindOuter (difference_type_1 d) const;

        /** bind the innermost dimension to a certain index.
            this reduces the dimensionality of the image to
            max { 1, N-1 }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 2D array by fixing the innermost index (i.e. index 0) to 23
            MultiArrayView <2, double, StridedArrayTag> array2 = array3.bindInner(23);
            \endcode
        */
    MultiArrayView <N-1, T, StridedArrayTag> bindInner (difference_type_1 d) const;

        /** bind dimension m to index d.
            this reduces the dimensionality of the image to
            max { 1, N-1 }.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a 2D array by fixing index 2 to 15
            MultiArrayView <2, double, StridedArrayTag> array2 = array3.bindAt(2, 15);
            \endcode
         */
    MultiArrayView <N-1, T, StridedArrayTag>
    bindAt (difference_type_1 m, difference_type_1 d) const;
    
        /** Create a view to channel 'i' of a vector-like value type. Possible value types
            (of the original array) are: \ref TinyVector, \ref RGBValue, \ref FFTWComplex, 
            and <tt>std::complex</tt>. The list can be extended to any type whose memory
            layout is equivalent to a fixed-size C array, by specializing 
            <tt>ExpandElementResult</tt>.

            <b>Usage:</b>
            \code
                MultiArray<2, RGBValue<float> > rgb_image(Shape2(w, h));
                
                MultiArrayView<2, float, StridedArrayTag> red   = rgb_image.bindElementChannel(0);
                MultiArrayView<2, float, StridedArrayTag> green = rgb_image.bindElementChannel(1);
                MultiArrayView<2, float, StridedArrayTag> blue  = rgb_image.bindElementChannel(2);
            \endcode
        */
    MultiArrayView <N, typename ExpandElementResult<T>::type, StridedArrayTag> 
    bindElementChannel(difference_type_1 i) const
    {
        vigra_precondition(0 <= i && i < ExpandElementResult<T>::size,
              "MultiArrayView::bindElementChannel(i): 'i' out of range.");
        return expandElements(0).bindInner(i);
    }

        /** Create a view where a vector-like element type is expanded into a new 
            array dimension. The new dimension is inserted at index position 'd',
            which must be between 0 and N inclusive.
            
            Possible value types of the original array are: \ref TinyVector, \ref RGBValue, 
            \ref FFTWComplex, <tt>std::complex</tt>, and the built-in number types (in this 
            case, <tt>expandElements</tt> is equivalent to <tt>insertSingletonDimension</tt>). 
            The list of supported types can be extended to any type whose memory
            layout is equivalent to a fixed-size C array, by specializing 
            <tt>ExpandElementResult</tt>.

            <b>Usage:</b>
            \code
                MultiArray<2, RGBValue<float> > rgb_image(Shape2(w, h));
                
                MultiArrayView<3, float, StridedArrayTag> multiband_image = rgb_image.expandElements(2);
            \endcode
        */
    MultiArrayView <N+1, typename ExpandElementResult<T>::type, StridedArrayTag> 
    expandElements(difference_type_1 d) const;
    
        /** Add a singleton dimension (dimension of length 1).

            Singleton dimensions don't change the size of the data, but introduce
            a new index that can only take the value 0. This is mainly useful for
            the 'reduce mode' of transformMultiArray() and combineTwoMultiArrays(),
            because these functions require the source and destination arrays to
            have the same number of dimensions.

            The range of \a i must be <tt>0 <= i <= N</tt>. The new dimension will become
            the i'th index, and the old indices from i upwards will shift one
            place to the right.

            <b>Usage:</b>

            Suppose we want have a 2D array and want to create a 1D array that contains
            the row average of the first array.
            \code
            typedef MultiArrayShape<2>::type Shape2;
            MultiArray<2, double> original(Shape2(40, 30));

            typedef MultiArrayShape<1>::type Shape1;
            MultiArray<1, double> rowAverages(Shape1(30));

            // temporarily add a singleton dimension to the destination array
            transformMultiArray(srcMultiArrayRange(original),
                                destMultiArrayRange(rowAverages.insertSingletonDimension(0)),
                                FindAverage<double>());
            \endcode
         */
    MultiArrayView <N+1, T, StrideTag>
    insertSingletonDimension (difference_type_1 i) const;
    
        /** create a multiband view for this array.

            The type <tt>MultiArrayView<N, Multiband<T> ></tt> tells VIGRA
            algorithms which recognize the <tt>Multiband</tt> modifier to
            interpret the outermost (last) dimension as a channel dimension. 
            In effect, these algorithms will treat the data as a set of 
            (N-1)-dimensional arrays instead of a single N-dimensional array.
        */
    MultiArrayView<N, Multiband<value_type>, StrideTag> multiband() const
    {
        return MultiArrayView<N, Multiband<value_type>, StrideTag>(*this);
    }

        /** Create a view to the diagonal elements of the array.
        
            This produces a 1D array view whose size equals the size
            of the shortest dimension of the original array.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a view to the diagonal elements
            MultiArrayView <1, double, StridedArrayTag> diagonal = array3.diagonal();
            assert(diagonal.shape(0) == 20);
            \endcode
        */
    MultiArrayView<1, T, StridedArrayTag> diagonal() const
    {
        return MultiArrayView<1, T, StridedArrayTag>(Shape1(vigra::min(m_shape)), 
                                                     Shape1(vigra::sum(m_stride)), m_ptr);
    }

        /** create a rectangular subarray that spans between the
            points p and q, where p is in the subarray, q not. 
            If an element of p or q is negative, it is subtracted
            from the correspongng shape.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a subarray set is smaller by one element at all sides
            MultiArrayView <3, double> subarray  = array3.subarray(Shape(1,1,1), Shape(39, 29, 19));
            
            // specifying the end point with a vector of '-1' is equivalent
            MultiArrayView <3, double> subarray2 = array3.subarray(Shape(1,1,1), Shape(-1, -1, -1));
            \endcode
        */
    MultiArrayView subarray (difference_type p, difference_type q) const
    {
        detail::RelativeToAbsoluteCoordinate<actual_dimension-1>::exec(shape(), p);
        detail::RelativeToAbsoluteCoordinate<actual_dimension-1>::exec(shape(), q);
        const difference_type_1 offset = dot (m_stride, p);
        return MultiArrayView (q - p, m_stride, m_ptr + offset);
    }

        /** apply an additional striding to the image, thereby reducing
            the shape of the array.
            for example, multiplying the stride of dimension one by three
            turns an appropriately laid out (interleaved) rgb image into
            a single band image.
        */
    MultiArrayView <N, T, StridedArrayTag>
    stridearray (const difference_type &s) const
    {
        difference_type shape = m_shape;
        for (unsigned int i = 0; i < actual_dimension; ++i)
            shape[i] = (shape[i] + s[i] - 1) / s[i];
        return MultiArrayView <N, T, StridedArrayTag>(shape, m_stride * s, m_ptr);
    }

        /** Transpose an array. If N==2, this implements the usual matrix transposition.
            For N > 2, it reverses the order of the indices.

            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);

            MultiArrayView<2, double, StridedArrayTag> transposed = array.transpose();

            for(int i=0; i<array.shape(0), ++i)
                for(int j=0; j<array.shape(1); ++j)
                    assert(array(i, j) == transposed(j, i));
            \endcode
        */
    MultiArrayView <N, T, StridedArrayTag>
    transpose () const
    {
        difference_type shape(m_shape.begin(), difference_type::ReverseCopy),
                        stride(m_stride.begin(), difference_type::ReverseCopy);
        return MultiArrayView <N, T, StridedArrayTag>(shape, stride, m_ptr);
    }

        /** Permute the dimensions of the array.
            The function exchanges the orer of the array's axes without copying the data.
            Argument\a permutation specifies the desired order such that 
            <tt>permutation[k] = j</tt> means that axis <tt>j</tt> in the original array
            becomes axis <tt>k</tt> in the transposed array. 

            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);

            MultiArrayView<2, double, StridedArrayTag> transposed = array.transpose(Shape(1,0));

            for(int i=0; i<array.shape(0), ++i)
                for(int j=0; j<array.shape(1); ++j)
                    assert(array(i, j) == transposed(j, i));
            \endcode
        */
    MultiArrayView <N, T, StridedArrayTag>
    transpose(const difference_type &permutation) const
    {
        return permuteDimensions(permutation);
    }

    MultiArrayView <N, T, StridedArrayTag>
    permuteDimensions (const difference_type &s) const;

        /** Permute the dimensions of the array so that the strides are in ascending order.
            Determines the appropriate permutation and then calls permuteDimensions().
        */
    MultiArrayView <N, T, StridedArrayTag>
    permuteStridesAscending() const;
    
        /** Permute the dimensions of the array so that the strides are in descending order.
            Determines the appropriate permutation and then calls permuteDimensions().
        */
    MultiArrayView <N, T, StridedArrayTag>
    permuteStridesDescending() const;
    
        /** Compute the ordering of the strides in this array.
            The result is describes the current permutation of the axes relative 
            to the standard ascending stride order.
        */
    difference_type strideOrdering() const
    {
        return strideOrdering(m_stride);
    }
    
        /** Compute the ordering of the given strides.
            The result is describes the current permutation of the axes relative 
            to the standard ascending stride order.
        */
    static difference_type strideOrdering(difference_type strides);

        /** number of the elements in the array.
         */
    difference_type_1 elementCount () const
    {
        difference_type_1 ret = m_shape[0];
        for(int i = 1; i < actual_dimension; ++i)
            ret *= m_shape[i];
        return ret;
    }

        /** number of the elements in the array.
            Same as <tt>elementCount()</tt>. Mostly useful to support the std::vector interface.
         */
    difference_type_1 size () const
    {
        return elementCount();
    }

        /** return the array's shape.
         */
    const difference_type & shape () const
    {
        return m_shape;
    }

        /** return the array's size at a certain dimension.
         */
    difference_type_1 size (difference_type_1 n) const
    {
        return m_shape [n];
    }

        /** return the array's shape at a certain dimension
            (same as <tt>size(n)</tt>).
         */
    difference_type_1 shape (difference_type_1 n) const
    {
        return m_shape [n];
    }

        /** return the array's width (same as <tt>shape(0)</tt>).
         */
    difference_type_1 width() const
    {
        return m_shape [0];
    }

        /** return the array's height (same as <tt>shape(1)</tt>).
         */
    difference_type_1 height() const
    {
        return m_shape [1];
    }

        /** return the array's stride for every dimension.
         */
    const difference_type & stride () const
    {
        return m_stride;
    }

        /** return the array's stride at a certain dimension.
         */
    difference_type_1 stride (int n) const
    {
        return m_stride [n];
    }

        /** check whether two arrays are elementwise equal.
         */
    template <class U, class C1>
    bool operator==(MultiArrayView<N, U, C1> const & rhs) const
    {
        if(this->shape() != rhs.shape())
            return false;
        return detail::equalityOfMultiArrays(traverser_begin(), shape(), rhs.traverser_begin(), MetaInt<actual_dimension-1>());
    }

        /** check whether two arrays are not elementwise equal.
            Also true when the two arrays have different shapes.
         */
    template <class U, class C1>
    bool operator!=(MultiArrayView<N, U, C1> const & rhs) const
    {
        return !operator==(rhs);
    }

        /** check whether the given point is in the array range.
         */
    bool isInside (difference_type const & p) const
    {
        for(int d=0; d<actual_dimension; ++d)
            if(p[d] < 0 || p[d] >= shape(d))
                return false;
        return true;
    }
        /** check whether the given point is not in the array range.
         */
    bool isOutside (difference_type const & p) const
    {
        for(int d=0; d<actual_dimension; ++d)
            if(p[d] < 0 || p[d] >= shape(d))
                return true;
        return false;
    }

        /** Check if the array contains only non-zero elements (or if all elements
            are 'true' if the value type is 'bool').
         */
    bool all() const
    {
        bool res = true;
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::AllTrueReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
    }

        /** Check if the array contains a non-zero element (or an element
            that is 'true' if the value type is 'bool').
         */
    bool any() const
    {
        bool res = false;
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::AnyTrueReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
    }

        /** Find the minimum and maximum element in this array. 
            See \ref FeatureAccumulators for a general feature 
            extraction framework.
         */
    void minmax(T * minimum, T * maximum) const
    {
        std::pair<T, T> res(NumericTraits<T>::max(), NumericTraits<T>::min());
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::MinmaxReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        *minimum = res.first;
        *maximum = res.second;
    }

        /** Compute the mean and variance of the values in this array. 
            See \ref FeatureAccumulators for a general feature 
            extraction framework.
         */
    template <class U>
    void meanVariance(U * mean, U * variance) const
    {
        typedef typename NumericTraits<U>::RealPromote R;
        R zero = R();
        triple<double, R, R> res(0.0, zero, zero);
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::MeanVarianceReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        *mean     = res.second;
        *variance = res.third / res.first;
    }

        /** Compute the sum of the array elements.

            You must provide the type of the result by an explicit template parameter:
            \code
            MultiArray<2, UInt8> A(width, height);
            
            double sum = A.sum<double>();
            \endcode
         */
    template <class U>
    U sum() const
    {
        U res = NumericTraits<U>::zero();
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::SumReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
    }

        /** Compute the sum of the array elements over selected axes.
        
            \arg sums must have the same shape as this array, except for the
            axes along which the sum is to be accumulated. These axes must be 
            singletons. Note that you must include <tt>multi_pointoperators.hxx</tt>
            for this function to work.

            <b>Usage:</b>
            \code
            #include <vigra/multi_array.hxx>
            #include <vigra/multi_pointoperators.hxx>
            
            MultiArray<2, double> A(Shape2(rows, cols));
            ... // fill A
            
            // make the first axis a singleton to sum over the first index
            MultiArray<2, double> rowSums(Shape2(1, cols));
            A.sum(rowSums);
            
            // this is equivalent to
            transformMultiArray(srcMultiArrayRange(A),
                                destMultiArrayRange(rowSums),
                                FindSum<double>());
            \endcode
         */
    template <class U, class S>
    void sum(MultiArrayView<N, U, S> sums) const
    {
        transformMultiArray(srcMultiArrayRange(*this),
                            destMultiArrayRange(sums),
                            FindSum<U>());
    }

        /** Compute the product of the array elements.

            You must provide the type of the result by an explicit template parameter:
            \code
            MultiArray<2, UInt8> A(width, height);
            
            double prod = A.product<double>();
            \endcode
         */
    template <class U>
    U product() const
    {
        U res = NumericTraits<U>::one();
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::ProdReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
    }

        /** Compute the squared Euclidean norm of the array (sum of squares of the array elements).
         */
    typename NormTraits<MultiArrayView>::SquaredNormType 
    squaredNorm() const
    {
        typedef typename NormTraits<MultiArrayView>::SquaredNormType SquaredNormType;
        SquaredNormType res = NumericTraits<SquaredNormType>::zero();
        detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     res, 
                                     detail::SquaredL2NormReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
    }

        /** Compute various norms of the array.
            The norm is determined by parameter \a type:

            <ul>
            <li> type == 0: maximum norm (L-infinity): maximum of absolute values of the array elements
            <li> type == 1: Manhattan norm (L1): sum of absolute values of the array elements
            <li> type == 2: Euclidean norm (L2): square root of <tt>squaredNorm()</tt> when \a useSquaredNorm is <tt>true</tt>,<br>
                 or direct algorithm that avoids underflow/overflow otherwise.
            </ul>

            Parameter \a useSquaredNorm has no effect when \a type != 2. Defaults: compute L2 norm as square root of
            <tt>squaredNorm()</tt>.
         */
    typename NormTraits<MultiArrayView>::NormType 
    norm(int type = 2, bool useSquaredNorm = true) const;

        /** return the pointer to the image data
         */
    pointer data () const
    {
        return m_ptr;
    }
    
    pointer & unsafePtr()
    {
        return m_ptr;
    }

        /**
         * returns true iff this view refers to valid data,
         * i.e. data() is not a NULL pointer.  (this is false after
         * default construction.)
         */
    bool hasData () const
    {
        return m_ptr != 0;
    }

        /** returns a scan-order iterator pointing
            to the first array element.
        */
    iterator begin()
    {
        return iterator(*this);
    }

        /** returns a const scan-order iterator pointing
            to the first array element.
        */
    const_iterator begin() const
    {
        return const_iterator(*this);
    }

        /** returns a scan-order iterator pointing
            beyond the last array element.
        */
    iterator end()
    {
        return begin().getEndIterator();
    }

        /** returns a const scan-order iterator pointing
            beyond the last array element.
        */
    const_iterator end() const
    {
        return begin().getEndIterator();
    }

        /** returns the N-dimensional MultiIterator pointing
            to the first element in every dimension.
        */
    traverser traverser_begin ()
    {
        traverser ret (m_ptr, m_stride.begin (), m_shape.begin ());
        return ret;
    }

        /** returns the N-dimensional MultiIterator pointing
            to the const first element in every dimension.
        */
    const_traverser traverser_begin () const
    {
        const_traverser ret (m_ptr, m_stride.begin (), m_shape.begin ());
        return ret;
    }

        /** returns the N-dimensional MultiIterator pointing
            beyond the last element in dimension N, and to the
            first element in every other dimension.
        */
    traverser traverser_end ()
    {
        traverser ret (m_ptr, m_stride.begin (), m_shape.begin ());
        ret += m_shape [actual_dimension-1];
        return ret;
    }

        /** returns the N-dimensional const MultiIterator pointing
            beyond the last element in dimension N, and to the
            first element in every other dimension.
        */
    const_traverser traverser_end () const
    {
        const_traverser ret (m_ptr, m_stride.begin (), m_shape.begin ());
        ret += m_shape [actual_dimension-1];
        return ret;
    }

    view_type view () const
    {
        return *this;
    }
};

template <unsigned int N, class T, class StrideTag>
class MultiArrayView<N, Multiband<T>, StrideTag>
: public MultiArrayView<N, T, StrideTag>
{
  public:
    MultiArrayView(MultiArrayView<N, T, StrideTag> const & v)
    : MultiArrayView<N, T, StrideTag>(v)
    {}
};


template <unsigned int N, class T, class Stride1>
template <class Stride2>
void
MultiArrayView <N, T, Stride1>::assignImpl(MultiArrayView<N, T, Stride2> const & rhs)
{
    if(m_ptr == 0)
    {
        vigra_precondition(rhs.checkInnerStride(Stride1()),
            "MultiArrayView<..., UnstridedArrayTag>::operator=(MultiArrayView const &): cannot create unstrided view from strided array.");
                           
        m_shape  = rhs.shape();
        m_stride = rhs.stride();
        m_ptr    = rhs.data();
    }
    else
    {
        vigra_precondition(this->shape() == rhs.shape(),
            "MultiArrayView::operator=(MultiArrayView const &): shape mismatch.");
        this->copyImpl(rhs);
    }
}

template <unsigned int N, class T, class StrideTag>
template <class CN>
bool
MultiArrayView <N, T, StrideTag>::arraysOverlap(const MultiArrayView <N, T, CN>& rhs) const
{
    vigra_precondition (shape () == rhs.shape (),
        "MultiArrayView::arraysOverlap(): shape mismatch.");
    const_pointer first_element = this->m_ptr,
                  last_element = first_element + dot(this->m_shape - difference_type(1), this->m_stride);
    typename MultiArrayView <N, T, CN>::const_pointer
           rhs_first_element = rhs.data(),
           rhs_last_element = rhs_first_element + dot(rhs.shape() - difference_type(1), rhs.stride());
    return !(last_element < rhs_first_element || rhs_last_element < first_element);
}

template <unsigned int N, class T, class StrideTag>
template <class U, class CN>
void
MultiArrayView <N, T, StrideTag>::copyImpl(const MultiArrayView <N, U, CN>& rhs)
{
    if(!arraysOverlap(rhs))
    {
        // no overlap -- can copy directly
        detail::copyMultiArrayData(rhs.traverser_begin(), shape(), traverser_begin(), MetaInt<actual_dimension-1>());
    }
    else
    {
        // overlap: we got different views to the same data -- copy to intermediate memory in order to avoid
        // overwriting elements that are still needed on the rhs.
        MultiArray<N, T> tmp(rhs);
        detail::copyMultiArrayData(tmp.traverser_begin(), shape(), traverser_begin(), MetaInt<actual_dimension-1>());
    }
}

#define VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT(name, op) \
template <unsigned int N, class T, class StrideTag> \
template<class U, class C1> \
MultiArrayView<N, T, StrideTag> &  \
MultiArrayView <N, T, StrideTag>::operator op(MultiArrayView<N, U, C1> const & rhs) \
{ \
    vigra_precondition(this->shape() == rhs.shape(), "MultiArrayView::operator" #op "() size mismatch."); \
    if(!arraysOverlap(rhs)) \
    { \
        detail::name##MultiArrayData(rhs.traverser_begin(), shape(), traverser_begin(), MetaInt<actual_dimension-1>()); \
    } \
    else \
    { \
        MultiArray<N, T> tmp(rhs); \
        detail::name##MultiArrayData(tmp.traverser_begin(), shape(), traverser_begin(), MetaInt<actual_dimension-1>()); \
    } \
    return *this; \
}

VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT(copyAdd, +=)
VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT(copySub, -=)
VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT(copyMul, *=)
VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT(copyDiv, /=)

#undef VIGRA_MULTI_ARRAY_COMPUTED_ASSIGNMENT

template <unsigned int N, class T, class StrideTag>
template <class U, class CN>
void
MultiArrayView <N, T, StrideTag>::swapDataImpl(MultiArrayView <N, U, CN> rhs)
{
    vigra_precondition (shape () == rhs.shape (),
        "MultiArrayView::swapData(): shape mismatch.");

    // check for overlap of this and rhs
    const_pointer first_element = this->m_ptr,
                  last_element = first_element + dot(this->m_shape - difference_type(1), this->m_stride);
    typename MultiArrayView <N, U, CN>::const_pointer
           rhs_first_element = rhs.data(),
           rhs_last_element = rhs_first_element + dot(rhs.shape() - difference_type(1), rhs.stride());
    if(last_element < rhs_first_element || rhs_last_element < first_element)
    {
        // no overlap -- can swap directly
        detail::swapDataImpl(traverser_begin(), shape(), rhs.traverser_begin(), MetaInt<actual_dimension-1>());
    }
    else
    {
        // overlap: we got different views to the same data -- copy to intermediate memory in order to avoid
        // overwriting elements that are still needed.
        MultiArray<N, T> tmp(*this);
        copy(rhs);
        rhs.copy(tmp);
    }
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::permuteDimensions (const difference_type &s) const
{
    difference_type shape, stride, check((typename difference_type::value_type)0);
    for (unsigned int i = 0; i < actual_dimension; ++i)
    {
        shape[i]  = m_shape[s[i]];
        stride[i] = m_stride[s[i]];
        ++check[s[i]];
    }
    vigra_precondition(check == difference_type(1),
       "MultiArrayView::transpose(): every dimension must occur exactly once.");
    return MultiArrayView <N, T, StridedArrayTag>(shape, stride, m_ptr);
}

template <unsigned int N, class T, class StrideTag>
typename MultiArrayView <N, T, StrideTag>::difference_type 
MultiArrayView <N, T, StrideTag>::strideOrdering(difference_type stride)
{
    difference_type permutation;
    for(int k=0; k<(int)N; ++k)
        permutation[k] = k;
    for(int k=0; k<(int)N-1; ++k)
    {
        int smallest = k;
        for(int j=k+1; j<(int)N; ++j)
        {
            if(stride[j] < stride[smallest])
                smallest = j;
        }
        if(smallest != k)
        {
            std::swap(stride[k], stride[smallest]);
            std::swap(permutation[k], permutation[smallest]);
        }
    }
    difference_type ordering;
    for(unsigned int k=0; k<N; ++k)
        ordering[permutation[k]] = k;
    return ordering;
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::permuteStridesAscending() const
{
    difference_type ordering(strideOrdering(m_stride)), permutation;
    for(MultiArrayIndex k=0; k<N; ++k)
        permutation[ordering[k]] = k;
    return permuteDimensions(permutation);
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::permuteStridesDescending() const
{
    difference_type ordering(strideOrdering(m_stride)), permutation;
    for(MultiArrayIndex k=0; k<N; ++k)
        permutation[N-1-ordering[k]] = k;
    return permuteDimensions(permutation);
}

template <unsigned int N, class T, class StrideTag>
template <int M, class Index>
MultiArrayView <N-M, T, StrideTag>
MultiArrayView <N, T, StrideTag>::bindOuter (const TinyVector <Index, M> &d) const
{
    TinyVector <MultiArrayIndex, M> stride;
    stride.init (m_stride.begin () + N-M, m_stride.end ());
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <MultiArrayIndex, NNew> inner_shape, inner_stride;
    if (N-M == 0)
    {
        inner_shape [0] = 1;
        inner_stride [0] = 1;
    }
    else
    {
        inner_shape.init (m_shape.begin (), m_shape.end () - M);
        inner_stride.init (m_stride.begin (), m_stride.end () - M);
    }
    return MultiArrayView <N-M, T, StrideTag> (inner_shape, inner_stride, ptr);
}

template <unsigned int N, class T, class StrideTag>
template <int M, class Index>
MultiArrayView <N - M, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::bindInner (const TinyVector <Index, M> &d) const
{
    TinyVector <MultiArrayIndex, M> stride;
    stride.init (m_stride.begin (), m_stride.end () - N + M);
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <MultiArrayIndex, NNew> outer_shape, outer_stride;
    if (N-M == 0)
    {
        outer_shape [0] = 1;
        outer_stride [0] = 1;
    }
    else
    {
        outer_shape.init (m_shape.begin () + M, m_shape.end ());
        outer_stride.init (m_stride.begin () + M, m_stride.end ());
    }
    return MultiArrayView <N-M, T, StridedArrayTag>
        (outer_shape, outer_stride, ptr);
}

template <unsigned int N, class T, class StrideTag>
template <unsigned int M>
MultiArrayView <N-1, T, typename detail::MaybeStrided<StrideTag, M>::type >
MultiArrayView <N, T, StrideTag>::bind (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> shape, stride;
    // the remaining dimensions are 0..n-1,n+1..N-1
    if (N-1 == 0)
    {
        shape[0] = 1;
        stride[0] = 1;
    }
    else
    {
        std::copy (m_shape.begin (), m_shape.begin () + M, shape.begin ());
        std::copy (m_shape.begin () + M+1, m_shape.end (),
                   shape.begin () + M);
        std::copy (m_stride.begin (), m_stride.begin () + M, stride.begin ());
        std::copy (m_stride.begin () + M+1, m_stride.end (),
                   stride.begin () + M);
    }
    return MultiArrayView <N-1, T, typename detail::MaybeStrided<StrideTag, M>::type>
        (shape, stride, m_ptr + d * m_stride[M]);
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N - 1, T, StrideTag>
MultiArrayView <N, T, StrideTag>::bindOuter (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> inner_shape, inner_stride;
    if (N-1 == 0)
    {
        inner_shape [0] = 1;
        inner_stride [0] = 1;
    }
    else
    {
        inner_shape.init (m_shape.begin (), m_shape.end () - 1);
        inner_stride.init (m_stride.begin (), m_stride.end () - 1);
    }
    return MultiArrayView <N-1, T, StrideTag> (inner_shape, inner_stride,
                                       m_ptr + d * m_stride [N-1]);
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N - 1, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::bindInner (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> outer_shape, outer_stride;
    if (N-1 == 0)
    {
        outer_shape [0] = 1;
        outer_stride [0] = 1;
    }
    else
    {
        outer_shape.init (m_shape.begin () + 1, m_shape.end ());
        outer_stride.init (m_stride.begin () + 1, m_stride.end ());
    }
    return MultiArrayView <N-1, T, StridedArrayTag>
        (outer_shape, outer_stride, m_ptr + d * m_stride [0]);
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N - 1, T, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::bindAt (difference_type_1 n, difference_type_1 d) const
{
    vigra_precondition (
        n < static_cast <int> (N),
        "MultiArrayView <N, T, StrideTag>::bindAt(): dimension out of range.");
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> shape, stride;
    // the remaining dimensions are 0..n-1,n+1..N-1
    if (N-1 == 0)
    {
        shape [0] = 1;
        stride [0] = 1;
    }
    else
    {
        std::copy (m_shape.begin (), m_shape.begin () + n, shape.begin ());
        std::copy (m_shape.begin () + n+1, m_shape.end (),
                   shape.begin () + n);
        std::copy (m_stride.begin (), m_stride.begin () + n, stride.begin ());
        std::copy (m_stride.begin () + n+1, m_stride.end (),
                   stride.begin () + n);
    }
    return MultiArrayView <N-1, T, StridedArrayTag>
        (shape, stride, m_ptr + d * m_stride[n]);
}


template <unsigned int N, class T, class StrideTag>
MultiArrayView <N+1, typename ExpandElementResult<T>::type, StridedArrayTag>
MultiArrayView <N, T, StrideTag>::expandElements(difference_type_1 d) const
{
    vigra_precondition(0 <= d && d <= static_cast <difference_type_1> (N),
          "MultiArrayView<N, ...>::expandElements(d): 0 <= 'd' <= N required.");
    
    int elementSize = ExpandElementResult<T>::size;
    typename MultiArrayShape<N+1>::type newShape, newStrides;
    for(int k=0; k<d; ++k)
    {
        newShape[k] = m_shape[k];
        newStrides[k] = m_stride[k]*elementSize;
    }   
    
    newShape[d] = elementSize;
    newStrides[d] = 1;
    
    for(int k=d; k<N; ++k)
    {
        newShape[k+1] = m_shape[k];
        newStrides[k+1] = m_stride[k]*elementSize;
    }   
    
    typedef typename ExpandElementResult<T>::type U;     
    return MultiArrayView<N+1, U, StridedArrayTag>(
                    newShape, newStrides, reinterpret_cast<U*>(m_ptr));
}

template <unsigned int N, class T, class StrideTag>
MultiArrayView <N+1, T, StrideTag>
MultiArrayView <N, T, StrideTag>::insertSingletonDimension (difference_type_1 i) const
{
    vigra_precondition (
        0 <= i && i <= static_cast <difference_type_1> (N),
        "MultiArrayView <N, T, StrideTag>::insertSingletonDimension(): index out of range.");
    TinyVector <MultiArrayIndex, N+1> shape, stride;
    std::copy (m_shape.begin (), m_shape.begin () + i, shape.begin ());
    std::copy (m_shape.begin () + i, m_shape.end (), shape.begin () + i + 1);
    std::copy (m_stride.begin (), m_stride.begin () + i, stride.begin ());
    std::copy (m_stride.begin () + i, m_stride.end (), stride.begin () + i + 1);
    shape[i] = 1;
    stride[i] = 1;

    return MultiArrayView <N+1, T, StrideTag>(shape, stride, m_ptr);
}

template <unsigned int N, class T, class StrideTag>
typename NormTraits<MultiArrayView <N, T, StrideTag> >::NormType
MultiArrayView <N, T, StrideTag>::norm(int type, bool useSquaredNorm) const
{
    typedef typename NormTraits<MultiArrayView>::NormType NormType;

    switch(type)
    {
      case 0:
      {
        NormType res = NumericTraits<NormType>::zero();
        detail::reduceOverMultiArray(traverser_begin(), shape(), 
                                     res, 
                                     detail::MaxNormReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
      }
      case 1:
      {
        NormType res = NumericTraits<NormType>::zero();
        detail::reduceOverMultiArray(traverser_begin(), shape(), 
                                     res, 
                                     detail::L1NormReduceFunctor(),
                                     MetaInt<actual_dimension-1>());
        return res;
      }
      case 2:
      {
        if(useSquaredNorm)
        {
            return sqrt((NormType)squaredNorm());
        }
        else
        {
            NormType normMax = NumericTraits<NormType>::zero();
            detail::reduceOverMultiArray(traverser_begin(), shape(), 
                                        normMax, 
                                        detail::MaxNormReduceFunctor(),
                                        MetaInt<actual_dimension-1>());
            if(normMax == NumericTraits<NormType>::zero())
                return normMax;
            NormType res  = NumericTraits<NormType>::zero();
            detail::reduceOverMultiArray(traverser_begin(), shape(), 
                                         res, 
                                         detail::WeightedL2NormReduceFunctor<NormType>(1.0/normMax),
                                         MetaInt<actual_dimension-1>());
            return sqrt(res)*normMax;
        }
      }
      default:
        vigra_precondition(false, "MultiArrayView::norm(): Unknown norm type.");
        return NumericTraits<NormType>::zero(); // unreachable
    }
}


/********************************************************/
/*                                                      */
/*                          norm                        */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class StrideTag>
inline typename NormTraits<MultiArrayView <N, T, StrideTag> >::SquaredNormType
squaredNorm(MultiArrayView <N, T, StrideTag> const & a)
{
    return a.squaredNorm();
}

template <unsigned int N, class T, class StrideTag>
inline typename NormTraits<MultiArrayView <N, T, StrideTag> >::NormType
norm(MultiArrayView <N, T, StrideTag> const & a)
{
    return a.norm();
}

/********************************************************/
/*                                                      */
/*                       MultiArray                     */
/*                                                      */
/********************************************************/

/** \brief Main <TT>MultiArray</TT> class containing the memory
    management.

This class inherits the interface of MultiArrayView, and implements
the memory ownership.
MultiArray's are always unstrided, striding them creates a MultiArrayView.


The template parameters are as follows
\code
    N: the array dimension

    T: the type of the array elements

    A: the allocator used for internal storage management
       (default: std::allocator<T>)
\endcode

<b>\#include</b> \<vigra/multi_array.hxx\> <br/>
Namespace: vigra
*/
template <unsigned int N, class T, class A /* default already declared above */>
class MultiArray 
: public MultiArrayView <N, typename vigra::detail::ResolveMultiband<T>::type, 
                            typename vigra::detail::ResolveMultiband<T>::Stride>
{
  public:
    typedef typename vigra::detail::ResolveMultiband<T>::Stride actual_stride;

        /** the view type associated with this array.
         */
    typedef MultiArrayView <N, typename vigra::detail::ResolveMultiband<T>::type, 
                               typename vigra::detail::ResolveMultiband<T>::Stride> view_type;
    
    using view_type::actual_dimension;

        /** the allocator type used to allocate the memory
         */
    typedef A allocator_type;

        /** the matrix type associated with this array.
         */
    typedef MultiArray <N, T, A> matrix_type;

        /** the array's value type
         */
    typedef typename view_type::value_type value_type;

        /** pointer type
         */
    typedef typename view_type::pointer pointer;

        /** const pointer type
         */
    typedef typename view_type::const_pointer const_pointer;

        /** reference type (result of operator[])
         */
    typedef typename view_type::reference reference;

        /** const reference type (result of operator[] const)
         */
    typedef typename view_type::const_reference const_reference;

        /** size type
         */
    typedef typename view_type::size_type size_type;

        /** difference type (used for multi-dimensional offsets and indices)
         */
    typedef typename view_type::difference_type difference_type;

        /** difference and index type for a single dimension
         */
    typedef typename view_type::difference_type_1 difference_type_1;

        /** traverser type
         */
    typedef typename view_type::traverser traverser;

        /** traverser type to const data
         */
    typedef typename view_type::const_traverser  const_traverser;

        // /** sequential (random access) iterator type
         // */
    // typedef typename vigra::detail::MultiIteratorChooser<actual_stride>::template Iterator<N, value_type, reference, pointer>::type
    // iterator;

        // /** sequential (random access) const iterator type
         // */
    // typedef typename vigra::detail::MultiIteratorChooser<actual_stride>::template Iterator<N, value_type, const_reference, const_pointer>::type
    // const_iterator;

        /** sequential (random access) iterator type
         */
    typedef typename view_type::iterator iterator;

        /** sequential (random access) const iterator type
         */
    typedef typename view_type::const_iterator const_iterator;

protected:

    typedef typename difference_type::value_type diff_zero_t;

        /** the allocator used to allocate the memory
         */
    allocator_type m_alloc;

        /** allocate memory for s pixels, write its address into the given
            pointer and initialize the pixels with init.
        */
    void allocate (pointer &ptr, difference_type_1 s, const_reference init);

        /** allocate memory for s pixels, write its address into the given
            pointer and initialize the linearized pixels to the values of init.
        */
    template <class U>
    void allocate (pointer &ptr, difference_type_1 s, U const * init);

        /** allocate memory, write its address into the given
            pointer and initialize it by copying the data from the given MultiArrayView.
        */
    template <class U, class StrideTag>
    void allocate (pointer &ptr, MultiArrayView<N, U, StrideTag> const & init);

        /** deallocate the memory (of length s) starting at the given address.
         */
    void deallocate (pointer &ptr, difference_type_1 s);

    template <class U, class StrideTag>
    void copyOrReshape (const MultiArrayView<N, U, StrideTag> &rhs);
public:
        /** default constructor
         */
    MultiArray ()
    : view_type (difference_type (diff_zero_t(0)),
                 difference_type (diff_zero_t(0)), 0)
    {}

        /** construct with given allocator
         */
    MultiArray (allocator_type const & alloc)
    : view_type(difference_type (diff_zero_t(0)),
                difference_type (diff_zero_t(0)), 0),
      m_alloc(alloc)
    {}

        /** construct with given length
        
            Use only for 1-dimensional arrays (<tt>N==1</tt>).
         */
    explicit MultiArray (difference_type_1 length,
                         allocator_type const & alloc = allocator_type());


        /** construct with given width and height
        
            Use only for 2-dimensional arrays (<tt>N==2</tt>).
         */
    MultiArray (difference_type_1 width, difference_type_1 height,
                         allocator_type const & alloc = allocator_type());

        /** construct with given shape
         */
    explicit MultiArray (const difference_type &shape,
                         allocator_type const & alloc = allocator_type());

        /** construct from shape with an initial value
         */
    MultiArray (const difference_type &shape, const_reference init,
                allocator_type const & alloc = allocator_type());

        /** construct from shape and initialize with a linear sequence in scan order 
            (i.e. first pixel gets value 0, second on gets value 1 and so on).
         */
    MultiArray (const difference_type &shape, MultiArrayInitializationTag init,
                allocator_type const & alloc = allocator_type());

        /** construct from shape and copy values from the given array
         */
    MultiArray (const difference_type &shape, const_pointer init,
                         allocator_type const & alloc = allocator_type());

        /** copy constructor
         */
    MultiArray (const MultiArray &rhs)
    : view_type(rhs.m_shape, rhs.m_stride, 0),
      m_alloc (rhs.m_alloc)
    {
        allocate (this->m_ptr, this->elementCount (), rhs.data ());
    }

        /** constructor from an array expression
         */
    template<class Expression>
    MultiArray (multi_math::MultiMathOperand<Expression> const & rhs,
                allocator_type const & alloc = allocator_type())
    : view_type(difference_type (diff_zero_t(0)),
                difference_type (diff_zero_t(0)), 0),
      m_alloc (alloc)
    {
        multi_math::math_detail::assignOrResize(*this, rhs);
    }

        /** construct by copying from a MultiArrayView
         */
    template <class U, class StrideTag>
    MultiArray (const MultiArrayView<N, U, StrideTag>  &rhs,
                allocator_type const & alloc = allocator_type());

        /** assignment.<br>
            If the size of \a rhs is the same as the left-hand side arrays's old size, only
            the data are copied. Otherwise, new storage is allocated, which invalidates all
            objects (array views, iterators) depending on the lhs array.
         */
    MultiArray & operator= (const MultiArray &rhs)
    {
        if (this != &rhs)
            this->copyOrReshape(rhs);
        return *this;
    }

        /** assignment from arbitrary MultiArrayView.<br>
            If the size of \a rhs is the same as the left-hand side arrays's old size, only
            the data are copied. Otherwise, new storage is allocated, which invalidates all
            objects (array views, iterators) depending on the lhs array.
         */
    template <class U, class StrideTag>
    MultiArray &operator= (const MultiArrayView<N, U, StrideTag> &rhs)
    {
        this->copyOrReshape(rhs);
        return *this;
    }

        /** assignment from scalar.<br>
            Equivalent to MultiArray::init(v).
         */
    MultiArray & operator=(value_type const & v)
    {
        return this->init(v);
    }

        /** Add-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
            If the left array has no data (hasData() is false), this function is 
            equivalent to a normal assignment (i.e. an empty
            array is interpreted as a zero-array of appropriate size).
         */
    template <class U, class StrideTag>
    MultiArray &operator+= (const MultiArrayView<N, U, StrideTag> &rhs)
    {
        if(this->hasData())
            view_type::operator+=(rhs);
        else
            *this = rhs;
        return *this;
    }

        /** Subtract-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
            If the left array has no data (hasData() is false), this function is 
            equivalent to an assignment of the negated rhs (i.e. an empty
            array is interpreted as a zero-array of appropriate size).
         */
    template <class U, class StrideTag>
    MultiArray &operator-= (const MultiArrayView<N, U, StrideTag> &rhs)
    {
        if(!this->hasData())
            this->reshape(rhs.shape());
        view_type::operator-=(rhs);
        return *this;
    }

        /** Multiply-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
            If the left array has no data (hasData() is false), this function is 
            equivalent to reshape(rhs.shape()) with zero initialisation (i.e. an empty
            array is interpreted as a zero-array of appropriate size).
         */
    template <class U, class StrideTag>
    MultiArray &operator*= (const MultiArrayView<N, U, StrideTag> &rhs)
    {
        if(this->hasData())
            view_type::operator*=(rhs);
        else
            this->reshape(rhs.shape());
        return *this;
    }

        /** Divide-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
            If the left array has no data (hasData() is false), this function is 
            equivalent to reshape(rhs.shape()) with zero initialisation (i.e. an empty
            array is interpreted as a zero-array of appropriate size).
         */
    template <class U, class StrideTag>
    MultiArray &operator/= (const MultiArrayView<N, U, StrideTag> &rhs)
    {
        if(this->hasData())
            view_type::operator/=(rhs);
        else
            this->reshape(rhs.shape());
        return *this;
    }

        /** Add-assignment of a scalar.
         */
    MultiArray &operator+= (const T &rhs)
    {
        view_type::operator+=(rhs);
        return *this;
    }

        /** Subtract-assignment of a scalar.
         */
    MultiArray &operator-= (const T &rhs)
    {
        view_type::operator-=(rhs);
        return *this;
    }

        /** Multiply-assignment of a scalar.
         */
    MultiArray &operator*= (const T &rhs)
    {
        view_type::operator*=(rhs);
        return *this;
    }

        /** Divide-assignment of a scalar.
         */
    MultiArray &operator/= (const T &rhs)
    {
        view_type::operator/=(rhs);
        return *this;
    }
        /** Assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::assignOrResize(*this, rhs);
        return *this;
    }

        /** Add-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator+=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::plusAssignOrResize(*this, rhs);
        return *this;
    }

        /** Subtract-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator-=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::minusAssignOrResize(*this, rhs);
        return *this;
    }

        /** Multiply-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator*=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::multiplyAssignOrResize(*this, rhs);
        return *this;
    }

        /** Divide-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator/=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::math_detail::divideAssignOrResize(*this, rhs);
        return *this;
    }

        /** destructor
         */
    ~MultiArray ()
    {
        deallocate (this->m_ptr, this->elementCount ());
    }


        /** init elements with a constant
         */
    template <class U>
    MultiArray & init(const U & init)
    {
        view_type::init(init);
        return *this;
    }

        /** Allocate new memory with the given shape and initialize with zeros.<br>
            <em>Note:</em> this operation invalidates all dependent objects
            (array views and iterators)
         */
    void reshape (const difference_type &shape)
    {
        reshape (shape, value_type());
    }

        /** Allocate new memory with the given shape and initialize it
            with the given value.<br>
            <em>Note:</em> this operation invalidates all dependent objects
            (array views and iterators)
         */
    void reshape (const difference_type &shape, const_reference init);

        /** Swap the contents with another MultiArray. This is fast,
            because no data are copied, but only pointers and shapes swapped.
            <em>Note:</em> this operation invalidates all dependent objects
            (array views and iterators)
         */
    void swap (MultiArray & other);

        // /** sequential iterator pointing to the first array element.
         // */
    // iterator begin ()
    // {
        // return vigra::detail::MultiIteratorChooser<actual_stride>::template constructIterator<iterator>((view_type *)this);
    // }

        // /** sequential iterator pointing beyond the last array element.
         // */
    // iterator end ()
    // {
        // return begin() + this->elementCount();
    // }

        // /** sequential const iterator pointing to the first array element.
         // */
    // const_iterator begin () const
    // {
        // return vigra::detail::MultiIteratorChooser<actual_stride>::template constructIterator<iterator>((view_type const *)this);
    // }

        // /** sequential const iterator pointing beyond the last array element.
         // */
    // const_iterator end () const
    // {
        // return begin() + this->elementCount();
    // }

        /** get the allocator.
         */
    allocator_type const & allocator () const
    {
        return m_alloc;
    }
    
    static difference_type defaultStride(difference_type const & shape)
    {
        return vigra::detail::ResolveMultiband<T>::defaultStride(shape);
    }
};

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (difference_type_1 length,
                                  allocator_type const & alloc)
: view_type(difference_type(length),
            defaultStride(difference_type(length)),
            0),
  m_alloc(alloc)
{
    allocate (this->m_ptr, this->elementCount (), value_type());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (difference_type_1 width, difference_type_1 height,
                                  allocator_type const & alloc)
: view_type(difference_type(width, height),
            defaultStride(difference_type(width, height)),
            0),
  m_alloc(alloc)
{
    allocate (this->m_ptr, this->elementCount (), value_type());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  allocator_type const & alloc)
: view_type(shape,
            defaultStride(shape),
            0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 1;
    }
    allocate (this->m_ptr, this->elementCount (), value_type());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, const_reference init,
                                  allocator_type const & alloc)
: view_type(shape,
            defaultStride(shape),
            0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 1;
    }
    allocate (this->m_ptr, this->elementCount (), init);
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, MultiArrayInitializationTag init,
                                  allocator_type const & alloc)
: view_type(shape,
            defaultStride(shape),
            0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 1;
    }
    allocate (this->m_ptr, this->elementCount (), value_type());
    switch(init)
    {
      case LinearSequence:
        linearSequence(this->begin(), this->end());
        break;
      default:
        vigra_precondition(false,
            "MultiArray(): invalid MultiArrayInitializationTag.");
    }
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, const_pointer init,
                                  allocator_type const & alloc)
: view_type(shape,
            defaultStride(shape),
            0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 1;
    }
    allocate (this->m_ptr, this->elementCount (), init);
}

template <unsigned int N, class T, class A>
template <class U, class StrideTag>
MultiArray <N, T, A>::MultiArray(const MultiArrayView<N, U, StrideTag>  &rhs,
                                 allocator_type const & alloc)
: view_type(rhs.shape(),
            defaultStride(rhs.shape()),
            0),
  m_alloc (alloc)
{
    allocate (this->m_ptr, rhs);
}

template <unsigned int N, class T, class A>
template <class U, class StrideTag>
void
MultiArray <N, T, A>::copyOrReshape(const MultiArrayView<N, U, StrideTag> &rhs)
{
    if (this->shape() == rhs.shape())
        this->copy(rhs);
    else
    {
        MultiArray t(rhs);
        this->swap(t);
    }
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::reshape (const difference_type & new_shape,
                                    const_reference initial)
{
    if (N == 0)
    {
        return;
    }
    else if(new_shape == this->shape())
    {
        this->init(initial);
    }
    else
    {
        difference_type new_stride = defaultStride(new_shape);
        difference_type_1 new_size = prod(new_shape);
        pointer new_ptr = pointer();
        allocate (new_ptr, new_size, initial);
        deallocate (this->m_ptr, this->elementCount ());
        this->m_ptr = new_ptr;
        this->m_shape = new_shape;
        this->m_stride = new_stride;
    }
}


template <unsigned int N, class T, class A>
inline void
MultiArray <N, T, A>::swap (MultiArray & other)
{
    if (this == &other)
        return;
    this->view_type::swap(other);
    std::swap(this->m_alloc,  other.m_alloc);
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::allocate (pointer & ptr, difference_type_1 s,
                                     const_reference init)
{
    if(s == 0)
    {
        ptr = 0;
        return;
    }
    ptr = m_alloc.allocate ((typename A::size_type)s);
    difference_type_1 i = 0;
    try {
        for (; i < s; ++i)
            m_alloc.construct (ptr + i, init);
    }
    catch (...) {
        for (difference_type_1 j = 0; j < i; ++j)
            m_alloc.destroy (ptr + j);
        m_alloc.deallocate (ptr, (typename A::size_type)s);
        throw;
    }
}

template <unsigned int N, class T, class A>
template <class U>
void MultiArray <N, T, A>::allocate (pointer & ptr, difference_type_1 s,
                                     U const * init)
{
    if(s == 0)
    {
        ptr = 0;
        return;
    }
    ptr = m_alloc.allocate ((typename A::size_type)s);
    difference_type_1 i = 0;
    try {
        for (; i < s; ++i, ++init)
            m_alloc.construct (ptr + i, *init);
    }
    catch (...) {
        for (difference_type_1 j = 0; j < i; ++j)
            m_alloc.destroy (ptr + j);
        m_alloc.deallocate (ptr, (typename A::size_type)s);
        throw;
    }
}

template <unsigned int N, class T, class A>
template <class U, class StrideTag>
void MultiArray <N, T, A>::allocate (pointer & ptr, MultiArrayView<N, U, StrideTag> const & init)
{
    difference_type_1 s = init.elementCount();
    if(s == 0)
    {
        ptr = 0;
        return;
    }
    ptr = m_alloc.allocate ((typename A::size_type)s);
    pointer p = ptr;
    try {
        detail::uninitializedCopyMultiArrayData(init.traverser_begin(), init.shape(),
                                                p, m_alloc, MetaInt<actual_dimension-1>());
    }
    catch (...) {
        for (pointer pp = ptr; pp < p; ++pp)
            m_alloc.destroy (pp);
        m_alloc.deallocate (ptr, (typename A::size_type)s);
        throw;
    }
}

template <unsigned int N, class T, class A>
inline void MultiArray <N, T, A>::deallocate (pointer & ptr, difference_type_1 s)
{
    if (ptr == 0)
        return;
    for (difference_type_1 i = 0; i < s; ++i)
        m_alloc.destroy (ptr + i);
    m_alloc.deallocate (ptr, (typename A::size_type)s);
    ptr = 0;
}

/********************************************************/
/*                                                      */
/*              argument object factories               */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class StrideTag>
inline triple<typename MultiArrayView<N,T,StrideTag>::const_traverser,
              typename MultiArrayView<N,T,StrideTag>::difference_type,
              typename AccessorTraits<T>::default_const_accessor >
srcMultiArrayRange( MultiArrayView<N,T,StrideTag> const & array )
{
    return triple<typename MultiArrayView<N,T,StrideTag>::const_traverser,
                  typename MultiArrayView<N,T,StrideTag>::difference_type,
                  typename AccessorTraits<T>::default_const_accessor >
      ( array.traverser_begin(),
        array.shape(),
        typename AccessorTraits<T>::default_const_accessor() );
}

template <unsigned int N, class T, class StrideTag, class Accessor>
inline triple<typename MultiArrayView<N,T,StrideTag>::const_traverser,
              typename MultiArrayView<N,T,StrideTag>::difference_type,
              Accessor >
srcMultiArrayRange( MultiArrayView<N,T,StrideTag> const & array, Accessor a )
{
    return triple<typename MultiArrayView<N,T,StrideTag>::const_traverser,
                  typename MultiArrayView<N,T,StrideTag>::difference_type,
                  Accessor >
      ( array.traverser_begin(),
        array.shape(),
        a);
}

template <unsigned int N, class T, class StrideTag>
inline pair<typename MultiArrayView<N,T,StrideTag>::const_traverser,
            typename AccessorTraits<T>::default_const_accessor >
srcMultiArray( MultiArrayView<N,T,StrideTag> const & array )
{
    return pair<typename MultiArrayView<N,T,StrideTag>::const_traverser,
                typename AccessorTraits<T>::default_const_accessor >
      ( array.traverser_begin(),
        typename AccessorTraits<T>::default_const_accessor() );
}

template <unsigned int N, class T, class StrideTag, class Accessor>
inline pair<typename MultiArrayView<N,T,StrideTag>::const_traverser,
            Accessor >
srcMultiArray( MultiArrayView<N,T,StrideTag> const & array, Accessor a )
{
    return pair<typename MultiArrayView<N,T,StrideTag>::const_traverser,
                Accessor >
      ( array.traverser_begin(), a );
}

template <unsigned int N, class T, class StrideTag>
inline triple<typename MultiArrayView<N,T,StrideTag>::traverser,
              typename MultiArrayView<N,T,StrideTag>::difference_type,
              typename AccessorTraits<T>::default_accessor >
destMultiArrayRange( MultiArrayView<N,T,StrideTag> & array )
{
    return triple<typename MultiArrayView<N,T,StrideTag>::traverser,
                  typename MultiArrayView<N,T,StrideTag>::difference_type,
                  typename AccessorTraits<T>::default_accessor >
      ( array.traverser_begin(),
        array.shape(),
        typename AccessorTraits<T>::default_accessor() );
}

template <unsigned int N, class T, class StrideTag, class Accessor>
inline triple<typename MultiArrayView<N,T,StrideTag>::traverser,
              typename MultiArrayView<N,T,StrideTag>::difference_type,
              Accessor >
destMultiArrayRange( MultiArrayView<N,T,StrideTag> & array, Accessor a )
{
    return triple<typename MultiArrayView<N,T,StrideTag>::traverser,
                  typename MultiArrayView<N,T,StrideTag>::difference_type,
                  Accessor >
      ( array.traverser_begin(),
        array.shape(),
        a );
}

template <unsigned int N, class T, class StrideTag>
inline pair<typename MultiArrayView<N,T,StrideTag>::traverser,
            typename AccessorTraits<T>::default_accessor >
destMultiArray( MultiArrayView<N,T,StrideTag> & array )
{
    return pair<typename MultiArrayView<N,T,StrideTag>::traverser,
                typename AccessorTraits<T>::default_accessor >
        ( array.traverser_begin(),
          typename AccessorTraits<T>::default_accessor() );
}

template <unsigned int N, class T, class StrideTag, class Accessor>
inline pair<typename MultiArrayView<N,T,StrideTag>::traverser,
            Accessor >
destMultiArray( MultiArrayView<N,T,StrideTag> & array, Accessor a )
{
    return pair<typename MultiArrayView<N,T,StrideTag>::traverser,
                Accessor >
        ( array.traverser_begin(), a );
}

/********************************************************************/

template <class PixelType, class Accessor>
inline triple<ConstStridedImageIterator<PixelType>,
              ConstStridedImageIterator<PixelType>, Accessor>
srcImageRange(const MultiArrayView<2, PixelType, StridedArrayTag> & img, Accessor a)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return triple<ConstStridedImageIterator<PixelType>,
                  ConstStridedImageIterator<PixelType>,
                  Accessor>(
                      ul, ul + Size2D(img.shape(0), img.shape(1)), a);
}

template <class PixelType, class Accessor>
inline pair<ConstStridedImageIterator<PixelType>, Accessor>
srcImage(const MultiArrayView<2, PixelType, StridedArrayTag> & img, Accessor a)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return pair<ConstStridedImageIterator<PixelType>, Accessor>
        (ul, a);
}

template <class PixelType, class Accessor>
inline triple<StridedImageIterator<PixelType>,
              StridedImageIterator<PixelType>, Accessor>
destImageRange(MultiArrayView<2, PixelType, StridedArrayTag> & img, Accessor a)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return triple<StridedImageIterator<PixelType>,
                  StridedImageIterator<PixelType>,
                  Accessor>(
                      ul, ul + Size2D(img.shape(0), img.shape(1)), a);
}

template <class PixelType, class Accessor>
inline pair<StridedImageIterator<PixelType>, Accessor>
destImage(MultiArrayView<2, PixelType, StridedArrayTag> & img, Accessor a)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return pair<StridedImageIterator<PixelType>, Accessor>
        (ul, a);
}

template <class PixelType, class Accessor>
inline pair<StridedImageIterator<PixelType>, Accessor>
maskImage(MultiArrayView<2, PixelType, StridedArrayTag> & img, Accessor a)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return pair<StridedImageIterator<PixelType>, Accessor>
        (ul, a);
}

// -------------------------------------------------------------------

template <class PixelType>
inline triple<ConstStridedImageIterator<PixelType>,
              ConstStridedImageIterator<PixelType>,
              typename AccessorTraits<PixelType>::default_const_accessor>
srcImageRange(MultiArrayView<2, PixelType, StridedArrayTag> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_const_accessor Accessor;
    return triple<ConstStridedImageIterator<PixelType>,
                  ConstStridedImageIterator<PixelType>,
                  Accessor>
        (ul, ul + Size2D(img.shape(0), img.shape(1)), Accessor());
}

template <class PixelType>
inline triple<ConstImageIterator<PixelType>,
              ConstImageIterator<PixelType>,
              typename AccessorTraits<PixelType>::default_const_accessor>
srcImageRange(MultiArrayView<2, PixelType, UnstridedArrayTag> const & img)
{
    ConstImageIterator<PixelType>
        ul(img.data(), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_const_accessor Accessor;
    return triple<ConstImageIterator<PixelType>,
                  ConstImageIterator<PixelType>,
                  Accessor>
        (ul, ul + Size2D(img.shape(0), img.shape(1)), Accessor());
}

template <class PixelType>
inline pair< ConstStridedImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_const_accessor>
srcImage(MultiArrayView<2, PixelType, StridedArrayTag> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_const_accessor Accessor;
    return pair<ConstStridedImageIterator<PixelType>,
                Accessor>
        (ul, Accessor());
}

template <class PixelType>
inline pair< ConstImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_const_accessor>
srcImage(MultiArrayView<2, PixelType, UnstridedArrayTag> const & img)
{
    ConstImageIterator<PixelType>
        ul(img.data(), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_const_accessor Accessor;
    return pair<ConstImageIterator<PixelType>,
                Accessor>
        (ul, Accessor());
}

template <class PixelType>
inline triple< StridedImageIterator<PixelType>,
               StridedImageIterator<PixelType>,
               typename AccessorTraits<PixelType>::default_accessor>
destImageRange(MultiArrayView<2, PixelType, StridedArrayTag> & img)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return triple<StridedImageIterator<PixelType>,
                  StridedImageIterator<PixelType>,
                  Accessor>
        (ul, ul + Size2D(img.shape(0), img.shape(1)), Accessor());
}

template <class PixelType>
inline triple< ImageIterator<PixelType>,
               ImageIterator<PixelType>,
               typename AccessorTraits<PixelType>::default_accessor>
destImageRange(MultiArrayView<2, PixelType, UnstridedArrayTag> & img)
{
    ImageIterator<PixelType>
        ul(img.data(), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return triple<ImageIterator<PixelType>,
                  ImageIterator<PixelType>,
                  Accessor>
        (ul, ul + Size2D(img.shape(0), img.shape(1)), Accessor());
}

template <class PixelType>
inline pair< StridedImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_accessor>
destImage(MultiArrayView<2, PixelType, StridedArrayTag> & img)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return pair<StridedImageIterator<PixelType>, Accessor>
        (ul, Accessor());
}

template <class PixelType>
inline pair< ImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_accessor>
destImage(MultiArrayView<2, PixelType, UnstridedArrayTag> & img)
{
    ImageIterator<PixelType> ul(img.data(), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return pair<ImageIterator<PixelType>, Accessor>(ul, Accessor());
}

template <class PixelType>
inline pair< ConstStridedImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_accessor>
maskImage(MultiArrayView<2, PixelType, StridedArrayTag> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return pair<ConstStridedImageIterator<PixelType>, Accessor>
        (ul, Accessor());
}

template <class PixelType>
inline pair< ConstImageIterator<PixelType>,
             typename AccessorTraits<PixelType>::default_accessor>
maskImage(MultiArrayView<2, PixelType, UnstridedArrayTag> const & img)
{
    ConstImageIterator<PixelType>
        ul(img.data(), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return pair<ConstImageIterator<PixelType>, Accessor>
        (ul, Accessor());
}

/********************************************************/
/*                                                      */
/*                  makeBasicImageView                  */
/*                                                      */
/********************************************************/

/** \addtogroup MultiArrayToImage Create BasicImageView from MultiArrayViews
  
    Some convenience functions for wrapping a \ref vigra::MultiArrayView's
    data in a \ref vigra::BasicImageView. 
*/
//@{
/** Create a \ref vigra::BasicImageView from an unstrided 2-dimensional
    \ref vigra::MultiArrayView.

    The \ref vigra::BasicImageView will have the same <tt>value_type </tt>
    as the original \ref vigra::MultiArrayView.
*/
template <class T, class Stride>
BasicImageView <T>
makeBasicImageView (MultiArrayView <2, T, Stride> const &array)
{
    vigra_precondition(array.isUnstrided(),
       "makeBasicImageView(array): array must be unstrided (i.e. array.isUnstrided() == true).");
    return BasicImageView <T> (array.data (), array.shape (0),
                               array.shape (1), array.stride(1));
}

/** Create a \ref vigra::BasicImageView from a 3-dimensional
    \ref vigra::MultiArray.

    This wrapper flattens the two innermost dimensions of the array
    into single rows of the resulting image.
    The \ref vigra::BasicImageView will have the same <tt>value_type </tt>
    as the original \ref vigra::MultiArray.
*/
template <class T>
BasicImageView <T>
makeBasicImageView (MultiArray <3, T> const &array)
{
    vigra_precondition(array.stride(1) == array.shape(0),
               "makeBasicImageView(): cannot join strided dimensions");
    return BasicImageView <T> (array.data (),
                               array.shape (0)*array.shape (1), array.shape (2), array.stride(2));
}

/** Create a \ref vigra::BasicImageView from a 3-dimensional
    \ref vigra::MultiArray.

    This wrapper only works if <tt>T</tt> is a scalar type and the
    array's innermost dimension has size 3. It then re-interprets
    the data array as a 2-dimensional array with value_type
    <tt>RGBValue<T></tt>.
*/
template <class T, class Stride>
BasicImageView <RGBValue<T> >
makeRGBImageView (MultiArrayView<3, T, Stride> const &array)
{
    vigra_precondition(array.shape (0) == 3, 
       "makeRGBImageView(): array.shape(0) must be 3.");
    vigra_precondition(array.isUnstrided(),
       "makeRGBImageView(array): array must be unstrided (i.e. array.isUnstrided() == true).");
    return BasicImageView <RGBValue<T> > (
        reinterpret_cast <RGBValue <T> *> (array.data ()),
        array.shape (1), array.shape (2));
}

//@}

} // namespace vigra

#undef VIGRA_ASSERT_INSIDE

#endif // VIGRA_MULTI_ARRAY_HXX
