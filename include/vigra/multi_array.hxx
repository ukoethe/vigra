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
#include "basicimageview.hxx"
#include "imageiterator.hxx"
#include "numerictraits.hxx"
#include "multi_iterator.hxx"
#include "metaprogramming.hxx"
#include "mathutil.hxx"

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
/*                    defaultStride                     */
/*                                                      */
/********************************************************/

/* generates the stride for a gapless shape.

    Namespace: vigra::detail
*/
template <unsigned int N>
inline TinyVector <MultiArrayIndex, N>
defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
{
    TinyVector <MultiArrayIndex, N> ret;
    ret [0] = 1;
    for (int i = 1; i < (int)N; ++i)
        ret [i] = ret [i-1] * shape [i-1];
    return ret;
}

/********************************************************/
/*                                                      */
/*                 ScanOrderToOffset                    */
/*                                                      */
/********************************************************/

/* transforms an index in scan order sense to a pointer offset in a possibly
   strided, multi-dimensional array.

    Namespace: vigra::detail
*/

template <int K>
struct ScanOrderToOffset
{
    template <int N>
    static MultiArrayIndex
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> &shape,
         const TinyVector <MultiArrayIndex, N> & stride)
    {
        return stride[N-K] * (d % shape[N-K]) +
               ScanOrderToOffset<K-1>::exec(d / shape[N-K], shape, stride);
    }
};

template <>
struct ScanOrderToOffset<1>
{
    template <int N>
    static MultiArrayIndex
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> & /*shape*/,
         const TinyVector <MultiArrayIndex, N> & stride)
    {
        return stride[N-1] * d;
    }
};

template <int K>
struct ScanOrderToCoordinate
{
    template <int N>
    static void
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> &shape,
         TinyVector <MultiArrayIndex, N> & result)
    {
        result[N-K] = (d % shape[N-K]);
        ScanOrderToCoordinate<K-1>::exec(d / shape[N-K], shape, result);
    }
};

template <>
struct ScanOrderToCoordinate<1>
{
    template <int N>
    static void
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> & /*shape*/,
         TinyVector <MultiArrayIndex, N> & result)
    {
        result[N-1] = d;
    }
};

template <int K>
struct CoordinateToScanOrder
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> &shape,
         const TinyVector <MultiArrayIndex, N> & coordinate)
    {
        return coordinate[N-K] + shape[N-K] * CoordinateToScanOrder<K-1>::exec(shape, coordinate);
    }
};

template <>
struct CoordinateToScanOrder<1>
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & /*shape*/,
         const TinyVector <MultiArrayIndex, N> & coordinate)
    {
        return coordinate[N-1];
    }
};


template <class C>
struct CoordinatesToOffest
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x)
    {
        return stride[0] * x;
    }
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x, MultiArrayIndex y)
    {
        return stride[0] * x + stride[1] * y;
    }
};

template <>
struct CoordinatesToOffest<UnstridedArrayTag>
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & /*stride*/, MultiArrayIndex x)
    {
        return x;
    }
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x, MultiArrayIndex y)
    {
        return x + stride[1] * y;
    }
};

/********************************************************/
/*                                                      */
/*                     MaybeStrided                     */
/*                                                      */
/********************************************************/

/* metatag implementing a test for marking MultiArrays that were
    indexed at the zero'th dimension as strided, and all others as
    unstrided.

<b>\#include</b>
\<vigra/multi_array.hxx\>

Namespace: vigra::detail
*/
template <unsigned int N>
struct MaybeStrided
{
    typedef UnstridedArrayTag type;
};

template <>
struct MaybeStrided <0>
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

<b>\#include</b>
\<vigra/multi_array.hxx\>

Namespace: vigra::detail
*/
template <class O>
struct MultiIteratorChooser
{
    struct Nil {};

    template <unsigned int N, class T, class REFERENCE, class POINTER>
    struct Traverser
    {
        typedef Nil type;
    };
};

/********************************************************/
/*                                                      */
/*       MultiIteratorChooser <StridedArrayTag>         */
/*                                                      */
/********************************************************/

/* specialization of the MultiIteratorChooser for strided arrays.

<b>\#include</b>
\<vigra/multi_array.hxx\>

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
};

/********************************************************/
/*                                                      */
/*      MultiIteratorChooser <UnstridedArrayTag>        */
/*                                                      */
/********************************************************/

/* specialization of the MultiIteratorChooser for unstrided arrays.

<b>\#include</b>
\<vigra/multi_array.hxx\>

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

#define VIGRA_COPY_MULTI_ARRAY_DATA(name, op) \
template <class SrcIterator, class Shape, class DestIterator> \
inline void \
name##MultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>) \
{     \
    SrcIterator send = s + shape[0]; \
    for(; s < send; ++s, ++d) \
    { \
        *d op detail::RequiresExplicitCast<typename DestIterator::value_type>::cast(*s); \
    } \
} \
 \
template <class SrcIterator, class Shape, class DestIterator, int N> \
void \
name##MultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<N>) \
{ \
    SrcIterator send = s + shape[N]; \
    for(; s < send; ++s, ++d) \
    { \
        name##MultiArrayData(s.begin(), shape, d.begin(), MetaInt<N-1>()); \
    } \
} \
\
template <class DestIterator, class Shape, class T> \
inline void \
name##ScalarMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<0>) \
{     \
    DestIterator dend = d + shape[0]; \
    for(; d < dend; ++d) \
    { \
        *d op detail::RequiresExplicitCast<typename DestIterator::value_type>::cast(init); \
    } \
} \
 \
template <class DestIterator, class Shape, class T, int N> \
void \
name##ScalarMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<N>) \
{     \
    DestIterator dend = d + shape[N]; \
    for(; d < dend; ++d) \
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

template <unsigned int N, class T, class C = UnstridedArrayTag>
class MultiArrayView;
template <unsigned int N, class T, class A = std::allocator<T> >
class MultiArray;

namespace multi_math {

template <class T>
struct MultiMathOperand;

namespace detail {

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

} // namespace detail

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
: public NormTraits<MultiArrayView<N, T, UnstridedArrayTag> >
{
    typedef NormTraits<MultiArrayView<N, T, UnstridedArrayTag> > BaseType;
    typedef MultiArray<N, T, A>                                  Type;
    typedef typename BaseType::SquaredNormType                   SquaredNormType;
    typedef typename BaseType::NormType                          NormType;
};

/** \brief Base class for, and view to, \ref vigra::MultiArray.

This class implements the interface of both MultiArray and
MultiArrayView.  By default, MultiArrayViews are tagged as
unstrided. If necessary, strided arrays are constructed automatically
by calls to a variant of the bind...() function.

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
       Possible values: UnstridedArrayTag (default), StridedArrayTag
\endcode

<b>\#include</b>
\<vigra/multi_array.hxx\>

Namespace: vigra
*/
template <unsigned int N, class T, class C>
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
        C>::template Traverser <actual_dimension, T, T &, T *>::type traverser;

        /** const traverser (MultiIterator) type
         */
    typedef typename vigra::detail::MultiIteratorChooser <
        C>::template Traverser <actual_dimension, T, T const &, T const *>::type const_traverser;

        /** the view type associated with this array.
         */
    typedef MultiArrayView <N, T, C> view_type;

        /** the matrix type associated with this array.
         */
    typedef MultiArray <N, T> matrix_type;

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

        /** construct from shape and pointer
         */
    MultiArrayView (const difference_type &shape, pointer ptr)
    : m_shape (shape),
      m_stride (detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape)),
      m_ptr (ptr)
    {}

        /** Construct from shape, strides (offset of a sample to the
            next) for every dimension, and pointer.  (Note that
            strides are not given in bytes, but in offsets of the
            respective pointer type.)
         */
    MultiArrayView (const difference_type &shape,
                    const difference_type &stride,
                    pointer ptr)
    : m_shape (shape),
      m_stride (stride),
      m_ptr (ptr)
    {}


        /** Assignment. There are 3 cases:

            <ul>
            <li> When this <tt>MultiArrayView</tt> does not point to valid data
                 (e.g. after default construction), it becomes a copy of \a rhs.
            <li> When the shapes of the two arrays match, the array contents are copied.
            <li> Otherwise, a <tt>PreconditionViolation</tt> exception is thrown.
            </ul>
         */
    MultiArrayView & operator=(MultiArrayView const & rhs);

        /** Assignment of a differently typed MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class U, class C1>
    MultiArrayView & operator=(MultiArrayView<N, U, C1> const & rhs)
    {
        vigra_precondition(this->shape() == rhs.shape(),
            "MultiArrayView::operator=() size mismatch.");
        this->copyImpl(rhs);
        return *this;
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
        multi_math::detail::assign(*this, rhs);
        return *this;
    }

        /** Add-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator+=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::plusAssign(*this, rhs);
        return *this;
    }

        /** Subtract-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator-=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::minusAssign(*this, rhs);
        return *this;
    }

        /** Multiply-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator*=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::multiplyAssign(*this, rhs);
        return *this;
    }

        /** Divide-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArrayView & operator/=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::divideAssign(*this, rhs);
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

        /** equvalent to bindInner(), when M < N.
         */
    template <unsigned int M>
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
        return m_ptr [detail::CoordinatesToOffest<C>::exec(m_stride, x)];
    }

        /** 2D array access. Use only if N == 2.
         */
    reference operator() (difference_type_1 x, difference_type_1 y)
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y));
        return m_ptr [detail::CoordinatesToOffest<C>::exec(m_stride, x, y)];
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
        return m_ptr [detail::CoordinatesToOffest<C>::exec(m_stride, x)];
    }

        /** 2D const array access. Use only if N == 2.
         */
    const_reference operator() (difference_type_1 x, difference_type_1 y) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(x, y));
        return m_ptr [detail::CoordinatesToOffest<C>::exec(m_stride, x, y)];
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
        difference_type p = shape() - difference_type(1);
        for(unsigned int k = dimension+1; k < N; ++k)
            p[k] = 0;
        return (&operator[](p) - m_ptr) == coordinateToScanOrderIndex(p);
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
    template <unsigned int M>
    MultiArrayView <N-M, T, C> bindOuter (const TinyVector <MultiArrayIndex, M> &d) const;

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
    template <unsigned int M>
    MultiArrayView <N-M, T, StridedArrayTag>
    bindInner (const TinyVector <MultiArrayIndex, M> &d) const;

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
    MultiArrayView <N-1, T, typename vigra::detail::MaybeStrided <M>::type >
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
    MultiArrayView <N-1, T, C> bindOuter (difference_type_1 d) const;

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
            layout is equivalent to a fixed-size C array, by specialzing 
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
            array dimension. The new deimension is inserted at index position 'd',
            which must be between 0 and N inclusive.
            
            Possible value types of the original array are: \ref TinyVector, \ref RGBValue, 
            \ref FFTWComplex, <tt>std::complex</tt>, and the built-in number types (in this 
            case, <tt>expandElements</tt> is equivalent to <tt>insertSingletonDimension</tt>). 
            The list of supported types can be extended to any type whose memory
            layout is equivalent to a fixed-size C array, by specialzing 
            <tt>ExpandElementResult</tt>.

            <b>Usage:</b>
            \code
                MultiArray<2, RGBValue<float> > rgb_image(Shape2(w, h));
                
                MultiArrayView<3, float, StridedArrayTag> multiband_image = rgb_image.expandElements(2);
            \endcode
        */
    MultiArrayView <N+1, typename ExpandElementResult<T>::type, StridedArrayTag> 
    expandElements(difference_type_1 d) const;
    
        /** Add a singleton dimension (dimension of legth 1).

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
    MultiArrayView <N+1, T, C>
    insertSingletonDimension (difference_type_1 i) const;

        /** create a rectangular subarray that spans between the
            points p and q, where p is in the subarray, q not.

            <b>Usage:</b>
            \code
            // create a 3D array of size 40x30x20
            typedef MultiArray<3, double>::difference_type Shape;
            MultiArray<3, double> array3(Shape(40, 30, 20));

            // get a subarray set is smaller by one element at all sides
            MultiArrayView <3, double> subarray = array3.subarray(Shape(1,1,1), Shape(39, 29, 19));
            \endcode
        */
    MultiArrayView subarray (const difference_type &p,
                             const difference_type &q) const
    {
        const difference_type_1 offset = dot (m_stride, p);
        return MultiArrayView (q - p, m_stride, m_ptr + offset);
    }

        /** apply an additional striding to the image, thereby reducing
            the shape of the array.
            for example, multiplying the stride of dimension one by three
            turns an appropriately layed out (interleaved) rgb image into
            a single band image.
        */
    MultiArrayView <N, T, StridedArrayTag>
    stridearray (const difference_type &s) const
    {
        difference_type shape = m_shape;
        for (unsigned int i = 0; i < actual_dimension; ++i)
            shape [i] /= s [i];
        return MultiArrayView <N, T, StridedArrayTag>(shape, m_stride * s, m_ptr);
    }

        /** Transpose an array. If N==2, this implements the usual matrix transposition.
            For N > 2, it reverses the order of the indices.

            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);

            MultiArray<2, double, StridedArrayTag> transposed = array.transpose();

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

        /** permute the dimensions of the array.
            The function exchanges the meaning of the dimensions without copying the data.
            In case of a 2-dimensional array, this is simply array transposition. In higher dimensions,
            there are more possibilities.

            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);

            MultiArray<2, double, StridedArrayTag> transposed = array.permuteDimensions(Shape(1,0));

            for(int i=0; i<array.shape(0), ++i)
                for(int j=0; j<array.shape(1); ++j)
                    assert(array(i, j) == transposed(j, i));
            \endcode
        */
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
         */
    template <class U>
    void meanVariance(U * mean, U * variance) const
    {
        typedef typename NumericTraits<U>::RealPromote R;
        triple<R, R, R> res(0.0, 0.0, 0.0);
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
            
            MultiArray<2, double> A(rows, cols);
            ... // fill A
            
            // make the first axis a singleton to sum over the first index
            MultiArray<2, double> rowSums(1, cols);
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
        return iterator(m_ptr, m_shape, m_stride);
    }

        /** returns a const scan-order iterator pointing
            to the first array element.
        */
    const_iterator begin() const
    {
        return const_iterator(m_ptr, m_shape, m_stride);
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

    view_type view ()
    {
        return *this;
    }
};

template <unsigned int N, class T, class C>
MultiArrayView<N, T, C> &
MultiArrayView <N, T, C>::operator=(MultiArrayView const & rhs)
{
    if(this == &rhs)
        return *this;
    vigra_precondition(this->shape() == rhs.shape() || m_ptr == 0,
        "MultiArrayView::operator=(MultiArrayView const &) size mismatch.");
    if(m_ptr == 0)
    {
        m_shape  = rhs.m_shape;
        m_stride = rhs.m_stride;
        m_ptr    = rhs.m_ptr;
    }
    else
        this->copyImpl(rhs);
    return *this;
}

template <unsigned int N, class T, class C>
template <class CN>
bool
MultiArrayView <N, T, C>::arraysOverlap(const MultiArrayView <N, T, CN>& rhs) const
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

template <unsigned int N, class T, class C>
template <class U, class CN>
void
MultiArrayView <N, T, C>::copyImpl(const MultiArrayView <N, U, CN>& rhs)
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
template <unsigned int N, class T, class C> \
template<class U, class C1> \
MultiArrayView<N, T, C> &  \
MultiArrayView <N, T, C>::operator op(MultiArrayView<N, U, C1> const & rhs) \
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

template <unsigned int N, class T, class C>
template <class U, class CN>
void
MultiArrayView <N, T, C>::swapDataImpl(MultiArrayView <N, U, CN> rhs)
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

template <unsigned int N, class T, class C>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, C>::permuteDimensions (const difference_type &s) const
{
    difference_type shape, stride, check((typename difference_type::value_type)0);
    for (unsigned int i = 0; i < actual_dimension; ++i)
    {
        shape[i]  = m_shape[s[i]];
        stride[i] = m_stride[s[i]];
        ++check[s[i]];
    }
    vigra_precondition(check == difference_type(1),
       "MultiArrayView::permuteDimensions(): every dimension must occur exactly once.");
    return MultiArrayView <N, T, StridedArrayTag>(shape, stride, m_ptr);
}

template <unsigned int N, class T, class C>
typename MultiArrayView <N, T, C>::difference_type 
MultiArrayView <N, T, C>::strideOrdering(difference_type stride)
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

template <unsigned int N, class T, class C>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, C>::permuteStridesAscending() const
{
    difference_type ordering(strideOrdering(m_stride)), permutation;
    for(MultiArrayIndex k=0; k<N; ++k)
        permutation[ordering[k]] = k;
    return permuteDimensions(permutation);
}

template <unsigned int N, class T, class C>
MultiArrayView <N, T, StridedArrayTag>
MultiArrayView <N, T, C>::permuteStridesDescending() const
{
    difference_type ordering(strideOrdering(m_stride)), permutation;
    for(MultiArrayIndex k=0; k<N; ++k)
        permutation[ordering[N-1-k]] = k;
    return permuteDimensions(permutation);
}

template <unsigned int N, class T, class C>
template <unsigned int M>
MultiArrayView <N-M, T, C>
MultiArrayView <N, T, C>::bindOuter (const TinyVector <MultiArrayIndex, M> &d) const
{
    TinyVector <MultiArrayIndex, M> stride;
    stride.init (m_stride.begin () + N-M, m_stride.end ());
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <MultiArrayIndex, NNew> inner_shape, inner_stride;
    if (N-M == 0)
    {
        inner_shape [0] = 1;
        inner_stride [0] = 0;
    }
    else
    {
        inner_shape.init (m_shape.begin (), m_shape.end () - M);
        inner_stride.init (m_stride.begin (), m_stride.end () - M);
    }
    return MultiArrayView <N-M, T, C> (inner_shape, inner_stride, ptr);
}

template <unsigned int N, class T, class C>
template <unsigned int M>
MultiArrayView <N - M, T, StridedArrayTag>
MultiArrayView <N, T, C>::bindInner (const TinyVector <MultiArrayIndex, M> &d) const
{
    TinyVector <MultiArrayIndex, M> stride;
    stride.init (m_stride.begin (), m_stride.end () - N + M);
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <MultiArrayIndex, NNew> outer_shape, outer_stride;
    if (N-M == 0)
    {
        outer_shape [0] = 1;
        outer_stride [0] = 0;
    }
    else
    {
        outer_shape.init (m_shape.begin () + M, m_shape.end ());
        outer_stride.init (m_stride.begin () + M, m_stride.end ());
    }
    return MultiArrayView <N-M, T, StridedArrayTag>
        (outer_shape, outer_stride, ptr);
}

template <unsigned int N, class T, class C>
template <unsigned int M>
MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type >
MultiArrayView <N, T, C>::bind (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> shape, stride;
    // the remaining dimensions are 0..n-1,n+1..N-1
    if (N-1 == 0)
    {
        shape[0] = 1;
        stride[0] = 0;
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
    return MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type>
        (shape, stride, m_ptr + d * m_stride[M]);
}

template <unsigned int N, class T, class C>
MultiArrayView <N - 1, T, C>
MultiArrayView <N, T, C>::bindOuter (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> inner_shape, inner_stride;
    if (N-1 == 0)
    {
        inner_shape [0] = 1;
        inner_stride [0] = 0;
    }
    else
    {
        inner_shape.init (m_shape.begin (), m_shape.end () - 1);
        inner_stride.init (m_stride.begin (), m_stride.end () - 1);
    }
    return MultiArrayView <N-1, T, C> (inner_shape, inner_stride,
                                       m_ptr + d * m_stride [N-1]);
}

template <unsigned int N, class T, class C>
MultiArrayView <N - 1, T, StridedArrayTag>
MultiArrayView <N, T, C>::bindInner (difference_type_1 d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> outer_shape, outer_stride;
    if (N-1 == 0)
    {
        outer_shape [0] = 1;
        outer_stride [0] = 0;
    }
    else
    {
        outer_shape.init (m_shape.begin () + 1, m_shape.end ());
        outer_stride.init (m_stride.begin () + 1, m_stride.end ());
    }
    return MultiArrayView <N-1, T, StridedArrayTag>
        (outer_shape, outer_stride, m_ptr + d * m_stride [0]);
}

template <unsigned int N, class T, class C>
MultiArrayView <N - 1, T, StridedArrayTag>
MultiArrayView <N, T, C>::bindAt (difference_type_1 n, difference_type_1 d) const
{
    vigra_precondition (
        n < static_cast <int> (N),
        "MultiArrayView <N, T, C>::bindAt(): dimension out of range.");
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <MultiArrayIndex, NNew> shape, stride;
    // the remaining dimensions are 0..n-1,n+1..N-1
    if (N-1 == 0)
    {
        shape [0] = 1;
        stride [0] = 0;
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


template <unsigned int N, class T, class C>
MultiArrayView <N+1, typename ExpandElementResult<T>::type, StridedArrayTag>
MultiArrayView <N, T, C>::expandElements(difference_type_1 d) const
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

template <unsigned int N, class T, class C>
MultiArrayView <N + 1, T, C>
MultiArrayView <N, T, C>::insertSingletonDimension (difference_type_1 i) const
{
    vigra_precondition (
        0 <= i && i <= static_cast <difference_type_1> (N),
        "MultiArrayView <N, T, C>::insertSingletonDimension(): index out of range.");
    TinyVector <MultiArrayIndex, N+1> shape, stride;
    std::copy (m_shape.begin (), m_shape.begin () + i, shape.begin ());
    std::copy (m_shape.begin () + i, m_shape.end (), shape.begin () + i + 1);
    std::copy (m_stride.begin (), m_stride.begin () + i, stride.begin ());
    std::copy (m_stride.begin () + i, m_stride.end (), stride.begin () + i + 1);
    shape[i] = 1;
    stride[i] = 1;

    return MultiArrayView <N+1, T, C>(shape, stride, m_ptr);
}

template <unsigned int N, class T, class C>
typename NormTraits<MultiArrayView <N, T, C> >::NormType
MultiArrayView <N, T, C>::norm(int type, bool useSquaredNorm) const
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

template <unsigned int N, class T, class C>
inline typename NormTraits<MultiArrayView <N, T, C> >::SquaredNormType
squaredNorm(MultiArrayView <N, T, C> const & a)
{
    return a.squaredNorm();
}

template <unsigned int N, class T, class C>
inline typename NormTraits<MultiArrayView <N, T, C> >::NormType
norm(MultiArrayView <N, T, C> const & a)
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

<b>\#include</b>
\<vigra/multi_array.hxx\>

Namespace: vigra
*/
template <unsigned int N, class T, class A /* default already declared above */>
class MultiArray : public MultiArrayView <N, T>
{

public:
    using MultiArrayView <N, T>::actual_dimension;

        /** the allocator type used to allocate the memory
         */
    typedef A allocator_type;

        /** the view type associated with this array.
         */
    typedef MultiArrayView <N, T> view_type;

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
    typedef typename vigra::detail::MultiIteratorChooser <
        UnstridedArrayTag>::template Traverser <N, T, T &, T *>::type
    traverser;

        /** traverser type to const data
         */
    typedef typename vigra::detail::MultiIteratorChooser <
        UnstridedArrayTag>::template Traverser <N, T, T const &, T const *>::type
    const_traverser;

        /** sequential (random access) iterator type
         */
    typedef T * iterator;

        /** sequential (random access) const iterator type
         */
    typedef T * const_iterator;

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
    template <class U, class C>
    void allocate (pointer &ptr, MultiArrayView<N, U, C> const & init);

        /** deallocate the memory (of length s) starting at the given address.
         */
    void deallocate (pointer &ptr, difference_type_1 s);

    template <class U, class C>
    void copyOrReshape (const MultiArrayView<N, U, C> &rhs);
public:
        /** default constructor
         */
    MultiArray ()
    : MultiArrayView <N, T> (difference_type (diff_zero_t(0)),
                             difference_type (diff_zero_t(0)), 0)
    {}

        /** construct with given allocator
         */
    MultiArray (allocator_type const & alloc)
    : MultiArrayView <N, T> (difference_type (diff_zero_t(0)),
                             difference_type (diff_zero_t(0)), 0),
      m_alloc(alloc)
    {}

        /** construct with given shape
         */
    explicit MultiArray (const difference_type &shape,
                         allocator_type const & alloc = allocator_type());

        /** construct from shape with an initial value
         */
    MultiArray (const difference_type &shape, const_reference init,
                allocator_type const & alloc = allocator_type());

        /** construct from shape and copy values from the given array
         */
    MultiArray (const difference_type &shape, const_pointer init,
                         allocator_type const & alloc = allocator_type());

        /** copy constructor
         */
    MultiArray (const MultiArray &rhs)
    : MultiArrayView <N, T> (rhs.m_shape, rhs.m_stride, 0),
      m_alloc (rhs.m_alloc)
    {
        allocate (this->m_ptr, this->elementCount (), rhs.data ());
    }

        /** constructor from an array expression
         */
    template<class Expression>
    MultiArray (multi_math::MultiMathOperand<Expression> const & rhs,
                allocator_type const & alloc = allocator_type())
    : MultiArrayView <N, T> (difference_type (diff_zero_t(0)),
                             difference_type (diff_zero_t(0)), 0),
      m_alloc (alloc)
    {
        multi_math::detail::assignOrResize(*this, rhs);
    }

        /** construct by copying from a MultiArrayView
         */
    template <class U, class C>
    MultiArray (const MultiArrayView<N, U, C>  &rhs,
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
    template <class U, class C>
    MultiArray &operator= (const MultiArrayView<N, U, C> &rhs)
    {
        this->copyOrReshape(rhs);
        return *this;
    }

        /** Add-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template <class U, class C>
    MultiArray &operator+= (const MultiArrayView<N, U, C> &rhs)
    {
        view_type::operator+=(rhs);
        return *this;
    }

        /** Subtract-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template <class U, class C>
    MultiArray &operator-= (const MultiArrayView<N, U, C> &rhs)
    {
        view_type::operator-=(rhs);
        return *this;
    }

        /** Multiply-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template <class U, class C>
    MultiArray &operator*= (const MultiArrayView<N, U, C> &rhs)
    {
        view_type::operator*=(rhs);
        return *this;
    }

        /** Divide-assignment from arbitrary MultiArrayView. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template <class U, class C>
    MultiArray &operator/= (const MultiArrayView<N, U, C> &rhs)
    {
        view_type::operator/=(rhs);
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
        multi_math::detail::assignOrResize(*this, rhs);
        return *this;
    }

        /** Add-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator+=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::plusAssignOrResize(*this, rhs);
        return *this;
    }

        /** Subtract-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator-=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::minusAssignOrResize(*this, rhs);
        return *this;
    }

        /** Multiply-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator*=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::multiplyAssignOrResize(*this, rhs);
        return *this;
    }

        /** Divide-assignment of an array expression. Fails with
            <tt>PreconditionViolation</tt> exception when the shapes do not match.
         */
    template<class Expression>
    MultiArray & operator/=(multi_math::MultiMathOperand<Expression> const & rhs)
    {
        multi_math::detail::divideAssignOrResize(*this, rhs);
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
        reshape (shape, T());
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

        /** sequential iterator pointing to the first array element.
         */
    iterator begin ()
    {
        return this->data();
    }

        /** sequential iterator pointing beyond the last array element.
         */
    iterator end ()
    {
        return this->data() + this->elementCount();
    }

        /** sequential const iterator pointing to the first array element.
         */
    const_iterator begin () const
    {
        return this->data();
    }

        /** sequential const iterator pointing beyond the last array element.
         */
    const_iterator end () const
    {
        return this->data() + this->elementCount();
    }

        /** get the allocator.
         */
    allocator_type const & allocator () const
    {
        return m_alloc;
    }
};

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  allocator_type const & alloc)
: MultiArrayView <N, T> (shape,
                         detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape),
                         0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 0;
    }
    allocate (this->m_ptr, this->elementCount (), T());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, const_reference init,
                                  allocator_type const & alloc)
: MultiArrayView <N, T> (shape,
                         detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape),
                         0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 0;
    }
    allocate (this->m_ptr, this->elementCount (), init);
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, const_pointer init,
                                  allocator_type const & alloc)
: MultiArrayView <N, T> (shape,
                         detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape),
                         0),
  m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 0;
    }
    allocate (this->m_ptr, this->elementCount (), init);
}

template <unsigned int N, class T, class A>
template <class U, class C>
MultiArray <N, T, A>::MultiArray(const MultiArrayView<N, U, C>  &rhs,
                                 allocator_type const & alloc)
: MultiArrayView <N, T> (rhs.shape(),
                         detail::defaultStride <MultiArrayView<N,T>::actual_dimension>(rhs.shape()),
                         0),
  m_alloc (alloc)
{
    allocate (this->m_ptr, rhs);
}

template <unsigned int N, class T, class A>
template <class U, class C>
void
MultiArray <N, T, A>::copyOrReshape(const MultiArrayView<N, U, C> &rhs)
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
    if (N== 0)
    {
        return;
    }
    else if(new_shape == this->shape())
    {
        this->init(initial);
    }
    else
    {
        difference_type new_stride = detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (new_shape);
        difference_type_1 new_size = new_shape [MultiArrayView<N,T>::actual_dimension-1] * new_stride [MultiArrayView<N,T>::actual_dimension-1];
        T *new_ptr;
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
    std::swap(this->m_shape,  other.m_shape);
    std::swap(this->m_stride, other.m_stride);
    std::swap(this->m_ptr,    other.m_ptr);
    std::swap(this->m_alloc,  other.m_alloc);
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::allocate (pointer & ptr, difference_type_1 s,
                                     const_reference init)
{
    ptr = m_alloc.allocate ((typename A::size_type)s);
    difference_type_1 i;
    try {
        for (i = 0; i < s; ++i)
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
    ptr = m_alloc.allocate ((typename A::size_type)s);
    difference_type_1 i;
    try {
        for (i = 0; i < s; ++i, ++init)
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
template <class U, class C>
void MultiArray <N, T, A>::allocate (pointer & ptr, MultiArrayView<N, U, C> const & init)
{
    difference_type_1 s = init.elementCount();
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

template <unsigned int N, class T, class C>
inline triple<typename MultiArrayView<N,T,C>::const_traverser,
              typename MultiArrayView<N,T,C>::difference_type,
              typename AccessorTraits<T>::default_const_accessor >
srcMultiArrayRange( MultiArrayView<N,T,C> const & array )
{
    return triple<typename MultiArrayView<N,T,C>::const_traverser,
                  typename MultiArrayView<N,T,C>::difference_type,
                  typename AccessorTraits<T>::default_const_accessor >
      ( array.traverser_begin(),
        array.shape(),
        typename AccessorTraits<T>::default_const_accessor() );
}

template <unsigned int N, class T, class C, class Accessor>
inline triple<typename MultiArrayView<N,T,C>::const_traverser,
              typename MultiArrayView<N,T,C>::difference_type,
              Accessor >
srcMultiArrayRange( MultiArrayView<N,T,C> const & array, Accessor a )
{
    return triple<typename MultiArrayView<N,T,C>::const_traverser,
                  typename MultiArrayView<N,T,C>::difference_type,
                  Accessor >
      ( array.traverser_begin(),
        array.shape(),
        a);
}

template <unsigned int N, class T, class C>
inline pair<typename MultiArrayView<N,T,C>::const_traverser,
            typename AccessorTraits<T>::default_const_accessor >
srcMultiArray( MultiArrayView<N,T,C> const & array )
{
    return pair<typename MultiArrayView<N,T,C>::const_traverser,
                typename AccessorTraits<T>::default_const_accessor >
      ( array.traverser_begin(),
        typename AccessorTraits<T>::default_const_accessor() );
}

template <unsigned int N, class T, class C, class Accessor>
inline pair<typename MultiArrayView<N,T,C>::const_traverser,
            Accessor >
srcMultiArray( MultiArrayView<N,T,C> const & array, Accessor a )
{
    return pair<typename MultiArrayView<N,T,C>::const_traverser,
                Accessor >
      ( array.traverser_begin(), a );
}

template <unsigned int N, class T, class C>
inline triple<typename MultiArrayView<N,T,C>::traverser,
              typename MultiArrayView<N,T,C>::difference_type,
              typename AccessorTraits<T>::default_accessor >
destMultiArrayRange( MultiArrayView<N,T,C> & array )
{
    return triple<typename MultiArrayView<N,T,C>::traverser,
                  typename MultiArrayView<N,T,C>::difference_type,
                  typename AccessorTraits<T>::default_accessor >
      ( array.traverser_begin(),
        array.shape(),
        typename AccessorTraits<T>::default_accessor() );
}

template <unsigned int N, class T, class C, class Accessor>
inline triple<typename MultiArrayView<N,T,C>::traverser,
              typename MultiArrayView<N,T,C>::difference_type,
              Accessor >
destMultiArrayRange( MultiArrayView<N,T,C> & array, Accessor a )
{
    return triple<typename MultiArrayView<N,T,C>::traverser,
                  typename MultiArrayView<N,T,C>::difference_type,
                  Accessor >
      ( array.traverser_begin(),
        array.shape(),
        a );
}

template <unsigned int N, class T, class C>
inline pair<typename MultiArrayView<N,T,C>::traverser,
            typename AccessorTraits<T>::default_accessor >
destMultiArray( MultiArrayView<N,T,C> & array )
{
    return pair<typename MultiArrayView<N,T,C>::traverser,
                typename AccessorTraits<T>::default_accessor >
        ( array.traverser_begin(),
          typename AccessorTraits<T>::default_accessor() );
}

template <unsigned int N, class T, class C, class Accessor>
inline pair<typename MultiArrayView<N,T,C>::traverser,
            Accessor >
destMultiArray( MultiArrayView<N,T,C> & array, Accessor a )
{
    return pair<typename MultiArrayView<N,T,C>::traverser,
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

/** \addtogroup MultiArrayToImage Wrap a \ref vigra::MultiArrayView in
                                  a \ref vigra::BasicImageView
*/
//@{
/** Create a \ref vigra::BasicImageView from an unstrided 2-dimensional
    \ref vigra::MultiArrayView.

    The \ref vigra::BasicImageView will have the same <tt>value_type </tt>
    as the original \ref vigra::MultiArrayView.
*/
template <class T>
BasicImageView <T>
makeBasicImageView (MultiArrayView <2, T, UnstridedArrayTag> const &array)
{
    return BasicImageView <T> (array.data (), array.shape (0),
                               array.shape (1));
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
    return BasicImageView <T> (array.data (),
                               array.shape (0)*array.shape (1), array.shape (2));
}

/** Create a \ref vigra::BasicImageView from a 3-dimensional
    \ref vigra::MultiArray.

    This wrapper only works if <tt>T</tt> is a scalar type and the
    array's innermost dimension has size 3. It then re-interprets
    the data array as a 2-dimensional array with value_type
    <tt>RGBValue<T></tt>.
*/
template <class T>
BasicImageView <RGBValue<T> >
makeRGBImageView (MultiArray<3, T> const &array)
{
    vigra_precondition (
        array.shape (0) == 3, "makeRGBImageView(): array.shape(0) must be 3.");
    return BasicImageView <RGBValue<T> > (
        reinterpret_cast <RGBValue <T> *> (array.data ()),
        array.shape (1), array.shape (2));
}

//@}

} // namespace vigra
#undef VIGRA_ASSERT_INSIDE
#endif // VIGRA_MULTI_ARRAY_HXX
