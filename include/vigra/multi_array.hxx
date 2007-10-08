/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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
TinyVector <ptrdiff_t, N> defaultStride(const TinyVector <ptrdiff_t, N> &shape)
{
    TinyVector <ptrdiff_t, N> ret;
    ret [0] = 1;
    for (unsigned int i = 1; i < N; ++i)
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

template <unsigned int K>
struct ScanOrderToOffset
{
    template <int N>
    static ptrdiff_t 
    exec(int d, const TinyVector <ptrdiff_t, N> &shape, const TinyVector <ptrdiff_t, N> & stride)
    {
        return stride[N-K] * (d % shape[N-K]) + 
               ScanOrderToOffset<K-1>::exec(d / shape[N-K], shape, stride);
    }
};

template <>
struct ScanOrderToOffset<1u>
{
    template <int N>
    static ptrdiff_t 
    exec(int d, const TinyVector <ptrdiff_t, N> &shape, const TinyVector <ptrdiff_t, N> & stride)
    {
        return stride[N-1] * d;
    }
};

template <>
struct ScanOrderToOffset<0u>
{
    template <int N>
    static ptrdiff_t 
    exec(int d, const TinyVector <ptrdiff_t, N> &shape, const TinyVector <ptrdiff_t, N> & stride)
    {
        return 0;
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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
void
initMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<0>)
{    
    DestIterator dend = d + shape[0];
    for(; d != dend; ++d)
    {
        *d = init;
    }
}

template <class DestIterator, class Shape, class T, int N>
void
initMultiArrayData(DestIterator d, Shape const & shape, T const & init, MetaInt<N>)
{    
    DestIterator dend = d + shape[N];
    for(; d != dend; ++d)
    {
        initMultiArrayData(d.begin(), shape, init, MetaInt<N-1>());
    }
}

template <class SrcIterator, class Shape, class DestIterator>
void
copyMultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>)
{    
    SrcIterator send = s + shape[0];
    for(; s != send; ++s, ++d)
    {
        *d = *s;
    }
}

template <class SrcIterator, class Shape, class DestIterator, int N>
void
copyMultiArrayData(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<N>)
{    
    SrcIterator send = s + shape[N];
    for(; s != send; ++s, ++d)
    {
        copyMultiArrayData(s.begin(), shape, d.begin(), MetaInt<N-1>());
    }
}

template <class SrcIterator, class Shape, class T, class ALLOC>
void
uninitializedCopyMultiArrayData(SrcIterator s, Shape const & shape, T * & d, ALLOC & a, MetaInt<0>)
{    
    SrcIterator send = s + shape[0];
    for(; s != send; ++s, ++d)
    {
        a.construct(d, static_cast<T const &>(*s));
    }
}

template <class SrcIterator, class Shape, class T, class ALLOC, int N>
void
uninitializedCopyMultiArrayData(SrcIterator s, Shape const & shape, T * & d, ALLOC & a, MetaInt<N>)
{    
    SrcIterator send = s + shape[N];
    for(; s != send; ++s)
    {
        uninitializedCopyMultiArrayData(s.begin(), shape, d, a, MetaInt<N-1>());
    }
}

template <class SrcIterator, class Shape, class T>
void
squaredNormOfMultiArray(SrcIterator s, Shape const & shape, T & result, MetaInt<0>)
{    
    SrcIterator send = s + shape[0];
    for(; s != send; ++s)
    {
        result += *s * *s;
    }
}

template <class SrcIterator, class Shape, class T, int N>
void
squaredNormOfMultiArray(SrcIterator s, Shape const & shape, T & result, MetaInt<N>)
{    
    SrcIterator send = s + shape[N];
    for(; s != send; ++s)
    {
        squaredNormOfMultiArray(s.begin(), shape, result, MetaInt<N-1>());
    }
}

template <class SrcIterator, class Shape, class DestIterator>
bool
equalityOfMultiArrays(SrcIterator s, Shape const & shape, DestIterator d, MetaInt<0>)
{    
    SrcIterator send = s + shape[0];
    for(; s != send; ++s, ++d)
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
    for(; s != send; ++s, ++d)
    {
        if(!equalityOfMultiArrays(s.begin(), shape, d.begin(), MetaInt<N-1>()))
            return false;
    }
    return true;
}


} // namespace detail

/********************************************************/
/*                                                      */
/*                     MultiArrayView                   */
/*                                                      */
/********************************************************/

// forward declaration
template <unsigned int N, class T, class C = UnstridedArrayTag>
class MultiArrayView;
template <unsigned int N, class T, class A = std::allocator<T> >
class MultiArray;

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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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

        /** difference type (used for offsetting)
         */
    typedef TinyVector <ptrdiff_t, actual_dimension> difference_type;

        /** size type
         */
    typedef difference_type size_type;

        /** traverser (MultiIterator) type
         */
    typedef typename detail::MultiIteratorChooser <
        C>::template Traverser <actual_dimension, T, T &, T *>::type traverser;

        /** const traverser (MultiIterator) type
         */
    typedef typename detail::MultiIteratorChooser <
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

public:

        /** default constructor: create an empty image of size 0.
         */
    MultiArrayView ()
        : m_shape (diff_zero_t(0)), m_stride (diff_zero_t(0)), m_ptr (0)
    {}

        /** construct from shape and pointer
         */
    MultiArrayView (const difference_type &shape, pointer ptr);

        /** construct from shape, strides (offset of a sample to the next)
            for every dimension) and pointer
         */
    MultiArrayView (const difference_type &shape,
                    const difference_type &stride,
                    pointer ptr);
                    
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
                    
        /** reset the view to point to the same data as the rhs.
         */
    void reset(MultiArrayView const & rhs)
    {
        m_shape  = rhs.m_shape;
        m_stride = rhs.m_stride;
        m_ptr    = rhs.m_ptr;
    }

        /** array access.
         */
    reference operator[] (const difference_type &d)
    {
        return m_ptr [dot (d, m_stride)];
    }

        /** array access.
         */
    const_reference operator[] (const difference_type &d) const
    {
        return m_ptr [dot (d, m_stride)];
    }
    
        /** array access in scan-order sense.
            Mostly useful to support standard indexing for 1-dimensional multi-arrays,
            but works for any N.
         */
    reference operator[](int d)
    {
        return m_ptr [detail::ScanOrderToOffset<N>::exec(d, m_shape, m_stride)];
    }
    
        /** array access in scan-order sense.
            Mostly useful to support standard indexing for 1-dimensional multi-arrays,
            but works for any N.
         */
    const_reference operator[](int d) const
    {
        return m_ptr [detail::ScanOrderToOffset<N>::exec(d, m_shape, m_stride)];
    }

        /** 1D array access. Use only if N == 1.
         */
    reference operator() (int x)
    {
        return m_ptr [m_stride[0]*x];
    }

        /** 2D array access. Use only if N == 2.
         */
    reference operator() (int x, int y)
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y];
    }

        /** 3D array access. Use only if N == 3.
         */
    reference operator() (int x, int y, int z)
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    }

        /** 4D array access. Use only if N == 4.
         */
    reference operator() (int x, int y, int z, int u)
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    }

        /** 5D array access. Use only if N == 5.
         */
    reference operator() (int x, int y, int z, int u, int v)
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    }

        /** 1D const array access. Use only if N == 1.
         */
    const_reference operator() (int x) const
    {
        return m_ptr [m_stride[0]*x];
    }

        /** 2D const array access. Use only if N == 2.
         */
    const_reference operator() (int x, int y) const
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y];
    }

        /** 3D const array access. Use only if N == 3.
         */
    const_reference operator() (int x, int y, int z) const
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    }

        /** 4D const array access. Use only if N == 4.
         */
    const_reference operator() (int x, int y, int z, int u) const
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    }

        /** 5D const array access. Use only if N == 5.
         */
    const_reference operator() (int x, int y, int z, int u, int v) const
    {
        return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    }

        /** Init with a constant.
         */
    template <class U>
    void init(const U & init);

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

        /** bind the M outmost dimensions to certain indices.
            this reduces the dimensionality of the image to
            max { 1, N-M }
        */
    template <unsigned int M>
    MultiArrayView <N-M, T, C> bindOuter (const TinyVector <ptrdiff_t, M> &d) const;

        /** bind the M innermost dimensions to certain indices.
            this reduces the dimensionality of the image to
            max { 1, N-M }
        */
    template <unsigned int M>
    MultiArrayView <N-M, T, StridedArrayTag>
    bindInner (const TinyVector <ptrdiff_t, M> &d) const;

        /** bind dimension M to index d.
            this reduces the dimensionality of the image to
            max { 1, N-1 }
         */
    template <unsigned int M>
    MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type >
    bind (int d) const;

        /** bind the outmost dimension to a certain index.
            this reduces the dimensionality of the image to
            max { 1, N-1 }
        */
    MultiArrayView <N-1, T, C> bindOuter (int d) const;

        /** bind the innermost dimension to a certain index.
            this reduces the dimensionality of the image to
            max { 1, N-1 }
        */
    MultiArrayView <N-1, T, StridedArrayTag> bindInner (int d) const;

        /** bind dimension m to index d.
            this reduces the dimensionality of the image to
            max { 1, N-1 }
         */
    MultiArrayView <N-1, T, StridedArrayTag>
    bindAt (int m, int d) const;

        /** create a rectangular subarray that spans between the
            points p and q, where p is in the subarray, q not.
        */
    MultiArrayView subarray (const difference_type &p,
                             const difference_type &q) const
    {
        const int offset = dot (m_stride, p);
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
        return MultiArrayView <N, T, StridedArrayTag>
            (shape, m_stride * s, m_ptr);
    }

        /** permute the dimensions of the array.
            The function exchanges the meaning of the dimensions without copying the data. 
            In case of 2-dimensional array, this is simply array transposition. In higher dimensions,
            there are more posibilities.
            
            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);
            
            MultiArray<2, double, StridedArrayTag> transposed = array.permuteDimensions(Shape(1,0));
            
            for(unsigned int i=0; i<array.shape(0), ++i)
                for(unsigned int j=0; j<array.shape(1); ++j)
                    assert(array(i, j) == transposed(j, i));
            \endcode
        */
    MultiArrayView <N, T, StridedArrayTag>
    permuteDimensions (const difference_type &s) const
    {
        difference_type shape, stride, check(0);
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

        /** transpose a 2-dimensional array. Use only if N==2.
            
            <b>Usage:</b><br>
            \code
            typedef MultiArray<2, double>::difference_type Shape;
            MultiArray<2, double> array(10, 20);
            
            MultiArray<2, double, StridedArrayTag> transposed = array.transpose();
            
            for(unsigned int i=0; i<array.shape(0), ++i)
                for(unsigned int j=0; j<array.shape(1); ++j)
                    assert(array(i, j) == transposed(j, i));
            \endcode
        */
    MultiArrayView <2, T, StridedArrayTag>
    transpose () const
    {
        difference_type shape(m_shape[1], m_shape[0]), 
                        stride(m_stride[1], m_stride[0]);
        return MultiArrayView <2, T, StridedArrayTag>(shape, stride, m_ptr);
    }

        /** number of the elements in the array.
         */
    std::size_t elementCount () const
    {
        std::size_t ret = m_shape[0];
        for(unsigned int i = 1; i < actual_dimension; ++i)
            ret *= m_shape[i];
        return ret;
    }

        /** number of the elements in the array.
            Same as <tt>elementCount()</tt>. Mostly useful to support the std::vector interface.
         */
    std::size_t size () const
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
    int size (int n) const
    {
        return m_shape [n];
    }

        /** return the array's shape at a certain dimension 
            (same as <tt>size(n)</tt>).
         */
    int shape (int n) const
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
    int stride (int n) const
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

        /** return the squared norm of the array (sum of squares of the array elements).
         */
    typename NormTraits<MultiArrayView>::SquaredNormType squaredNorm() const
    {
        typedef typename NormTraits<MultiArrayView>::SquaredNormType SquaredNormType;
        SquaredNormType res = NumericTraits<SquaredNormType>::zero();  
        detail::squaredNormOfMultiArray(traverser_begin(), shape(), res, MetaInt<actual_dimension-1>());
        return res;
    }

        /** return the norm of the array (equals <tt>sqrt(array.squaredNorm())</tt>).
         */
    typename NormTraits<MultiArrayView>::NormType norm() const
    {
        typedef typename NormTraits<MultiArrayView>::SquaredNormType SquaredNormType;
        return sqrt(static_cast<typename SquareRootTraits<SquaredNormType>::SquareRootArgument>(this->squaredNorm()));
    }

        /** return the pointer to the image data
         */
    pointer data () const
    {
        return m_ptr;
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
MultiArrayView <N, T, C>::MultiArrayView (const difference_type &shape,
                                          pointer ptr)
    : m_shape (shape), m_stride (detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape)), m_ptr (ptr)
{}

template <unsigned int N, class T, class C>
MultiArrayView <N, T, C>::MultiArrayView
(const difference_type &shape, const difference_type &stride, pointer ptr)
    : m_shape (shape), m_stride (stride), m_ptr (ptr)
{}

template <unsigned int N, class T, class C>
MultiArrayView <N, T, C> &
MultiArrayView <N, T, C>::operator=(MultiArrayView const & rhs)
{
    if(this == &rhs)
        return *this;
    vigra_precondition(this->shape() == rhs.shape() || m_ptr == 0,
        "MultiArrayView::operator=(MultiArrayView const &) size mismatch - use MultiArrayView::reset().");
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
template <class U>
void 
MultiArrayView <N, T, C>::init(const U & init)
{
    detail::initMultiArrayData(traverser_begin(), shape(), init, MetaInt<actual_dimension-1>());
}

template <unsigned int N, class T, class C>
template <class U, class CN>
void 
MultiArrayView <N, T, C>::copyImpl(const MultiArrayView <N, U, CN>& rhs)
{
    vigra_precondition (shape () == rhs.shape (),
        "MultiArrayView::copy(): shape mismatch.");
    // check for overlap of this and rhs
    const_pointer first_element = this->m_ptr,
                  last_element = first_element + dot(this->m_shape - difference_type(1), this->m_stride);
    typename MultiArrayView <N, U, CN>::const_pointer 
           rhs_first_element = rhs.data(),
           rhs_last_element = rhs_first_element + dot(rhs.shape() - difference_type(1), rhs.stride());
    if(last_element < rhs_first_element || rhs_last_element < first_element)
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

template <unsigned int N, class T, class C>
template <unsigned int M>
MultiArrayView <N-M, T, C>
MultiArrayView <N, T, C>::bindOuter (const TinyVector <ptrdiff_t, M> &d) const
{
    TinyVector <ptrdiff_t, M> stride;
    stride.init (m_stride.begin () + N-M, m_stride.end ());
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <ptrdiff_t, NNew> inner_shape, inner_stride;
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
MultiArrayView <N, T, C>::bindInner (const TinyVector <ptrdiff_t, M> &d) const
{
    TinyVector <ptrdiff_t, M> stride;
    stride.init (m_stride.begin (), m_stride.end () - N + M);
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <ptrdiff_t, NNew> outer_shape, outer_stride;
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
MultiArrayView <N, T, C>::bind (int d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <ptrdiff_t, NNew> shape, stride;
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
MultiArrayView <N, T, C>::bindOuter (int d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <ptrdiff_t, NNew> inner_shape, inner_stride;
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
MultiArrayView <N, T, C>::bindInner (int d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <ptrdiff_t, NNew> outer_shape, outer_stride;
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
MultiArrayView <N, T, C>::bindAt (int n, int d) const
{
    vigra_precondition (
        n < static_cast <int> (N),
        "MultiArrayView <N, T, C>::bindAt(): dimension out of range.");
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <ptrdiff_t, NNew> shape, stride;
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
"<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>"

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

        /** difference type (used for offsetting)
         */
    typedef typename view_type::difference_type difference_type;

        /** traverser type
         */
    typedef typename detail::MultiIteratorChooser <
        UnstridedArrayTag>::template Traverser <N, T, T &, T *>::type
    traverser;

        /** traverser type to const data
         */
    typedef typename detail::MultiIteratorChooser <
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
    void allocate (pointer &ptr, std::size_t s, const_reference init);

        /** allocate memory for s pixels, write its address into the given
            pointer and initialize the linearized pixels to the values of init.
        */
    template <class U>
    void allocate (pointer &ptr, std::size_t s, U const * init);

        /** allocate memory, write its address into the given
            pointer and initialize it by copying the data from the given MultiArrayView.
        */
    template <class U, class C>
    void allocate (pointer &ptr, MultiArrayView<N, U, C> const & init);

        /** deallocate the memory (of length s) starting at the given address.
         */
    void deallocate (pointer &ptr, std::size_t s);

    template <class U, class C>
    void copyOrReshape (const MultiArrayView<N, U, C> &rhs);
public:

        /** default constructor
         */
    MultiArray ();

        /** construct with given allocator
         */
    MultiArray (allocator_type const & alloc);

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
    MultiArray (const MultiArray &rhs);

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
    MultiArray &operator= (const MultiArray &rhs)
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

        /** destructor
         */
   ~MultiArray ();


        /** change the shape and allocate new memory.<br>
            <em>Note:</em> this operation invalidates all dependent objects
            (array views and iterators)
         */
    void reshape (const difference_type &shape)
    {
        reshape (shape, NumericTraits <T>::zero ());
    }

        /** change the shape, allocate new memory and initialize it
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
MultiArray <N, T, A>::MultiArray ()
    : MultiArrayView <N, T> (difference_type (diff_zero_t(0)), 
                             difference_type (diff_zero_t(0)), 0)
{}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (allocator_type const & alloc)
    : MultiArrayView <N, T> (difference_type (diff_zero_t(0)), 
                             difference_type (diff_zero_t(0)), 0),
      m_alloc(alloc)
{}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape), 0),
      m_alloc(alloc)
{
    if (N == 0)
    {
        this->m_shape [0] = 1;
        this->m_stride [0] = 0;
    }
    allocate (this->m_ptr, this->elementCount (), NumericTraits<T>::zero ());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  const_reference init,
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape), 0),
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
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  const_pointer init,
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (shape), 0),
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
MultiArray <N, T, A>::MultiArray (const MultiArray &rhs)
    : MultiArrayView <N, T> (rhs.m_shape, rhs.m_stride, 0),
    m_alloc (rhs.m_alloc)
{
    allocate (this->m_ptr, this->elementCount (), rhs.data ());
}

template <unsigned int N, class T, class A>
template <class U, class C>
MultiArray <N, T, A>::MultiArray (const MultiArrayView<N, U, C> &rhs,
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (rhs.shape(), 
                             detail::defaultStride <MultiArrayView<N,T>::actual_dimension> (rhs.shape()), 0),
    m_alloc (alloc)
{
    allocate (this->m_ptr, rhs); 
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::~MultiArray ()
{
    deallocate (this->m_ptr, this->elementCount ());
}

template <unsigned int N, class T, class A>
template <class U, class C>
void
MultiArray <N, T, A>::copyOrReshape(const MultiArrayView<N, U, C> &rhs)
{
    if (this->shape() == rhs.shape())
    {
        this->copy(rhs);
    }
    else
    {
        pointer new_ptr;
        allocate (new_ptr, rhs);
        deallocate (this->m_ptr, this->elementCount ());
        this->m_shape = rhs.shape();
        this->m_stride = rhs.stride();
        this->m_ptr = new_ptr;
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
        std::size_t new_size = new_shape [MultiArrayView<N,T>::actual_dimension-1] * new_stride [MultiArrayView<N,T>::actual_dimension-1];
        T *new_ptr;
        allocate (new_ptr, new_size, initial);
        deallocate (this->m_ptr, this->elementCount ());
        this->m_ptr = new_ptr;
        this->m_shape = new_shape;
        this->m_stride = new_stride;
    }
}


template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::swap (MultiArray <N, T, A> & other)
{
    if (this == &other)
        return;
    std::swap(this->m_shape,  other.m_shape);
    std::swap(this->m_stride, other.m_stride);
    std::swap(this->m_ptr,    other.m_ptr);
    std::swap(this->m_alloc,  other.m_alloc);
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::allocate (pointer & ptr, std::size_t s,
                                     const_reference init)
{
    ptr = m_alloc.allocate (s);
    std::size_t i;
    try {
        for (i = 0; i < s; ++i)
            m_alloc.construct (ptr + i, init);
    }
    catch (...) {
        for (std::size_t j = 0; j < i; ++j)
            m_alloc.destroy (ptr + j);
        m_alloc.deallocate (ptr, s);
        throw;
    }
}

template <unsigned int N, class T, class A>
template <class U>
void MultiArray <N, T, A>::allocate (pointer & ptr, std::size_t s,
                                     U const * init)
{
    ptr = m_alloc.allocate (s);
    std::size_t i;
    try {
        for (i = 0; i < s; ++i, ++init)
            m_alloc.construct (ptr + i, *init);
    }
    catch (...) {
        for (std::size_t j = 0; j < i; ++j)
            m_alloc.destroy (ptr + j);
        m_alloc.deallocate (ptr, s);
        throw;
    }
}

template <unsigned int N, class T, class A>
template <class U, class C>
void MultiArray <N, T, A>::allocate (pointer & ptr, MultiArrayView<N, U, C> const & init)
{
    std::size_t s = init.elementCount();
    ptr = m_alloc.allocate (s);
    pointer p = ptr;
    try {
        detail::uninitializedCopyMultiArrayData(init.traverser_begin(), init.shape(), 
                                                p, m_alloc, MetaInt<actual_dimension-1>());
    }
    catch (...) {
        for (pointer pp = ptr; pp < p; ++pp)
            m_alloc.destroy (pp);
        m_alloc.deallocate (ptr, s);
        throw;
    }
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::deallocate (pointer & ptr, std::size_t s)
{
    if (ptr == 0)
        return;
    for (std::size_t i = 0; i < s; ++i)
        m_alloc.destroy (ptr + i);
    m_alloc.deallocate (ptr, s);
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
    <tt>RGBValue&lt; T &gt;</tt>.
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

#endif // VIGRA_MULTI_ARRAY_HXX
