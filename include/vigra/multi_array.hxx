/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_MULTI_ARRAY_HXX
#define VIGRA_MULTI_ARRAY_HXX

#include <memory>
#include <algorithm>
#include "vigra/accessor.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/multi_iterator.hxx"
#include "vigra/metaprogramming.hxx"

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
TinyVector <int, N> defaultStride (const TinyVector <int, N> &shape)
{
    TinyVector <int, N> ret;
    ret [0] = 1;
    for (unsigned int i = 1; i < N; ++i)
        ret [i] = ret [i-1] * shape [i-1];
    return ret;
}

/********************************************************/
/*                                                      */
/*                     MaybeStrided                     */
/*                                                      */
/********************************************************/

/* metatag implementing a test for marking MultiArrays that were
    indexed at the zero'th dimension as strided, and all others as
    unstrided.

<b>\#include</b>
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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

/** \brief Base class for, and view to, \ref vigra::MultiArray.

This class implements the interface of both MultiArray and MultiArrayView.
By default, MultiArrayViews are tagged as unstrided. I necessary, strided arrays are
constructed automatically by calls to a variant of the bind...() function.

<b>\#include</b>
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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

        /** size type
         */
    typedef TinyVector <int, actual_dimension> size_type;

        /** difference type (used for offsetting)
         */
    typedef TinyVector <int, actual_dimension> difference_type;

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

private:

        /** traverse an array element-by-element by using a difference_type
         *  object as access coordinates.
         */
    bool inc_navigator (difference_type &nav)
    {
        // essentially, this is nothing more than incrementing a number
        // in radix representation up to a certain limit, the shape.
        int carry = 1;
        for (unsigned int d = 0; d < actual_dimension; ++d) {
            if (nav [d] + carry < m_shape [d]) {
                nav [d] += carry;
                return true;
            } else {
                nav [d] = 0;
            }
        }
        return false;
    }

public:

        /** default constructor: create an empty image of size 0.
         */
    MultiArrayView ()
        : m_shape (0), m_stride (0), m_ptr (0)
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


#if 0
        /** shape assignment.
         */
    template <class CN>
    reference operator= (const MultiArrayView <N, T, CN>& rhs)
    {
        assert (shape () == rhs.shape ());
        difference_type nav;
        if (shape () != nav)
            (*this) [nav] = rhs [nav];
        while (inc_navigator (nav))
            (*this) [nav] = rhs [nav];
    }
#endif /* #if 0 */

    
#ifndef NO_OUT_OF_LINE_MEMBER_TEMPLATES
        /** bind the M outmost dimensions to certain indices.
            this reduces the dimensionality of the image to
            min { 1, N-M }
        */
    template <unsigned int M>
    MultiArrayView <N-M, T, C> bindOuter (const TinyVector <int, M> &d) const;

        /** bind the M innermost dimensions to certain indices.
            this reduces the dimensionality of the image to
            min { 1, N-M }
        */
    template <unsigned int M>
    MultiArrayView <N-M, T, StridedArrayTag>
    bindInner (const TinyVector <int, M> &d) const;
    
        /** bind dimension M to index d.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
         */
    template <unsigned int M>
    MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type >
    bind (int d) const;
#else
    template <unsigned int M>
    MultiArrayView <N-M, T, C> bindOuter (const TinyVector <int, M> &d) const
    {
        return bindOuterImpl(*this, d);
    }
 
    template <unsigned int M>
    MultiArrayView <N-M, T, StridedArrayTag>
    bindInner (const TinyVector <int, M> &d) const
    {
        return bindInnerImpl(*this, d);
    }

    template <unsigned int M>
    MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type >
    bind (int d) const
    {
        return bindImpl<M>(*this, d);
    }
#endif

        /** bind the outmost dimension to a certain index.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
        */
    MultiArrayView <N-1, T, C> bindOuter (int d) const;
    
        /** bind the innermost dimension to a certain index.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
        */
    MultiArrayView <N-1, T, StridedArrayTag> bindInner (int d) const;

        /** bind dimension m to index d.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
         */
    MultiArrayView <N-1, T, StridedArrayTag>
    bindAt (int m, int d) const;

        /** bind dimension m to index d.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
         */
    MultiArrayView <N-1, T, StridedArrayTag>
    bindRow (int d) const;

        /** bind dimension m to index d.
            this reduces the dimensionality of the image to
            min { 1, N-1 }
         */
    MultiArrayView <N-1, T, C>
    bindColumn (int d) const;

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

        /** number of the elements in the array.
         */
    std::size_t elementCount () const
    {
        return m_shape [actual_dimension-1] * m_stride [actual_dimension-1];
    }

        /** return the array's size (same as the shape).
         */
    const size_type & size () const
    {
        return m_shape;
    }

        /** return the array's shape.
         */
    const difference_type & shape () const
    {
        return m_shape;
    }

        /** return the array's shape at a certain dimension.
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
    : m_shape (shape), m_stride (defaultStride (shape)), m_ptr (ptr)
{}

template <unsigned int N, class T, class C>
MultiArrayView <N, T, C>::MultiArrayView
(const difference_type &shape, const difference_type &stride, pointer ptr)
    : m_shape (shape), m_stride (stride), m_ptr (ptr)
{}

#ifndef NO_OUT_OF_LINE_MEMBER_TEMPLATES
template <unsigned int N, class T, class C>
template <unsigned int M>
MultiArrayView <N-M, T, C>
MultiArrayView <N, T, C>::bindOuter (const TinyVector <int, M> &d) const
{
    TinyVector <int, M> stride;
    stride.init (m_stride.begin () + N-M, m_stride.end ());
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <int, NNew> inner_shape, inner_stride;
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
MultiArrayView <N, T, C>::bindInner (const TinyVector <int, M> &d) const
{
    TinyVector <int, M> stride;
    stride.init (m_stride.begin (), m_stride.end () - N + M);
    pointer ptr = m_ptr + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <int, NNew> outer_shape, outer_stride;
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
    TinyVector <int, NNew> shape, stride;
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

#else // NO_OUT_OF_LINE_MEMBER_TEMPLATES
template <unsigned int N, class T, class C, unsigned int M>
MultiArrayView <N-M, T, C>
bindOuterImpl(MultiArrayView <N, T, C> const & self, const TinyVector <int, M> &d)
{
    TinyVector <int, M> stride;
    stride.init (self.stride().begin () + N-M, self.stride().end ());
    typename MultiArrayView <N, T, C>::pointer ptr = self.data() + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <int, NNew> inner_shape, inner_stride;
    if (N-M == 0)
    {
        inner_shape [0] = 1;
        inner_stride [0] = 0;
    }
    else
    {
        inner_shape.init (self.shape().begin (), self.shape().end () - M);
        inner_stride.init (self.stride().begin (), self.stride().end () - M);
    }
    return MultiArrayView <N-M, T, C> (inner_shape, inner_stride, ptr);
}

template <unsigned int N, class T, class C, unsigned int M>
MultiArrayView <N - M, T, StridedArrayTag>
bindInnerImpl(MultiArrayView <N, T, C> const & self, const TinyVector <int, M> &d)
{
    TinyVector <int, M> stride;
    stride.init (self.stride().begin (), self.stride().end () - N + M);
    typename MultiArrayView <N, T, C>::pointer ptr = self.data() + dot (d, stride);
    static const int NNew = (N-M == 0) ? 1 : N-M;
    TinyVector <int, NNew> outer_shape, outer_stride;
    if (N-M == 0)
    {
        outer_shape [0] = 1;
        outer_stride [0] = 0;
    }
    else
    {
        outer_shape.init (self.shape().begin () + M, self.shape().end ());
        outer_stride.init (self.stride().begin () + M, self.stride().end ());
    }
    return MultiArrayView <N-M, T, StridedArrayTag>
        (outer_shape, outer_stride, ptr);
}

template <unsigned int M, unsigned int N, class T, class C>
MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type >
bindImpl(MultiArrayView <N, T, C> const & self, int d)
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <int, NNew> shape, stride;
    // the remaining dimensions are 0..n-1,n+1..N-1
    if (N-1 == 0)
    {
        shape[0] = 1;
        stride[0] = 0;
    }
    else
    {
        std::copy (self.shape().begin (), self.shape().begin () + M, shape.begin ());
        std::copy (self.shape().begin () + M+1, self.shape().end (),
                   shape.begin () + M);
        std::copy (self.stride().begin (), self.stride().begin () + M, stride.begin ());
        std::copy (self.stride().begin () + M+1, self.stride().end (),
                   stride.begin () + M);
    }
    return MultiArrayView <N-1, T, typename detail::MaybeStrided <M>::type>
        (shape, stride, self.data() + d * self.stride(M));
}

#endif // NO_OUT_OF_LINE_MEMBER_TEMPLATES

template <unsigned int N, class T, class C>
MultiArrayView <N - 1, T, C>
MultiArrayView <N, T, C>::bindOuter (int d) const
{
    static const int NNew = (N-1 == 0) ? 1 : N-1;
    TinyVector <int, NNew> inner_shape, inner_stride;
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
    TinyVector <int, NNew> outer_shape, outer_stride;
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
    TinyVector <int, NNew> shape, stride;
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
MultiArrayView <N-1, T, StridedArrayTag>
MultiArrayView <N, T, C>::bindRow (int d) const
{
    return this->bindInner (d);
}

template <unsigned int N, class T, class C>
MultiArrayView <N-1, T, C>
MultiArrayView <N, T, C>::bindColumn (int d) const
{
    return this->bindOuter (d);
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

    <b>\#include</b>
    "<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

    Namespace: vigra
*/
template <unsigned int N, class T, class A /* default already declared above */>
class MultiArray : public MultiArrayView <N, T>
{

public:

        /** the allocator type used to allocate the memory
         */
    typedef A allocator_type;

        /** the view type associated with this array.
         */
    typedef MultiArrayView <N, T> view_type;

        /** the matrix type associated with this array.
         */
    typedef MultiArray <N, T> matrix_type;

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
    void allocate (pointer &ptr, std::size_t s, const_pointer init);

        /** deallocate the memory (of length s) starting at the given address.
         */
    void deallocate (pointer &ptr, std::size_t s);

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
 
        /** assignment.<br>
            <em>Note:</em> this operation invalidates all dependent objects
            (array views and iterators)
         */
    MultiArray &operator= (const MultiArray &rhs);

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

        /** sequential iterator pointing to the first array element.
         */
    iterator begin ()
    {
        return data();
    }

        /** sequential iterator pointing beyond the last array element.
         */
    iterator end ()
    {
        return data() + elementCount();
    }

        /** sequential const iterator pointing to the first array element.
         */
    const_iterator begin () const 
    {
        return data();
    }

        /** sequential const iterator pointing beyond the last array element.
         */
    const_iterator end () const
    {
        return data() + elementCount();
    }
};

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray ()
    : MultiArrayView <N, T> (difference_type (0), difference_type (0), 0)
{}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (allocator_type const & alloc)
    : MultiArrayView <N, T> (difference_type (0), difference_type (0), 0),
      m_alloc(alloc)
{}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape, 
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <actual_dimension> (shape), 0),
      m_alloc(alloc)
{
    if (N == 0)
    {
        m_shape [0] = 1;
        m_stride [0] = 0;
    }
    allocate (m_ptr, elementCount (), NumericTraits<T>::zero ());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  const_reference init, 
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <actual_dimension> (shape), 0),
      m_alloc(alloc)
{
    if (N == 0)
    {
        m_shape [0] = 1;
        m_stride [0] = 0;
    }
    allocate (m_ptr, elementCount (), init);
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const difference_type &shape,
                                  const_pointer init, 
                                  allocator_type const & alloc)
    : MultiArrayView <N, T> (shape, detail::defaultStride <actual_dimension> (shape), 0),
      m_alloc(alloc)
{
    if (N == 0)
    {
        m_shape [0] = 1;
        m_stride [0] = 0;
    }
    allocate (m_ptr, elementCount (), init);
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::MultiArray (const MultiArray &rhs)
    : MultiArrayView <N, T> (rhs.m_shape, rhs.m_stride, 0),
    m_alloc (rhs.m_alloc)
{
    allocate (m_ptr, elementCount (), rhs.data ());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A>::~MultiArray ()
{
    deallocate (m_ptr, elementCount ());
}

template <unsigned int N, class T, class A>
MultiArray <N, T, A> &
MultiArray <N, T, A>::operator= (const MultiArray &rhs)
{
    if (this == &rhs)
        return *this;
    pointer new_ptr;
    allocate (new_ptr, rhs.elementCount (), rhs.data ());
    deallocate (m_ptr, elementCount ());
    m_alloc = rhs.m_alloc;
    m_shape = rhs.m_shape;
    m_stride = rhs.m_stride;
    m_ptr = new_ptr;
    return *this;
}

template <unsigned int N, class T, class A>
void MultiArray <N, T, A>::reshape (const difference_type & new_shape,
                                    const_reference init)
{
    if (N== 0)
        return;

    difference_type new_stride = detail::defaultStride <actual_dimension> (new_shape);
    std::size_t new_size = new_shape [actual_dimension-1] * new_stride [actual_dimension-1];
    T *new_ptr;
    allocate (new_ptr, new_size, init);
    deallocate (m_ptr, elementCount ());
    m_ptr = new_ptr;
    m_shape = new_shape;
    m_stride = new_stride;
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
void MultiArray <N, T, A>::allocate (pointer & ptr, std::size_t s,
                                     const_pointer init)
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
void MultiArray <N, T, A>::deallocate (pointer & ptr, std::size_t s)
{
    if (ptr == 0)
        return;
    for (std::size_t i = 0; i < s; ++i)
        m_alloc.destroy (ptr + i);
    m_alloc.deallocate (ptr, s);
    ptr = 0;
}

template <class T>
BasicImageView <T>
makeBasicImageView (MultiArrayView <2, T, UnstridedArrayTag> const &array)
{
    return BasicImageView <T> (array.data (), array.shape (0),
                               array.shape (1));
}

template <class T>
BasicImageView <T>
makeBasicImageView (MultiArray <3, T> const &array)
{
    vigra_precondition (
        array.shape (0) == 1, "makeBasicImageView(): array.shape(0) must be 1.");
    return BasicImageView <T> (array.data (), 
                               array.shape (1), array.shape (2));
}

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

} // namespace vigra

#endif // VIGRA_MULTI_ARRAY_HXX
