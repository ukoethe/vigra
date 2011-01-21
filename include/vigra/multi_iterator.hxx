/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2008 by Gunnar Kedenburg and Ullrich Koethe       */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.3.0, Sep 10 2004 )                                    */
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


#ifndef VIGRA_MULTI_ITERATOR_HXX
#define VIGRA_MULTI_ITERATOR_HXX

#include <sys/types.h>
#include "tinyvector.hxx"
#include "iteratortags.hxx"

namespace vigra {


template <unsigned int N, class T, 
          class REFERENCE = T &, class POINTER = T *> class MultiIterator;
template <unsigned int N, class T, 
          class REFERENCE = T &, class POINTER = T *> class StridedMultiIterator;

/** \page MultiIteratorPage  Multi-dimensional Array Iterators

    General iterators for arrays of arbitrary dimension.


<p>
<UL style="list-style-image:url(documents/bullet.gif)">
<LI> \ref vigra::MultiArrayShape
     <BR>&nbsp;&nbsp;&nbsp;<em>Difference type for \ref vigra::MultiArrayView or \ref vigra::MultiIterator</em>
<LI> \ref vigra::MultiIterator
     <BR>&nbsp;&nbsp;&nbsp;<em>Iterator for unstrided \ref vigra::MultiArrayView</em>
<LI> \ref vigra::StridedMultiIterator
     <BR>&nbsp;&nbsp;&nbsp;<em>Iterator for strided \ref vigra::MultiArrayView</em>
</UL>
</p>

<p>
    The Multidimensional Iterator concept allows navigation on arrays
    of arbitrary dimension. It provides two modes of iteration: 
    <em>direct traveral</em>, and <em>hierarchical traversal</em>.
    In general, hierarchical traversal will be faster, while only 
    direct traversal allows for true random access in all dimensions.
    Via the <tt>dim<K>()</tt> function, operations applying to a particular
    dimension can be used in the direct traversal mode. In contrast,
    direct traversal functions should not be used in the hierarchical mode
    because the hierarchical functions are only well-defined if the
    iterator points to element 0 in all dimensions below its current dimension.
    The current dimension of a <tt>MultiIterator<N, ...></tt> is <tt>N-1</tt>.
</p>
<h3>Gerneral Requirements for MultiIterator</h3>
<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr><th colspan=2>
    Local Types
    </th><th>
    Meaning
    </th>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::value_type</tt></td><td>the underlying arrays's pixel type</td>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::reference</tt></td>
    <td>the iterator's reference type (return type of <TT>*iter</TT>). Will be
    <tt>value_type &</tt> for a mutable iterator, and convertible to
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::pointer</tt></td>
    <td>the iterator's pointer type (return type of <TT>iter.operator->()</TT>). Will be
    <tt>value_type *</tt> for a mutable iterator, and convertible to
    <tt>value_type const *</tt> for a const iterator.</td>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::iterator_category</tt></td>
    <td>the iterator tag (<tt>vigra::multi_dimensional_traverser_tag</tt>)</td>
</tr>
<tr><th>
    Operation
    </th><th>
    Result
    </th><th>
    Semantics
    </th>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator k;</tt></td><td>default constructor</td>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator k(i);</tt></td><td>copy constructor</td>
</tr>
<tr>
    <td><tt>k = i</tt></td>
    <td><tt>MultiIterator &</tt></td><td>assignment</td>
</tr>
<tr>
    <td><tt>i == j</tt></td><td><tt>bool</tt></td>
    <td>equality (iterators point to the same element)</td>
</tr>
<tr>
    <td><tt>i != j</tt></td><td><tt>bool</tt></td>
    <td>inequality (iterators don't point to the same element)</td>
</tr>
<tr>
    <td><tt>*i</tt></td><td><tt>MultiIterator::reference</tt></td>
    <td>access the current element</td>
</tr>
<tr>
    <td><tt>i->member()</tt></td><td>depends on operation</td>
    <td>call member function of underlying pixel type via <tt>operator-></tt> of iterator</td>
</tr>
</table>
</p>
<h3>Requirements for Direct Traversal</h3>
<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr><th colspan=2>
    Local Types
    </th><th>
    Meaning
    </th>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::multi_difference_type</tt></td>
    <td>the iterator's multi-dimensional difference type (<TT>TinyVector<MultiArrayIndex, N></TT>)</td>
</tr>
<tr><th>
    Operation
    </th><th>
    Result
    </th><th>
    Semantics
    </th>
</tr>
<tr>
    <td><tt>i += diff</tt></td><td><tt>MultiIterator &</tt></td>
    <td>add offset to current position</td>
</tr>
<tr>
    <td><tt>i -= diff</tt></td><td><tt>MultiIterator &</tt></td>
    <td>subtract offset from current position</td>
</tr>
<tr>
    <td><tt>i + diff</tt></td><td><tt>MultiIterator</tt></td>
    <td>create traverser by adding offset</td>
</tr>
<tr>
    <td><tt>i - diff</tt></td><td><tt>MultiIterator</tt></td>
    <td>create traverser by subtracting offset</td>
</tr>
<tr>
    <td><tt>i[diff]</tt></td><td><tt>MultiIterator::reference</tt></td>
    <td>access element at offset <tt>diff</tt></td>
</tr>
<tr>
    <td><tt>i.dim<K>()</tt></td><td><tt>MultiIterator<K+1, T, ...></tt></td>
    <td>Access the traverser with the current dimension set to K. Typically used to call
        navigation functions referring to a particular dimension.<br>
        Example (assuming <tt>i, j</tt> are 3-dimensional):<br>
        \code
        i.dim<0>()++;   // increment dimension 0 
        i.dim<1>()++;   // increment dimension 1 
        i.dim<2>()++;   // increment dimension 2 
        
        j += MultiIterator::multi_difference_type(1,1,1);    // same effect
        \endcode
    </td>
</tr>
<tr><td colspan=3>
       <tt>i, j</tt> are of type <tt>MultiIterator</tt><br>
       <tt>diff</tt> is of type <tt>MultiIterator::multi_difference_type</tt><br>
       <tt>K</tt> is an integer compile-time constant
    </td>
</tr>
</table>
</p>
<p>
Note that it is impossible to support an <tt>operator-</tt> between two iterators which returns
a <tt>MultiIterator::multi_difference_type</tt> because it is impossible to decide to which
dimension a difference applies. Consider for example, a 2-dimensional iterator <tt>i</tt>, and
let <tt>j = i + multi_difference_type(width, 0)</tt>, <tt>k = i + multi_difference_type(0,1)</tt>, 
where <tt>width</tt> is the array's total width. In general, <tt>j</tt> and <tt>k</tt> point to 
the same memory location, so that the two cases cannot easily be distinguished (it is possible,
but iterator performance will suffer significantly, as is experienced with 
\ref vigra::ImageIterator where differencing is allowed).
</p>

<h3>Requirements for Hierarchical Traversal</h3>
<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr><th colspan=2>
    Local Types
    </th><th>
    Meaning
    </th>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::difference_type</tt></td>
    <td>the iterator's difference type (<TT>MultiArrayIndex</TT>)</td>
</tr>
<tr><td colspan=2>
    <tt>MultiIterator::next_type</tt></td><td>type of the next iterator
       (referring to the next lower dimension) in the hierarchy</td>
</tr>
<tr><th>
    Operation
    </th><th>
    Result
    </th><th>
    Semantics
    </th>
</tr>
<tr>
    <td><tt>++i</tt></td><td><tt>MultiIterator &</tt></td>
    <td>pre-increment iterator in its current dimension</td>
</tr>
<tr>
    <td><tt>i++</tt></td><td><tt>MultiIterator</tt></td>
    <td>post-increment iterator in its current dimension</td>
</tr>
<tr>
    <td><tt>--i</tt></td><td><tt>MultiIterator &</tt></td>
    <td>pre-decrement iterator in its current dimension</td>
</tr>
<tr>
    <td><tt>i--</tt></td><td><tt>MultiIterator</tt></td>
    <td>post-decrement iterator in its current dimension</td>
</tr>
<tr>
    <td><tt>i += d</tt></td><td><tt>MultiIterator &</tt></td>
    <td>add <tt>d</tt> in current dimension</td>
</tr>
<tr>
    <td><tt>i -= d</tt></td><td><tt>MultiIterator &</tt></td>
    <td>subtract <tt>d</tt> in from dimension</td>
</tr>
<tr>
    <td><tt>i + d</tt></td><td><tt>MultiIterator</tt></td>
    <td>create new iterator by adding <tt>d</tt> in current dimension</td>
</tr>
<tr>
    <td><tt>i - d</tt></td><td><tt>MultiIterator</tt></td>
    <td>create new iterator by subtracting <tt>d</tt> in current dimension</td>
</tr>
<tr>
    <td><tt>i - j</tt></td><td><tt>difference_type</tt></td>
    <td>difference of <tt>i</tt> and <tt>j</tt> in the current dimension<br>
    <em>Note:</em> The result of this operation is undefined if the iterator
    doesn't point to element 0 in all dimensions below its current dimension.</td>
</tr>
<tr>
    <td><tt>i < j</tt></td><td><tt>bool</tt></td>
    <td><tt>i - j < 0</tt><br>
    <em>Note:</em> The result of this operation is undefined if the iterator
    doesn't point to element 0 in all dimensions below its current dimension.</td>
</tr>
<tr>
    <td><tt>i[d]</tt></td><td><tt>MultiIterator::reference</tt></td>
    <td>access element by adding offset <tt>d</tt> in current dimension</td>
</tr>
<tr>
    <td><tt>i.begin()</tt></td><td><tt>next_type</tt></td>
    <td>create the hierarchical iterator poiting to the first element in the 
    next lower dimension.<br>
    <em>Note:</em> The result of this operation is undefined if the iterator
    doesn't point to element 0 in all dimensions below its current dimension.<br>
    Usage:<br>
    \code
    MultiIterator<3, int> i3 = ..., end3 = ...;
    for(; i3 != end3; ++i3)
    {
        MultiIterator<3, int>::next_type i2 = i3.begin(), end2 = i3.end();
        for(; i2 != end2; ++i2)
        {
            MultiIterator<3, int>::next_type::next_type i1 = i2.begin(), end1 = i2.end();
            for(; i1 != end1; ++i1)
            {
                ... // do something with the current element
            }
        }
    }
    
    \endcode
    </td>
</tr>
<tr>
    <td><tt>i.end()</tt></td><td><tt>next_type</tt></td>
    <td>create the hierarchical iterator poiting to the past-the-end location in the 
    next lower dimension.<br>
    <em>Note:</em> The result of this operation is undefined if the iterator
    doesn't point to element 0 in all dimensions below its current dimension.</td>
</tr>
<tr><td colspan=3>
       <tt>i, j</tt> are of type <tt>MultiIterator</tt><br>
       <tt>d</tt> is of type <tt>MultiIterator::difference_type</tt>
    </td>
</tr>
</table>
</p>

*/

/** \addtogroup MultiIteratorGroup  Multi-dimensional Array Iterators

    \brief General iterators for arrays of arbitrary dimension.
*/
//@{

    /** Index type for a single dimension of a MultiArrayView or
        MultiArray.
    */
typedef std::ptrdiff_t MultiArrayIndex;

    /** Traits class for the difference type of all MultiIterator, MultiArrayView, and
        MultiArray variants.
    */
template <unsigned int N>
class MultiArrayShape
{
  public:
        /** The difference type of all MultiIterator, MultiArrayView, and
            MultiArray variants.
        */
    typedef TinyVector<MultiArrayIndex, N> type;
};

/********************************************************/
/*                                                      */
/*                      MultiIterator                   */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class REFERENCE, class POINTER>
class MultiIterator;

/********************************************************/
/*                                                      */
/*                   MultiIterator<1>                   */
/*                                                      */
/********************************************************/

//
template <class T, class REFERENCE, class POINTER>
class MultiIterator<1, T, REFERENCE, POINTER>
{
  public:
    enum { level = 0 };
    typedef T value_type;
    typedef REFERENCE reference;
    typedef const value_type &const_reference;
    typedef POINTER pointer;
    typedef const value_type *const_pointer;
    typedef typename MultiArrayShape<1>::type multi_difference_type;
    typedef MultiArrayIndex difference_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef std::random_access_iterator_tag iterator_category;

  protected:
    pointer m_ptr;

  public:
    MultiIterator ()
        : m_ptr (0)
    {}

    MultiIterator (pointer ptr,
          const difference_type *,
          const difference_type *)
        : m_ptr (ptr)
    {}

    void operator++ ()
    {
        ++m_ptr;
    }

    void operator-- ()
    {
        --m_ptr;
    }

    MultiIterator operator++ (int)
    {
        MultiIterator ret = *this;
        ++(*this);
        return ret;
    }

    MultiIterator operator-- (int)
    {
        MultiIterator ret = *this;
        --(*this);
        return ret;
    }

    MultiIterator &operator+= (difference_type n)
    {
        m_ptr += n;
        return *this;
    }

    MultiIterator & operator+= (multi_difference_type const & d)
    {
        m_ptr += d[level];
        return *this;
    }

    MultiIterator &operator-= (difference_type n)
    {
        m_ptr -= n;
        return *this;
    }

    MultiIterator & operator-= (multi_difference_type const & d)
    {
        m_ptr -= d[level];
        return *this;
    }

    MultiIterator operator+ (difference_type n) const
    {
        MultiIterator ret = *this;
        ret += n;
        return ret;
    }

    MultiIterator operator+ (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret += d;
        return ret;
    }

    difference_type operator- (MultiIterator const & d) const
    {
        return (m_ptr - d.m_ptr);
    }

    MultiIterator operator- (difference_type n) const
    {
        MultiIterator ret = *this;
        ret -= n;
        return ret;
    }

    MultiIterator operator- (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret -= d;
        return ret;
    }

    reference operator[] (difference_type n) const
    {
        return m_ptr [n];
    }

    reference operator[] (multi_difference_type const & d) const
    {
        return m_ptr [d[level]];
    }

    reference operator* () const
    {
        return *m_ptr;
    }

    pointer get () const
    {
        return m_ptr;
    }

    pointer operator->() const
    {
        return &(operator*());
    }

    bool operator!= (const MultiIterator &rhs) const
    {
        return m_ptr != rhs.m_ptr;
    }

    bool operator== (const MultiIterator &rhs) const
    {
        return m_ptr == rhs.m_ptr;
    }

    bool operator< (const MultiIterator &rhs) const
    {
        return m_ptr < rhs.m_ptr;
    }

    bool operator<= (const MultiIterator &rhs) const
    {
        return m_ptr <= rhs.m_ptr;
    }

    bool operator> (const MultiIterator &rhs) const
    {
        return m_ptr > rhs.m_ptr;
    }

    bool operator>= (const MultiIterator &rhs) const
    {
        return m_ptr >= rhs.m_ptr;
    }

    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d == 0,
            "MultiIterator<1>::iteratorForDimension(d): d == 0 required");
        const difference_type stride = 1;
        return iterator(m_ptr, &stride, 0);
    }

    template <unsigned int K>
    MultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    MultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level];
    }
};

/********************************************************/
/*                                                      */
/*                   MultiIterator<2>                   */
/*                                                      */
/********************************************************/

//
template <class T, class REFERENCE, class POINTER>
class MultiIterator<2, T, REFERENCE, POINTER>
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public MultiIterator<1, T, REFERENCE, POINTER>
#endif
{
  public:

    typedef MultiIterator<1, T, REFERENCE, POINTER> base_type;
    enum { level = 1 };
    typedef T value_type;
    typedef REFERENCE reference;
    typedef const value_type &const_reference;
    typedef POINTER pointer;
    typedef const value_type *const_pointer;
    typedef typename MultiArrayShape<2>::type multi_difference_type;
    typedef MultiArrayIndex difference_type;
    typedef base_type next_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef multi_dimensional_traverser_tag iterator_category;

  protected:
    const difference_type *m_stride;
    const difference_type *m_shape;

  public:
    /* use default copy constructor and assignment operator */

    MultiIterator ()
        : base_type (),
          m_stride (0), m_shape (0)
    {}

    MultiIterator (pointer ptr,
          const difference_type *stride,
          const difference_type *shape)
        : base_type (ptr, stride, shape),
          m_stride (stride), m_shape (shape)
    {}

    void operator++ ()
    {
        this->m_ptr += m_stride [level];
    }

    void operator-- ()
    {
        this->m_ptr -= m_stride [level];
    }

    MultiIterator operator++ (int)
    {
        MultiIterator ret = *this;
        ++(*this);
        return ret;
    }

    MultiIterator operator-- (int)
    {
        MultiIterator ret = *this;
        --(*this);
        return ret;
    }

    MultiIterator & operator+= (difference_type n)
    {
        this->m_ptr += n * m_stride [level];
        return *this;
    }

    MultiIterator & operator+= (multi_difference_type const & d)
    {
        this->m_ptr += total_stride(d.begin());
        return *this;
    }

    MultiIterator  &operator-= (difference_type n)
    {
        this->m_ptr -= n * m_stride [level];
        return *this;
    }

    MultiIterator & operator-= (multi_difference_type const & d)
    {
        this->m_ptr -= total_stride(d.begin());
        return *this;
    }

    MultiIterator operator+ (difference_type n) const
    {
        MultiIterator ret = *this;
        ret += n;
        return ret;
    }

    MultiIterator operator+ (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret += d;
        return ret;
    }

    difference_type operator- (MultiIterator const & d) const
    {
        return (this->m_ptr - d.m_ptr) / this->m_stride[level];
    }

    MultiIterator operator- (difference_type n) const
    {
        MultiIterator ret = *this;
        ret -= n;
        return ret;
    }

    MultiIterator operator- (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret -= d;
        return ret;
    }

    reference operator[] (difference_type n) const
    {
        return this->m_ptr [n*m_stride [level]];
    }

    reference operator[] (multi_difference_type const & d) const
    {
        return this->m_ptr [total_stride(d.begin())];
    }

    next_type begin () const
    {
        return *this;
    }

    next_type end () const
    {
        next_type ret = *this;
        ret += m_shape [level-1];
        return ret;
    }

    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d <= level,
            "MultiIterator<N>::iteratorForDimension(d): d < N required");
        return iterator(this->m_ptr, &m_stride [d], 0);
    }

    template <unsigned int K>
    MultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    MultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }
    MultiIterator<2, T, REFERENCE, POINTER> &
    dim1() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level]*m_stride[level] + base_type::total_stride(d);
    }
};

/********************************************************/
/*                                                      */
/*                   MultiIterator<N>                   */
/*                                                      */
/********************************************************/

/** \brief A multi-dimensional hierarchical iterator to be used with 
           \ref vigra::MultiArrayView if it is not strided.

     See \ref MultiIteratorPage for further documentation.

<b>\#include</b> \<vigra/multi_iterator.hxx\>

Namespace: vigra
*/
template <unsigned int N, class T, class REFERENCE, class POINTER>
class MultiIterator 
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public MultiIterator<N-1, T, REFERENCE, POINTER>
#endif
{
public:

        /** the type of the parent in the inheritance hierarchy.
         */
    typedef MultiIterator<N-1, T, REFERENCE, POINTER> base_type;
        
        /** the iterator's level in the dimension hierarchy
         */
    enum { level = N-1 };

        /** the iterator's value type
         */
    typedef T value_type;

        /** reference type (result of operator[])
         */
    typedef REFERENCE reference;

        /** const reference type (result of operator[] const)
         */
    typedef const value_type &const_reference;

        /** pointer type
         */
    typedef POINTER pointer;

        /** const pointer type
         */
    typedef const value_type *const_pointer;

        /** multi difference type 
            (used for offsetting along all axes simultaneously)
         */
    typedef typename MultiArrayShape<N>::type multi_difference_type;
    
        /** difference type (used for offsetting)
         */
    typedef MultiArrayIndex difference_type;

        /** the MultiIterator for the next lower dimension.
         */
    typedef base_type next_type;

        /** the 1-dimensional iterator for this iterator hierarchy
            (result of iteratorForDimension()).
        */
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;

        /** the iterator tag (image traverser)
        */
    typedef multi_dimensional_traverser_tag iterator_category;
            
    /* use default copy constructor and assignment operator */

        /** default constructor.
        */
    MultiIterator ()
    {}

        /** construct from pointer, strides (offset of a sample to the
            next) for every dimension, and the shape.
        */
    MultiIterator (pointer ptr,
                   const difference_type *stride,
                   const difference_type *shape)
        : base_type (ptr, stride, shape)
    {}


        /** prefix-increment the iterator in it's current dimension
         */
    void operator++ ()
    {
        this->m_ptr += this->m_stride [level];
    }

        /** prefix-decrement the iterator in it's current dimension
         */
    void operator-- ()
    {
        this->m_ptr -= this->m_stride [level];
    }

        /** postfix-increment the iterator in it's current dimension
         */
    MultiIterator operator++ (int)
    {
        MultiIterator ret = *this;
        ++(*this);
        return ret;
    }

        /** postfix-decrement the iterator in it's current dimension
         */
    MultiIterator operator-- (int)
    {
        MultiIterator ret = *this;
        --(*this);
        return ret;
    }

        /** increment the iterator in it's current dimension
            by the given value.
        */
    MultiIterator & operator+= (difference_type n)
    {
        this->m_ptr += n * this->m_stride [level];
        return *this;
    }

        /** increment the iterator in all dimensions
            by the given offset.
        */
    MultiIterator & operator+= (multi_difference_type const & d)
    {
        this->m_ptr += total_stride(d.begin());
        return *this;
    }

        /** decrement the iterator in it's current dimension
            by the given value.
        */
    MultiIterator & operator-= (difference_type n)
    {
        this->m_ptr -= n * this->m_stride [level];
        return *this;
    }

        /** decrement the iterator in all dimensions
            by the given offset.
        */
    MultiIterator & operator-= (multi_difference_type const & d)
    {
        this->m_ptr -= total_stride(d.begin());
        return *this;
    }

        /** addition within current dimension
         */
    MultiIterator operator+ (difference_type n) const
    {
        MultiIterator ret = *this;
        ret += n;
        return ret;
    }

        /** addition along all dimensions
         */
    MultiIterator operator+ (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret += d;
        return ret;
    }

        /** difference of two iterators in the current dimension.
            The result of this operation is undefined if the iterator
            doesn't point to element 0 in all dimensions below its current dimension.
        */
    difference_type operator- (MultiIterator const & d) const
    {
        return (this->m_ptr - d.m_ptr) / this->m_stride[level];
    }

        /** subtraction within current dimension
         */
    MultiIterator operator- (difference_type n) const
    {
        MultiIterator ret = *this;
        ret -= n;
        return ret;
    }

        /** subtraction along all dimensions
         */
    MultiIterator operator- (multi_difference_type const & d) const
    {
        MultiIterator ret = *this;
        ret -= d;
        return ret;
    }

#ifdef DOXYGEN /* documentation only: operators *, ->, ==, !=, <, <=, >, >= are inherited */
        /** derefenrence item
         */
    reference operator* () const;

        /** get address of current item
         */
    pointer get () const;

        /** call method of current item
         */
    pointer operator->() const;

        /** inequality. True if iterators reference different items.
         */
    bool operator!= (const MultiIterator &rhs) const;

        /** equality. True if iterators reference the same items.
         */
    bool operator== (const MultiIterator &rhs) const;

        /** less than.
         */
    bool operator< (const MultiIterator &rhs) const;

        /** less or equal.
         */
    bool operator<= (const MultiIterator &rhs) const;

        /** greater than.
         */
    bool operator> (const MultiIterator &rhs) const;
    
        /** greater or equal.
         */
    bool operator>= (const MultiIterator &rhs) const;
#endif

        /** access the array element at the given offset in 
        the current dimension.
        */
    reference operator[] (difference_type n) const
    {
        return this->m_ptr [n* this->m_stride [level]];
    }

        /** access the array element at the given offset.
        */
    reference operator[] (multi_difference_type const & d) const
    {
        return this->m_ptr [total_stride(d.begin())];
    }

        /** Return the (N-1)-dimensional multi-iterator that points to 
            the first (N-1)-dimensional subarray of the 
            N-dimensional array this iterator is referring to.
            The result is only valid if this iterator refers to location
            0 in <em>all</em> dimensions below its current dimension N,
            otherwise it is undefined. Usage:

            \code

            MultiIterator<2, int> outer = ...;  // this iterator

            MultiIterator<2, int>::next_type inner = outer.begin();
            for(; inner != outer.end(); ++inner)
            {
                // manipulate current 1D subimage
            }
            \endcode
        */
    next_type begin () const
    {
        return *this;
    }

        /** Return the (N-1)-dimensional multi-iterator that points beyond 
            the last (N-1)-dimensional subarray of the 
            N-dimensional array this iterator is referring to.
            The result is only valid if this iterator refers to location
            0 in <em>all</em> dimensions below its current dimension N,
            otherwise it is undefined.
        */
    next_type end () const
    {
        next_type ret = *this;
        ret += this->m_shape [level-1];
        return ret;
    }

        /** Get a 1-dimensional, STL-compatible iterator for the
            given dimension, pointing to the current element of <TT>this</TT>.
            Usage:

            \code

            MultiIterator<3, int> outer = ...;  // this iterator

            MultiIterator<3, int>::iterator i = outer.iteratorForDimension(1);
            MultiIterator<3, int>::iterator end = i + height;
            for(; i != end; ++i)
            {
                // go down the current column starting at the location of 'outer'
            }
            \endcode            
        */
    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d <= level,
            "MultiIterator<N>::iteratorForDimension(d): d < N required");
        return iterator(this->m_ptr, &this->m_stride [d], 0);
    }
        /** Return the multi-iterator that operates on dimension K in order
            to manipulate this dimension directly. Usage:
               
            \code
                
            MultiIterator<3, int> i3 = ...;
                
            i3.template dim<2>()++;  // increment outer dimension
            i3.template dim<0>()++;  // increment inner dimension
            \endcode
            
            For convenience, the same functionality is also available
            as <tt>dim0()</tt>, <tt>dim1()</tt> etc. up to <tt>dim4()</tt>:
            
            \code
                
            MultiIterator<3, int> i3 = ...;
                
            i3.dim2()++;  // increment outer dimension
            i3.dim0()++;  // increment inner dimension
            \endcode            
        */
    template <unsigned int K>
    MultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    MultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }
    MultiIterator<2, T, REFERENCE, POINTER> &
    dim1() { return *this; }
    MultiIterator<3, T, REFERENCE, POINTER> &
    dim2() { return *this; }
    MultiIterator<4, T, REFERENCE, POINTER> &
    dim3() { return *this; }
    MultiIterator<5, T, REFERENCE, POINTER> &
    dim4() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level]*this->m_stride[level] + base_type::total_stride(d);
    }

};

/********************************************************/
/*                                                      */
/*                      StridedMultiIterator            */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class REFERENCE, class POINTER>
class StridedMultiIterator;

/********************************************************/
/*                                                      */
/*                   StridedMultiIterator<1>            */
/*                                                      */
/********************************************************/

//
template <class T, class REFERENCE, class POINTER>
class StridedMultiIterator<1, T, REFERENCE, POINTER>
{
  public:
    enum { level = 0 };
    typedef T value_type;
    typedef REFERENCE reference;
    typedef const value_type &const_reference;
    typedef POINTER pointer;
    typedef const value_type *const_pointer;
    typedef typename MultiArrayShape<1>::type multi_difference_type;
    typedef MultiArrayIndex difference_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef std::random_access_iterator_tag iterator_category;

  protected:
    pointer m_ptr;
    difference_type m_stride;

    /* use default copy constructor and assignment operator */

  public:
    StridedMultiIterator ()
        : m_ptr (0), m_stride (0)
    {}

    StridedMultiIterator (pointer ptr,
          const difference_type *stride,
          const difference_type *)
        : m_ptr (ptr), m_stride (stride [level])
    {}

    void operator++ ()
    {
        m_ptr += m_stride;
    }

    void operator-- ()
    {
        m_ptr -= m_stride;
    }

    StridedMultiIterator operator++ (int)
    {
        StridedMultiIterator ret = *this;
        ++(*this);
        return ret;
    }

    StridedMultiIterator operator-- (int)
    {
        StridedMultiIterator ret = *this;
        --(*this);
        return ret;
    }

    StridedMultiIterator &operator+= (difference_type n)
    {
        m_ptr += n * m_stride;
        return *this;
    }

    StridedMultiIterator & operator+= (multi_difference_type const & d)
    {
        m_ptr += d[level] * m_stride;
        return *this;
    }

    StridedMultiIterator &operator-= (difference_type n)
    {
        m_ptr -= n * m_stride;
        return *this;
    }

    StridedMultiIterator & operator-= (multi_difference_type const & d)
    {
        m_ptr -= d[level] * m_stride;
        return *this;
    }

    StridedMultiIterator operator+ (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret += n;
        return ret;
    }

    StridedMultiIterator operator+ (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret += d;
        return ret;
    }

    difference_type operator- (StridedMultiIterator const & d) const
    {
        return (m_ptr - d.m_ptr) / m_stride;
    }

    StridedMultiIterator operator- (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret -= n;
        return ret;
    }

    StridedMultiIterator operator- (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret -= d;
        return ret;
    }

    reference operator[] (difference_type n) const
    {
        return m_ptr [n*m_stride];
    }

    reference operator[] (multi_difference_type const & d) const
    {
        return m_ptr [d[level]*m_stride];
    }

    reference operator* () const
    {
        return *m_ptr;
    }

    pointer get () const
    {
        return m_ptr;
    }

    pointer operator->() const
    {
        return &(operator*());
    }

    bool operator!= (const StridedMultiIterator &rhs) const
    {
        return m_ptr != rhs.m_ptr;
    }

    bool operator== (const StridedMultiIterator &rhs) const
    {
        return m_ptr == rhs.m_ptr;
    }

    bool operator< (const StridedMultiIterator &rhs) const
    {
        return m_ptr < rhs.m_ptr;
    }

    bool operator<= (const StridedMultiIterator &rhs) const
    {
        return m_ptr <= rhs.m_ptr;
    }

    bool operator> (const StridedMultiIterator &rhs) const
    {
        return m_ptr > rhs.m_ptr;
    }

    bool operator>= (const StridedMultiIterator &rhs) const
    {
        return m_ptr >= rhs.m_ptr;
    }

    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d == 0,
            "StridedMultiIterator<1>::iteratorForDimension(d): d == 0 required");
        const difference_type stride = 1;
        return iterator(m_ptr, &stride, 0);
    }

    template <unsigned int K>
    StridedMultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    StridedMultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level] * m_stride;
    }
};

/********************************************************/
/*                                                      */
/*                   StridedMultiIterator<2>            */
/*                                                      */
/********************************************************/

//
template <class T, class REFERENCE, class POINTER>
class StridedMultiIterator<2, T, REFERENCE, POINTER>
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public StridedMultiIterator<1, T, REFERENCE, POINTER>
#endif
{
  public:

    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> base_type;
    enum { level = 1 };
    typedef T value_type;
    typedef REFERENCE reference;
    typedef const value_type &const_reference;
    typedef POINTER pointer;
    typedef const value_type *const_pointer;
    typedef typename MultiArrayShape<2>::type multi_difference_type;
    typedef MultiArrayIndex difference_type;
    typedef base_type next_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef multi_dimensional_traverser_tag iterator_category;

  protected:
    const difference_type *m_stride;
    const difference_type *m_shape;

  public:
    /* use default copy constructor and assignment operator */

    StridedMultiIterator ()
        : base_type (),
          m_stride (0), m_shape (0)
    {}

    StridedMultiIterator (pointer ptr,
          const difference_type *stride,
          const difference_type *shape)
        : base_type (ptr, stride, shape),
          m_stride (stride), m_shape (shape)
    {}

    void operator++ ()
    {
        this->m_ptr += m_stride [level];
    }

    void operator-- ()
    {
        this->m_ptr -= m_stride [level];
    }

    StridedMultiIterator operator++ (int)
    {
        StridedMultiIterator ret = *this;
        ++(*this);
        return ret;
    }

    StridedMultiIterator operator-- (int)
    {
        StridedMultiIterator ret = *this;
        --(*this);
        return ret;
    }

    StridedMultiIterator & operator+= (difference_type n)
    {
        this->m_ptr += n * m_stride [level];
        return *this;
    }

    StridedMultiIterator & operator+= (multi_difference_type const & d)
    {
        this->m_ptr += total_stride(d.begin());
        return *this;
    }

    StridedMultiIterator  &operator-= (difference_type n)
    {
        this->m_ptr -= n * m_stride [level];
        return *this;
    }

    StridedMultiIterator & operator-= (multi_difference_type const & d)
    {
        this->m_ptr -= total_stride(d.begin());
        return *this;
    }

    StridedMultiIterator operator+ (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret += n;
        return ret;
    }

    StridedMultiIterator operator+ (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret += d;
        return ret;
    }

    difference_type operator- (StridedMultiIterator const & d) const
    {
        return (this->m_ptr - d.m_ptr) / this->m_stride[level];
    }

    StridedMultiIterator operator- (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret -= n;
        return ret;
    }

    StridedMultiIterator operator- (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret -= d;
        return ret;
    }

    reference operator[] (difference_type n) const
    {
        return this->m_ptr [n*m_stride [level]];
    }

    reference operator[] (multi_difference_type const & d) const
    {
        return this->m_ptr [total_stride(d.begin())];
    }

    next_type begin () const
    {
        return *this;
    }

    next_type end () const
    {
        next_type ret = *this;
        ret += m_shape [level-1];
        return ret;
    }

    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d <= level,
            "StridedMultiIterator<N>::iteratorForDimension(d): d < N required");
        return iterator(this->m_ptr, &m_stride [d], 0);
    }

    template <unsigned int K>
    StridedMultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    StridedMultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }
    StridedMultiIterator<2, T, REFERENCE, POINTER> &
    dim1() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level]*m_stride[level] + base_type::total_stride(d);
    }
};

/********************************************************/
/*                                                      */
/*                   StridedMultiIterator<N>            */
/*                                                      */
/********************************************************/

/** \brief A multi-dimensional hierarchical iterator to be used with 
           \ref vigra::MultiArrayView if it is not strided.

     See \ref MultiIteratorPage for further documentation.

<b>\#include</b> \<vigra/multi_iterator.hxx\>

Namespace: vigra
*/
template <unsigned int N, class T, class REFERENCE, class POINTER>
class StridedMultiIterator 
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public StridedMultiIterator<N-1, T, REFERENCE, POINTER>
#endif
{
public:

        /** the type of the parent in the inheritance hierarchy.
         */
    typedef StridedMultiIterator<N-1, T, REFERENCE, POINTER> base_type;
        
        /** the iterator's level in the dimension hierarchy
         */
    enum { level = N-1 };

        /** the iterator's value type
         */
    typedef T value_type;

        /** reference type (result of operator[])
         */
    typedef REFERENCE reference;

        /** const reference type (result of operator[] const)
         */
    typedef const value_type &const_reference;

        /** pointer type
         */
    typedef POINTER pointer;

        /** const pointer type
         */
    typedef const value_type *const_pointer;

        /** multi difference type 
            (used for offsetting along all axes simultaneously)
         */
    typedef typename MultiArrayShape<N>::type multi_difference_type;

        /** difference type (used for offsetting)
         */
    typedef MultiArrayIndex difference_type;
    
        /** the StridedMultiIterator for the next lower dimension.
         */
    typedef base_type next_type;

        /** the 1-dimensional iterator for this iterator hierarchy
            (result of iteratorForDimension()).
        */
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;

        /** the iterator tag (image traverser)
        */
    typedef multi_dimensional_traverser_tag iterator_category;
            
    /* use default copy constructor and assignment operator */

        /** default constructor.
        */
    StridedMultiIterator ()
    {}

        /** construct from pointer, strides (offset of a sample to the
            next) for every dimension, and the shape.
        */
    StridedMultiIterator (pointer ptr,
                   const difference_type *stride,
                   const difference_type *shape)
        : base_type (ptr, stride, shape)
    {}


        /** prefix-increment the iterator in it's current dimension
         */
    void operator++ ()
    {
        this->m_ptr += this->m_stride [level];
    }

        /** prefix-decrement the iterator in it's current dimension
         */
    void operator-- ()
    {
        this->m_ptr -= this->m_stride [level];
    }

        /** postfix-increment the iterator in it's current dimension
         */
    StridedMultiIterator operator++ (int)
    {
        StridedMultiIterator ret = *this;
        ++(*this);
        return ret;
    }

        /** postfix-decrement the iterator in it's current dimension
         */
    StridedMultiIterator operator-- (int)
    {
        StridedMultiIterator ret = *this;
        --(*this);
        return ret;
    }

        /** increment the iterator in it's current dimension
            by the given value.
        */
    StridedMultiIterator & operator+= (difference_type n)
    {
        this->m_ptr += n * this->m_stride [level];
        return *this;
    }

        /** increment the iterator in all dimensions
            by the given offset.
        */
    StridedMultiIterator & operator+= (multi_difference_type const & d)
    {
        this->m_ptr += total_stride(d.begin());
        return *this;
    }

        /** decrement the iterator in it's current dimension
            by the given value.
        */
    StridedMultiIterator & operator-= (difference_type n)
    {
        this->m_ptr -= n * this->m_stride [level];
        return *this;
    }

        /** decrement the iterator in all dimensions
            by the given offset.
        */
    StridedMultiIterator & operator-= (multi_difference_type const & d)
    {
        this->m_ptr -= total_stride(d.begin());
        return *this;
    }

        /** addition within current dimension
         */
    StridedMultiIterator operator+ (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret += n;
        return ret;
    }

        /** addition along all dimensions
         */
    StridedMultiIterator operator+ (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret += d;
        return ret;
    }

        /** difference of two iterators in the current dimension.
            The result of this operation is undefined if the iterator
            doesn't point to element 0 in all dimensions below its current dimension.
        */
    difference_type operator- (StridedMultiIterator const & d) const
    {
        return (this->m_ptr - d.m_ptr) / this->m_stride[level];
    }

        /** subtraction within current dimension
         */
    StridedMultiIterator operator- (difference_type n) const
    {
        StridedMultiIterator ret = *this;
        ret -= n;
        return ret;
    }

        /** subtraction along all dimensions
         */
    StridedMultiIterator operator- (multi_difference_type const & d) const
    {
        StridedMultiIterator ret = *this;
        ret -= d;
        return ret;
    }

#ifdef DOXYGEN /* documentation only: operators *, ->, ==, !=, <, <=, >, >= are inherited */
        /** derefenrence item
         */
    reference operator* () const;

        /** get address of current item
         */
    pointer get () const;

        /** call method of current item
         */
    pointer operator->() const;

        /** inequality. True if iterators reference different items.
         */
    bool operator!= (const StridedMultiIterator &rhs) const;

        /** equality. True if iterators reference the same items.
         */
    bool operator== (const StridedMultiIterator &rhs) const;

        /** less than.
         */
    bool operator< (const StridedMultiIterator &rhs) const;

        /** less or equal.
         */
    bool operator<= (const StridedMultiIterator &rhs) const;

        /** greater than.
         */
    bool operator> (const StridedMultiIterator &rhs) const;
    
        /** greater or equal.
         */
    bool operator>= (const StridedMultiIterator &rhs) const;
#endif

        /** access the array element at the given offset in 
        the current dimension.
        */
    reference operator[] (difference_type n) const
    {
        return this->m_ptr [n* this->m_stride [level]];
    }

        /** access the array element at the given offset.
        */
    reference operator[] (multi_difference_type const & d) const
    {
        return this->m_ptr [total_stride(d.begin())];
    }

        /** Return the (N-1)-dimensional multi-iterator that points to 
            the first (N-1)-dimensional subarray of the 
            N-dimensional array this iterator is referring to.
            The result is only valid if this iterator refers to location
            0 in <em>all</em> dimensions below its current dimension N,
            otherwise it is undefined. Usage:

            \code

            StridedMultiIterator<2, int> outer = ...;  // this iterator

            StridedMultiIterator<2, int>::next_type inner = outer.begin();
            for(; inner != outer.end(); ++inner)
            {
                // manipulate current 1D subimage
            }
            \endcode
        */
    next_type begin () const
    {
        return *this;
    }

        /** Return the (N-1)-dimensional multi-iterator that points beyond 
            the last (N-1)-dimensional subarray of the 
            N-dimensional array this iterator is referring to.
            The result is only valid if this iterator refers to location
            0 in <em>all</em> dimensions below its current dimension N,
            otherwise it is undefined.
        */
    next_type end () const
    {
        next_type ret = *this;
        ret += this->m_shape [level-1];
        return ret;
    }

        /** Get a 1-dimensional, STL-compatible iterator for the
            given dimension, pointing to the current element of <TT>this</TT>.
            Usage:

            \code

            StridedMultiIterator<3, int> outer = ...;  // this iterator

            StridedMultiIterator<3, int>::iterator i = outer.iteratorForDimension(1);
            StridedMultiIterator<3, int>::iterator end = i + height;
            for(; i != end; ++i)
            {
                // go down the current column starting at the location of 'outer'
            }
            \endcode            
        */
    iterator iteratorForDimension(unsigned int d) const
    {
        vigra_precondition(d <= level,
            "StridedMultiIterator<N>::iteratorForDimension(d): d < N required");
        return iterator(this->m_ptr, &this->m_stride [d], 0);
    }
        /** Return the multi-iterator that operates on dimension K in order
            to manipulate this dimension directly. Usage:
               
            \code
                
            StridedMultiIterator<3, int> i3 = ...;
                
            i3.template dim<2>()++;  // increment outer dimension
            i3.template dim<0>()++;  // increment inner dimension
            \endcode
            
            For convenience, the same functionality is also available
            as <tt>dim0()</tt>, <tt>dim1()</tt> etc. up to <tt>dim4()</tt>:
            
            \code
                
            StridedMultiIterator<3, int> i3 = ...;
                
            i3.dim2()++;  // increment outer dimension
            i3.dim0()++;  // increment inner dimension
            \endcode            
        */
    template <unsigned int K>
    StridedMultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }

    StridedMultiIterator<1, T, REFERENCE, POINTER> &
    dim0() { return *this; }
    StridedMultiIterator<2, T, REFERENCE, POINTER> &
    dim1() { return *this; }
    StridedMultiIterator<3, T, REFERENCE, POINTER> &
    dim2() { return *this; }
    StridedMultiIterator<4, T, REFERENCE, POINTER> &
    dim3() { return *this; }
    StridedMultiIterator<5, T, REFERENCE, POINTER> &
    dim4() { return *this; }

  protected:

    difference_type 
    total_stride(typename multi_difference_type::const_iterator d) const
    {
        return d[level]*this->m_stride[level] + base_type::total_stride(d);
    }

};

namespace detail {

template <unsigned int M>
struct MoveToScanOrderId
{
    template <class Shape, class Ptr>
    static void 
    exec(MultiArrayIndex & newId, MultiArrayIndex shapeStride, Ptr & p,
         Shape & point, Shape const & shape, Shape const & strides)
    {
        enum { N = Shape::static_size };
        MoveToScanOrderId<M-1>::exec(newId, shapeStride*shape[N-1-M], p, point, shape, strides);
        MultiArrayIndex newPos = newId / shapeStride;
        p += (newPos - point[N-1-M]) * strides[N-1-M];
        point[N-1-M] = newPos;
        newId %= shapeStride;
    }
};

template <>
struct MoveToScanOrderId<0>
{
    template <class Shape, class Ptr>
    static void 
    exec(MultiArrayIndex & newId, MultiArrayIndex shapeStride, Ptr & p,
         Shape & point, Shape const & shape, Shape const & strides)
    {
        enum { N = Shape::static_size }; 
        MultiArrayIndex newPos = newId / shapeStride;
        p += (newPos - point[N-1]) * strides[N-1];
        point[N-1] = newPos;
        newId %= shapeStride;
    }
};

}



template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M = N>
class StridedScanOrderIterator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: protected StridedScanOrderIterator<N, T, REFERENCE, POINTER, M-1>
#endif
{
    typedef StridedScanOrderIterator<N, T, REFERENCE, POINTER, M-1> base_type;
    enum { level = M-1 };

  public:

    typedef typename base_type::value_type value_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::shape_type shape_type;
    typedef MultiArrayIndex difference_type;
    typedef StridedScanOrderIterator iterator;
    typedef std::random_access_iterator_tag iterator_category;
    
    StridedScanOrderIterator(pointer i, 
                             shape_type const & shape, shape_type const & strides)
    : base_type(i, shape, strides)
    {}

    StridedScanOrderIterator & operator++()
    {
        base_type::operator++();
        if(this->point_[level-1] == this->shape_[level-1] && 
           this->point_[level] < this->shape_[level]) 
        {
            base_type::reset();
            this->i_ += this->strides_[level];
            ++this->point_[level];
        }
        return *this;
    }

        /** Advance to next starting location.
         */
    StridedScanOrderIterator operator++(int)
    {
        StridedScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    StridedScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToScanOrderId(this->id_+i);
        return *this;
    }

    StridedScanOrderIterator & operator--()
    {
        base_type::operator--();
        if(this->point_[level-1] == -1 && this->point_[level] > 0) 
        {
            base_type::inverseReset();
            this->i_ -= this->strides_[level];
            --this->point_[level];
        }
        return *this;
    }

        /** Advance to next starting location.
         */
    StridedScanOrderIterator operator--(int)
    {
        StridedScanOrderIterator res(*this);
        --*this;
        return res;
    }

    StridedScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    StridedScanOrderIterator getEndIterator() const
    {
        StridedScanOrderIterator res(*this);
        res.makeEndIterator();
        return res;
    }

    bool atBorder() const
    {
        return base_type::atBorder() || point_[level] == 0 || point_[level] == shape_[level] - 1;
    }
    
    using base_type::point;
    using base_type::id;
    using base_type::operator*;
    using base_type::operator->;
    using base_type::operator[];

  protected:
    void reset()
    {
        this->i_ -= this->shape_[level]*this->strides_[level];
        this->point_[level] = 0;
    }

    void inverseReset()
    {
        this->i_ += this->shape_[level]*this->strides_[level];
        this->point_[level] = this->shape_[level]-1;
    }
};

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline StridedScanOrderIterator<N, T, REFERENCE, POINTER, M>
operator+(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> i, MultiArrayIndex d)
{
    return i += d;
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline StridedScanOrderIterator<N, T, REFERENCE, POINTER, M>
operator-(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> i, MultiArrayIndex d)
{
    return i -= d;
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline MultiArrayIndex
operator-(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
           StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() - r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator==(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
            StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() == r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator!=(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
            StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() != r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator<(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
            StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() < r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator<=(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
            StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() <= r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator>(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
           StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() > r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER, unsigned int M>
inline bool
operator>=(StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & l, 
           StridedScanOrderIterator<N, T, REFERENCE, POINTER, M> const & r)
{
    return l.id() >= r.id();
}

template <unsigned int N, class T, class REFERENCE, class POINTER>
class StridedScanOrderIterator<N, T, REFERENCE, POINTER, 1>
{
    enum { level = 0 };

  public:

    typedef T value_type;
    typedef POINTER pointer;
    typedef T const * const_pointer;
    typedef REFERENCE reference;
    typedef T const & const_reference;
    typedef typename MultiArrayShape<N>::type shape_type;
    typedef MultiArrayIndex difference_type;
    typedef StridedScanOrderIterator iterator;
    typedef std::random_access_iterator_tag iterator_category;

    StridedScanOrderIterator(pointer i, 
                             shape_type const & shape, shape_type const & strides)
    : i_(i),
      shape_(shape),
      strides_(strides),
      id_(0)
    {}

    StridedScanOrderIterator & operator++()
    {
        i_ += strides_[level];
        ++point_[level];
        ++id_;
        return *this;
    }

    StridedScanOrderIterator operator++(int)
    {
        StridedScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    StridedScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToScanOrderId(id_+i);
        return *this;
    }
    
    StridedScanOrderIterator & operator--()
    {
        i_ -= strides_[level];
        --point_[level];
        --id_;
        return *this;
    }

    StridedScanOrderIterator operator--(int)
    {
        StridedScanOrderIterator res(*this);
        --this;
        return res;
    }

    StridedScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }
    
    reference operator*()
    {
        return *i_;
    }

    const_reference operator*() const
    {
        return *i_;
    }
    
    pointer operator->()
    {
        return i_;
    }

    const_pointer operator->() const
    {
        return i_;
    }

    reference operator[](MultiArrayIndex i)
    {
        StridedScanOrderIterator t(*this);
        t.moveToScanOrderId(id_+i);
        return *t;
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        StridedScanOrderIterator t(*this);
        t.moveToScanOrderId(id_+i);
        return *t;
    }

    bool atBorder() const
    {
        return point_[level] == 0 || point_[level] == shape_[level] - 1;
    }
    
    MultiArrayIndex id() const
    {
        return id_;
    }
    
    shape_type const & point() const
    {
        return point_;
    }
    
    StridedScanOrderIterator getEndIterator() const
    {
        StridedScanOrderIterator res(*this);
        res.makeEndIterator();
        return res;
    }

  protected:
    void reset()
    {
        i_ -= shape_[level]*strides_[level];
        point_[level] = 0;
    }
    
    void inverseReset()
    {
        i_ += shape_[level]*strides_[level];
        point_[level] = shape_[level] - 1;
    }
    
    void makeEndIterator()
    {
        moveToScanOrderId(prod(shape_));
    }
    
    void moveToScanOrderId(MultiArrayIndex newId)
    {
        id_ = newId;
        detail::MoveToScanOrderId<N-1>::exec(newId, 1, i_, point_, shape_, strides_);
    }

    pointer i_;
    shape_type point_, shape_, strides_;
    MultiArrayIndex id_;
};


//@}

} // namespace vigra

#endif // VIGRA_MULTI_ITERATOR_HXX
