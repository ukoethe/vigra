/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2012 by Gunnar Kedenburg and Ullrich Koethe       */
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
#include "iteratortags.hxx"
#include "multi_iterator_coupled.hxx"

namespace vigra {

/** \addtogroup MultiIteratorGroup
*/
//@{

    /** \brief Iterate over a virtual array where each element contains its coordinate.

        MultiCoordinateIterator behaves like a read-only random access iterator. 
        It moves accross the given region of interest in scan-order (with the first
        index changing most rapidly), and dereferencing the iterator returns the 
        coordinate (i.e. multi-dimensional index) of the current array element. 
        The functionality is thus similar to a meshgrid in Matlab or numpy. 
        
        Internally, it is just a wrapper of a \ref CoupledScanOrderIterator that
        has been created without any array and whose reference type is not a 
        \ref CoupledHandle, but the coordinate itself.
                
        The iterator supports all functions listed in the STL documentation for 
        <a href="http://www.sgi.com/tech/stl/RandomAccessIterator.html">Random Access Iterators</a>.

        <b>Usage:</b>

        <b>\#include</b> \<vigra/multi_iterator.hxx\><br/>
        Namespace: vigra
        
        \code
        MultiCoordinateIterator<3> i(Shape3(3,2,1)), end = i.getEndIterator();
        
        for(; i != end; ++i)
            std::cout << *i << "\n";
            
        // Output:
        // (0, 0, 0)
        // (1, 0, 0)
        // (2, 0, 0)
        // (0, 1, 0)
        // (1, 1, 0)
        // (2, 1, 0)
        \endcode
    */
template<unsigned int N>
class MultiCoordinateIterator
    : public CoupledScanOrderIterator<N>
{
  public:
    typedef CoupledScanOrderIterator<N>            base_type;

    typedef typename base_type::shape_type         shape_type;
    typedef typename base_type::difference_type    difference_type;
    typedef MultiCoordinateIterator                iterator;
    typedef std::random_access_iterator_tag        iterator_category;

    typedef typename base_type::value_type         handle_type;
    typedef typename handle_type::value_type       value_type;
    typedef typename handle_type::reference        reference;
    typedef typename handle_type::const_reference  const_reference;
    typedef typename handle_type::pointer          pointer;
    typedef typename handle_type::const_pointer    const_pointer;

    MultiCoordinateIterator() 
        : base_type(handle_type())
    {}

    MultiCoordinateIterator(shape_type const & shape) 
        : base_type(handle_type(shape))
    {}

    // dereferencing the iterator yields the coordinate object
    // (used as vertex_descriptor)
    reference operator*()
    {
        return this->template get<0>();
    }
    
    const_reference operator*() const
    {
        return this->template get<0>();
    }

    pointer operator->()
    {
        return &this->template get<0>();
    }
    
    const_pointer operator->() const
    {
        return &this->template get<0>();
    }

    reference operator[](MultiArrayIndex i)
    {
        return *(MultiCoordinateIterator(*this) += i);
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        return *(MultiCoordinateIterator(*this) += i);
    }

    MultiCoordinateIterator & operator++()
    {
        base_type::operator++();
        return *this;
    }
    
    MultiCoordinateIterator operator++(int)
    {
        MultiCoordinateIterator res(*this);
        ++*this;
        return res;
    }

    MultiCoordinateIterator & operator+=(MultiArrayIndex i)
    {
        base_type::operator+=(i);
        return *this;
    }

    MultiCoordinateIterator & operator+=(const shape_type &coordOffset)
    {
        base_type::operator+=(coordOffset);
        return *this;
    }

    MultiCoordinateIterator & operator--()
    {
        base_type::operator--();
        return *this;
    }

    MultiCoordinateIterator operator--(int)
    {
        MultiCoordinateIterator res(*this);
        --*this;
        return res;
    }

    MultiCoordinateIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    MultiCoordinateIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

    MultiCoordinateIterator getEndIterator() const
    {
        return MultiCoordinateIterator(base_type::getEndIterator());
    }

    MultiCoordinateIterator operator+(MultiArrayIndex d) const
    {
        return MultiCoordinateIterator(*this) += d;
    }

    MultiCoordinateIterator operator-(MultiArrayIndex d) const
    {
        return MultiCoordinateIterator(*this) -= d;
    }

    MultiCoordinateIterator operator+(const shape_type &coordOffset) const
    {
        return MultiCoordinateIterator(*this) += coordOffset;
    }

    MultiCoordinateIterator operator-(const shape_type &coordOffset) const
    {
        return MultiCoordinateIterator(*this) -= coordOffset;
    }

    MultiArrayIndex operator-(const MultiCoordinateIterator & other) const
    {
        return base_type::operator-(other);
    }
    
  protected:
    MultiCoordinateIterator(base_type const & base) 
        : base_type(base)
    {}
};

    /** \brief Sequential iterator for MultiArrayView.
        
        This iterator provides STL-compatible random access iterator functionality for arbitrary 
        \ref MultiArrayView instances, regardless of their shapes and strides. The
        class uses an implementation that minimizes speed penalties that could result from 
        non-trivial strides. The <i>scan-order</i> is defined such that dimensions are iterated 
        from front to back (first to last).
        
        You normally construct instances of this class by calling \ref MultiArrayView::begin() 
        and \ref MultiArrayView::end(). 
        
        The iterator supports all functions listed in the STL documentation for 
        <a href="http://www.sgi.com/tech/stl/RandomAccessIterator.html">Random Access Iterators</a>.
        
        <b>\#include</b> \<vigra/multi_iterator.hxx\><br/>
        Namespace: vigra
    */
template <unsigned int N, class T, class REFERENCE, class POINTER>
class StridedScanOrderIterator
    : public CoupledIteratorType<N, T>::type
{
  public:
    typedef typename CoupledIteratorType<N, T>::type  base_type;

    typedef typename base_type::shape_type         shape_type;
    typedef typename base_type::difference_type    difference_type;
    typedef StridedScanOrderIterator              iterator;
    typedef std::random_access_iterator_tag        iterator_category;

    typedef T                                      value_type;
    typedef REFERENCE                              reference;
    typedef T const &                              const_reference;
    typedef POINTER                                pointer;
    typedef T const *                              const_pointer;

    StridedScanOrderIterator() 
        : base_type()
    {}

    template <class S>
    StridedScanOrderIterator(MultiArrayView<N, T, S> const & view) 
        : base_type(createCoupledIterator(view))
    {}

    StridedScanOrderIterator(POINTER p, shape_type const & shape, shape_type const & strides) 
        : base_type(createCoupledIterator(MultiArrayView<N, T, StridedArrayTag>(shape, strides, const_cast<T *>(p))))
    {}

    reference operator*()
    {
        return this->template get<1>();
    }
    
    const_reference operator*() const
    {
        return this->template get<1>();
    }

    pointer operator->()
    {
        return &this->template get<1>();
    }
    
    const_pointer operator->() const
    {
        return &this->template get<1>();
    }

    reference operator[](MultiArrayIndex i)
    {
        return *(StridedScanOrderIterator(*this) += i);
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        return *(StridedScanOrderIterator(*this) += i);
    }

    StridedScanOrderIterator & operator++()
    {
        base_type::operator++();
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
        base_type::operator+=(i);
        return *this;
    }

    StridedScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        base_type::operator+=(coordOffset);
        return *this;
    }

    StridedScanOrderIterator & operator--()
    {
        base_type::operator--();
        return *this;
    }

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

    StridedScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

    StridedScanOrderIterator getEndIterator() const
    {
        return StridedScanOrderIterator(base_type::getEndIterator());
    }

    StridedScanOrderIterator operator+(MultiArrayIndex d) const
    {
        return StridedScanOrderIterator(*this) += d;
    }

    StridedScanOrderIterator operator-(MultiArrayIndex d) const
    {
        return StridedScanOrderIterator(*this) -= d;
    }

    MultiArrayIndex operator-(StridedScanOrderIterator const & other) const
    {
        return base_type::operator-(other);
    }

    StridedScanOrderIterator operator+(const shape_type &coordOffset) const
    {
        return StridedScanOrderIterator(*this) += coordOffset;
    }

    StridedScanOrderIterator operator-(const shape_type &coordOffset) const
    {
        return StridedScanOrderIterator(*this) -= coordOffset;
    }
    
    MultiArrayIndex index() const
    {
        return this->scanOrderIndex();
    }
    
  protected:
    StridedScanOrderIterator(base_type const & base) 
        : base_type(base)
    {}
};

//@}

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
<LI> \ref vigra::StridedScanOrderIterator
     <BR>&nbsp;&nbsp;&nbsp;<em>STL-compatible random access iterator for \ref vigra::MultiArrayView</em>
<LI> \ref vigra::CoupledScanOrderIterator
     <BR>&nbsp;&nbsp;&nbsp;<em>Iterate over multiple images simultaneously in scan order</em>
</UL>
</p>

<p>
    The Multidimensional Iterator concept allows navigation on arrays
    of arbitrary dimension. It provides two modes of iteration: 
    <em>direct traversal</em>, and <em>hierarchical traversal</em>.
    In general, hierarchical traversal will be faster, while only 
    direct traversal allows for true random access in all dimensions.
    Via the <tt>dim<K>()</tt> function, operations applying to a particular
    dimension can be used in the direct traversal mode. In contrast,
    direct traversal functions should not be used in the hierarchical mode
    because the hierarchical functions are only well-defined if the
    iterator points to element 0 in all dimensions below its current dimension.
    The current dimension of a <tt>MultiIterator<N, ...></tt> is <tt>N-1</tt>.
</p>
<h3>General Requirements for MultiIterator</h3>
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
    <td>create the hierarchical iterator pointing to the first element in the 
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
    <td>create the hierarchical iterator pointing to the past-the-end location in the 
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

/** \addtogroup MultiIteratorGroup  
*/
//@{

template <class POINTER>
struct MultiIteratorStrideTraits
{
    typedef MultiArrayIndex    stride_type;
    typedef const stride_type* stride_array_type;
    typedef stride_array_type  shape_array_type;
    static stride_array_type shift(stride_array_type s, unsigned d)
    {
        return s + d;
    }
};

template <unsigned int N, class T, class REFERENCE = T &, class POINTER = T *>
class MultiIterator;

template <unsigned int N, class T, class REFERENCE = T &, class POINTER = T *>
class StridedMultiIterator;

/********************************************************/
/*                                                      */
/*                      MultiIterator                   */
/*                                                      */
/********************************************************/

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
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
    typedef typename stride_traits::shape_array_type shape_array_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef std::random_access_iterator_tag iterator_category;

  protected:
    pointer m_ptr;

  public:
    MultiIterator ()
        : m_ptr (0)
    {}

    MultiIterator (pointer ptr,
          const difference_array_type &,
          const shape_array_type &)
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
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
    typedef typename stride_traits::shape_array_type shape_array_type;
    typedef base_type next_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef multi_dimensional_traverser_tag iterator_category;

  protected:
    difference_array_type m_stride;
    shape_array_type      m_shape;

  public:
    /* use default copy constructor and assignment operator */

    MultiIterator ()
        : base_type (),
          m_stride (0), m_shape (0)
    {}

    MultiIterator (pointer ptr,
          const difference_array_type & stride,
          const shape_array_type & shape)
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
        return iterator(this->m_ptr, stride_traits::shift(m_stride, d), 0);
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
#ifndef DOXYGEN
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
    typedef typename stride_traits::shape_array_type shape_array_type;
#else
    typedef MultiArrayIndex difference_type;
#endif

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
                   const difference_array_type & stride,
                   const shape_array_type & shape)
        : base_type (ptr, stride, shape)
    {}


        /** prefix-increment the iterator in its current dimension
         */
    void operator++ ()
    {
        this->m_ptr += this->m_stride [level];
    }

        /** prefix-decrement the iterator in its current dimension
         */
    void operator-- ()
    {
        this->m_ptr -= this->m_stride [level];
    }

        /** postfix-increment the iterator in its current dimension
         */
    MultiIterator operator++ (int)
    {
        MultiIterator ret = *this;
        ++(*this);
        return ret;
    }

        /** postfix-decrement the iterator in its current dimension
         */
    MultiIterator operator-- (int)
    {
        MultiIterator ret = *this;
        --(*this);
        return ret;
    }

        /** increment the iterator in its current dimension
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

        /** decrement the iterator in its current dimension
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
        return iterator(this->m_ptr, stride_traits::shift(this->m_stride, d),0);
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
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
    typedef typename stride_traits::shape_array_type shape_array_type;
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
          const difference_array_type & stride,
          const shape_array_type &)
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
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
    typedef typename stride_traits::shape_array_type shape_array_type;
    typedef base_type next_type;
    typedef StridedMultiIterator<1, T, REFERENCE, POINTER> iterator;
    typedef multi_dimensional_traverser_tag iterator_category;

  protected:
    difference_array_type m_stride;
    shape_array_type m_shape;

  public:
    /* use default copy constructor and assignment operator */

    StridedMultiIterator ()
        : base_type (),
          m_stride (0), m_shape (0)
    {}

    StridedMultiIterator (pointer ptr,
          const difference_array_type & stride,
          const shape_array_type & shape)
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
        return iterator(this->m_ptr, stride_traits::shift(m_stride, d), 0);
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
#ifndef DOXYGEN
    typedef MultiIteratorStrideTraits<POINTER> stride_traits;
    typedef typename stride_traits::stride_type difference_type;
    typedef typename stride_traits::stride_array_type difference_array_type;
#else
    typedef MultiArrayIndex difference_type;
#endif
    
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
                   const difference_array_type & stride,
                   const difference_array_type & shape)
        : base_type (ptr, stride, shape)
    {}


        /** prefix-increment the iterator in its current dimension
         */
    void operator++ ()
    {
        this->m_ptr += this->m_stride [level];
    }

        /** prefix-decrement the iterator in its current dimension
         */
    void operator-- ()
    {
        this->m_ptr -= this->m_stride [level];
    }

        /** postfix-increment the iterator in its current dimension
         */
    StridedMultiIterator operator++ (int)
    {
        StridedMultiIterator ret = *this;
        ++(*this);
        return ret;
    }

        /** postfix-decrement the iterator in its current dimension
         */
    StridedMultiIterator operator-- (int)
    {
        StridedMultiIterator ret = *this;
        --(*this);
        return ret;
    }

        /** increment the iterator in its current dimension
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

        /** decrement the iterator in its current dimension
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
        return iterator(this->m_ptr, stride_traits::shift(this->m_stride, d),0);
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

//@}

} // namespace vigra

#endif // VIGRA_MULTI_ITERATOR_HXX
