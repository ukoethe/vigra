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


#ifndef VIGRA_MULTI_ITERATOR_HXX
#define VIGRA_MULTI_ITERATOR_HXX

#include <sys/types.h>
#include "vigra/tinyvector.hxx"
#include "vigra/iteratortags.hxx"

namespace vigra {


template <unsigned int N, class T, 
          class REFERENCE = T &, class POINTER = T *> class MultiIterator;
template <unsigned int N, class T, 
          class REFERENCE = T &, class POINTER = T *> class StridedMultiIterator;

/** \page MultiIteratorPage  Multi-dimensional Array Iterators

    General iterators for arrays of arbitrary dimension.


<p>
<DL>
<DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
    \ref vigra::MultiIterator
    <DD> <em>Iterator for unstrided \ref vigra::MultiArrayView</em>
<DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
    \ref vigra::MultiIteratorBase::type
    <DD> <em>Inner class implementing most of the functionality of \ref vigra::MultiIterator</em>
<DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
    \ref vigra::StridedMultiIterator
    <DD> <em>Iterator for strided \ref vigra::MultiArrayView</em>
<DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
    \ref vigra::StridedMultiIteratorBase::type
    <DD> <em>Inner class implementing most of the functionality of \ref vigra::StridedMultiIterator</em>
</DL>
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
    The current dimension of a <tt>MultiIterator&lt;N, ..&gt;</tt> is <tt>N-1</tt>.
</p>
<h3>Gerneral Requirements for MultiIterator</h3>
<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr><td>
    \htmlonly
    <th colspan=2>
    \endhtmlonly
    Local Types
    \htmlonly
    </th><th>
    \endhtmlonly
    Meaning
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
<tr><td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::value_type</tt></td><td>the underlying arrays's pixel type</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::reference</tt></td>
    <td>the iterator's reference type (return type of <TT>*iter</TT>). Will be
    <tt>value_type &</tt> for a mutable iterator, and convertible to
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::pointer</tt></td>
    <td>the iterator's pointer type (return type of <TT>iter.operator->()</TT>). Will be
    <tt>value_type *</tt> for a mutable iterator, and convertible to
    <tt>value_type const *</tt> for a const iterator.</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::iterator_category</tt></td>
    <td>the iterator tag (<tt>vigra::multi_dimensional_traverser_tag</tt>)</td>
</tr>
<tr><td>
    \htmlonly
    <th>
    \endhtmlonly
    Operation
    \htmlonly
    </th><th>
    \endhtmlonly
    Result
    \htmlonly
    </th><th>
    \endhtmlonly
    Semantics
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
<tr><td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator k;</tt></td><td>default constructor</td>
</tr>
<tr><td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
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
<tr><td>
    \htmlonly
    <th colspan=2>
    \endhtmlonly
    Local Types
    \htmlonly
    </th><th>
    \endhtmlonly
    Meaning
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
<tr><td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::multi_difference_type</tt></td>
    <td>the iterator's multi-dimensional difference type (<TT>TinyVector&lt;ptrdiff_t, N&gt;</TT>)</td>
</tr>
<tr><td>
    \htmlonly
    <th>
    \endhtmlonly
    Operation
    \htmlonly
    </th><th>
    \endhtmlonly
    Result
    \htmlonly
    </th><th>
    \endhtmlonly
    Semantics
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
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
    <td><tt>i.dim<K>()</tt></td><td><tt>MultiIterator&lt;K+1, T, ...&gt;</tt></td>
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
<tr><td>
    \htmlonly
    <td colspan=3>
    \endhtmlonly
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
<tr><td>
    \htmlonly
    <th colspan=2>
    \endhtmlonly
    Local Types
    \htmlonly
    </th><th>
    \endhtmlonly
    Meaning
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::difference_type</tt></td>
    <td>the iterator's difference type (<TT>ptrdiff_t</TT>)</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>MultiIterator::next_type</tt></td><td>type of the next iterator
       (referring to the next lower dimension) in the hierarchy</td>
</tr>
<tr><td>
    \htmlonly
    <th>
    \endhtmlonly
    Operation
    \htmlonly
    </th><th>
    \endhtmlonly
    Result
    \htmlonly
    </th><th>
    \endhtmlonly
    Semantics
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
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
<tr>
    <td>
    \htmlonly
    <td colspan=3>
    \endhtmlonly
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


/********************************************************/
/*                                                      */
/*                    MultiIteratorBase                 */
/*                                                      */
/********************************************************/

/** \brief Enclosing class for \ref vigra::MultiIterator base classes.

This design is necessary for compilers that do not support partial
specialization (otherwise, MultiIterator could be specialized directly).

<b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

Namespace: vigra
*/
template <unsigned int N>
class MultiIteratorBase
{
  public:
        /** \brief Base class for \ref vigra::MultiIterator.

        This class implements the multi-iterator by means of the enclosed template 
        class <tt>type</tt>. This design is necessary for compilers that do not support partial
        specialization (otherwise, MultiIterator could be specialized directly).

        <b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

        Namespace: vigra
        */
    template <class T, class REFERENCE, class POINTER>
    class type : public MultiIterator <N-1, T, REFERENCE, POINTER>
    {
    public:
            /** the type of the parent in the inheritance hierarchy.
             */
        typedef MultiIterator <N-1, T, REFERENCE, POINTER> base_type;

    public:
        
            /** the iterator's level in the dimension hierarchy
             */
        enum { level = N-1 };

            /** the iterator's value type
             */
        typedef T value_type;

            /** reference type (result of operator[] and operator*())
             */
        typedef REFERENCE reference;

            /** pointer type (result of operator->())
             */
        typedef POINTER pointer;

            /** difference type (used for offsetting along one axis)
             */
        typedef ptrdiff_t difference_type;

            /** multi difference type 
                (used for offsetting along all axes simultaneously)
             */
        typedef TinyVector<difference_type, N> multi_difference_type;

            /** the next type, this is a non-standard typedef denoting the
                type of the multi-iterator with the next-lower dimension.
            */
        typedef MultiIterator <level, T, REFERENCE, POINTER> next_type;

            /** the 1-dimensional iterator for this iterator hierarchy
                (result of iteratorForDimension()).
            */
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;

            /** the iterator tag (image traverser)
            */
        typedef multi_dimensional_traverser_tag iterator_category;

            
        /* use default copy constructor and assignment operator */
            
            /** default constructor.
            */
        type ()
        {}

            /** construct from pointer, strides (offset of a sample to the
                next) for every dimension, and the shape.
            */
        type (pointer ptr,
              const difference_type *stride,
              const difference_type *shape)
            : base_type (ptr, stride, shape)
        {}

            /** prefix-increment the iterator in it's current dimension
             */
        void operator++ ()
        {
            m_ptr += m_stride [level];
        }

            /** prefix-decrement the iterator in it's current dimension
             */
        void operator-- ()
        {
            m_ptr -= m_stride [level];
        }

            /** postfix-increment the iterator in it's current dimension
             */
        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

            /** postfix-decrement the iterator in it's current dimension
             */
        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

            /** increment the iterator in it's current dimension
                by the given value.
            */
        type & operator+= (difference_type n)
        {
            m_ptr += n * m_stride [level];
            return *this;
        }
        
            /** increment the iterator in all dimensions
                by the given offset.
            */
        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += total_stride(d.begin());
            return *this;
        }

            /** decrement the iterator in it's current dimension
                by the given value.
            */
        type & operator-= (difference_type n)
        {
            m_ptr -= n * m_stride [level];
            return *this;
        }
        
            /** decrement the iterator in all dimensions
                by the given offset.
            */
        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= total_stride(d.begin());
            return *this;
        }
        
            /** difference of two iterators in the current dimension.
                The result of this operation is undefined if the iterator
                doesn't point to element 0 in all dimensions below its current dimension.
            */
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr) / m_stride[level];
        }

        /* operators *, ->, ==, !=, < inherited */

            /** access the array element at the given offset in 
            the current dimension.
            */
        reference operator[] (difference_type n) const
        {
            return m_ptr [n* m_stride [level]];
        }

            /** access the array element at the given offset.
            */
        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [total_stride(d.begin())];
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
            ret += m_shape [level-1];
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
            return iterator(m_ptr, &m_stride [d], 0);
        }
        
      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level]*m_stride[level] + base_type::total_stride(d);
        }
    };
};

/********************************************************/
/*                                                      */
/*                      MultiIteratorBase <2>           */
/*                                                      */
/********************************************************/

//
template <>
class MultiIteratorBase <2>
{
  public:
    template <class T, class REFERENCE, class POINTER>
    class type : public MultiIterator <1, T, REFERENCE, POINTER>
    {

    public:

        enum { level = 1 };
        typedef MultiIterator <1, T, REFERENCE, POINTER> base_type;
        typedef T value_type;
        typedef REFERENCE reference;
        typedef const value_type &const_reference;
        typedef POINTER pointer;
        typedef const value_type *const_pointer;
        typedef ptrdiff_t difference_type;
        typedef MultiIterator <1, T, REFERENCE, POINTER> next_type;
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;
        typedef TinyVector<difference_type, 2> multi_difference_type;
        typedef std::random_access_iterator_tag iterator_category;

        const difference_type *m_stride;
        const difference_type *m_shape;
            
        /* use default copy constructor and assignment operator */
            
        type ()
            : base_type (),
              m_stride (0), m_shape (0)
        {}

        type (pointer ptr,
              const difference_type *stride,
              const difference_type *shape)
            : base_type (ptr, stride, shape),
              m_stride (stride), m_shape (shape)
        {}

        void operator++ ()
        {
            m_ptr += m_stride [level];
        }

        void operator-- ()
        {
            m_ptr -= m_stride [level];
        }

        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

        type & operator+= (difference_type n)
        {
            m_ptr += n * m_stride [level];
            return *this;
        }
        
        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += total_stride(d.begin());
            return *this;
        }

        type  &operator-= (difference_type n)
        {
            m_ptr -= n * m_stride [level];
            return *this;
        }
        
        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= total_stride(d.begin());
            return *this;
        }
        
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr) / m_stride[level];
        }

        reference operator[] (difference_type n) const
        {
            return m_ptr [n*m_stride [level]];
        }

        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [total_stride(d.begin())];
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
            return iterator(m_ptr, &m_stride [d], 0);
        }
        
      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level]*m_stride[level] + base_type::total_stride(d);
        }
    };
};

/********************************************************/
/*                                                      */
/*                   MultiIteratorBase <1>              */
/*                                                      */
/********************************************************/

//
template <>
class MultiIteratorBase <1>
{
  public:
    template <class T, class REFERENCE, class POINTER>
    class type
    {
      public:
        enum { level = 0 };
        typedef T value_type;
        typedef REFERENCE reference;
        typedef const value_type &const_reference;
        typedef POINTER pointer;
        typedef const value_type *const_pointer;
        typedef ptrdiff_t difference_type;
        typedef void next_type;
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;
        typedef TinyVector<difference_type, 1> multi_difference_type;
        typedef multi_dimensional_traverser_tag iterator_category;

        pointer m_ptr;

            
        /* use default copy constructor and assignment operator */
            
        type ()
            : m_ptr (0)
        {}
        
        type (pointer ptr,
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

        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

        type &operator+= (difference_type n)
        {
            m_ptr += n;
            return *this;
        }

        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += d[level];
            return *this;
        }

        type &operator-= (difference_type n)
        {
            m_ptr -= n;
            return *this;
        }

        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= d[level];
            return *this;
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

        reference operator[] (difference_type n) const
        {
            return m_ptr [n];
        }

        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [d[level]];
        }
        
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr);
        }

        bool operator!= (const type &rhs) const
        {
            return m_ptr != rhs.m_ptr;
        }

        bool operator== (const type &rhs) const
        {
            return m_ptr == rhs.m_ptr;
        }

        bool operator< (const type &rhs) const
        {
            return m_ptr < rhs.m_ptr;
        }

        bool operator<= (const type &rhs) const
        {
            return m_ptr <= rhs.m_ptr;
        }

        iterator iteratorForDimension(unsigned int d) const
        {
            vigra_precondition(d == 0,
                "MultiIterator<1>::iteratorForDimension(d): d == 0 required");
            const difference_type stride = 1;
            return iterator(m_ptr, &stride, 0);
        }
        
      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level];
        }
    };
};


/********************************************************/
/*                                                      */
/*                      MultiIterator                   */
/*                                                      */
/********************************************************/

/** \brief A multi-dimensional hierarchical iterator to be used with 
           \ref vigra::MultiArrayView if it is not strided.

This class wraps the MultiIteratorBase in a template of arity two.

<b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

Namespace: vigra
*/
template <unsigned int N, class T, class REFERENCE, class POINTER>
class MultiIterator 
: public MultiIteratorBase <N>::template type <T, REFERENCE, POINTER>
{
public:

        /** the type of the parent in the inheritance hierarchy.
         */
    typedef typename MultiIteratorBase <N>::template type <T, REFERENCE, POINTER> base_type;

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

        /** difference type (used for offsetting)
         */
    typedef ptrdiff_t difference_type;

        /** multi difference type 
            (used for offsetting along all axes simultaneously)
         */
    typedef TinyVector<difference_type, N> multi_difference_type;
    
        /** the MultiIterator for the next lower dimension.
         */
    typedef typename base_type::next_type next_type;

        /** the 1-dimensional iterator for this iterator hierarchy
            (result of iteratorForDimension()).
        */
    typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;

        /** the iterator tag (image traverser)
        */
    typedef typename base_type::iterator_category iterator_category;
            
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
        return base_type::operator-(d);
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

        /** Return the multi-iterator that operates on dimension K in order
            to manipulate this dimension directly. Usage:
               
            \code
                
            MultiIterator<3, int> i3 = ...;
                
            i3.dim<2>()++;  // increment outer dimension
            i3.dim<0>()++;  // increment inner dimension
            \endcode
        */
    template <unsigned int K>
    MultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*               StridedMultiIteratorBase               */
/*                                                      */
/********************************************************/

/** \brief Encloses the base class for \ref vigra::StridedMultiIterator.

This design is necessary for compilers that do not support partial
specialization (otherwise, StridedMultiIterator could be specialized directly).

<b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

Namespace: vigra
*/
template <unsigned int N>
class StridedMultiIteratorBase
{
  public:
        /** \brief Base class for \ref vigra::StridedMultiIterator.

        This class implements the multi-iterator for strided arrays
        by means of the enclosed template 
        class <tt>type</tt>. This design is necessary for compilers that do not support partial
        specialization (otherwise, MultiIterator could be specialized directly).

        <b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

        Namespace: vigra
        */
    template <class T, class REFERENCE, class POINTER>
    class type : public StridedMultiIterator <N-1, T, REFERENCE, POINTER>
    {
    public:

            /** the type of the parent in the inheritance hierarchy.
             */
        typedef StridedMultiIterator <N-1, T, REFERENCE, POINTER> base_type;

    public:
        
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

            /** difference type (used for offsetting)
             */
        typedef ptrdiff_t difference_type;

            /** multi difference type 
                (used for offsetting along all axes simultaneously)
             */
        typedef TinyVector<difference_type, N> multi_difference_type;

            /** the next type, this is a non-standard typedef denoting the
                type of the multi-iterator with the next-lower dimension.
            */
        typedef StridedMultiIterator <level, T, REFERENCE, POINTER> next_type;

            /** the 1-dimensional iterator for this iterator hierarchy
                (result of iteratorForDimension()).
            */
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;

            /** the iterator tag (image traverser)
            */
        typedef multi_dimensional_traverser_tag iterator_category;
            
        /* use default copy constructor and assignment operator */
            
            /** default constructor.
            */
        type ()
        {}

            /** construct from pointer, strides (offset of a sample to the
                next) for every dimension, and the shape.
            */
        type (pointer ptr,
              const difference_type *stride,
              const difference_type *shape)
            : base_type (ptr, stride, shape)
        {}

            /** prefix-increment the iterator in it's current dimension
             */
        void operator++ ()
        {
            m_ptr += m_stride [level];
        }

            /** prefix-decrement the iterator in it's current dimension
             */
        void operator-- ()
        {
            m_ptr -= m_stride [level];
        }

            /** postfix-increment the iterator in it's current dimension
             */
        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

            /** postfix-decrement the iterator in it's current dimension
             */
        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

            /** increment the iterator in it's current dimension
                by the given value.
            */
        type &operator+= (difference_type n)
        {
            m_ptr += n * m_stride [level];
            return *this;
        }

            /** increment the iterator in all dimensions
                by the given offset.
            */
        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += total_stride(d.begin());
            return *this;
        }

            /** decrement the iterator in it's current dimension
                by the given value.
            */
        type &operator-= (difference_type n)
        {
            m_ptr -= n * m_stride [level];
            return *this;
        }

            /** decrement the iterator in all dimensions
                by the given offset.
            */
        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= total_stride(d.begin());
            return *this;
        }
        
            /** difference of two iterators in the current dimension.
                The result of this operation is undefined if the iterator
                doesn't point to element 0 in all dimensions below its current dimension.
            */
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr) / m_stride[level];
        }

        /* operators *, ->, ==, !=, < inherited */

            /** access the array element at the given offset
                in the iterator's current dimension.
            */
        reference operator[] (difference_type n) const
        {
            return m_ptr [n* m_stride [level]];
        }

            /** access the array element at the given offset.
            */
        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [total_stride(d.begin())];
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
                otherwise it is undefined. Usage:
            */
        next_type end () const
        {
            next_type ret = *this;
            ret += m_shape [level-1];
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
            vigra_precondition(d <= N,
                "StridedMultiIterator<N>::iteratorForDimension(d): d <= N required");
            return iterator(m_ptr, &m_stride [d], 0);
        }

      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level]*m_stride[level] + base_type::total_stride(d);
        }
    };
};

/********************************************************/
/*                                                      */
/*              StridedMultiIteratorBase <2>            */
/*                                                      */
/********************************************************/

//
template <>
class StridedMultiIteratorBase <2>
{
  public:
    template <class T, class REFERENCE, class POINTER>
    class type : public StridedMultiIterator <1, T, REFERENCE, POINTER>
    {
    public:
        enum { level = 1 };
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> base_type;
        typedef T value_type;
        typedef REFERENCE reference;
        typedef const value_type &const_reference;
        typedef POINTER pointer;
        typedef const value_type *const_pointer;
        typedef ptrdiff_t difference_type;
        typedef TinyVector<difference_type, 2>      multi_difference_type;
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> next_type;
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;
        typedef multi_dimensional_traverser_tag iterator_category;

        const difference_type *m_stride;
        const difference_type *m_shape;
            
        /* use default copy constructor and assignment operator */
            
        type ()
            : base_type (),
              m_stride (0), m_shape (0)
        {}

        type (pointer ptr,
              const difference_type *stride,
              const difference_type *shape)
            : base_type (ptr, stride, shape),
              m_stride (stride), m_shape (shape)
        {}

        void operator++ ()
        {
            m_ptr += m_stride [level];
        }

        void operator-- ()
        {
            m_ptr -= m_stride [level];
        }

        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

        type &operator+= (int n)
        {
            m_ptr += n * m_stride [level];
            return *this;
        }

        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += total_stride(d.begin());
            return *this;
        }

        type  &operator-= (difference_type n)
        {
            m_ptr -= n * m_stride [level];
            return *this;
        }
        
        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= total_stride(d.begin());
            return *this;
        }

        reference operator[] (difference_type n) const
        {
            return m_ptr [n*m_stride [level]];
        }
        
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr) / m_stride[level];
        }

        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [total_stride(d.begin())];
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
            vigra_precondition(d <= N,
                "StridedMultiIterator<N>::iteratorForDimension(d): d <= N required");
            return iterator(m_ptr, &m_stride [d], 0);
        }

      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level]*m_stride[level] + base_type::total_stride(d);
        }
    };
};

/********************************************************/
/*                                                      */
/*               StridedMultiIteratorBase <1>           */
/*                                                      */
/********************************************************/

//
template <>
class StridedMultiIteratorBase <1>
{
  public:
    template <class T, class REFERENCE, class POINTER>
    class type
    {
    public:

        enum { level = 0 };
        typedef T value_type;
        typedef REFERENCE reference;
        typedef const value_type &const_reference;
        typedef POINTER pointer;
        typedef const value_type *const_pointer;
        typedef ptrdiff_t difference_type;
        typedef TinyVector<difference_type, 1> multi_difference_type;
        typedef void next_type;
        typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;
        typedef std::random_access_iterator_tag iterator_category;
        
        pointer m_ptr;
        difference_type m_stride;
            
        /* use default copy constructor and assignment operator */
            
        type ()
            : m_ptr (0), m_stride (0)
        {}

        type (pointer ptr,
              const difference_type *stride,
              const difference_type *)
            : m_ptr (ptr), m_stride (stride [level])
        {}

        reference operator* () const
        {
            return *m_ptr;
        }

        pointer get () const
        {
            return m_ptr;
        }

        pointer operator-> () const
        {
            return &(operator*());
        }

        void operator++ ()
        {
            m_ptr += m_stride;
        }

        void operator-- ()
        {
            m_ptr -= m_stride;
        }

        type operator++ (int)
        {
            type ret = *this;
            ++(*this);
            return ret;
        }

        type operator-- (int)
        {
            type ret = *this;
            --(*this);
            return ret;
        }

        type &operator+= (difference_type n)
        {
            m_ptr += n * m_stride;
            return *this;
        }

        type & operator+= (multi_difference_type const & d)
        {
            m_ptr += d[level] * m_stride;
            return *this;
        }

        type &operator-= (difference_type n)
        {
            m_ptr -= n * m_stride;
            return *this;
        }
        
        type & operator-= (multi_difference_type const & d)
        {
            m_ptr -= d[level] * m_stride;
            return *this;
        }
        
        difference_type operator- (type const & d) const
        {
            return (d.m_ptr - m_ptr) / m_stride;
        }

        reference operator[] (difference_type n) const
        {
            return m_ptr [n*m_stride];
        }

        reference operator[] (multi_difference_type const & d) const
        {
            return m_ptr [d[level]*m_stride];
        }

        bool operator!= (const type &rhs) const
        {
            return m_ptr != rhs.m_ptr;
        }

        bool operator== (const type &rhs) const
        {
            return m_ptr == rhs.m_ptr;
        }

        bool operator< (const type &rhs) const
        {
            return m_ptr < rhs.m_ptr;
        }

        bool operator<= (const type &rhs) const
        {
            return m_ptr <= rhs.m_ptr;
        }
        
        iterator iteratorForDimension(unsigned int d) const
        {
            vigra_precondition(d == 0,
                "StridedMultiIterator<1>::iteratorForDimension(d): d == 0 required");
            return *this;
        }

      protected:
      
        difference_type 
        total_stride(typename multi_difference_type::const_iterator d) const
        {
            return d[level]*m_stride;
        }
    };
};

/********************************************************/
/*                                                      */
/*                 StridedMultiIterator                 */
/*                                                      */
/********************************************************/

/** \brief A multi-dimensional hierarchical iterator to be used with 
           \ref vigra::MultiArrayView is it is strided.

This class wraps the StridedMultiIteratorBase in a template of arity two.

<b>\#include</b> "<a href="multi_iterator_8hxx-source.html">vigra/multi_iterator.hxx</a>"

Namespace: vigra
*/
template <unsigned int N, class T, class REFERENCE, class POINTER>
class StridedMultiIterator
    : public StridedMultiIteratorBase <N>::template type <T, REFERENCE, POINTER>
{
public:

        /** the type of the parent in the inheritance hierarchy.
         */
    typedef typename StridedMultiIteratorBase <
        N>::template type <T, REFERENCE, POINTER> base_type;

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

        /** difference type (used for offsetting)
         */
    typedef ptrdiff_t difference_type;

        /** multi difference type 
            (used for offsetting along all axes simultaneously)
         */
    typedef TinyVector<difference_type, N> multi_difference_type;
    
        /** the StridedMultiIterator for the next lower dimension.
         */
    typedef typename base_type::next_type next_type;

        /** the 1-dimensional iterator for this iterator hierarchy
            (result of iteratorForDimension()).
        */
    typedef StridedMultiIterator <1, T, REFERENCE, POINTER> iterator;

        /** the iterator tag (image traverser)
        */
    typedef typename base_type::iterator_category iterator_category;
            
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
        return base_type::operator-(d);
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

        /** Return the multi-iterator that operates on dimension K in order
            to manipulate this dimension directly. Usage:
                
            \code
            
            StridedMultiIterator<3, int> i3 = ...;
                
            i3.dim<2>()++;  // increment outer dimension
            i3.dim<0>()++;  // increment inner dimension
            \endcode
        */
    template <unsigned int K>
    StridedMultiIterator<K+1, T, REFERENCE, POINTER> &
    dim()
    {
        return *this;
    }
};


//@}

} // namespace vigra

#endif // VIGRA_MULTI_ITERATOR_HXX
