/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_ARRAY_VECTOR_HXX
#define VIGRA_ARRAY_VECTOR_HXX

#include "error.hxx"
#include "memory.hxx"
#include "numerictraits.hxx"
#include <memory>
#include <algorithm>
#include <iosfwd>

#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(diff >= 0, "Index out of bounds");\
  vigra_precondition(diff < (difference_type)size_, "Index out of bounds");
#else
#define VIGRA_ASSERT_INSIDE(diff)
#endif

namespace vigra
{

template <class T, class Alloc = std::allocator<T> >
class ArrayVector;

/** Provide STL conforming interface for C-arrays.

    This template implements much of the functionality of <a href="http://www.sgi.com/tech/stl/Vector.html">std::vector</a>
    on top of a C-array. <tt>ArrayVectorView</tt> does not manage the memory
    it refers to (i.e. it does not allocate or deallocate any memory).
    Thus, if the underlying memory changes, all dependent <tt>ArrayVectorView</tt>
    objects are invalidated. This is especially important when <tt>ArrayVectorView</tt>
    is used as a base class for <tt>ArrayVector</tt>, where several functions
    (e.g. resize(), insert()) can allocate new memory and thus invalidate the
    dependent views. The rules what operations invalidate view objects are the
    same as the rules concerning standard iterators.

    <b>\#include</b> \<vigra/array_vector.hxx\><br>
    Namespace: vigra
*/
template <class T>
class ArrayVectorView
{
    typedef ArrayVectorView<T> this_type;

public:
        /** default constructor
        */
    typedef T value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef value_type * iterator;
    typedef value_type const * const_iterator;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

public:
        /** default constructor.
            View contains NULL pointer.
        */
    ArrayVectorView()
    : size_(0),
      data_(0)
    {}

        /** Construct for given array \a data of length \a size.
            <tt>data, data+size</tt> must form a valid range.
        */
    ArrayVectorView( size_type size, pointer const & data)
    : size_(size),
      data_(data)
    {}

        /** Copy constructor.
        */
    ArrayVectorView( this_type const & rhs )
    : size_(rhs.size_),
      data_(rhs.data_)
    {}

        /** Copy assignment. There are 3 cases:

            <ul>
            <li> When this <tt>ArrayVectorView</tt> does not point to valid data
                 (e.g. after default construction), it becomes a copy of \a rhs.
            <li> When the shapes of the two arrays match, the array contents
                 (not the pointers) are copied.
            <li> Otherwise, a <tt>PreconditionViolation</tt> exception is thrown.
            </ul>
        */
    ArrayVectorView & operator=( ArrayVectorView const & rhs );

        /** Copy assignment.
            When the shapes of the two arrays match, the array contents
            (not the pointers) are copied. Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    template <class U>
    this_type & operator=( ArrayVectorView<U> const & rhs )
    {
        copyImpl(rhs);
        return *this;
    }

        /** Overwrite all array elements with the value \a initial.
        */
    template <class U>
    void init(U const & initial)
    {
        std::fill(begin(), end(), initial);
    }

        /** Copy array elements.
            When the shapes of the two arrays match, the array contents
            (not the pointers) are copied. Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    void copy( this_type const & rhs )
    {
        if(data_ != rhs.data_)
            copyImpl(rhs);
    }

        /** Copy array elements.
            When the shapes of the two arrays match, the array contents
            (not the pointers) are copied. Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    template <class U>
    void copy( ArrayVectorView<U> const & rhs )
    {
        copyImpl(rhs);
    }

        /** Swap array elements.
            When the shapes of the two arrays match, the array contents
            (not the pointers) are swapped. Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    void swapData(this_type rhs)
    {
        if(data_ != rhs.data_)
            swapDataImpl(rhs);
    }

        /** Swap array elements.
            When the shapes of the two arrays match, the array contents
            (not the pointers) are swapped. Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    template <class U>
    void swapData(ArrayVectorView<U> rhs)
    {
        swapDataImpl(rhs);
    }

        /** Construct <tt>ArrayVectorView</tt> referring to a subarray.
            \a begin and \a end must be a valid sub-range of the current array.
            Otherwise, a <tt>PreconditionViolation</tt>
            exception is thrown.
        */
    this_type subarray (size_type begin, size_type end) const
    {
        vigra_precondition(begin <= end && end <= size_,
                "ArrayVectorView::subarray(): Limits out of range.");
        return this_type(end-begin, data_ + begin);
    }

        /** Get contained const pointer to the data.
        */
    inline const_pointer data() const
    {
        return data_;
    }

        /** Get contained pointer to the data.
        */
    inline pointer data()
    {
        return data_;
    }

        /** Get const iterator referring to the first array element.
        */
    inline const_iterator begin() const
    {
        return data();
    }

        /** Get iterator referring to the first array element.
        */
    inline iterator begin()
    {
        return data();
    }

        /** Get const iterator pointing beyond the last array element.
        */
    inline const_iterator end() const
    {
        return data() + size();
    }

        /** Get iterator pointing beyond the last array element.
        */
    inline iterator end()
    {
        return data() + size();
    }

        /** Get const iterator referring to the first array element.
        */
    inline const_iterator cbegin() const
    {
        return data();
    }

        /** Get const iterator pointing beyond the last array element.
        */
    inline const_iterator cend() const
    {
        return data() + size();
    }

        /** Get reverse iterator referring to the last array element.
        */
    inline reverse_iterator rbegin()
    {
        return (reverse_iterator(end()));
    }

        /** Get const reverse iterator referring to the last array element.
        */
    inline const_reverse_iterator rbegin() const
    {
        return (const_reverse_iterator(end()));
    }

        /** Get reverse iterator pointing before the first array element.
        */
    inline reverse_iterator rend()
    {
        return (reverse_iterator(begin()));
    }

        /** Get const reverse iterator pointing before the first array element.
        */
    inline const_reverse_iterator rend() const
    {
        return (const_reverse_iterator(begin()));
    }

        /** Get const reverse iterator referring to the last array element.
        */
    inline const_reverse_iterator crbegin() const
    {
        return (const_reverse_iterator(end()));
    }

        /** Get const reverse iterator pointing before the first array element.
        */
    inline const_reverse_iterator crend() const
    {
        return (const_reverse_iterator(begin()));
    }

        /** Access first array element.
        */
    reference front()
    {
        return *data_;
    }

        /** Read first array element.
        */
    const_reference front() const
    {
        return *data_;
    }

        /** Access last array element.
        */
    reference back()
    {
        return data_[size_-1];
    }

        /** Read last array element.
        */
    const_reference back() const
    {
        return data_[size_-1];
    }

        /** Access array element \a i.
        */
    reference operator[]( difference_type i )
    {
        VIGRA_ASSERT_INSIDE(i);
        return data()[i];
    }

        /** Read array element \a i.
        */
    const_reference operator[]( difference_type i ) const
    {
        VIGRA_ASSERT_INSIDE(i);
        return data()[i];
    }

        /** Equivalent to <tt>size() == 0</tt>.
        */
    bool empty() const
    {
        return size_ == 0;
    }

        /** Number of elements in the array.
        */
    size_type size() const
    {
        return size_;
    }

        /** Check for element-wise equality of two array.
            Also returns <tt>false</tt> if the two arrays have different sizes.
        */
    template <class U>
    bool operator==(ArrayVectorView<U> const & rhs) const;

        /** check whether two arrays are not elementwise equal.
            Also returns <tt>true</tt> if the two arrays have different sizes.
         */
    template <class U>
    bool operator!=(ArrayVectorView<U> const & rhs) const
    {
        return !operator==(rhs);
    }

        /** check whether the given point is in the array range.
         */
    bool isInside (difference_type const & p) const
    {
        return p >= 0 && p < size_;
    }

  protected:

    template <class U>
    void copyImpl(const ArrayVectorView <U>& rhs);

    void copyImpl(const ArrayVectorView & rhs);

    template <class U>
    void swapDataImpl(const ArrayVectorView <U>& rhs);

    size_type size_;
    pointer data_;
};

template <class T>
ArrayVectorView<T> & ArrayVectorView<T>::operator=( ArrayVectorView<T> const & rhs )
{
    if(data_ == 0)
    {
        size_ = rhs.size_;
        data_ = rhs.data_;
    }
    else if(data_ != rhs.data_)
        copyImpl(rhs);
    return *this;
}

template <class T>
template <class U>
bool ArrayVectorView<T>::operator==(ArrayVectorView<U> const & rhs) const
{
    if(size() != rhs.size())
        return false;
    for(size_type k=0; k<size(); ++k)
        if(data_[k] != rhs[k])
            return false;
    return true;
}

template <class T>
void
ArrayVectorView <T>::copyImpl(const ArrayVectorView & rhs)
{
    vigra_precondition (size() == rhs.size(),
        "ArrayVectorView::copy(): shape mismatch.");
    if(size() == 0)  // needed because MSVC debug assertions in std::copy() may fire
        return;      // "invalid address: data_ == NULL" even when nothing is to be copied
    // use copy() or copy_backward() according to possible overlap of this and rhs
    if(data_ <= rhs.data())
    {
        std::copy(rhs.begin(), rhs.end(), begin());
    }
    else
    {
        std::copy_backward(rhs.begin(), rhs.end(), end());
    }
}

template <class T>
template <class U>
void
ArrayVectorView <T>::copyImpl(const ArrayVectorView <U>& rhs)
{
    vigra_precondition (size() == rhs.size(),
        "ArrayVectorView::copy(): shape mismatch.");
    std::copy(rhs.begin(), rhs.end(), begin());
}

template <class T>
template <class U>
void
ArrayVectorView <T>::swapDataImpl(const ArrayVectorView <U>& rhs)
{
    vigra_precondition (size () == rhs.size() (),
        "ArrayVectorView::swapData(): size mismatch.");

    // check for overlap
    if(data_ + size_ <= rhs.data_ || rhs.data_ + size_ <= data_)
    {
        for(size_type k=0; k<size_; ++k)
            std::swap(data_[k], rhs.data_[k]);
    }
    else
    {
        ArrayVector<T> t(*this);
        copyImpl(rhs);
        rhs.copyImpl(t);
    }
}


/** Replacement for <tt>std::vector</tt>.

    This template implements the same functionality as <tt><a href="http://www.sgi.com/tech/stl/Vector.html">std::vector</a></tt> (see there for detailed documentation).
    However, it gives two useful guarantees, that <tt>std::vector</tt> fails
    to provide:

    <ul>
    <li>The memory is always allocated as one contiguous piece.</li>
    <li>The iterator is always a <TT>T *</TT> </li>
    </ul>

    This means that memory managed by <tt>ArrayVector</tt> can be passed
    to algorithms that expect raw memory. This is especially important
    when legacy or C code has to be called, but it is also useful for certain
    optimizations.

    Moreover, <tt>ArrayVector</tt> is derived from <tt>ArrayVectorView</tt> so that one
    can create views of the array (in particular, subarrays). This implies another
    important difference to <tt>std::vector</tt>: the indexing operator
    (<tt>ArrayVector::operator[]</tt>) takes <tt>signed</tt> indices. In this way,
    an <tt>ArrayVectorView</tt> can be used with negative indices:

    \code
    ArrayVector<int> data(100);
    ArrayVectorView<int> view = data.subarray(50, 100);

    view[-50] = 1; // valid access
    \endcode

    Refer to the documentation of <tt>std::vector</tt> for a detailed
    description of <tt>ArrayVector</tt> functionality.

    <b>\#include</b> \<vigra/array_vector.hxx\><br>
    Namespace: vigra
*/
template <class T, class Alloc /* = std::allocator<T> */ >
class ArrayVector
: public ArrayVectorView<T>
{
    typedef ArrayVector<T, Alloc> this_type;
    enum { minimumCapacity = 2, resizeFactor = 2 };

public:
    typedef ArrayVectorView<T> view_type;
    typedef typename view_type::value_type value_type;
    typedef typename view_type::reference reference;
    typedef typename view_type::const_reference const_reference;
    typedef typename view_type::pointer pointer;
    typedef typename view_type::const_pointer const_pointer;
    typedef typename view_type::iterator iterator;
    typedef typename view_type::const_iterator const_iterator;
    typedef typename view_type::size_type size_type;
    typedef typename view_type::difference_type difference_type;
    typedef typename view_type::reverse_iterator reverse_iterator;
    typedef typename view_type::const_reverse_iterator const_reverse_iterator;
    typedef Alloc        allocator_type;

public:
    ArrayVector()
    : view_type(),
      capacity_(minimumCapacity),
      alloc_(Alloc())
    {
        this->data_ = reserve_raw(capacity_);
    }

    explicit ArrayVector(Alloc const & alloc)
    : view_type(),
      capacity_(minimumCapacity),
      alloc_(alloc)
    {
        this->data_ = reserve_raw(capacity_);
    }

    explicit ArrayVector( size_type size, Alloc const & alloc = Alloc())
    : view_type(),
      alloc_(alloc)
    {
        initImpl(size, value_type(), VigraTrueType());
    }

    ArrayVector( size_type size, value_type const & initial, Alloc const & alloc = Alloc())
    : view_type(),
      alloc_(alloc)
    {
        initImpl(size, initial, VigraTrueType());
    }


    ArrayVector( this_type const & rhs )
    : view_type(),
      alloc_(rhs.alloc_)
    {
        initImpl(rhs.begin(), rhs.end(), VigraFalseType());
    }

    template <class U>
    explicit ArrayVector( ArrayVectorView<U> const & rhs, Alloc const & alloc = Alloc() )
    : view_type(),
      alloc_(alloc)
    {
        initImpl(rhs.begin(), rhs.end(), VigraFalseType());
    }

    template <class InputIterator>
    ArrayVector(InputIterator i, InputIterator end)
    {
        initImpl(i, end, typename NumericTraits<InputIterator>::isIntegral());
    }

    template <class InputIterator>
    ArrayVector(InputIterator i, InputIterator end, Alloc const & alloc)
    : alloc_(alloc)
    {
        initImpl(i, end, typename NumericTraits<InputIterator>::isIntegral());
    }

    this_type & operator=( this_type const & rhs )
    {
        if(this == &rhs)
            return *this;
        if(this->size_ == rhs.size_)
            this->copyImpl(rhs);
        else
        {
            ArrayVector t(rhs);
            this->swap(t);
        }
        return *this;
    }

    template <class U>
    this_type & operator=( ArrayVectorView<U> const & rhs);

    ~ArrayVector()
    {
        deallocate(this->data_, this->size_, this->capacity_);
    }

    void pop_back();

    void push_back( value_type const & t );

    iterator insert(iterator p, value_type const & v);

    iterator insert(iterator p, size_type n, value_type const & v);

    template <class InputIterator>
    iterator insert(iterator p, InputIterator i, InputIterator iend);

    iterator erase(iterator p);

    iterator erase(iterator p, iterator q);

    void clear();

    pointer reserveImpl( bool dealloc, size_type new_capacity );

    pointer reserveImpl( bool dealloc);

    void reserve()
    {
        reserveImpl(true);
    }

    void reserve( size_type new_capacity )
    {
        reserveImpl(true, new_capacity);
    }

    void resize( size_type new_size, value_type const & initial );

    void resize( size_type new_size )
    {
        resize(new_size, value_type());
    }

    size_type capacity() const
    {
        return capacity_;
    }

    void swap(this_type & rhs);

  private:

    void deallocate(pointer data, size_type size, size_type capacity);

    pointer reserve_raw(size_type capacity);

    void initImpl( size_type size, value_type const & initial, VigraTrueType /*isIntegral*/);

    template <class Iter>
    void initImpl( Iter i, Iter end, VigraFalseType /*isIntegral*/);

    template <class Iter>
    void initImpl( Iter i, Iter end, Error_NumericTraits_not_specialized_for_this_case)
    {
        initImpl(i, end, VigraFalseType());
    }

    size_type capacity_;
    Alloc alloc_;
};

template <class T, class Alloc>
template <class U>
ArrayVector<T, Alloc> & ArrayVector<T, Alloc>::operator=( ArrayVectorView<U> const & rhs )
{
    if(this->size_ == rhs.size())
        this->copyImpl(rhs);
    else
    {
        ArrayVector t(rhs);
        this->swap(t);
    }
    return *this;
}

template <class T, class Alloc>
inline void ArrayVector<T, Alloc>::pop_back()
{
    --this->size_;
    alloc_.destroy(this->data_ + this->size_);
}

template <class T, class Alloc>
inline void ArrayVector<T, Alloc>::push_back( value_type const & t )
{
    pointer old_data = reserveImpl(false);
    size_type old_capacity = this->capacity_;
    alloc_.construct(this->data_ + this->size_, t);
    // deallocate old data _after_ construction of new element, so that
    // 't' can refer to the old data as in 'push_back(front())'
    deallocate(old_data, this->size_, old_capacity);
    ++this->size_;
}

template <class T, class Alloc>
inline void ArrayVector<T, Alloc>::clear()
{
    detail::destroy_n(this->data_, this->size_);
    this->size_ = 0;
}

template <class T, class Alloc>
typename ArrayVector<T, Alloc>::iterator
ArrayVector<T, Alloc>::insert(iterator p, value_type const & v)
{
    difference_type pos = p - this->begin();
    if(p == this->end())
    {
        push_back(v);
        p = this->begin() + pos;
    }
    else
    {
        T lastElement = this->back();
        push_back(lastElement);
        p = this->begin() + pos;
        std::copy_backward(p, this->end() - 2, this->end() - 1);
        *p = v;
    }
    return p;
}

template <class T, class Alloc>
typename ArrayVector<T, Alloc>::iterator
ArrayVector<T, Alloc>::insert(iterator p, size_type n, value_type const & v)
{
    difference_type pos = p - this->begin();
    size_type new_size = this->size() + n;
    if(new_size > capacity_)
    {
        size_type new_capacity = std::max(new_size, resizeFactor*capacity_);
        pointer new_data = reserve_raw(new_capacity);
        try
        {
            std::uninitialized_copy(this->begin(), p, new_data);
            std::uninitialized_fill(new_data + pos, new_data + pos + n, v);
            std::uninitialized_copy(p, this->end(), new_data + pos + n);
        }
        catch(...)
        {
            alloc_.deallocate(new_data, new_capacity);
            throw;
        }
        deallocate(this->data_, this->size_, this->capacity_);
        capacity_ = new_capacity;
        this->data_ = new_data;
    }
    else if(pos + n > this->size_)
    {
        size_type diff = pos + n - this->size_;
        std::uninitialized_copy(p, this->end(), this->end() + diff);
        std::uninitialized_fill(this->end(), this->end() + diff, v);
        std::fill(p, this->end(), v);
    }
    else
    {
        size_type diff = this->size_ - (pos + n);
        std::uninitialized_copy(this->end() - n, this->end(), this->end());
        std::copy_backward(p, p + diff, this->end());
        std::fill(p, p + n, v);
    }
    this->size_ = new_size;
    return this->begin() + pos;
}

template <class T, class Alloc>
template <class InputIterator>
typename ArrayVector<T, Alloc>::iterator
ArrayVector<T, Alloc>::insert(iterator p, InputIterator i, InputIterator iend)
{
    size_type n = std::distance(i, iend);
    size_type pos = p - this->begin();
    size_type new_size = this->size() + n;
    if(new_size > capacity_)
    {
        size_type new_capacity = std::max(new_size, resizeFactor*capacity_);
        pointer new_data = reserve_raw(new_capacity);
        try
        {
            std::uninitialized_copy(this->begin(), p, new_data);
            std::uninitialized_copy(i, iend, new_data + pos);
            std::uninitialized_copy(p, this->end(), new_data + pos + n);
        }
        catch(...)
        {
            alloc_.deallocate(new_data, new_capacity);
            throw;
        }
        deallocate(this->data_, this->size_, this->capacity_);
        capacity_ = new_capacity;
        this->data_ = new_data;
    }
    else if(pos + n > this->size_)
    {
        size_type diff = pos + n - this->size_;
        std::uninitialized_copy(p, this->end(), this->end() + diff);
        InputIterator split = i;
        std::advance(split, n - diff);
        std::uninitialized_copy(split, iend, this->end());
        std::copy(i, split, p);
    }
    else
    {
        size_type diff = this->size_ - (pos + n);
        std::uninitialized_copy(this->end() - n, this->end(), this->end());
        std::copy_backward(p, p + diff, this->end());
        std::copy(i, iend, p);
    }
    this->size_ = new_size;
    return this->begin() + pos;
}

template <class T, class Alloc>
typename ArrayVector<T, Alloc>::iterator
ArrayVector<T, Alloc>::erase(iterator p)
{
    std::copy(p+1, this->end(), p);
    pop_back();
    return p;
}

template <class T, class Alloc>
typename ArrayVector<T, Alloc>::iterator
ArrayVector<T, Alloc>::erase(iterator p, iterator q)
{
    std::copy(q, this->end(), p);
    difference_type eraseCount = q - p;
    detail::destroy_n(this->end() - eraseCount, eraseCount);
    this->size_ -= eraseCount;
    return p;
}

template <class T, class Alloc>
typename ArrayVector<T, Alloc>::pointer
ArrayVector<T, Alloc>::reserveImpl( bool dealloc, size_type new_capacity)
{
    if(new_capacity <= capacity_)
        return 0;
    pointer new_data = reserve_raw(new_capacity),
            old_data = this->data_;
    if(this->size_ > 0)
        std::uninitialized_copy(old_data, old_data+this->size_, new_data);
    this->data_ = new_data;
    if(!dealloc)
    {
        this->capacity_ = new_capacity;
        return old_data;
    }
    deallocate(old_data, this->size_, this->capacity_);
    this->capacity_ = new_capacity;
    return 0;
}

template <class T, class Alloc>
inline typename ArrayVector<T, Alloc>::pointer
ArrayVector<T, Alloc>::reserveImpl(bool dealloc)
{
    if(capacity_ == 0)
        return reserveImpl(dealloc, minimumCapacity);
    else if(this->size_ == capacity_)
        return reserveImpl(dealloc, resizeFactor*capacity_);
    else
        return 0;
}

template <class T, class Alloc>
inline void
ArrayVector<T, Alloc>::resize( size_type new_size, value_type const & initial)
{
    if(new_size < this->size_)
        erase(this->begin() + new_size, this->end());
    else if(this->size_ < new_size)
    {
        insert(this->end(), new_size - this->size(), initial);
    }
}

template <class T, class Alloc>
inline void
ArrayVector<T, Alloc>::initImpl( size_type size, value_type const & initial, VigraTrueType /*isIntegral*/)
{
    this->size_ = size;
    capacity_ = size;
    this->data_ = reserve_raw(capacity_);
    if(this->size_ > 0)
        std::uninitialized_fill(this->data_, this->data_+this->size_, initial);
}

template <class T, class Alloc>
template <class Iter>
inline void
ArrayVector<T, Alloc>::initImpl( Iter i, Iter end, VigraFalseType /*isIntegral*/)
{
    this->size_ = std::distance(i, end);
    capacity_ = this->size_;
    this->data_ = reserve_raw(capacity_);
    if(this->size_ > 0)
        detail::uninitializedCopy(i, end, this->data_);
}

template <class T, class Alloc>
inline void
ArrayVector<T, Alloc>::swap(this_type & rhs)
{
    std::swap(this->size_, rhs.size_);
    std::swap(capacity_, rhs.capacity_);
    std::swap(this->data_, rhs.data_);
}

template <class T, class Alloc>
inline void
ArrayVector<T, Alloc>::deallocate(pointer data, size_type size, size_type capacity)
{
    if(data)
    {
        detail::destroy_n(data, size);
        alloc_.deallocate(data, capacity);
    }
}

template <class T, class Alloc>
inline typename ArrayVector<T, Alloc>::pointer
ArrayVector<T, Alloc>::reserve_raw(size_type capacity)
{
    pointer data = 0;
    if(capacity)
    {
        data = alloc_.allocate(capacity);
    }
    return data;
}

} // namespace vigra

namespace std {

template <class T>
ostream & operator<<(ostream & s, vigra::ArrayVectorView<T> const & a)
{
    for(std::size_t k=0; k<a.size()-1; ++k)
        s << a[k] << ", ";
    if(a.size())
            s << a.back();
    return s;
}

} // namespace std

#undef VIGRA_ASSERT_INSIDE
#endif /* VIGRA_ARRAY_VECTOR_HXX */
