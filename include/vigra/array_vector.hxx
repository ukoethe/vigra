/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2003 by Ullrich Koethe                  */
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

#ifndef VIGRA_ARRAY_VECTOR_HXX
#define VIGRA_ARRAY_VECTOR_HXX

#include <memory>
#include <algorithm>
#include <vigra/memory.hxx>

namespace vigra
{

template <class T>
class ArrayVector
{
    typedef ArrayVector<T> this_type;

public:
    typedef T value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef value_type * iterator;
    typedef value_type const * const_iterator;
    typedef unsigned int size_type;
    typedef int          difference_type;

public:
    ArrayVector();

    ArrayVector( size_type size);

    ArrayVector( size_type size, value_type const & initial);

    ArrayVector( this_type const & rhs );

    template <class InputIterator>
    ArrayVector(InputIterator i, InputIterator end);

    this_type & operator=( this_type const & rhs );

    ~ArrayVector();

    inline const_pointer data() const
    {
        return data_;
    }

    inline pointer data()
    {
        return data_;
    }

    inline const_iterator begin() const
    {
        return data();
    }

    inline iterator begin()
    {
        return data();
    }

    inline const_iterator end() const
    {
        return data() + size();
    }

    inline iterator end()
    {
        return data() + size();
    }

    reference front()
    {
        return *data_;
    }

    const_reference front() const
    {
        return *data_;
    }

    reference back()
    {
        return data_[size_-1];
    }

    const_reference back() const
    {
        return data_[size_-1];
    }

    reference operator[]( size_type i )
    {
        return data()[i];
    }

    const_reference operator[]( size_type i ) const
    {
        return data()[i];
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

    void reserve( size_type new_capacity );

    void reserve();

    void resize( size_type new_size, value_type const & initial );

    void resize( size_type new_size )
    {
        resize(new_size, value_type());
    }

    bool empty() const
    {
        return size_ == 0;
    }

    size_type size() const
    {
        return size_;
    }

    size_type capacity() const
    {
        return capacity_;
    }

    void swap(this_type & rhs);

  private:

    static void deallocate(pointer data, size_type size);

    static pointer reserve_raw(size_type capacity);

    size_type size_, capacity_;
    pointer data_;
};

template <class T>
ArrayVector<T>::ArrayVector()
: size_(0),
  capacity_(5),
  data_(reserve_raw(5))
{}

template <class T>
ArrayVector<T>::ArrayVector( size_type size)
: size_(size),
  capacity_(size),
  data_(reserve_raw(size))
{
    if(size_ > 0)
        std::uninitialized_fill(data_, data_+size_, value_type());
}

template <class T>
ArrayVector<T>::ArrayVector( size_type size, value_type const & initial)
: size_(size),
  capacity_(size),
  data_(reserve_raw(size))
{
    if(size_ > 0)
        std::uninitialized_fill(data_, data_+size_, initial);
}

template <class T>
ArrayVector<T>::ArrayVector( this_type const & rhs )
: size_(rhs.size_),
  capacity_(rhs.capacity_),
  data_(reserve_raw(rhs.capacity_))
{
    if(size_ > 0)
        std::uninitialized_copy(rhs.data_, rhs.data_+size_, data_);
}

template <class T>
template <class InputIterator>
ArrayVector<T>::ArrayVector(InputIterator i, InputIterator end)
: size_(std::distance(i, end)),
  capacity_(size_),
  data_(reserve_raw(size_))
{
    std::uninitialized_copy(i, end, data_);
}


template <class T>
ArrayVector<T> & ArrayVector<T>::operator=( this_type const & rhs )
{
    if(this == &rhs)
        return *this;
    ArrayVector new_vector(rhs);
    swap(new_vector);
    return *this;
}

template <class T>
ArrayVector<T>::~ArrayVector()
{
    deallocate(data_, size_);
}

template <class T>
void ArrayVector<T>::pop_back()
{
    --size_;
    detail::destroy(data_ + size_);
}

template <class T>
void ArrayVector<T>::push_back( value_type const & t )
{
    reserve();
    new (static_cast<void*>(data_ + size_)) T(t);
    ++size_;
}

template <class T>
void ArrayVector<T>::clear()
{
    detail::destroy_n(data_, size_);
    size_ = 0;
}

template <class T>
typename ArrayVector<T>::iterator
ArrayVector<T>::insert(iterator p, value_type const & v)
{
    difference_type pos = p - begin();
    if(p == end())
    {
        push_back(v);
        p = begin() + pos;
    }
    else
    {
        push_back(back());
        p = begin() + pos;
        std::copy_backward(p, end() - 2, end() - 1);
        *p = v;
    }
    return p;
}

template <class T>
typename ArrayVector<T>::iterator
ArrayVector<T>::insert(iterator p, size_type n, value_type const & v)
{
    difference_type pos = p - begin();
    size_type new_size = size() + n;
    if(new_size >= capacity_)
    {
        pointer new_data = reserve_raw(new_size);
        std::uninitialized_copy(begin(), p, new_data);
        std::uninitialized_fill(new_data + pos, new_data + pos + n, v);
        std::uninitialized_copy(p, end(), new_data + pos + n);
        deallocate(data_, size_);
        capacity_ = new_size;
        data_ = new_data;
    }
    else if(pos + n >= size_)
    {
        size_type diff = pos + n - size_;
        std::uninitialized_copy(p, end(), end() + diff);
        std::uninitialized_fill(end(), end() + diff, v);
        std::fill(p, end(), v);
    }
    else
    {
        size_type diff = size_ - (pos + n);
        std::uninitialized_copy(end() - n, end(), end());
        std::copy_backward(p, p + diff, end());
        std::fill(p, p + n, v);
    }
    size_ = new_size;
    return begin() + pos;
}

template <class T>
template <class InputIterator>
typename ArrayVector<T>::iterator
ArrayVector<T>::insert(iterator p, InputIterator i, InputIterator iend)
{
    difference_type n = iend - i;
    difference_type pos = p - begin();
    size_type new_size = size() + n;
    if(new_size >= capacity_)
    {
        pointer new_data = reserve_raw(new_size);
        std::uninitialized_copy(begin(), p, new_data);
        std::uninitialized_copy(i, iend, new_data + pos);
        std::uninitialized_copy(p, end(), new_data + pos + n);
        std::deallocate(data_, size_);
        capacity_ = new_size;
        data_ = new_data;
    }
    else if(pos + n >= size_)
    {
        size_type diff = pos + n - size_;
        std::uninitialized_copy(p, end(), end() + diff);
        std::uninitialized_copy(iend - diff, iend, end());
        std::copy(i, iend - diff, p);
    }
    else
    {
        size_type diff = size_ - (pos + n);
        std::uninitialized_copy(end() - n, end(), end());
        std::copy_backward(p, p + diff, end());
        std::copy(i, iend, p);
    }
    size_ = new_size;
    return begin() + pos;
}

template <class T>
typename ArrayVector<T>::iterator
ArrayVector<T>::erase(iterator p)
{
    std::copy(p+1, end(), p);
    pop_back();
    return p;
}

template <class T>
typename ArrayVector<T>::iterator
ArrayVector<T>::erase(iterator p, iterator q)
{
    std::copy(q, end(), p);
    size_type eraseCount = q - p;
    detail::destroy_n(end() - eraseCount, eraseCount);
    size_ -= eraseCount;
    return p;
}

template <class T>
void ArrayVector<T>::reserve( size_type new_capacity )
{
    if(new_capacity <= capacity_)
        return;
    pointer new_data = reserve_raw(new_capacity);
    if(size_ > 0)
        std::uninitialized_copy(data_, data_+size_, new_data);
    deallocate(data_, size_);
    data_ = new_data;
    capacity_ = new_capacity;
}

template <class T>
void ArrayVector<T>::reserve()
{
    if(size_ == capacity_)
        reserve(2*capacity_);
}

template <class T>
void ArrayVector<T>::resize( size_type new_size, value_type const & initial)
{
    if(new_size < size_)
        erase(begin() + new_size, end());
    else if(size_ < new_size)
    {
        insert(end(), new_size - size(), initial);
    }
}

template <class T>
void ArrayVector<T>::swap(this_type & rhs)
{
    std::swap(size_, rhs.size_);
    std::swap(capacity_, rhs.capacity_);
    std::swap(data_, rhs.data_);
}

template <class T>
void ArrayVector<T>::deallocate(pointer data, size_type size)
{
    if(data)
    {
        detail::destroy_n(data, size);
        ::operator delete(data);
    }
}

template <class T>
typename ArrayVector<T>::pointer
ArrayVector<T>::reserve_raw(size_type capacity)
{
    pointer data = 0;
    if(capacity)
        data = static_cast<pointer>(::operator new(capacity*sizeof(value_type)));
    return data;
}

} // namespace vigra


#endif /* VIGRA_ARRAY_VECTOR_HXX */
