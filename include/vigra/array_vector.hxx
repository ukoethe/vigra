/************************************************************************/
/*                                                                      */
/*               Copyright 200-2003 by Ullrich Koethe                   */
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

namespace vigra
{

template< class T >
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
    ArrayVector()
    : size_(0),
      capacity_(5),
      data_(reserve_raw(5)),
    {}

    ArrayVector( size_type size )
    : size_(size),
      capacity_(size),
      data_(reserve_raw(size)),
    {
        if(size_ > 0)
            std::uninitialized_fill(data_, data_+size_, value_type());
    }

    ArrayVector( this_type const & rhs )
    : size_(rhs.size_),
      capacity_(rhs.capacity_),
      data_(reserve_raw(rhs.capacity_)),
    {
        if(size_ > 0)
            std::uninitialized_copy(rhs.data_, rhs.data_+size_, data_);
    }

    this_type & operator=( this_type const & rhs )
    {
        if(this == &rhs)
            return *this;
        ArrayVector new_vector(rhs);
        this->swap(new_vector);
        return *this;
    }

    ~ArrayVector()
    {
        deallocate(data_, size_);
    }

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

    void push_back( value_type const & t )
    {
        if ( size_ == capacity_ ) reserve();
        data()[m_size++] = t;
    }

    reference operator[]( size_type i )
    {
        return data()[i];
    }

    const_reference operator[]( size_type i ) const
    {
        return data()[i];
    }

    void reserve( size_type new_capacity )
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

    void reserve()
    {
        reserve(2*capacity_);
    }

    void resize( size_type new_size )
    {
        if(new_size < size_)
            std::destroy(data+new_size, data+size_, (value_type*)0);
        else if(size_ < new_size && new_size <= capacity_)
            std::uninitialized_fill(data_+size_, data_+new_size, value_type());
        else if(capacity_ < new_size)
        {
            reserve(new_size);
            std::uninitialized_fill(data_+size_, data_+new_size, value_type());
        }
        size_ = new_size;
    }

    size_type size() const
    {
        return size_;
    }

    size_type capacity() const
    {
        return capacity_;
    }
    
    void swap(this_type & rhs)
    {
        std::swap(size_, rhs.size_);
        std::swap(capacity_, rhs.capacity_);
        std::swap(data_, rhs.data_);
    }
    
  private:
    static void deallocate(pointer data, size_type size)
    {
        if(data)
        {
            std::destroy(data, data+size, (value_type*)0);
            ::operator delete(data);
        }
    }
  
    static pointer reserve_raw(size_type capacity)
    {
        pointer data = 0;
        if(capacity)
            data = ::operator new(capacity*sizeof(value_type));
        return data;
    }
  
    size_type size_, capacity_;
    pointer data_;
};



} // namespace vigra


#endif /* VIGRA_ARRAY_VECTOR_HXX */
