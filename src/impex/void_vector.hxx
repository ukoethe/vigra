/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Gunnar Kedenburg                     */
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

#ifndef VIGRA_IMPEX_VOIDVECTOR
#define VIGRA_IMPEX_VOIDVECTOR

// use this only for pod types

#include <algorithm>
#include <cstring>

namespace vigra
{

    // void_vector_base does not support iteration or subscripting
    // (because the payload is not typed)
    class void_vector_base
    {
        friend void swap_void_vector( void_vector_base &, void_vector_base & );

        typedef void_vector_base this_type;

    protected:

	typedef size_t size_type;
        void * m_data;
        size_type m_size, m_capacity;

    public:
        void_vector_base()
            : m_data(0), m_size(0), m_capacity(0)
        {}

        void_vector_base( size_type size )
            : m_data(0), m_size(size), m_capacity(size)
        {
            if (size)
                m_data = ::operator new(size);
        }
    
        void_vector_base( this_type const & rhs )
            : m_size(rhs.m_size), m_capacity(rhs.m_capacity)
        {
            m_data = ::operator new(rhs.m_capacity);
            std::memcpy( m_data, rhs.m_data, m_size );
        }
    
        ~void_vector_base()
        {
            if(m_data)
                ::operator delete(m_data);
        }
    
        this_type & operator=( this_type const & rhs )
        {
            if ( this == &rhs ) return *this;
            ::operator delete(m_data);
            m_capacity = rhs.m_capacity;
            m_size = rhs.m_size;
            m_data = ::operator new(m_capacity);
            std::memcpy( m_data, rhs.m_data, m_size );
            return *this;
        }
    
        void append( this_type const & rhs )
        {
            size_type new_size = m_size + rhs.m_size;
            if ( m_capacity < new_size )
                reserve( m_capacity + rhs.m_capacity );
            std::memcpy( static_cast< unsigned char * >
                         (m_data) + m_size, rhs.m_data, rhs.m_size );
            m_size = new_size;
        }

        inline void reserve()
        {
            reserve( 2 * m_capacity );
        }

        void reserve( size_type new_capacity )
        {
            if ( new_capacity <= m_capacity ) return;
            void * new_data = ::operator new(new_capacity);
            std::memcpy( new_data, m_data, m_size );
            ::operator delete(m_data);
            m_data = new_data;
            m_capacity = new_capacity;
        }

        inline void resize( size_type new_size )
        {
            reserve(new_size);
            m_size = m_capacity;
        }

        inline const void * data() const
        {
            return m_data;
        }

        inline void * data()
        {
            return m_data;
        }

        inline size_type size() const
        {
            return m_size;
        }
    
        inline size_type capacity() const
        {
            return m_capacity;
        }
    
        inline void clear()
        {
            m_size = 0;
        }
    };


    template< class T >
    class void_vector : public void_vector_base
    {
        typedef void_vector<T> this_type;

    public:
        typedef T value_type;
        typedef value_type * iterator;
        typedef value_type const * const_iterator;
    
    public:
        void_vector()
            : void_vector_base( 20 * sizeof(value_type) )
        {}

        void_vector( size_type size )
            : void_vector_base( size * sizeof(value_type) )
        {}

        this_type & operator=( this_type const & rhs )
        {
            void_vector_base::operator=(rhs);
            return *this;
        }

        inline const value_type * data() const
        {
            return static_cast< const value_type * >
                ( void_vector_base::data() );
        }

        inline value_type * data()
        {
            return static_cast< value_type * >
                ( void_vector_base::data() );
        }

        inline const iterator begin() const
        {
            return data();
        }

        inline iterator begin()
        {
            return data();
        }

        inline const iterator end() const
        {
            return data() + size();
        }
    
        inline iterator end()
        {
            return data() + size();
        }
    
        void push_back( value_type const & t )
        {
            if ( size() == capacity() ) resize();
            data()[m_size++] = t;
        }
    
        value_type & operator[]( size_type i )
        {
            return data()[i];
        }
    
        value_type const & operator[]( size_type i ) const
        {
            return data()[i];
        }

        void reserve( size_type new_capacity )
        {
            void_vector_base::reserve( sizeof(value_type) * new_capacity );
        }

        void resize( size_type new_size )
        {
            void_vector_base::resize( sizeof(value_type) * new_size );
        }

        size_type size() const
        {
            return m_size / sizeof(value_type);
        }

        size_type capacity() const
        {
            return m_capacity / sizeof(value_type);
        }
    };

    void swap_void_vector( void_vector_base & left, void_vector_base & right );
}

#endif // VIGRA_IMPEX_VOIDVECTOR
