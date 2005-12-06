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

#ifndef VIGRA_BYTEORDER_HXX
#define VIGRA_BYTEORDER_HXX

#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <algorithm>
#include "vigra/sized_int.hxx"

namespace vigra
{
    class byteorder
    {
        class host
        {
            std::string m_string;

        public:
            
            // ctor, dtor

            host();

            // methods

            const std::string & get() const;
        };

        // attributes

        static const host m_host;

        // data byte order, can be changed

        std::string m_string;
        bool native;

        // delegation methods

        template< class T >
        void reversebytes( T & x ) const
        {
            const unsigned int n = sizeof(T);
            UInt8 t[n];
            UInt8 * c = reinterpret_cast< UInt8 * >(&x);
            unsigned int i;
            for( i = 0; i < n; ++i )
                t[i] = c[ n - 1 - i ];
            for( i = 0; i < n; ++i )
                c[i] = t[i];
        }

    public:

        // ctor, dtor

        byteorder(); // uses the host byteorder
        byteorder( const std::string & );

        // methods

        void set( const std::string & );
        const std::string & get() const;
        const std::string & get_host_byteorder() const;

        template< class T >
        void convert_to_host( T & x ) const
        {
            if (!native)
                reversebytes(x);
        }

        template< class T >
        void convert_to_host( T * x, unsigned int num ) const
        {
            if (!native)
                for( unsigned int i = 0; i < num; ++i )
                    reversebytes(x[i]);
        }

        template< class T >
        void convert_from_host( T & x ) const
        {
            if (!native)
                reversebytes(x);
        }

        template< class T >
        void convert_from_host( T * x, unsigned int num ) const
        {
            if (!native)
                for( unsigned int i = 0; i < num; ++i )
                    reversebytes(x[i]);
        }

        void convert_to_host( char & x ) const {}
        void convert_to_host( Int8 & x ) const {}
        void convert_to_host( UInt8 & x ) const {}

        void convert_to_host( char * x , unsigned int) const {}
        void convert_to_host( Int8 * x, unsigned int) const {}
        void convert_to_host( UInt8 * x, unsigned int) const {}

        void convert_from_host( char & x ) const {}
        void convert_from_host( Int8 & x ) const {}
        void convert_from_host( UInt8 & x ) const {}

        void convert_from_host( char * x , unsigned int) const {}
        void convert_from_host( Int8 * x, unsigned int) const {}
        void convert_from_host( UInt8 * x, unsigned int) const {}
    };

    template< class T >
    void read_field( std::ifstream & stream, const byteorder & bo, T & x )
    {
        stream.read( reinterpret_cast< char * >(&x), sizeof(T) );
        bo.convert_to_host(x);
    }

    template< class T >
    void read_array( std::ifstream & stream, const byteorder & bo, T * x,
                     unsigned int num )
    {
        stream.read( reinterpret_cast< char * >(x), sizeof(T) * num );
        bo.convert_to_host( x, num );
    }

    template< class T >
    void write_field( std::ofstream & stream, const byteorder & bo, T t )
    {
        bo.convert_from_host(t);
        stream.write( reinterpret_cast< char * >(&t), sizeof(T) );
    }

    template< class T >
    void write_array( std::ofstream & stream, const byteorder & bo,
                      const T * x, unsigned int num )
    {
        for( unsigned int i = 0; i < num; ++i )
            write_field( stream, bo, x[i] );
    }

} // namespace vigra

#endif // VIGRA_BYTEORDER_HXX
