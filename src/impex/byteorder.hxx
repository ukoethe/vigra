/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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
            const size_t n = sizeof(T);
            UInt8 t[n];
            UInt8 * c = reinterpret_cast< UInt8 * >(&x);
            size_t i;
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
        void convert_to_host( T * x, size_t num ) const
        {
            if (!native)
                for( size_t i = 0; i < num; ++i )
                    reversebytes(x[i]);
        }

        template< class T >
        void convert_from_host( T & x ) const
        {
            if (!native)
                reversebytes(x);
        }

        template< class T >
        void convert_from_host( T * x, size_t num ) const
        {
            if (!native)
                for( size_t i = 0; i < num; ++i )
                    reversebytes(x[i]);
        }

        void convert_to_host( char & ) const {}
        void convert_to_host( Int8 & ) const {}
        void convert_to_host( UInt8 & ) const {}

        void convert_to_host( char * , size_t) const {}
        void convert_to_host( Int8 * , size_t) const {}
        void convert_to_host( UInt8 * , size_t) const {}

        void convert_from_host( char & ) const {}
        void convert_from_host( Int8 & ) const {}
        void convert_from_host( UInt8 & ) const {}

        void convert_from_host( char * , size_t) const {}
        void convert_from_host( Int8 * , size_t) const {}
        void convert_from_host( UInt8 * , size_t) const {}
    };

    template< class T >
    void read_field( std::ifstream & stream, const byteorder & bo, T & x )
    {
        stream.read( reinterpret_cast< char * >(&x), sizeof(T) );
        bo.convert_to_host(x);
    }

    template< class T >
    void read_array( std::ifstream & stream, const byteorder & bo, T * x,
                     size_t num )
    {
        stream.read( reinterpret_cast< char * >(x), static_cast<std::streamsize>(sizeof(T) * num) );
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
                      const T * x, size_t num )
    {
        for( size_t i = 0; i < num; ++i )
            write_field( stream, bo, x[i] );
    }

} // namespace vigra

#endif // VIGRA_BYTEORDER_HXX
