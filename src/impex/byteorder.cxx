/************************************************************************/
/*                                                                      */
/*      Copyright 2002-2004 by Ullrich Koethe and Gunnar Kedenburg      */
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

#include <algorithm>
#include <stdexcept>
#include "byteorder.hxx"
#include "vigra/sized_int.hxx"

namespace vigra
{
    const byteorder::host byteorder::m_host;

    byteorder::host::host()
    {
        // byteorder check: if the first byte is the least significant,
        // we have little endian byteorder.
        UIntBiggest testint = 0x01;
        UInt8 * testchar = reinterpret_cast< UInt8 * >(&testint);
        if ( testchar[0] == 0x01 )
            m_string = "little endian";
        else
            m_string = "big endian";
    }

    const std::string & byteorder::host::get() const
    {
        return m_string;
    }

    const std::string & byteorder::get_host_byteorder() const
    {
        return m_host.get();
    }

    byteorder::byteorder()
    {
        set( m_host.get() );
    }

    byteorder::byteorder( const std::string & s )
    {
        set(s);
    }

    const std::string & byteorder::get() const
    {
        return m_string;
    }

    void byteorder::set( const std::string & s )
    {
        m_string = s;
        native = get_host_byteorder() == s;
    }
}
