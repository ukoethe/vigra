/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_BASICS_HXX
#define VIGRA_BASICS_HXX

#include "config.hxx"
#include "error.hxx"
#include "metaprogramming.hxx"
#include "tuple.hxx"
#include "diff2d.hxx"
#include "mathutil.hxx"
#include <string>
#include <sstream>

namespace vigra {

/** Convert a value to a string. Available for integral and floating point types
    and void *.
*/
doxygen_overloaded_function(template <class T> std::string asString(T t))

#define VIGRA_AS_STRING(T) \
inline std::string asString(T t) \
{ \
    std::stringstream s; \
    s << t; \
    return s.str(); \
}

VIGRA_AS_STRING(bool)
VIGRA_AS_STRING(signed char)
VIGRA_AS_STRING(unsigned char)
VIGRA_AS_STRING(signed short)
VIGRA_AS_STRING(unsigned short)
VIGRA_AS_STRING(signed long)
VIGRA_AS_STRING(unsigned long)
VIGRA_AS_STRING(signed long long)
VIGRA_AS_STRING(unsigned long long)
VIGRA_AS_STRING(signed int)
VIGRA_AS_STRING(unsigned int)
VIGRA_AS_STRING(float)
VIGRA_AS_STRING(double)
VIGRA_AS_STRING(long double)
VIGRA_AS_STRING(void *)

#undef VIGRA_AS_STRING

template <class T>
std::string & operator<<(std::string & s, T const & t)
{
    std::stringstream ss;
    ss << t; 
    return s += ss.str();
}

} // namespace vigra

/*! \page Utilities Utilities
    Basic helper functionality needed throughout.

    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref vigra::ArrayVector
         <BR>&nbsp;&nbsp;&nbsp;<em>replacement for std::vector (always uses consecutive memory)</em>
    <LI> \ref vigra::BucketQueue and \ref vigra::MappedBucketQueue
         <BR>&nbsp;&nbsp;&nbsp;<em>efficient priority queues for integer priorities</em>
    <LI> \ref RangesAndPoints
         <BR>&nbsp;&nbsp;&nbsp;<em>2-D and N-D positions, extents, and boxes</em>
    <LI> \ref PixelNeighborhood
         <BR>&nbsp;&nbsp;&nbsp;<em>4- and 8-neighborhood definitions and circulators</em>
    <LI> \ref VoxelNeighborhood
         <BR>&nbsp;&nbsp;&nbsp;<em>6- and 26-neighborhood definitions and circulators</em>
    <LI> \ref vigra::IteratorAdaptor
         <BR>&nbsp;&nbsp;&nbsp;<em>Quickly create STL-compatible 1D iterator adaptors</em>
    <LI> \ref TupleTypes
         <BR>&nbsp;&nbsp;&nbsp;<em>pair, triple, tuple4, tuple5</em>
    <LI> \ref MathConstants
         <BR>&nbsp;&nbsp;&nbsp;<em>M_PI, M_SQRT2</em>
    <LI> \ref TimingMacros
         <BR>&nbsp;&nbsp;&nbsp;<em>Macros for taking execution speed measurements</em>
    </UL>
*/

#endif // VIGRA_BASICS_HXX
