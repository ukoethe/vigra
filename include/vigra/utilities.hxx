/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#include "vigra/config.hxx"
#include "vigra/error.hxx"
#include "vigra/metaprogramming.hxx"
#include "vigra/tuple.hxx"
#include "vigra/diff2d.hxx"
#include "vigra/mathutil.hxx"

/*! \page Utilities Utilities
    Basic helper functionality needed throughout.

    <DL>
    <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref vigra::ArrayVector
     <DD><em>replacement for std::vector</em>
    <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref RangesAndPoints
     <DD><em>2-dimensioanl positions, extents, amd rectangles</em>
    <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref PixelNeighborhood
     <DD><em>4- and 8-neighborhood definitions and circulators</em>
    <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref vigra::IteratorAdaptor
     <DD><em>Quickly create STL-compatible 1D iterator adaptors</em>
     <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref TupleTypes
     <DD><em>pair, triple, tuple4, tuple5</em>
      <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
     \ref MathConstants
     <DD><em>M_PI, M_SQRT2</em>
    </DL>
*/

#endif // VIGRA_BASICS_HXX
