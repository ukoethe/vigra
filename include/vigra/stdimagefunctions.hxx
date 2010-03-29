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
 
 
#ifndef VIGRA_STDIMAGEFUNCTIONS_HXX
#define VIGRA_STDIMAGEFUNCTIONS_HXX

/** \page PointOperators Point Operators 

    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref InitAlgo
         <BR>&nbsp;&nbsp;&nbsp; <em>init images or image borders </em>
    <LI> \ref InspectAlgo
         <BR>&nbsp;&nbsp;&nbsp; <em>Apply read-only functor to every pixel</em>
    <LI> \ref InspectFunctor
         <BR>&nbsp;&nbsp;&nbsp; <em>Functors which report image statistics</em>
    <LI> \ref CopyAlgo
         <BR>&nbsp;&nbsp;&nbsp; <em>Copy images or regions</em>
    <LI> \ref TransformAlgo
         <BR>&nbsp;&nbsp;&nbsp; <em>apply functor to calculate a pixelwise transformation of one image</em>
    <LI> \ref TransformFunctor
         <BR>&nbsp;&nbsp;&nbsp; <em>frequently used unary transformation functors</em>
    <LI> \ref CombineAlgo
         <BR>&nbsp;&nbsp;&nbsp; <em>apply functor to calculate a pixelwise transformation from several image</em>
    <LI> \ref CombineFunctor
         <BR>&nbsp;&nbsp;&nbsp; <em>frequently used binary transformations functors</em>
    <LI> \ref MultiPointoperators
         <BR>&nbsp;&nbsp;&nbsp; <em>Point operators on multi-dimensional arrays</em>
    </UL>
    
    <b>\#include</b> \<<a href="stdimagefunctions_8hxx-source.html">vigra/stdimagefunctions.hxx</a>\><br>
    Namespace: vigra
        
    see also: \ref FunctorExpressions "Automatic Functor Creation"
*/

#include "initimage.hxx"
#include "inspectimage.hxx"
#include "copyimage.hxx"
#include "transformimage.hxx"
#include "combineimages.hxx"
#include "resizeimage.hxx"

#endif // VIGRA_STDIMAGEFUNCTIONS_HXX
