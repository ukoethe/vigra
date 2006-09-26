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
 
 
#ifndef VIGRA_STDIMAGEFUNCTIONS_HXX
#define VIGRA_STDIMAGEFUNCTIONS_HXX

/** \page PointOperators Point Operators 

    <DL>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref InitAlgo
        <DD><em>init images or image borders </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref InspectAlgo
        <DD> <em>Apply read-only functor to every pixel</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref InspectFunctor
        <DD><em>Functors which report image statistics</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref CopyAlgo
        <DD> <em>Copy images or regions</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref TransformAlgo
        <DD><em>apply functor to calculate a pixelwise transformation of one image</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref TransformFunctor
        <DD> <em>frequently used unary transformation functors</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref CombineAlgo
        <DD><em>apply functor to calculate a pixelwise transformation from several image</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref CombineFunctor
        <DD> <em>frequently used binary transformations functors</em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref MultiPointoperators
        <DD> <em>Point operators on multi-dimensional arrays</em>
    </DL>
    
    <b>\#include</b> "<a href="stdimagefunctions_8hxx-source.html">vigra/stdimagefunctions.hxx</a>"<br>
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
