/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
    </DL>
    
    <b>\#include</b> "<a href="stdimagefunctions_8hxx-source.html">vigra/stdimagefunctions.hxx</a>"<br>
    Namespace: vigra
        
    see also: \ref FunctorExpressions "Automatic Functor Creation"
*/

#include "vigra/initimage.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/resizeimage.hxx"

#endif // VIGRA_STDIMAGEFUNCTIONS_HXX
