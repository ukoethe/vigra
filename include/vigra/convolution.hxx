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
 

#ifndef VIGRA_CONVOLUTION_HXX
#define VIGRA_CONVOLUTION_HXX

/** \page Convolution Functions to Convolve Images and Signals 

    1D and 2D filters, including separable and recursive convolution, and non-linear diffusion
    
    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"<br>
    Namespace: vigra

    <DL>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref BorderTreatmentMode
        <DD><em>Choose between different border treatment modes </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref StandardConvolution
        <DD><em>2D non-separable convolution, with and without ROI mask </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref vigra::Kernel2D
        <DD><em>Generic 2-dimensional convolution kernel </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref SeparableConvolution
        <DD> <em>1D convolution and separable filters in 2 dimensions </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref vigra::Kernel1D
        <DD> <em>Generic 1-dimensional convolution kernel </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref RecursiveConvolution
        <DD> <em>Recursive implementation of the exponential filter and its derivatives </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref NonLinearDiffusion
        <DD> <em>Edge-preserving smoothing </em>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
        \ref KernelArgumentObjectFactories
        <DD> <em>Factory functions to create argument objects to simplify passing kernels</em>
    </DL>
*/

#include "vigra/stdconvolution.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/nonlineardiffusion.hxx"

/** \page KernelArgumentObjectFactories Kernel Argument Object Factories

    These factory functions allow to create argument objects for 1D
    and 2D convolution kernel analogously to 
    \ref ArgumentObjectFactories for images. 

    \section Kernel1dFactory kernel1d()
        
	Pass a \ref vigra::Kernel1D to a 1D or separable convolution algorithm.
        
        These factories can be used to create argument objects when we 
	are given instances or subclasses of \ref vigra::Kernel1D
	(analogous to the \ref ArgumentObjectFactories for images).
	These factory functions access <TT>kernel.center()</TT>,
        <TT>kernel.left()</TT>, <TT>kernel.right()</TT>, <TT>kernel.accessor()</TT>,
	and  <TT>kernel.borderTreatment()</TT> to obtain the necessary 
	information. The following factory functions are provided:
	
	<table>
        <tr><td>
            \htmlonly
            <th bgcolor="#f0e0c0" colspan=2 align=left>
            \endhtmlonly
            <TT>\ref vigra::Kernel1D "vigra::Kernel1D<SomeType>" kernel;</TT>
            \htmlonly
            </th>
            \endhtmlonly
        </td></tr>
        <tr><td>
	<TT>kernel1d(kernel)</TT>
        </td><td>
	    create argument object from information provided by 
	    kernel
	    
        </td></tr>
        <tr><td>
	<TT>kernel1d(kernel, vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
	    create argument object from information provided by 
	    kernel, but use given border treatment mode
	    
        </td></tr>
        <tr><td>
	<TT>kernel1d(kerneliterator, kernelaccessor,</TT><br>
	<TT>                kernelleft, kernelright,</TT><br>
	<TT>                vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
	    create argument object from explicitly given iterator
	    (pointing to the center of th kernel), accessor,
	    left and right boundaries, and border treatment mode

	</table>
	
	For usage examples see 
	\ref SeparableConvolution "one-dimensional and separable convolution functions".

    \section Kernel2dFactory kernel2d()
        
	Pass a \ref vigra::Kernel2D to a 2D (non-separable) convolution algorithm.
                
	These factories can be used to create argument objects when we 
	are given instances or subclasses of \ref vigra::Kernel2D
	(analogous to the \ref ArgumentObjectFactories for images).
	These factory functions access <TT>kernel.center()</TT>,
        <TT>kernel.upperLeft()</TT>, <TT>kernel.lowerRight()</TT>, <TT>kernel.accessor()</TT>,
	and  <TT>kernel.borderTreatment()</TT> to obtain the necessary 
	information. The following factory functions are provided:
	
	<table>
        <tr><td>
            \htmlonly
            <th bgcolor="#f0e0c0" colspan=2 align=left>
            \endhtmlonly
            <TT>\ref vigra::Kernel2D "vigra::Kernel2D<SomeType>" kernel;</TT>
            \htmlonly
            </th>
            \endhtmlonly
        </td></tr>
        <tr><td>
	<TT>kernel2d(kernel)</TT>
        </td><td>
	    create argument object from information provided by 
	    kernel
	    
        </td></tr>
        <tr><td>
	<TT>kernel2d(kernel, vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
	    create argument object from information provided by 
	    kernel, but use given border treatment mode
	    
        </td></tr>
        <tr><td>
	<TT>kernel2d(kerneliterator, kernelaccessor,</TT>
	<TT>                upperleft, lowerright,</TT>
	<TT>                vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
	    create argument object from explicitly given iterator
	    (pointing to the center of th kernel), accessor,
	    upper left and lower right corners, and border treatment mode

	</table>
	
	For usage examples see \ref StandardConvolution "two-dimensional convolution functions".
    */


#endif // VIGRA_CONVOLUTION_HXX
