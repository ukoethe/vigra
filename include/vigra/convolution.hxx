/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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

/** @heading Functions to Convolve Images and Signals 

    Include-File:
    \URL[vigra/convolution.hxx]{../include/vigra/convolution.hxx}\\
    Namespace: vigra
    
    @memo including standard and separable convolution, and recursive filters.
*/
//@{

//@Include: bordertreatment.hxx stdconvolution.hxx 
//@Include: separableconvolution.hxx recursiveconvolution.hxx

#include "vigra/stdconvolution.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/recursiveconvolution.hxx"

/** @heading Kernel Argument Object Factories

    These factory functions allow to create argument objects for 1D
    and 2D convolution kernel analogously to 
    \Ref{Argument Object Factories} for images. 
    
    @memo Factory functions to create argument objects to simplify passing kernels
*/
//@{
    /** @heading kernel1d()
        
	These factories can be used to create argument objects when we 
	are given instances or subclasses of \Ref{Kernel1D}
	(analogous to the \Ref{Argument Object Factories} for images).
	These factory functions access #kernel.center()#,
        #kernel.left()#, #kernel.right()#, #kernel.accessor()#,
	and  #kernel.borderTreatment()# to obtain the necessary 
	information. The following factory functions are provided:
	
	\begin{tabular}{ll}
	#vigra::Kernel1D<SomeType> kernel;# & 
	
	    \\
	    
	#kernel1d(kernel)# &
	    create argument object from information provided by 
	    kernel
	    
	    \\
	    
	#kernel1d(kernel, vigra::BORDER_TREATMENT_CLIP)# &
	    create argument object from information provided by 
	    kernel, but use given border treatment mode
	    
	    \\
	    
	#kernel1d(kerneliterator, kernelaccessor,#
	#                kernelleft, kernelright,#
	#                vigra::BORDER_TREATMENT_CLIP)# &
	    create argument object from explicitly given iterator
	    (pointing to the center of th kernel), accessor,
	    left and right boundaries, and border treatment mode

	    \\
	\end{tabular}
	
	For usage examples see 
	\Ref{One-dimensional and separable convolution functions}.
	
	@memo Argument Object Factories for \Ref{Kernel1D}
    */
    /** @heading kernel2d()
        
	These factories can be used to create argument objects when we 
	are given instances or subclasses of \Ref{Kernel2D}
	(analogous to the \Ref{Argument Object Factories} for images).
	These factory functions access #kernel.center()#,
        #kernel.upperLeft()#, #kernel.lowerRight()#, #kernel.accessor()#,
	and  #kernel.borderTreatment()# to obtain the necessary 
	information. The following factory functions are provided:
	
	\begin{tabular}{ll}
	#vigra::Kernel2D<SomeType> kernel;# & 
	
	    \\
	    
	#kernel2d(kernel)# &
	    create argument object from information provided by 
	    kernel
	    
	    \\
	    
	#kernel2d(kernel, vigra::BORDER_TREATMENT_CLIP)# &
	    create argument object from information provided by 
	    kernel, but use given border treatment mode
	    
	    \\
	    
	#kernel2d(kerneliterator, kernelaccessor,#
	#                upperleft, lowerright,#
	#                vigra::BORDER_TREATMENT_CLIP)# &
	    create argument object from explicitly given iterator
	    (pointing to the center of th kernel), accessor,
	    upper left and lower right corners, and border treatment mode

	    \\
	\end{tabular}
	
	For usage examples see \Ref{Two-dimensional convolution functions}.
	
	@memo Argument Object Factories for \Ref{Kernel2D}
    */
//@}

//@}



#endif // VIGRA_CONVOLUTION_HXX
