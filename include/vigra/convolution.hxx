/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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

#include <functional>
#include "stdconvolution.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/nonlineardiffusion.hxx"
#include "vigra/combineimages.hxx"

/** \page Convolution Functions to Convolve Images and Signals

    1D and 2D filters, including separable and recursive convolution, and non-linear diffusion

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"<br>
    Namespace: vigra

    <DL>
    <DT>
        <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
        \ref CommonConvolutionFilters
        <DD><em>Short-hands for the most common convolution filters</em>
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

namespace vigra {



/********************************************************/
/*                                                      */
/*             Common convolution filters               */
/*                                                      */
/********************************************************/

/** \addtogroup CommonConvolutionFilters Common Filters

    These functions calculate common filters by appropriate sequences of calls 
    to \link SeparableConvolution#separableConvolveX separableConvolveX\endlink()
    and \link SeparableConvolution#separableConvolveY separableConvolveY\endlink().
*/
//@{

/********************************************************/
/*                                                      */
/*                    convolveImage                     */
/*                                                      */
/********************************************************/

/** \brief Apply two separable filters successively, the first in x-direction, 
           the second in y-direction.

    This function is a shorthand for the concatenation of a call to
    \link SeparableConvolution#separableConvolveX separableConvolveX\endlink()
    and \link SeparableConvolution#separableConvolveY separableConvolveY\endlink() 
    with the given kernels.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class T>
        void convolveImage(SrcIterator supperleft,
                           SrcIterator slowerright, SrcAccessor sa,
                           DestIterator dupperleft, DestAccessor da,
                           Kernel1D<T> const & kx, Kernel1D<T> const & ky);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class T>
        inline void
        convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      pair<DestIterator, DestAccessor> dest,
                      Kernel1D<T> const & kx, Kernel1D<T> const & ky);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // implement sobel filter in x-direction
    Kernel1D<double> kx, ky;
    kx.initSymmetricGradient();
    ky.initBinomial(1);
    
    vigra::convolveImage(srcImageRange(src), destImage(dest), kx, ky);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void convolveImage(SrcIterator supperleft,
                   SrcIterator slowerright, SrcAccessor sa,
                   DestIterator dupperleft, DestAccessor da,
                   Kernel1D<T> const & kx, Kernel1D<T> const & ky)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(kx));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleft, da), kernel1d(ky));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
inline void
convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              Kernel1D<T> const & kx, Kernel1D<T> const & ky)
{
    convolveImage(src.first, src.second, src.third,
                  dest.first, dest.second, kx, ky);
}

/********************************************************/
/*                                                      */
/*                    simpleSharpening                  */
/*                                                      */
/********************************************************/

/** \brief Perform simple sharpening function.

    This function use \link StandardConvolution#convolveImage convolveImage\endlink( ) with following filter:
    
    \code
    -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0,
    -sharpening_factor/8.0,   1.0+sharpening_factor*0.75,  -sharpening_factor/8.0,
    -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0;    
    \endcode
    
    and use <TT>BORDER_TREATMENT_REFLECT</TT> as border treatment mode.

    <b> Preconditions:</b>
    \code  
    1. sharpening_factor >= 0
    2. scale >= 0
    \endcode

    <b> Declarations:</b>

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor,
                class DestIterator, class DestAccessor>
      void simpleSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                            DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor)

    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor, 
	        class DestIterator, class DestAccessor>
      inline
      void simpleSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
    	    		    pair<DestIterator, DestAccessor> dest, double sharpening_factor)
      {
	  simpleSharpening(src.first, src.second, src.third,
		  	   dest.first, dest.second, sharpening_factor);
      }

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 0.1
    vigra::simpleSharpening(srcImageRange(src), destImage(dest), 0.1);

    \endcode

*/    
template <class SrcIterator, class SrcAccessor,
	  class DestIterator, class DestAccessor>
void simpleSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                    DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor)
{

    vigra_precondition(sharpening_factor >= 0.0,
		       "simpleSharpening(): amount of sharpening must be >= 0.");

    Kernel2D<double> kernel;

    kernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0,
                                                        -sharpening_factor/8.0,   1.0+sharpening_factor*0.75,  -sharpening_factor/8.0,
                                                        -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0;

    convolveImage(src_ul, src_lr, src_acc, dest_ul, dest_acc, 
		  kernel.center(), kernel.accessor(), 
		  kernel.upperLeft(), kernel.lowerRight() , BORDER_TREATMENT_REFLECT );
}

template <class SrcIterator, class SrcAccessor, 
	  class DestIterator, class DestAccessor>
inline
void simpleSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
		    pair<DestIterator, DestAccessor> dest, double sharpening_factor)
{
    simpleSharpening(src.first, src.second, src.third,
		     dest.first, dest.second, sharpening_factor);
}


/********************************************************/
/*                                                      */
/*                    gaussianSharpening                */
/*                                                      */
/********************************************************/

/** \brief Perform sharpening function with gaussian filter.


    This function use the 
    \link vigra::gaussianSmoothing gaussianSmoothing \endlink()
    at first and scale the source image 
    (\code src \endcode) with the \code scale \endcode
    factor in an temporary image (\code tmp \endcode). At second the new 
    pixel in the destination image will be with following
    formel calculate:
    \code
    dest = (1 + sharpening_factor)*src - sharpening_factor*tmp
    \endcode

    <b> Preconditions:</b>
    \code  
    1. sharpening_factor >= 0
    2. scale >= 0
    \endcode

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor,
	        class DestIterator, class DestAccessor>
      void gaussianSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
	  		      DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor, 
			      double scale)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor,
	        class DestIterator, class DestAccessor>
      void gaussianSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
	 		      pair<DestIterator, DestAccessor> dest, double sharpening_factor, 
			      double scale)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 3.0
    // smoothing with scale = 0.5
    vigra::gaussianSmoothing(srcImageRange(src), destImage(dest), 3.0, 0.5);

    \endcode

*/    
template <class SrcIterator, class SrcAccessor,
	  class DestIterator, class DestAccessor>
void gaussianSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
			DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor, 
			double scale)
{
    vigra_precondition(sharpening_factor >= 0.0,
		       "gaussianSharpening(): amount of sharpening must be >= 0");
    vigra_precondition(scale >= 0.0,
		       "gaussianSharpening(): scale parameter should be >= 0.");

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote ValueType;

    BasicImage<ValueType> tmp(src_lr - src_ul);
    typename BasicImage<ValueType>::Accessor tmp_acc = tmp.accessor();

    gaussianSmoothing(src_ul, src_lr, src_acc, tmp.upperLeft(), tmp_acc, scale);

    SrcIterator i_src = src_ul;
    DestIterator i_dest = dest_ul;
    typename BasicImage<ValueType>::traverser tmp_ul = tmp.upperLeft();
    typename BasicImage<ValueType>::traverser i_tmp = tmp_ul;

    for(; i_src.y != src_lr.y ; i_src.y++, i_dest.y++, i_tmp.y++ )
    {
	for (;i_src.x != src_lr.x ; i_src.x++, i_dest.x++, i_tmp.x++ )
        {
	    dest_acc.set((1.0 + sharpening_factor)*src_acc(i_src) - sharpening_factor*tmp_acc(i_tmp), i_dest);
	}
	i_src.x = src_ul.x;
	i_dest.x = dest_ul.x;
	i_tmp.x = tmp_ul.x;
    }
}

template <class SrcIterator, class SrcAccessor,
	  class DestIterator, class DestAccessor>
void gaussianSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
			pair<DestIterator, DestAccessor> dest, double sharpening_factor, 
			double scale)
{
    gaussianSharpening(src.first, src.second, src.third,
		       dest.first, dest.second,
		       sharpening_factor, scale);
}



/********************************************************/
/*                                                      */
/*                    gaussianSmoothing                 */
/*                                                      */
/********************************************************/

/** \brief Perform isotropic Gaussian convolution.

    This function is a shorthand for the concatenation of a call to
    \link SeparableConvolution#separableConvolveX separableConvolveX\endlink()
    and \link SeparableConvolution#separableConvolveY separableConvolveY\endlink() with a
    Gaussian kernel of the given scale. The function uses 
    <TT>BORDER_TREATMENT_REFLECT</TT>.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void gaussianSmoothing(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIterator dupperleft, DestAccessor da,
                                double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        inline void
        gaussianSmoothing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // smooth with scale = 3.0
    vigra::gaussianSmoothing(srcImageRange(src), destImage(dest), 3.0);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void gaussianSmoothing(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft);

    Kernel1D<double> smooth;
    smooth.initGaussian(scale);
    smooth.setBorderTreatment(BORDER_TREATMENT_REFLECT);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleft, da), kernel1d(smooth));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  double scale)
{
    gaussianSmoothing(src.first, src.second, src.third,
                 dest.first, dest.second, scale);
}

/********************************************************/
/*                                                      */
/*                     gaussianGradient                 */
/*                                                      */
/********************************************************/

/** \brief Calculate the gradient vector by means of a 1st derivatives of
    Gaussian filter.

    This function is a shorthand for the concatenation of a call to
    \link SeparableConvolution#separableConvolveX separableConvolveX\endlink()
    and \link SeparableConvolution#separableConvolveY separableConvolveY\endlink() with the
    appropriate kernels at the given scale. Note that this function produces
    <i>two</i> result images.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorY, class DestAccessorY>
        void gaussianGradient(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIteratorX dupperleftx, DestAccessorX dax,
                                DestIteratorY dupperlefty, DestAccessorY day,
                                double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorY, class DestAccessorY>
        inline void
        gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorX, DestAccessorX> destx,
                          pair<DestIteratorY, DestAccessorY> desty,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), gradx(w,h), grady(w,h);
    ...

    // calculate gradient vector at scale = 3.0
    vigra::gaussianGradient(srcImageRange(src),
                             destImage(gradx), destImage(grady), 3.0);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorY, class DestAccessorY>
void gaussianGradient(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft);

    Kernel1D<double> smooth, grad;
    smooth.initGaussian(scale);
    grad.initGaussianDerivative(scale, 1);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(grad));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftx, dax), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperlefty, day), kernel1d(grad));
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorY, class DestAccessorY>
inline void
gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorX, DestAccessorX> destx,
                  pair<DestIteratorY, DestAccessorY> desty,
                  double scale)
{
    gaussianGradient(src.first, src.second, src.third,
                 destx.first, destx.second, desty.first, desty.second, scale);
}

/********************************************************/
/*                                                      */
/*                 laplacianOfGaussian                  */
/*                                                      */
/********************************************************/

/** \brief Filter image with the Laplacian of Gaussian operator
    at the given scale.

    This function calls \link SeparableConvolution#separableConvolveX separableConvolveX\endlink() and
    \link SeparableConvolution#separableConvolveY separableConvolveY\endlink() with the appropriate 2nd derivative
    of Gaussian kernels in x- and y-direction and then sums the results
    to get the Laplacian.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void laplacianOfGaussian(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIterator dupperleft, DestAccessor da,
                                double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        inline void
        laplacianOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // calculate Laplacian of Gaussian at scale = 3.0
    vigra::laplacianOfGaussian(srcImageRange(src), destImage(dest), 3.0);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void laplacianOfGaussian(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft),
                        tmpx(slowerright - supperleft),
                        tmpy(slowerright - supperleft);

    Kernel1D<double> smooth, deriv;
    smooth.initGaussian(scale);
    deriv.initGaussianDerivative(scale, 2);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv));
    separableConvolveY(srcImageRange(tmp),
                       destImage(tmpx), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destImage(tmpy), kernel1d(deriv));
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpy),
                       destIter(dupperleft, da), std::plus<TmpType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
laplacianOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  double scale)
{
    laplacianOfGaussian(src.first, src.second, src.third,
                 dest.first, dest.second, scale);
}

/********************************************************/
/*                                                      */
/*               hessianMatrixOfGaussian                */
/*                                                      */
/********************************************************/

/** \brief Filter image with the 2nd derivatives of the Gaussian
    at the given scale to get the Hessian matrix.

    The Hessian matrix is a symmetric matrix defined as:

	    \f[
                \mbox{\rm Hessian}(I) = \left(
                \begin{array}{cc}
                G_{xx} \ast I & G_{xy} \ast I \\
                G_{xy} \ast I & G_{yy} \ast I
                \end{array} \right)
	    \f]

    where \f$G_{xx}, G_{xy}, G_{yy}\f$ denote 2nd derivatives of Gaussians
    at the given scale, and
    \f$\ast\f$ is the convolution symbol. This function calls
    \link SeparableConvolution#separableConvolveX separableConvolveX\endlink() and
    \link SeparableConvolution#separableConvolveY separableConvolveY\endlink()
    with the appropriate 2nd derivative
    of Gaussian kernels and puts the results in
    the three destination images. The first destination image will
    contain the second derivative in x-direction, the second one the mixed
    derivative, and the third one holds the derivative in y-direction.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void hessianMatrixOfGaussian(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIteratorX dupperleftx, DestAccessorX dax,
                                DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                                DestIteratorY dupperlefty, DestAccessorY day,
                                double scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        inline void
        hessianMatrixOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorX, DestAccessorX> destx,
                          pair<DestIteratorXY, DestAccessorXY> destxy,
                          pair<DestIteratorY, DestAccessorY> desty,
                          double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), hxx(w,h), hxy(w,h), hyy(w,h);
    ...

    // calculate Hessian of Gaussian at scale = 3.0
    vigra::hessianMatrixOfGaussian(srcImageRange(src),
        destImage(hxx), destImage(hxy), destImage(hyy), 3.0);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
void hessianMatrixOfGaussian(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft);

    Kernel1D<double> smooth, deriv1, deriv2;
    smooth.initGaussian(scale);
    deriv1.initGaussianDerivative(scale, 1);
    deriv2.initGaussianDerivative(scale, 2);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv2));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftx, dax), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperlefty, day), kernel1d(deriv2));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv1));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftxy, daxy), kernel1d(deriv1));
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
inline void
hessianMatrixOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorX, DestAccessorX> destx,
                  pair<DestIteratorXY, DestAccessorXY> destxy,
                  pair<DestIteratorY, DestAccessorY> desty,
                  double scale)
{
    hessianMatrixOfGaussian(src.first, src.second, src.third,
                 destx.first, destx.second,
                 destxy.first, destxy.second,
                 desty.first, desty.second,
                 scale);
}

/********************************************************/
/*                                                      */
/*                   structureTensor                    */
/*                                                      */
/********************************************************/

/** \brief Calculate the Structure Tensor for each pixel of
 and image, using Gaussian (derivative) filters.

    The Structure Tensor is is a smoothed version of the Euclidean product
    of the gradient vector with itself. I.e. it's a symmetric matrix defined as:

	    \f[
                \mbox{\rm StructurTensor}(I) = \left(
                \begin{array}{cc}
                G \ast (I_x I_x) & G \ast (I_x I_y) \\
                G \ast (I_x I_y) & G \ast (I_y I_y)
                \end{array} \right) = \left(
                \begin{array}{cc}
                A & C \\
                C & B
                \end{array} \right)
	    \f]

    where \f$G\f$ denotes Gaussian smoothing at the <i>outer scale</i>,
    \f$I_x, I_y\f$ are the gradient components taken at the <i>inner scale</i>,
    \f$\ast\f$ is the convolution symbol, and \f$I_x I_x\f$ etc. are pixelwise
    products of the 1st derivative images. This function calls
    \link SeparableConvolution#separableConvolveX separableConvolveX\endlink()
    and \link SeparableConvolution#separableConvolveY separableConvolveY\endlink() with the
    appropriate Gaussian kernels and puts the results in
    the three destination images. The first destination image will
    contain \f$G \ast (I_x I_x)\f$, the second one \f$G \ast (I_x I_y)\f$, and the
    third one holds \f$G \ast (I_y I_y)\f$.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void structureTensor(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIteratorX dupperleftx, DestAccessorX dax,
                                DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                                DestIteratorY dupperlefty, DestAccessorY day,
                                double inner_scale, double outer_scale);
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        inline void
        structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorX, DestAccessorX> destx,
                          pair<DestIteratorXY, DestAccessorXY> destxy,
                          pair<DestIteratorY, DestAccessorY> desty,
                          double nner_scale, double outer_scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="convolution_8hxx-source.html">vigra/convolution.hxx</a>"


    \code
    vigra::FImage src(w,h), stxx(w,h), stxy(w,h), styy(w,h);
    ...

    // calculate Structure Tensor at inner scale = 1.0 and outer scale = 3.0
    vigra::structureTensor(srcImageRange(src),
        destImage(stxx), destImage(stxy), destImage(styy), 1.0, 3.0);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
void structureTensor(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double inner_scale, double outer_scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft),
                        tmpx(slowerright - supperleft),
                        tmpy(slowerright - supperleft);

    gaussianGradient(srcIterRange(supperleft, slowerright, sa),
           destImage(tmpx), destImage(tmpy), inner_scale);
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpx),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperleftx, dax), outer_scale);
    combineTwoImages(srcImageRange(tmpy), srcImage(tmpy),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperlefty, day), outer_scale);
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpy),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperleftxy, daxy), outer_scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
inline void
structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorX, DestAccessorX> destx,
                  pair<DestIteratorXY, DestAccessorXY> destxy,
                  pair<DestIteratorY, DestAccessorY> desty,
                  double inner_scale, double outer_scale)
{
    structureTensor(src.first, src.second, src.third,
                 destx.first, destx.second,
                 destxy.first, destxy.second,
                 desty.first, desty.second,
                 inner_scale, outer_scale);
}

//@}

} // namespace vigra

#endif // VIGRA_CONVOLUTION_HXX
