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
 
 
#ifndef VIGRA_CORNERDETECTION_HXX
#define VIGRA_CORNERDETECTION_HXX

#include <vigra/utilities.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/convolution.hxx>

namespace vigra {

template <class SrcType, class DestType>
struct InternalCornerResponseFunctor
{
    DestType operator()(SrcType const & a1, 
                        SrcType const & a2, SrcType const & a3) const
    {
	return (a1*a2 - a3*a3) - 0.04 * (a1 + a2) * (a1 + a2);
    }

    SrcType operator()(SrcType const & a1, SrcType const & a2) const
    {
	return (a1*a2);
    }

    SrcType operator()(SrcType const & a) const
    {
	return (a*a);
    }
};

/** \addtogroup CornerDetection Corner Detection
    Measure the 'cornerness' at each pixel
*/
//@{ 
                                    
/********************************************************/
/*                                                      */
/*                 cornerResponseFunction               */
/*                                                      */
/********************************************************/

/** \brief Find corners in an image.

    This algorithm implements the so called 'corner response function'
    to measure the 'cornerness' of each pixel in the image, according to
    [C.G. Harris and M.J. Stevens: <em> "A Combined Corner and Edge Detector"</em>,
    Proc. of 4th Alvey Vision Conference, 1988]. We have found this to be a
    very robust corner detector, although it moves the corners somewhat into one
    region, depending on the scale.
    
    The algorithm proceeds as follows: At the given <TT>scale</TT>, it calculates the 
    structure tensor (which is the smoothed matrix of gradient products) at each 
    pixel like this:
    
    \f[ 
        \left( 
	\begin{array}{cc}
	G^s * (g_x^s g_x^s) & G^s * (g_x^s g_y^s) \\
	G^s * (g_x^s g_y^s) & G^s * (g_y^s g_y^s) \\
	\end{array} 
        \right) = 
        \left( 
	\begin{array}{cc}
	A & C \\
	C & B \\
	\end{array} 
        \right)
    \f]
    
    where <math>G<sup>s</sup> *</math> denotes convolution with a Gaussian filter at 
    scale s, <math>g<sub>x</sub><sup>s</sup></math> and 
    <math>g<sub>y</sub><sup>s</sup></math> are first Gaussian 
    derivatives in x and y directions at scale s, and products mean pixelwise
    multiplication. Then the corner response is calculated pixelwise as 
    
    \f[
        CornerResponse = A B - C^2 - 0.04 (A + B)^2
    \f]
    
    The local maxima of the corner response denote the corners in the gray level 
    function.
    
    The source value type must be an linaer algebra, i.e. addition, subtraction, and
    multiplication with itself, multiplication with doubles and 
    \ref NumericTraits "NumericTraits" must 
    be defined. (Sorry, only the built-in types are currently supported.)
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        void
        cornerResponseFunction(SrcIterator sul, SrcIterator slr, SrcAccessor as,
			       DestIterator dul, DestAccessor ad,
			       double scale)
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        inline 
        void cornerResponseFunction(
	           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	           pair<DestIterator, DestAccessor> dest,
	           double scale)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="cornerdetection_8hxx-source.html">vigra/cornerdetection.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), corners(w,h);
    vigra::FImage corner_response(w,h);
    
    // empty corner image
    corners = 0;
    ...
    
    // find corner response at scale 1.0
    vigra::cornerResponseFunction(srcImageRange(src), destImage(corner_response), 
                           1.0);
    
    // find local maxima of corner response, mark with 1
    vigra::localMaxima(srcImageRange(corner_response), destImage(corners));
    
    // threshold corner response to keep only strong corners (above 400.0)
    transformImage(srcImageRange(corner_response), destImage(corner_response),
          vigra::Threshold<double, double>(
               400.0, std::numeric_limits<double>::max(), 0.0, 1.0)); 

    // combine thresholding and local maxima
    vigra::combineTwoImages(srcImageRange(corners), srcImage(corner_response),
                     destImage(corners), std::multiplies<float>());
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    double d;
    
    u = u + u
    u = u - u
    u = u * u
    u = d * u
    
    dest_accessor.set(u, dest_upperleft);
    \endcode
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
cornerResponseFunction(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                       DestIterator dul, DestAccessor ad,
		       double scale)
{
    vigra_precondition(scale > 0.0,
                 "cornerResponseFunction(): Scale must be > 0");
		 
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    if(w <= 0 || h <= 0) return;
    
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
	
    typedef BasicImage<TmpType> TmpImage;
    
    TmpImage tmp(w,h);
    TmpImage gx(w,h);
    TmpImage gy(w,h);
    TmpImage gxy(w,h);
    
    Kernel1D<double> smooth, grad;
    smooth.initGaussian(scale);
    grad.initGaussianDerivative(scale, 1);
    
    separableConvolveX(srcIterRange(sul, slr, as), destImage(tmp), kernel1d(grad));
    separableConvolveY(srcImageRange(tmp), destImage(gx), kernel1d(smooth));
    
    separableConvolveY(srcIterRange(sul, slr, as), destImage(tmp), kernel1d(grad));
    separableConvolveX(srcImageRange(tmp), destImage(gy), kernel1d(smooth));
    
    InternalCornerResponseFunctor<TmpType, 
                typename DestAccessor::value_type > CRF;

    combineTwoImages(srcImageRange(gx), srcImage(gy), destImage(gxy), CRF);
		  
    transformImage(srcImageRange(gx), destImage(gx), CRF);
    transformImage(srcImageRange(gy), destImage(gy), CRF);
    
    separableConvolveX(srcImageRange(gx), destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp), destImage(gx), kernel1d(smooth));

    separableConvolveX(srcImageRange(gy), destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp), destImage(gy), kernel1d(smooth));

    separableConvolveX(srcImageRange(gxy), destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp), destImage(gxy), kernel1d(smooth));

    combineThreeImages(srcImageRange(gx), srcImage(gy), srcImage(gxy), 
                       destIter(dul, ad), CRF );
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void cornerResponseFunction(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	   pair<DestIterator, DestAccessor> dest,
	   double scale)
{
    cornerResponseFunction(src.first, src.second, src.third,
                            dest.first, dest.second,
			    scale);
}

//@}

} // namespace vigra

#endif // VIGRA_CORNERDETECTION_HXX
