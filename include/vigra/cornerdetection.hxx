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
 
 
#ifndef VIGRA_CORNERDETECTION_HXX
#define VIGRA_CORNERDETECTION_HXX

#include <vigra/utilities.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/convolution.hxx>

namespace vigra {

template <class SrcType>
struct CornerResponseFunctor
{
    typedef typename NumericTraits<SrcType>::RealPromote argument_type;
    typedef argument_type result_type;
    
    result_type operator()(argument_type a1, 
                        argument_type a2, argument_type a3) const
    {
	return (a1*a2 - a3*a3) - 0.04 * (a1 + a2) * (a1 + a2);
    }
};

template <class SrcType>
struct FoerstnerCornerFunctor
{
    typedef typename NumericTraits<SrcType>::RealPromote argument_type;
    typedef argument_type result_type;
    
    result_type operator()(argument_type a1, 
                           argument_type a2, argument_type a3) const
    {
	return (a1*a2 - a3*a3) / (a1 + a2);
    }
};

template <class SrcType>
struct RohrCornerFunctor
{
    typedef typename NumericTraits<SrcType>::RealPromote argument_type;
    typedef argument_type result_type;
    
    result_type operator()(argument_type a1, 
                        argument_type a2, argument_type a3) const
    {
	return (a1*a2 - a3*a3);
    }
};

template <class SrcType>
struct BeaudetCornerFunctor
{
    typedef typename NumericTraits<SrcType>::RealPromote argument_type;
    typedef argument_type result_type;
    
    result_type operator()(argument_type a1, 
                        argument_type a2, argument_type a3) const
    {
	return (a3*a3 - a1*a2);
    }
};

/** \addtogroup CornerDetection Corner Detection
    Measure the 'cornerness' at each pixel.
    Note: The Kitchen-Rosenfeld detector is not implemented because of its
    inferior performance. The SUSAN detector is missing because it's patented.
*/
//@{ 
                                    
/********************************************************/
/*                                                      */
/*                 cornerResponseFunction               */
/*                                                      */
/********************************************************/

/** \brief Find corners in an image (1).

    This algorithm implements the so called 'corner response function'
    to measure the 'cornerness' of each pixel in the image, according to
    [C.G. Harris and M.J. Stevens: <em> "A Combined Corner and Edge Detector"</em>,
    Proc. of 4th Alvey Vision Conference, 1988]. Several studies have found this to be a
    very robust corner detector, although it moves the corners somewhat into one
    region, depending on the scale.
    
    The algorithm first determines the structure tensor at each pixel by calling
    \link CommonConvolutionFilters#structureTensor structureTensor\endlink(). 
    Then the entries of the structure tensor are combined as 
    
    \f[
        \mbox{\rm CornerResponse} = \mbox{\rm det(StructureTensor)} - 0.04 \mbox{\rm tr(StructureTensor)}^2
        = A B - C^2 - 0.04 (A + B)^2
    \f]
    
    The local maxima of the corner response denote the corners in the gray level 
    image.
    
    The source value type must be a linaer algebra, i.e. addition, subtraction, and
    multiplication with itself, multiplication with doubles and 
    \ref NumericTraits "NumericTraits" must 
    be defined. 
    
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
    
    TmpImage gx(w,h);
    TmpImage gy(w,h);
    TmpImage gxy(w,h);

    structureTensor(srcIterRange(sul, slr, as), 
                    destImage(gx), destImage(gxy), destImage(gy), 
                    scale, scale);
    CornerResponseFunctor<typename SrcAccessor::value_type > cf;
                    
    combineThreeImages(srcImageRange(gx), srcImage(gy), srcImage(gxy), 
                       destIter(dul, ad), cf );
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

/********************************************************/
/*                                                      */
/*               foerstnerCornerDetector                */
/*                                                      */
/********************************************************/

/** \brief Find corners in an image (2).

    This algorithm implements the so called 'Foerstner Corner Detector'
    to measure the 'cornerness' of each pixel in the image, according to
    [W. F&ouml;rstner: <em> "A feature based correspondence algorithms for image
    matching"</em>, Intl. Arch. Photogrammetry and Remote Sensing, vol. 24, pp 160-166, 
    1986]. It is also known as the "Plessey Detector" by Harris. However, it should not 
    be confused with the
    "\link CornerDetection#cornerResponseFunction Corner Repsonse Function\endlink ",
    another detector invented by Harris.
    
    The algorithm first determines the structure tensor at each pixel by calling
    \link CommonConvolutionFilters#structureTensor structureTensor\endlink(). 
    Then the entries of the structure tensor are combined as 
    
    \f[
        \mbox{\rm FoerstnerCornerStrength} = \frac{\mbox{\rm det(StructureTensor)}}{\mbox{\rm tr(StructureTensor)}} = 
        \frac{A B - C^2}{A + B}
    \f]
    
    The local maxima of the corner strength denote the corners in the gray level 
    image. Its performance is similar to the 
    \link CornerDetection#cornerResponseFunction cornerResponseFunction\endlink().
    
    The source value type must be a division algebra, i.e. addition, subtraction,
    multiplication, and division with itself, multiplication with doubles and 
    \ref NumericTraits "NumericTraits" must 
    be defined.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        void
        foerstnerCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
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
        void foerstnerCornerDetector(
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
    vigra::FImage foerstner_corner_strength(w,h);
    
    // empty corner image
    corners = 0;
    ...
    
    // find corner response at scale 1.0
    vigra::foerstnerCornerDetector(srcImageRange(src), destImage(foerstner_corner_strength), 
                                   1.0);
    
    // find local maxima of corner response, mark with 1
    vigra::localMaxima(srcImageRange(foerstner_corner_strength), destImage(corners));
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
    u = u / u
    u = d * u
    
    dest_accessor.set(u, dest_upperleft);
    \endcode
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
foerstnerCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                       DestIterator dul, DestAccessor ad,
		       double scale)
{
    vigra_precondition(scale > 0.0,
                 "foerstnerCornerDetector(): Scale must be > 0");
		 
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    if(w <= 0 || h <= 0) return;
    
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
	
    typedef BasicImage<TmpType> TmpImage;
    
    TmpImage gx(w,h);
    TmpImage gy(w,h);
    TmpImage gxy(w,h);

    structureTensor(srcIterRange(sul, slr, as), 
                    destImage(gx), destImage(gxy), destImage(gy), 
                    scale, scale);
    FoerstnerCornerFunctor<typename SrcAccessor::value_type > cf;
                    
    combineThreeImages(srcImageRange(gx), srcImage(gy), srcImage(gxy), 
                       destIter(dul, ad), cf );
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void foerstnerCornerDetector(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	   pair<DestIterator, DestAccessor> dest,
	   double scale)
{
    foerstnerCornerDetector(src.first, src.second, src.third,
                            dest.first, dest.second,
			    scale);
}

/********************************************************/
/*                                                      */
/*                   rohrCornerDetector                 */
/*                                                      */
/********************************************************/

/** \brief Find corners in an image (3).

    This algorithm implements yet another structure tensor-based corner detector, 
    according to [K. Rohr: <em> "Untersuchung von grauwertabh&auml;ngigen 
    Transformationen zur Ermittlung der optischen Flusses in Bildfolgen"</em>, 
    Diploma thesis, Inst. f&uuml;r Nachrichtensysteme, Univ. Karlsruhe, 1987]. 
    
    The algorithm first determines the structure tensor at each pixel by calling
    \link CommonConvolutionFilters#structureTensor structureTensor\endlink(). 
    Then the entries of the structure tensor are combined as 
    
    \f[
        \mbox{\rm RohrCornerStrength} = \mbox{\rm det(StructureTensor)} = A B - C^2
    \f]
    
    The local maxima of the corner strength denote the corners in the gray level 
    image. Its performance is similar to the 
    \link CornerDetection#cornerResponseFunction cornerResponseFunction\endlink().
    
    The source value type must be a linear algebra, i.e. addition, subtraction, and
    multiplication with itself, multiplication with doubles and 
    \ref NumericTraits "NumericTraits" must 
    be defined.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        void
        rohrCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
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
        void rohrCornerDetector(
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
    vigra::FImage rohr_corner_strength(w,h);
    
    // empty corner image
    corners = 0;
    ...
    
    // find corner response at scale 1.0
    vigra::rohrCornerDetector(srcImageRange(src), destImage(rohr_corner_strength), 
                              1.0);
    
    // find local maxima of corner response, mark with 1
    vigra::localMaxima(srcImageRange(rohr_corner_strength), destImage(corners));
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
rohrCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                       DestIterator dul, DestAccessor ad,
		       double scale)
{
    vigra_precondition(scale > 0.0,
                 "rohrCornerDetector(): Scale must be > 0");
		 
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    if(w <= 0 || h <= 0) return;
    
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
	
    typedef BasicImage<TmpType> TmpImage;
    
    TmpImage gx(w,h);
    TmpImage gy(w,h);
    TmpImage gxy(w,h);

    structureTensor(srcIterRange(sul, slr, as), 
                    destImage(gx), destImage(gxy), destImage(gy), 
                    scale, scale);
    RohrCornerFunctor<typename SrcAccessor::value_type > cf;
                    
    combineThreeImages(srcImageRange(gx), srcImage(gy), srcImage(gxy), 
                       destIter(dul, ad), cf );
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void rohrCornerDetector(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	   pair<DestIterator, DestAccessor> dest,
	   double scale)
{
    rohrCornerDetector(src.first, src.second, src.third,
                            dest.first, dest.second,
			    scale);
}

/********************************************************/
/*                                                      */
/*                 beaudetCornerDetector                */
/*                                                      */
/********************************************************/

/** \brief Find corners in an image (4).

    This algorithm implements a corner detector  
    according to [P.R. Beaudet: <em> "Rotationally Invariant Image Operators"</em>, 
    Proc. Intl. Joint Conf. on Pattern Recognition, Kyoto, Japan, 1978, pp. 579-583]. 
    
    The algorithm calculates the corner strength as the negative determinant of the 
    \link CommonConvolutionFilters#hessianMatrixOfGaussian Hessian Matrix\endlink. 
    The local maxima of the corner strength denote the corners in the gray level 
    image. 
    
    The source value type must be a linear algebra, i.e. addition, subtraction, and
    multiplication with itself, multiplication with doubles and 
    \ref NumericTraits "NumericTraits" must 
    be defined.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        void
        beaudetCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
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
        void beaudetCornerDetector(
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
    vigra::FImage beaudet_corner_strength(w,h);
    
    // empty corner image
    corners = 0;
    ...
    
    // find corner response at scale 1.0
    vigra::beaudetCornerDetector(srcImageRange(src), destImage(beaudet_corner_strength), 
                              1.0);
    
    // find local maxima of corner response, mark with 1
    vigra::localMaxima(srcImageRange(beaudet_corner_strength), destImage(corners));
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
beaudetCornerDetector(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                       DestIterator dul, DestAccessor ad,
		       double scale)
{
    vigra_precondition(scale > 0.0,
                 "beaudetCornerDetector(): Scale must be > 0");
		 
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    if(w <= 0 || h <= 0) return;
    
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
	
    typedef BasicImage<TmpType> TmpImage;
    
    TmpImage gx(w,h);
    TmpImage gy(w,h);
    TmpImage gxy(w,h);

    hessianMatrixOfGaussian(srcIterRange(sul, slr, as), 
                    destImage(gx), destImage(gxy), destImage(gy), 
                    scale);
    BeaudetCornerFunctor<typename SrcAccessor::value_type > cf;
                    
    combineThreeImages(srcImageRange(gx), srcImage(gy), srcImage(gxy), 
                       destIter(dul, ad), cf );
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void beaudetCornerDetector(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	   pair<DestIterator, DestAccessor> dest,
	   double scale)
{
    beaudetCornerDetector(src.first, src.second, src.third,
                            dest.first, dest.second,
			    scale);
}


//@}

} // namespace vigra

#endif // VIGRA_CORNERDETECTION_HXX
