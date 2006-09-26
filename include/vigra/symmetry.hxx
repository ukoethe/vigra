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

#ifndef VIGRA_SYMMETRY_HXX
#define VIGRA_SYMMETRY_HXX

#include "utilities.hxx"
#include "numerictraits.hxx"
#include "stdimage.hxx"
#include "convolution.hxx"

namespace vigra {

/** \addtogroup SymmetryDetection Symmetry Detection
    Measure the local symmetry at each pixel.
*/
//@{

/********************************************************/
/*                                                      */
/*                 radialSymmetryTransform              */
/*                                                      */
/********************************************************/

/** \brief Find centers of radial symmetry in an image.

    This algorithm implements the Fast Radial Symmetry Transform according to
    [G. Loy, A. Zelinsky: <em> "A Fast Radial Symmetry Transform for Detecting
    Points of Interest"</em>, in: A. Heyden et al. (Eds.): Proc. of 7th European
    Conf. on Computer Vision, Part 1, pp. 358-368, Springer LNCS 2350, 2002].
    Minima of the algorithm response mark dark blobs, maxima correspond to light blobs.
    The "radial strictness parameter" is fixed at <TT>alpha</tt> = 2.0, the
    spatial spreading of the raw response is done by a Gaussian convolution
    at <tt>0.25*scale</TT> (these values are recommendations from the paper).
    Loy and Zelinsky additionally propose to add the operator response from several
    scales (see usage example below).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        void
        radialSymmetryTransform(SrcIterator sul, SrcIterator slr, SrcAccessor as,
			       DestIterator dul, DestAccessor ad,
			       double scale)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
	          class DestIterator, class DestAccessor>
        inline
        void radialSymmetryTransform(
	           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	           pair<DestIterator, DestAccessor> dest,
	           double scale)
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="symmetry_8hxx-source.html">vigra/symmetry.hxx</a>"<br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), centers(w,h);
    vigra::FImage symmetry(w,h);

    // empty result image
    centers.init(128);
    symmetry.init(0.0);

    // input width of edge detection filter
    for(double scale = 2.0; scale <= 8.0; scale *= 2.0)
    {
        vigra::FImage tmp(w,h);

        // find centers of symmetry
        radialSymmetryTransform(srcImageRange(src), destImage(tmp), scale);

        combineTwoImages(srcImageRange(symmetry), srcImage(tmp), destImage(symmetry),
                         std::plus<float>());
    }

    localMinima(srcImageRange(symmetry), destImage(centers), 0);
    localMaxima(srcImageRange(symmetry), destImage(centers), 255);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    // SrcAccessor::value_type must be a built-in type
    SrcAccessor::value_type u = src_accessor(src_upperleft);

    dest_accessor.set(u, dest_upperleft);
    \endcode
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
radialSymmetryTransform(SrcIterator sul, SrcIterator slr, SrcAccessor as,
               DestIterator dul, DestAccessor ad,
		double scale)
{
    vigra_precondition(scale > 0.0,
                 "radialSymmetryTransform(): Scale must be > 0");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    if(w <= 0 || h <= 0) return;

    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;

    typedef BasicImage<TmpType> TmpImage;
    typedef typename TmpImage::Iterator TmpIterator;

    TmpImage gx(w,h);
    TmpImage gy(w,h);
    IImage   orientationCounter(w,h);
    TmpImage magnitudeAccumulator(w,h);

    gaussianGradient(srcIterRange(sul, slr, as),
                     destImage(gx), destImage(gy),
                     scale);

    orientationCounter.init(0);
    magnitudeAccumulator.init(NumericTraits<TmpType>::zero());

    TmpIterator gxi = gx.upperLeft();
    TmpIterator gyi = gy.upperLeft();
    int y;
    for(y=0; y<h; ++y, ++gxi.y, ++gyi.y)
    {
        typename TmpIterator::row_iterator gxr = gxi.rowIterator();
        typename TmpIterator::row_iterator gyr = gyi.rowIterator();

        for(int x = 0; x<w; ++x, ++gxr, ++gyr)
        {
            double angle = VIGRA_CSTD::atan2(-*gyr, *gxr);
            double magnitude = VIGRA_CSTD::sqrt(*gxr * *gxr + *gyr * *gyr);

            if(magnitude < NumericTraits<TmpType>::epsilon()*10.0)
                continue;

            int dx = NumericTraits<int>::fromRealPromote(scale * VIGRA_CSTD::cos(angle));
            int dy = NumericTraits<int>::fromRealPromote(scale * VIGRA_CSTD::sin(angle));

            int xx = x + dx;
            int yy = y - dy;

            if(xx >= 0 && xx < w && yy >= 0 && yy < h)
            {
                orientationCounter(xx, yy) += 1;
                magnitudeAccumulator(xx, yy) += magnitude;
            }

            xx = x - dx;
            yy = y + dy;

            if(xx >= 0 && xx < w && yy >= 0 && yy < h)
            {
                orientationCounter(xx, yy) -= 1;
                magnitudeAccumulator(xx, yy) -= magnitude;
            }
        }
    }

    int maxOrientation = 0;
    TmpType maxMagnitude = NumericTraits<TmpType>::zero();

    for(y=0; y<h; ++y)
    {
        for(int x = 0; x<w; ++x)
        {
            int o = VIGRA_CSTD::abs(orientationCounter(x,y));

            if(o > maxOrientation)
                maxOrientation = o;

            TmpType m = VIGRA_CSTD::abs(magnitudeAccumulator(x,y));

            if(m > maxMagnitude)
                maxMagnitude = m;
        }
    }

    for(y=0; y<h; ++y)
    {
        for(int x = 0; x<w; ++x)
        {
            double o = (double)orientationCounter(x, y) / maxOrientation;
            magnitudeAccumulator(x, y) = o * o * magnitudeAccumulator(x, y) / maxMagnitude;
        }
    }

    gaussianSmoothing(srcImageRange(magnitudeAccumulator), destIter(dul, ad), 0.25*scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void radialSymmetryTransform(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
	   pair<DestIterator, DestAccessor> dest,
	   double scale)
{
    radialSymmetryTransform(src.first, src.second, src.third,
                            dest.first, dest.second,
			    scale);
}


//@}

} // namespace vigra


#endif /* VIGRA_SYMMETRY_HXX */
