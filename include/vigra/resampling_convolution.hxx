/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_RESAMPLING_CONVOLUTION_HXX
#define VIGRA_RESAMPLING_CONVOLUTION_HXX

#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/rational.hxx"

namespace vigra {

namespace resampling_detail
{

struct MapTargetToSourceCoordinate
{
    MapTargetToSourceCoordinate(Rational<int> const & samplingRatio, 
                                Rational<int> const & offset)
    : a(samplingRatio.denominator()*offset.denominator()),
      b(samplingRatio.numerator()*offset.numerator()),
      c(samplingRatio.numerator()*offset.denominator())
    {}
    
//        the following funcions are more efficient realizations of:
//             rational_cast<T>(i / samplingRatio + offset);
//        we need efficiency because this may be called in the inner loop

    int operator()(int i) const
    {
        return (i * a + b) / c;
    }
    
    double toDouble(int i) const
    {
        return double(i * a + b) / c;
    }

    Rational<int> toRational(int i) const
    {
        return Rational<int>(i * a + b, c);
    }

    int a, b, c;
};

} // namespace resampling_detail

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class KernelArray,
          class Functor>
void 
resamplingConvolveLine(SrcIter s, SrcIter send, SrcAcc src,
                       DestIter d, DestIter dend, DestAcc dest,
                       KernelArray const & kernels,
                       Functor mapTargetToSourceCoordinate)
{
    typedef typename 
        NumericTraits<typename SrcAcc::value_type>::RealPromote
        TmpType;
    typedef typename KernelArray::value_type Kernel;
    typedef typename Kernel::const_iterator KernelIter;
    
    int wo = send - s;
    int wn = dend - d;
    int wo2 = 2*wo - 2;
    
    int i;
    typename KernelArray::const_iterator kernel = kernels.begin();
    for(i=0; i<wn; ++i, ++d, ++kernel)
    {
        // use the kernels periodically
        if(kernel == kernels.end())
            kernel = kernels.begin();
        
        // calculate current target point into source location
        int is = mapTargetToSourceCoordinate(i);
        
        TmpType sum = NumericTraits<TmpType>::zero();

        int lbound = is - kernel->right(), 
            hbound = is - kernel->left();
                    
        KernelIter k = kernel->center() + kernel->right();
        if(lbound < 0 || hbound >= wo)
        {    
            vigra_precondition(-lbound < wo && wo2 - hbound >= 0,
                "resamplingConvolveLine(): kernel or offset larger than image.");
            for(int m=lbound; m <= hbound; ++m, --k)
            {
                int mm = (m < 0) ?
                            -m :
                            (m >= wo) ?
                                wo2 - m :
                                m;
                sum += *k * src(s, mm);
            }
        }
        else
        {
            SrcIter ss = s + lbound;
            SrcIter ssend = s + hbound;
            
            for(; ss <= ssend; ++ss, --k)
            {
                sum += *k * src(ss);
            }
        }
        
        dest.set(sum, d);
    }
}

template <class Kernel, class MapCoordinate, class KernelArray>
void
createResamplingKernels(Kernel const & kernel, 
             MapCoordinate const & mapCoordinate, KernelArray & kernels)
{
    for(unsigned int idest = 0; idest < kernels.size(); ++idest)
    {
        int isrc = mapCoordinate(idest);
        double idsrc = mapCoordinate.toDouble(idest);
        double offset = idsrc - isrc;
        double radius = kernel.radius();
        int left = int(ceil(-radius - offset));
        int right = int(floor(radius - offset));
        kernels[idest].initExplicitly(left, right);
        
        double x = left + offset;
        for(int i = left; i <= right; ++i, ++x)
            kernels[idest][i] = kernel(x);
        kernels[idest].normalize(1.0, kernel.derivativeOrder(), offset);
    }
}

/** \addtogroup ResamplingConvolutionFilters Resampling Convolution Filters

    These functions implement the convolution operation when the source and target images
    have different sizes. This is realized by accessing a continous kernel at the
    appropriate non-integer positions. The technique is, for example, described in
    D. Schumacher: <i>General Filtered Image Rescaling</i>, in: Graphics Gems III, 
    Academic Press, 1992.
*/
//@{

/********************************************************/
/*                                                      */
/*                  resamplingConvolveX                 */
/*                                                      */
/********************************************************/

/** \brief Apply a resampling filter in the x-direction.

    This function implements a convolution operation in x-direction
    (i.e. applies a 1D filter to every row) where the width of the source
    and destination images differ. This is typically used to avoid aliasing if 
    the image is scaled down, or to interpolate smoothly if the image is scaled up.
    The target coordinates are transformed into source coordinates by
    
    \code
    xsource = (xtarget - offset) / samplingRatio
    \endcode
    
    The <tt>samplingRatio</tt> and <tt>offset</tt> must be given as \ref vigra::Rational 
    in order to avoid rounding errors in this transformation. It is required that for all
    pixels of the target image, <tt>xsource</tt> remains within the range of the source
    image (i.e. <tt>0 <= xsource <= sourceWidth-1</tt>. Since <tt>xsource</tt> is
    in general not an integer, the <tt>kernel</tt> must be a functor that can be accessed at 
    arbitrary (<tt>double</tt>) coordinates. It must also provide a member function <tt>radius()</tt>
    which specifies the support (non-zero interval) of the kernel. VIGRA already
    provides a number of suitable functors, e.g. \ref vigra::Gaussian, \ref vigra::BSpline
    \ref vigra::CatmullRomSpline, and \ref vigra::CoscotFunction. The function
    \ref resizeImageSplineInterpolation() is implemented by means resamplingConvolveX() and 
    resamplingConvolveY().

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIter, class SrcAcc,
                  class DestIter, class DestAcc,
                  class Kernel>
        void 
        resamplingConvolveX(SrcIter sul, SrcIter slr, SrcAcc src,
                            DestIter dul, DestIter dlr, DestAcc dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIter, class SrcAcc,
                  class DestIter, class DestAcc,
                  class Kernel>
        void 
        resamplingConvolveX(triple<SrcIter, SrcIter, SrcAcc> src,
                            triple<DestIter, DestIter, DestAcc> dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="resampling_convolution_8hxx-source.html">vigra/resampling_convolution.hxx</a>"


    \code
    Rational<int> ratio(2), offset(0);

    FImage src(w,h), 
           dest(rational_cast<int>(ratio*w), h);

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simpultaneously enlarge and smooth source image
    resamplingConvolveX(srcImageRange(src), destImageRange(dest), 
                        smooth, ratio, offset);
    \endcode

    <b> Required Interface:</b>
    
    \code
    Kernel kernel;
    int kernelRadius = kernel.radius();
    double x = ...;  // must be <= radius()
    double value = kernel(x);
    \endcode
*/
template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
void 
resamplingConvolveX(SrcIter sul, SrcIter slr, SrcAcc src,
                    DestIter dul, DestIter dlr, DestAcc dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    int wold = slr.x - sul.x;
    int wnew = dlr.x - dul.x;
    
    vigra_precondition(!samplingRatio.is_inf() && samplingRatio > 0,
                "resamplingConvolveX(): sampling ratio must be > 0 and < infinity");
    vigra_precondition(!offset.is_inf(),
                "resamplingConvolveX(): offset must be < infinity");

    int period = lcm(samplingRatio.numerator(), samplingRatio.denominator());
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(period);
    
    createResamplingKernels(kernel, mapCoordinate, kernels);

    for(; sul.y < slr.y; ++sul.y, ++dul.y)
    {
        typename SrcIter::row_iterator sr = sul.rowIterator();
        typename DestIter::row_iterator dr = dul.rowIterator();
        resamplingConvolveLine(sr, sr+wold, src, dr, dr+wnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
inline void 
resamplingConvolveX(triple<SrcIter, SrcIter, SrcAcc> src,
                    triple<DestIter, DestIter, DestAcc> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveX(src.first, src.second, src.third,
                        dest.first, dest.second, dest.third,
                        kernel, samplingRatio, offset);
}

/********************************************************/
/*                                                      */
/*                  resamplingConvolveY                 */
/*                                                      */
/********************************************************/

/** \brief Apply a resampling filter in the y-direction.

    This function implements a convolution operation in y-direction
    (i.e. applies a 1D filter to every column) where the height of the source
    and destination images differ. This is typically used to avoid aliasing if 
    the image is scaled down, or to interpolate smoothly if the image is scaled up.
    The target coordinates are transformed into source coordinates by
    
    \code
    ysource = (ytarget - offset) / samplingRatio
    \endcode
    
    The <tt>samplingRatio</tt> and <tt>offset</tt> must be given as \ref vigra::Rational 
    in order to avoid rounding errors in this transformation. It is required that for all
    pixels of the target image, <tt>ysource</tt> remains within the range of the source
    image (i.e. <tt>0 <= ysource <= sourceHeight-1</tt>. Since <tt>ysource</tt> is
    in general not an integer, the <tt>kernel</tt> must be a functor that can be accessed at 
    arbitrary (<tt>double</tt>) coordinates. It must also provide a member function <tt>radius()</tt>
    which specifies the support (non-zero interval) of the kernel. VIGRA already
    provides a number of suitable functors, e.g. \ref vigra::Gaussian, \ref vigra::BSpline
    \ref vigra::CatmullRomSpline, and \ref vigra::CoscotFunction. The function
    \ref resizeImageSplineInterpolation() is implemented by means resamplingConvolveX() and 
    resamplingConvolveY().

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIter, class SrcAcc,
                  class DestIter, class DestAcc,
                  class Kernel>
        void 
        resamplingConvolveY(SrcIter sul, SrcIter slr, SrcAcc src,
                            DestIter dul, DestIter dlr, DestAcc dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIter, class SrcAcc,
                  class DestIter, class DestAcc,
                  class Kernel>
        void 
        resamplingConvolveY(triple<SrcIter, SrcIter, SrcAcc> src,
                            triple<DestIter, DestIter, DestAcc> dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="resampling_convolution_8hxx-source.html">vigra/resampling_convolution.hxx</a>"


    \code
    Rational<int> ratio(2), offset(0);

    FImage src(w,h), 
           dest(w, rational_cast<int>(ratio*h));

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simpultaneously enlarge and smooth source image
    resamplingConvolveY(srcImageRange(src), destImageRange(dest), 
                        smooth, ratio, offset);
    \endcode

    <b> Required Interface:</b>
    
    \code
    Kernel kernel;
    int kernelRadius = kernel.radius();
    double y = ...;  // must be <= radius()
    double value = kernel(y);
    \endcode
*/
template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
void 
resamplingConvolveY(SrcIter sul, SrcIter slr, SrcAcc src,
                    DestIter dul, DestIter dlr, DestAcc dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    int hold = slr.y - sul.y;
    int hnew = dlr.y - dul.y;
    
    vigra_precondition(!samplingRatio.is_inf() && samplingRatio > 0,
                "resamplingConvolveY(): sampling ratio must be > 0 and < infinity");
    vigra_precondition(!offset.is_inf(),
                "resamplingConvolveY(): offset must be < infinity");

    int period = lcm(samplingRatio.numerator(), samplingRatio.denominator());
    
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(period);
    
    createResamplingKernels(kernel, mapCoordinate, kernels);

    for(; sul.x < slr.x; ++sul.x, ++dul.x)
    {
        typename SrcIter::column_iterator sc = sul.columnIterator();
        typename DestIter::column_iterator dc = dul.columnIterator();
        resamplingConvolveLine(sc, sc+hold, src, dc, dc+hnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class Kernel>
inline void 
resamplingConvolveY(triple<SrcIter, SrcIter, SrcAcc> src,
                    triple<DestIter, DestIter, DestAcc> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveY(src.first, src.second, src.third,
                        dest.first, dest.second, dest.third,
                        kernel, samplingRatio, offset);
}

/********************************************************/
/*                                                      */
/*               resamplingConvolveImage                */
/*                                                      */
/********************************************************/

/** \brief Apply two separable resampling filters successively, the first in x-direction, 
           the second in y-direction.

    This function is a shorthand for the concatenation of a call to
    \link ResamplingConvolutionFilters#resamplingConvolveX resamplingConvolveX\endlink()
    and \link ResamplingConvolutionFilters#resamplingConvolveY resamplingConvolveY\endlink() 
    with the given kernels. See there for detailed documentation.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelX, class KernelY>
        void resamplingConvolveImage(SrcIterator sul,SrcIterator slr, SrcAccessor src,
                           DestIterator dul, DestIterator dlr, DestAccessor dest,
                           KernelX const & kx, 
                           Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                           KernelY const & ky, 
                           Rational<int> const & samplingRatioY, Rational<int> const & offsetY);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelX, class KernelY>
        void
        resamplingConvolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           triple<DestIterator, DestIterator, DestAccessor> dest,
                           KernelX const & kx, 
                           Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                           KernelY const & ky, 
                           Rational<int> const & samplingRatioY, Rational<int> const & offsetY);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="resampling_convolution_8hxx-source.html">vigra/resampling_convolution.hxx</a>"


    \code
    Rational<int> xratio(2), yratio(3), offset(0);

    FImage src(w,h), 
           dest(rational_cast<int>(xratio*w), rational_cast<int>(yratio*h));

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simpultaneously enlarge and smooth source image
    resamplingConvolveImage(srcImageRange(src), destImageRange(dest), 
                            smooth, xratio, offset,
                            smooth, yratio, offset);

    \endcode

*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelX, class KernelY>
void resamplingConvolveImage(SrcIterator sul,SrcIterator slr, SrcAccessor src,
                   DestIterator dul, DestIterator dlr, DestAccessor dest,
                   KernelX const & kx, 
                   Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                   KernelY const & ky, 
                   Rational<int> const & samplingRatioY, Rational<int> const & offsetY)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    
    BasicImage<TmpType> tmp(dlr.x - dul.x, slr.y - sul.y);

    resamplingConvolveX(srcIterRange(sul, slr, src),
                        destImageRange(tmp), 
                        kx, samplingRatioX, offsetX);
    resamplingConvolveY(srcImageRange(tmp),
                        destIterRange(dul, dlr, dest), 
                        ky, samplingRatioY, offsetY);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelX, class KernelY>
inline void
resamplingConvolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   triple<DestIterator, DestIterator, DestAccessor> dest,
                   KernelX const & kx, 
                   Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                   KernelY const & ky, 
                   Rational<int> const & samplingRatioY, Rational<int> const & offsetY)
{
    resamplingConvolveImage(src.first, src.second, src.third,
                            dest.first, dest.second, dest.third,
                            kx, samplingRatioX, offsetX,
                            ky, samplingRatioY, offsetY);
}

} // namespace vigra 


#endif /* VIGRA_RESAMPLING_CONVOLUTION_HXX */
