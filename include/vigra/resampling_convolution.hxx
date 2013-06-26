/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_RESAMPLING_CONVOLUTION_HXX
#define VIGRA_RESAMPLING_CONVOLUTION_HXX

#include <cmath>
#include "stdimage.hxx"
#include "array_vector.hxx"
#include "rational.hxx"
#include "functortraits.hxx"
#include "functorexpression.hxx"
#include "transformimage.hxx"
#include "imagecontainer.hxx"
#include "multi_shape.hxx"

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

//        the following functions are more efficient realizations of:
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
    
    bool isExpand2() const
    {
        return a == 1 && b == 0 && c == 2;
    }
    
    bool isReduce2() const
    {
        return a == 2 && b == 0 && c == 1;
    }

    int a, b, c;
};

} // namespace resampling_detail

template <>
class FunctorTraits<resampling_detail::MapTargetToSourceCoordinate>
: public FunctorTraitsBase<resampling_detail::MapTargetToSourceCoordinate>
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class KernelArray>
void
resamplingExpandLine2(SrcIter s, SrcIter send, SrcAcc src,
                       DestIter d, DestIter dend, DestAcc dest,
                       KernelArray const & kernels)
{
    typedef typename KernelArray::value_type Kernel;
    typedef typename KernelArray::const_reference KernelRef;
    typedef typename Kernel::const_iterator KernelIter;

    typedef typename
        PromoteTraits<typename SrcAcc::value_type, typename Kernel::value_type>::Promote
        TmpType;

    int wo = send - s;
    int wn = dend - d;
    int wo2 = 2*wo - 2;
    
    int ileft = std::max(kernels[0].right(), kernels[1].right());
    int iright = wo + std::min(kernels[0].left(), kernels[1].left()) - 1;
    for(int i = 0; i < wn; ++i, ++d)
    {
        int is = i / 2;
        KernelRef kernel = kernels[i & 1];
        KernelIter k = kernel.center() + kernel.right();
        TmpType sum = NumericTraits<TmpType>::zero();        
        if(is < ileft)
        {
            for(int m=is-kernel.right(); m <= is-kernel.left(); ++m, --k)
            {
                int mm = (m < 0) 
                        ? -m 
                        : m;
                sum += *k * src(s, mm);
            }        
        }
        else if(is > iright)
        {
            for(int m=is-kernel.right(); m <= is-kernel.left(); ++m, --k)
            {
                int mm =  (m >= wo) 
                            ? wo2 - m
                            : m;
                sum += *k * src(s, mm);
            }        
        }
        else
        {
            SrcIter ss = s + is - kernel.right();
            for(int m = 0; m < kernel.size(); ++m, --k, ++ss)
            {
                sum += *k * src(ss);
            }        
        }
        dest.set(sum, d);
    }
}

template <class SrcIter, class SrcAcc,
          class DestIter, class DestAcc,
          class KernelArray>
void
resamplingReduceLine2(SrcIter s, SrcIter send, SrcAcc src,
                       DestIter d, DestIter dend, DestAcc dest,
                       KernelArray const & kernels)
{
    typedef typename KernelArray::value_type Kernel;
    typedef typename KernelArray::const_reference KernelRef;
    typedef typename Kernel::const_iterator KernelIter;

    KernelRef kernel = kernels[0];
    KernelIter kbegin = kernel.center() + kernel.right();

    typedef typename
        PromoteTraits<typename SrcAcc::value_type, typename Kernel::value_type>::Promote
        TmpType;

    int wo = send - s;
    int wn = dend - d;
    int wo2 = 2*wo - 2;
    
    int ileft = kernel.right();
    int iright = wo + kernel.left() - 1;
    for(int i = 0; i < wn; ++i, ++d)
    {
        int is = 2 * i;
        KernelIter k = kbegin;
        TmpType sum = NumericTraits<TmpType>::zero();        
        if(is < ileft)
        {
            for(int m=is-kernel.right(); m <= is-kernel.left(); ++m, --k)
            {
                int mm = (m < 0) 
                        ? -m 
                        : m;
                sum += *k * src(s, mm);
            }        
        }
        else if(is > iright)
        {
            for(int m=is-kernel.right(); m <= is-kernel.left(); ++m, --k)
            {
                int mm =  (m >= wo) 
                            ? wo2 - m
                            : m;
                sum += *k * src(s, mm);
            }        
        }
        else
        {
            SrcIter ss = s + is - kernel.right();
            for(int m = 0; m < kernel.size(); ++m, --k, ++ss)
            {
                sum += *k * src(ss);
            }        
        }
        dest.set(sum, d);
    }
}

/** \addtogroup ResamplingConvolutionFilters Resampling Convolution Filters

    These functions implement the convolution operation when the source and target images
    have different sizes. This is realized by accessing a continuous kernel at the
    appropriate non-integer positions. The technique is, for example, described in
    D. Schumacher: <i>General Filtered Image Rescaling</i>, in: Graphics Gems III,
    Academic Press, 1992.
*/
//@{

/********************************************************/
/*                                                      */
/*                resamplingConvolveLine                */
/*                                                      */
/********************************************************/

/** \brief Performs a 1-dimensional resampling convolution of the source signal using the given
    set of kernels.

    This function is mainly used internally: It is called for each dimension of a 
    higher dimensional array in order to perform a separable resize operation.

    <b> Declaration:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br/>
    Namespace: vigra

    \code
    namespace vigra {
        template <class SrcIter, class SrcAcc,
                  class DestIter, class DestAcc,
                  class KernelArray,
                  class Functor>
        void
        resamplingConvolveLine(SrcIter s, SrcIter send, SrcAcc src,
                               DestIter d, DestIter dend, DestAcc dest,
                               KernelArray const & kernels,
                               Functor mapTargetToSourceCoordinate)    
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void resamplingConvolveLine)

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
    if(mapTargetToSourceCoordinate.isExpand2())
    {
        resamplingExpandLine2(s, send, src, d, dend, dest, kernels);
        return;
    }
    if(mapTargetToSourceCoordinate.isReduce2())
    {
        resamplingReduceLine2(s, send, src, d, dend, dest, kernels);
        return;
    }
    
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
                sum = TmpType(sum + *k * src(s, mm));
            }
        }
        else
        {
            SrcIter ss = s + lbound;
            SrcIter ssend = s + hbound;

            for(; ss <= ssend; ++ss, --k)
            {
                sum = TmpType(sum + *k * src(ss));
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
        int left = std::min(0, int(ceil(-radius - offset)));
        int right = std::max(0, int(floor(radius - offset)));
        kernels[idest].initExplicitly(left, right);

        double x = left + offset;
        for(int i = left; i <= right; ++i, ++x)
            kernels[idest][i] = kernel(x);
        kernels[idest].normalize(1.0, kernel.derivativeOrder(), offset);
    }
}

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
    \ref resizeImageSplineInterpolation() is implemented by means of resamplingConvolveX() and
    resamplingConvolveY().

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class Kernel>
        void
        resamplingConvolveX(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode

    \deprecatedAPI{resamplingConvolveX}
    pass \ref ImageIterators and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br/>
    Namespace: vigra

    \code
    Rational<int> ratio(2), offset(0);

    MultiArray<2, float> src(w,h),
                         dest(rational_cast<int>(ratio*w), h);

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simultaneously enlarge and smooth source image
    resamplingConvolveX(src, dest,
                        smooth, ratio, offset);
    \endcode

    \deprecatedUsage{resamplingConvolveX}
    \code
    Rational<int> ratio(2), offset(0);

    FImage src(w,h),
           dest(rational_cast<int>(ratio*w), h);

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simultaneously enlarge and smooth source image
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void resamplingConvolveX)

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

template <class T1, class S1,
          class T2, class S2,
          class Kernel>
inline void
resamplingConvolveX(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveX(srcImageRange(src),
                        destImageRange(dest),
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
    \ref resizeImageSplineInterpolation() is implemented by means of resamplingConvolveX() and
    resamplingConvolveY().

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class Kernel>
        void
        resamplingConvolveY(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Kernel const & kernel,
                            Rational<int> const & samplingRatio, Rational<int> const & offset);
    }
    \endcode

    \deprecatedAPI{resamplingConvolveY}
    pass \ref ImageIterators and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br/>
    Namespace: vigra

    \code
    Rational<int> ratio(2), offset(0);

    MultiArray<2, float> src(w,h),
                         dest(w, rational_cast<int>(ratio*h));

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simultaneously enlarge and smooth source image
    resamplingConvolveY(src, dest,
                        smooth, ratio, offset);
    \endcode

    \deprecatedUsage{resamplingConvolveY}
    \code
    Rational<int> ratio(2), offset(0);

    FImage src(w,h),
           dest(w, rational_cast<int>(ratio*h));

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simultaneously enlarge and smooth source image
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void resamplingConvolveY)

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

template <class SrcIter, class SrcAccessor,
          class DestIter, class DestAccessor,
          class Kernel>
inline void
resamplingConvolveY(triple<SrcIter, SrcIter, SrcAccessor> src,
                    triple<DestIter, DestIter, DestAccessor> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveY(src.first, src.second, src.third,
                        dest.first, dest.second, dest.third,
                        kernel, samplingRatio, offset);
}

template <class T1, class S1,
          class T2, class S2,
          class Kernel>
inline void
resamplingConvolveY(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    Kernel const & kernel,
                    Rational<int> const & samplingRatio, Rational<int> const & offset)
{
    resamplingConvolveY(srcImageRange(src),
                        destImageRange(dest),
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
    \ref resamplingConvolveX() and \ref resamplingConvolveY()
    with the given kernels. See there for detailed documentation.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class KernelX, class KernelY>
        void
        resamplingConvolveImage(MultiArrayView<2, T1, S1> const & src,
                                MultiArrayView<2, T2, S2> dest,
                                KernelX const & kx,
                                Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                                KernelY const & ky,
                                Rational<int> const & samplingRatioY, Rational<int> const & offsetY);
    }
    \endcode

    \deprecatedAPI{resamplingConvolveImage}
    pass \ref ImageIterators and \ref DataAccessors :
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
    use argument objects in conjunction with \ref ArgumentObjectFactories :
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
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br/>
    Namespace: vigra

    \code
    Rational<int> xratio(2), yratio(3), offset(0);

    MultiArray<2, float> src(w,h),
                         dest(rational_cast<int>(xratio*w), rational_cast<int>(yratio*h));

    float sigma = 2.0;
    Gaussian<float> smooth(sigma);
    ...

    // simultaneously enlarge and smooth source image
    resamplingConvolveImage(src, dest,
                            smooth, xratio, offset,
                            smooth, yratio, offset);
    \endcode
*/
doxygen_overloaded_function(template <...> void resamplingConvolveImage)

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

template <class T1, class S1,
          class T2, class S2,
          class KernelX, class KernelY>
inline void
resamplingConvolveImage(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                        KernelX const & kx,
                        Rational<int> const & samplingRatioX, Rational<int> const & offsetX,
                        KernelY const & ky,
                        Rational<int> const & samplingRatioY, Rational<int> const & offsetY)
{
    resamplingConvolveImage(srcImageRange(src),
                            destImageRange(dest),
                            kx, samplingRatioX, offsetX,
                            ky, samplingRatioY, offsetY);
}

/** \brief Two-fold down-sampling for image pyramid construction.

    Sorry, no \ref detailedDocumentation() available yet.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidReduceBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                     DestIterator dul, DestIterator dlr, DestAccessor dest,
                                     double centerValue = 0.4);
    }
    \endcode

    \deprecatedAPI{pyramidReduceBurtFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidReduceBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                     DestIterator dul, DestIterator dlr, DestAccessor dest,
                                     double centerValue = 0.4);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidReduceBurtFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     triple<DestIterator, DestIterator, DestAccessor> dest,
                                     double centerValue = 0.4);
    }
    \endcode
    \deprecatedEnd

    use a \ref vigra::ImagePyramid :
    \code
    namespace vigra {
        template <class Image, class Alloc>
        void pyramidReduceBurtFilter(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                                     double centerValue = 0.4);
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void pyramidReduceBurtFilter)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void pyramidReduceBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                             DestIterator dul, DestIterator dlr, DestAccessor dest,
                             double centerValue = 0.4)
{
    vigra_precondition(0.25 <= centerValue && centerValue <= 0.5,
             "pyramidReduceBurtFilter(): centerValue must be between 0.25 and 0.5.");
             
    int wold = slr.x - sul.x;
    int wnew = dlr.x - dul.x;
    int hold = slr.y - sul.y;
    int hnew = dlr.y - dul.y;
    
    vigra_precondition(wnew == (wold + 1) / 2 && hnew == (hold + 1) / 2,
       "pyramidReduceBurtFilter(): oldSize = ceil(newSize / 2) required.");
    
    vigra_precondition(wnew == (wold + 1) / 2 && hnew == (hold + 1) / 2,
       "pyramidReduceBurtFilter(): oldSize = ceil(newSize / 2) required.");
    
    Rational<int> samplingRatio(1,2), offset(0);
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(1);
    kernels[0].initExplicitly(-2, 2) = 0.25 - centerValue / 2.0, 0.25, centerValue, 0.25, 0.25 - centerValue / 2.0;
   
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    typedef BasicImage<TmpType> TmpImage;
    typedef typename TmpImage::traverser TmpIterator;
    
    BasicImage<TmpType> tmp(wnew, hold);
    
    TmpIterator tul = tmp.upperLeft();

    for(; sul.y < slr.y; ++sul.y, ++tul.y)
    {
        typename SrcIterator::row_iterator sr = sul.rowIterator();
        typename TmpIterator::row_iterator tr = tul.rowIterator();
        // FIXME: replace with reduceLineBurtFilter()
        resamplingConvolveLine(sr, sr+wold, src, tr, tr+wnew, tmp.accessor(),
                               kernels, mapCoordinate);
    }
    
    tul  = tmp.upperLeft();

    for(; dul.x < dlr.x; ++dul.x, ++tul.x)
    {
        typename DestIterator::column_iterator dc = dul.columnIterator();
        typename TmpIterator::column_iterator tc = tul.columnIterator();
        resamplingConvolveLine(tc, tc+hold, tmp.accessor(), dc, dc+hnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void pyramidReduceBurtFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                             triple<DestIterator, DestIterator, DestAccessor> dest,
                             double centerValue = 0.4)
{
    pyramidReduceBurtFilter(src.first, src.second, src.third, 
                            dest.first, dest.second, dest.third, centerValue);
}

template <class Image, class Alloc>
inline
void pyramidReduceBurtFilter(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                             double centerValue = 0.4)
{
    vigra_precondition(fromLevel  < toLevel,
       "pyramidReduceBurtFilter(): fromLevel must be smaller than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= fromLevel && toLevel <= pyramid.highestLevel(),
       "pyramidReduceBurtFilter(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel+1; i <= toLevel; ++i)
        pyramidReduceBurtFilter(srcImageRange(pyramid[i-1]), destImageRange(pyramid[i]), centerValue);
}

/** \brief Two-fold up-sampling for image pyramid reconstruction.

    Sorry, no \ref detailedDocumentation() available yet.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidExpandBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                     DestIterator dul, DestIterator dlr, DestAccessor dest,
                                     double centerValue = 0.4);
    }
    \endcode

    \deprecatedAPI{pyramidExpandBurtFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidExpandBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                     DestIterator dul, DestIterator dlr, DestAccessor dest,
                                     double centerValue = 0.4);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void pyramidExpandBurtFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     triple<DestIterator, DestIterator, DestAccessor> dest,
                                     double centerValue = 0.4);
    }
    \endcode
    \deprecatedEnd

    use a \ref vigra::ImagePyramid :
    \code
    namespace vigra {
        template <class Image, class Alloc>
        void pyramidExpandBurtFilter(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                                     double centerValue = 0.4);
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void pyramidExpandBurtFilter)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void pyramidExpandBurtFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                             DestIterator dul, DestIterator dlr, DestAccessor dest,
                             double centerValue = 0.4)
{
    vigra_precondition(0.25 <= centerValue && centerValue <= 0.5,
             "pyramidExpandBurtFilter(): centerValue must be between 0.25 and 0.5.");
             
    int wold = slr.x - sul.x;
    int wnew = dlr.x - dul.x;
    int hold = slr.y - sul.y;
    int hnew = dlr.y - dul.y;
    
    vigra_precondition(wold == (wnew + 1) / 2 && hold == (hnew + 1) / 2,
       "pyramidExpandBurtFilter(): oldSize = ceil(newSize / 2) required.");
    
    vigra_precondition(wold == (wnew + 1) / 2 && hold == (hnew + 1) / 2,
       "pyramidExpandBurtFilter(): oldSize = ceil(newSize / 2) required.");
    
    Rational<int> samplingRatio(2), offset(0);
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(samplingRatio, offset);
    
    ArrayVector<Kernel1D<double> > kernels(2);
    kernels[0].initExplicitly(-1, 1) = 0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue;
    kernels[1].initExplicitly(-1, 0) = 0.5, 0.5;
   
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    typedef BasicImage<TmpType> TmpImage;
    typedef typename TmpImage::traverser TmpIterator;
    
    BasicImage<TmpType> tmp(wnew, hold);
    
    TmpIterator tul = tmp.upperLeft();

    for(; sul.y < slr.y; ++sul.y, ++tul.y)
    {
        typename SrcIterator::row_iterator sr = sul.rowIterator();
        typename TmpIterator::row_iterator tr = tul.rowIterator();
        // FIXME: replace with expandLineBurtFilter()
        resamplingConvolveLine(sr, sr+wold, src, tr, tr+wnew, tmp.accessor(),
                               kernels, mapCoordinate);
    }
    
    tul  = tmp.upperLeft();

    for(; dul.x < dlr.x; ++dul.x, ++tul.x)
    {
        typename DestIterator::column_iterator dc = dul.columnIterator();
        typename TmpIterator::column_iterator tc = tul.columnIterator();
        resamplingConvolveLine(tc, tc+hold, tmp.accessor(), dc, dc+hnew, dest,
                               kernels, mapCoordinate);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void pyramidExpandBurtFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                             triple<DestIterator, DestIterator, DestAccessor> dest,
                             double centerValue = 0.4)
{
    pyramidExpandBurtFilter(src.first, src.second, src.third, 
                            dest.first, dest.second, dest.third, centerValue);
}


template <class Image, class Alloc>
inline
void pyramidExpandBurtFilter(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                             double centerValue = 0.4)
{
    vigra_precondition(fromLevel  > toLevel,
       "pyramidExpandBurtFilter(): fromLevel must be larger than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= toLevel && fromLevel <= pyramid.highestLevel(),
       "pyramidExpandBurtFilter(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel-1; i >= toLevel; --i)
        pyramidExpandBurtFilter(srcImageRange(pyramid[i+1]), destImageRange(pyramid[i]), centerValue);
}

/** \brief Create a Laplacian pyramid.

    Sorry, no \ref detailedDocumentation() available yet.

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br>
    Namespace: vigra
*/
template <class Image, class Alloc>
inline void
pyramidReduceBurtLaplacian(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                           double centerValue = 0.4)
{
    using namespace functor;
    
    pyramidReduceBurtFilter(pyramid, fromLevel, toLevel, centerValue);
    for(int i=fromLevel; i < toLevel; ++i)
    {
        typename ImagePyramid<Image, Alloc>::value_type tmpImage(pyramid[i].size());
        pyramidExpandBurtFilter(srcImageRange(pyramid[i+1]), destImageRange(tmpImage), centerValue);
        combineTwoImages(srcImageRange(tmpImage), srcImage(pyramid[i]), destImage(pyramid[i]),
                       Arg1() - Arg2()); 
    }
}

/** \brief Reconstruct a Laplacian pyramid.

    Sorry, no \ref detailedDocumentation() available yet.

    <b>\#include</b> \<vigra/resampling_convolution.hxx\><br>
    Namespace: vigra
*/
template <class Image, class Alloc>
inline void
pyramidExpandBurtLaplacian(ImagePyramid<Image, Alloc> & pyramid, int fromLevel, int toLevel,
                           double centerValue = 0.4)
{
    using namespace functor;
    
    vigra_precondition(fromLevel  > toLevel,
       "pyramidExpandBurtLaplacian(): fromLevel must be larger than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= toLevel && fromLevel <= pyramid.highestLevel(),
       "pyramidExpandBurtLaplacian(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel-1; i >= toLevel; --i)
    {
        typename ImagePyramid<Image, Alloc>::value_type tmpImage(pyramid[i].size());
        pyramidExpandBurtFilter(srcImageRange(pyramid[i+1]), destImageRange(tmpImage), centerValue);
        combineTwoImages(srcImageRange(tmpImage), srcImage(pyramid[i]), destImage(pyramid[i]),
                       Arg1() - Arg2()); 
    }
}

//@}

} // namespace vigra


#endif /* VIGRA_RESAMPLING_CONVOLUTION_HXX */
