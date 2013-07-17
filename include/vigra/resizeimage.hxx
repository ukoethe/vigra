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


#ifndef VIGRA_RESIZEIMAGE_HXX
#define VIGRA_RESIZEIMAGE_HXX

#include <vector>
#include "utilities.hxx"
#include "numerictraits.hxx"
#include "stdimage.hxx"
#include "recursiveconvolution.hxx"
#include "separableconvolution.hxx"
#include "resampling_convolution.hxx"
#include "splines.hxx"
#include "multi_shape.hxx"

namespace vigra {

/*****************************************************************/
/*                                                               */
/*                         CoscotFunction                        */
/*                                                               */
/*****************************************************************/

/** The Coscot interpolation function.

    Implements the Coscot interpolation function proposed by Maria Magnusson Seger
    (maria@isy.liu.se) in the context of tomographic reconstruction. It provides a fast
    transition between the pass- and stop-bands and minimal ripple outside the transition
    region. Both properties are important for this application and can be tuned by the parameters
    <i>m</i> and <i>h</i> (with defaults 3 and 0.5). The function is defined by

    \f[ f_{m,h}(x) = \left\{ \begin{array}{ll}
                                   \frac{1}{2m}\sin(\pi x)\cot(\pi x / (2 m))(h + (1-h)\cos(\pi x/m)) & |x| \leq m \\
                                  0 & \mbox{otherwise}
                        \end{array}\right.
    \f]

    It can be used as a functor, and as a kernel for
    \ref resamplingConvolveImage() to create a differentiable interpolant
    of an image.

    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra

    \ingroup MathFunctions
*/
template <class T>
class CoscotFunction
{
  public:

        /** the kernel's value type
        */
    typedef T            value_type;
        /** the unary functor's argument type
        */
    typedef T            argument_type;
        /** the splines polynomial order
        */
    typedef T            result_type;

    CoscotFunction(unsigned int m = 3, double h = 0.5)
    : m_(m),
      h_(h)
    {}

        /** function (functor) call
        */
    result_type operator()(argument_type x) const
    {
        return x == 0.0 ?
                    1.0
                  : abs(x) < m_ ?
                        VIGRA_CSTD::sin(M_PI*x) / VIGRA_CSTD::tan(M_PI * x / 2.0 / m_) *
                             (h_ + (1.0 - h_) * VIGRA_CSTD::cos(M_PI * x / m_)) / 2.0 / m_
                      : 0.0;
    }

        /** index operator -- same as operator()
        */
    value_type operator[](value_type x) const
        { return operator()(x); }

        /** Radius of the function's support.
            Needed for  \ref resamplingConvolveImage(), equals m.
        */
    double radius() const
        { return m_; }

        /** Derivative order of the function: always 0.
        */
    unsigned int derivativeOrder() const
        { return 0; }

        /** Prefilter coefficients for compatibility with \ref vigra::BSpline.
            (array has zero length, since prefiltering is not necessary).
        */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }
    
  protected:
    static ArrayVector<double> prefilterCoefficients_;
    unsigned int m_;
    double h_;
};

template <class T>
ArrayVector<double> CoscotFunction<T>::prefilterCoefficients_;



/** \addtogroup GeometricTransformations Geometric Transformations
    Zoom up and down by repeating pixels, or using various interpolation schemes.

    See also: \ref resamplingConvolveImage(), \ref resampleImage(), \ref resizeMultiArraySplineInterpolation()

    <b>\#include</b> \<vigra/stdimagefunctions.hxx\><br>
    <b>or</b><br>
    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
*/
//@{

/********************************************************/
/*                                                      */
/*               resizeLineNoInterpolation              */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeLineNoInterpolation(SrcIterator i1, SrcIterator iend, SrcAccessor as,
                           DestIterator id, DestIterator idend, DestAccessor ad)
{
    int wold = iend - i1;
    int wnew = idend - id;

    if(wnew == 1)
    {
        ad.set(as(i1), id);
        return;
    }
    
    double dx = (double)(wold - 1) / (wnew - 1);
    double x = 0.5;
    for(; id != idend; ++id, x += dx)
    {
        int ix = (int)x;
        ad.set(as(i1, ix), id);
    }
}

/********************************************************/
/*                                                      */
/*              resizeImageNoInterpolation              */
/*                                                      */
/********************************************************/

/** \brief Resize image by repeating the nearest pixel values.

    This algorithm is very fast and does not require any arithmetic on
    the pixel types.

    The range of both the input and output images (resp. regions) must
    be given. Both images must have a size of at least 2x2 pixels. The
    scaling factors are then calculated accordingly. Destination
    pixels are directly copied from the appropriate source pixels.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        resizeImageNoInterpolation(MultiArrayView<2, T1, S1> const & src,
                                   MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{resizeImageNoInterpolation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void
        resizeImageNoInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
              DestImageIterator id, DestImageIterator idend, DestAccessor da)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void
        resizeImageNoInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w, h);
    MultiArray<2, float>         dest(w_new, h_new);
    
    resizeImageNoInterpolation(src, dest);
    \endcode

    \deprecatedUsage{resizeImageNoInterpolation}
    \code
    vigra::resizeImageNoInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(),
               dest.upperLeft(), dest.lowerRight(), dest.accessor());

    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    dest_accessor.set(src_accessor(src_upperleft), dest_upperleft);

    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    Source and destination must have at least 2 pixels along each axis.
*/
doxygen_overloaded_function(template <...> void resizeImageNoInterpolation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageNoInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;

    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;

    vigra_precondition((w > 1) && (h > 1),
                 "resizeImageNoInterpolation(): "
                 "Source image too small.\n");
    vigra_precondition((wnew > 1) && (hnew > 1),
                 "resizeImageNoInterpolation(): "
                 "Destination image too small.\n");

    typedef BasicImage<typename SrcAccessor::value_type> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    TmpImage tmp(w, hnew);

    TmpImageIterator yt = tmp.upperLeft();

    for(int x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        typename SrcIterator::column_iterator c1 = is.columnIterator();
        typename TmpImageIterator::column_iterator ct = yt.columnIterator();

        resizeLineNoInterpolation(c1, c1 + h, sa, ct, ct + hnew, tmp.accessor());
    }

    yt = tmp.upperLeft();

    for(int y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        typename TmpImageIterator::row_iterator rt = yt.rowIterator();

        resizeLineNoInterpolation(rt, rt + w, tmp.accessor(), rd, rd + wnew, da);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageNoInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageNoInterpolation(src.first, src.second, src.third,
                               dest.first, dest.second, dest.third);
}

template <class T1, class S1,
          class T2, class S2>
inline void
resizeImageNoInterpolation(MultiArrayView<2, T1, S1> const & src,
                           MultiArrayView<2, T2, S2> dest)
{
    resizeImageNoInterpolation(srcImageRange(src),
                               destImageRange(dest));
}

/********************************************************/
/*                                                      */
/*             resizeLineLinearInterpolation            */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeLineLinearInterpolation(SrcIterator i1, SrcIterator iend, SrcAccessor as,
                           DestIterator id, DestIterator idend, DestAccessor ad)
{
    int wold = iend - i1;
    int wnew = idend - id;

    if((wold <= 1) || (wnew <= 1)) return; // oder error ?

    typedef
        NumericTraits<typename DestAccessor::value_type> DestTraits;
    typedef typename DestTraits::RealPromote RealPromote;

    ad.set(DestTraits::fromRealPromote(as(i1)), id);
    ++id;

    --iend, --idend;
    ad.set(DestTraits::fromRealPromote(as(iend)), idend);

    double dx = (double)(wold - 1) / (wnew - 1);
    double x = dx;

    for(; id != idend; ++id, x += dx)
    {
        if(x >= 1.0)
        {
            int xx = (int)x;
            i1 += xx;
            x -= (double)xx;
        }
        double x1 = 1.0 - x;

        ad.set(DestTraits::fromRealPromote(RealPromote(x1 * as(i1) + x * as(i1, 1))), id);
    }
}

/********************************************************/
/*                                                      */
/*           resizeImageLinearInterpolation             */
/*                                                      */
/********************************************************/

/** \brief Resize image using linear interpolation.

    The function uses the standard separable bilinear interpolation algorithm to
    obtain a good compromise between quality and speed.

    The range must of both the input and output images (resp. regions)
    must be given. Both images must have a size of at
    least 2x2. The scaling factors are then calculated
    accordingly. If the source image is larger than the destination, it
    is smoothed (band limited) using a recursive
    exponential filter. The source value_type (SrcAccessor::value_type) must
    be a linear space, i.e. it must support addition, multiplication
    with a scalar real number and \ref NumericTraits "NumericTraits".

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        resizeImageLinearInterpolation(MultiArrayView<2, T1, S1> const & src,
                                       MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{resizeImageLinearInterpolation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void
        resizeImageLinearInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
              DestImageIterator id, DestImageIterator idend, DestAccessor da)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void
        resizeImageLinearInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              triple<DestImageIterator, DestImageIterator, DestAccessor> dest)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w, h);
    MultiArray<2, float>         dest(w_new, h_new);
    
    resizeImageLinearInterpolation(src, dest);
    \endcode

    \deprecatedUsage{resizeImageLinearInterpolation}
    \code
    vigra::resizeImageLinearInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(),
               dest.upperLeft(), dest.lowerRight(), dest.accessor());

    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    NumericTraits<SrcAccessor::value_type>::RealPromote
                             u = src_accessor(src_upperleft),
                 v = src_accessor(src_upperleft, 1);
    double d;

    u = d * v;
    u = u + v;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(u),
    dest_upperleft);

    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    Source and destination must have at least 2 pixels along each axis.
*/
doxygen_overloaded_function(template <...> void resizeImageLinearInterpolation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageLinearInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    int w = iend.x - is.x;
    int h = iend.y - is.y;

    int wnew = idend.x - id.x;
    int hnew = idend.y - id.y;

    vigra_precondition((w > 1) && (h > 1),
                 "resizeImageLinearInterpolation(): "
                 "Source image too small.\n");
    vigra_precondition((wnew > 1) && (hnew > 1),
                 "resizeImageLinearInterpolation(): "
                 "Destination image too small.\n");

    double const scale = 2.0;

    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    typedef BasicImage<TMPTYPE> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<TMPTYPE> tmp(w, hnew);
    BasicImage<TMPTYPE> line((h > w) ? h : w, 1);

    int x,y;

    typename BasicImage<TMPTYPE>::Iterator yt = tmp.upperLeft();
    typename TmpImageIterator::row_iterator lt = line.upperLeft().rowIterator();

    for(x=0; x<w; ++x, ++is.x, ++yt.x)
    {
        typename SrcIterator::column_iterator c1 = is.columnIterator();
        typename TmpImageIterator::column_iterator ct = yt.columnIterator();

        if(hnew < h)
        {
            recursiveSmoothLine(c1, c1 + h, sa,
                 lt, line.accessor(), (double)h/hnew/scale);

            resizeLineLinearInterpolation(lt, lt + h, line.accessor(),
                                          ct, ct + hnew, tmp.accessor());
        }
        else
        {
            resizeLineLinearInterpolation(c1, c1 + h, sa,
                                          ct, ct + hnew, tmp.accessor());
        }
    }

    yt = tmp.upperLeft();

    for(y=0; y < hnew; ++y, ++yt.y, ++id.y)
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        typename TmpImageIterator::row_iterator rt = yt.rowIterator();

        if(wnew < w)
        {
            recursiveSmoothLine(rt, rt + w, tmp.accessor(),
                              lt, line.accessor(), (double)w/wnew/scale);

            resizeLineLinearInterpolation(lt, lt + w, line.accessor(),
                                          rd, rd + wnew, da);
        }
        else
        {
            resizeLineLinearInterpolation(rt, rt + w, tmp.accessor(),
                                          rd, rd + wnew, da);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageLinearInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageLinearInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

template <class T1, class S1,
          class T2, class S2>
inline void
resizeImageLinearInterpolation(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest)
{
    resizeImageLinearInterpolation(srcImageRange(src),
                                   destImageRange(dest));
}

/***************************************************************/
/*                                                             */
/*                resizeImageSplineInterpolation               */
/*                                                             */
/***************************************************************/

/** \brief Resize image using B-spline interpolation.

    The function implements separable spline interpolation algorithm described in

    M. Unser, A. Aldroubi, M. Eden, <i>"B-Spline Signal Processing"</i>
    IEEE Transactions on Signal Processing, vol. 41, no. 2, pp. 821-833 (part I),
    pp. 834-848 (part II), 1993.

    to obtain optimal interpolation quality and speed. You may pass the function
    a spline of arbitrary order (e.g. <TT>BSpline<ORDER, double></tt> or
    <TT>CatmullRomSpline<double></tt>). The default is a third order spline
    which gives a twice continuously differentiable interpolant.
    The implementation ensures that image values are interpolated rather
    than smoothed by first calling a recursive (sharpening) prefilter as
    described in the above paper. Then the actual interpolation is done
    using \ref resamplingConvolveLine().

    The range of both the input and output images (resp. regions)
    must be given. The input image must have a size of at
    least 4x4, the destination of at least 2x2. The scaling factors are then calculated
    accordingly. If the source image is larger than the destination, it
    is smoothed (band limited) using a recursive
    exponential filter. The source value_type (SrcAccessor::value_type) must
    be a linear algebra, i.e. it must support addition, subtraction,
    and multiplication (+, -, *), multiplication with a scalar
    real number and \ref NumericTraits "NumericTraits".
    The function uses accessors.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class SPLINE>
        void
        resizeImageSplineInterpolation(MultiArrayView<2, T1, S1> const & src,
                                       MultiArrayView<2, T2, S2> dest,
                                       SPLINE const & spline = BSpline<3, double>());
    }
    \endcode

    \deprecatedAPI{resizeImageSplineInterpolation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class SPLINE>
        void
        resizeImageSplineInterpolation(
              SrcImageIterator is, SrcImageIterator iend, SrcAccessor sa,
              DestImageIterator id, DestImageIterator idend, DestAccessor da,
              SPLINE spline = BSpline<3, double>())
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class SPLINE>
        void
        resizeImageSplineInterpolation(
              triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              triple<DestImageIterator, DestImageIterator, DestAccessor> dest,
              SPLINE spline = BSpline<3, double>())
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, unsigned char> src(w, h);
    MultiArray<2, float>         dest(w_new, h_new);
    
    // use default cubic spline interpolator
    resizeImageSplineInterpolation(src, dest);
    
    // use 5th-order spline interpolator
    resizeImageSplineInterpolation(src, dest, BSpline<5, double>());
    \endcode

    \deprecatedUsage{resizeImageSplineInterpolation}
    \code
    vigra::resizeImageSplineInterpolation(
               src.upperLeft(), src.lowerRight(), src.accessor(),
               dest.upperLeft(), dest.lowerRight(), dest.accessor());

    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft, src_lowerright;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    NumericTraits<SrcAccessor::value_type>::RealPromote
                             u = src_accessor(src_upperleft),
                 v = src_accessor(src_upperleft, 1);
    double d;

    u = d * v;
    u = u + v;
    u = u - v;
    u = u * v;
    u += v;
    u -= v;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(u),
    dest_upperleft);

    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    Source and destination must have at least 2 pixels along each axis.
*/
doxygen_overloaded_function(template <...> void resizeImageSplineInterpolation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class SPLINE>
void
resizeImageSplineInterpolation(
    SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
    DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc,
    SPLINE const & spline)
{

    int width_old = src_iter_end.x - src_iter.x;
    int height_old = src_iter_end.y - src_iter.y;

    int width_new = dest_iter_end.x - dest_iter.x;
    int height_new = dest_iter_end.y - dest_iter.y;

    vigra_precondition((width_old > 1) && (height_old > 1),
                 "resizeImageSplineInterpolation(): "
                 "Source image too small.\n");

    vigra_precondition((width_new > 1) && (height_new > 1),
                 "resizeImageSplineInterpolation(): "
                 "Destination image too small.\n");

    Rational<int> xratio(width_new - 1, width_old - 1);
    Rational<int> yratio(height_new - 1, height_old - 1);
    Rational<int> offset(0);
    resampling_detail::MapTargetToSourceCoordinate xmapCoordinate(xratio, offset);
    resampling_detail::MapTargetToSourceCoordinate ymapCoordinate(yratio, offset);
    int xperiod = lcm(xratio.numerator(), xratio.denominator());
    int yperiod = lcm(yratio.numerator(), yratio.denominator());

    double const scale = 2.0;

    typedef typename SrcAccessor::value_type SRCVT;
    typedef typename NumericTraits<SRCVT>::RealPromote TMPTYPE;
    typedef BasicImage<TMPTYPE> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<TMPTYPE> tmp(width_old, height_new);

    BasicImage<TMPTYPE> line((height_old > width_old) ? height_old : width_old, 1);
    typename BasicImage<TMPTYPE>::Accessor tmp_acc = tmp.accessor();
    ArrayVector<double> const & prefilterCoeffs = spline.prefilterCoefficients();

    int x,y;

    ArrayVector<Kernel1D<double> > kernels(yperiod);
    createResamplingKernels(spline, ymapCoordinate, kernels);

    typename BasicImage<TMPTYPE>::Iterator y_tmp = tmp.upperLeft();
    typename TmpImageIterator::row_iterator line_tmp = line.upperLeft().rowIterator();

    for(x=0; x<width_old; ++x, ++src_iter.x, ++y_tmp.x)
    {

        typename SrcIterator::column_iterator c_src = src_iter.columnIterator();
        typename TmpImageIterator::column_iterator c_tmp = y_tmp.columnIterator();

        if(prefilterCoeffs.size() == 0)
        {
            if(height_new >= height_old)
            {
                resamplingConvolveLine(c_src, c_src + height_old, src_acc,
                                       c_tmp, c_tmp + height_new, tmp_acc,
                                       kernels, ymapCoordinate);
            }
            else
            {
                recursiveSmoothLine(c_src, c_src + height_old, src_acc,
                     line_tmp, line.accessor(), (double)height_old/height_new/scale);
                resamplingConvolveLine(line_tmp, line_tmp + height_old, line.accessor(),
                                       c_tmp, c_tmp + height_new, tmp_acc,
                                       kernels, ymapCoordinate);
            }
        }
        else
        {
            recursiveFilterLine(c_src, c_src + height_old, src_acc,
                                line_tmp, line.accessor(),
                                prefilterCoeffs[0], BORDER_TREATMENT_REFLECT);
            for(unsigned int b = 1; b < prefilterCoeffs.size(); ++b)
            {
                recursiveFilterLine(line_tmp, line_tmp + height_old, line.accessor(),
                                    line_tmp, line.accessor(),
                                    prefilterCoeffs[b], BORDER_TREATMENT_REFLECT);
            }
            if(height_new < height_old)
            {
                recursiveSmoothLine(line_tmp, line_tmp + height_old, line.accessor(),
                     line_tmp, line.accessor(), (double)height_old/height_new/scale);
            }
            resamplingConvolveLine(line_tmp, line_tmp + height_old, line.accessor(),
                                   c_tmp, c_tmp + height_new, tmp_acc,
                                   kernels, ymapCoordinate);
        }
    }

    y_tmp = tmp.upperLeft();

    kernels.resize(xperiod);
    createResamplingKernels(spline, xmapCoordinate, kernels);

    for(y=0; y < height_new; ++y, ++y_tmp.y, ++dest_iter.y)
    {
        typename DestIterator::row_iterator r_dest = dest_iter.rowIterator();
        typename TmpImageIterator::row_iterator r_tmp = y_tmp.rowIterator();

        if(prefilterCoeffs.size() == 0)
        {
            if(width_new >= width_old)
            {
                resamplingConvolveLine(r_tmp, r_tmp + width_old, tmp.accessor(),
                                       r_dest, r_dest + width_new, dest_acc,
                                       kernels, xmapCoordinate);
            }
            else
            {
                recursiveSmoothLine(r_tmp, r_tmp + width_old, tmp.accessor(),
                                  line_tmp, line.accessor(), (double)width_old/width_new/scale);
                resamplingConvolveLine(line_tmp, line_tmp + width_old, line.accessor(),
                                       r_dest, r_dest + width_new, dest_acc,
                                       kernels, xmapCoordinate);
            }
        }
        else
        {
            recursiveFilterLine(r_tmp, r_tmp + width_old, tmp.accessor(),
                                line_tmp, line.accessor(),
                                prefilterCoeffs[0], BORDER_TREATMENT_REFLECT);
            for(unsigned int b = 1; b < prefilterCoeffs.size(); ++b)
            {
                recursiveFilterLine(line_tmp, line_tmp + width_old, line.accessor(),
                                    line_tmp, line.accessor(),
                                    prefilterCoeffs[b], BORDER_TREATMENT_REFLECT);
            }
            if(width_new < width_old)
            {
                recursiveSmoothLine(line_tmp, line_tmp + width_old, line.accessor(),
                                    line_tmp, line.accessor(), (double)width_old/width_new/scale);
            }
            resamplingConvolveLine(line_tmp, line_tmp + width_old, line.accessor(),
                                   r_dest, r_dest + width_new, dest_acc,
                                   kernels, xmapCoordinate);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageSplineInterpolation(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestIterator idend, DestAccessor da)
{
    resizeImageSplineInterpolation(is, iend, sa, id, idend, da, BSpline<3, double>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class SPLINE>
inline void
resizeImageSplineInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               triple<DestIterator, DestIterator, DestAccessor> dest,
                               SPLINE const & spline)
{
    resizeImageSplineInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third, spline);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageSplineInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageSplineInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

template <class T1, class S1,
          class T2, class S2,
          class SPLINE>
inline void
resizeImageSplineInterpolation(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest,
                               SPLINE const & spline)
{
    resizeImageSplineInterpolation(srcImageRange(src),
                                   destImageRange(dest), spline);
}

template <class T1, class S1,
          class T2, class S2>
inline void
resizeImageSplineInterpolation(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest)
{
    resizeImageSplineInterpolation(srcImageRange(src),
                                   destImageRange(dest));
}

/*****************************************************************/
/*                                                               */
/*              resizeImageCatmullRomInterpolation               */
/*                                                               */
/*****************************************************************/

/** \brief Resize image using the Catmull/Rom interpolation function.

    The function calls like \ref resizeImageSplineInterpolation() with
    \ref vigra::CatmullRomSpline as an interpolation kernel.
    The interpolated function has one continuous derivative.
    (See \ref resizeImageSplineInterpolation() for more documentation)

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        resizeImageCatmullRomInterpolation(MultiArrayView<2, T1, S1> const & src,
                                           MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{resizeImageCatmullRomInterpolation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        resizeImageCatmullRomInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                              DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        resizeImageCatmullRomInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                              triple<DestIterator, DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd

    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra
    
    \code
    MultiArray<2, unsigned char> src(w, h);
    MultiArray<2, float>         dest(w_new, h_new);
    
    resizeImageCatmullRomInterpolation(src, dest);
    \endcode
*/
doxygen_overloaded_function(template <...> void resizeImageCatmullRomInterpolation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageCatmullRomInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                      DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{
    resizeImageSplineInterpolation(src_iter, src_iter_end, src_acc, dest_iter, dest_iter_end, dest_acc,
                                  CatmullRomSpline<double>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageCatmullRomInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageCatmullRomInterpolation(src.first, src.second, src.third,
                                       dest.first, dest.second, dest.third);
}

template <class T1, class S1,
          class T2, class S2>
inline void
resizeImageCatmullRomInterpolation(MultiArrayView<2, T1, S1> const & src,
                                   MultiArrayView<2, T2, S2> dest)
{
    resizeImageCatmullRomInterpolation(srcImageRange(src),
                                       destImageRange(dest));
}

/*****************************************************************/
/*                                                               */
/*              resizeImageCoscotInterpolation                   */
/*                                                               */
/*****************************************************************/

/** \brief Resize image using the Coscot interpolation function.

    The function calls \ref resizeImageSplineInterpolation() with
    \ref vigra::CoscotFunction as an interpolation kernel.
    The interpolated function has one continuous derivative.
    (See \ref resizeImageSplineInterpolation() for more documentation)

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        resizeImageCoscotInterpolation(MultiArrayView<2, T1, S1> const & src,
                                       MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{resizeImageCoscotInterpolation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        resizeImageCoscotInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                              DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        resizeImageCoscotInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                              triple<DestIterator, DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd


    <b>\#include</b> \<vigra/resizeimage.hxx\><br>
    Namespace: vigra
    
    \code
    MultiArray<2, unsigned char> src(w, h);
    MultiArray<2, float>         dest(w_new, h_new);
    
    resizeImageCoscotInterpolation(src, dest);
    \endcode
*/
doxygen_overloaded_function(template <...> void resizeImageCoscotInterpolation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
resizeImageCoscotInterpolation(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                      DestIterator dest_iter, DestIterator dest_iter_end, DestAccessor dest_acc)
{
    resizeImageSplineInterpolation(src_iter, src_iter_end, src_acc, dest_iter, dest_iter_end, dest_acc,
                                   CoscotFunction<double>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeImageCoscotInterpolation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               triple<DestIterator, DestIterator, DestAccessor> dest)
{
    resizeImageCoscotInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

template <class T1, class S1,
          class T2, class S2>
inline void
resizeImageCoscotInterpolation(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest)
{
    resizeImageCoscotInterpolation(srcImageRange(src),
                                   destImageRange(dest));
}

//@}

} // namespace vigra

#endif // VIGRA_RESIZEIMAGE_HXX
