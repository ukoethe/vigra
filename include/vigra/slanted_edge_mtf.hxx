/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2006 by Ullrich Koethe                  */
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


#ifndef VIGRA_SLANTED_EDGE_MTF_HXX
#define VIGRA_SLANTED_EDGE_MTF_HXX

#include <algorithm>
#include "array_vector.hxx"
#include "basicgeometry.hxx"
#include "edgedetection.hxx"
#include "fftw3.hxx"
#include "functorexpression.hxx"
#include "linear_solve.hxx"
#include "mathutil.hxx"
#include "numerictraits.hxx"
#include "separableconvolution.hxx"
#include "static_assert.hxx"
#include "stdimage.hxx"
#include "transformimage.hxx"
#include "utilities.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup SlantedEdgeMTF Camera MTF Estimation
    Determine the magnitude transfer function (MTF) of a camera using the slanted edge method.
*/
//@{

/********************************************************/
/*                                                      */
/*                  SlantedEdgeMTFOptions               */
/*                                                      */
/********************************************************/

/** \brief Pass options to one of the \ref slantedEdgeMTF() functions.

    <tt>SlantedEdgeMTFOptions</tt>  is an argument objects that holds various optional
    parameters used by the \ref slantedEdgeMTF() functions. If a parameter is not explicitly
    set, a suitable default will be used. Changing the defaults is only necessary if you can't
    obtain good input data, but absolutely need an MTF estimate.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/slanted_edge_mtf.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h);
    std::vector<vigra::TinyVector<double, 2> > mtf;

    ...
    slantedEdgeMTF(src, mtf,
                   SlantedEdgeMTFOptions().mtfSmoothingScale(1.0));

    // print the frequency / attenuation pairs found
    for(int k=0; k<result.size(); ++k)
        std::cout << "frequency: " << mtf[k][0] << ", estimated attenuation: " << mtf[k][1] << std::endl;
    \endcode
*/

class SlantedEdgeMTFOptions
{
  public:
        /** Initialize all options with default values.
        */
    SlantedEdgeMTFOptions()
    : minimum_number_of_lines(20),
      desired_edge_width(10),
      minimum_edge_width(5),
      mtf_smoothing_scale(2.0)
    {}

        /** Minimum number of pixels the edge must cross.

            The longer the edge the more accurate the resulting MTF estimate. If you don't have good
            data, but absolutely have to compute an MTF, you may force a lower value here.<br>
            Default: 20
        */
    SlantedEdgeMTFOptions & minimumNumberOfLines(unsigned int n)
    {
        minimum_number_of_lines = n;
        return *this;
    }

        /** Desired number of pixels perpendicular to the edge.

            The larger the regions to either side of the edge,
            the more accurate the resulting MTF estimate. If you don't have good
            data, but absolutely have to compute an MTF, you may force a lower value here.<br>
            Default: 10
        */
    SlantedEdgeMTFOptions & desiredEdgeWidth(unsigned int n)
    {
        desired_edge_width = n;
        return *this;
    }

        /** Minimum acceptable number of pixels perpendicular to the edge.

            The larger the regions to either side of the edge,
            the more accurate the resulting MTF estimate. If you don't have good
            data, but absolutely have to compute an MTF, you may force a lower value here.<br>
            Default: 5
        */
    SlantedEdgeMTFOptions & minimumEdgeWidth(unsigned int n)
    {
        minimum_edge_width = n;
        return *this;
    }

        /** Amount of smoothing of the computed MTF.

            If the data is noisy, so will be the MTF. Thus, some smoothing is useful.<br>
            Default: 2.0
        */
    SlantedEdgeMTFOptions & mtfSmoothingScale(double scale)
    {
        vigra_precondition(scale >= 0.0,
            "SlantedEdgeMTFOptions: MTF smoothing scale must not be < 0.");
        mtf_smoothing_scale = scale;
        return *this;
    }

    unsigned int minimum_number_of_lines, desired_edge_width, minimum_edge_width;
    double mtf_smoothing_scale;
};

//@}

namespace detail {

struct SortEdgelsByStrength
{
    bool operator()(Edgel const & l, Edgel const & r) const
    {
        return l.strength > r.strength;
    }
};

    /* Make sure that the edge runs vertically, intersects the top and bottom border
       of the image, and white is on the left.
    */
template <class SrcIterator, class SrcAccessor, class DestImage>
ptrdiff_t
prepareSlantedEdgeInput(SrcIterator sul, SrcIterator slr, SrcAccessor src, DestImage & res,
                        SlantedEdgeMTFOptions const & options)
{
    ptrdiff_t w = slr.x - sul.x;
    ptrdiff_t h = slr.y - sul.y;

    // find rough estimate of the edge
    ArrayVector<Edgel> edgels;
    cannyEdgelList(sul, slr, src, edgels, 2.0);
    std::sort(edgels.begin(), edgels.end(), SortEdgelsByStrength());

    double x = 0.0, y = 0.0, x2 = 0.0, y2 = 0.0, xy = 0.0;
    ptrdiff_t c = std::min(w, h);

    for(ptrdiff_t k = 0; k < c; ++k)
    {
        x += edgels[k].x;
        y += edgels[k].y;
        x2 += sq(edgels[k].x);
        xy += edgels[k].x*edgels[k].y;
        y2 += sq(edgels[k].y);
    }
    double xc = x / c, yc = y / c;
    x2 = x2 / c - sq(xc);
    xy = xy / c - xc*yc;
    y2 = y2 / c - sq(yc);
    double angle = 0.5*VIGRA_CSTD::atan2(2*xy, x2 - y2);

    DestImage tmp;
    // rotate image when slope is less than +-45 degrees
    if(VIGRA_CSTD::fabs(angle) < M_PI / 4.0)
    {
        xc = yc;
        yc = w - xc - 1;
        std::swap(w, h);
        tmp.resize(w, h);
        rotateImage(srcIterRange(sul, slr, src), destImage(tmp), 90);
        angle += M_PI / 2.0;
    }
    else
    {
        tmp.resize(w, h);
        copyImage(srcIterRange(sul, slr, src), destImage(tmp));
        if(angle < 0.0)
            angle += M_PI;
    }
    // angle is now between pi/4 and 3*pi/4
    double slope = VIGRA_CSTD::tan(M_PI/2.0 - angle);
    vigra_precondition(slope != 0.0,
          "slantedEdgeMTF(): Input edge is not slanted");

    // trim image so that the edge only intersects the top and bottom border
    ptrdiff_t minimumNumberOfLines = options.minimum_number_of_lines, //20,
                 edgeWidth = options.desired_edge_width, // 10
                 minimumEdgeWidth = options.minimum_edge_width; // 5

    ptrdiff_t y0 = 0, y1 = h;
    for(; edgeWidth >= minimumEdgeWidth; --edgeWidth)
    {
        y0 = ptrdiff_t(VIGRA_CSTD::floor((edgeWidth - xc) / slope + yc + 0.5));
        y1 = ptrdiff_t(VIGRA_CSTD::floor((w - edgeWidth - 1 - xc) / slope + yc + 0.5));
        if(slope < 0.0)
            std::swap(y0, y1);
        if(y1 - y0 >= (ptrdiff_t)minimumNumberOfLines)
            break;
    }

    vigra_precondition(edgeWidth >= minimumEdgeWidth,
        "slantedEdgeMTF(): Input edge is too slanted or image is too small");

    y0 = std::max(y0, ptrdiff_t(0));
    y1 = std::min(y1+1, h);

    res.resize(w, y1-y0);

    // ensure that white is on the left
    if(tmp(0,0) < tmp(w-1, h-1))
    {
        rotateImage(srcIterRange(tmp.upperLeft() + Diff2D(0, y0), tmp.upperLeft() + Diff2D(w, y1), tmp.accessor()),
                    destImage(res), 180);
    }
    else
    {
        copyImage(srcImageRange(tmp), destImage(res));
    }
    return edgeWidth;
}

template <class Image>
void slantedEdgeShadingCorrection(Image & i, ptrdiff_t edgeWidth)
{
    using namespace functor;

    // after prepareSlantedEdgeInput(), the white region is on the left
    // find a plane that approximates the logarithm of the white ROI

    transformImage(srcImageRange(i), destImage(i), log(Arg1() + Param(1.0)));

    ptrdiff_t w = i.width(),
              h = i.height();

    Matrix<double> m(3,3), r(3, 1), l(3, 1);
    for(ptrdiff_t y = 0; y < h; ++y)
    {
        for(ptrdiff_t x = 0; x < edgeWidth; ++x)
        {
            l(0,0) = x;
            l(1,0) = y;
            l(2,0) = 1.0;
            m += outer(l);
            r += i(x,y)*l;
        }
    }

    linearSolve(m, r, l);
    double a = l(0,0),
           b = l(1,0),
           c = l(2,0);

    // subtract the plane and go back to the non-logarithmic representation
    for(ptrdiff_t y = 0; y < h; ++y)
    {
        for(ptrdiff_t x = 0; x < w; ++x)
        {
            i(x, y) = VIGRA_CSTD::exp(i(x,y) - a*x - b*y - c);
        }
    }
}

template <class Image, class BackInsertable>
void slantedEdgeSubpixelShift(Image const & img, BackInsertable & centers, double & angle,
                              SlantedEdgeMTFOptions const & options)
{
    ptrdiff_t w = img.width();
    ptrdiff_t h = img.height();

    Image grad(w, h);
    Kernel1D<double> kgrad;
    kgrad.initGaussianDerivative(1.0, 1);

    separableConvolveX(srcImageRange(img), destImage(grad), kernel1d(kgrad));

    ptrdiff_t desiredEdgeWidth = (ptrdiff_t)options.desired_edge_width;
    double sy = 0.0, sx = 0.0, syy = 0.0, sxy = 0.0;
    for(ptrdiff_t y = 0; y < h; ++y)
    {
        double a = 0.0,
               b = 0.0;
        for(ptrdiff_t x = 0; x < w; ++x)
        {
            a += x*grad(x,y);
            b += grad(x,y);
        }
        ptrdiff_t c = ptrdiff_t(a / b);
        // c is biased because the numbers of black and white pixels differ
        // repeat the analysis with a symmetric window around the edge
        a = 0.0;
        b = 0.0;
        ptrdiff_t ew = desiredEdgeWidth;
        if(c-desiredEdgeWidth < 0)
            ew = c;
        if(c + ew + 1 >= w)
            ew = w - c - 1;
        for(ptrdiff_t x = c-ew; x < c+ew+1; ++x)
        {
            a += x*grad(x,y);
            b += grad(x,y);
        }
        double sc = a / b;
        sy += y;
        sx += sc;
        syy += sq(y);
        sxy += sc*y;
    }
    // fit a line to the subpixel locations
    double a = (h * sxy - sx * sy) / (h * syy - sq(sy));
    double b = (sx * syy - sxy * sy) / (h * syy - sq(sy));

    // compute the regularized subpixel values of the edge location
    angle = VIGRA_CSTD::atan(a);
    for(ptrdiff_t y = 0; y < h; ++y)
    {
        centers.push_back(a * y + b);
    }
}

template <class Image, class Vector>
void slantedEdgeCreateOversampledLine(Image const & img, Vector const & centers,
                                      Image & result)
{
    ptrdiff_t w = img.width();
    ptrdiff_t h = img.height();
    ptrdiff_t w2 = std::min(std::min(ptrdiff_t(centers[0]), ptrdiff_t(centers[h-1])),
                            std::min(ptrdiff_t(w - centers[0]) - 1,
                                     ptrdiff_t(w - centers[h-1]) - 1));
    ptrdiff_t ww = 8*w2;

    Image r(ww, 1), s(ww, 1);
    for(ptrdiff_t y = 0; y < h; ++y)
    {
        ptrdiff_t x0 = ptrdiff_t(centers[y]) - w2;
        ptrdiff_t x1 = ptrdiff_t((VIGRA_CSTD::ceil(centers[y]) - centers[y])*4);
        for(; x1 < ww; x1 += 4)
        {
            r(x1, 0) += img(x0, y);
            ++s(x1, 0);
            ++x0;
        }
    }

    for(ptrdiff_t x = 0; x < ww; ++x)
    {
        vigra_precondition(s(x,0) > 0.0,
           "slantedEdgeMTF(): Input edge is not slanted enough");
        r(x,0) /= s(x,0);
    }

    result.resize(ww-1, 1);
    for(ptrdiff_t x = 0; x < ww-1; ++x)
    {
        result(x,0) = r(x+1,0) - r(x,0);
    }
}

template <class Image, class BackInsertable>
void slantedEdgeMTFImpl(Image const & i, BackInsertable & mtf, double angle,
                        SlantedEdgeMTFOptions const & options)
{
    ptrdiff_t w = i.width();
    ptrdiff_t h = i.height();

    double slantCorrection = VIGRA_CSTD::cos(angle);
    ptrdiff_t desiredEdgeWidth = 4*options.desired_edge_width;

    Image magnitude;

    if(w - 2*desiredEdgeWidth < 64)
    {
        FFTWComplexImage otf(w, h);
        fourierTransform(srcImageRange(i), destImage(otf));

        magnitude.resize(w, h);
        for(ptrdiff_t y = 0; y < h; ++y)
        {
            for(ptrdiff_t x = 0; x < w; ++x)
            {
                magnitude(x, y) = norm(otf(x, y));
            }
        }
    }
    else
    {
        w -= 2*desiredEdgeWidth;
        FFTWComplexImage otf(w, h);
        fourierTransform(srcImageRange(i, Rect2D(Point2D(desiredEdgeWidth, 0), Size2D(w, h))),
                         destImage(otf));

        // create an image where the edge is skipped - presumably it only contains the
        // noise which can then be subtracted
        Image noise(w,h);
        copyImage(srcImageRange(i, Rect2D(Point2D(0,0), Size2D(w/2, h))),
                  destImage(noise));
        copyImage(srcImageRange(i, Rect2D(Point2D(i.width()-w/2, 0), Size2D(w/2, h))),
                  destImage(noise, Point2D(w/2, 0)));
        FFTWComplexImage fnoise(w, h);
        fourierTransform(srcImageRange(noise), destImage(fnoise));

        // subtract the noise power spectrum from the total power spectrum
        magnitude.resize(w, h);
        for(ptrdiff_t y = 0; y < h; ++y)
        {
            for(ptrdiff_t x = 0; x < w; ++x)
            {
                magnitude(x, y) = VIGRA_CSTD::sqrt(std::max(0.0, squaredNorm(otf(x, y))-squaredNorm(fnoise(x, y))));
            }
        }
    }

    Kernel1D<double> gauss;
    gauss.initGaussian(options.mtf_smoothing_scale);
    Image smooth(w,h);
    separableConvolveX(srcImageRange(magnitude), destImage(smooth), kernel1d(gauss));

    ptrdiff_t ww = w/4;
    double maxVal = smooth(0,0),
           minVal = maxVal;
    for(ptrdiff_t k = 1; k < ww; ++k)
    {
        if(smooth(k,0) >= 0.0 && smooth(k,0) < minVal)
            minVal = smooth(k,0);
    }
    double norm = maxVal-minVal;

    typedef typename BackInsertable::value_type Result;
    mtf.push_back(Result(0.0, 1.0));
    for(ptrdiff_t k = 1; k < ww; ++k)
    {
        double n = (smooth(k,0) - minVal)/norm*sq(M_PI*k/w/VIGRA_CSTD::sin(M_PI*k/w));
        double xx = 4.0*k/w/slantCorrection;
        if(n < 0.0 || xx > 1.0)
            break;
        mtf.push_back(Result(xx, n));
    }
}

} // namespace detail

/** \addtogroup SlantedEdgeMTF Camera MTF Estimation
    Determine the magnitude transfer function (MTF) of a camera using the slanted edge method.
*/
//@{

/********************************************************/
/*                                                      */
/*                     slantedEdgeMTF                   */
/*                                                      */
/********************************************************/

/** \brief Determine the magnitude transfer function of the camera.

    This operator estimates the magnitude transfer function (MTF) of a camera by means of the
    slanted edge method described in:

    ISO Standard No. 12233: <i>"Photography - Electronic still picture cameras - Resolution measurements"</i>, 2000

    The input must be an image that contains a single step edge with bright pixels on one side and dark pixels on
    the other. However, the intensity values must be neither saturated nor zero. The algorithms computes the MTF
    from the Fourier transform of the edge's derivative. Thus, if the actual MTF is anisotropic, the estimated
    MTF does actually only apply in the direction perpendicular to the edge - several edges at different
    orientations are required to estimate an anisotropic MTF.

    The algorithm returns a sequence of frequency / attenuation pairs. The frequency axis is normalized so that the
    Nyquist frequency of the original image is 0.5. Since the edge's derivative is computed with subpixel accuracy,
    the attenuation can usually be computed for frequencies significantly above the Nyquist frequency as well. The
    MTF estimate ends at either the first zero crossing of the MTF or at frequency 1, whichever comes earlier.

    The present implementation improves the original slanted edge algorithm according to ISO 12233 in a number of
    ways:

    <ul>
    <li> The edge is not required to run nearly vertically or horizontally (i.e. with a slant of approximately 5 degrees).
         The algorithm will automatically compute the edge's actual angle and adjust estimates accordingly.
         However, it is still necessary for the edge to be somewhat slanted, because subpixel-accurate estimation
         of the derivative is impossible otherwise (i.e. the edge position perpendicular to the edge direction must
         differ by at least 1 pixel between the two ends of the edge).

    <li> Our implementation uses a more accurate subpixel derivative algorithm. In addition, we first perform a shading
         correction in order to reduce possible derivative bias due to nonuniform illumination.

    <li> If the input image is large enough (i.e. there are at least 20 pixels on either side of the edge over
         the edge's entire length), our algorithm attempts to subtract the estimated noise power spectrum
         from the estimated MTF.
    </ul>

    The source value type <TT>T1</TT> must be a scalar type which is convertible to <tt>double</tt>.
    The result is written into the \a result sequence, which must be back-insertable (supports <tt>push_back()</tt>)
    and whose <tt>value_type</tt> must be constructible
    from two <tt>double</tt> values. Algorithm options can be set via the \a options object
    (see \ref vigra::NoiseNormalizationOptions for details).

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1, class BackInsertable>
        void
        slantedEdgeMTF(MultiArrayView<2, T1, S1> const & src, BackInsertable & mtf,
                       SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions());
    }
    \endcode

    \deprecatedAPI{slantedEdgeMTF}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        slantedEdgeMTF(SrcIterator sul, SrcIterator slr, SrcAccessor src, BackInsertable & mtf,
                       SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions());
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        slantedEdgeMTF(triple<SrcIterator, SrcIterator, SrcAccessor> src, BackInsertable & mtf,
                       SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions())
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/slanted_edge_mtf.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h);
    std::vector<vigra::TinyVector<double, 2> > mtf;

    ...
    // keep all options at their default values
    slantedEdgeMTF(src, mtf);

    // print the frequency / attenuation pairs found
    for(int k=0; k<result.size(); ++k)
        std::cout << "frequency: " << mtf[k][0] << ", estimated attenuation: " << mtf[k][1] << std::endl;
    \endcode

    \deprecatedUsage{slantedEdgeMTF}
    \code
    vigra::BImage src(w,h);
    std::vector<vigra::TinyVector<double, 2> > mtf;

    ...
    vigra::slantedEdgeMTF(srcImageRange(src), mtf);

    // print the frequency / attenuation pairs found
    for(int k=0; k<result.size(); ++k)
        std::cout << "frequency: " << mtf[k][0] << ", estimated attenuation: " << mtf[k][1] << std::endl;
    \endcode
    <b> Required Interface:</b>
    \code
    SrcIterator upperleft, lowerright;
    SrcAccessor src;

    typedef SrcAccessor::value_type SrcType;
    typedef NumericTraits<SrcType>::isScalar isScalar;
    assert(isScalar::asBool == true);

    double value = src(uperleft);

    BackInsertable result;
    typedef BackInsertable::value_type ResultType;
    double intensity, variance;
    result.push_back(ResultType(intensity, variance));
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void slantedEdgeMTF)

template <class SrcIterator, class SrcAccessor, class BackInsertable>
void
slantedEdgeMTF(SrcIterator sul, SrcIterator slr, SrcAccessor src, BackInsertable & mtf,
               SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions())
{
    DImage preparedInput;
    ptrdiff_t edgeWidth = detail::prepareSlantedEdgeInput(sul, slr, src, preparedInput, options);
    detail::slantedEdgeShadingCorrection(preparedInput, edgeWidth);

    ArrayVector<double> centers;
    double angle = 0.0;
    detail::slantedEdgeSubpixelShift(preparedInput, centers, angle, options);

    DImage oversampledLine;
    detail::slantedEdgeCreateOversampledLine(preparedInput, centers, oversampledLine);

    detail::slantedEdgeMTFImpl(oversampledLine, mtf, angle, options);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
inline void
slantedEdgeMTF(triple<SrcIterator, SrcIterator, SrcAccessor> src, BackInsertable & mtf,
               SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions())
{
    slantedEdgeMTF(src.first, src.second, src.third, mtf, options);
}

template <class T1, class S1, class BackInsertable>
inline void
slantedEdgeMTF(MultiArrayView<2, T1, S1> const & src, BackInsertable & mtf,
               SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions())
{
    slantedEdgeMTF(srcImageRange(src), mtf, options);
}

/********************************************************/
/*                                                      */
/*                     mtfFitGaussian                   */
/*                                                      */
/********************************************************/

/** \brief Fit a Gaussian function to a given MTF.

    This function expects a sequence of frequency / attenuation pairs as produced by \ref slantedEdgeMTF()
    and finds the best fitting Gaussian point spread function (Gaussian functions are good approximations
    of the PSF of many real cameras). It returns the standard deviation (scale) of this function. The algorithm
    computes the standard deviation by means of a linear least square on the logarithm of the MTF, i.e.
    an algebraic fit rather than a Euclidean fit - thus, the resulting Gaussian may not be the one that
    intuitively fits the data optimally.

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <class Vector>
        double mtfFitGaussian(Vector const & mtf);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/slanted_edge_mtf.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h);
    std::vector<vigra::TinyVector<double, 2> > mtf;

    ...
    slantedEdgeMTF(src, mtf);
    double scale = vigra::mtfFitGaussian(mtf)

    std::cout << "The camera PSF is approximately a Gaussian at scale " << scale << std::endl;
    \endcode

    <b> Required Interface:</b>

    \code
    Vector mtf;
    int numberOfMeasurements = mtf.size()

    double frequency = mtf[0][0];
    double attenuation = mtf[0][1];
    \endcode
*/
template <class Vector>
double mtfFitGaussian(Vector const & mtf)
{
    double minVal = mtf[0][1];
    for(ptrdiff_t k = 1; k < (ptrdiff_t)mtf.size(); ++k)
    {
        if(mtf[k][1] < minVal)
            minVal = mtf[k][1];
    }
    double x2 = 0.0,
           xy = 0.0;
    for(ptrdiff_t k = 1; k < (ptrdiff_t)mtf.size(); ++k)
    {
        if(mtf[k][1] <= 0.0)
            break;
        double x = mtf[k][0],
               y = VIGRA_CSTD::sqrt(-VIGRA_CSTD::log(mtf[k][1])/2.0)/M_PI;
        x2 += vigra::sq(x);
        xy += x*y;
        if(mtf[k][1] == minVal)
            break;
    }
    return xy / x2;
}

//@}

} // namespace vigra

#endif // VIGRA_SLANTED_EDGE_MTF_HXX
