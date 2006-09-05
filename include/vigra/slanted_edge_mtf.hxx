/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2006 by Ullrich Koethe                  */
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
 

#ifndef VIGRA_SLANTED_EDGE_MTF_HXX
#define VIGRA_SLANTED_EDGE_MTF_HXX

#include <algorithm>
#include "vigra/utilities.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/basicgeometry.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/linear_solve.hxx"
#include "vigra/fftw3.hxx"
#include "vigra/static_assert.hxx"

namespace vigra {

class SlantedEdgeMTFOptions;

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
unsigned int 
prepareSlantedEdgeInput(SrcIterator sul, SrcIterator slr, SrcAccessor src, DestImage & res,
                        SlantedEdgeMTFOptions const & options)
{
    unsigned int w = slr.x - sul.x;
    unsigned int h = slr.y - sul.y;

    // find rough estimate of the edge
    ArrayVector<Edgel> edgels;
    cannyEdgelList(sul, slr, src, edgels, 2.0);
    std::sort(edgels.begin(), edgels.end(), SortEdgelsByStrength());
    
    double x = 0.0, y = 0.0, x2 = 0.0, y2 = 0.0, xy = 0.0;
    unsigned int c = std::min(w, h);
    
    for(unsigned int k = 0; k < c; ++k)
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
    unsigned int minimumNumberOfLines = options.minimum_number_of_lines, //20,
                 edgeWidth = options.desired_edge_width, // 10
                 minimumEdgeWidth = options.minimum_edge_width; // 5
    
    int y0, y1;
    for(; edgeWidth >= minimumEdgeWidth; --edgeWidth)
    {
        y0 = int(VIGRA_CSTD::floor((edgeWidth - xc) / slope + yc + 0.5));
        y1 = int(VIGRA_CSTD::floor((w - edgeWidth - 1 - xc) / slope + yc + 0.5));
        if(slope < 0.0)
            std::swap(y0, y1);
        if(y1 - y0 >= (int)minimumNumberOfLines)
            break;
    }
    
    vigra_precondition(edgeWidth >= minimumEdgeWidth,
        "slantedEdgeMTF(): Input edge is too slanted or image is too small");

    y0 = std::max(y0, 0);
    y1 = std::min(y1+1, (int)h);
    
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
void slantedEdgeShadingCorrection(Image & i, unsigned int edgeWidth)
{
    using namespace vigra::functor;
    
    // after prepareSlantedEdgeInput(), the white region is on the left
    // find a plane that approximates the logarithm of the white ROI
    
    transformImage(srcImageRange(i), destImage(i), log(Arg1() + Param(1.0)));
    
    unsigned int w = i.width(),
                 h = i.height(),
                 s = edgeWidth*h;
                 
    Matrix<double> m(3,3), r(3, 1), l(3, 1);
    for(unsigned int y = 0; y < h; ++y)
    {
        for(unsigned int x = 0; x < edgeWidth; ++x)
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
    for(unsigned int y = 0; y < h; ++y)
    {
        for(unsigned int x = 0; x < w; ++x)
        {
            i(x, y) = VIGRA_CSTD::exp(i(x,y) - a*x - b*y - c);
        }
    }
}

template <class Image, class BackInsertable>
void slantedEdgeSubpixelShift(Image const & img, BackInsertable & centers, double & angle,
                              SlantedEdgeMTFOptions const & options)
{
    unsigned int w = img.width();
    unsigned int h = img.height();

    Image grad(w, h);    
    Kernel1D<double> kgrad;
    kgrad.initGaussianDerivative(1.0, 1);
    
    separableConvolveX(srcImageRange(img), destImage(grad), kernel1d(kgrad));
    
    int desiredEdgeWidth = (int)options.desired_edge_width;
    double sy = 0.0, sx = 0.0, syy = 0.0, sxy = 0.0;
    for(unsigned int y = 0; y < h; ++y)
    {
        double a = 0.0, 
               b = 0.0;
        for(unsigned int x = 0; x < w; ++x)
        {
            a += x*grad(x,y);
            b += grad(x,y);
        }
        int c = int(a / b);
        // c is biased because the numbers of black and white pixels differ
        // repeat the analysis with a symmetric window around the edge
        a = 0.0;
        b = 0.0;
        int ew = desiredEdgeWidth;
        if(c-desiredEdgeWidth < 0)
            ew = c;
        if(c + ew + 1 >= (int)w)
            ew = w - c - 1;
        for(int x = c-ew; x < c+ew+1; ++x)
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
    for(unsigned int y = 0; y < h; ++y)
    {
        centers.push_back(a * y + b);
    }
}

template <class Image, class Vector>
void slantedEdgeCreateOversampledLine(Image const & img, Vector const & centers, 
                                      Image & result)
{
    unsigned int w = img.width();
    unsigned int h = img.height();
    unsigned int w2 = std::min(std::min(int(centers[0]), int(centers[h-1])), 
                               std::min(int(w - centers[0]) - 1, int(w - centers[h-1]) - 1));
    unsigned int ww = 8*w2;
    
    Image r(ww, 1), s(ww, 1);
    for(unsigned int y = 0; y < h; ++y)
    {
        int x0 = int(centers[y]) - w2;
        int x1 = int((VIGRA_CSTD::ceil(centers[y]) - centers[y])*4);
        for(; x1 < (int)ww; x1 += 4)
        {
            r(x1, 0) += img(x0, y);
            ++s(x1, 0);
            ++x0;
        }
    }

    for(unsigned int x = 0; x < ww; ++x)
    {
        vigra_precondition(s(x,0) > 0.0,
           "slantedEdgeMTF(): Input edge is not slanted enough");
        r(x,0) /= s(x,0);
    }
    
    result.resize(ww-1, 1);
    for(unsigned int x = 0; x < ww-1; ++x)
    {
        result(x,0) = r(x+1,0) - r(x,0);
    }
}

template <class Image, class BackInsertable>
void slantedEdgeMTFImpl(Image const & i, BackInsertable & mtf, double angle, 
                        SlantedEdgeMTFOptions const & options)
{
    unsigned int w = i.width();
    unsigned int h = i.height();
    
    double slantCorrection = VIGRA_CSTD::cos(angle);
    int desiredEdgeWidth = 4*options.desired_edge_width;
    
    Image magnitude;
    
    if(w - 2*desiredEdgeWidth < 64)
    {
        FFTWComplexImage otf(w, h);
        fourierTransform(srcImageRange(i), destImage(otf));
        
        magnitude.resize(w, h);
        for(unsigned int y = 0; y < h; ++y)
        {
            for(unsigned int x = 0; x < w; ++x)
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
        for(unsigned int y = 0; y < h; ++y)
        {
            for(unsigned int x = 0; x < w; ++x)
            {
                magnitude(x, y) = VIGRA_CSTD::sqrt(std::max(0.0, squaredNorm(otf(x, y))-squaredNorm(fnoise(x, y))));
            }
        }
    }
    
    Kernel1D<double> gauss;
    gauss.initGaussian(options.mtf_smoothing_scale);
    Image smooth(w,h);
    separableConvolveX(srcImageRange(magnitude), destImage(smooth), kernel1d(gauss));
    
    unsigned int ww = w/4;
    double maxVal = smooth(0,0),
           minVal = maxVal;
    for(unsigned int k = 1; k < ww; ++k)
    {
        if(smooth(k,0) >= 0.0 && smooth(k,0) < minVal)
            minVal = smooth(k,0);
    }
    double norm = maxVal-minVal;
    
    typedef typename BackInsertable::value_type Result;
    mtf.push_back(Result(0.0, 1.0));
    for(unsigned int k = 1; k < ww; ++k)
    {
        double n = (smooth(k,0) - minVal)/norm*sq(M_PI*k/w/VIGRA_CSTD::sin(M_PI*k/w));
        double xx = 4.0*k/w/slantCorrection;
        if(n < 0.0 || xx > 1.0)
            break;
        mtf.push_back(Result(xx, n));
    }
}

} // namespace detail 

class SlantedEdgeMTFOptions
{
  public:
    SlantedEdgeMTFOptions()
    : minimum_number_of_lines(20),
      desired_edge_width(10),
      minimum_edge_width(5),
      mtf_smoothing_scale(2.0)
    {}

    SlantedEdgeMTFOptions & minimumNumberOfLines(unsigned int n)
    {
        minimum_number_of_lines = n;
        return *this;
    }

    SlantedEdgeMTFOptions & desiredEdgeWidth(unsigned int n)
    {
        desired_edge_width = n;
        return *this;
    }

    SlantedEdgeMTFOptions & minimumEdgeWidth(unsigned int n)
    {
        minimum_edge_width = n;
        return *this;
    }

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

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
void
slantedEdgeMTF(SrcIterator sul, SrcIterator slr, SrcAccessor src, BackInsertable & mtf,
               SlantedEdgeMTFOptions const & options = SlantedEdgeMTFOptions())
{
    DImage preparedInput;
    unsigned int edgeWidth = detail::prepareSlantedEdgeInput(sul, slr, src, preparedInput, options);
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

} // namespace vigra

template <class Vector>
double mtfFitGaussian(Vector const & mtf)
{
    double minVal = mtf[0].second;
    for(unsigned int k = 1; k < mtf.size(); ++k)
    {
        if(mtf[k].second < minVal)
            minVal = mtf[k].second;
    }
    double x2 = 0.0,
           xy = 0.0;
    for(unsigned int k = 1; k < mtf.size(); ++k)
    {
        if(mtf[k].second <= 0.0)
            break;
        double x = mtf[k].first,
               y = VIGRA_CSTD::sqrt(-VIGRA_CSTD::log(mtf[k].second)/2.0)/M_PI;
        x2 += sq(x);
        xy += x*y;
        if(mtf[k].second == minVal)
            break;
    }
    return xy / x2;
}

#endif // VIGRA_SLANTED_EDGE_MTF_HXX
