/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_ORIENTEDTENSORFILTERS_HXX
#define VIGRA_ORIENTEDTENSORFILTERS_HXX

#include <cmath>
#include "utilities.hxx"
#include "initimage.hxx"
#include "stdconvolution.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup TensorImaging Tensor Image Processing
*/
//@{

/********************************************************/
/*                                                      */
/*                     hourGlassFilter                  */
/*                                                      */
/********************************************************/

/** \brief Anisotropic tensor smoothing with the hourglass filter.

    This function implements anisotropic tensor smoothing by an
    hourglass-shaped filters as described in
    
    U. K&ouml;the: <a href="http://hci.iwr.uni-heidelberg.de/people/ukoethe/papers/index.php#cite_structureTensor">
    <i>"Edge and Junction Detection with an Improved Structure Tensor"</i></a>, 
     in: Proc. of 25th DAGM Symposium, Magdeburg 2003, Lecture Notes in Computer Science 2781, 
     pp. 25-32, Heidelberg: Springer, 2003
     
    It is closely related to the structure tensor (see \ref structureTensor()), but
    replaces the linear tensor smoothing with a smoothing along edges only. 
    Smoothing across edges is largely suppressed. This means that the
    image structure is preserved much better because nearby features
    such as parallel edges are not blended into each other. 
    
    The hourglass filter is typically applied to a gradient tensor, i.e. the 
    Euclidean product of the gradient with itself, which can be obtained by a
    gradient operator followed with \ref vectorToTensor(), see example below. 
    The hourglass shape of the filter can be interpreted as indicating the likely 
    continuations of a local edge element. The parameter <tt>sigma</tt> determines
    the radius of the hourglass (i.e. how far the influence of the edge element 
    reaches), and <tt>rho</tt> controls its opening angle (i.e. how narrow the 
    edge orientation os followed). Recommended values are <tt>sigma = 1.4</tt>
    (or, more generally, two to three times the scale of the gradient operator
    used in the first step), and <tt>rho = 0.4</tt> which corresponds to an 
    opening angle of 22.5 degrees to either side of the edge.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        hourGlassFilter(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                        double sigma, double rho);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void hourGlassFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                             DestIterator dul, DestAccessor dest,
                             double sigma, double rho);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        inline
        void hourGlassFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                             pair<DestIterator, DestAccessor> d,
                             double sigma, double rho);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/orientedtensorfilters.hxx\>

    \code
    FImage img(w,h);
    FVector2Image gradient(w,h);
    FVector3Image tensor(w,h), smoothedTensor(w,h);
    
    gaussianGradient(srcImageRange(img), destImage(gradient), 1.0);
    vectorToTensor(srcImageRange(gradient), destImage(tensor));
    hourGlassFilter(srcImageRange(tensor), destImage(smoothedTensor), 2.0, 0.4);
    \endcode
    
    \see vectorToTensor()
*/
doxygen_overloaded_function(template <...> void hourGlassFilter)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void hourGlassFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                     DestIterator dul, DestAccessor dest,
                     double sigma, double rho)
{
    vigra_precondition(sigma >= 0.0 && rho >= 0.0,
                       "hourGlassFilter(): sigma and rho must be >= 0.0");
    vigra_precondition(src.size(sul) == 3,
                       "hourGlassFilter(): input image must have 3 bands.");
    vigra_precondition(dest.size(dul) == 3,
                       "hourGlassFilter(): output image must have 3 bands.");

    // TODO: normalization

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    double radius = VIGRA_CSTD::floor(3.0*sigma + 0.5);
    double sigma2 = -0.5 / sigma / sigma;
    double rho2 = -0.5 / rho / rho;
    double norm = 1.0 / (2.0 * M_PI * sigma * sigma);

    initImage(dul, dul+Diff2D(w,h), dest, NumericTraits<typename DestAccessor::value_type>::zero());

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator s = sul;
        DestIterator d = dul;
        for(int x=0; x<w; ++x, ++s.x, ++d.x)
        {
            double phi = 0.5 * VIGRA_CSTD::atan2(
                                     2.0*src.getComponent(s,1),
                                     (double)src.getComponent(s,0) - src.getComponent(s,2));
            double u = VIGRA_CSTD::sin(phi);
            double v = VIGRA_CSTD::cos(phi);

            double x0 = x - radius < 0 ? -x : -radius;
            double y0 = y - radius < 0 ? -y : -radius;
            double x1 = x + radius >= w ? w - x - 1 : radius;
            double y1 = y + radius >= h ? h - y - 1 : radius;

            DestIterator dwul = d + Diff2D((int)x0, (int)y0);

            for(double yy=y0; yy <= y1; ++yy, ++dwul.y)
            {
                typename DestIterator::row_iterator dw = dwul.rowIterator();
                for(double xx=x0; xx <= x1; ++xx, ++dw)
                {
                    double r2 = xx*xx + yy*yy;
                    double p  = u*xx - v*yy;
                    double q  = v*xx + u*yy;
                    double kernel = (p == 0.0) ?
                                      (q == 0.0) ?
                                       norm :
                                       0.0 :
                                       norm * VIGRA_CSTD::exp(sigma2*r2 + rho2*q*q/p/p);
                    dest.set(dest(dw) + kernel*src(s), dw);
                }
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
hourGlassFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                pair<DestIterator, DestAccessor> d,
                double sigma, double rho)
{
    hourGlassFilter(s.first, s.second, s.third, d.first, d.second, sigma, rho);
}

template <class T1, class S1,
          class T2, class S2>
inline void
hourGlassFilter(MultiArrayView<2, T1, S1> const & src,
                MultiArrayView<2, T2, S2> dest,
                double sigma, double rho)
{
    vigra_precondition(src.shape() == dest.shape(),
        "hourGlassFilter(): shape mismatch between input and output.");
    hourGlassFilter(srcImageRange(src), destImage(dest), sigma, rho);
}

/********************************************************/
/*                                                      */
/*                    ellipticGaussian                  */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void ellipticGaussian(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                      DestIterator dul, DestAccessor dest,
                      double sigmax, double sigmin)
{
    vigra_precondition(sigmax >= sigmin && sigmin >= 0.0,
                       "ellipticGaussian(): "
                       "sigmax >= sigmin and sigmin >= 0.0 required");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    double radius = VIGRA_CSTD::floor(3.0*sigmax + 0.5);
    double sigmin2 = -0.5 / sigmin / sigmin;
    double norm = 1.0 / (2.0 * M_PI * sigmin * sigmax);

    initImage(dul, dul+Diff2D(w,h), dest, NumericTraits<typename DestAccessor::value_type>::zero());

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator s = sul;
        DestIterator d = dul;
        for(int x=0; x<w; ++x, ++s.x, ++d.x)
        {
            typedef typename 
               NumericTraits<typename SrcAccessor::component_type>::RealPromote TmpType;
            TmpType d1 = src.getComponent(s,0) + src.getComponent(s,2);
            TmpType d2 = src.getComponent(s,0) - src.getComponent(s,2);
            TmpType d3 = 2.0 * src.getComponent(s,1);
            TmpType d4 = VIGRA_CSTD::sqrt(sq(d2) + sq(d3));
            TmpType excentricity = 1.0 - (d1 - d4) / (d1 + d4);
            double sigmax2 = -0.5 / sq((sigmax - sigmin)*excentricity + sigmin);
            
            double phi = 0.5 * VIGRA_CSTD::atan2(d3, d2);
            double u = VIGRA_CSTD::sin(phi);
            double v = VIGRA_CSTD::cos(phi);

            double x0 = x - radius < 0 ? -x : -radius;
            double y0 = y - radius < 0 ? -y : -radius;
            double x1 = x + radius >= w ? w - x - 1 : radius;
            double y1 = y + radius >= h ? h - y - 1 : radius;

            DestIterator dwul = d + Diff2D((int)x0, (int)y0);

            for(double yy=y0; yy <= y1; ++yy, ++dwul.y)
            {
                typename DestIterator::row_iterator dw = dwul.rowIterator();
                for(double xx=x0; xx <= x1; ++xx, ++dw)
                {
                    double p  = u*xx - v*yy;
                    double q  = v*xx + u*yy;
                    double kernel = norm * VIGRA_CSTD::exp(sigmax2*p*p + sigmin2*q*q);
                    dest.set(dest(dw) + kernel*src(s), dw);
                }
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
ellipticGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                 pair<DestIterator, DestAccessor> dest,
                 double sigmax, double sigmin)
{
    ellipticGaussian(src.first, src.second, src.third, dest.first, dest.second, sigmax, sigmin);
}

template <class T1, class S1,
          class T2, class S2>
inline void
ellipticGaussian(MultiArrayView<2, T1, S1> const & src,
                 MultiArrayView<2, T2, S2> dest,
                 double sigmax, double sigmin)
{
    vigra_precondition(src.shape() == dest.shape(),
        "ellipticGaussian(): shape mismatch between input and output.");
    ellipticGaussian(srcImageRange(src), destImage(dest), sigmax, sigmin);
}

/********************************************************/
/*                                                      */
/*         kernels for orientedTrigonometricFilter      */
/*                                                      */
/********************************************************/

class FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
    typedef VectorType::value_type    ValueType;
  
    FoerstnerKernelBase(double scale, bool ringShaped = false)
    : radius_((int)(3.0*scale+0.5)),
      weights_(2*radius_+1, 2*radius_+1),
      vectors_(2*radius_+1, 2*radius_+1)
    {
        double norm = 1.0 / (2.0 * M_PI * scale * scale);
        double s2 = -0.5 / scale / scale;
        
        for(int y = -radius_; y <= radius_; ++y)
        {
            for(int x = -radius_; x <= radius_; ++x)
            {
                double d2 = x*x + y*y;
                double d = VIGRA_CSTD::sqrt(d2);
                vectors_(x+radius_,y+radius_) = d != 0.0 ?
                                                  VectorType(x/d, -y/d) :
                                                  VectorType(ValueType(0), ValueType(0));
                weights_(x+radius_,y+radius_) = ringShaped ? 
                                       norm * d2 * VIGRA_CSTD::exp(d2 * s2) :
                                       norm * VIGRA_CSTD::exp(d2 * s2);
            }
        }
    }   
    
    ResultType operator()(int x, int y, VectorType const &) const
    {
        // isotropic filtering
        return weights_(radius_, radius_);
    }

    int radius_;
    DImage weights_;
    DVector2Image vectors_;
};

class FoerstnerRingKernelBase
: public FoerstnerKernelBase
{
  public:
    FoerstnerRingKernelBase(double scale)
    : FoerstnerKernelBase(scale, true)
    {}
};

class Cos2RingKernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Cos2RingKernel(double scale)
    : FoerstnerKernelBase(scale, true)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return d * d * weights_(x+radius_, y+radius_);
    }
};

class Cos2Kernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Cos2Kernel(double scale)
    : FoerstnerKernelBase(scale, false)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return d * d * weights_(x+radius_, y+radius_);
    }
};

class Sin2RingKernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Sin2RingKernel(double scale)
    : FoerstnerKernelBase(scale, true)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return (1.0 - d * d) * weights_(x+radius_, y+radius_);
    }
};

class Sin2Kernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Sin2Kernel(double scale)
    : FoerstnerKernelBase(scale, false)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return (1.0 - d * d) * weights_(x+radius_, y+radius_);
    }
};

class Sin6RingKernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Sin6RingKernel(double scale)
    : FoerstnerKernelBase(scale, true)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return VIGRA_CSTD::pow(1.0 - d * d, 3) * weights_(x+radius_, y+radius_);
    }
};

class Sin6Kernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Sin6Kernel(double scale)
    : FoerstnerKernelBase(scale, false)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return VIGRA_CSTD::pow(1.0 - d * d, 3) * weights_(x+radius_, y+radius_);
    }
};

class Cos6RingKernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Cos6RingKernel(double scale)
    : FoerstnerKernelBase(scale, true)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return (1.0 - VIGRA_CSTD::pow(1.0 - d * d, 3)) * weights_(x+radius_, y+radius_);
    }
};

class Cos6Kernel
: public FoerstnerKernelBase
{
  public:
    typedef double ResultType;
    typedef double WeightType;
    typedef DVector2Image::value_type VectorType;
  
    Cos6Kernel(double scale)
    : FoerstnerKernelBase(scale, false)
    {}
    
    ResultType operator()(int x, int y, VectorType const & v) const
    {
        if(x == 0 && y == 0)
            return weights_(radius_, radius_);
        double d = dot(vectors_(x+radius_, y+radius_), v);
        return (1.0 - VIGRA_CSTD::pow(1.0 - d * d, 3)) * weights_(x+radius_, y+radius_);
    }
};

/********************************************************/
/*                                                      */
/*              orientedTrigonometricFilter             */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Kernel>
void orientedTrigonometricFilter(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestAccessor dest,
                    Kernel const & kernel)
{
    vigra_precondition(src.size(sul) == 2,
                       "orientedTrigonometricFilter(): input image must have 2 bands.");
    vigra_precondition(dest.size(dul) == 3,
                       "orientedTrigonometricFilter(): output image must have 3 bands.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int radius = kernel.radius_;
    
    typedef typename SrcAccessor::value_type VectorType;
    typedef typename DestAccessor::value_type TensorType;

    initImage(dul, dul+Diff2D(w,h), dest, NumericTraits<TensorType>::zero());

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator s = sul;
        DestIterator d = dul;
        for(int x=0; x<w; ++x, ++s.x, ++d.x)
        {
            int x0 = x - radius < 0 ? -x : -radius;
            int y0 = y - radius < 0 ? -y : -radius;
            int x1 = x + radius >= w ? w - x - 1 : radius;
            int y1 = y + radius >= h ? h - y - 1 : radius;

            VectorType v(src(s));
            TensorType t(sq(v[0]), v[0]*v[1], sq(v[1]));
            double sqMag = t[0] + t[2];
            double mag = VIGRA_CSTD::sqrt(sqMag);
            if(mag != 0.0)
                v /= mag;
            else
                v *= 0.0;
            Diff2D dd;
            for(dd.y = y0; dd.y <= y1; ++dd.y)
            {
                for(dd.x = x0; dd.x <= x1; ++dd.x)
                {
                    dest.set(dest(d, dd) + kernel(dd.x, dd.y, v) * t, d, dd);
                }
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Kernel>
inline void
orientedTrigonometricFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            Kernel const & kernel)
{
    orientedTrigonometricFilter(src.first, src.second, src.third, dest.first, dest.second, kernel);
}

template <class T1, class S1,
          class T2, class S2,
          class Kernel>
inline void
orientedTrigonometricFilter(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Kernel const & kernel)
{
    vigra_precondition(src.shape() == dest.shape(),
        "orientedTrigonometricFilter(): shape mismatch between input and output.");
    orientedTrigonometricFilter(srcImageRange(src), destImage(dest), kernel);
}

//@}

} // namespace vigra

#endif /* VIGRA_ORIENTEDTENSORFILTERS_HXX */
