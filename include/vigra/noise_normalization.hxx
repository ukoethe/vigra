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
 

#ifndef VIGRA_NOISE_NORMALIZATION_HXX
#define VIGRA_NOISE_NORMALIZATION_HXX

#include <algorithm>
#include "vigra/utilities.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/linear_solve.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/static_assert.hxx"

namespace vigra {

class NoiseNormalizationOptions;

template <class ArgumentType, class ResultType>
class NonparametricNoiseNormalizationFunctor
{
    struct Segment
    {
        double lower, a, b, shift;
    };
    
    ArrayVector<Segment> segments_;
    
    template <class T>
    double exec(unsigned int k, T t) const
    {
        if(segments_[k].a == 0.0)
        {
            return t / VIGRA_CSTD::sqrt(segments_[k].b);
        }
        else
        {
            return 2.0 / segments_[k].a * VIGRA_CSTD::sqrt(std::max(0.0, segments_[k].a * t + segments_[k].b));
        }
    }
  
  public:
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  
    template <class Vector>
    NonparametricNoiseNormalizationFunctor(Vector const & clusters)
    : segments_(clusters.size()-1)
    {
        for(unsigned int k = 0; k<segments_.size(); ++k)
        {
            segments_[k].lower = clusters[k][0];
            segments_[k].a = (clusters[k+1][1] - clusters[k][1]) / (clusters[k+1][0] - clusters[k][0]);
            segments_[k].b = clusters[k][1] - segments_[k].a * clusters[k][0];
            // FIXME: set a to zero if it is very small 
            //          - determine what 'very small' means
            //          - shouldn't the two formulas (for a == 0, a != 0) be equal in the limit a -> 0 ?
            
            if(k == 0)
            {
                segments_[k].shift = segments_[k].lower - exec(k, segments_[k].lower);
            }
            else
            {
                segments_[k].shift = exec(k-1, segments_[k].lower) - exec(k, segments_[k].lower) + segments_[k-1].shift;
            }
        }
    }
    
    result_type operator()(argument_type t) const
    {
        // find the segment
        unsigned int k = 0;
        for(; k < segments_.size(); ++k)
            if(t < segments_[k].lower)
                break;
        if(k > 0)
            --k;
        return detail::RequiresExplicitCast<ResultType>::cast(exec(k, t) + segments_[k].shift);
    }
};

template <class ArgumentType, class ResultType>
class QuadraticNoiseNormalizationFunctor
{
    double a, b, c, d, f, o;

    void init(double ia, double ib, double ic, double xmin)
    {
        a = ia;
        b = ib;
        c = ic;
        d = VIGRA_CSTD::sqrt(VIGRA_CSTD::fabs(c));
        if(c > 0.0)
        {
            o = VIGRA_CSTD::log(VIGRA_CSTD::fabs((2.0*c*xmin + b)/d + 2*VIGRA_CSTD::sqrt(c*sq(xmin) +b*xmin + a)))/d;
            f = 0.0;
        }
        else
        {
            f = VIGRA_CSTD::sqrt(b*b - 4.0*a*c);
            o = -VIGRA_CSTD::asin((2.0*c*xmin+b)/f)/d;
        }
    }
        
  public:
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  
    template <class Vector>
    QuadraticNoiseNormalizationFunctor(Vector const & clusters)
    {
        double xmin = NumericTraits<double>::max();
        Matrix<double> m(3,3), r(3, 1), l(3, 1);
        for(unsigned int k = 0; k<clusters.size(); ++k)
        {
            l(0,0) = 1.0;
            l(1,0) = clusters[k][0];
            l(2,0) = sq(clusters[k][0]);
            m += outer(l);
            r += clusters[k][1]*l;
            if(clusters[k][0] < xmin)
                xmin = clusters[k][0];
        }

        linearSolve(m, r, l);
        init(l(0,0), l(1,0), l(2,0), xmin);
    }
    
    result_type operator()(argument_type t) const
    {
        double r;
        if(c > 0.0)
            r = VIGRA_CSTD::log(VIGRA_CSTD::fabs((2.0*c*t + b)/d + 2.0*VIGRA_CSTD::sqrt(c*t*t +b*t + a)))/d-o;
        else
            r = -VIGRA_CSTD::asin((2.0*c*t+b)/f)/d-o;
        return detail::RequiresExplicitCast<ResultType>::cast(r);
    }
};

template <class ArgumentType, class ResultType>
class LinearNoiseNormalizationFunctor
{
    double a, b, o;

    void init(double ia, double ib, double xmin)
    {
        a = ia;
        b = ib;
        if(b != 0.0)
        {
            o = xmin - 2.0 / b * VIGRA_CSTD::sqrt(a + b * xmin);
        }
        else
        {
            o = xmin - xmin / VIGRA_CSTD::sqrt(a);
        }
    }
        
  public:
    typedef ArgumentType argument_type;
    typedef ResultType result_type;
  
    template <class Vector>
    LinearNoiseNormalizationFunctor(Vector const & clusters)
    {
        double xmin = NumericTraits<double>::max();
        Matrix<double> m(2,2), r(2, 1), l(2, 1);
        for(unsigned int k = 0; k<clusters.size(); ++k)
        {
            l(0,0) = 1.0;
            l(1,0) = clusters[k][0];
            m += outer(l);
            r += clusters[k][1]*l;
            if(clusters[k][0] < xmin)
                xmin = clusters[k][0];
        }

        linearSolve(m, r, l);
        init(l(0,0), l(1,0), xmin);
    }
    
    result_type operator()(argument_type t) const
    {
        double r;
        if(b != 0.0)
            r = 2.0 / b * VIGRA_CSTD::sqrt(a + b*t) + o;
        else
            r =  t / VIGRA_CSTD::sqrt(a) + o;
        return detail::RequiresExplicitCast<ResultType>::cast(r);
    }
};

#define VIGRA_NoiseNormalizationFunctor(name, type, size) \
template <class ResultType> \
class name<type, ResultType> \
{ \
    ResultType lut_[size]; \
    \
  public: \
    typedef type argument_type; \
    typedef ResultType result_type; \
     \
    template <class Vector> \
    name(Vector const & clusters) \
    { \
        name<double, ResultType> f(clusters); \
         \
        for(unsigned int k = 0; k < size; ++k) \
        { \
            lut_[k] = f(k); \
        } \
    } \
     \
    result_type operator()(argument_type t) const \
    { \
        return lut_[t]; \
    } \
};

VIGRA_NoiseNormalizationFunctor(NonparametricNoiseNormalizationFunctor, UInt8, 256)
VIGRA_NoiseNormalizationFunctor(NonparametricNoiseNormalizationFunctor, UInt16, 65536)
VIGRA_NoiseNormalizationFunctor(QuadraticNoiseNormalizationFunctor, UInt8, 256)
VIGRA_NoiseNormalizationFunctor(QuadraticNoiseNormalizationFunctor, UInt16, 65536)
VIGRA_NoiseNormalizationFunctor(LinearNoiseNormalizationFunctor, UInt8, 256)
VIGRA_NoiseNormalizationFunctor(LinearNoiseNormalizationFunctor, UInt16, 65536)

#undef VIGRA_NoiseNormalizationFunctor

namespace detail {

template <class SrcIterator, class SrcAcessor,
          class GradIterator>
bool 
iterativeNoiseEstimationChi2(SrcIterator s, SrcAcessor src, GradIterator g, 
                         double & mean, double & variance, 
                         double robustnessThreshold, int windowRadius) 
{
    double l2 = sq(robustnessThreshold);
    double countThreshold = 1.0 - VIGRA_CSTD::exp(-l2);
    double f = (1.0 - VIGRA_CSTD::exp(-l2)) / (1.0 - (1.0 + l2)*VIGRA_CSTD::exp(-l2));

    Diff2D ul(-windowRadius, -windowRadius);
    int r2 = sq(windowRadius);
    
    for(int iter=0; iter<100 ; ++iter) // maximum iteration 100 only for terminating
                                       // if something is wrong
    {
        double sum=0.0; 
        double gsum=0.0; 
        unsigned int count = 0;
        unsigned int tcount = 0;
        
        SrcIterator siy = s + ul;
        GradIterator giy = g + ul;
        for(int y=-windowRadius; y <= windowRadius; y++, ++siy.y, ++giy.y)
        {
            typename SrcIterator::row_iterator six = siy.rowIterator();   
            typename GradIterator::row_iterator gix = giy.rowIterator();   
            for(int x=-windowRadius; x <= windowRadius; x++, ++six, ++gix)
            {  
                if (sq(x) + sq(y) > r2) 
                    continue;

                ++tcount;
                if (*gix < l2*variance) 
                {
                    sum += src(six);
                    gsum += *gix;
                    ++count;
                }
            }
        }
        if (count==0) // not homogeneous enough
            return false;

        double oldvariance = variance;
        variance= f * gsum / count;
        mean = sum / count;
        
        if ( closeAtTolerance(oldvariance - variance, 0.0, 1e-10)) 
            return (count >= tcount * countThreshold / 2.0); // sufficiently many valid points
    }
    return false; // no convergence
}

template <class SrcIterator, class SrcAcessor,
          class GradIterator>
bool 
iterativeNoiseEstimationGauss(SrcIterator s, SrcAcessor src, GradIterator, 
                         double & mean, double & variance, 
                         double robustnessThreshold, int windowRadius) 
{
    double l2 = sq(robustnessThreshold);
    double countThreshold = erf(VIGRA_CSTD::sqrt(0.5 * l2));
    double f = countThreshold / (countThreshold - VIGRA_CSTD::sqrt(2.0/M_PI*l2)*VIGRA_CSTD::exp(-l2/2.0));

    mean = src(s);

    Diff2D ul(-windowRadius, -windowRadius);
    int r2 = sq(windowRadius);
    
    for(int iter=0; iter<100 ; ++iter) // maximum iteration 100 only for terminating
                                       // if something is wrong
    {
        double sum = 0.0; 
        double sum2 = 0.0; 
        unsigned int count = 0;
        unsigned int tcount = 0;
        
        SrcIterator siy = s + ul;
        for(int y=-windowRadius; y <= windowRadius; y++, ++siy.y)
        {
            typename SrcIterator::row_iterator six = siy.rowIterator();   
            for(int x=-windowRadius; x <= windowRadius; x++, ++six)
            {  
                if (sq(x) + sq(y) > r2) 
                    continue;

                ++tcount;
                if (sq(src(six) - mean) < l2*variance) 
                {
                    sum += src(six);
                    sum2 += sq(src(six));
                    ++count;
                }
            }
        }
        if (count==0) // not homogeneous enough
            return false;

        double oldmean = mean;
        double oldvariance = variance;
        mean = sum / count;
        variance= f * (sum2 / count - sq(mean));
        
        if ( closeAtTolerance(oldmean - mean, 0.0, 1e-10) &&
             closeAtTolerance(oldvariance - variance, 0.0, 1e-10)) 
            return (count >= tcount * countThreshold / 2.0); // sufficiently many valid points
    }
    return false; // no convergence
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
symmetricDifferenceSquaredMagnitude(
     SrcIterator sul, SrcIterator slr, SrcAccessor src,
     DestIterator dul, DestAccessor dest)
{
    using namespace vigra::functor;
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TmpType> TmpImage;

    Kernel1D<double> mask;
    mask.initSymmetricGradient();
    mask.setBorderTreatment(BORDER_TREATMENT_REFLECT);

    TmpImage dx(w, h), dy(w, h);
    separableConvolveX(srcIterRange(sul, slr, src), destImage(dx),  kernel1d(mask));
    separableConvolveY(srcIterRange(sul, slr, src), destImage(dy),  kernel1d(mask));
    combineTwoImages(srcImageRange(dx), srcImage(dy), destIter(dul, dest), Arg1()*Arg1() + Arg2()*Arg2());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
findHomogeneousRegionsFoerstner(
     SrcIterator sul, SrcIterator slr, SrcAccessor src,
     DestIterator dul, DestAccessor dest,
     unsigned int windowRadius = 6, double homogeneityThreshold = 40.0)
{
    using namespace vigra::functor;
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TmpType> TmpImage;

    BImage btmp(w, h);
    transformImage(srcIterRange(sul, slr, src), destImage(btmp), 
                    ifThenElse(Arg1() <= Param(homogeneityThreshold), Param(1), Param(0)));

    // Erosion
    discErosion(srcImageRange(btmp), destIter(dul, dest), windowRadius);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
findHomogeneousRegions(
     SrcIterator sul, SrcIterator slr, SrcAccessor src,
     DestIterator dul, DestAccessor dest)
{
    localMinima(sul, slr, src, dul, dest);
}

template <class Vector1, class Vector2>
void noiseVarianceListMedianCut(Vector1 const & noise, Vector2 & clusters,
                                unsigned int maxClusterCount)
{
    typedef typename Vector2::value_type Result;
    
    clusters.push_back(Result(0, noise.size()));
    
    while(clusters.size() <= maxClusterCount)
    {
        // find biggest cluster
        unsigned int kMax;
        double diffMax = 0.0;
        for(unsigned int k=0; k < clusters.size(); ++k)
        {
            double diff = noise[clusters[k][1]-1][0] - noise[clusters[k][0]][0];
            if(diff > diffMax)
            {
                diffMax = diff;
                kMax = k;
            }
        }
        
        if(diffMax == 0.0)
            return; // all clusters have only one value

        unsigned int k1 = clusters[kMax][0],
                     k2 = clusters[kMax][1];
        unsigned int kSplit = k1 + (k2 - k1) / 2;
        clusters[kMax][1] = kSplit;
        clusters.push_back(Result(kSplit, k2));
    }
}

struct SortNoiseByMean
{
    template <class T>
    bool operator()(T const & l, T const & r) const
    {
        return l[0] < r[0];
    }
};

struct SortNoiseByVariance
{
    template <class T>
    bool operator()(T const & l, T const & r) const
    {
        return l[1] < r[1];
    }
};

template <class Vector1, class Vector2, class Vector3>
void noiseVarianceClusterAveraging(Vector1 & noise, Vector2 & clusters,
                                   Vector3 & result, double quantile)
{
    typedef typename Vector1::iterator Iter;
    typedef typename Vector3::value_type Result;
    
    for(unsigned int k=0; k<clusters.size(); ++k)
    {
        Iter i1 = noise.begin() + clusters[k][0];
        Iter i2 = noise.begin() + clusters[k][1];
        
        std::sort(i1, i2, SortNoiseByVariance());
        
        std::size_t size = std::max<std::size_t>(1, 
                                    std::min<std::size_t>(i2 - i1, 
                                                          static_cast<std::size_t>(VIGRA_CSTD::ceil(quantile*(i2 - i1)))));
        i2 = i1 + size;
        
        double mean = 0.0,
               variance = 0.0;
        for(; i1 < i2; ++i1)
        {
            mean += (*i1)[0];
            variance += (*i1)[1];
        }
        
        result.push_back(Result(mean / size, variance / size));
    }
}

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
void noiseVarianceEstimationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                           BackInsertable & result,
                           NoiseNormalizationOptions const & options)
{
    typedef typename BackInsertable::value_type ResultType;
    
    unsigned int w = slr.x - sul.x;
    unsigned int h = slr.y - sul.y;

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TmpType> TmpImage;

    TmpImage gradient(w, h);
    symmetricDifferenceSquaredMagnitude(sul, slr, src, gradient.upperLeft(), gradient.accessor());

    BImage homogeneous(w, h);
    findHomogeneousRegions(gradient.upperLeft(), gradient.lowerRight(), gradient.accessor(),
                                   homogeneous.upperLeft(), homogeneous.accessor());
                                       
    // Generate noise of each of the remaining pixels == centers of homogenous areas (border is not used)
    unsigned int windowRadius = options.window_radius;
    for(unsigned int y=windowRadius; y<h-windowRadius; ++y)
    {
        for(unsigned int x=windowRadius; x<w-windowRadius; ++x)
        {
            if (! homogeneous(x, y))
                continue;

            Diff2D center(x, y);
            double mean = 0.0, variance = options.noise_variance_initial_guess;
            
            bool success;
            
            if(options.use_gradient)
            {
                success = iterativeNoiseEstimationChi2(sul + center, src,
                              gradient.upperLeft() + center, mean, variance, 
                              options.noise_estimation_quantile, windowRadius);
            }
            else
            {
                success = iterativeNoiseEstimationGauss(sul + center, src,
                              gradient.upperLeft() + center, mean, variance, 
                              options.noise_estimation_quantile, windowRadius);
            }
            if (success)
            {
                result.push_back(ResultType(mean, variance));
            }
        }
    }
}

template <class Vector, class BackInsertable>  
void noiseVarianceClusteringImpl(Vector & noise, BackInsertable & result,
                           unsigned int clusterCount, double quantile)
{
    std::sort(noise.begin(), noise.end(), detail::SortNoiseByMean());
    
    ArrayVector<TinyVector<unsigned int, 2> > clusters;
    detail::noiseVarianceListMedianCut(noise, clusters, clusterCount);

    std::sort(clusters.begin(), clusters.end(), detail::SortNoiseByMean());

    detail::noiseVarianceClusterAveraging(noise, clusters, result, quantile);
}

template <class Functor,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
noiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                       DestIterator dul, DestAccessor dest,
                       NoiseNormalizationOptions const & options)
{
    ArrayVector<TinyVector<double, 2> > noiseData;
    noiseVarianceEstimationImpl(sul, slr, src, noiseData, options);
    
    if(noiseData.size() < 10)
        return false;
        
    std::sort(noiseData.begin(), noiseData.end(), SortNoiseByMean());

    ArrayVector<TinyVector<double, 2> > noiseClusters;

    noiseVarianceClusteringImpl(noiseData, noiseClusters, 
                                  options.cluster_count, options.averaging_quantile);
    
    transformImage(sul, slr, src, dul, dest, Functor(noiseClusters));
    
    return true;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
nonparametricNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                    DestIterator dul, DestAccessor dest,
                                    NoiseNormalizationOptions const & options,
                                    VigraTrueType /* isScalar */)
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    return noiseNormalizationImpl<NonparametricNoiseNormalizationFunctor<SrcType, DestType> >
                                                         (sul, slr, src, dul, dest, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
nonparametricNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              NoiseNormalizationOptions const & options,
                              VigraFalseType /* isScalar */)
{
    int bands = src.size(sul);
    for(int b=0; b<bands; ++b)
    {
        VectorElementAccessor<SrcAccessor> sband(b, src);
        VectorElementAccessor<DestAccessor> dband(b, dest);
        typedef typename VectorElementAccessor<SrcAccessor>::value_type SrcType;
        typedef typename VectorElementAccessor<DestAccessor>::value_type DestType;

        if(!noiseNormalizationImpl<NonparametricNoiseNormalizationFunctor<SrcType, DestType> >
                                                           (sul, slr, sband, dul, dband, options))
            return false;
    }
    return true;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
quadraticNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                    DestIterator dul, DestAccessor dest,
                                    NoiseNormalizationOptions const & options,
                                    VigraTrueType /* isScalar */)
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    return noiseNormalizationImpl<QuadraticNoiseNormalizationFunctor<SrcType, DestType> >
                                                         (sul, slr, src, dul, dest, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
quadraticNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              NoiseNormalizationOptions const & options,
                              VigraFalseType /* isScalar */)
{
    int bands = src.size(sul);
    for(int b=0; b<bands; ++b)
    {
        VectorElementAccessor<SrcAccessor> sband(b, src);
        VectorElementAccessor<DestAccessor> dband(b, dest);
        typedef typename VectorElementAccessor<SrcAccessor>::value_type SrcType;
        typedef typename VectorElementAccessor<DestAccessor>::value_type DestType;

        if(!noiseNormalizationImpl<QuadraticNoiseNormalizationFunctor<SrcType, DestType> >
                                                           (sul, slr, sband, dul, dband, options))
            return false;
    }
    return true;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
void 
quadraticNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              double a0, double a1, double a2,
                              VigraTrueType /* isScalar */)
{
    ArrayVector<TinyVector<double, 2> > noiseClusters;
    noiseClusters.push_back(TinyVector<double, 2>(0.0, a0));
    noiseClusters.push_back(TinyVector<double, 2>(1.0, a0 + a1 + a2));
    noiseClusters.push_back(TinyVector<double, 2>(2.0, a0 + 2.0*a1 + 4.0*a2));
    transformImage(sul, slr, src, dul, dest, 
                   QuadraticNoiseNormalizationFunctor<typename SrcAccessor::value_type,
                                                   typename DestAccessor::value_type>(noiseClusters));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
void 
quadraticNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              double a0, double a1, double a2,
                              VigraFalseType /* isScalar */)
{
    int bands = src.size(sul);
    for(int b=0; b<bands; ++b)
    {
        VectorElementAccessor<SrcAccessor> sband(b, src);
        VectorElementAccessor<DestAccessor> dband(b, dest);
        quadraticNoiseNormalizationImpl(sul, slr, sband, dul, dband, a0, a1, a2, VigraTrueType());
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
linearNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                    DestIterator dul, DestAccessor dest,
                                    NoiseNormalizationOptions const & options,
                                    VigraTrueType /* isScalar */)
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    return noiseNormalizationImpl<LinearNoiseNormalizationFunctor<SrcType, DestType> >
                                                         (sul, slr, src, dul, dest, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
bool 
linearNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              NoiseNormalizationOptions const & options,
                              VigraFalseType /* isScalar */)
{
    int bands = src.size(sul);
    for(int b=0; b<bands; ++b)
    {
        VectorElementAccessor<SrcAccessor> sband(b, src);
        VectorElementAccessor<DestAccessor> dband(b, dest);
        typedef typename VectorElementAccessor<SrcAccessor>::value_type SrcType;
        typedef typename VectorElementAccessor<DestAccessor>::value_type DestType;

        if(!noiseNormalizationImpl<LinearNoiseNormalizationFunctor<SrcType, DestType> >
                                                           (sul, slr, sband, dul, dband, options))
            return false;
    }
    return true;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
void 
linearNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              double a0, double a1,
                              VigraTrueType /* isScalar */)
{
    ArrayVector<TinyVector<double, 2> > noiseClusters;
    noiseClusters.push_back(TinyVector<double, 2>(0.0, a0));
    noiseClusters.push_back(TinyVector<double, 2>(1.0, a0 + a1));
    transformImage(sul, slr, src, dul, dest, 
                   LinearNoiseNormalizationFunctor<typename SrcAccessor::value_type,
                                                   typename DestAccessor::value_type>(noiseClusters));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
void 
linearNoiseNormalizationImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              double a0, double a1,
                              VigraFalseType /* isScalar */)
{
    int bands = src.size(sul);
    for(int b=0; b<bands; ++b)
    {
        VectorElementAccessor<SrcAccessor> sband(b, src);
        VectorElementAccessor<DestAccessor> dband(b, dest);
        linearNoiseNormalizationImpl(sul, slr, sband, dul, dband, a0, a1, VigraTrueType());
    }
}

} // namespace detail 

class NoiseNormalizationOptions
{
  public:
    NoiseNormalizationOptions()
    : window_radius(6),
      cluster_count(10),
      noise_estimation_quantile(1.5),
      averaging_quantile(0.8),
      noise_variance_initial_guess(10.0),
      use_gradient(true)
    {}
    
    NoiseNormalizationOptions & useGradient(bool r)
    {
        use_gradient = r;
        return *this;
    }
    
    NoiseNormalizationOptions & windowRadius(unsigned int r)
    {
        vigra_precondition(r > 0,
            "NoiseNormalizationOptions: window radius must be > 0.");
        window_radius = r;
        return *this;
    }
    
    NoiseNormalizationOptions & clusterCount(unsigned int c)
    {
        vigra_precondition(c > 0,
            "NoiseNormalizationOptions: cluster count must be > 0.");
        cluster_count = c;
        return *this;
    }
    
    NoiseNormalizationOptions & averagingQuantile(double quantile)
    {
        vigra_precondition(quantile > 0.0 && quantile <= 1.0,
            "NoiseNormalizationOptions: averaging quantile must be between 0 and 1.");
        averaging_quantile = quantile;
        return *this;
    }
    
    NoiseNormalizationOptions & noiseEstimationQuantile(double quantile)
    {
        vigra_precondition(quantile > 0.0,
            "NoiseNormalizationOptions: noise estimation quantile must be > 0.");
        noise_estimation_quantile = quantile;
        return *this;
    }
    
    NoiseNormalizationOptions & noiseVarianceInitialGuess(double guess)
    {
        vigra_precondition(guess > 0.0,
            "NoiseNormalizationOptions: noise variance initial guess must be > 0.");
        noise_variance_initial_guess = guess;
        return *this;
    }
    
    unsigned int window_radius, cluster_count;
    double noise_estimation_quantile, averaging_quantile, noise_variance_initial_guess;
    bool use_gradient;
};

template <bool P>
struct noiseVarianceEstimation_can_only_work_on_scalar_images
: vigra::staticAssert::AssertBool<P>
{};

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
inline
void noiseVarianceEstimation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                           BackInsertable & result,
                           NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    typedef typename BackInsertable::value_type ResultType;
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename NumericTraits<SrcType>::isScalar isScalar;
    
    VIGRA_STATIC_ASSERT((
        noiseVarianceEstimation_can_only_work_on_scalar_images<(isScalar::asBool)>));
    
    detail::noiseVarianceEstimationImpl(sul, slr, src, result, options);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
inline 
void noiseVarianceEstimation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           BackInsertable & result,
                           NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    noiseVarianceEstimation(src.first, src.second, src.third, result, options);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
inline
void noiseVarianceClustering(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                           BackInsertable & result,
                           NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    ArrayVector<TinyVector<double, 2> > variance;
    noiseVarianceEstimation(sul, slr, src, variance, options);
    detail::noiseVarianceClusteringImpl(variance, result, options.cluster_count, options.averaging_quantile);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable> 
inline 
void noiseVarianceClustering(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           BackInsertable & result,
                           NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    noiseVarianceClustering(src.first, src.second, src.third, result, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
inline bool 
nonparametricNoiseNormalization(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                DestIterator dul, DestAccessor dest,
                                NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    typedef typename SrcAccessor::value_type SrcType;
    
    return detail::nonparametricNoiseNormalizationImpl(sul, slr, src, dul, dest, options,
                                         typename NumericTraits<SrcType>::isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
bool nonparametricNoiseNormalization(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<DestIterator, DestAccessor> dest,
                                     NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    return nonparametricNoiseNormalization(src.first, src.second, src.third, dest.first, dest.second, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
inline bool 
quadraticNoiseNormalization(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                DestIterator dul, DestAccessor dest,
                                NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    typedef typename SrcAccessor::value_type SrcType;
    
    return detail::quadraticNoiseNormalizationImpl(sul, slr, src, dul, dest, options,
                                         typename NumericTraits<SrcType>::isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
bool quadraticNoiseNormalization(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<DestIterator, DestAccessor> dest,
                                     NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    return quadraticNoiseNormalization(src.first, src.second, src.third, dest.first, dest.second, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
inline void 
quadraticNoiseNormalization(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestAccessor dest,
                            double a0, double a1, double a2)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::quadraticNoiseNormalizationImpl(sul, slr, src, dul, dest, a0, a1, a2,
                                         typename NumericTraits<SrcType>::isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void quadraticNoiseNormalization(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 double a0, double a1, double a2)
{
    quadraticNoiseNormalization(src.first, src.second, src.third, dest.first, dest.second, a0, a1, a2);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline bool 
linearNoiseNormalization(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                DestIterator dul, DestAccessor dest,
                                NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    typedef typename SrcAccessor::value_type SrcType;
    
    return detail::linearNoiseNormalizationImpl(sul, slr, src, dul, dest, options,
                                         typename NumericTraits<SrcType>::isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
bool linearNoiseNormalization(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                     pair<DestIterator, DestAccessor> dest,
                                     NoiseNormalizationOptions const & options = NoiseNormalizationOptions())
{
    return linearNoiseNormalization(src.first, src.second, src.third, dest.first, dest.second, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor> 
inline
void linearNoiseNormalization(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                              DestIterator dul, DestAccessor dest,
                              double a0, double a1)
{
    typedef typename SrcAccessor::value_type SrcType;
    
    detail::linearNoiseNormalizationImpl(sul, slr, src, dul, dest, a0, a1, 
                                         typename NumericTraits<SrcType>::isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void linearNoiseNormalization(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                              pair<DestIterator, DestAccessor> dest,
                              double a0, double a1)
{
    linearNoiseNormalization(src.first, src.second, src.third, dest.first, dest.second, a0, a1);
}



} // namespace vigra

#endif // VIGRA_NOISE_NORMALIZATION_HXX
