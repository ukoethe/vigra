/************************************************************************/
/*                                                                      */
/*               Copyright 2005-2006 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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
 
#ifndef VIGRA_AFFINE_REGISTRATION_HXX
#define VIGRA_AFFINE_REGISTRATION_HXX

#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"
#include <cmath>

namespace vigra {

/** \addtogroup Registration Image Registration
*/
//@{

/********************************************************/
/*                                                      */
/*       affineMatrix2DFromCorrespondingPoints          */
/*                                                      */
/********************************************************/

/** \brief Create homogeneous matrix that maps corresponding points onto each other.
 
    For use with \ref affineWarpImage(). Since only two corresponding points are given,
    the matrix will not use a full affine transform, but only a similarity transform 
    (translation, rotation, and uniform scaling). See \
*/
template <class SrcIterator, class DestIterator>
linalg::TemporaryMatrix<double> 
affineMatrix2DFromCorrespondingPoints(SrcIterator s, SrcIterator send, DestIterator d)
{
    int size = send - s;
    
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));

    if(size == 1)
    {
        ret(0,2) = (*d)[0] - (*s)[0];
        ret(1,2) = (*d)[1] - (*s)[1];
    }
    else if(size == 2)
    {
        Matrix<double> m(4,4), r(4,1), so(4,1);
        
        for(int k=0; k<size; ++k, ++s, ++d)
        {
            m(2*k,0) = (*s)[0];
            m(2*k,1) = -(*s)[1];
            m(2*k,2) = 1.0;
            m(2*k,3) = 0.0;
            r(2*k,0) = (*d)[0];
            
            m(2*k+1,0) = (*s)[1];
            m(2*k+1,1) = (*s)[0];
            m(2*k+1,2) = 0.0;
            m(2*k+1,3) = 1.0;
            r(2*k+1,0) = (*d)[1];
        }
    
        if(!linearSolve(m, r, so))
            vigra_fail("affineMatrix2DFromCorrespondingPoints(): singular solution matrix.");
        
        ret(0,0) = so(0,0);
        ret(1,1) = so(0,0);
        ret(0,1) = -so(1,0);
        ret(1,0) = so(1,0);
        ret(0,2) = so(2,0);
        ret(1,2) = so(3,0);
    }
    else if(size >= 3)
    {
        Matrix<double> m(3,3),  rx(3,1), sx(3,1), ry(3,1), sy(3,1), c(3,1);
        c(2,0) = 1.0;
        for(int k=0; k<size; ++k, ++s, ++d)
        {
            c(0,0) = (*s)[0];
            c(1,0) = (*s)[1];
            
            m  += outer(c);
            rx += (*d)[0]*c;
            ry += (*d)[1]*c;
        }
    
        if(!linearSolve(m, rx, sx) || !linearSolve(m, ry, sy))
            vigra_fail("affineMatrix2DFromCorrespondingPoints(): singular solution matrix.");
        
        ret(0,0) = sx(0,0);
        ret(0,1) = sx(1,0);
        ret(0,2) = sx(2,0);
        ret(1,0) = sy(0,0);
        ret(1,1) = sy(1,0);
        ret(1,2) = sy(2,0);
    }

    return ret;
}

template <int SPLINEORDER = 2>
class AffineMotionEstimationOptions
{
  public:
    double burt_filter_strength;
    int highest_level, iterations_per_level;
    bool use_laplacian_pyramid;
    
    AffineMotionEstimationOptions()
    : burt_filter_strength(0.4),
      highest_level(4),
      iterations_per_level(4),
      use_laplacian_pyramid(false)
    {}
    
    template <int ORDER>
    AffineMotionEstimationOptions(AffineMotionEstimationOptions<ORDER>  const & other)
    : burt_filter_strength(other.burt_filter_strength),
      highest_level(other.highest_level),
      iterations_per_level(other.iterations_per_level),
      use_laplacian_pyramid(other.use_laplacian_pyramid)
    {}
    
    template <int NEWORDER>
    AffineMotionEstimationOptions<NEWORDER> splineOrder() const
    {
        return AffineMotionEstimationOptions<NEWORDER>(*this);
    }
    
    AffineMotionEstimationOptions & burtFilterStrength(double strength)
    {
        vigra_precondition(0.25 <= strength && strength <= 0.5,
          "AffineMotionEstimationOptions::burtFilterStrength(): strength must be between 0.25 and 0.5 (inclusive).");
        burt_filter_strength = strength;
        return *this;
    }
    
    AffineMotionEstimationOptions & highestPyramidLevel(unsigned int level)
    {
        highest_level = (int)level;
        return *this;
    }
    
    AffineMotionEstimationOptions & iterationsPerLevel(unsigned int iter)
    {
        vigra_precondition(0 < iter,
          "AffineMotionEstimationOptions::iterationsPerLevel(): must do at least one iteration per level.");
        iterations_per_level = (int)iter;
        return *this;
    }
    
    AffineMotionEstimationOptions & useGaussianPyramid(bool f = true)
    {
        use_laplacian_pyramid = !f;
        return *this;
    }
    
    AffineMotionEstimationOptions & useLaplacianPyramid(bool f = true)
    {
        use_laplacian_pyramid = f;
        return *this;
    }
};

namespace detail {

struct TranslationEstimationFunctor
{
    template <class SplineImage, class Image>
    void operator()(SplineImage const & src, Image const & dest, Matrix<double> & matrix) const
    {
        int w = dest.width();
        int h = dest.height();
        
        Matrix<double> grad(2,1), m(2,2), r(2,1), s(2,1);
        double dx = matrix(0,0), dy = matrix(1,0);
        
        for(int y = 0; y < h; ++y)
        {
            double sx = matrix(0,1)*y + matrix(0,2);
            double sy = matrix(1,1)*y + matrix(1,2);
            for(int x = 0; x < w; ++x, sx += dx, sy += dy)
            {
                if(!src.isInside(sx, sy))
                    continue;
                  
                grad(0,0) = src.dx(sx, sy);
                grad(1,0) = src.dy(sx, sy);
                double diff = dest(x, y) - src(sx, sy);
                
                m += outer(grad);
                r -= diff*grad;
            }
        }
        
        linearSolve(m, r, s);
        
        matrix(0,2) -= s(0,0);
        matrix(1,2) -= s(1,0);
    }
};

struct SimilarityTransformEstimationFunctor
{
    template <class SplineImage, class Image>
    void operator()(SplineImage const & src, Image const & dest, Matrix<double> & matrix) const
    {
        int w = dest.width();
        int h = dest.height();
        
        Matrix<double> grad(2,1), coord(4, 2), c(4, 1), m(4, 4), r(4,1), s(4,1);
        coord(0,0) = 1.0;
        coord(1,1) = 1.0;
        double dx = matrix(0,0), dy = matrix(1,0);
        
        for(int y = 0; y < h; ++y)
        {
            double sx = matrix(0,1)*y + matrix(0,2);
            double sy = matrix(1,1)*y + matrix(1,2);
            for(int x = 0; x < w; ++x, sx += dx, sy += dy)
            {
                if(!src.isInside(sx, sy))
                    continue;
                  
                grad(0,0) = src.dx(sx, sy);
                grad(1,0) = src.dy(sx, sy);
                coord(2,0) = (double)x;
                coord(3,1) = (double)x;
                coord(3,0) = -(double)y;
                coord(2,1) = (double)y;
                double diff = dest(x, y) - src(sx, sy);
                
                c = coord * grad;
                m += outer(c);
                r -= diff*c;
            }
        }
        
        linearSolve(m, r, s);
        
        matrix(0,2) -= s(0,0);
        matrix(1,2) -= s(1,0);
        matrix(0,0) -= s(2,0);
        matrix(1,1) -= s(2,0);
        matrix(0,1) += s(3,0);
        matrix(1,0) -= s(3,0);
    }
};

struct AffineTransformEstimationFunctor
{
    template <class SplineImage, class Image>
    void operator()(SplineImage const & src, Image const & dest, Matrix<double> & matrix) const
    {
        int w = dest.width();
        int h = dest.height();
        
        Matrix<double> grad(2,1), coord(6, 2), c(6, 1), m(6,6), r(6,1), s(6,1);
        coord(0,0) = 1.0;
        coord(1,1) = 1.0;
        double dx = matrix(0,0), dy = matrix(1,0);
        
        for(int y = 0; y < h; ++y)
        {
            double sx = matrix(0,1)*y + matrix(0,2);
            double sy = matrix(1,1)*y + matrix(1,2);
            for(int x = 0; x < w; ++x, sx += dx, sy += dy)
            {
                if(!src.isInside(sx, sy))
                    continue;
                  
                grad(0,0) = src.dx(sx, sy);
                grad(1,0) = src.dy(sx, sy);
                coord(2,0) = (double)x;
                coord(4,1) = (double)x;
                coord(3,0) = (double)y;
                coord(5,1) = (double)y;
                double diff = dest(x, y) - src(sx, sy);
                
                c = coord * grad;
                m += outer(c);
                r -= diff*c;
            }
        }
        
        linearSolve(m, r, s);
        
        matrix(0,2) -= s(0,0);
        matrix(1,2) -= s(1,0);
        matrix(0,0) -= s(2,0);
        matrix(0,1) -= s(3,0);
        matrix(1,0) -= s(4,0);
        matrix(1,1) -= s(5,0);
    }
};

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, 
          int SPLINEORDER, class Functor>
void 
estimateAffineMotionImpl(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                         DestIterator dul, DestIterator dlr, DestAccessor dest,
                         Matrix<double> & affineMatrix, 
                         AffineMotionEstimationOptions<SPLINEORDER> const & options,
                         Functor motionModel)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote STmpType;
    typedef BasicImage<STmpType> STmpImage;
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote DTmpType;
    typedef BasicImage<DTmpType> DTmpImage;
    
    int toplevel = options.highest_level;
    ImagePyramid<STmpImage> srcPyramid(0, toplevel, sul, slr, src);
    ImagePyramid<DTmpImage> destPyramid(0, toplevel, dul, dlr, dest);
    
    if(options.use_laplacian_pyramid)
    {
        pyramidReduceBurtLaplacian(srcPyramid, 0, toplevel, options.burt_filter_strength);
        pyramidReduceBurtLaplacian(destPyramid, 0, toplevel, options.burt_filter_strength);
    }
    else
    {
        pyramidReduceBurtFilter(srcPyramid, 0, toplevel, options.burt_filter_strength);
        pyramidReduceBurtFilter(destPyramid, 0, toplevel, options.burt_filter_strength);
    }
        
    Matrix<double> currentMatrix(affineMatrix(2,2) == 0.0 
                                    ? identityMatrix<double>(3)
                                    : affineMatrix);
    currentMatrix(0,2) /= std::pow(2.0, toplevel);
    currentMatrix(1,2) /= std::pow(2.0, toplevel);
    
    for(int level = toplevel; level >= 0; --level)
    {
        SplineImageView<SPLINEORDER, STmpType> sp(srcImageRange(srcPyramid[level]));
        
        for(int iter = 0; iter < options.iterations_per_level; ++iter)
        {
            motionModel(sp, destPyramid[level], currentMatrix);
        }
        
        if(level > 0)
        {
            currentMatrix(0,2) *= 2.0;
            currentMatrix(1,2) *= 2.0;
        }
    }
    
    affineMatrix = currentMatrix;
}

} // namespace detail 

/********************************************************/
/*                                                      */
/*                 estimateTranslation                  */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to a translation model.

    Sorry, no \ref detailedDocumentation() available yet.

    <b> Declarations:</b>

    <b>\#include</b> \<<a href="affine__registration_8hxx-source.html">vigra/affine_registration.hxx</a>\><br>
    Namespace: vigra

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateTranslation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options = 
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            triple<DestIterator, DestIterator, DestAccessor> dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options =
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void estimateTranslation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          int SPLINEORDER>
inline void
estimateTranslation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestIterator dlr, DestAccessor dest,
                    Matrix<double> & affineMatrix,
                    AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    detail::estimateAffineMotionImpl(sul, slr, src, dul, dlr, dest, affineMatrix,
                                     options, detail::TranslationEstimationFunctor());
}
             
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
estimateTranslation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestIterator dlr, DestAccessor dest,
                    Matrix<double> & affineMatrix)
{
    estimateTranslation(sul, slr, src, dul, dlr, dest,
                        affineMatrix, AffineMotionEstimationOptions<>());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, 
          int SPLINEORDER>
inline void 
estimateTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    triple<DestIterator, DestIterator, DestAccessor> dest,
                    Matrix<double> & affineMatrix, 
                    AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    estimateTranslation(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                        affineMatrix, options);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
estimateTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    triple<DestIterator, DestIterator, DestAccessor> dest,
                    Matrix<double> & affineMatrix)
{
    estimateTranslation(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                        affineMatrix, AffineMotionEstimationOptions<>());
}

/********************************************************/
/*                                                      */
/*              estimateSimilarityTransform             */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to a similarity transform model
           (e.g. translation, rotation, and uniform scaling).

    Sorry, no \ref detailedDocumentation() available yet.

    <b> Declarations:</b>

    <b>\#include</b> \<<a href="affine__registration_8hxx-source.html">vigra/affine_registration.hxx</a>\><br>
    Namespace: vigra

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateSimilarityTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options = 
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateSimilarityTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            triple<DestIterator, DestIterator, DestAccessor> dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options =
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void estimateSimilarityTransform)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          int SPLINEORDER>
inline void 
estimateSimilarityTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest,
                            Matrix<double> & affineMatrix, 
                            AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    detail::estimateAffineMotionImpl(sul, slr, src, dul, dlr, dest, affineMatrix,
                                     options, detail::SimilarityTransformEstimationFunctor());
}
             
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
estimateSimilarityTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest,
                            Matrix<double> & affineMatrix)
{
    estimateSimilarityTransform(sul, slr, src, dul, dlr, dest,
                                affineMatrix, AffineMotionEstimationOptions<>());
}
             
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, 
          int SPLINEORDER>
inline void 
estimateSimilarityTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            triple<DestIterator, DestIterator, DestAccessor> dest,
                            Matrix<double> & affineMatrix, 
                            AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    estimateSimilarityTransform(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                                affineMatrix, options);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
estimateSimilarityTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            triple<DestIterator, DestIterator, DestAccessor> dest,
                            Matrix<double> & affineMatrix)
{
    estimateSimilarityTransform(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                                affineMatrix, AffineMotionEstimationOptions<>());
}

/********************************************************/
/*                                                      */
/*                estimateAffineTransform               */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to an affine transform model
           (e.g. translation, rotation, non-uniform scaling, and shearing).

    Sorry, no \ref detailedDocumentation() available yet.

    <b> Declarations:</b>

    <b>\#include</b> \<<a href="affine__registration_8hxx-source.html">vigra/affine_registration.hxx</a>\><br>
    Namespace: vigra

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateAffineTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options = 
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  int SPLINEORDER = 2>
        void
        estimateAffineTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            triple<DestIterator, DestIterator, DestAccessor> dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options =
                                                        AffineMotionEstimationOptions<>())
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void estimateAffineTransform)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          int SPLINEORDER>
inline void 
estimateAffineTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                        DestIterator dul, DestIterator dlr, DestAccessor dest,
                        Matrix<double> & affineMatrix, 
                        AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    detail::estimateAffineMotionImpl(sul, slr, src, dul, dlr, dest, affineMatrix,
                                     options, detail::AffineTransformEstimationFunctor());
}
             
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
estimateAffineTransform(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                        DestIterator dul, DestIterator dlr, DestAccessor dest,
                        Matrix<double> & affineMatrix)
{
    estimateAffineTransform(sul, slr, src, dul, dlr, dest,
                            affineMatrix, AffineMotionEstimationOptions<>());
}
             
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, 
          int SPLINEORDER>
inline void 
estimateAffineTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        triple<DestIterator, DestIterator, DestAccessor> dest,
                        Matrix<double> & affineMatrix, 
                        AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    estimateAffineTransform(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                            affineMatrix, options);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
estimateAffineTransform(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        triple<DestIterator, DestIterator, DestAccessor> dest,
                        Matrix<double> & affineMatrix)
{
    estimateAffineTransform(src.first, src.second, src.third, dest.first, dest.second, dest.third,
                            affineMatrix, AffineMotionEstimationOptions<>());
}

//@}

} // namespace vigra


#endif /* VIGRA_AFFINE_REGISTRATION_HXX */
