/************************************************************************/
/*                                                                      */
/*                 Copyright 2005-2006 by Ullrich Koethe                */
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

#ifndef VIGRA_AFFINE_REGISTRATION_HXX
#define VIGRA_AFFINE_REGISTRATION_HXX

#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"
#include "imagecontainer.hxx"
#include "multi_shape.hxx"
#include "affinegeometry.hxx"

#include <cmath>

namespace vigra {

/** \addtogroup Registration Image Registration

    Transform different image into a common coordinate system.
*/
//@{

/********************************************************/
/*                                                      */
/*         affineMatrix2DFromCorrespondingPoints        */
/*                                                      */
/********************************************************/

/** \brief Create homogeneous matrix that maps corresponding points onto each other.

    For use with \ref affineWarpImage(). When only two corresponding points are given,
    the matrix will only represent a similarity transform
    (translation, rotation, and uniform scaling). When only one point pair is given,
    the result will be a pure translation.
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
        Matrix<double> m(3,3),    rx(3,1), sx(3,1), ry(3,1), sy(3,1), c(3,1);
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

    /** \brief Option object for affine registration functions.

        The template parameter <tt>SPLINEORDER</tt> (default: 2) specifies
        the order of interpolation for the intensities at non-integer image
        coordinates.
    */
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
    AffineMotionEstimationOptions(AffineMotionEstimationOptions<ORDER>    const & other)
    : burt_filter_strength(other.burt_filter_strength),
      highest_level(other.highest_level),
      iterations_per_level(other.iterations_per_level),
      use_laplacian_pyramid(other.use_laplacian_pyramid)
    {}

        /** \brief Change the spline order for intensity interpolation.

            Usage:
            \code
            // use linear interpolation
            AffineMotionEstimationOptions<>().splineOrder<1>();
            \endcode

            Default: order = 2 (quadratic interpolation)
        */
    template <int NEWORDER>
    AffineMotionEstimationOptions<NEWORDER> splineOrder() const
    {
        return AffineMotionEstimationOptions<NEWORDER>(*this);
    }

        /** \brief Define amount of smoothing before subsampling to the next pyramid level.

            Pyramids are created with the Burt filter:
            \code
            [0.25 - center / 2.0, 0.25, center, 0.25, 0.25 - center / 2.0]
            \endcode
            \a center must thus be between 0.25 and 0.5, and the smaller it is,
            the more smoothing is applied.

            Default: 0.4 (moderate smoothing)
        */
    AffineMotionEstimationOptions & burtFilterCenterStrength(double center)
    {
        vigra_precondition(0.25 <= center && center <= 0.5,
          "AffineMotionEstimationOptions::burtFilterCenterStrength(): center must be between 0.25 and 0.5 (inclusive).");
        burt_filter_strength = center;
        return *this;
    }

        /** \brief Define the highest level of the image pyramid.

            The original image is at level 0, and each level downsamples
            the image by 1/2.

            Default: 4 (16-fold downsampling)
        */
    AffineMotionEstimationOptions & highestPyramidLevel(unsigned int level)
    {
        highest_level = (int)level;
        return *this;
    }

        /** \brief Number of Gauss-Newton iterations per level.

            Default: 4
        */
    AffineMotionEstimationOptions & iterationsPerLevel(unsigned int iter)
    {
        vigra_precondition(0 < iter,
          "AffineMotionEstimationOptions::iterationsPerLevel(): must do at least one iteration per level.");
        iterations_per_level = (int)iter;
        return *this;
    }

        /** \brief Base registration on a Gaussian pyramid.

            Images are registered such that the similarity in intensity is
            maximized.

            Default: true
        */
    AffineMotionEstimationOptions & useGaussianPyramid(bool f = true)
    {
        use_laplacian_pyramid = !f;
        return *this;
    }

        /** \brief Base registration on a Laplacian pyramid.

            Images are registered such that the similarity in second
            derivatives (=Laplacian operator results) is maximized.

            Default: false
        */
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
/*                   estimateTranslation                */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to a translation model.

    This function applies the same algorithm as \ref estimateAffineTransform()
    with the additional constraint that the motion model must be a translation
    rather than affine.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/affine_registration.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  int SPLINEORDER>
        void
        estimateTranslation(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options = AffineMotionEstimationOptions<SPLINEORDER>());
    }
    \endcode

    \deprecatedAPI{estimateTranslation}
    pass \ref ImageIterators and \ref DataAccessors :
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
    \deprecatedEnd
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

template <class T1, class S1,
          class T2, class S2,
          int SPLINEORDER>
inline void
estimateTranslation(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    Matrix<double> & affineMatrix,
                    AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    estimateTranslation(srcImageRange(src), destImageRange(dest),
                        affineMatrix, options);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateTranslation(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    Matrix<double> & affineMatrix)
{
    estimateTranslation(srcImageRange(src), destImageRange(dest),
                        affineMatrix, AffineMotionEstimationOptions<>());
}

/********************************************************/
/*                                                      */
/*                estimateSimilarityTransform           */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to a similarity transform model
           (e.g. translation, rotation, and uniform scaling).

    This function applies the same algorithm as \ref estimateAffineTransform()
    with the additional constraint that the motion model must be a similarity
    transform rather than affine.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/affine_registration.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  int SPLINEORDER>
        void
        estimateSimilarityTransform(MultiArrayView<2, T1, S1> const & src,
                                    MultiArrayView<2, T2, S2> dest,
                                    Matrix<double> & affineMatrix,
                                    AffineMotionEstimationOptions<SPLINEORDER> const & options =
                                                        AffineMotionEstimationOptions<SPLINEORDER>());
    }
    \endcode

    \deprecatedAPI{estimateSimilarityTransform}
    pass \ref ImageIterators and \ref DataAccessors :
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
    \deprecatedEnd
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

template <class T1, class S1,
          class T2, class S2,
          int SPLINEORDER>
inline void
estimateSimilarityTransform(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Matrix<double> & affineMatrix,
                            AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    estimateSimilarityTransform(srcImageRange(src), destImageRange(dest),
                                affineMatrix, options);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateSimilarityTransform(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Matrix<double> & affineMatrix)
{
    estimateSimilarityTransform(srcImageRange(src), destImageRange(dest),
                                affineMatrix, AffineMotionEstimationOptions<>());
}

/********************************************************/
/*                                                      */
/*                  estimateAffineTransform             */
/*                                                      */
/********************************************************/

/** \brief Estimate the optical flow between two images according to an affine transform model
           (e.g. translation, rotation, non-uniform scaling, and shearing).

    This function implements the algorithm described in

    J.R. Bergen, P. Anandan, K.J. Hanna, R. Hingorani: <i>"Hierarchical model-based motion estimation"</i>, ECCV 1992

    Specifically, it minimizes the squared loss between the images I at two consecutive time
    points t-1 and t:
    \f[ \min_{\theta} \sum_{\mathbf{x}} \left(I(\mathbf{x}, t) - I(\mathbf{x} - \mathbf{u}_{\theta}(\mathbf{x}), t-1)\right)^2
    \f]
    where \f$\mathbf{x}\f$ are the pixel coordinates and \f$\mathbf{u}_{\theta}(\mathbf{x})\f$
    is an affine motion model parameterized by \f$\theta\f$. Since the objective is
    non-linear, it is linearized by first-order Taylor expansion w.r.t. \f$\theta\f$,
    and a local optimum is determined iteratively by the Gauss-Newton method. To handle
    larger displacements, the algorithm employs a coarse-to-fine strategy, where the
    motion is first estimated on downsampled versions of the images and then refined at
    consecutively higher resolutions.

    The algorithm's parameters can be controlled by the option object
    \ref vigra::AffineMotionEstimationOptions. In particular, one can determine if
    <ul>
    <li> I in the objective refers to the original images (default) or their second
    derivatives (<tt>options.useLaplacianPyramid()</tt> -- makes motion
    estimation invariant against additive intensity offsets);</li>
    <li> the highest pyramid level to be used
    (<tt>options.highestPyramidLevel(h)</tt> -- images are downsampled to 2<sup>-h</sup> times their original size, default: h=4);</li>
    <li> the number of Gauss-Newton iterations per resolution level
    (<tt>options.iterationsPerLevel(i)</tt>, default: i=4);</li>
    <li> the interpolation order to compute subpixel intensities \f$I(\mathbf{x} - \mathbf{u}_{\theta}(\mathbf{x}), t-1)\f$ (default: 2)
    </ul>

    The resulting affine model is stored in parameter <tt>affineMatrix</tt>, which
    can be used by \ref affineWarpImage() to apply the transformation to time frame t-1.
    See documentation there for the precise meaning of the matrix elements.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/affine_registration.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  int SPLINEORDER>
        void
        estimateAffineTransform(MultiArrayView<2, T1, S1> const & src,
                                MultiArrayView<2, T2, S2> dest,
                                Matrix<double> & affineMatrix,
                                AffineMotionEstimationOptions<SPLINEORDER> const & options =
                                                   AffineMotionEstimationOptions<SPLINEORDER>());
    }
    \endcode

    \deprecatedAPI{estimateAffineTransform}
    pass \ref ImageIterators and \ref DataAccessors :
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
    \deprecatedEnd
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

template <class T1, class S1,
          class T2, class S2,
          int SPLINEORDER>
inline void
estimateAffineTransform(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                        Matrix<double> & affineMatrix,
                        AffineMotionEstimationOptions<SPLINEORDER> const & options)
{
    vigra_precondition(src.shape() == dest.shape(),
        "estimateAffineTransform(): shape mismatch between input and output.");
    estimateAffineTransform(srcImageRange(src), destImageRange(dest),
                            affineMatrix, options);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateAffineTransform(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                        Matrix<double> & affineMatrix)
{
    vigra_precondition(src.shape() == dest.shape(),
        "estimateAffineTransform(): shape mismatch between input and output.");
    estimateAffineTransform(srcImageRange(src), destImageRange(dest),
                            affineMatrix, AffineMotionEstimationOptions<>());
}

//@}

} // namespace vigra


#endif /* VIGRA_AFFINE_REGISTRATION_HXX */
