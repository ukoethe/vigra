/************************************************************************/
/*                                                                      */
/*               Copyright 2005-2006 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_AFFINEGEOMETRY_HXX
#define VIGRA_AFFINEGEOMETRY_HXX

#include "mathutil.hxx"
#include "matrix.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"
#include <cmath>

namespace vigra {

/** \addtogroup GeometricTransformations Geometric Transformations
*/
//@{

/********************************************************/
/*                                                      */
/*                create affine matrices                */
/*                                                      */
/********************************************************/

/** \brief Create homogeneous matrix representing a 2D translation.
 
    For use with \ref affineWarpImage().
*/
inline
linalg::TemporaryMatrix<double> translationMatrix2D(TinyVector<double, 2> const & shift)
{
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));
    ret(0,2) = shift[0];
    ret(1,2) = shift[1];
    return ret;
}

/** \brief Create homogeneous matrix representing a 2D uniform scaling about the coordinate origin.
 
    For use with \ref affineWarpImage().
*/
inline
linalg::TemporaryMatrix<double> scalingMatrix2D(double scalingFactor)
{
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));
    ret(0,0) = scalingFactor;
    ret(1,1) = scalingFactor;
    return ret;
}

/** \brief Create homogeneous matrix representing a 2D non-uniform scaling about the coordinate origin.
 
    For use with \ref affineWarpImage().
*/
inline
linalg::TemporaryMatrix<double> scalingMatrix2D(double sx, double sy)
{
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));
    ret(0,0) = sx;
    ret(1,1) = sy;
    return ret;
}

/** \brief Create homogeneous matrix representing a 2D shearing.
 
    For use with \ref affineWarpImage().
*/
inline
linalg::TemporaryMatrix<double> shearMatrix2D(double s01, double s10)
{
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));
    ret(0,1) = s01;
    ret(1,0) = s10;
    return ret;
}

/** \brief Create homogeneous matrix representing a 2D rotation about the coordinate origin.
 
    For use with \ref affineWarpImage(). Angle must be in radians.
*/
inline
linalg::TemporaryMatrix<double> rotationMatrix2DRadians(double angle)
{
    linalg::TemporaryMatrix<double> ret(identityMatrix<double>(3));
    double s = std::sin(angle);
    double c = std::cos(angle);
    ret(0,0) = c;
    ret(1,1) = c;
    ret(0,1) = -s;
    ret(1,0) = s;
    return ret;
}

/** \brief Create homogeneous matrix representing a 2D rotation about the coordinate origin.
 
    For use with \ref affineWarpImage(). Angle must be in degrees.
*/
inline
linalg::TemporaryMatrix<double> rotationMatrix2DDegrees(double angle)
{
    return rotationMatrix2DRadians(angle*M_PI/180.0);
}

/** \brief Create homogeneous matrix representing a 2D rotation about the given point.
 
    For use with \ref affineWarpImage(). Angle must be in radians.
*/
inline
linalg::TemporaryMatrix<double> rotationMatrix2DRadians(double angle, TinyVector<double, 2> const & center)
{
    return translationMatrix2D(center) * rotationMatrix2DRadians(angle) * translationMatrix2D(-center);
}

/** \brief Create homogeneous matrix representing a 2D rotation about the given point.
 
    For use with \ref affineWarpImage(). Angle must be in degrees.
*/
inline
linalg::TemporaryMatrix<double> rotationMatrix2DDegrees(double angle, TinyVector<double, 2> const & center)
{
    return rotationMatrix2DRadians(angle*M_PI/180.0, center);
}

/********************************************************/
/*                                                      */
/*                      rotateImage                     */
/*                                                      */
/********************************************************/

/** \brief Rotate image by an arbitrary angle.

    The algorithm performs a rotation about the given center point (the image center by default)
    using the given SplineImageView for interpolation. The destination image must have the same size
    as the source SplineImageView. The rotation is counter-clockwise, and the angle must be given in degrees.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        // rotate about given center point
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor>
        void rotateImage(SplineImageView<ORDER, T> const & src,
                         DestIterator id, DestAccessor dest, 
                         double angleInDegree, TinyVector<double, 2> const & center);
                         
        // rotate about image center
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor>
        void 
        rotateImage(SplineImageView<ORDER, T> const & src,
                    DestIterator id, DestAccessor dest, 
                    double angleInDegree)
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // rotate about given center point
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor>
        void 
        rotateImage(SplineImageView<ORDER, T> const & src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double angleInDegree, TinyVector<double, 2> const & center);

        // rotate about image center
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor>
        void 
        rotateImage(SplineImageView<ORDER, T> const & src,
                    pair<DestImageIterator, DestAccessor> dest, 
                    double angleInDegree);
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> \<vigra/affinegeometry.hxx\><br>
        Namespace: vigra
    
    \code
    
    Image src(width, height);
    vigra::SplineImageView<3, Image::value_type> spline(srcImageRange(src));
    
    Image dest(width, height);
    
    vigra::rotateImage(spline, destImage(dest), 38.5);
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    DestImageIterator dest_upperleft;
    
    double x = ..., y = ...;
    
    if (spline.isInside(x,y))
        dest_accessor.set(spline(x, y), dest_upperleft);

    \endcode
*/
doxygen_overloaded_function(template <...> void rotateImage)

template <int ORDER, class T, 
          class DestIterator, class DestAccessor>
void rotateImage(SplineImageView<ORDER, T> const & src,
                 DestIterator id, DestAccessor dest, 
                 double angleInDegree, TinyVector<double, 2> const & center)
{
    int w = src.width();
    int h = src.height();
    
    double angle = angleInDegree*M_PI/180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);
    
    // avoid round-off errors for simple rotations
    if(closeAtTolerance(std::fmod(angleInDegree, 45.0), 0.0))
    {
        // convert angle into a multiple of pi/4
        int ai = roundi(angleInDegree / 45.0) % 8;
        if(ai < 0)
            ai += 8;
        
        static double sqrt2_2 = 0.5*M_SQRT2;
        static double ss[8] = {0.0, sqrt2_2, 1.0,  sqrt2_2,  0.0, -sqrt2_2, -1.0, -sqrt2_2};
        static double cc[8] = {1.0, sqrt2_2, 0.0, -sqrt2_2, -1.0, -sqrt2_2,  0.0,  sqrt2_2};
        
        s = ss[ai];
        c = cc[ai];
    }
    
    for(int y = 0; y < h; ++y, ++id.y)
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        double sy =  (y - center[1])*c - center[0]*s + center[1];
        double sx = -(y - center[1])*s - center[0]*c + center[0];
        for(int x=0; x < w; ++x, ++rd, sx += c, sy += s)
        {
            if(src.isInside(sx, sy))
                dest.set(src(sx, sy), rd);
        }
    }
}

template <int ORDER, class T, 
          class DestIterator, class DestAccessor>
inline void 
rotateImage(SplineImageView<ORDER, T> const & src,
            pair<DestIterator, DestAccessor> dest, 
            double angleInDegree, TinyVector<double, 2> const & center)
{
    rotateImage(src, dest.first, dest.second, angleInDegree, center);
}

template <int ORDER, class T, 
          class DestIterator, class DestAccessor>
inline void 
rotateImage(SplineImageView<ORDER, T> const & src,
            DestIterator id, DestAccessor dest, 
            double angleInDegree)
{
    TinyVector<double, 2> center((src.width()-1.0) / 2.0, (src.height()-1.0) / 2.0);
    rotateImage(src, id, dest, angleInDegree, center);
}

template <int ORDER, class T, 
          class DestIterator, class DestAccessor>
inline void 
rotateImage(SplineImageView<ORDER, T> const & src,
            pair<DestIterator, DestAccessor> dest, 
            double angleInDegree)
{
    TinyVector<double, 2> center((src.width()-1.0) / 2.0, (src.height()-1.0) / 2.0);
    rotateImage(src, dest.first, dest.second, angleInDegree, center);
}

/********************************************************/
/*                                                      */
/*                  affineWarpImage                     */
/*                                                      */
/********************************************************/

/** \brief Warp an image according to an affine transformation.

    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <int ORDER, class T, 
                class DestIterator, class DestAccessor,
                class C>
        void affineWarpImage(SplineImageView<ORDER, T> const & src,
                            DestIterator dul, DestIterator dlr, DestAccessor dest, 
                            MultiArrayView<2, double, C> const & affineMatrix);
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <int ORDER, class T, 
                class DestIterator, class DestAccessor,
                class C>
        void affineWarpImage(SplineImageView<ORDER, T> const & src,
                            triple<DestIterator, DestIterator, DestAccessor> dest, 
                            MultiArrayView<2, double, C> const & affineMatrix);
    }
    \endcode
    
    The algorithm applies the given \a affineMatrix to the <i>destination coordinates</i> and copies
    the image value from the resulting source coordinates, using the given SplineImageView \a src for interpolation. 
    If the resulting coordinate is outside the source image, nothing will be writen at that destination point.
    
    \code
        for all dest pixels:
            currentSrcCoordinate = affineMatrix * currentDestCoordinate;
            if src.isInside(currentSrcCoordinate):
                dest[currentDestCoordinate] = src[currentSrcCoordinate]; // copy an interpolated value
    \endcode
    
    The matrix represent a 2-dimensional affine transform by means of homogeneous coordinates,
    i.e. it must be a 3x3 matrix whose last row is (0,0,1).
    
    <b> Usage:</b>
    
        <b>\#include</b> \<vigra/affinegeometry.hxx\><br>
        Namespace: vigra
    
    \code
    
    Image src(width, height);
    vigra::SplineImageView<3, Image::value_type> spline(srcImageRange(src));
    
    Image dest1(width, height);
    
    // equivalent (up to round-off errors) with 
    //     rotateImage(spline, destImage(dest1), 45.0);
    TinyVector<double, 2> center((width-1.0)/2.0, (height-1.0)/2.0);
    affineWarpImage(spline, destImageRange(dest1), rotationMatrix2DDegrees(45.0, center));
    
    Image dest2(2*width-1, 2*height-1);
    
    // equivalent (up to round-off errors) with 
    //     resizeImageSplineInterpolation(srcImageRange(img), destImageRange(dest2));
    // note that scaleFactor = 0.5, because we must pass the transformation from destination to source
    affineWarpImage(spline, destImageRange(dest2), scalingMatrix2D(0.5));
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    DestImageIterator dest_upperleft;
    
    double x = ..., y = ...;
    
    if (spline.isInside(x,y))
        dest_accessor.set(spline(x, y), dest_upperleft);

    \endcode
    
    <b>See also:</b> Functions to specify affine transformation: \ref translationMatrix2D(), \ref scalingMatrix2D(), 
                    \ref shearMatrix2D(), \ref rotationMatrix2DRadians(), \ref rotationMatrix2DDegrees()
*/
doxygen_overloaded_function(template <...> void affineWarpImage)

template <int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
void affineWarpImage(SplineImageView<ORDER, T> const & src,
                     DestIterator dul, DestIterator dlr, DestAccessor dest, 
                     MultiArrayView<2, double, C> const & affineMatrix)
{
    vigra_precondition(rowCount(affineMatrix) == 3 && columnCount(affineMatrix) == 3 && 
                       affineMatrix(2,0) == 0.0 && affineMatrix(2,1) == 0.0 && affineMatrix(2,2) == 1.0,
        "affineWarpImage(): matrix doesn't represent an affine transformation with homogeneous 2D coordinates.");
         
    
    double w = dlr.x - dul.x;
    double h = dlr.y - dul.y;
    
    for(double y = 0.0; y < h; ++y, ++dul.y)
    {
        typename DestIterator::row_iterator rd = dul.rowIterator();
        for(double x=0.0; x < w; ++x, ++rd)
        {
            double sx = x*affineMatrix(0,0) + y*affineMatrix(0,1) + affineMatrix(0,2);
            double sy = x*affineMatrix(1,0) + y*affineMatrix(1,1) + affineMatrix(1,2);
            if(src.isInside(sx, sy))
                dest.set(src(sx, sy), rd);
        }
    }
}

template <int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
inline
void affineWarpImage(SplineImageView<ORDER, T> const & src,
                     triple<DestIterator, DestIterator, DestAccessor> dest, 
                     MultiArrayView<2, double, C> const & affineMatrix)
{
    affineWarpImage(src, dest.first, dest.second, dest.third, affineMatrix);
}


//@}

} // namespace vigra


#endif /* VIGRA_AFFINEGEOMETRY_HXX */
