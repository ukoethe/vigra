/************************************************************************/
/*                                                                      */
/*               Copyright 2005-2006 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_AFFINEGEOMETRY_HXX
#define VIGRA_AFFINEGEOMETRY_HXX

#include "vigra/mathutil.hxx"
#include "vigra/splineimageview.hxx"
#include <cmath>

namespace vigra {

/** \addtogroup GeometricTransformations Geometric Transformations
*/
//@{

/********************************************************/
/*                                                      */
/*                      rotateImage                     */
/*                                                      */
/********************************************************/

/** \brief Rotate image by an arbitrary angle.

    The algorithm performs a rotation about the given center point (the image center by default)
    using the given SplineImageView for interpolation. The destination imae must have the same size
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
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
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
    
        <b>\#include</b> "<a href="affinegeometry_8hxx-source.html">vigra/affinegeometry.hxx</a>"<br>
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

//@}

} // namespace vigra


#endif /* VIGRA_AFFINEGEOMETRY_HXX */
