/************************************************************************/
/*                                                                      */
/*               Copyright 2013-2017 by Benjamin Seppke                 */
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

#ifndef VIGRA_PIECEWISEAFFINE_REGISTRATION_HXX
#define VIGRA_PIECEWISEAFFINE_REGISTRATION_HXX

#ifndef WITH_LEMON
    #error "Should only be included with flag \"WITH_LEMON\""
#endif

#include <lemon/list_graph.h>
#include <lemon/dim2.h>

#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"
#include "affine_registration.hxx"
#include "delaunay.hxx"

namespace vigra
{

/** \addtogroup Registration Image Registration
 */
//@{

//used vigra types for points and triangles
typedef vigra::TinyVector<double,2> PointType;
typedef vigra::triple<PointType, PointType, PointType> TriangleType;
typedef std::pair<TriangleType, vigra::Matrix<double> > TriangleTransformationType;
    
/**
 * Tiny helper to check whether a point is inside a triangl.
 */
inline
bool isInside(const PointType & point, const TriangleType & triangle, bool allowBorder=true)
{
    // Compute direction vectors
    PointType v0 = triangle.third  - triangle.first;
    PointType v1 = triangle.second - triangle.first;
    PointType v2 = point -           triangle.first;
    
    // Compute dot products
    double    dot00 = dot(v0, v0),
            dot01 = dot(v0, v1),
            dot02 = dot(v0, v2),
            dot11 = dot(v1, v1),
            dot12 = dot(v1, v2);
    
    // Compute barycentric coordinates
    double    invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01),
                   u = (dot11 * dot02 - dot01 * dot12) * invDenom,
                   v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    // Check if point is in triangle
    if (allowBorder)
    {
        return (u >= 0) && (v >= 0) && (u + v <= 1);
    }
    else
    {
        return (u > 0) && (v > 0) && (u + v < 1);
    }
}

/********************************************************/
/*                                                      */
/*  triangleTransformations2DFromCorrespondingPoints    */
/*                                                      */
/********************************************************/

/** \brief Create a list of triangles of the dest image together with an affine transform for each.

    For use with \ref piecewiseAffineWarpImage().

    This function performs a delaunay triangulation on the source image an computes the
    affine transform from the dest image triangles to the first one.
 
    For the delaunay triangulation, LEMON is neccessary.
*/
template <class SrcPointIterator, class DestPointIterator>
std::vector<TriangleTransformationType>
triangleTransformations2DFromCorrespondingPoints(SrcPointIterator s, SrcPointIterator s_end, DestPointIterator d)
{
    using namespace lemon;
    typedef dim2::Point<double> Point;
    
    int point_count = s_end - s;
    
    //Create a new graph with all the points as nodes
    ListGraph g;
    ListGraph::NodeMap<Point> coords(g);
    ListGraph::NodeMap<int> idxs(g);
    
    //Add nodes and their geometric embedding
    for (int i=0; i<point_count;++i)
    {
        ListGraph::Node n = g.addNode();
        coords[n] = Point(s[i][0], s[i][1]);
        idxs[n] = i;
    }
    
    //Perform a delaunay triangulation on the graph
    delaunay(g, coords);
    
    //Get the Triangles by means of triples of node indices
    std::vector<vigra::triple<int,int,int> > triangles = trianglesFromDelaunay(g, idxs);
    
    //Fill the transformation matrixes
    std::vector<PointType> s_points(3), d_points(3);
    std::vector<TriangleTransformationType> result;
    for (auto iter = triangles.begin(); iter!=triangles.end(); iter++)
    {
        int idx_a = iter->first,
            idx_b = iter->second,
            idx_c = iter->third;
        
        s_points[0] = s[idx_a]; s_points[1] = s[idx_b]; s_points[2] = s[idx_c];
        d_points[0] = d[idx_a]; d_points[1] = d[idx_b]; d_points[2] = d[idx_c];
        
        result.push_back(TriangleTransformationType(TriangleType(d_points[0],d_points[1],d_points[2]),
                                                    affineMatrix2DFromCorrespondingPoints(d_points.begin(), d_points.end(), s_points.begin())));
    }
    
    return result;
}

/********************************************************/
/*                                                      */
/*                piecewiseAffineWarpImage              */
/*                                                      */
/********************************************************/

/** \brief Warp an image according to a piecewise affine transformation.

    To get more information about the structure of the triangle+matrix
    vector, see \ref polynomialMatrix2DFromCorrespondingPoints().

    <b>\#include</b> \<vigra/piecewiseaffine_registration.hxx\><br>
    Namespace: vigra

    pass 2D array views:
    \code
    namespace vigra {
        template <int ORDER, class T,
                  class T2, class S2,
                  class C>
        void
        piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                                 MultiArrayView<2, T2, S2> dest,
                                 const std::vector<TriangleTransformationType> & tri_trans);
    }
    \endcode

    \deprecatedAPI{piecewiseAffineWarpImage}

    pass arguments explicitly:
    \code
    namespace vigra {
        template <int ORDER, class T,
                  class DestIterator, class DestAccessor,
                  class C>
        void piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                                 DestIterator dul, DestIterator dlr, DestAccessor dest,
                                 const std::vector<TriangleTransformationType> & tri_trans);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <int ORDER, class T,
                  class DestIterator, class DestAccessor,
                  class C>
        void piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                                 triple<DestIterator, DestIterator, DestAccessor> dest,
                                 const std::vector<TriangleTransformationType> & tri_trans);
    }
    \endcode
    \deprecatedEnd
 */
doxygen_overloaded_function(template <...> void piecewiseAffineWarpImage)

template <int ORDER, class T,
          class DestIterator, class DestAccessor>
void piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                              DestIterator dul, DestIterator dlr, DestAccessor dest,
                              const std::vector<TriangleTransformationType> & tri_trans)
{
    double w = dlr.x - dul.x;
    double h = dlr.y - dul.y;
    
    for(double y = 0.0; y < h; ++y, ++dul.y)
    {
        typename DestIterator::row_iterator rd = dul.rowIterator();
        for(double x=0.0; x < w; ++x, ++rd)
        {
            for(std::vector<TriangleTransformationType>::const_iterator iter=tri_trans.begin(); iter!=tri_trans.end(); ++iter)
            {
                if(isInside(PointType(x,y), iter->first))
                {
                    vigra::Matrix<double> transformation = iter->second;
                    double sx = transformation(0,0)*x + transformation(0,1)*y + transformation(0,2);
                    double sy = transformation(1,0)*x + transformation(1,1)*y + transformation(1,2);
                    
                    if(src.isInside(sx, sy))
                        dest.set(src(sx, sy), rd);
                }
            }
        }
    }
}

template <int ORDER, class T,
          class DestIterator, class DestAccessor>
inline
void piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                              triple<DestIterator, DestIterator, DestAccessor> dest,
                              const std::vector<TriangleTransformationType> & tri_trans)
{
    piecewiseAffineWarpImage(src, dest.first, dest.second, dest.third, tri_trans);
}


template <int ORDER, class T,
          class T2, class S2>
inline
void piecewiseAffineWarpImage(SplineImageView<ORDER, T> const & src,
                              MultiArrayView<2, T2, S2> dest,
                              const std::vector<TriangleTransformationType> & tri_trans)
{
    piecewiseAffineWarpImage(src, destImageRange(dest), tri_trans);
}


//@}

} // namespace vigra


#endif /* VIGRA_PIECEWISEAFFINE_REGISTRATION_HXX */
