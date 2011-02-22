/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2010 by Ullrich Koethe                  */
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

#ifndef VIGRA_POLYGON_HXX
#define VIGRA_POLYGON_HXX

#include <cmath>
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"

namespace vigra {

/** \addtogroup MathFunctions
*/
//@{

namespace detail {

template < class Point >    
bool sortPoints(Point const & p1, Point const & p2) 
{
    return (p1[0]<p2[0]) || (p1[0] == p2[0] && p1[1] < p2[1]);
}

template < class Point >    
bool orderedClockwise(const Point &O, const Point &A, const Point &B)
{
    return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]) <= 0;
}

} // namespace detail


/** \brief Compute convex hull of a 2D polygon.

    The input array \a points contains a (not necessarily ordered) set of 2D points
    whose convex hull is to be computed. The array's <tt>value_type</tt> (i.e. the point type)
    must be compatible with std::vector (in particular, it must support indexing, 
    copying, and have <tt>size() == 2</tt>). The points of the convex hull will be appended
    to the output array \a convex_hull (which must support <tt>std::back_inserter(convex_hull)</tt>). 
    Since the convex hull is a closed polygon, the first and last point of the output will 
    be the same (i.e. the first point will simply be inserted at the end again). The points
    of the convex hull will be ordered counter-clockwise, starting with the leftmost point
    of the imput.
*/
template<class PointArray1, class PointArray2>
void convexHull(const PointArray1 &points, PointArray2 & convex_hull)
{
    vigra_precondition(points.size() >= 2,
                       "convexHull(): at least two input points are needed.");
    vigra_precondition(points[0].size() == 2,
                       "convexHull(): 2-dimensional points required.");
    
    typedef typename PointArray1::value_type Point;
    
    ArrayVector<Point> ordered(points.begin(), points.end());
    std::sort(ordered.begin(), ordered.end(), detail::sortPoints<Point>);
    
    ArrayVector<Point> H(points.size()*2);
    
    int n = points.size(), k=0;
    
    // Build lower hull
    for (int i = 0; i < n; i++) 
    {
        while (k >= 2 && detail::orderedClockwise(H[k-2], H[k-1], ordered[i])) 
        {
            k--;
        }
        H[k++] = ordered[i];
    }
    
    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) 
    {
        while (k >= t && detail::orderedClockwise(H[k-2], H[k-1], ordered[i])) 
        {
            k--;
        }
        H[k++] = ordered[i];
    }
    
    convex_hull.resize(k);
    for (int i = 0; i < k; ++i) 
    {
        convex_hull[i] = H[i];
    }
}

//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
