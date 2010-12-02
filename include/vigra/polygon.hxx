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

template<class Point>
struct CCWCompare
{
    Point p0_;
    CCWCompare(const Point &p0)
    : p0_(p0)
    {}

    bool operator()(const Point &a, const Point &b) const
    {
        return (a[1]-p0_[1])*(b[0]-p0_[0]) - (a[0]-p0_[0])*(b[1]-p0_[1]) < 0;
    }
};

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
void convexHull(
    const PointArray1 &points, PointArray2 &convex_hull)
{
    vigra_precondition(points.size() >= 2,
                       "convexHull(): at least two input points are needed.");
    vigra_precondition(points[0].size() == 2,
                       "convexHull(): 2-dimensional points required.");

    typedef typename PointArray1::value_type Point;
    typedef typename Point::value_type Coordinate;

    // find extremal point (min. x, then min. y):
    unsigned int i0 = 0;
    Point p0 = points[0];
    for(unsigned int i = 1; i < points.size(); ++i)
    {
        Coordinate xDiff = points[i][0] - p0[0];
        if(xDiff < 0 || (xDiff == 0 && points[i][1] < p0[1]))
        {
            p0 = points[i];
            i0 = i;
        }
    }

    // sort other points by angle from p0:
    ArrayVector<Point> other(points.begin(), points.begin() + i0);
    other.insert(other.end(), points.begin()+i0+1, points.end());
    
    // the current definition of CCWCompare ensures that points identical to p0
    // end up at the end of the list => those duplicates will be removed during 
    // Graham scan
    std::sort(other.begin(), other.end(), detail::CCWCompare<Point>(p0));
    
    ArrayVector<Point> result(points.size()+1);
    result[0] = p0;
    result[1] = other[0];

    typename ArrayVector<Point>::iterator currentEnd = result.begin() + 1;

    // Graham's scan:
    Point endSegment = *currentEnd - currentEnd[-1];
    Coordinate sa2;
    for(unsigned int i = 1; i < other.size(); ++i)
    {
        if(other[i] == other[i-1] || other[i] == p0) // skip duplicate points
            continue;
        do
        {
            Point diff = other[i] - currentEnd[-1];
            sa2 = diff[0]*endSegment[1] - endSegment[0]*diff[1];
            if(sa2 < 0)
            {
                // point is to the left, add to convex hull:
                *(++currentEnd) = other[i];
                endSegment = other[i] - currentEnd[-1];
            }
            else if(sa2 == 0)
            {
                // points are collinear, keep far one:
                if(diff.squaredMagnitude() > endSegment.squaredMagnitude())
                {
                    *currentEnd = other[i];
                    endSegment = diff;
                }
            }
            else
            {
                // point is to the right, backtracking needed:
                --currentEnd;
                endSegment = *currentEnd - currentEnd[-1];
            }
        }
        while(sa2 > 0);
    }

    // return closed Polygon:
    *(++currentEnd) = p0;
    ++currentEnd;
    std::copy(result.begin(), currentEnd, std::back_inserter(convex_hull));
}

//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
