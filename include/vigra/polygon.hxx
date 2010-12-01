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
#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"

namespace vigra {

/** \addtogroup MathFunctionsons and functors.
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

template<class PointArray>
void convexHull(
    const PointArray &poly, PointArray &chull)
{
    vigra_precondition(poly.size() >= 2,
                       "convexHull(): at least two input points are needed.");
    vigra_precondition(poly[0].size() == 2,
                       "convexHull(): 2-dimensional points required.");

    typedef typename PointArray::value_type Point;
    typedef typename Point::value_type Coordinate;

    // find extremal point (min. x, then min. y):
    unsigned int i0 = 0;
    Point p0 = poly[0];
    for(unsigned int i = 1; i < poly.size(); ++i)
    {
        Coordinate xDiff = poly[i][0] - p0[0];
        if(xDiff < 0 || (xDiff == 0 && poly[i][1] < p0[1]))
        {
            p0 = poly[i];
            i0 = i;
        }
    }

    // sort other points by angle from p0:
    ArrayVector<Point> other(poly.begin(), poly.begin() + i0);
    other.insert(other.end(), poly.begin()+i0+1, poly.end());
    std::sort(other.begin(), other.end(), detail::CCWCompare<Point>(p0));
    
    ArrayVector<Point> result(poly.size()+1);
    result[0] = p0;
    result[1] = other[0];
    typename ArrayVector<Point>::iterator currentEnd = result.begin() + 1;

    // Graham's scan:
    Point endSegment = *currentEnd - currentEnd[-1];
    Coordinate sa2;
    for(unsigned int i = 1; i < other.size(); ++i)
    {
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
    chull.insert(chull.end(), result.begin(), currentEnd);
}

//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
