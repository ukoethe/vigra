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
#include "tinyvector.hxx"
#include "array_vector.hxx"

namespace vigra {

/** \addtogroup MathFunctions
*/
//@{


template<class POINT=TinyVector<double, 2> >
class Polygon 
: protected ArrayVector<POINT>
{
  public:
    typedef ArrayVector<POINT> Base;

    typedef typename POINT                        Point;
    typedef typename Base::value_type             value_type;
    typedef typename Base::reference              reference;
    typedef typename Base::const_reference        const_reference;
    typedef typename Base::pointer                pointer;
    typedef typename Base::const_pointer          const_pointer;
    typedef typename Base::iterator               iterator;
    typedef typename Base::const_iterator         const_iterator;
    typedef typename Base::reverse_iterator       reverse_iterator;
    typedef typename Base::const_reverse_iterator const_reverse_iterator;
    typedef typename Base::size_type              size_type;
    typedef typename Base::difference_type        difference_type;
    typedef typename POINT::value_type            coordinate_type;
    
    using Base::size;
    using Base::empty;
    using Base::front;
    using Base::back;
    using Base::begin;
    using Base::end;
    using Base::rbegin;
    using Base::rend;
    using Base::operator[];

    Polygon(Polygon const & points)
    : Base(points),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    Polygon()
    : length_(0.0),
      lengthValid_(true),
      partialArea_(0.0),
      partialAreaValid_(true)
    {}

    Polygon(size_type n)
    : Base(n),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    template <class InputIterator>
    Polygon(InputIterator b, InputIterator e)
    : Base(b, e),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    void invalidateProperties()
    {
        lengthValid_ = false;
        partialAreaValid_ = false;
    }

    double length() const
    {
        if(!lengthValid_)
        {
            length_ = 0.0;
            for(unsigned int i = 1; i < size(); ++i)
                length_ += ((*this)[i] - (*this)[i-1]).magnitude();
            lengthValid_ = true;
        }
        return length_;
    }

    double partialArea() const
    {
        if(!partialAreaValid_)
        {
            partialArea_ = 0.0;
            for(unsigned int i = 1; i < size(); ++i)
                partialArea_ += ((*this)[i][0]*(*this)[i-1][1] -
                                 (*this)[i][1]*(*this)[i-1][0]);
            partialArea_ /= 2.0;
            partialAreaValid_ = true;
        }
        return partialArea_;
    }

        /// Returns true iff the last and first points are equal.
    bool closed() const
    {
        return size() > 0 && back() == front();
    }

        /** Linearly interpolate at <tt>offset</tt> between knots 
            <tt>index</tt> and <tt>index+1</tt>.
            
            Preconditions: <tt>0 <= index < size()-1</tt> and <tt>0 <= offset <= 1</tt>.
        */
    Point interpolate(unsigned int index, double offset) const
    {
        return (1.0 - offset) * (*this)[index] + offset * (*this)[index+1];
    }

        /**
         * Tests whether the given point lies within this polygon.
         * Requires that this polygon is closed.

         * Points which lie directly on the polylines or coincide with a knot
         * are considered inside (this behavior is consistent with fillPolygon()).
         * Parameter \a tolerance (interpreted as an absolute error bound)
         * controls the numerical accuracy of this test.
         */
    bool contains(const_reference point, 
                  coordinate_type tolerance=2.0*NumericTraits<coordinate_type>::epsilon()) const;

    void push_back(const_reference v)
    {
        if(size())
        {
            if(lengthValid_)
                length_ += (v - back()).magnitude();
            if(partialAreaValid_)
                partialArea_ += (v[0]*back()[1] - v[1]*back()[0]);
        }
        Base::push_back(v);
    }

    void extend(const Polygon &other)
    {
        if(!other.size())
            return;

        const_iterator otherBegin(other.begin());
        if(size())
        {
            if(*otherBegin == back())
            {
                // don't copy first pixel
                ++otherBegin;
            }
            else
            {
                if(lengthValid_)
                    length_ += (other.front() - back()).magnitude();
                if(partialAreaValid_)
                    partialArea_ += (other.front()[0]*back()[1] -
                                     other.front()[1]*back()[0]);
            }
        }
        if(lengthValid_)
            length_ += other.length();
        if(partialAreaValid_)
            partialArea_ += other.partialArea();
        Base::insert(end(), otherBegin, other.end());
    }

    void setPoint(unsigned int pos, const_reference x)
    {
        if(lengthValid_)
        {
            if(pos > 0)
            {
                length_ += (x - (*this)[pos-1]).magnitude() -
                           ((*this)[pos] - (*this)[pos-1]).magnitude();
            }
            if(pos < size() - 1)
            {
                length_ += (x - (*this)[pos+1]).magnitude() -
                           ((*this)[pos] - (*this)[pos+1]).magnitude();
            }
        }
        partialAreaValid_ = false;
        (*this)[pos] = x;
    }

    void erase(iterator pos)
    {
        invalidateProperties();
        Base::erase(pos);
    }

    iterator insert(iterator pos, const_reference x)
    {
        if(lengthValid_)
        {
            if(pos > begin())
                length_ += (x - pos[-1]).magnitude();
            if(end() - pos >= 1)
            {
                length_ += (x - *pos).magnitude();
                if(pos > begin())
                    length_ -= (*pos - pos[-1]).magnitude();
            }
        }
        partialAreaValid_ = false;
        return Base::insert(pos, x);
    }

    template <class InputIterator>
    iterator insert(iterator pos, InputIterator i, InputIterator end)
    {
        partialAreaValid_ = false;
        lengthValid_ = false;
        return Base::insert(pos, i, end);
    }

    Polygon split(unsigned int pos)
    {
        if(pos == 0)
        {
            Polygon result(1);
            result[0] = (*this)[0];
            swap(result);
            return result;
        }

        Polygon result(begin() + pos, end());
        Base::erase(begin() + pos + 1, end());

#if 0 // FIXME: somehow this does not work!
        if(pos > size() * 2 / 3)
        {
             // heuristic: when splitting off only a "small part",
             // re-use existing information
            if(lengthValid_)
                length_ -= result.length();
            if(partialAreaValid_)
                partialArea_ -= result.partialArea();
        }
        else
#endif
            invalidateProperties();

        return result;
    }

    template <class Sequence>
    void arcLengthList(Sequence & arcLengths) const
    {
        double length = 0.0;
        arcLengths.push_back(0.0);
        for(unsigned int i = 1; i < size(); ++i)
        {
            length += ((*this)[i] - (*this)[i-1]).magnitude();
            arcLengths.push_back(length);
        }
    }

    void swap(Polygon &rhs)
    {
        Base::swap(rhs);
        std::swap(length_, rhs.length_);
        std::swap(lengthValid_, rhs.lengthValid_);
        std::swap(partialArea_, rhs.partialArea_);
        std::swap(partialAreaValid_, rhs.partialAreaValid_);
    }

    void reverse()
    {
        Base::reverse();
        if(partialAreaValid_)
            partialArea_ = -partialArea_;
    }

    POINT nearestPoint(const_reference p) const;

    Polygon operator+(const Point &offset) const
    {
        Polygon result(size());
        for(unsigned int i = 0; i < size(); ++i)
            result[i] = (*this)[i] + offset;
        return result;
    }

    Polygon operator-(const Point &offset) const
    {
        return operator+(-offset);
    }

    Polygon operator*(double scale) const
    {
        Polygon result(size());
        for(unsigned int i = 0; i < size(); ++i)
            result[i] = (*this)[i] * scale;
        return result;
    }
    
  protected:
  
    Base & points()
    {
        return (Base &)*this;
    }
  
    Base const & points() const
    {
        return (Base const &)*this;
    }
    
    mutable double length_;
    mutable bool lengthValid_;
    mutable double partialArea_;
    mutable bool partialAreaValid_;
};

template <class POINT>
POINT Polygon<POINT>::nearestPoint(const_reference p) const
{
    double dist = NumericTraits<double>::max();
    POINT r;
    for(unsigned int k=1; k<size(); ++k)
    {
        POINT dp = (*this)[k] - (*this)[k-1];
        POINT dc = p - (*this)[k-1];
        double t = dot(dp, dc);
        if(t != 0.0)
            t /= squaredNorm(dp);
        if(t > 1.0)
        {
            double d = norm((*this)[k]-p);
            if (d < dist)
            {
                dist = d;
                r = (*this)[k];
            }
        }
        else if(t < 0.0)
        {
            double d = norm((*this)[k-1]-p);
            if (d < dist)
            {
                dist = d;
                r = (*this)[k-1];
            }
        }
        else
        {
            POINT pp = (*this)[k-1] + t*dp;
            double d = norm(pp-p);
            if (d < dist)
            {
                dist = d;
                r = pp;
            }
        }
    }
    return r;
}

template <class POINT>
bool 
Polygon<POINT>::contains(const_reference point, 
                         coordinate_type tolerance) const
{
    vigra_precondition(closed(),
                       "Polygon::contains() requires polygon to be closed!");
    
    // NOTE: the following code is very similar to detail::createScanIntervals()
    //       to ensure consistent results.
    
    Polygon p = (*this) - point; // shift the polygon so that we only need to test
                                 // for intersections with scanline 0
    int n = p.size();
    for(int k=0; k<n; ++k)
        if(closeAtTolerance(p[k][1], 0.0, tolerance))
            p[k][1] = 0.0;
        
    int result = 0;
    bool drop_next_start_point = false;
    int first_point_maybe_dropped = -1;
    for(int k=0; k<n-1; ++k)
    {
        Point const & p1 = p[k];
        Point const & p2 = p[k+1];
        
        if(p1[1] == p2[1]) // ignore horizontal lines
            continue;

        double t = (p2[0] - p1[0]) / (p2[1] - p1[1]);
        double y, yend, dy;
        if(p1[1] < p2[1])
        {
            y = ceil(p1[1]);
            yend = floor(p2[1]);
            dy = 1.0;
        }
        else
        {
            y = floor(p1[1]);
            yend = ceil(p2[1]);
            dy = -1.0;
        }
        if(yend != p2[1])
            yend += dy;
        if(drop_next_start_point)
        {
            y += dy;
            drop_next_start_point = false;
        }
        if(first_point_maybe_dropped == -1)
        {
            if(y == 0.0 && p1[0] - p1[1]*t < 0.0)
                first_point_maybe_dropped = 1;
            else
                first_point_maybe_dropped = 0;
        }
        if(y*dy <= 0.0 && yend*dy > 0.0)  // intersects scanline 0
        {
            double x = p1[0] - p1[1]*t;
            if(closeAtTolerance(x, 0.0, tolerance))
                return true;
            if(x < 0.0)
                ++result;
        }
        else if(p2[1] == 0.0) // degenerate case
        {
            int j = (k+2)%n;
            bool convex = detail::orderedClockwise(p1, p2, p[j]);
            if(convex) 
            {
                double x = p2[0] - p2[1]*t;
                if(closeAtTolerance(x, 0.0, tolerance))
                    return true;
                if(x < 0.0)
                    ++result;
            }
            for(; j != k+1; j = (j+1)%n)
            {
                double bend = dy*(p[j][1] - yend);
                if(bend == 0.0)
                    continue;
                // Drop startpoint of next segment when the polygon after a convex 
                // degenerate knot eventually crosses the scanline, or when it 
                // returns to the original side of the scanline after a concave 
                // degenerate knot.
                if((convex && bend > 0.0) || (!convex && bend < 0.0))
                    drop_next_start_point = true;
                break;
            }
        }
    }
    
    if(drop_next_start_point && first_point_maybe_dropped == 1)
        --result;

    return (result % 2) != 0;
}

template <class POINT>
Polygon<Shape2> roundi(Polygon<POINT> const & p) 
{
    Polygon<Shape2> result(p.size());
    for(unsigned int i = 0; i < p.size(); ++i)
    {
        result[i] = round(p[i]);
    }
    return result;
}

} // namespace vigra

/********************************************************************/

namespace std {

template<class T>
void swap(vigra::Polygon<T> &a, vigra::Polygon<T> &b)
{
    a.swap(b);
}

} // namespace std

/********************************************************************/

namespace vigra {

/** \brief Create a polygon from the interpixel contour of a labeled region.

    The point \a anchor_point must be in the region whose contour we want to extract,
    and must be adjacent to the contour. The algorithm uses the 'left hand on the wall'
    algorithm to trace the connected component whose label equals the label of the
    \a anchor_point. The contour is returned in \a countour_points as a closed polygon 
    that circles the region counter-clockwise in the image coordinate system (i.e. the
    coordinate system where x points to the right and y points downwards). Since the
    resulting polygon represents the interpixel contour, all points will have one integer 
    and one half-integer coordinate.
*/
template<class T, class S, class PointArray>
void 
extractContour(MultiArrayView<2, T, S> const &label_image,
               Shape2 const & anchor_point,
               PointArray & contour_points) 
{
    typedef typename PointArray::value_type Point;
    
    Shape2 step[4] = { Shape2(0, -1), Shape2(1, 0), Shape2(0, 1), Shape2(-1, 0) };
    Point contour_offsets[4] = { Point(-0.5, 0), Point(0, -0.5), Point(0.5, 0), Point(0, 0.5) };
    
    T foreground = label_image[anchor_point];

    int direction;
    Shape2 position;
    // find a position outside the object from which we place the hand
    for(direction = 3; direction >= 0; --direction)
    {
        position = anchor_point + step[(direction + 1) % 4];
        if(!label_image.isInside(position) || label_image[position] != foreground)
            break;
    }
    
    vigra_precondition(direction >= 0,
        "extractContour(): the anchor point must be at the region border.");
        
    int initial_direction = direction;
    Shape2 initial_position = position;

    // go around the object
    do 
    {
        contour_points.push_back(position + contour_offsets[direction]);

        Shape2 next_position = position + step[direction];

        if(label_image.isInside(next_position) && 
           label_image[next_position] == foreground)
        {
            // we have bumped into a wall => turn right to touch the wall again
            direction = (direction + 1) % 4;
        }
        else
        {
            position = next_position;
            int next_direction = (direction + 3) % 4;
            next_position += step[next_direction];
            if(!label_image.isInside(next_position) || 
               label_image[next_position] != foreground)
            {
                // we have lost the wall => turn left and move forward to touch the wall again
                direction = next_direction;
                position = next_position;
            }
        }
    } 
    while (position != initial_position || direction != initial_direction);
    
    contour_points.push_back(contour_points.front()); // make it a closed polygon
}

namespace detail {

template < class Point >    
bool pointYXOrdering(Point const & p1, Point const & p2) 
{
    return (p1[1]<p2[1]) || (p1[1] == p2[1] && p1[0] < p2[0]);
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
    of the input. The function implements Andrew's Monotone Chain algorithm.
*/
template<class PointArray1, class PointArray2>
void convexHull(const PointArray1 &points, PointArray2 & convex_hull)
{
    vigra_precondition(points.size() >= 2,
                       "convexHull(): at least two input points are needed.");
    vigra_precondition(points[0].size() == 2,
                       "convexHull(): 2-dimensional points required.");
    
    typedef typename PointArray1::value_type Point;
    
    auto begin = points.begin();
    if(points.front() == points.back()) // closed polygon
        ++begin;                        // => remove redundant start point
    ArrayVector<Point> ordered(begin, points.end());
    std::sort(ordered.begin(), ordered.end(), detail::pointYXOrdering<Point>);
    
    ArrayVector<Point> H;
    
    int n = ordered.size(), k=0;
    
    // Build lower hull
    for (int i = 0; i < n; i++) 
    {
        while (k >= 2 && detail::orderedClockwise(H[k-2], H[k-1], ordered[i])) 
        {
            H.pop_back();
            --k;
        }
        H.push_back(ordered[i]);
        ++k;
    }
    
    // Build upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) 
    {
        while (k >= t && detail::orderedClockwise(H[k-2], H[k-1], ordered[i])) 
        {
            H.pop_back();
            --k;
        }
        H.push_back(ordered[i]);
        ++k;
    }
    
    for(int i=k-1; i>=0; --i)
        convex_hull.push_back(H[i]);
}

/********************************************************************/
/*                                                                  */
/*                         polygon drawing                          */
/*                                                                  */
/********************************************************************/

namespace detail {

/*
 * Find and sort all intersection points of the polygon with scanlines.
 * Polygons are considered as closed set, i.e. pixels on the polygon 
 * contour are included. The function handles degenerate cases (i.e.
 * knots on scanlines) correctly.
 */
template<class Point, class Array>
void createScanIntervals(Polygon<Point> const &p, Array & result) 
{
    bool drop_next_start_point = false;
    int n = p.size();
    for(int k=0; k<n-1; ++k)
    {
        Point const & p1 = p[k];
        Point const & p2 = p[k+1];
        
        if(p1[1] == p2[1]) // ignore horizontal lines
            continue;

        double t = (p2[0] - p1[0]) / (p2[1] - p1[1]);
        double y, yend, dy;
        if(p1[1] < p2[1])
        {
            y = ceil(p1[1]);
            yend = floor(p2[1]);
            dy = 1.0;
        }
        else
        {
            y = floor(p1[1]);
            yend = ceil(p2[1]);
            dy = -1.0;
        }
        if(yend != p2[1]) // in general don't include the segment's endpoint
            yend += dy;   // (since it is also the start point of the next segment)
        if(drop_next_start_point) // handle degeneracy from previous iteration
        {
            y += dy;
            drop_next_start_point = false;
        }
        for(; (y-yend)*dy < 0.0; y += dy)  // compute scanline intersections
        {
            double x = p1[0] + (y - p1[1])*t;
            result.push_back(Point(x,y));
        }
        if(yend == p2[1]) // degenerate case: p2 is exactly on a scanline (yend is integer)
        {
            int j = (k+2)%n;
            bool convex = detail::orderedClockwise(p1, p2, p[j]);
            if(convex) // include the segment's endpoint p2 when it is a convex knot
            {
                result.push_back(p2);
            }
            for(; j != k+1; j = (j+1)%n)
            {
                double bend = dy*(p[j][1] - yend);
                if(bend == 0.0)
                    continue;
                // Drop startpoint of next segment when the polygon after a convex 
                // degenerate knot eventually crosses the scanline, or when it 
                // returns to the original side of the scanline after a concave 
                // degenerate knot.
                if((convex && bend > 0.0) || (!convex && bend < 0.0))
                    drop_next_start_point = true;
                break;
            }
        }
    }
    
    if(drop_next_start_point)
        result.erase(result.begin());
    
    vigra_invariant((result.size() & 1) == 0,
        "createScanIntervals(): internal error - should return an even number of points.");
    sort(result.begin(), result.end(), pointYXOrdering<Point>);
}


} // namespace detail

/** \brief Render closed polygon \a p into the image \a output_image. 

    All pixels on the polygon's contour and in its interior are 
    set to the given \a value. Parts of the polygon outside the image
    region are clipped. The function uses a robust X-intersection array 
    algorithm that is able to handle all corner cases (concave and 
    self-intersecting polygons, knots on integer coordinates).
 */
template<class Point, class T, class S>
void fillPolygon(Polygon<Point> const &p,
                 MultiArrayView<2, T, S> &output_image, 
                 T value) 
{
    vigra_precondition(p.closed(),
        "fillPolygon(): polygon must be closed (i.e. first point == last point).");
        
    std::vector<Point> scan_intervals;
    detail::createScanIntervals(p, scan_intervals);

    for(unsigned int k=0; k < scan_intervals.size(); k+=2)
    {
        MultiArrayIndex x    = (MultiArrayIndex)ceil(scan_intervals[k][0]),
                        y    = (MultiArrayIndex)scan_intervals[k][1],
                        xend = (MultiArrayIndex)floor(scan_intervals[k+1][0]) + 1;
        vigra_invariant(y == scan_intervals[k+1][1],
            "fillPolygon(): internal error - scan interval should have same y value.");
        // clipping
        if(y < 0)
            continue;
        if(y >= output_image.shape(1))
            break;
        if(x < 0)
            x = 0;
        if(xend > output_image.shape(0))
            xend = output_image.shape(0);
        // drawing
        for(; x < xend; ++x)
            output_image(x,y) = value;
    }
}

//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
