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
#include "gaussians.hxx"
#include "splines.hxx"
#include "linear_solve.hxx"

namespace vigra {

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


/** \addtogroup Geometry
*/
//@{

    /** Polygons in two and higher dimenions.
    
    */
template<class POINT=TinyVector<double, 2> >
class Polygon 
: protected ArrayVector<POINT>
{
  public:
    typedef ArrayVector<POINT> Base;

    typedef POINT                                 Point;
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
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
    using Base::rbegin;
    using Base::rend;
    using Base::crbegin;
    using Base::crend;

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
    : Base(n, Point()),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    template <class InputIterator>
    Polygon(InputIterator b, InputIterator e)
    : Base(b, e),
      lengthValid_(false),
      partialAreaValid_(false)
    {}
    
    void clear()
    {
        invalidateProperties();
        Base::clear();
    }

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
            partialArea_ *= 0.5;
            partialAreaValid_ = true;
        }
        return partialArea_;
    }

    double area() const
    {
        vigra_precondition(closed(),
                           "Polygon::area() requires polygon to be closed!");
        return abs(partialArea());
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
    bool contains(const_reference point, coordinate_type tolerance) const;

    bool contains(const_reference point) const
    {
        return contains(point, 2.0*NumericTraits<coordinate_type>::epsilon());
    }
    
    void push_back(const_reference v)
    {
        if(size())
        {
            if(lengthValid_)
                length_ += (v - back()).magnitude();
            if(partialAreaValid_)
                partialArea_ += 0.5*(v[0]*back()[1] - v[1]*back()[0]);
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
            if(*otherBegin == Base::back())
            {
                // don't copy first pixel
                ++otherBegin;
            }
            else
            {
                if(lengthValid_)
                    length_ += (other.front() - Base::back()).magnitude();
                if(partialAreaValid_)
                    partialArea_ += (other.front()[0]*Base::back()[1] -
                                     other.front()[1]*Base::back()[0]);
            }
        }
        if(lengthValid_)
            length_ += other.length();
        if(partialAreaValid_)
            partialArea_ += other.partialArea();
        Base::insert(Base::end(), otherBegin, other.end());
    }

    void setPoint(unsigned int pos, const_reference x)
    {
        invalidateProperties();
        Base::operator[](pos) = x;
    }

        // doesn't call invalidateProperties()
    void setPointUnsafe(unsigned int pos, const_reference x)
    {
        Base::operator[](pos) = x;
    }

        // alternative, but it will also invalidate if the caller only reads 
        // reads the return value.
    // reference operator[](unsigned int pos)
    // {
        // invalidateProperties();
        // return Base::operator[](pos);
    // }

    const_reference operator[](unsigned int pos) const
    {
        return Base::operator[](pos);
    }

    const_reference front() const
    {
        return Base::front();
    }
    
    const_reference back() const
    {
        return Base::back();
    }
    
    iterator begin()
    {
        invalidateProperties();
        return Base::begin();
    }
    
    iterator end()
    {
        invalidateProperties();
        return Base::end();
    }
    
    reverse_iterator rbegin()
    {
        invalidateProperties();
        return Base::rbegin();
    }
    
    reverse_iterator rend()
    {
        invalidateProperties();
        return Base::rend();
    }

    void erase(iterator pos)
    {
        invalidateProperties();
        Base::erase(pos);
    }

    void erase(iterator pos, iterator end)
    {
        invalidateProperties();
        Base::erase(pos, end);
    }

    iterator insert(iterator pos, const_reference x)
    {
        invalidateProperties();
        return Base::insert(pos, x);
    }

    template <class InputIterator>
    iterator insert(iterator pos, InputIterator i, InputIterator end)
    {
        invalidateProperties();
        return Base::insert(pos, i, end);
    }

    Polygon split(unsigned int pos)
    {
        Polygon result;
        if(pos == 0)
        {
            swap(result);
        }
        else if(pos < size())
        {
            result.insert(result.begin(), begin() + pos, end());
            erase(begin() + pos, end());
        }
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
        std::reverse(Base::begin(), Base::end());
        if(partialAreaValid_)
            partialArea_ = -partialArea_;
    }

    POINT nearestPoint(const_reference p) const;
    
    Polygon & operator+=(POINT const & offset)
    {
        if(!closed())
            partialAreaValid_ = false;
        for(unsigned int i = 0; i < size(); ++i)
            Base::operator[](i) += offset;
        return *this;
    }

    Polygon & operator-=(POINT const & offset)
    {
        if(!closed())
            partialAreaValid_ = false;
        for(unsigned int i = 0; i < size(); ++i)
            Base::operator[](i) -= offset;
        return *this;
    }

    Polygon & operator*=(double scale)
    {
        partialArea_ *= sq(scale);
        length_ *= scale;
        for(unsigned int i = 0; i < size(); ++i)
            Base::operator[](i) *= scale;
        return *this;
    }

    Polygon & operator/=(double scale)
    {
        partialArea_ /= sq(scale);
        length_ /= scale;
        for(unsigned int i = 0; i < size(); ++i)
            Base::operator[](i) /= scale;
        return *this;
    }

    bool operator==(Polygon const & rhs) const
    {
        if(size() != rhs.size())
            return false;
        for(size_type k=0; k<size(); ++k)
            if((*this)[k] != rhs[k])
                return false;
        return true;
    }

    bool operator!=(Polygon const & rhs) const
    {
        return !((*this) == rhs);
    }
    
  protected:
    
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
    typedef typename Polygon<POINT>::Base Base;
    vigra_precondition(closed(),
                       "Polygon::contains() requires polygon to be closed!");
    
    // NOTE: the following code is very similar to detail::createScanIntervals()
    //       to ensure consistent results.
    
    Polygon p = (*this) - point; // shift the polygon so that we only need to test
                                 // for intersections with scanline 0
    int n = p.size();
    for(int k=0; k<n; ++k)
        if(closeAtTolerance(p[k][1], 0.0, tolerance))
            ((Base&)p)[k][1] = 0.0;
        
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
inline Polygon<POINT> round(Polygon<POINT> const & p) 
{
    Polygon<POINT> result(p.size());
    for(unsigned int i = 0; i < p.size(); ++i)
    {
        result.setPointUnsafe(i, round(p[i]));
    }
    return result;
}

template <class POINT>
inline Polygon<TinyVector<std::ptrdiff_t, 2> > roundi(Polygon<POINT> const & p) 
{
    Polygon<TinyVector<std::ptrdiff_t, 2> > result(p.size());
    for(unsigned int i = 0; i < p.size(); ++i)
    {
        result.setPointUnsafe(i,roundi(p[i]));
    }
    return result;
}

template <class POINT>
inline Polygon<POINT> 
operator+(Polygon<POINT> const & p, POINT const & offset)
{
    return Polygon<POINT>(p) += offset;
}

template <class POINT>
inline Polygon<POINT> 
operator+(POINT const & offset, Polygon<POINT> const & p)
{
    return Polygon<POINT>(p) += offset;
}

template <class POINT>
inline Polygon<POINT> 
operator-(Polygon<POINT> const & p)
{
    Polygon<POINT> result(p.size());
    for(unsigned int i = 0; i < p.size(); ++i)
        result.setPointUnsafe(i, -p[i]);
    return result;
}

template <class POINT>
inline Polygon<POINT> 
operator-(Polygon<POINT> const & p, POINT const & offset)
{
    return Polygon<POINT>(p) -= offset;
}

template <class POINT>
inline Polygon<POINT> 
operator*(Polygon<POINT> const & p, double scale)
{
    return Polygon<POINT>(p) *= scale;
}

template <class POINT>
inline Polygon<POINT> 
operator*(double scale, Polygon<POINT> const & p)
{
    return Polygon<POINT>(p) *= scale;
}


template <class POINT>
inline Polygon<POINT> 
operator/(Polygon<POINT> const & p, double scale)
{
    return Polygon<POINT>(p) /= scale;
}

template <class POINT>
inline Polygon<POINT> 
transpose(Polygon<POINT> const & p)
{
    Polygon<POINT> result(p.size());
    for(unsigned int i = 0; i < p.size(); ++i)
    {
        result.setPointUnsafe(i, POINT(p[i][1], p[i][0]));
    }
    return result;
}

template<class Point>
Point centroid(const Polygon<Point> &polygon)
{
    vigra_precondition(polygon.closed(),
                       "centroid() expects a closed polygon");
    double a = 0.0;
    TinyVector<double, 2> result;
    for(unsigned int i = 1; i < polygon.size(); ++i)
    {
        double pa = (polygon[i-1][0]*polygon[i][1] -
                     polygon[i-1][1]*polygon[i][0]);
        a += pa;
        result += (polygon[i-1] + polygon[i])*pa;
    }
    return result / (3.0*a);
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
    
    typename PointArray1::const_iterator begin = points.begin();
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
template<class Point, class T, class S, class Value>
void fillPolygon(Polygon<Point> const &p,
                 MultiArrayView<2, T, S> &output_image, 
                 Value value) 
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

#if 0

// the following sophisticated polygon resampling functions have no tests yet

/********************************************************************/

template<bool useMaxStep, class PointIterator, class TargetArray>
void simplifyPolygonHelper(
    const PointIterator &polyBegin, const PointIterator &polyEnd,
    TargetArray &simple, double epsilon,
    double maxStep2 = vigra::NumericTraits<double>::max())
{
    if(polyEnd - polyBegin <= 2)
        return; // no splitpoint necessary / possible

    PointIterator splitPos(polyEnd), lastPoint(polyEnd);
    --lastPoint;

    double maxDist = epsilon;

    // calculate normal of straight end point connection
    typename TargetArray::value_type
        straight(*lastPoint - *polyBegin);
    double straightLength2 = straight.squaredMagnitude();

    // search splitpoint
    if(straightLength2 > 1e-16)
    {
        typename TargetArray::value_type
            normal(straight[1], -straight[0]);

        normal /= std::sqrt(straightLength2);

        // iterate over points *between* first and last point:
        PointIterator it(polyBegin);
        for(++it; it != lastPoint; ++it)
        {
            double dist = fabs(dot(*it - *polyBegin, normal));
            if(dist > maxDist)
            {
                splitPos = it;
                maxDist = dist;
            }
        }
    }
    else
    {
        // start- and end-points identical?! -> look for most distant point
        PointIterator it(polyBegin);
        for(++it; it != lastPoint; ++it)
        {
            double dist = (*it - *polyBegin).magnitude();
            if(dist > maxDist)
            {
                splitPos = it;
                maxDist = dist;
            }
        }
    }

    if(useMaxStep && (straightLength2 > maxStep2) && (splitPos == polyEnd))
    {
        PointIterator it(polyBegin);
        ++it;
        double bestD2D = std::fabs((*it - *polyBegin).squaredMagnitude()
                                   - (*it - *lastPoint).squaredMagnitude());
        splitPos = it;
        for(++it; it != lastPoint; ++it)
        {
            double dist2Diff = std::fabs((*it - *polyBegin).squaredMagnitude()
                                         - (*it - *lastPoint).squaredMagnitude());
            if(dist2Diff < bestD2D)
            {
                bestD2D = dist2Diff;
                splitPos = it;
            }
        }
    }

    if(splitPos != polyEnd)
    {
        simplifyPolygonHelper<useMaxStep>(
            polyBegin, splitPos + 1, simple, epsilon, maxStep2);
        simple.push_back(*splitPos);
        simplifyPolygonHelper<useMaxStep>(
            splitPos, polyEnd, simple, epsilon, maxStep2);
    }
}

template<class PointArray>
void simplifyPolygon(
    const PointArray &poly, PointArray &simple, double epsilon)
{
    simple.push_back(poly[0]);
    simplifyPolygonHelper<false>(poly.begin(), poly.end(), simple, epsilon);
    simple.push_back(poly[poly.size()-1]);
}

template<class PointArray>
void simplifyPolygon(
    const PointArray &poly, PointArray &simple, double epsilon, double maxStep)
{
    simple.push_back(poly[0]);
    simplifyPolygonHelper<true>(poly.begin(), poly.end(),
                                simple, epsilon, maxStep*maxStep);
    simple.push_back(poly[poly.size()-1]);
}

/********************************************************************/

namespace detail {

template <class Point>
Point digitalLineIntersection(TinyVector<double, 3> const & l1, TinyVector<double, 3> const & l2)
{
    double d = l1[0]*l2[1] - l1[1]*l2[0];
    return Point((l1[1]*l2[2] - l1[2]*l2[1]) / d, (l1[2]*l2[0] - l1[0]*l2[2]) / d);
}

} // namespace detail

template<class PointArray>
void simplifyPolygonDigitalLine(
    const PointArray &poly, PointArray &simple, int connectivity)
{
    typedef typename PointArray::value_type Point;

    int size = poly.size();
    if(size <= 2)
    {
        simple = poly;
        return;
    }

    vigra_precondition(connectivity == 4 || connectivity == 8,
       "simplifyPolygonDigitalLine(): connectivity must be 4 or 8.");

    bool isOpenPolygon = poly[size-1] != poly[0];

    ArrayVector<TinyVector<double, 3> > lines;
    Point l1 = poly[0],
          r1 = l1,
          l2 = poly[1],
          r2 = l2;
    double a = l2[1] - l1[1],
           b = l1[0] - l2[0],
           c = -a*l2[0] - b*l2[1];
    for(int k=2; k < size; ++k)
    {
        double ab = (connectivity == 4)
                        ? std::fabs(a) + std::fabs(b)
                        : std::max(std::fabs(a), std::fabs(b));
        double d = a*poly[k][0] + b*poly[k][1] + c;
        if(d < -1.0 || d > ab)
        {
            // finish current segment
            c = (c - a*r2[0] - b*r2[1]) / 2.0;
            lines.push_back(TinyVector<double, 3>(a, b, c));
            // initialize new segment
            l1 = poly[k-1];
            r1 = l1;
            l2 = poly[k];
            r2 = l2;
            a = l2[1] - l1[1];
            b = l1[0] - l2[0];
            c = -a*l2[0] - b*l2[1];
        }
        else if(d <= 0.0)
        {
            l2 = poly[k];
            if(d < 0.0)
            {
                r1 = r2;
                a = l2[1] - l1[1];
                b = l1[0] - l2[0];
                c = -a*l2[0] - b*l2[1];
            }
        }
        else if(d >= ab - 1.0)
        {
            r2 = poly[k];
            if(d > ab - 1.0)
            {
                l1 = l2;
                a = r2[1] - r1[1];
                b = r1[0] - r2[0];
                c = -a*l2[0] - b*l2[1];
            }
        }
    }

    c = (c - a*r2[0] - b*r2[1]) / 2.0;
    lines.push_back(TinyVector<double, 3>(a, b, c));
    int segments = lines.size();

    if(isOpenPolygon)
        simple.push_back(poly[0]);
    else
        simple.push_back(detail::digitalLineIntersection<Point>(lines[0], lines[segments-1]));

    for(int k=1; k<segments; ++k)
        simple.push_back(detail::digitalLineIntersection<Point>(lines[k-1], lines[k]));

    if(isOpenPolygon)
        simple.push_back(poly[size-1]);
    else
        simple.push_back(detail::digitalLineIntersection<Point>(lines[0], lines[segments-1]));
}

/********************************************************************/

template<class PointArray>
void resamplePolygon(
    const PointArray &poly, PointArray &simple, double desiredPointDistance)
{
    typedef typename PointArray::value_type Point;

    int size = poly.size();
    bool isOpenPolygon = !poly.closed();

    ArrayVector<double> arcLength;
    poly.arcLengthList(arcLength);
    int segmentCount = int(std::ceil(arcLength[size-1] / desiredPointDistance));
    if(segmentCount < 2)
    {
        simple.push_back(poly[0]);
        if(!isOpenPolygon)
        {
            Point p = poly[1];
            double dist = (p - poly[0]).magnitude();
            for(int k=2; k < size-1; ++k)
            {
                double d = (poly[k] - poly[0]).magnitude();
                if(d > dist)
                {
                    dist = d;
                    p = poly[k];
                }
            }
            simple.push_back(p);
        }
        simple.push_back(poly[size-1]);
        return;
    }

    for(int k=0; k<size; ++k)
        arcLength[k] *= segmentCount / arcLength[size-1];

    ArrayVector<Point> integrals(segmentCount+1, Point(0.0, 0.0));
    Point p1 = poly[0];
    double t1 = 0.0;
    int l = 1;
    for(int k=1; k<size; ++k)
    {
        double d = arcLength[k];
        while(d >= l && l <= segmentCount)
        {
            double t2 = 1.0;
            double dt = t2 - t1;
            Point p2 = poly.interpolate(k-1, (l - arcLength[k-1]) / (d - arcLength[k-1]));
            Point sum1 = 0.5 * dt * (p1 + p2);
            Point sumt = dt / 6.0 * (p1*(2.0*t1+t2) + p2*(t1+2.0*t2));
            integrals[l-1] += sum1 - sumt;
            integrals[l] += sumt;
            if(isOpenPolygon && l==1)
            {
                integrals[0] +=  poly[0] - integrals[1];
            }
            p1 = p2;
            t1 = 0.0;
            ++l;
            if(isOpenPolygon && l==segmentCount)
            {
                integrals[segmentCount] += integrals[segmentCount-1];
            }
        }
        if(d < l && l <= segmentCount)
        {
            double t2 = std::fmod(d, 1.0);
            double dt = t2 - t1;
            Point p2 = poly[k];
            Point sum1 = 0.5 * dt * (p1 + p2);
            Point sumt = dt / 6.0 * (p1*(2.0*t1+t2) + p2*(t1+2.0*t2));
            integrals[l-1] += sum1 - sumt;
            integrals[l] += sumt;
            p1 = p2;
            t1 = t2;
        }
    }

    if(isOpenPolygon)
    {
        integrals[segmentCount] += poly[size-1] - integrals[segmentCount-1];
        integrals[1] -= integrals[0] / 6.0;
        integrals[segmentCount-1] -= integrals[segmentCount] / 6.0;

        ArrayVector<double> g(segmentCount);
        double b = 2.0 / 3.0;
        simple.push_back(poly[0]);
        simple.push_back(integrals[1] / b);
        for(int k=2; k<segmentCount; ++k)
        {
            g[k] = 1.0 / 6.0 / b;
            b = 2.0 / 3.0 - g[k] / 6.0;
            simple.push_back((integrals[k] - simple[k-1] / 6.0) / b);
        }
        for(int k = segmentCount-2; k >= 1; --k)
        {
            simple[k] -= g[k+1] * simple[k+1];
        }

        simple.push_back(poly[size-1]);
    }
    else
    {
        integrals[0] += integrals[segmentCount];

        int initializationSteps = std::min(segmentCount, 5);
        ArrayVector<Point> p(segmentCount+2*initializationSteps);
        double b = 0.6220084679281461,
               g = 0.26794919243112275;
        p[0] = integrals[0] / b;

        for(int k=1; k<segmentCount+2*initializationSteps; ++k)
        {
            p[k] = (integrals[k % segmentCount] - p[k-1] / 6.0) / b;
        }
        for(int k = segmentCount+2*initializationSteps-2; k >= initializationSteps; --k)
        {
            p[k] -= g * p[k+1];
        }

        for(int k=segmentCount; k<segmentCount+initializationSteps; ++k)
            simple.push_back(p[k]);
        for(int k=initializationSteps; k<=segmentCount; ++k)
            simple.push_back(p[k]);
    }
}

/********************************************************************/

template<class PointArray>
void resamplePolygonLinearInterpolation(
    const PointArray &poly, PointArray &simple, double desiredPointDistance)
{
    int size = poly.size();
    if(size <= 2)
    {
        simple = poly;
        return;
    }

    ArrayVector<double> arcLengths;
    poly.arcLengthList(arcLengths);

    int steps = int(std::ceil(arcLengths[size-1] / desiredPointDistance));
    double newStep = arcLengths[size-1] / steps;

    simple.push_back(poly[0]);
    int l = 1;
    for(int k=1; k<steps; ++k)
    {
        double currentArcLength = k*newStep;
        for(; l < size; ++l)
        {
            if(arcLengths[l] >= currentArcLength)
                break;
        }
        double o = (arcLengths[l] - currentArcLength) / (arcLengths[l] - arcLengths[l-1]);
        simple.push_back(o*poly[l-1] + (1.0-o)*poly[l]);
    }
    simple.push_back(poly[size-1]);
}

/********************************************************************/

template<class PointArray>
void resamplePolygonExponentialFilter(
    const PointArray &poly, PointArray &simple, double scale, double desiredPointDistance)
{
    int size = poly.size();
    if(size <= 2)
    {
        simple = poly;
        return;
    }

    bool isOpenPolygon = !poly.closed();

    typedef typename PointArray::value_type Point;
    ArrayVector<Point> pforward(size), pbackward(size);
    ArrayVector<double> wforward(size), wbackward(size), weights(size-1);
    for(int k=0; k < size - 1; ++k)
        weights[k] = std::exp(-(poly[k] - poly[k+1]).magnitude()/scale);

    // init recursion with cyclic boundary conditions
    Point p = poly[0];
    double w = 1.0;
    for(int k=1; k < size; ++k)
    {
        p = poly[k] + weights[k-1]*p;
        w = 1.0 + weights[k-1]*w;
    }
    pforward[0] = p;
    wforward[0] = w;

    p = poly[size-1];
    w = 1.0;
    for(int k=size-2; k>=0; --k)
    {
        p = poly[k] + weights[k]*p;
        w = 1.0 + weights[k]*w;
    }
    pbackward[size-1] = p;
    wbackward[size-1] = w;

    if(isOpenPolygon)
    {
        // change initialization into anti-reflective boundary conditions for open polygons
        std::swap(wbackward[size-1], wforward[0]);
        std::swap(pbackward[size-1], pforward[0]);
        pforward[0] = 2.0*wforward[0]*poly[0] - pforward[0];
        pbackward[size-1] = 2.0*wbackward[size-1]*poly[size-1] - pbackward[size-1];
    }

    // forward and backward pass of the recursive filter
    for(int k=1; k < size; ++k)
    {
        pforward[k] = poly[k] + weights[k-1]*pforward[k-1];
        wforward[k] = 1.0 + weights[k-1]*wforward[k-1];
    }
    for(int k=size-2; k >= 0; --k)
    {
        pbackward[k] = poly[k] + weights[k]*pbackward[k+1];
        wbackward[k] = 1.0 + weights[k]*wbackward[k+1];
    }

    // measure the arc length of the new polygon (after possible shrinkage)
    p = (pforward[0]+weights[0]*pbackward[1]) / (wforward[0] + weights[0]*wbackward[1]);
    simple.push_back(p);

    Point pend = isOpenPolygon
                   ? (weights[size-2]*pforward[size-2]+pbackward[size-1]) /
                     (weights[size-2]*wforward[size-2] + wbackward[size-1])
                   : p;

    ArrayVector<double> arcLength;
    double length = 0.0;
    arcLength.push_back(length);
    for(int k=1; k<size-1; ++k)
    {
        Point pc = (pforward[k]+weights[k]*pbackward[k+1]) / (wforward[k] + weights[k]*wbackward[k+1]);
        length += (pc - p).magnitude();
        arcLength.push_back(length);
        p = pc;
    }
    length += (p-pend).magnitude();
    arcLength.push_back(length);

//    alternative: use the arc lenth of the original polygon
//    poly.arcLengthList(arcLength);

    int steps = int(std::floor(arcLength[size-1] / desiredPointDistance+0.5));
    double newStep = arcLength[size-1] / steps;

    int l = 1;
    for(int k=1; k < steps; ++k)
    {
        double currentArcLength = k*newStep;
        for(; l < size; ++l)
        {
            if(arcLength[l] >= currentArcLength)
                break;
        }
        double w = weights[l-1];
        double o = (arcLength[l] - currentArcLength) / (arcLength[l] - arcLength[l-1]);
        double wl = std::pow(w, 1.0-o);
        double wr = std::pow(w, o);
        simple.push_back((wl*pforward[l-1]+wr*pbackward[l]) / (wl*wforward[l-1] + wr*wbackward[l]));
    }
    simple.push_back(pend);
}

/********************************************************************/

namespace detail {

template<class ArcLengthList, class PointList>
typename PointList::value_type
singleGaussianConvolvePolygonReflective(
    const ArcLengthList &arcLengthList,
    const PointList &pointList,
    int i, double arcLengthPos,
    const Gaussian<double> &g)
{
    typedef typename PointList::value_type ValueType;

    ValueType sum(vigra::NumericTraits<ValueType>::zero());
    double norm = 0.0;

    int size = arcLengthList.size(),
        lastIndex = size - 1;
    double reflectLength = 2.0*arcLengthList[lastIndex];

    ValueType reflectPoint = 2.0*pointList[lastIndex];
    for(int j = i; true; ++j)
    {
        int k = (j > lastIndex)
                ? 2*lastIndex - j
                : j;
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j > lastIndex)
        {
            pos = reflectLength - pos;
            point = reflectPoint - point;
        }
        double diff = pos - arcLengthPos;
        if(diff > g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    reflectPoint = 2.0*pointList[0];
    for(int j = i - 1; true; --j)
    {
        int k = std::abs(j);
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j < 0)
        {
            pos = -pos;
            point = reflectPoint - point;
        }
        double diff = pos - arcLengthPos;
        if(diff < -g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    return sum / norm;
}

template<class ArcLengthList, class PointList>
typename PointList::value_type
singleGaussianConvolvePolygonCyclic(
    const ArcLengthList &arcLengthList,
    const PointList &pointList,
    int i, double arcLengthPos,
    const Gaussian<double> &g)
{
    typedef typename PointList::value_type ValueType;

    ValueType sum(vigra::NumericTraits<ValueType>::zero());
    double norm = 0.0;

    int size = arcLengthList.size() - 1,
        lastIndex = size - 1;
    double totalLength = arcLengthList[size];

    for(int j = i; true; ++j)
    {
        int k = j % size;
        double pos = j > lastIndex
                        ? arcLengthList[k] + totalLength
                        : arcLengthList[k];
        ValueType point = pointList[k];
        double diff = pos - arcLengthPos;
        if(diff > g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    for(int j = i - 1; true; --j)
    {
        int k = (j + size) % size;
        double pos = j < 0
                       ? arcLengthList[k] - totalLength
                       : arcLengthList[k];
        ValueType point = pointList[k];
        double diff = pos - arcLengthPos;
        if(diff < -g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    return sum / norm;
}

} // namespace detail

template<class PointArray>
void resamplePolygonGaussianFilter(
    const PointArray &poly, PointArray &simple, double scale, double desiredPointDistance)
{
    int size = poly.size();
    if(size <= 2)
    {
        simple = poly;
        return;
    }

    ArrayVector<double> arcLengths;
    poly.arcLengthList(arcLengths);

    Gaussian<double> g(scale);

    vigra_precondition(arcLengths[size-1] > g.radius(),
        "resamplePolygonGaussianFilter(): Filter longer than polygon.");

    bool isOpenPolygon = !poly.closed();

    int steps = int(std::ceil(arcLengths[size-1] / desiredPointDistance));
    double newStep = arcLengths[size-1] / steps;

    int l = 0;
    for(int k=0; k<steps; ++k)
    {
        double currentArcLength = k*newStep;
        for(; l < size; ++l)
        {
            if(arcLengths[l] >= currentArcLength)
                break;
        }
        if(isOpenPolygon)
            simple.push_back(detail::singleGaussianConvolvePolygonReflective(arcLengths, poly, l, currentArcLength, g));
        else
            simple.push_back(detail::singleGaussianConvolvePolygonCyclic(arcLengths, poly, l, currentArcLength, g));
    }
    if(isOpenPolygon)
        simple.push_back(detail::singleGaussianConvolvePolygonReflective(arcLengths, poly, size-1, arcLengths[size-1], g));
    else
        simple.push_back(simple[0]);
}

/********************************************************************/

namespace detail {

template<class Point>
Point spline3Integral(Point const & p1, Point const & p2, double t1, double t2)
{
    StaticPolynomial<5, double> p[2];
    p[0][0] = p[1][0] = 0.0;
    if(t1 >= 1.0)
    {
        return (t1 - t2) / 120.0 *
               (p1 * (-80.0 + t1*(80.0 + 2.0*t2*(t2 - 10.0) + t1*(3.0*t2 - 30.0 + 4.0*t1)) + t2*(40.0 + t2*(t2 - 10.0))) +
                p2 * (-80.0 + t1*(40.0 + t2*(3.0*t2 - 20.0) + t1*(2.0*t2 - 10.0 + t1)) + t2*(80.0 + t2*(4.0*t2 - 30.0))));
    }
    else
    {
        return (t2 - t1) / 120.0 *
               (p1 * (40.0 + t1*(2.0*t2*(3.0*t2 - 10.0) + t1*(9.0*t2 - 30.0 + 12.0*t1)) + t2*t2*(3.0*t2 - 10.0)) +
                p2 * (40.0 + t1*(t2*(9.0*t2 - 20.0) + t1*(6.0*t2 - 10.0 + 3.0*t1)) + t2*t2*(12.0*t2 - 30.0)));
    }
}

template<class ArcLengthList, class PointList>
typename PointList::value_type
singleSpline3ConvolvePolygon(
    const ArcLengthList &arcLengthList,
    const PointList &pointList,
    int left, int center, int right)
{
    typedef typename PointList::value_type ValueType;

    ValueType sum(vigra::NumericTraits<ValueType>::zero());
    double arcLengthPos = arcLengthList[center];
    for(int j = center + 1; j <= right; ++j)
    {
        double t1 = arcLengthList[j-1] - arcLengthPos,
               t2 = arcLengthList[j] - arcLengthPos;
        sum += spline3Integral(pointList[j-1], pointList[j], t1, t2);
    }
    for(int j = center - 1; j >= left; --j)
    {
        double t1 = arcLengthPos - arcLengthList[j+1],
               t2 = arcLengthPos - arcLengthList[j];
        sum -= spline3Integral(-pointList[j+1], -pointList[j], t1, t2);
    }

    return sum;
}

} // namespace detail

template<class PointArray>
void polygonSplineControlPoints(
    const PointArray &poly, PointArray &splinePoints, int segmentCount)
{
    typedef typename PointArray::value_type Point;

    int size = poly.size();
    vigra_precondition(size >= 4,
        "polygonSplineControlPoints(): Polygon must have at least 4 points.");

    bool isOpenPolygon = !poly.closed();

    ArrayVector<double> arcLength;
    poly.arcLengthList(arcLength);
    double totalLength = segmentCount / arcLength[size-1];
    for(int k=0; k<size; ++k)
        arcLength[k] *= totalLength;

    PointArray augmentedPoly;
    augmentedPoly.push_back(poly[0]);

    ArrayVector<double> augmentedArcLength;
    augmentedArcLength.push_back(0.0);

    ArrayVector<int> splineIndices(segmentCount + 1);
    splineIndices[0] = 0;
    int l = 1;
    for(int k=1; k<size-1; ++k)
    {
        double d = arcLength[k];
        while(d > l)
        {
            augmentedPoly.push_back(poly.interpolate(k-1, (l - arcLength[k-1]) / (d - arcLength[k-1])));
            augmentedArcLength.push_back(l);
            splineIndices[l] = augmentedPoly.size()-1;
            ++l;
        }
        augmentedPoly.push_back(poly[k]);
        augmentedArcLength.push_back(d);
        if(d == l)
        {
            splineIndices[l] = augmentedPoly.size()-1;
            ++l;
        }
    }
    augmentedPoly.push_back(poly[size-1]);
    augmentedArcLength.push_back(segmentCount);
    splineIndices[segmentCount] = augmentedPoly.size()-1;
    size = augmentedPoly.size();

    ArrayVector<Point> integrals(segmentCount+1);
    if(isOpenPolygon)
    {
        integrals[0] = augmentedPoly[0];
        PointArray reflectedPoly;
        Point reflectPoint = 2.0*poly[0];
        ArrayVector<double> reflectedArcLength;
        for(int k=-splineIndices[1]; k <= splineIndices[3]; ++k)
        {
            if(k < 0)
            {
                reflectedPoly.push_back(reflectPoint - augmentedPoly[-k]);
                reflectedArcLength.push_back(-augmentedArcLength[-k]);
            }
            else
            {
                reflectedPoly.push_back(augmentedPoly[k]);
                reflectedArcLength.push_back(augmentedArcLength[k]);
            }
        }
        integrals[1] = detail::singleSpline3ConvolvePolygon(reflectedArcLength, reflectedPoly,
                                    0, 2*splineIndices[1], splineIndices[1] + splineIndices[3]);

        reflectPoint = 2.0*augmentedPoly[size-1];
        for(int k=size-2; k>=splineIndices[segmentCount-1]; --k)
        {
            augmentedPoly.push_back(reflectPoint - augmentedPoly[k]);
            augmentedArcLength.push_back(2.0*segmentCount - augmentedArcLength[k]);
        }
        integrals[segmentCount-1] = detail::singleSpline3ConvolvePolygon(augmentedArcLength, augmentedPoly,
               splineIndices[segmentCount-3], splineIndices[segmentCount-1],
               2*splineIndices[segmentCount] - splineIndices[segmentCount-1]);
        integrals[segmentCount] = augmentedPoly[size-1];
    }
    else
    {
        PointArray wrappedPoly;
        ArrayVector<double> wrappedArcLength;
        for(int k=splineIndices[segmentCount-1]; k < splineIndices[segmentCount]; ++k)
        {
            wrappedPoly.push_back(augmentedPoly[k]);
            wrappedArcLength.push_back(augmentedArcLength[k] - segmentCount);
        }
        int indexShift = wrappedPoly.size();
        for(int k=0; k <= splineIndices[3]; ++k)
        {
            wrappedPoly.push_back(augmentedPoly[k]);
            wrappedArcLength.push_back(augmentedArcLength[k]);
        }
        integrals[1] = detail::singleSpline3ConvolvePolygon(wrappedArcLength, wrappedPoly,
                             0, splineIndices[1] + indexShift, splineIndices[3] + indexShift);

        for(int k=1; k <= splineIndices[2]; ++k)
        {
            augmentedPoly.push_back(augmentedPoly[k]);
            augmentedArcLength.push_back(segmentCount + augmentedArcLength[k]);
        }
        integrals[segmentCount-1] = detail::singleSpline3ConvolvePolygon(augmentedArcLength, augmentedPoly,
               splineIndices[segmentCount-3], splineIndices[segmentCount-1],
               splineIndices[segmentCount] + splineIndices[1]);
        integrals[0] = detail::singleSpline3ConvolvePolygon(augmentedArcLength, augmentedPoly,
               splineIndices[segmentCount-2], splineIndices[segmentCount],
               splineIndices[segmentCount] + splineIndices[2]);
    }

    for(int k=2; k <= segmentCount-2; ++k)
        integrals[k] = detail::singleSpline3ConvolvePolygon(augmentedArcLength, augmentedPoly,
                                    splineIndices[k-2], splineIndices[k], splineIndices[k+2]);

    BSpline<7, double> spline7;
    if(isOpenPolygon)
    {
        int solutionSize = segmentCount + 1;
        Matrix<double> m(solutionSize, solutionSize),
                    rhs(solutionSize, 2),
                    solution(solutionSize, 2);
        for(int k=0; k<solutionSize; ++k)
        {
            for(int l=-3; l<=3; ++l)
            {
                if(k + l < 0)
                {
                    m(k, 0) += 2.0*spline7(l);
                    m(k, abs(k+l)) -= spline7(l);
                }
                else if(k + l >= solutionSize)
                {
                    m(k, solutionSize - 1) += 2.0*spline7(l);
                    m(k, 2*solutionSize - k - l - 2) -= spline7(l);
                }
                else
                    m(k,k+l) += spline7(l);
            }
            rhs(k, 0) = integrals[k][0];
            rhs(k, 1) = integrals[k][1];
        }

        linearSolve(m, rhs, solution);

        for(int k=0; k<solutionSize; ++k)
        {
            splinePoints.push_back(Point(solution(k,0), solution(k,1)));
        }
        splinePoints.push_back(2.0*splinePoints[solutionSize-1] - splinePoints[solutionSize-2]);
        splinePoints.insert(splinePoints.begin(), 2.0*splinePoints[0] - splinePoints[1]);
    }
    else
    {
        int solutionSize = segmentCount;
        Matrix<double> m(solutionSize, solutionSize),
                    rhs(solutionSize, 2),
                    solution(solutionSize, 2);
        for(int k=0; k<solutionSize; ++k)
        {
            for(int l=-3; l<=3; ++l)
                m(k, (k+l+solutionSize) % solutionSize) = spline7(l);
            rhs(k, 0) = integrals[k][0];
            rhs(k, 1) = integrals[k][1];
        }
        linearSolve(m, rhs, solution);

        for(int k=0; k<solutionSize; ++k)
        {
            splinePoints.push_back(Point(solution(k,0), solution(k,1)));
        }
        splinePoints.push_back(splinePoints[0]);
    }
}

#endif

//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
