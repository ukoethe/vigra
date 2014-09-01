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

        /**
         * Tests whether the given point lies within this polygon.
         * Requires that this polygon is closed.

         * The result of testing points which lie directly on the
         * polylines (or are incident with the support points) is
         * undefined.  (ATM, the implementation uses half-open
         * intervals, so points on the left/top border are included,
         * in contrast to the ones on the right/bottom.)
         */
    bool contains(const_reference point) const
    {
        vigra_precondition(closed(),
                           "Polygon::contains() requires polygon to be closed!");
        int result = 0;
        bool above = (*this)[0][1] < point[1];
        for(unsigned int i = 1; i < size(); ++i)
        {
            bool now = (*this)[i][1] < point[1];
            if(now != above)
            {
                typename Point::value_type intersectX =
                    (*this)[i-1][0] +
                    ((*this)[i][0] - (*this)[i-1][0]) *
                    (point[1]      - (*this)[i-1][1]) /
                    ((*this)[i][1] - (*this)[i-1][1]);
                if(intersectX < point[0])
                    ++result;

                above = now;
            }
        }
        return (result % 2) != 0;
    }

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

/********************************************************************/

/** \brief Create a polygon from the interpixel contour of a labeled region.

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
/*
 * Left hand on the wall contour extraction
 * the label of the anchor point marks object, everything else is background
 * anchor_point is the first point inside the object when traversing in scan-order
 * contour_points contains a half-integer point for each section of the wall
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

/*
 * Scan line filling, works on both convex and non convex polygons
 * The contour points are outside the region to draw
 * Adjacent polygon points must be next to each other in p too.
 */
template<class PointArray, class T, class S>
void fillPolygon(PointArray const &p,
                 MultiArrayView<2, T, S> &output_image, 
                 T value) 
{
    typedef typename PointArray::value_type Point;
    
    std::vector<Point> contour_points;

    auto ip = p.begin();

    for (; ip < p.end() - 1; ++ip) {
        contour_points.push_back(*ip);
        pushLinePoints(*ip, *(ip + 1), contour_points);
    }
    contour_points.push_back(*(p.end() - 1));
    pushLinePoints(p.back(), p.front(), contour_points);

    sort(contour_points.begin(), contour_points.end(), detail::pointYXOrdering<Point>);

    auto points_iterator = contour_points.begin();

    float min_y = (contour_points.front())[1];
    float max_y = (contour_points.back())[1];

    while (points_iterator != contour_points.end()) {
        float y = (*points_iterator)[1];
        if ((y > min_y) && (y < max_y)) { // We are inside the polygon
            Point current = *points_iterator;
            float min_x = current[0];
            current[0] = ceil(current[0]);
            current[1] = ceil(current[1]);

            bool drawing = true;

            while ((*points_iterator)[1] == y) {
                ++points_iterator;

                Point endPoint = *points_iterator;
                float max_x = endPoint[0];

                // Draw the scan line

                for (; current[0] < endPoint[0]; ++current[0]) {
                    if ((current[0] > min_x) && (current[0] < max_x)) {
                        if (drawing) {
                            output_image(current[0], current[1]) = value;
                        }
                    }
                }

                drawing = !drawing;
            }

        } else { // We are on an edge, just continue
            ++points_iterator;
        }
    }
}


//@}

} // namespace vigra

#endif /* VIGRA_POLYGON_HXX */
