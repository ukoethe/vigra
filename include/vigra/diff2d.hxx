/************************************************************************/
/*                                                                      */
/*                  Copyright 1998-2003 by Hans Meine                   */
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

#ifndef VIGRA_DIFF2D_HXX
#define VIGRA_DIFF2D_HXX

#include <cmath> // for sqrt()
#include <iosfwd>
#include "vigra/config.hxx"
#include "vigra/iteratortags.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/iteratoradapter.hxx"
#include "vigra/tuple.hxx"

namespace vigra {

template <class Diff>
class Diff2DConstRowIteratorPolicy
{
  public:
    typedef Diff                            BaseType;
    typedef Diff                            value_type;
    typedef typename Diff::MoveX            difference_type;
    typedef Diff const &                    reference;
    typedef Diff                            index_reference;
    typedef Diff const *                    pointer;
    typedef std::random_access_iterator_tag iterator_category;

    static void initialize(BaseType &) {}

    static reference dereference(BaseType const & d)
        { return d; }

    static index_reference dereference(BaseType d, difference_type n)
    {
        d.x += n;
        return d;
    }

    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.x == d2.x; }

    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.x < d2.x; }

    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return d1.x - d2.x; }

    static void increment(BaseType & d)
        { ++d.x; }

    static void decrement(BaseType & d)
        { --d.x; }

    static void advance(BaseType & d, difference_type n)
        { d.x += n; }
};

template <class Diff>
class Diff2DConstColumnIteratorPolicy
{
  public:
    typedef Diff                            BaseType;
    typedef Diff                            value_type;
    typedef typename Diff::MoveY            difference_type;
    typedef Diff const &                    reference;
    typedef Diff                            index_reference;
    typedef Diff const *                    pointer;
    typedef std::random_access_iterator_tag iterator_category;

    static void initialize(BaseType & /*d*/) {}

    static reference dereference(BaseType const & d)
        { return d; }

    static index_reference dereference(BaseType d, difference_type n)
    {
        d.y += n;
        return d;
    }

    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.y == d2.y; }

    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.y < d2.y; }

    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return d1.y - d2.y; }

    static void increment(BaseType & d)
        { ++d.y; }

    static void decrement(BaseType & d)
        { --d.y; }

    static void advance(BaseType & d, difference_type n)
        { d.y += n; }
};

/** \addtogroup RangesAndPoints Two-dimensional Ranges and Points

    Specify a 2D position, extent, or rectangle.
*/
//@{

/********************************************************/
/*                                                      */
/*                      Diff2D                          */
/*                                                      */
/********************************************************/

/** \brief Two dimensional difference vector.

    This class acts primarily as a difference vector for specifying
    pixel coordinates and region sizes. In addition, Diff2D fulfills
    the requirements of an \ref ImageIterator, so that it can be used to
    simulate an image whose pixels' values equal their coordinates. This
    secondary usage is explained on page \ref CoordinateIterator.

    Standard usage as a difference vector is mainly needed in the context
    of images. For example, Diff2D may be used as an index for <TT>operator[]</TT>:

    \code
    vigra::Diff2D location(...);

    value = image[location];
    \endcode

    This is especially important in connection with accessors, where the
    offset variant of <TT>operator()</TT> takes only one offset object:

    \code
    // accessor(iterator, dx, dy); is not allowed
    value = accessor(iterator, vigra::Diff2D(dx, dy));
    \endcode


    Diff2D is also returned by <TT>image.size()</TT>, so that we can create
    new images by calculating their size using Diff2D's arithmetic
    functions:

    \code
    // create an image that is 10 pixels smaller in each direction
    Image new_image(old_image.size() - Diff2D(10,10));
    \endcode

    <b>\#include</b> "<a href="diff2d_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra
*/
class Diff2D
{
  public:
        /** The iterator's value type: a coordinate.
        */
    typedef Diff2D PixelType;

        /** The iterator's value type: a coordinate.
        */
    typedef Diff2D value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef Diff2D const &       reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef Diff2D const *       pointer;

        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               difference_type;

        /** the iterator tag (image traverser)
        */
    typedef image_traverser_tag  iterator_category;

        /** The associated row iterator.
        */
    typedef IteratorAdaptor<Diff2DConstRowIteratorPolicy<Diff2D> >    row_iterator;

        /** The associated column iterator.
        */
   typedef IteratorAdaptor<Diff2DConstColumnIteratorPolicy<Diff2D> > column_iterator;

        /** type of the iterator's x-navigator
        */
    typedef int MoveX;
        /** type of the iterator's y-navigator
        */
    typedef int MoveY;


        /** Default Constructor. Init iterator at position (0,0)
        */
    Diff2D()
    : x(0), y(0)
    {}

        /** Construct at given position.
        */
    Diff2D(int ax, int ay)
    : x(ax), y(ay)
    {}

        /** Copy Constructor.
        */
    Diff2D(Diff2D const & v)
    : x(v.x), y(v.y)
    {}

        /** Copy Assigment.
        */
    Diff2D & operator=(Diff2D const & v)
    {
        if(this != &v)
        {
            x = v.x;
            y = v.y;
        }
        return *this;
    }

        /** Unary negation.
        */
    Diff2D operator-() const
    {
        return Diff2D(-x, -y);
    }

        /** Increase coordinate by specified offset.
        */
    Diff2D & operator+=(Diff2D const & offset)
    {
        x += offset.x;
        y += offset.y;
        return *this;
    }

        /** Decrease coordinate by specified vector.
        */
    Diff2D & operator-=(Diff2D const & offset)
    {
        x -= offset.x;
        y -= offset.y;
        return *this;
    }

       /** Create vector by scaling by factor.
        */
    Diff2D & operator*=(int factor)
    {
        x *= factor;
        y *= factor;
        return *this;
    }

       /** Create vector by scaling by factor.
        */
    Diff2D & operator*=(double factor)
    {
        x = (int)(x * factor);
        y = (int)(y * factor);
        return *this;
    }

       /** Create vector by scaling by 1/factor.
        */
    Diff2D & operator/=(int factor)
    {
        x /= factor;
        y /= factor;
        return *this;
    }

       /** Create vector by scaling by 1/factor.
        */
    Diff2D & operator/=(double factor)
    {
        x = (int)(x / factor);
        y = (int)(y / factor);
        return *this;
    }

       /** Create vector by scaling by factor.
        */
    Diff2D operator*(int factor) const
    {
        return Diff2D(x * factor, y * factor);
    }

       /** Create vector by scaling by factor.
        */
    Diff2D operator*(double factor) const
    {
        return Diff2D((int)(x * factor), (int)(y * factor));
    }

       /** Create vector by scaling by 1/factor.
        */
    Diff2D operator/(int factor) const
    {
        return Diff2D(x / factor, y / factor);
    }

       /** Create vector by scaling by 1/factor.
        */
    Diff2D operator/(double factor) const
    {
        return Diff2D((int)(x / factor), (int)(y / factor));
    }

        /** Calculate length of difference vector.
        */
    int squaredMagnitude() const
    {
        return x*x + y*y;
    }

        /** Calculate length of difference vector.
        */
    double magnitude() const
    {
        return VIGRA_CSTD::sqrt((double)squaredMagnitude());
    }

        /** Equality.
        */
    bool operator==(Diff2D const & r) const
    {
        return (x == r.x) && (y == r.y);
    }

        /** Inequality.
        */
    bool operator!=(Diff2D const & r) const
    {
        return (x != r.x) || (y != r.y);
    }

        /** Used for both access to the current x-coordinate \em and
            to specify that an iterator navigation command is to be
            applied in x-direction. <br>
            usage:  <TT> x = diff2d.x </TT> (use \p Diff2D::x  as component of difference vector) <br>
            or <TT>&nbsp; ++diff.x &nbsp; </TT> (use Diff2D as iterator, move right)
         */
    int x;
        /** Used for both access to the current y-coordinate \em and
            to specify that an iterator navigation command is to be
            applied in y-direction. <br>
            usage:  <TT> y = diff2d.y </TT> (use \p Diff2D::y as component of difference vector) <br>
            or <TT>&nbsp; ++diff.y &nbsp; </TT> (use Diff2D as iterator, move right)
        */
    int y;

        /** Access current coordinate.
        */
    reference operator*() const
    {
        return *this;
    }

        /** Read coordinate at an offset.
        */
    index_reference operator()(int const & dx, int const & dy) const
    {
        return Diff2D(x + dx, y + dy);
    }

        /** Read coordinate at an offset.
        */
    index_reference operator[](Diff2D const & offset) const
    {
        return Diff2D(x + offset.x, y + offset.y);
    }

        /** Read vector components.
        */
    int operator[](int index) const
    {
        return (&x)[index];
    }

        /** Access current coordinate.
        */
    pointer operator->() const
    {
        return this;
    }

        /** Get a row iterator at the current position.
        */
    row_iterator rowIterator() const
        { return row_iterator(*this); }

        /** Get a column iterator at the current position.
        */
    column_iterator columnIterator() const
        { return column_iterator(*this); }
};


template <>
struct IteratorTraits<Diff2D >
{
    typedef Diff2D                               Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardConstValueAccessor<Diff2D>   DefaultAccessor;
    typedef StandardConstValueAccessor<Diff2D>   default_accessor;
    typedef VigraTrueType                        hasConstantStrides;

};


/********************************************************/
/*                                                      */
/*                      Size2D                          */
/*                                                      */
/********************************************************/

/** \brief Two dimensional size object.

    Specializes \ref Diff2D for the specification of a 2-dimensional
    extent, in contrast to a point or position (for the latter
    use \ref Point2D).

    \code
    // create an image that is 10 pixels squared
    Image new_image(Size2D(10,10));
    \endcode

    <b>\#include</b> "<a href="diff2d_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra
*/
class Size2D : public Diff2D
{
public:
        /** Default Constructor. Init point at position (0,0)
        */
    Size2D()
    {}

        /** Construct point at given position.
        */
    Size2D(int width, int height)
    : Diff2D(width, height)
    {}

        /** Copy Constructor.
        */
    Size2D(Size2D const & v)
    : Diff2D(v)
    {}

        /** Explicit conversion Constructor.
        */
    explicit Size2D(Diff2D const & v)
    : Diff2D(v)
    {}

        /** Query the width.
         */
    int width() const
    {
        return x;
    }

        /** Query the height.
         */
    int height() const
    {
        return y;
    }

        /** Returns width()*height(), the area of a rectangle of this size.
         */
    int area() const
    {
        return width()*height();
    }

        /** Copy Assigment.
        */
    Size2D & operator=(Diff2D const & v)
    {
        return static_cast<Size2D &>(Diff2D::operator=(v));
    }

        /** Unary negation.
        */
    Size2D operator-() const
    {
        return Size2D(-x, -y);
    }

        /** Increase size by specified offset.
        */
    Size2D & operator+=(Diff2D const & offset)
    {
        return static_cast<Size2D &>(Diff2D::operator+=(offset));
    }

        /** Decrease size by specified offset.
        */
    Size2D & operator-=(Diff2D const & offset)
    {
        return static_cast<Size2D &>(Diff2D::operator-=(offset));
    }
};

/********************************************************/
/*                                                      */
/*                     Point2D                          */
/*                                                      */
/********************************************************/

/** \brief Two dimensional point or position.

    Specializes \ref Diff2D for the specification of a 2-dimensional
    point or position, in contrast to an extent (for the latter
    use \ref Size2D).

    \code
    // access an image at a point
    value = image[Point2D(10, 20)];
    \endcode

    <b>\#include</b> "<a href="diff2d_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra
*/
class Point2D : public Diff2D
{
public:
        /** The iterator's value type: a coordinate.
        */
    typedef Point2D PixelType;

        /** The iterator's value type: a coordinate.
        */
    typedef Point2D value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef Point2D const & reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef Point2D         index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef Point2D const * pointer;

        /** Default Constructor. Init point at position (0,0)
        */
    Point2D()
    {}

        /** Construct point at given position.
        */
    Point2D(int x, int y)
    : Diff2D(x, y)
    {}

        /** Copy Constructor.
        */
    Point2D(Point2D const & v)
    : Diff2D(v)
    {}

        /** Explicit conversion Constructor.
        */
    explicit Point2D(Diff2D const & v)
    : Diff2D(v)
    {}

        /** Query the points' x coordinate
         */
    int px() const
    {
        return x;
    }

        /** Query the points' y coordinate
         */
    int py() const
    {
        return y;
    }

        /** Copy Assigment.
        */
    Point2D & operator=(Diff2D const & v)
    {
        return static_cast<Point2D &>(Diff2D::operator=(v));
    }

        /** Unary negation.
        */
    Point2D operator-() const
    {
        return Point2D(-x, -y);
    }

        /** Increase point coordinates by specified offset.
        */
    Point2D & operator+=(Diff2D const & offset)
    {
        return static_cast<Point2D &>(Diff2D::operator+=(offset));
    }

        /** Decrease point coordinates by specified offset.
        */
    Point2D & operator-=(Diff2D const & offset)
    {
        return static_cast<Point2D &>(Diff2D::operator-=(offset));
    }

        /** Access current point coordinate.
        */
    reference operator*() const
    {
        return *this;
    }

        /** Read point coordinate at an offset.
        */
    index_reference operator()(int const & dx, int const & dy) const
    {
        return Point2D(x + dx, y + dy);
    }

        /** Read point coordinate at an offset.
        */
    index_reference operator[](Diff2D const & offset) const
    {
        return Point2D(x + offset.x, y + offset.y);
    }

        /** Access current point coordinate.
        */
    pointer operator->() const
    {
        return this;
    }
};

/** Create vector by subtracting specified offset.
 */
inline Diff2D operator-(Diff2D const &a, Diff2D const &b)
{
    return Diff2D(a.x - b.x, a.y - b.y);
}

/** Create size by subtracting specified offset.
 */
inline Size2D operator-(Size2D const & s, Diff2D const &offset)
{
    return Size2D(s.x - offset.x, s.y - offset.y);
}

/** Calculate size of rect between two points.
 */
inline Point2D operator-(Point2D const & s, Diff2D const & offset)
{
    return Point2D(s.x - offset.x, s.y - offset.y);
}

/** The difference of two points is a size
 */
inline Size2D operator-(Point2D const & s, Point2D const & p)
{
    return Size2D(s.x - p.x, s.y - p.y);
}

/** Create vector by adding specified offset.
 */
inline Diff2D operator+(Diff2D const &a, Diff2D const &b)
{
    return Diff2D(a.x + b.x, a.y + b.y);
}

/** Create size by adding specified offset.
 */
inline Size2D operator+(Size2D const &a, Diff2D const &b)
{
    return Size2D(a.x + b.x, a.y + b.y);
}

/** Create point by adding specified offset.
 */
inline Point2D operator+(Point2D const &a, Diff2D const &b)
{
    return Point2D(a.x + b.x, a.y + b.y);
}

/** Add size and point
 */
inline Point2D operator+(Size2D const & s, Point2D const & p)
{
    return Point2D(s.x + p.x, s.y + p.y);
}

inline Point2D operator*(Point2D l, double r)
{
    l *= r;
    return l;
}

inline Point2D operator*(double l, Point2D r)
{
    r *= l;
    return r;
}

inline Size2D operator*(Size2D l, double r)
{
    l *= r;
    return l;
}

inline Size2D operator*(double l, Size2D r)
{
    r *= l;
    return r;
}

inline Point2D operator/(Point2D l, double r)
{
    l /= r;
    return l;
}

inline Size2D operator/(Size2D l, double r)
{
    l /= r;
    return l;
}

inline Point2D operator*(Point2D l, int r)
{
    l *= r;
    return l;
}

inline Point2D operator*(int l, Point2D r)
{
    r *= l;
    return r;
}

inline Size2D operator*(Size2D l, int r)
{
    l *= r;
    return l;
}

inline Size2D operator*(int l, Size2D r)
{
    r *= l;
    return r;
}

inline Point2D operator/(Point2D l, int r)
{
    l /= r;
    return l;
}

inline Size2D operator/(Size2D l, int r)
{
    l /= r;
    return l;
}


/********************************************************/
/*                                                      */
/*                      Rect2D                          */
/*                                                      */
/********************************************************/

/** \brief Two dimensional rectangle.

    This class stores a 2-dimensional rectangular range or region. Thus,
    it follows the VIGRA convention that the upper left corner is inside
    the rectangle, while the lower right is 1 pixel to the right and below the
    last pixel in the rectangle.

    A major advantage of this class is that it can be constructed from either
    a pair of \ref Point2D, or from a \ref Point2D and an extend
    (\ref Size2D). Rect2D overloads operators |=, &=, |, & to realize set
    union (in the sense of a minimal bounding rectangle) and set intersection.

    \code
    Rect2D r1(Point2D(0,0), Point2D(10, 20)),
           r2(Point2D(10, 15), Size2D(20, 20));
    Point2D p(0,100);

    Rect2D r3 =  r1 | r2; // upper left is (0,0), lower right is (30, 35)
    assert(r3.contains(r2));
    assert(!r3.contains(p));

    r3 |= p;       // lower right now (30,101) so that p is inside r3
    assert(r3.contains(p));
    \endcode

    <b>\#include</b> "<a href="diff2d_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra
*/
class Rect2D
{
    Point2D upperLeft_, lowerRight_;

public:
        /** Construct a null rectangle (isEmpty() will return true)
         */
    Rect2D()
    {}

        /** Construct a rectangle representing the given range
         * (lowerRight is considered to be outside the rectangle as
         * usual in the VIGRA)
         */
    Rect2D(Point2D const &upperLeft, Point2D const &lowerRight)
    : upperLeft_(upperLeft), lowerRight_(lowerRight)
    {}

        /** Construct a rectangle representing the given range
         */
    Rect2D(int left, int top, int right, int bottom)
    : upperLeft_(left, top), lowerRight_(right, bottom)
    {}

        /** Construct a rectangle of given position and size
         */
    Rect2D(Point2D const &upperLeft, Size2D const &size)
    : upperLeft_(upperLeft), lowerRight_(upperLeft + size)
    {}

        /** Construct a rectangle of given size at position (0,0)
         */
    explicit Rect2D(Size2D const &size)
    : lowerRight_(Point2D(size))
    {}

        /** Return the first point (scan-order wise) which is
         * considered to be "in" the rectangle.
         */
    Point2D const & upperLeft() const
    {
        return upperLeft_;
    }

        /** Return the first point to the right and below the
         * rectangle.
         */
    Point2D const & lowerRight() const
    {
        return lowerRight_;
    }

        /** Change upperLeft() without changing lowerRight(), which
         * will change the size most probably.
         */
    void setUpperLeft(Point2D const &ul)
    {
        upperLeft_ = ul;
    }

        /** Change lowerRight() without changing upperLeft(), which
         * will change the size most probably.
         */
    void setLowerRight(Point2D const &lr)
    {
        lowerRight_ = lr;
    }

        /** Move the whole rectangle so that the given point will be
         * upperLeft() afterwards.
         */
    void moveTo(Point2D const &newUpperLeft)
    {
        lowerRight_ += newUpperLeft - upperLeft_;
        upperLeft_ = newUpperLeft;
    }

        /** Move the whole rectangle so that upperLeft() will become
         * Point2D(left, top) afterwards.
         */
    void moveTo(int left, int top)
    {
        moveTo(Point2D(left, top));
    }

        /** Move the whole rectangle by the given 2D offset.
         */
    void moveBy(Diff2D const &offset)
    {
        upperLeft_ += offset;
        lowerRight_ += offset;
    }

        /** Move the whole rectangle by the given x- and y-offsets.
         */
    void moveBy(int xOffset, int yOffset)
    {
        moveBy(Diff2D(xOffset, yOffset));
    }

        /** Return the left coordinate of this rectangle.
         */
    int left() const
    {
        return upperLeft_.x;
    }

        /** Return the top coordinate of this rectangle.
         */
    int top() const
    {
        return upperLeft_.y;
    }

        /** Return the right coordinate of this rectangle. That is the
         * first column to the right of the rectangle.
         */
    int right() const
    {
        return lowerRight_.x;
    }

        /** Return the bottom coordinate of this rectangle. That is the
         * first row below the rectangle.
         */
    int bottom() const
    {
        return lowerRight_.y;
    }

        /** Determine and return the width of this rectangle. It might be
         * zero or even negative, and if so, isEmpty() will return true.
         */
    int width() const
    {
        return lowerRight_.x - upperLeft_.x;
    }

        /** Determine and return the height of this rectangle. It might be
         * zero or even negative, and if so, isEmpty() will return true.
         */
    int height() const
    {
        return lowerRight_.y - upperLeft_.y;
    }

        /** Determine and return the area of this rectangle. That is, if
         * this rect isEmpty(), returns zero, otherwise returns
         * width()*height().
         */
    int area() const
    {
        return isEmpty() ? 0 : width()*height();
    }

        /** Determine and return the size of this rectangle. The width
         * and/or height might be zero or even negative, and if so,
         * isEmpty() will return true.
         */
    Size2D size() const
    {
        return lowerRight_ - upperLeft_;
    }

        /** Resize this rectangle to the given extents. This will move
         * the lower right corner only.
         */
    void setSize(Size2D const &size)
    {
        lowerRight_ = upperLeft_ + size;
    }

        /** Resize this rectangle to the given extents. This will move
         * the lower right corner only.
         */
    void setSize(int width, int height)
    {
        lowerRight_ = upperLeft_ + Size2D(width, height);
    }

        /** Increase the size of the rectangle by the given offset. This
         * will move the lower right corner only. (If any of offset's
         * components is negative, the rectangle will get smaller
         * accordingly.)
         */
    void addSize(Size2D const &offset)
    {
        lowerRight_ += offset;
    }

        /** Adds a border of the given width around the rectangle. That
         * means, upperLeft()'s components are moved by -borderWidth
         * and lowerRight()'s by borderWidth. (If borderWidth is
         * negative, the rectangle will get smaller accordingly.)
         */
    void addBorder(int borderWidth)
    {
        upperLeft_ += Diff2D(-borderWidth, -borderWidth);
        lowerRight_ += Diff2D(borderWidth, borderWidth);
    }

        /** Adds a border with possibly different widths in x- and
         * y-directions around the rectangle. That means, each x
         * component is moved borderWidth pixels and each y component
         * is moved borderHeight pixels to the outside. (If
         * borderWidth is negative, the rectangle will get smaller
         * accordingly.)
         */
    void addBorder(int borderWidth, int borderHeight)
    {
        upperLeft_ += Diff2D(-borderWidth, -borderHeight);
        lowerRight_ += Diff2D(borderWidth, borderHeight);
    }

        /// equality check
    bool operator==(Rect2D const &r) const
    {
        return (upperLeft_ == r.upperLeft_) && (lowerRight_ == r.lowerRight_);
    }

        /// inequality check
    bool operator!=(Rect2D const &r) const
    {
        return (upperLeft_ != r.upperLeft_) || (lowerRight_ != r.lowerRight_);
    }

        /** Return whether this rectangle is considered empty. It is
         * non-empty if both coordinates of the lower right corner are
         * greater than the corresponding coordinate of the upper left
         * corner. Uniting an empty rectangle with something will return
         * the bounding rectangle of the 'something', intersecting with an
         * empty rectangle will yield again an empty rectangle.
         */
    bool isEmpty() const
    {
        return ((lowerRight_.x <= upperLeft_.x) ||
                (lowerRight_.y <= upperLeft_.y));
    }

        /** Return whether this rectangle contains the given point. That
         * is, if the point lies within the valid range of an
         * ImageIterator walking from upperLeft() to lowerRight()
         * (excluding the latter).
         */
    bool contains(Point2D const &p) const
    {
        return ((upperLeft_.x <= p.x) &&
                (upperLeft_.y <= p.y) &&
                (p.x < lowerRight_.x) &&
                (p.y < lowerRight_.y));
    }

        /** Return whether this rectangle contains the given
         * one. <tt>r1.contains(r2)</tt> returns the same as
         * <tt>r1 == (r1|r2)</tt> (but is of course more
         * efficient). That also means, a rectangle (even an empty one!)
         * contains() any empty rectangle.
         */
    bool contains(Rect2D const &r) const
    {
        return r.isEmpty() ||
            contains(r.upperLeft()) && contains(r.lowerRight()-Diff2D(1,1));
    }

        /** Return whether this rectangle overlaps with the given
         * one. <tt>r1.intersects(r2)</tt> returns the same as
         * <tt>!(r1&r2).isEmpty()</tt> (but is of course much more
         * efficient).
         */
    bool intersects(Rect2D const &r) const
    {
        return ((r.upperLeft_.x < lowerRight_.x) &&
                (upperLeft_.x < r.lowerRight_.x) &&
                (r.upperLeft_.y < lowerRight_.y) &&
                (upperLeft_.y < r.lowerRight_.y))
            && !r.isEmpty();
    }

        /** Modifies this rectangle by including the given point. The
         * result is the bounding rectangle of the rectangle and the
         * point. If isEmpty returns true, the union will be a
         * rectangle containing only the given point.
         */
    Rect2D &operator|=(Point2D const &p)
    {
        if(isEmpty())
        {
            upperLeft_ = p;
            lowerRight_ = p + Diff2D(1, 1);
        }
        else
        {
            if(p.x < upperLeft_.x)
                upperLeft_.x = p.x;
            if(p.y < upperLeft_.y)
                upperLeft_.y = p.y;
            if(lowerRight_.x <= p.x)
                lowerRight_.x = p.x + 1;
            if(lowerRight_.y <= p.y)
                lowerRight_.y = p.y + 1;
        }
        return *this;
    }

        /** Returns the union of this rectangle and the given
         * point. The result is the bounding rectangle of the
         * rectangle and the point. If isEmpty returns true, the union
         * will be a rectangle containing only the given point.
         */
    Rect2D operator|(Point2D const &p) const
    {
        Rect2D result(*this);
        result |= p;
        return result;
    }

        /** Modifies this rectangle by uniting it with the given
         * one. The result is the bounding rectangle of both
         * rectangles. If one of the rectangles isEmpty(), the union
         * will be the other one.
         */
    Rect2D &operator|=(Rect2D const &r)
    {
        if(r.isEmpty())
            return *this;
        if(isEmpty())
            return operator=(r);

        if(r.upperLeft_.x < upperLeft_.x)
            upperLeft_.x = r.upperLeft_.x;
        if(r.upperLeft_.y < upperLeft_.y)
            upperLeft_.y = r.upperLeft_.y;
        if(lowerRight_.x < r.lowerRight_.x)
            lowerRight_.x = r.lowerRight_.x;
        if(lowerRight_.y < r.lowerRight_.y)
            lowerRight_.y = r.lowerRight_.y;
        return *this;
    }

        /** Returns the union of this rectangle and the given one. The
         * result is the bounding rectangle of both rectangles. If one
         * of the rectangles isEmpty(), the union will be the other
         * one.
         */
    Rect2D operator|(Rect2D const &r) const
    {
        Rect2D result(*this);
        result |= r;
        return result;
    }

        /** Modifies this rectangle by intersecting it with the given
         * point. The result is the bounding rect of the point (with
         * width and height equal to 1) if it was contained in the
         * original rect, or an empty rect otherwise.
         */
    Rect2D &operator&=(Point2D const &p)
    {
        if(contains(p))
        {
            upperLeft_ = p;
            lowerRight_ = p + Diff2D(1, 1);
        }
        else
            lowerRight_ = upperLeft_;
        return *this;
    }

        /** Intersects this rectangle with the given point. The result
         * is the bounding rect of the point (with width and height
         * equal to 1) if it was contained in the original rect, or an
         * empty rect otherwise.
         */
    Rect2D operator&(Point2D const &p) const
    {
        Rect2D result(*this);
        result &= p;
        return result;
    }

        /** Modifies this rectangle by intersecting it with the given
         * one. The result is the maximal rectangle contained in both
         * original ones. Intersecting with an empty rectangle will
         * yield again an empty rectangle.
         */
    Rect2D &operator&=(Rect2D const &r)
    {
        if(isEmpty())
            return *this;
        if(r.isEmpty())
            return operator=(r);

        if(upperLeft_.x < r.upperLeft_.x)
            upperLeft_.x = r.upperLeft_.x;
        if(upperLeft_.y < r.upperLeft_.y)
            upperLeft_.y = r.upperLeft_.y;
        if(r.lowerRight_.x < lowerRight_.x)
            lowerRight_.x = r.lowerRight_.x;
        if(r.lowerRight_.y < lowerRight_.y)
            lowerRight_.y = r.lowerRight_.y;
        return *this;
    }

        /** Intersects this rectangle with the given one. The result
         * is the maximal rectangle contained in both original ones.
         * Intersecting with an empty rectangle will yield again an
         * empty rectangle.
         */
    Rect2D operator&(Rect2D const &r) const
    {
        Rect2D result(*this);
        result &= r;
        return result;
    }
};

/********************************************************/
/*                                                      */
/*                      Dist2D                          */
/*                                                      */
/********************************************************/

/** @deprecated use \ref vigra::Diff2D instead
*/
class Dist2D
{
  public:
    Dist2D(int the_width, int the_height)
    : width(the_width),
      height(the_height)
    {}

    Dist2D(Dist2D const & s)
    : width(s.width),
      height(s.height)
    {}

    Dist2D & operator=(Dist2D const & s)
    {
        if(this != &s)
        {
            width = s.width;
            height = s.height;
        }
        return *this;
    }

    Dist2D & operator+=(Dist2D const & s)
    {
        width += s.width;
        height += s.height;

        return *this;
    }

    Dist2D  operator+(Dist2D const & s) const
    {
        Dist2D ret(*this);
        ret += s;

        return ret;
    }

    operator Diff2D()
        { return Diff2D(width, height); }

    int width;
    int height;
 };

//@}

inline
std::ostream & operator<<(std::ostream & o, vigra::Diff2D const & d)
{
    o << '(' << d.x << ", " << d.y << ')';
    return o;
}

inline
std::ostream &operator <<(std::ostream &s, vigra::Size2D const &d)
{
    s << '(' << d.x << 'x' << d.y << ')';
    return s;
}

inline
std::ostream &operator <<(std::ostream &s, vigra::Rect2D const &r)
{
    s << "[" << r.upperLeft() << " to " << r.lowerRight()
      << " = " << r.size() << "]";
    return s;
}

} // namespace vigra

#endif // VIGRA_DIFF2D_HXX
