/************************************************************************/
/*                                                                      */
/*          Copyright 1998-2002 by Hans Meine, Ullrich Koethe           */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_PIXELNEIGHBORHOOD_HXX
#define VIGRA_PIXELNEIGHBORHOOD_HXX

#include <vigra/utilities.hxx>

namespace vigra {

/** \addtogroup PixelNeighborhood Utilities to manage pixel neighborhoods

    4- and 8-neighborhood definitions and circulators.

    <b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>

    <b>See also:</b> \ref vigra::NeighborhoodCirculator
 */
//@{

/********************************************************/
/*                                                      */
/*                    FourNeighborhood                  */
/*                                                      */
/********************************************************/

/** Utilities for 4-neighborhood. */
namespace FourNeighborhood
{

/** \brief Encapsulation of direction management for 4-neighborhood.

    This helper class allows the transformation between Freeman chain codes
    (East = 0, North = 1 etc.) and the corresponding Diff2D instances
    and back.

    You can either use the chain codes by explicit qualification:

    \code
    // the following three lines are equivalent
    FourNeighborhood::NeighborCode::Direction d = FourNeighborhood::NeighborCode::East;
    FourNeighborCode::Direction d = FourNeighborCode::East;
    FourNeighborhood::Direction d = FourNeighborhood::East;
    \endcode

    or you can fix 4-neighborhood by importing the entire namespace in
    your function:

    \code
    using namespace FourNeighborhood;

    Direction d = East;
    \endcode

    If you want to pass 4-neighborhood codes as a template parameter, use
    the class FourNeighborhood::NeighborCode.

    <b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>
    Namespace: vigra::FourNeighborhood
*/
class NeighborCode
{
  public:
        /** Freeman direction codes for 4-neighborhood.
            East = 0, North = 1 etc.
            DirectionCount may be used for portable loop termination conditions.
        */
    enum Direction {
        Error = -1,    ///< &nbsp;
        East = 0,      ///< &nbsp;
        North,         ///< &nbsp;
        West,          ///< &nbsp;
        South,         ///< &nbsp;
        DirectionCount ///< &nbsp;
    };

        /** Transform direction code into corresponding Diff2D offset.
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & diff(Direction code)
    {
        static Diff2D d[] = {
            Diff2D(1, 0), Diff2D(0, -1), Diff2D(-1, 0), Diff2D(0, 1)
        };
        return d[code];
    }

        /** Equivalent to <tt>diff(static_cast<Direction>(code))</tt>.
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & diff(int code) { return diff(static_cast<Direction>(code)); }

        /** Get the relative offset from one neighbor to the other.
            For example, <tt>relativeDiff(East, West) == Diff2D(-2,0)</tt>.
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & relativeDiff(Direction fromCode, Direction toCode)
    {
        static Diff2D d[][4] = {
            { Diff2D(0, 0), Diff2D(-1, -1), Diff2D(-2, 0), Diff2D(-1, 1) },
            { Diff2D(1, 1), Diff2D(0, 0), Diff2D(-1, 1), Diff2D(0, 2) },
            { Diff2D(2, 0), Diff2D(1, -1), Diff2D(0, 0), Diff2D(1, 1) },
            { Diff2D(1, -1), Diff2D(0, -2), Diff2D(-1, -1), Diff2D(0, 0) }
        };

        return d[fromCode][toCode];
    }

        /** Equivalent to relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode)).
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & relativeDiff(int fromCode, int toCode)
    {
        return relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode));
    }

        /**  X-component of diff() */
    static int dX(Direction code) { return diff(code).x; }
        /**  Y-component of diff() */
    static int dY(Direction code) { return diff(code).y; }
        /**  X-component of diff() */
    static int dX(int code) { return diff(code).x; }
        /**  Y-component of diff() */
    static int dY(int code) { return diff(code).y; }

        /** Transform Diff2D offset into corresponding direction code.
            The code <tt>Direction::Error</tt> will be returned if <tt>diff</tt>
            is not in the 4-neighborhood.
        */
    static Direction code(Diff2D const & diff)
    {
        switch(diff.x)
        {
            case  0:
            {
                switch(diff.y)
                {
                    case 1:
                        return South;
                    case -1:
                        return North;
                    default:
                        return Error;
                }
            }
            case -1:
            {
                return (diff.y == 0) ?
                            West :
                            Error;
            }
            case  1:
            {
                return (diff.y == 0) ?
                            East :
                            Error;
            }
        }
        return Error;
    }

        /** Check whether a code refers to a diagonal direction.
            Useful if you want to abstract the differences between 4- and 8-neighborhood.
            Always <tt>false</tt> for 4-neighborhood.
        */
    static bool isDiagonal(Direction) { return false; }

    static Diff2D const & right()        { return diff(East); }    /**<  Offset to the right neighbor */
    static Diff2D const & top()          { return diff(North); }   /**<  Offset to the top neighbor */
    static Diff2D const & left()         { return diff(West); }    /**<  Offset to the left neighbor */
    static Diff2D const & bottom()       { return diff(South); }   /**<  Offset to the bottom neighbor */

    static Diff2D const & east()       { return diff(East); }    /**<  Offset to the east neighbor */
    static Diff2D const & north()      { return diff(North); }   /**<  Offset to the north neighbor */
    static Diff2D const & west()       { return diff(West); }    /**<  Offset to the west neighbor */
    static Diff2D const & south()      { return diff(South); }   /**<  Offset to the south neighbor */
};

    /** Export NeighborCode::Direction into the scope of namespace FourNeighborhood.
    */
typedef NeighborCode::Direction Direction;

static const Direction East           = NeighborCode::East;           /**<  Export NeighborCode::East to namespace FourNeighborhood */
static const Direction North          = NeighborCode::North;          /**<  Export NeighborCode::North to namespace FourNeighborhood */
static const Direction West           = NeighborCode::West;           /**<  Export NeighborCode::West to namespace FourNeighborhood */
static const Direction South          = NeighborCode::South;          /**<  Export NeighborCode::South to namespace FourNeighborhood */
static const Direction DirectionCount = NeighborCode::DirectionCount; /**<  Export NeighborCode::DirectionCount to namespace FourNeighborhood */

} // namespace FourNeighborhood

    /** Export \ref vigra::FourNeighborhood::NeighborCode into the scope of namespace vigra.
    */
typedef FourNeighborhood::NeighborCode FourNeighborCode;

/********************************************************/
/*                                                      */
/*                   EightNeighborhood                  */
/*                                                      */
/********************************************************/

/** Utilities for 8-neighborhood. */
namespace EightNeighborhood
{
/** \brief Encapsulation of direction management for the 8-neighborhood.

    This helper class allows the transformation between Freeman chain codes
    (East = 0, NorthEast = 1 etc.) and the corresponding Diff2D instances
    and back.

    You can either use the chain codes by explicit qualification:

    \code
    // the following three lines are equivalent
    EightNeighborhood::NeighborCode::Direction d = EightNeighborhood::NeighborCode::East;
    EightNeighborCode::Direction d               = EightNeighborCode::East;
    EightNeighborhood::Direction d               = EightNeighborhood::East;
    \endcode

    or you can fix 8-neighborhood by importing the entire namespace in
    your function:

    \code
    using namespace EightNeighborhood;

    Direction d = East;
    \endcode

    If you want to pass 8-neighborhood codes as a template parameter, use
    the class EightNeighborhood::NeighborCode.

    <b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>
    Namespace: vigra::EightNeighborhood
*/
class NeighborCode
{
  public:
        /** Freeman direction codes for the 8-neighborhood.
            East = 0, NorthEast = 1 etc.
            DirectionCount may be used for portable loop termination conditions.
        */
    enum Direction {
        Error = -1,     ///< &nbsp;
        East = 0,       ///< &nbsp;
        NorthEast,      ///< &nbsp;
        North,          ///< &nbsp;
        NorthWest,      ///< &nbsp;
        West,           ///< &nbsp;
        SouthWest,      ///< &nbsp;
        South,          ///< &nbsp;
        SouthEast,      ///< &nbsp;
        DirectionCount  ///< &nbsp;
    };

        /** Transform direction code into corresponding Diff2D offset.
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & diff(Direction code)
    {
        static Diff2D d[] = {
            Diff2D(1, 0), Diff2D(1, -1), Diff2D(0, -1), Diff2D(-1, -1),
            Diff2D(-1, 0), Diff2D(-1, 1), Diff2D(0, 1), Diff2D(1, 1)
        };
        return d[code];
    }

        /** Equivalent to diff(static_cast<Direction>(code)).
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & diff(int code) { return diff(static_cast<Direction>(code)); }

        /** Get the relative offset from one neighbor to the other.
            For example, <tt>relativeDiff(East, West) == Diff2D(-2,0)</tt>.
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & relativeDiff(Direction fromCode, Direction toCode)
    {
        static Diff2D d[][8] = {
            { Diff2D(0, 0), Diff2D(0, -1), Diff2D(-1, -1), Diff2D(-2, -1),
              Diff2D(-2, 0), Diff2D(-2, 1), Diff2D(-1, 1), Diff2D(0, 1) },
            { Diff2D(0, 1), Diff2D(0, 0), Diff2D(-1, 0), Diff2D(-2, 0),
              Diff2D(-2, 1), Diff2D(-2, 2), Diff2D(-1, 2), Diff2D(0, 2) },
            { Diff2D(1, 1), Diff2D(1, 0), Diff2D(0, 0), Diff2D(-1, 0),
              Diff2D(-1, 1), Diff2D(-1, 2), Diff2D(0, 2), Diff2D(1, 2) },
            { Diff2D(2, 1), Diff2D(2, 0), Diff2D(1, 0), Diff2D(0, 0),
              Diff2D(0, 1), Diff2D(0, 2), Diff2D(1, 2), Diff2D(2, 2) },
            { Diff2D(2, 0), Diff2D(2, -1), Diff2D(1, -1), Diff2D(0, -1),
              Diff2D(0, 0), Diff2D(0, 1), Diff2D(1, 1), Diff2D(2, 1) },
            { Diff2D(2, -1), Diff2D(2, -2), Diff2D(1, -2), Diff2D(0, -2),
              Diff2D(0, -1), Diff2D(0, 0), Diff2D(1, 0), Diff2D(2, 0) },
            { Diff2D(1, -1), Diff2D(1, -2), Diff2D(0, -2), Diff2D(-1, -2),
              Diff2D(-1, -1), Diff2D(-1, 0), Diff2D(0, 0), Diff2D(1, 0) },
            { Diff2D(0, -1), Diff2D(0, -2), Diff2D(-1, -2), Diff2D(-2, -2),
              Diff2D(-2, -1), Diff2D(-2, 0), Diff2D(-1, 0), Diff2D(0, 0) }
        };

        return d[fromCode][toCode];
    }

        /** Equivalent to relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode)).
            (note: there is no bounds checking on the code you pass.)
        */
    static Diff2D const & relativeDiff(int fromCode, int toCode)
    {
        return relativeDiff(static_cast<Direction>(fromCode), static_cast<Direction>(toCode));
    }

        /**  X-component of diff() */
    static int dX(Direction code) { return diff(code).x; }
        /**  Y-component of diff() */
    static int dY(Direction code) { return diff(code).y; }
        /**  X-component of diff() */
    static int dX(int code) { return diff(code).x; }
        /**  Y-component of diff() */
    static int dY(int code) { return diff(code).y; }

        /** Transform 4-neighborhood code into 8-neighborhood code.
        */
    static Direction code(FourNeighborhood::Direction d)
        { return static_cast<Direction>(2*d); }

        /** Transform Diff2D offset into corresponding direction code.
            The code <tt>Direction::Error</tt> will be returned if <tt>diff</tt>
            is not in the 8-neighborhood.
        */
    static Direction code(Diff2D const & diff)
    {
        switch(diff.x)
        {
            case  0:
            {
                switch(diff.y)
                {
                    case 1:
                        return South;
                    case -1:
                        return North;
                    default:
                        return Error;
                }
            }
            case -1:
            {
                switch(diff.y)
                {
                    case 0:
                        return West;
                    case 1:
                        return SouthWest;
                    case -1:
                        return NorthWest;
                    default:
                        return Error;
                }
            }
            case  1:
            {
                switch(diff.y)
                {
                    case 0:
                        return East;
                    case 1:
                        return SouthEast;
                    case -1:
                        return NorthEast;
                    default:
                        return Error;
                }
            }
        }
        return Error;
    }

        /** Check whether a code refers to a diagonal direction.
            Useful if you want to abstract the differences between 4- and 8-neighborhood.
        */
    static bool isDiagonal(Direction code) { return (code % 2) != 0; }

    static Diff2D const & right()        { return diff(East); }        /**<  Offset to the right neighbor */
    static Diff2D const & topRight()     { return diff(NorthEast); }   /**<  Offset to the topRight neighbor */
    static Diff2D const & top()          { return diff(North); }       /**<  Offset to the top neighbor */
    static Diff2D const & topLeft()      { return diff(NorthWest); }   /**<  Offset to the topLeft neighbor */
    static Diff2D const & left()         { return diff(West); }        /**<  Offset to the left neighbor */
    static Diff2D const & bottomLeft()   { return diff(SouthWest); }   /**<  Offset to the bottomLeft neighbor */
    static Diff2D const & bottom()       { return diff(South); }       /**<  Offset to the bottom neighbor */
    static Diff2D const & bottomRight()  { return diff(SouthEast); }   /**<  Offset to the bottomRight neighbor */

    static Diff2D const & east()       { return diff(East); }        /**<  Offset to the east neighbor */
    static Diff2D const & northEast()  { return diff(NorthEast); }   /**<  Offset to the northEast neighbor */
    static Diff2D const & north()      { return diff(North); }       /**<  Offset to the north neighbor */
    static Diff2D const & northWest()  { return diff(NorthWest); }   /**<  Offset to the northWest neighbor */
    static Diff2D const & west()       { return diff(West); }        /**<  Offset to the west neighbor */
    static Diff2D const & southWest()  { return diff(SouthWest); }   /**<  Offset to the southWest neighbor */
    static Diff2D const & south()      { return diff(South); }       /**<  Offset to the south neighbor */
    static Diff2D const & southEast()  { return diff(SouthEast); }   /**<  Offset to the southEast neighbor */
};

    /** Export NeighborCode::Direction into the scope of namespace EightNeighborhood.
    */
typedef NeighborCode::Direction Direction;

static const Direction East           = NeighborCode::East;        /**<  Export NeighborCode::East to namespace EightNeighborhood */
static const Direction NorthEast      = NeighborCode::NorthEast;   /**<  Export NeighborCode::NorthEast to namespace EightNeighborhood */
static const Direction North          = NeighborCode::North;       /**<  Export NeighborCode::North to namespace EightNeighborhood */
static const Direction NorthWest      = NeighborCode::NorthWest;   /**<  Export NeighborCode::NorthWest to namespace EightNeighborhood */
static const Direction West           = NeighborCode::West;        /**<  Export NeighborCode::West to namespace EightNeighborhood */
static const Direction SouthWest      = NeighborCode::SouthWest;   /**<  Export NeighborCode::SouthWest to namespace EightNeighborhood */
static const Direction South          = NeighborCode::South;       /**<  Export NeighborCode::South to namespace EightNeighborhood */
static const Direction SouthEast      = NeighborCode::SouthEast;   /**<  Export NeighborCode::SouthEast to namespace EightNeighborhood */
static const Direction DirectionCount = NeighborCode::DirectionCount;   /**<  Export NeighborCode::DirectionCount to namespace EightNeighborhood */

} // namespace EightNeighborhood

    /** Export \ref vigra::EightNeighborhood::NeighborCode into the scope of namespace vigra.
    */
typedef EightNeighborhood::NeighborCode EightNeighborCode;

/********************************************************/
/*                                                      */
/*              NeighborOffsetCirculator                */
/*                                                      */
/********************************************************/

/** \brief Circulator that walks around a given location.

    The template parameter defines the kind of neighborhood used, e.g.

    \code
    NeighborOffsetCirculator<EightNeighborCode> eight_circulator;
    NeighborOffsetCirculator<FourNeighborCode>  four_circulator;
    \endcode

    Since this circulator doesn't now about the pixels in any particular image,
    you usually doesn't use it directly but rather as a base class or helper for
    neighborhood circulators refering to a particular image (e.g. NeighborhoodCirculator)

    <b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>
    Namespace: vigra
*/
template<class NEIGHBORCODE>
class NeighborOffsetCirculator
: public NEIGHBORCODE
{
public:
    typedef NEIGHBORCODE NeighborCode;

        /** return type of direction()
        */
    typedef typename NEIGHBORCODE::Direction Direction;

        /** the circulator's value type
        */
    typedef Diff2D value_type;

        /** the circulator's reference type (return type of <TT>*circ</TT>)
        */
    typedef Diff2D const & reference;

        /** the circulator's index reference type (return type of <TT>circ[n]</TT>)
        */
    typedef Diff2D const & index_reference;

        /** the circulator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef Diff2D const * pointer;

        /** the circulator's difference type (argument type of <TT>circ[diff]</TT>)
        */
    typedef int difference_type;

        /** the circulator tag (random access iterator)
        */
    typedef random_access_circulator_tag iterator_category;

protected:
    Direction direction_;

public:
        /** Create circulator refering to the given direction.
        */
    NeighborOffsetCirculator(Direction dir = East)
        : direction_(dir)
    {
    }

        /** pre-increment */
    NeighborOffsetCirculator & operator++()
    {
        direction_ = static_cast<Direction>((direction_+1) % DirectionCount);
        return *this;
    }

        /** pre-decrement */
    NeighborOffsetCirculator & operator--()
    {
        direction_ = static_cast<Direction>((direction_ + DirectionCount-1) % DirectionCount);
        return *this;
    }

        /** post-increment */
    NeighborOffsetCirculator operator++(int)
    {
        NeighborOffsetCirculator ret(*this);
        operator++();
        return ret;
    }

        /** post-decrement */
    NeighborOffsetCirculator operator--(int)
    {
        NeighborOffsetCirculator ret(*this);
        operator--();
        return ret;
    }

        /** add-assignment */
    NeighborOffsetCirculator & operator+=(difference_type d)
    {
        direction_ = static_cast<Direction>((direction_ + d) % DirectionCount);
        if(direction_ < 0)
            direction_ = static_cast<Direction>(direction_ + DirectionCount);
        return *this;
    }

        /** subtract-assignment */
    NeighborOffsetCirculator & operator-=(difference_type d)
    {
        direction_ = static_cast<Direction>((direction_ - d) % DirectionCount);
        if(direction_ < 0)
            direction_ = static_cast<Direction>(direction_ + DirectionCount);
        return *this;
    }

        /** addition */
    NeighborOffsetCirculator operator+(difference_type d) const
    {
        return NeighborOffsetCirculator(*this) += d;
    }

        /** subtraction */
    NeighborOffsetCirculator operator-(difference_type d) const
    {
        return NeighborOffsetCirculator(*this) -= d;
    }

        /** Move to the direction that is 'right' relative to the current direction.
            This is equivalent to <tt>four_circulator--</tt> and
            <tt>eight_circulator -= 2</tt> respectively.
        */
    NeighborOffsetCirculator & turnRight()
    {
        direction_ = static_cast<Direction>((direction_ + South) % DirectionCount);
        return *this;
    }

        /** Move to the direction that is 'left' relative to the current direction.
            This is equivalent to <tt>four_circulator++</tt> and
            <tt>eight_circulator += 2</tt> respectively.
        */
    NeighborOffsetCirculator & turnLeft()
    {
        direction_ = static_cast<Direction>((direction_ + North) % DirectionCount);
        return *this;
    }

        /** Move to the opposite direction of the current direction.
            This is equivalent to <tt>four_circulator += 2</tt> and
            <tt>eight_circulator += 4</tt> respectively.
        */
    NeighborOffsetCirculator & turnRound()
    {
        direction_ = static_cast<Direction>((direction_ + West) % DirectionCount);
        return *this;
    }

        /** Move to the given direction.
        */
    NeighborOffsetCirculator & turnTo(Direction d)
    {
        direction_ = d;
        return *this;
    }

        /** equality */
    bool operator==(NeighborOffsetCirculator const & o) const
    {
        return direction_ == o.direction_;
    }

        /** unequality */
    bool operator!=(NeighborOffsetCirculator const & o) const
    {
        return direction_ != o.direction_;
    }

        /** subtraction */
    difference_type operator-(NeighborOffsetCirculator const & o) const
    {
        return direction_ - o.direction_;
    }

        /** dereference */
    reference operator*() const
    {
        return diff();
    }

        /** index */
    index_reference operator[](difference_type d) const
    {
        return NEIGHBORCODE::diff(direction(d));
    }

        /** member access */
    pointer operator->() const
    {
        return &diff();
    }

        /** Get Diff2D offset from center to current neighbor.
        */
    Diff2D const & diff() const
    {
        return NEIGHBORCODE::diff(direction_);
    }

        /** Get Diff2D offset to given direction.
        */
    static Diff2D const & diff(Direction dir)
    {
        return NEIGHBORCODE::diff(dir);
    }

        /** Get relative distance (Diff2D) from current neighbor to neighbor
            at given offset.
        */
    Diff2D const &relativeDiff(difference_type offset) const
    {
        Direction toDir = static_cast<Direction>((direction_ + offset) % DirectionCount);
        if(toDir < 0)
            toDir = static_cast<Direction>(toDir + DirectionCount);
        return NEIGHBORCODE::relativeDiff(direction_, toDir);
    }

        /** X-component of diff()  */
    int dX() const
    {
        return NEIGHBORCODE::dX(direction_);
    }

        /** Y-component of diff() */
    int dY() const
    {
        return NEIGHBORCODE::dY(direction_);
    }

        /** Check whether current direction is a diagonal one.
        */
    bool isDiagonal() const
    {
        return NEIGHBORCODE::isDiagonal(direction_);
    }

        /** Get current direction.
        */
    Direction direction() const
    {
        return direction_;
    }

        /** Get direction code at offset of current direction.
        */
    Direction direction(difference_type offset) const
    {
        int result = (direction_ + offset) % DirectionCount;
        if(result < 0)
            result += DirectionCount;
        return static_cast<Direction>(result);
    }
};

/** Specialization of NeighborOffsetCirculator for 8-neighborhood.
*/
typedef NeighborOffsetCirculator<EightNeighborCode> EightNeighborOffsetCirculator;

/** Specialization of NeighborOffsetCirculator for 4-neighborhood.
*/
typedef NeighborOffsetCirculator<FourNeighborCode> FourNeighborOffsetCirculator;


//@}

/** \addtogroup ImageIteratorAdapters
 */
//@{

/********************************************************/
/*                                                      */
/*                NeighborhoodCirculator                */
/*                                                      */
/********************************************************/

/** \brief Circulator that walks around a given location in a given image.

    The template parameters define the kind of neighborhood used and the underlying
    image, e.g.

    \code
    NeighborhoodCirculator<BImage::traverser, EightNeighborCode> eight_circulator(image.upperLeft()+Diff2D(2,2));
    NeighborhoodCirculator<BImage::traverser, FourNeighborCode>  four_circulator(image.upperLeft()+Diff2D(2,2));
    \endcode

    The access functions return the value of the current neighbor pixel. Use <tt>center()</tt> to
    access the center pixel of the neighborhood. The center can be changed by calling
    <tt>moveCenterToNeighbor()</tt> or <tt>swapCenterNeighbor()</tt>.

    <b>\#include</b> "<a href="pixelneighborhood_8hxx-source.html">vigra/pixelneighborhood.hxx</a>"<br>
    Namespace: vigra
*/
template <class IMAGEITERATOR, class NEIGHBORCODE>
class NeighborhoodCirculator : private IMAGEITERATOR
{
    typedef NeighborOffsetCirculator<NEIGHBORCODE> NEIGHBOROFFSETCIRCULATOR;

public:
        /** type of the underlying image iterator
        */
    typedef IMAGEITERATOR base_type;

        /** type of the used neighbor code
        */
    typedef NEIGHBORCODE NeighborCode;

        /** the circulator's value type
        */
    typedef typename IMAGEITERATOR::value_type value_type;

        /** type of the direction code
        */
    typedef typename NEIGHBORCODE::Direction Direction;

        /** the circulator's reference type (return type of <TT>*circ</TT>)
        */
    typedef typename IMAGEITERATOR::reference reference;

        /** the circulator's index reference type (return type of <TT>circ[n]</TT>)
        */
    typedef typename IMAGEITERATOR::index_reference index_reference;

        /** the circulator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename IMAGEITERATOR::pointer pointer;

        /** the circulator's difference type (argument type of <TT>circ[diff]</TT>)
        */
    typedef typename NEIGHBOROFFSETCIRCULATOR::difference_type difference_type;

        /** the circulator tag (random_access_circulator_tag)
        */
    typedef typename NEIGHBOROFFSETCIRCULATOR::iterator_category iterator_category;

        /** Construct circulator with given <tt>center</tt> pixel, pointing to the neighbor
            at the given direction <tt>d</tt>.
        */
    NeighborhoodCirculator(IMAGEITERATOR const & center = IMAGEITERATOR(),
                           Direction d = NEIGHBOROFFSETCIRCULATOR::East)
        : IMAGEITERATOR(center), neighborCode_(d)
    {
        IMAGEITERATOR::operator+=(neighborCode_.diff());
    }

        /** pre-increment */
    NeighborhoodCirculator & operator++()
    {
        return operator+=(1);
    }

        /** pre-decrement */
    NeighborhoodCirculator operator++(int)
    {
        NeighborhoodCirculator ret(*this);
        operator++();
        return ret;
    }

        /** post-increment */
    NeighborhoodCirculator & operator--()
    {
        return operator+=(-1);
    }

        /** post-decrement */
    NeighborhoodCirculator operator--(int)
    {
        NeighborhoodCirculator ret(*this);
        operator--();
        return ret;
    }

        /** add-assignment */
    NeighborhoodCirculator & operator+=(difference_type d)
    {
        IMAGEITERATOR::operator+=(neighborCode_.relativeDiff(d));
        neighborCode_+= d;
        return *this;
    }

        /** subtract-assignment */
    NeighborhoodCirculator & operator-=(difference_type d)
    {
        return operator+=(-d);
    }

        /** addition */
    NeighborhoodCirculator operator+(difference_type d) const
    {
        NeighborhoodCirculator result(*this);
        result+= d;
        return result;
    }

        /** subtraction */
    NeighborhoodCirculator operator-(difference_type d) const
    {
        NeighborhoodCirculator result(*this);
        result-= d;
        return result;
    }

        /** Move to the direction that is 'right' relative to the current direction.
            This is equivalent to <tt>four_circulator--</tt> and
            <tt>eight_circulator -= 2</tt> respectively.
        */
    NeighborhoodCirculator & turnRight()
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnRight();
        IMAGEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

        /** Move to the direction that is 'left' relative to the current direction.
            This is equivalent to <tt>four_circulator++</tt> and
            <tt>eight_circulator += 2</tt> respectively.
        */
    NeighborhoodCirculator & turnLeft()
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnLeft();
        IMAGEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

        /** Move to the opposite direction of the current direction.
            This is equivalent to <tt>four_circulator += 2</tt> and
            <tt>eight_circulator += 4</tt> respectively.
        */
    NeighborhoodCirculator & turnRound()
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnRound();
        IMAGEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

        /** Move to the given direction.
        */
    NeighborhoodCirculator & turnTo(Direction d)
    {
        Direction oldDirection = neighborCode_.direction();
        neighborCode_.turnTo(d);
        IMAGEITERATOR::operator+=(NeighborCode::relativeDiff
                                  (oldDirection, neighborCode_.direction()));
        return *this;
    }

        /** Move the center in the current direction.
            The current neighbor becomes the new center, the direction does not change.
        */
    NeighborhoodCirculator & moveCenterToNeighbor()
    {
        IMAGEITERATOR::operator+=(neighborCode_.diff());
        return *this;
    }

        /** Exchange the center with the current neighbor.
            Equivalent to <tt>circ.moveCenterToNeighbor().turnRound()</tt>
            (but shorter and more efficient).
        */
    NeighborhoodCirculator & swapCenterNeighbor()
    {
        neighborCode_.turnRound();
        IMAGEITERATOR::operator+=(neighborCode_.diff());
        return *this;
    }

        /** equality */
    bool operator==(NeighborhoodCirculator const & rhs) const
    {
        return neighborCode_ == rhs.neighborCode_ &&
               IMAGEITERATOR::operator==(rhs);
    }

        /** inequality */
    bool operator!=(NeighborhoodCirculator const & rhs) const
    {
        return neighborCode_ != rhs.neighborCode_ ||
               IMAGEITERATOR::operator!=(rhs);
    }

        /** subtraction */
    difference_type operator-(NeighborhoodCirculator const & rhs) const
    {
        return neighborCode_ - rhs.neighborCode_;
    }

        /** dereference */
    reference operator*() const
    {
        return IMAGEITERATOR::operator*();
    }

        /** index */
    index_reference operator[](difference_type d) const
    {
        return IMAGEITERATOR::operator[](neighborCode_.relativeDiff(d));
    }

        /** member access */
    pointer operator->() const
    {
        return IMAGEITERATOR::operator->();
    }

        /** Get the base iterator for the current neighbor. */
    base_type const & base() const
    {
        return *this;
    }

        /** Get the base iterator for the center of the circulator. */
    base_type center() const
    {
        return (base_type)*this - neighborCode_.diff();
    }

        /** Get the current direction. */
    Direction direction() const
    {
        return neighborCode_.direction();
    }

        /** Get the difference vector (Diff2D) from the center to the current neighbor. */
    Diff2D const & diff() const
    {
        return neighborCode_.diff();
    }

        /** Is the current neighbor a diagonal neighbor? */
    bool isDiagonal() const
    {
        return neighborCode_.isDiagonal();
    }

private:
    NEIGHBOROFFSETCIRCULATOR neighborCode_;
};

//@}

} // namespace vigra

#endif /* VIGRA_PIXELNEIGHBORHOOD_HXX */
