/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_BASICS_HXX
#define VIGRA_BASICS_HXX

#include <cmath>     // for sqrt()
#include <utility>    // for pair

#include "vigra/config.hxx"
#include "vigra/error.hxx"

namespace vigra {

struct VigraTrueType 
{
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { asBool = true };
#else
    static const bool asBool = true;
#endif
};
struct VigraFalseType
{
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { asBool = false };
#else
    static const bool asBool = false;
#endif
};

/*! \page Utilities Utilities
    Basic helper functionality needed throughout.
     
    <DL>
    <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
     \ref vigra::Diff2D 
     <DD><em>Two dimensional difference vector</em>
     <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
     <b>vigra::Dist2D</b>
     <DD><em>Deprecated - use \ref vigra::Diff2D</em>
     <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
     \ref TupleTypes
     <DD><em>pair, triple, tuple4, tuple5</em>
      <DT>
    <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif"> 
     \ref MathConstants
     <DD><em>M_PI, M_SQRT2</em>
    </DL>
*/

/********************************************************/
/*                                                      */
/*                      Diff2D                          */
/*                                                      */
/********************************************************/

/*! \brief Two dimensional difference vector.
    
    This class acts primarily as a difference vector for specifying 
    pixel coordinates and region sizes. In addition, Diff2D fulfills 
    the requirements of an \ref ImageIterator, so that it can be used to
    simulate an image whose pixels' values equal their coordinates. This
    secondary usage is explained on page \ref CoordinateIterator.
    
    Standard usage as a difference vector is mainly needed in the context
    of images. For example, Diff2D may be used as an index for <TT>operator[]<TT>:
    
    \code
    vigra::Diff2D location(...);
    
    value = image[location];    
    \endcode
    
    This is especially important in connection with accessors, where the
    offset variant of <TT>operator()<TT> takes only one offset object:
    
    \code
    // accessor(iterator, dx, dy); is not allowed
    value = accessor(iterator, vigra::Diff2D(dx, dy));
    \endcode
    
    
    Diff2D is also returned by <TT>image.size()<TT>, so that we can create 
    new images by calculating their size using Diff2D's arithmetic 
    functions:
    
    \code
    // create an image that is 10 pixels smaller in each direction
    Image new_image(old_image.size() - Diff2D(10,10));  
    \endcode 
    
    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"
*/
class Diff2D
{
  public:
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

       /** Create vector by adding specified offset.
        */
    Diff2D operator+(Diff2D const & offset) const
    {
        return Diff2D(x + offset.x, y + offset.y);
    }
    
        /** Create vector by subtracting specified offset.
        */
    Diff2D operator-(Diff2D const & offset) const
    {
        return Diff2D(x - offset.x, y - offset.y);
    }
    
        /** Calculate length of difference vector.
        */
    double magnitude() const
    {
#ifndef CMATH_NOT_IN_STD
        return std::sqrt((double)(x*x + y*y));
#else
        return sqrt((double)(x*x + y*y));
#endif
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
    Diff2D & operator*()
    {
        return *this;
    }
    
        /** Read current coordinate.
        */
    Diff2D operator*() const
    {
        return *this;
    }
    
        /** Read coordinate at an offset.
        */
    Diff2D operator()(int const & dx, int const & dy) const
    {
        return Diff2D(x + dx, y + dy);
    }

        /** Read coordinate at an offset.
        */
    Diff2D operator[](Diff2D const & offset) const
    {
        return Diff2D(x + offset.x, y + offset.y);
    }
    
        /** the iterator's value type
        */
    typedef Diff2D value_type;
        /** the iterator's PixelType
        */
    typedef Diff2D PixelType;
        /** type of the x-navigator
        */
    typedef int MoveX;
        /** type of the y-navigator
        */
    typedef int MoveY;
};


/********************************************************/
/*                                                      */
/*                      Dist2D                          */
/*                                                      */
/********************************************************/

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


/*! \page TupleTypes Tuple Types 

    pair, triple, tuple4, tuple5

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"
    
    VIGRA defines tuple types \p vigra::triple, \p vigra::tuple4, \p vigra::tuple5. 
    In addition, \p std::pair is imported into namespace vigra from the C++ standard 
    library. All these types are defined similarly:
    
    <ul>
    
    <li> They are parameterized by the respective number of types. For each tuple,
    a constructor is defined that takes that many arguments, e.g.:
    \code
    template <class First, class Second, class Third>
    class Triple { ... };
    \endcode
    </li>
    <li> A number of \p typedef's tells the types stored in the tuple:
    
    \code
    typedef ... first_type; 
    typedef ... second_type; 
    typedef ... third_type;  // triple, tuple4, tuple5 only
    typedef ... forth_type;  // tuple4, tuple5 only
    typedef ... fifth_type;  // tuple5 only
    \endcode
    </li>
    <li> Items are stored in the following public attributes:
    
    \code
    
    first; 
    second; 
    third;  // triple, tuple4, tuple5 only
    forth;  // tuple4, tuple5 only
    fifth;  // tuple5 only
    
    \endcode
    </li>
    </ul>

    
*/

/********************************************************/
/*                                                      */
/*                          pair                        */
/*                                                      */
/********************************************************/

using std::pair;

/********************************************************/
/*                                                      */
/*                          triple                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3>
struct triple {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;
    triple() {}
    triple(const T1& a, const T2& b, const T3& c) 
    : first(a), second(b), third(c) {}
};

/********************************************************/
/*                                                      */
/*                          tuple4                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3, class T4>
struct tuple4 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    tuple4() {}
    tuple4(const T1& a, const T2& b, const T3& c, const T4& d) 
    : first(a), second(b), third(c), fourth(d) {}
};

/********************************************************/
/*                                                      */
/*                          tuple5                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3, class T4, class T5>
struct tuple5 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;
    typedef T5 fifth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    T5 fifth;
    tuple5() {}
    tuple5(const T1& a, const T2& b, const T3& c, const T4& d, const T5& e) 
    : first(a), second(b), third(c), fourth(d), fifth(e) {}
};


} // namespace vigra

/*! \page MathConstants Mathematical Constants 

    <TT>M_PI, M_SQRT2</TT>

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"
    
    Since <TT>M_PI</TT> and <TT>M_SQRT2</TT> are not officially standardized, 
    we provide definitions here for those compilers that don't support them.   
    
    \code
    #ifndef M_PI
    #    define M_PI     3.14159265358979323846
    #endif

    #ifndef M_SQRT2
    #    define M_SQRT2  1.41421356237309504880
    #endif
    \endcode
*/
#ifndef M_PI
#    define M_PI     3.14159265358979323846
#endif

#ifndef M_SQRT2
#    define M_SQRT2  1.41421356237309504880
#endif

#endif // VIGRA_BASICS_HXX
