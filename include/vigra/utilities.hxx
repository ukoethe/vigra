/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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

#include <math.h>     // for sqrt()
#include <utility>    // for pair

#include "vigra/config.hxx"
#include "vigra/error.hxx"

struct VigraTrueType 
{
    static const bool asBool = true;
};
struct VigraFalseType
{
    static const bool asBool = false;
};

/********************************************************/
/*                                                      */
/*                      Diff2D                          */
/*                                                      */
/********************************************************/

/** Basic helper functionality needed throughout.
    @name Utilities
*/
//@{

/** Two dimensional difference vector.
    
    This class acts primarily as a difference vector for specifying 
    pixel coordinates and region sizes. In addition, Diff2D fulfills 
    the requirements of an \Ref{ImageIterator}, so that it can be used to
    simulate an image whose pixels' values equal their coordinates. This
    secondary usage is explained on page \Ref{CoordinateIterator}.
    
    Standard usage as a difference vector is mainly needed in the context
    of images. For example, Diff2D may be used as an index for #operator[]#:
    
    \begin{verbatim}
    Diff2D location(...);
    
    value = image[location];    
    \end{verbatim}
    
    This is especially important in connection with accessors, where the
    offset variant of #operator()# takes only one offset object:
    
    \begin{verbatim}
    // accessor(iterator, dx, dy); is not allowed
    value = accessor(iterator, Diff2D(dx, dy));
    \end{verbatim}
    
    
    Diff2D is also returned by #image.size()#, so that we can create 
    new images by calculating their size using Diff2D's arithmetic 
    functions:
    
    \begin{verbatim}
    // create an image that is 10 pixels smaller in each direction
    Image new_image(old_image.size() - Diff2D(10,10));  
    \end{verbatim} 
    
    Include-File: \URL[vigra/utilities.hxx]{../include/vigra/utilities.hxx}
*/
class Diff2D
{
  public:
    /** @name Construction and Assignment */
    //@{
        /** Default Constructor. Init iterator at position (0,0)
	    @memo
	*/
    Diff2D()
    : x(0), y(0)
    {}
    
        /** Construct at given position.
	    @memo
	*/
    Diff2D(int ax, int ay)
    : x(ax), y(ay)
    {}
    
        /** Copy Constructor.
	    @memo
	*/
    Diff2D(Diff2D const & v)
    : x(v.x), y(v.y)
    {}
    
        /** Copy Assigment.
	    @memo
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
    
        /** Increase coordinate by specified offset.
	    @memo
	*/
    Diff2D & operator+=(Diff2D const & offset)
    {
        x += offset.x;
	y += offset.y;
        return *this;
    }
    
        /** Decrease coordinate by specified vector.
	    @memo
	*/
    Diff2D & operator-=(Diff2D const & offset)
    {
        x -= offset.x;
	y -= offset.y;
        return *this;
    }
    //@}
    
    /** @name Methods */
    //@{

        /** Create vector by adding specified offset.
	    @memo
	*/
    Diff2D operator+(Diff2D const & offset) const
    {
        return Diff2D(x + offset.x, y + offset.y);
    }
    
        /** Create vector by subtracting specified offset.
	    @memo
	*/
    Diff2D operator-(Diff2D const & offset) const
    {
        return Diff2D(x - offset.x, y - offset.y);
    }
    
        /** Calculate length of difference vector.
	    @memo
	*/
    double magnitude() const
    {
        return sqrt(x*x + y*y);
    }
    
        /** Equality.
	    @memo
	*/
    bool operator==(Diff2D const & r) const
    {
        return (x == r.x) && (y == r.y);
    }
    
        /** Inequality.
	    @memo
	*/
    bool operator!=(Diff2D const & r) const
    {
        return (x != r.x) || (y != r.y);
    }
    
    //@}
  
    /** @name Data members. */
    //@{
    
        /** Used for both access to the current x-coordinate {\em and}
            to specify that an iterator navigation command is to be
            applied in x-direction. \\
            usage:  # x = diff2d.x # (use Diff2D::x as component of difference vector) \\
            or #&nbsp; ++diff.x &nbsp;# (use Diff2D as iterator, move right)
            @memo
	*/
    int x;
        /** Used for both access to the current y-coordinate {\em and}
            to specify that an iterator navigation command is to be
            applied in y-direction. \\
            usage:  # y = diff2d.y # (use Diff2D::y as component of difference vector) \\
            or #&nbsp; ++diff.y &nbsp;# (use Diff2D as iterator, move down)
            @memo
	*/
    int y;
    
    //@}
    
    /** @name Functionality for use as ImageIterator. */
    //@{
  
        /** Access current coordinate.
	    @memo
	*/
    Diff2D & operator*()
    {
        return *this;
    }
    
        /** Read current coordinate.
	    @memo
	*/
    Diff2D operator*() const
    {
        return *this;
    }
    
        /** Read coordinate at an offset.
	    @memo
	*/
    Diff2D operator()(int const & dx, int const & dy) const
    {
        return Diff2D(x + dx, y + dy);
    }

        /** Read coordinate at an offset.
	    @memo
	*/
    Diff2D operator[](Diff2D const & offset) const
    {
        return Diff2D(x + offset.x, y + offset.y);
    }

    /** @name an ImageIterator's required embedded types
    */
    //@{
        /** the iterator's value type
	    @memo
	*/
    typedef Diff2D value_type;
    typedef Diff2D PixelType;
        /** type of the x-navigator
	    @memo
	*/
    typedef int MoveX;
        /** type of the y-navigator
	    @memo
	*/
    typedef int MoveY;
    
    //@}
    
    //@}
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



/********************************************************/
/*                                                      */
/*                          pair                        */
/*                                                      */
/********************************************************/

#if !defined(NO_NAMESPACE_STD)
using std::pair;
#endif

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

/** @name Tuple Types

    VIGRA defines tuple types #triple#, #tuple4#, #tuple5#. In addition,
    #pair# is imported from the C++ standard library. All these types are 
    defined similarly:
    
    \begin{itemize}
    
    \item They are parameterized by the respective number of types. For each tuple,
    a constructor is defined that takes that many arguments.
    
    \item A number of #typedef#s tells the types stored in the tuple:
    
    \begin{verbatim}
    
    typedef ... first_type; 
    typedef ... second_type; 
    typedef ... third_type;  // triple, tuple4, tuple5 only
    typedef ... forth_type;  // tuple4, tuple5 only
    typedef ... fifth_type;  // tuple5 only
    
    \end{verbatim}
    
    \item Items are stored in the following public attributes:
    
    \begin{verbatim}
    
    first; 
    second; 
    third;  // triple, tuple4, tuple5 only
    forth;  // tuple4, tuple5 only
    fifth;  // tuple5 only
    
    \end{verbatim}
    
    \end{itemize}

    Include-File: \URL[vigra/utilities.hxx]{../include/vigra/utilities.hxx}

    @memo pair, triple, tuple4, tuple5
*/

/** @name Dist2D

    This class is deprecated - use \Ref{Diff2D} instead.
*/

//@}

#endif // VIGRA_BASICS_HXX
