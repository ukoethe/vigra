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
/*                      Dist2D                          */
/*                                                      */
/********************************************************/

/** Basic helper functionality needed throughout.
    @name Utilities
*/
//@{

/** Two dimensional distance type.

    Include-File: \URL[vigra/utilities.hxx]{../include/vigra/utilities.hxx}
*/
class Dist2D
{
  public:
  /** @name Construction and Assignment */
  //@{
    /** init to specified width and height.
        @memo
    */
    Dist2D(int the_width, int the_height)
    : width(the_width),
      height(the_height)
    {}
    
    /** copy constructor
        @memo
    */
    Dist2D(Dist2D const & s)
    : width(s.width),
      height(s.height)
    {}
    
    /** copy assignment
        @memo
    */
    Dist2D & operator=(Dist2D const & s)
    {
        if(this != &s)
        {
            width = s.width;
            height = s.height;
        }
        return *this;
    }
    
    /** add-assign
        @memo
    */
    Dist2D & operator+=(Dist2D const & s)
    {
        width += s.width;
        height += s.height;
    
        return *this;
    }
    
    /** binary addition
        @memo
    */
    Dist2D  operator+(Dist2D const & s) const
    {
        Dist2D ret(*this);
        ret += s;
    
        return ret;
    }
  //@}
    
  /** @name Data members */
  //@{
    /** width of a range
        @memo
    */
    int width;
    /** height of a range
        @memo
    */
    int height;
  //@}
};

class Diff2D
{
  public:
  public:
        /** the iterator's PixelType
	    @memo
	*/
    typedef Diff2D value_type;
    typedef Diff2D PixelType;
    
	/** Let operations act in X direction
	*/
    typedef int MoveX;

	/** Let operations act in Y direction
	*/
    typedef int MoveY;
    
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
    
        /** Move iterator by specified distance.
	    @memo
	*/
    Diff2D & operator+=(Diff2D const & d)
    {
        x += d.x;
	y += d.y;
        return *this;
    }
    
        /** Move iterator by specified distance.
	    @memo
	*/
    Diff2D & operator-=(Diff2D const & d)
    {
        x -= d.x;
	y -= d.y;
        return *this;
    }
    //@}
    
    /** @name Methods */
    //@{

        /** Create vector by adding specified offset.
	    @memo
	*/
    Diff2D operator+(Diff2D const & r) const
    {
        return Diff2D(x + r.x, y + r.y);
    }
    
        /** Create vector by subtracting specified offset.
	    @memo
	*/
    Diff2D operator-(Diff2D const & r) const
    {
        return Diff2D(x - r.x, y - r.y);
    }
    
        /** Create vector by subtracting specified offset.
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
    
        /** Read coordinate at a distance.
	    @memo
	*/
    Diff2D operator()(int const & dx, int const & dy) const
    {
        return Diff2D(x + dx, y + dy);
    }

        /** Read coordinate at a distance.
	    @memo
	*/
    Diff2D operator[](Diff2D const & d) const
    {
        return Diff2D(x + d.x, y + d.y);
    }
    //@}
  
    /** @name Specify coordinate direction for navigation commands */
    //@{
        /// refer to x coordinate
    int x;
        /// refer to y coordinate
    int y;
    //@}
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


//@}

#endif // VIGRA_BASICS_HXX
