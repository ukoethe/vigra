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
#include <iterator>   // iterator tags

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
     \ref vigra::IteratorAdaptor 
     <DD><em>Quickly create STL-compatible 1D iterator adaptors</em>
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
/*                    IteratorAdaptor                   */
/*                                                      */
/********************************************************/

/*! \brief Quckly create 1-dimensional iterator adapters.
    
    This class supports the easy creation of 1D iterator adpaters out 
    of existing iterators. To use it, you must first implement a policy class
    that defines the iterator's behavior. The policy is used to 
    instantiate the IteratorAdapter template, which thus automatically
    obtains all required functions of an STL-compatible iterator.
    General information on how this works can be found on the 
    <a href="http://www.boost.org/libs/utility/iterator_adaptors.htm">Boost Iterator Adaptor</a>
    page, although there are some differences in the details of the
    boost and VIGRA implementations.
    Here is an example policy class that just exports the behaviour
    of the underlying iterator:
    
    \code
    template <class Iterator>
    class TrivialIteratorAdaptorPolicy
    {
      public:
        // the underlying iterator
        typedef Iterator                               BaseType;
        
        // the adpator's value type
        typedef typename Iterator::value_type          value_type;
        
        // the adpator's difference type (result of 'iter1 - iter2',
        //                                argument of 'iter[n]')
        typedef typename Iterator::difference_type     difference_type;
        
        // the adpator's reference type (result of '*iter')
        typedef typename Iterator::reference           reference;
        
        // the adpator's index_reference type (result of 'iter[n]')
        typedef typename Iterator::index_reference     index_reference;
        
        // the adpator's pointer type (result of 'iter.operator->()')
        typedef typename Iterator::pointer             pointer;
        
        // the adpator's iterator category
        typedef typename Iterator::iterator_category   iterator_category;

        // do some additional initialization in the adaptor's constructor
        static void initialize(BaseType & d) {}

        // called by '*iter', 'iter->'
        static reference dereference(BaseType const & d)
            { return *d; }

        // called by 'iter[n]'
        static index_reference dereference(BaseType d, difference_type n)
            { return d[n]; }

        // called by 'iter1 == iter2', 'iter1 != iter2'
        static bool equal(BaseType const & d1, BaseType const & d2)
            { return d1 == d2; }

        // called by 'iter1 < iter2', 'iter1 <= iter2', 'iter1 > iter2', 'iter1 >= iter2'
        static bool less(BaseType const & d1, BaseType const & d2)
            { return d1 < d2; }

        // called by 'iter1 - iter2'
        static difference_type difference(BaseType const & d1, BaseType const & d2)
            { return d1 - d2; }

        // called by '++iter', 'iter++'
        static void increment(BaseType & d)
            { ++d; }

        // called by '--iter', 'iter--'
        static void decrement(BaseType & d)
            { --d; }

        // called by 'iter += n', 'iter -= n'
        static void advance(BaseType & d, difference_type n)
            { d += n; }
    };
    \endcode
    
    This policy class is used like this:
    
    \code
    SomeIterator iter = ...;
    
    vigra::IteratorAdaptor<vigra::TrivialIteratorAdaptorPolicy<SomeIterator> > iter_adaptor(iter);
    \endcode
    
    By changing the definition of the policy members, a wide range of 
    adaptor behaviors can be achieved. If the base iterator isn't a
    random access iterator, just drop the functions that cannot be implemented.
    This simply means that some adaptor functions may not be called,
    as one would expect from an iterator that doesn't support random access.
    Note also that the <TT>BaseType</TT> needs not be an iterator -
    it can be any type that contains the information necessary for the
    adaptor to do it's work.
    
    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra

*/
template <class Policy>
class IteratorAdaptor
{
  public:
  
    typedef typename Policy::BaseType BaseType;
    typedef typename Policy::value_type        value_type;
    typedef typename Policy::difference_type   difference_type;
    typedef typename Policy::reference         reference;
    typedef typename Policy::index_reference   index_reference;
    typedef typename Policy::pointer           pointer;
    typedef typename Policy::iterator_category iterator_category;
    
    IteratorAdaptor()
    : adaptee_()
    {}
    
        /** Construct from an instance of the policy class' BaseType
            Note that the functions of the adaptor implement the
            interface of an random access iterator as defined in the
            C++ standard, so there is no need for explicit documentation.
        */
    explicit IteratorAdaptor(BaseType const & o)
    : adaptee_(o)
    {
        Policy::initialize(adaptee_);
    }
    
    IteratorAdaptor(IteratorAdaptor const & o)
    : adaptee_(o.adaptee_)
    {}
    
    IteratorAdaptor & operator=(BaseType const & o)
    {
        if(this != &o)
        {
            adaptee_ = o;
            Policy::initialize(adaptee_);
        }
        return *this;
    }
    
    IteratorAdaptor & operator=(IteratorAdaptor const & o)
    {
        if(this != &o)
            adaptee_ = o.adaptee_;
        return *this;
    }
    
    IteratorAdaptor & operator+=(difference_type d)
    {
        Policy::advance(adaptee_, d);
        return *this;
    }
    
    IteratorAdaptor operator+(difference_type d) const
    {
        return IteratorAdaptor(*this) += d;
    }
    
    IteratorAdaptor & operator-=(difference_type d)
    {
        Policy::advance(adaptee_, -d);
        return *this;
    }
    
    IteratorAdaptor operator-(difference_type d) const
    {
        return IteratorAdaptor(*this) -= d;
    }
    
    IteratorAdaptor & operator++()
    {
        Policy::increment(adaptee_);
        return *this;
    }
    
    IteratorAdaptor operator++(int)
    {
        IteratorAdaptor res(*this);
        Policy::increment(adaptee_);
        return res;
    }
    
    IteratorAdaptor & operator--()
    {
        Policy::decrement(adaptee_);
        return *this;
    }
    
    IteratorAdaptor operator--(int)
    {
        IteratorAdaptor res(*this);
        Policy::decrement(adaptee_);
        return res;
    }
    
    bool operator==(IteratorAdaptor const & o) const
    {
        return Policy::equal(adaptee_, o.adaptee_);
    }
    
    bool operator!=(IteratorAdaptor const & o) const
    {
        return !Policy::equal(adaptee_, o.adaptee_);
    }
    
    bool operator<(IteratorAdaptor const & o) const
    {
        return Policy::less(adaptee_, o.adaptee_);
    }
    
    bool operator<=(IteratorAdaptor const & o) const
    {
        return !Policy::less(o.adaptee_, adaptee_);
    }
    
    bool operator>(IteratorAdaptor const & o) const
    {
        return Policy::less(o.adaptee_, adaptee_);
    }
    
    bool operator>=(IteratorAdaptor const & o) const
    {
        return !Policy::less(adaptee_, o.adaptee_);
    }
    
    difference_type operator-(IteratorAdaptor const & o) const
    {
        return Policy::difference(adaptee_, o.adaptee_);
    }
    
    reference operator*() const
    {
        return Policy::dereference(adaptee_);
    }
    
    index_reference operator[](difference_type d) const
    {
        return Policy::dereference(adaptee_, d);
    }
    
    pointer operator->() const
    {
        return &Policy::dereference(adaptee_);
    }
    
  protected:

    BaseType adaptee_;
};

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
    
    static void initialize(BaseType & d) {}
    
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
    
    static void initialize(BaseType & d) {}
    
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

struct image_traverser_tag {};

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
    
    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"<br>
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

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra
    
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
