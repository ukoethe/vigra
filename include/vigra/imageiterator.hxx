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
 
 
#ifndef VIGRA_IMAGEITERATOR_HXX
#define VIGRA_IMAGEITERATOR_HXX

#include "vigra/utilities.hxx"
#include "vigra/accessor.hxx"
#include "vigra/iteratortraits.hxx"

namespace vigra { 

/** \addtogroup ImageIterators  Image Iterators

    \brief General image iterator implementations.
    
    These iterators may be used for any image data type that stores its
    data as a linear array of pixels. The array will be interpreted as a
    row-major matrix with a particular width.
*/
//@{
/********************************************************/
/*                                                      */
/*                      ImageIteratorBase               */
/*                                                      */
/********************************************************/

/** \brief Base class for 2D random access iterators.

    This class contains the nevigational part of the iterator.
    Use derived classes to specifiy the access to the pixels.

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"
    
    Namespace: vigra
    
    The usage examples assume that you constructed two iterators like
    this:

    \code
    vigra::ImageIteratorBase<SomePixelType> iterator(base, width);
    vigra::ImageIteratorBase<SomePixelType> iterator1(base, width);
    \endcode
    
    See the paper: U. Koethe: 
    <a href="documents/GenericProg2D.ps">Reusable Algorithms in Image Processing</a> 
    for a discussion of the concepts behind ImageIterators.
    
*/
template <class PIXELTYPE>
class ImageIteratorBase
{
  public:
    typedef PIXELTYPE PixelType;
    typedef PIXELTYPE value_type;
    
        /** Let operations act in X direction
        */
    class MoveX
    {
      public:
             /** Makes x coordinate identical to rhs, but lets y untouched.
                 Result is undefined if rhs points to different image. <br>
                 Usage: <TT> iterator.x = iterator1.x </TT>
             */
        MoveX & operator=(MoveX const & rhs)
        {
            current_ = rhs.current_;
            return *this;
        }
    
        /** @name Increment and Decrement X */
        //@{
            /** usage: <TT> ++iterator.x </TT> 
            */
        void operator++() {++current_;}
            /** usage: <TT> iterator.x++ </TT>
            */
        void operator++(int) {++current_;}
            /** usage: <TT> --iterator.x </TT>
            */
        void operator--() {--current_;}
            /** usage: <TT> iterator.x-- </TT>
            */
        void operator--(int) {--current_;}
    
        //@}

        /** @name Random Offset in X */
        //@{
            /** usage: <TT> iterator.x += 10 </TT>.
            */
        void operator+=(int dx) {current_ += dx; }
            /** usage: <TT> iterator.x -= 10 </TT>.
            */
        void operator-=(int dx) {current_ -= dx; }

        //@}
        /** @name Comparison of X Coordinates */
        //@{
            /** usage: <TT> iterator.x == iterator1.x </TT>
            */
        bool operator==(MoveX const & rhs) const 
         { return current_ == rhs.current_; }
            /** usage: <TT> iterator.x != iterator1.x </TT>
            */
        bool operator!=(MoveX const & rhs) const 
         { return current_ != rhs.current_; }
            /** usage: <TT> iterator.x < iterator1.x </TT>
            */
        bool operator<(MoveX const & rhs) const 
         { return current_ < rhs.current_; }

        //@}
        /** @name Difference of X Coordinates */
        //@{
            /** usage: <TT> int dx = iterator.x - iterator1.x </TT>
            */
        int operator-(MoveX const & rhs) const 
         { return current_ - rhs.current_; }
        //@}
      
      private:
        friend class ImageIteratorBase<PixelType>;

        MoveX(PixelType * base)
        : current_(base)
        {}

        MoveX(MoveX const & rhs)
        : current_(rhs.current_)
        {}

        PixelType * current_;
    };

        /** Let operations act in Y direction
        */
    class MoveY
    {
      public:
      
             /** Makes y coordinate identical to rhs, but lets x untouched.
                 Result is undefined if rhs points to different image. <br>
                 Usage: <TT> iterator.y = iterator1.y </TT>
             */
        MoveY & operator=(MoveY const & rhs)
        {
            width_ = rhs.width_;
            offset_ = rhs.offset_;
            return *this;
        }
    
        /** @name Increment and Decrement Y */
        //@{
            /** usage: <TT> ++iterator.y </TT>
            */
        void operator++() {offset_ += width_; }
            /** usage: <TT> iterator.y++ </TT>
            */
        void operator++(int) {offset_ += width_; }
            /** usage: <TT> --iterator.y </TT>
            */
        void operator--() {offset_ -= width_; }
            /** usage: <TT> iterator.y-- </TT>
            */
        void operator--(int) {offset_ -= width_; }
    
        //@}

        /** @name Random Offset in Y Direction */
        //@{
            /** usage: <TT> iterator.y += 20 </TT>
            */
        void operator+=(int dy) {offset_ += dy*width_; }
            /** usage: <TT> iterator.y -= 20 </TT>
            */
        void operator-=(int dy) {offset_ -= dy*width_; }

        //@}

        /** @name Comparison of Y Coordinates */
        //@{
            /** usage: <TT> iterator.y == iterator1.y </TT>
            */
        bool operator==(MoveY const & rhs) const 
         { return (offset_ == rhs.offset_); }
            /** usage: <TT> iterator.y != iterator1.y </TT>
            */
        bool operator!=(MoveY const & rhs) const 
         { return (offset_ != rhs.offset_); }
            /** usage: <TT> iterator.y < iterator1.y </TT>
            */
        bool operator<(MoveY const & rhs) const 
         { return (offset_ < rhs.offset_); }

        //@}

        /** @name Difference of Y Coordinates */
        //@{
            /** usage: <TT> int dy = iterator.y - iterator1.y </TT>
            */
        int operator-(MoveY const & rhs) const 
         { return (offset_ - rhs.offset_) / width_; }

        //@}

      private:
        friend class ImageIteratorBase<PixelType>;

            // only construct within ImageIteratorBase<PixelType>
        MoveY(int width)
        : width_(width),
          offset_(0)
        {}

            // only construct within ImageIteratorBase<PixelType>
        MoveY(MoveY const & rhs)
        : width_(rhs.width_),
          offset_(rhs.offset_)
        {}
    
        int offset_;
        int width_;
    };
    
    /** @name Comparison of Iterators */
    //@{
        /** usage: <TT> iterator == iterator1 </TT>
        */
    bool operator==(ImageIteratorBase const & rhs) const
    {
        return (x == rhs.x) && (y == rhs.y);
    }
    
        /** usage: <TT> iterator != iterator1 </TT>
        */
    bool operator!=(ImageIteratorBase const & rhs) const
    {
        return (x != rhs.x) || (y != rhs.y);
    }
    
        /** usage: <TT> Diff2D dist = iterator - iterator1 </TT>
        */
    Diff2D operator-(ImageIteratorBase const & rhs) const
    {
        return Diff2D(x - rhs.x, y - rhs.y);
    }
    
    //@}
    
    /** @name Specify coordinate to operate on */
    //@{
        /// refer to x coordinate
    MoveX x;
        /// refer to y coordinate
    MoveY y;
    //@}
    
  protected:
        /** Construct from raw memory with 'width' of physical lines.
        If the raw memory is encapsulated in an image object this
        object should have a fatory function that constructs the
        iterator.
        */
    ImageIteratorBase(PixelType * base, int width)
    : x(base),
      y(width)
    {}
    
        /** Copy constructor */
    ImageIteratorBase(ImageIteratorBase const & rhs)
    : x(rhs.x),
      y(rhs.y)
    {}
    
        /** Default constructor */
    ImageIteratorBase()
    : x(0),
      y(0)
    {}
    
        /** Copy assignment */
    ImageIteratorBase & operator=(ImageIteratorBase const & rhs)
    {
        if(this != &rhs)
        {
            x = rhs.x;
            y = rhs.y;
        }
        return *this;
    }
    
        /** Add offset via Diff2D
        */
    inline ImageIteratorBase & operator+=(Diff2D const & s)
    {
        x += s.x;
        y += s.y;
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    inline ImageIteratorBase & operator-=(Diff2D const & s)
    {
        x -= s.x;
        y -= s.y;
        return *this;
    }
    
        /** access current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current() </TT>
        */
    inline PixelType & current()
    {
        return *(x.current_ + y.offset_);
    }
    
        /** access pixel at offset (dx, dy) from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(dx, dy) </TT>
        */
    inline PixelType & current(int const & dx, int const & dy)
    {
        return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }
    
        /** access pixel at offset d from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(Diff2D(1,1)) </TT>
        */
    inline PixelType & current(Diff2D const & d)
    {
        return *(x.current_ + d.x + y.offset_ + y.width_ * d.y);
    }
    
        /** read current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current() </TT>
        */
    inline PixelType const & current() const
    {
        return *(x.current_ + y.offset_);
    }

        /** read pixel with offset (dx, dy) from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(dx, dy) </TT>
        */
    inline PixelType const & current(int const & dx, int const & dy) const
    {
        return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }

        /** read pixel with offset d from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(Diff2D(1,1)) </TT>
        */
    inline PixelType const & current(Diff2D const & d) const
    {
        return *(x.current_ + d.x + y.offset_ + y.width_ * d.y);
    }
    
    inline int width() const { return y.width_; }
}; 
    
/********************************************************/
/*                                                      */
/*                      ImageIterator                   */
/*                                                      */
/********************************************************/

/** \brief Standard 2D random access iterator for images that store the
    data in a linear array.
    
    Adds <TT>%operator*()</TT>, <TT>operator()(dx, dy)</TT>, <TT>operator[](dy)</TT> and <TT>operator[](Diff2D)</TT>
    to access the current
    pixel and the pixel at distance (dx, dy).
    The usage examples assume that you constructed two iterators like
    this:
    
    \code
    vigra::ImageIteratorBase<SomePixelType> iterator(base, width);
    vigra::ImageIteratorBase<SomePixelType> iterator1(base, width);
    \endcode
    
    See the paper: U. Koethe: 
    <a href="documents/GenericProg2D.ps">Reusable Algorithms in Image Processing</a> 
    for a discussion of the concepts behind ImageIterators.

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"
    
    Namespace: vigra
    
*/    
template <class PIXELTYPE>
class ImageIterator: public ImageIteratorBase<PIXELTYPE>
{
  public:

        /** The underlying image's pixel type.
        */
    typedef PIXELTYPE PixelType;
    
        /** Construct from raw memory with 'offset' between lines.
        If the raw memory is encapsulated in an image object this
        object should have a fatory function that constructs the
        iterator.
        */
    inline ImageIterator(PixelType * base, int offset)
    : ImageIteratorBase<PIXELTYPE>(base, offset)
    {}
    
        /** Copy constructor */
    inline ImageIterator(const ImageIterator & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Default constructor */
    inline ImageIterator()
    : ImageIteratorBase<PIXELTYPE>()
    {}    
    
        /** Copy assignment */
    inline ImageIterator & operator=(const ImageIterator & rhs)
    {
        if(this != &rhs)
        {
            ImageIteratorBase<PIXELTYPE>::operator=(rhs);
        }
        return *this;
    }

    /** @name Computed Assignment */
    //@{
        /** Add offset via Diff2D
        */
    inline ImageIterator & operator+=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator+=(s);
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    inline ImageIterator & operator-=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator-=(s);
        return *this;
    }
    
    //@}
   
    Diff2D operator-(ImageIterator const & rhs) const
    {
        return Diff2D(x - rhs.x, y - rhs.y);
    }

    /** @name Construct iterator at a distance */
    //@{
        /** Add a distance
        */
    inline ImageIterator operator+(Diff2D const & s) const
    {
        ImageIterator ret(*this);

        ret += s;

        return ret;
    }
        /** Subtract a distance
        */
    inline ImageIterator operator-(Diff2D const & s) const
    {
        ImageIterator ret(*this);

        ret -= s;

        return ret;
    }
   //@}

    /** @name Access the Pixels */
    //@{
        /** Access current pixel. <br>
            usage: <TT> SomePixelType value = *iterator </TT>
        */
    inline PixelType & operator*()
    {
        return this->current();
    }
    
        /** Read current pixel. <br>
            usage: <TT> SomePixelType value = *iterator </TT>
        */
    inline PixelType const & operator*() const
    {
        return this->current();
    }
    
        /** Call member of current pixel. <br>
            usage: <TT> iterator->pixelMemberFunction() </TT>
        */
    inline PixelType * operator->()
    {
        return &(this->current)();
    }
    
        /** Call const member of current pixel. <br>
            usage: <TT> iterator->constPixelMemberFunction() </TT>
            Use this operation in inner loop.
        */
    inline PixelType const * operator->() const
    {
        return &(this->current)();
    }
    
        /** Access pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    inline PixelType & operator[](Diff2D const & d)
    {
        return this->current(d);
    }
    
        /** Read pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    inline PixelType const & operator[](Diff2D const & d) const
    {
        return this->current(d);
    }
    
        /** Access pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    inline PixelType & operator()(int const & dx, int const & dy)
    {
        return this->current(dx, dy);
    }
    
        /** Read pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    inline PixelType const & operator()(int const & dx, int const & dy) const
    {
        return this->current(dx, dy);
    }
    
        /** Access pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    inline PixelType * operator[](int const & dy)
    {
        return &(this->current())+dy*width();
    }
    
        /** Read pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    inline PixelType const * operator[](int const & dy) const
    {
        return &(this->current())+dy*width();
    }
    //@}
    
};

/********************************************************/
/*                                                      */
/*                   ConstImageIterator                 */
/*                                                      */
/********************************************************/

/** \brief Standard 2D random access const iterator for images that
    store the data as a linear array.
    
    Adds <TT>operator*()</TT>, <TT>operator()(dx, dy)</TT>, and <TT>operator[](dy)</TT>
    to access the current
    pixel and the pixel at distance (dx, dy).
    The usage examples assume that you constructed two iterators like
    this:
    
    \code
    vigra::ConstImageIterator<SomePixelType> iterator(base, width) 
    vigra::ConstImageIterator<SomePixelType> iterator1(base, width)
    \endcode

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"
    
    Namespace: vigra
    
*/
template <class PIXELTYPE>
class ConstImageIterator: public ImageIteratorBase<PIXELTYPE>
{
  public:

        /** The underlying image's pixel type.
        */
    typedef PIXELTYPE PixelType;
    
        /** Construct from raw memory with 'offset' between lines.
        If the raw memory is encapsulated in an image object this
        object should have a fatory function that constructs the
        iterator.
        */
    inline ConstImageIterator(PixelType * base, int offset)
    : ImageIteratorBase<PIXELTYPE>(base, offset)
    {}
    
        /** Copy constructor */
    inline ConstImageIterator(const ConstImageIterator & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Defult constructor */
    inline ConstImageIterator()
    : ImageIteratorBase<PIXELTYPE>()
    {}    
    
        /** Constructor from mutable iterator */
    inline ConstImageIterator(const ImageIterator<PIXELTYPE> & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Copy assignment */
    inline ConstImageIterator & operator=(const ConstImageIterator & rhs)
    {
        if(this != &rhs)
        {
            ImageIteratorBase<PIXELTYPE>::operator=(rhs);
        }
        return *this;
    }
    
        /** Assignment from mutable iterator */
    inline ConstImageIterator & 
    operator=(const ImageIterator<PIXELTYPE> & rhs)
    {
        ImageIteratorBase<PIXELTYPE>::operator=(rhs);
        return *this;
    }

    /** @name Computed Assignment */
    //@{
        /** Add offset via Diff2D
        */
    inline ConstImageIterator & operator+=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator+=(s);
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    inline ConstImageIterator & operator-=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator-=(s);
        return *this;
    }
    
    //@}
    
    Diff2D operator-(ConstImageIterator const & rhs) const
    {
        return Diff2D(x - rhs.x, y - rhs.y);
    }

    /** @name Construct iterator at a distance */
    //@{
        /** Add a distance
        */
    inline ConstImageIterator operator+(Diff2D const & s) const
    {
        ConstImageIterator ret(*this);
    
        ret += s;
    
        return ret;
    }
        /** Subtract a distance
        */
    inline ConstImageIterator operator-(Diff2D const & s) const
    {
        ConstImageIterator ret(*this);

        ret -= s;

        return ret;
    }
    //@}

    /** @name Access the Pixels */
    //@{
        /** Read current pixel. <br>
            usage: <TT> SomePixelType value = *iterator </TT>
        */
    inline PixelType const & operator*() const
    {
       return this->current();
    }
    
        /** Call member of current pixel. <br>
            usage: <TT> iterator->pixelMemberFunction() </TT>
        */
    inline PixelType const * operator->() const
    {
        return &(this->current)();
    }
    
        /** Read pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    inline PixelType const & operator()(int const & dx, int const & dy) const
    {
       return this->current(dx, dy);
    }

        /** Read pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    inline PixelType const & operator[](Diff2D const & d) const
    {
        return this->current(d);
    }
    
        /** Read pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    inline PixelType const * operator[](int const & dy) const
    {
        return &(this->current())+dy*width();
    }
    //@}
    
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T> 
struct IteratorTraits<ImageIterator<T> > 
{
    typedef StandardAccessor<T> DefaultAccessor;
};

template <class T> 
struct IteratorTraits<ConstImageIterator<T> > 
{
    typedef StandardConstAccessor<T> DefaultAccessor;
};

#endif

/********************************************************/
/*                                                      */
/*                 ConstValueIterator                   */
/*                                                      */
/********************************************************/

/** \brief Iterator that always returns the constant specified in the 
    constructor. 
    
    This iterator can be used to simulate an image that
    does not actually exist.
    
    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"
    
    Namespace: vigra
    
*/
template <class PIXELTYPE>
class ConstValueIterator
{
  public:
        /** the iterator's PixelType
        */
    typedef PIXELTYPE PixelType;
    
        /** Let operations act in X direction
        */
    typedef int MoveX;

        /** Let operations act in Y direction
        */
    typedef int MoveY;
    
        /** Default Constructor. (the constant is set to
        <TT>NumericTraits<PIXELTYPE>::zero()</TT> )
        */
    ConstValueIterator()
    : value_(NumericTraits<PIXELTYPE>::zero()), x(0), y(0)
    {}
    
        /** Construct with given constant.
        */
    ConstValueIterator(PixelType const & v)
    : value_(v), x(0), y(0)
    {}
    
        /** Copy Constructor.
       */
    ConstValueIterator(ConstValueIterator const & v)
    : value_(v.value_), x(v.x), y(v.y)
    {}
    
        /** Copy Assigment.
        */
    ConstValueIterator & operator=(ConstValueIterator const & v)
    {
        if(this != &v)
    {
        value_ = v.value_;
        x = v.x;
        y = v.y;
    }
    return *this;
    }
    
        /** Move iterator by specified distance.
        */
    ConstValueIterator & operator+=(Diff2D const & d)
    {
        x += d.x;
        y += d.y;
        return *this;
    }
    
        /** Move iterator by specified distance.
        */
    ConstValueIterator & operator-=(Diff2D const & d)
    {
        x -= d.x;
        y -= d.y;
        return *this;
    }
        
        /** Create iterator at specified distance.
        */
    ConstValueIterator operator+(Diff2D const & d) const
    {
        ConstValueIterator ret(*this);
        ret += d;
        return ret;
    }
    
        /** Create iterator at specified distance.
        */
    ConstValueIterator operator-(Diff2D const & d) const
    {
        ConstValueIterator ret(*this);
        ret -= d;
        return ret;
    }
    
        /** Compute distance between two iterators
        */
    Diff2D operator-(ConstValueIterator const & r) const
    {
        return Diff2D(x - r.x, y - r.y);
    }
    
        /** Equality.
        */
    bool operator==(ConstValueIterator const & r) const
    {
        return (x == r.x) && (y == r.y);
    }
    
        /** Inequality.
        */
    bool operator!=(ConstValueIterator const & r) const
    {
        return (x != r.x) || (y != r.y);
    }
    
        /** Read current pixel (return specified constant).
        */
    PixelType const & operator*() const
    {
        return value_;
    }
    
        /** Read pixel at a distance (return specified constant).
        */
    PixelType const & operator()(int const &, int const &) const
    {
        return value_;
    }

        /** Read pixel at a distance (return specified constant).
        */
    PixelType const & operator[](Diff2D const &) const
    {
        return value_;
    }
  
    /** @name Specify coordinate direction for navigation commands */
    //@{
        /// refer to x coordinate
    int x;
        /// refer to y coordinate
    int y;
    //@}
    
  private:
  
    PixelType value_;
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T> 
struct IteratorTraits<ConstValueIterator<T> > 
{
    typedef StandardConstAccessor<T> DefaultAccessor;
};

#endif

typedef Diff2D CoordinateIterator;

/** \class CoordinateIterator 

    This used to be a separate class, 
    but has now become an alias for \ref vigra::Diff2D. This is possible because
    Diff2D now provides all the necessary functionality. 
    
    CoordinateIterator behaves like a read-only \ref ImageIterator for 
    an image in which each pixel contains its coordinate. This is useful for
    algorithms that need access to the current pixel's location.
    For example, you can use CoordinateIterator/Diff2D to 
    find the center of mass of an image region. To implement this,
    we first need a functor for center-of-mass calculations:

    \code
    
    struct CenterOfMassFunctor
    {
        CenterOfMassFunctor()
        : x(0.0), y(0.0), size(0)
        {}
        
        void operator()(Diff2d const& diff)
        {
            ++size;
            x += diff.x;
            y += diff.y;
        }
        
        float xCenter() const
        {   return x / size; }
        
        float yCenter() const
        {   return y / size; }       
        
        float x;
        float y;
        int size;
    };
    \endcode
    
    Using this functor, we find the center of mass like so:
    
    \code
    vigra::BImage img(w,h);
    ... // mark a region in the image with '1', background with '0'
    
    CenterOfMassFunctor center;
    
    vigra::inspectImageIf(
        srcIterRange(Diff2D(), Diff2D() + img.size()),
        srcImage(img),
        center);
        
    std::cout << "Center of mass: " << center.xCenter() << 
                                ", " << center.yCenter() << std::endl;
    \endcode    
    
    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"
    
    Namespace: vigra
    
    \brief Simulate an image where each pixel contains its coordinate
*/

//@}

} // namespace vigra

#endif // VIGRA_IMAGEITERATOR_HXX
