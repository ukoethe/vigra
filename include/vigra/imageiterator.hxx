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


/** Basic image iterator implementations 
    @name Image Iterators  
*/
//@{

/********************************************************/
/*                                                      */
/*                      ImageIteratorBase               */
/*                                                      */
/********************************************************/

/** Base class for 2D random access iterators for images containing the 
    data as a plain C array.
    This class contains the nevigational part of the iterator.
    Use derived classes to specifiy the access to the pixels.

    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
    The usage examples assume that you constructed two iterators like
    this:

    \begin{verbatim}
    ImageIteratorBase<SomePixelType> iterator(base, width);
    ImageIteratorBase<SomePixelType> iterator1(base, width);
    \end{verbatim}
    
    See the paper: U. Koethe: 
    \URL[Reusable Algorithms in Image Processing]{documents/GenericProg2D.ps} 
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
        /** @name Assignment (no public constructors) */
	//@{
		 /** Makes x coordinate identical to rhs, but lets y untouched.
		     Result is undefined if rhs points to different image. \\
		     Usage: # iterator.x = iterator1.x #
		     @memo
		 */
	MoveX & operator=(MoveX const & rhs)
	{
	    current_ = rhs.current_;
	    return *this;
	}
	
	//@}
	
        /** @name Increment and Decrement X */
	//@{
		/** usage: # ++iterator.x # 
		    @memo
		*/
        void operator++() {++current_;}
		/** usage: # iterator.x++ #
		    @memo
		*/
        void operator++(int) {++current_;}
		/** usage: # --iterator.x #
		    @memo
		*/
        void operator--() {--current_;}
		/** usage: # iterator.x-- #
		    @memo
		*/
        void operator--(int) {--current_;}
	
	//@}
	
	/** @name Random Offset in X */
	//@{
		/** usage: # iterator.x += 10 #.
		    @memo
		*/
	void operator+=(int dx) {current_ += dx; }
		/** usage: # iterator.x -= 10 #.
		    @memo
		*/
	void operator-=(int dx) {current_ -= dx; }
		
	//@}
	/** @name Comparision of X Coordinates */
	//@{
		/** usage: # iterator.x == iterator1.x #
		    @memo
		*/
	bool operator==(MoveX const & rhs) const 
	 { return current_ == rhs.current_; }
		/** usage: # iterator.x != iterator1.x #
		    @memo
		*/
	bool operator!=(MoveX const & rhs) const 
	 { return current_ != rhs.current_; }
		/** usage: # iterator.x < iterator1.x #
		    @memo
		*/
	bool operator<(MoveX const & rhs) const 
	 { return current_ < rhs.current_; }
	
	//@}
	/** @name Difference of X Coordinates */
	//@{
		/** usage: # int dx = iterator.x - iterator1.x #
		    @memo
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
      
        /** @name Assignment */
	//@{
		 /** Makes y coordinate identical to rhs, but lets x untouched.
		     Result is undefined if rhs points to different image. \\
		     Usage: # iterator.y = iterator1.y #
		     @memo
		 */
	    MoveY & operator=(MoveY const & rhs)
	    {
		width_ = rhs.width_;
		offset_ = rhs.offset_;
		return *this;
	    }
	//@}
	
        /** @name Increment and Decrement Y */
	//@{
		/** usage: # ++iterator.y #
		    @memo
		*/
        void operator++() {offset_ += width_; }
		/** usage: # iterator.y++ #
		    @memo
		*/
        void operator++(int) {offset_ += width_; }
		/** usage: # --iterator.y #
		    @memo
		*/
        void operator--() {offset_ -= width_; }
		/** usage: # iterator.y-- #
		    @memo
		*/
        void operator--(int) {offset_ -= width_; }
	
	//@}
	
	/** @name Random Offset in Y Direction */
	//@{
		/** usage: # iterator.y += 20 #
		    @memo
		*/
	void operator+=(int dy) {offset_ += dy*width_; }
		/** usage: # iterator.y -= 20 #
		    @memo
		*/
	void operator-=(int dy) {offset_ -= dy*width_; }

	//@}
	
	/** @name Comparision of Y Coordinates */
	//@{
		/** usage: # iterator.y == iterator1.y #
		    @memo
		*/
	bool operator==(MoveY const & rhs) const 
	 { return (offset_ == rhs.offset_); }
		/** usage: # iterator.y != iterator1.y #
		    @memo
		*/
	bool operator!=(MoveY const & rhs) const 
	 { return (offset_ != rhs.offset_); }
		/** usage: # iterator.y < iterator1.y #
		    @memo
		*/
	bool operator<(MoveY const & rhs) const 
	 { return (offset_ < rhs.offset_); }
    
	//@}
	
	/** @name Difference of Y Coordinates */
	//@{
		/** usage: # int dy = iterator.y - iterator1.y #
		    @memo
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
	    /** usage: # iterator == iterator1 #
		@memo
	    */
    bool operator==(ImageIteratorBase const & rhs) const
    {
        return (x == rhs.x) && (y == rhs.y);
    }
    
	    /** usage: # iterator != iterator1 #
		@memo
	    */
    bool operator!=(ImageIteratorBase const & rhs) const
    {
        return (x != rhs.x) || (y != rhs.y);
    }
    
	    /** usage: # Diff2D dist = iterator - iterator1 #
		@memo
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
	    @memo
	*/
    inline ImageIteratorBase & operator+=(Diff2D const & s)
    {
	x += s.x;
	y += s.y;
	return *this;
    }
        /** Subtract offset via Diff2D
	    @memo
	*/
    inline ImageIteratorBase & operator-=(Diff2D const & s)
    {
	x -= s.x;
	y -= s.y;
	return *this;
    }
	
    	    /** access current pixel. \\
	        usage: # SomePixelType value = iterator.current() #
		@memo
	    */
    inline PixelType & current()
    {
	return *(x.current_ + y.offset_);
    }
    
	    /** access pixel at offset (dx, dy) from current pixel. \\
	        usage: # SomePixelType value = iterator.current(dx, dy) #
		@memo
	    */
    inline PixelType & current(int const & dx, int const & dy)
    {
	return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }
    
	    /** access pixel at offset d from current pixel. \\
	        usage: # SomePixelType value = iterator.current(Diff2D(1,1)) #
		@memo
	    */
    inline PixelType & current(Diff2D const & d)
    {
	return *(x.current_ + d.x + y.offset_ + y.width_ * d.y);
    }
    
  	    /** read current pixel. \\
	        usage: # SomePixelType value = iterator.current() #
		@memo
	    */
    inline PixelType const & current() const
    {
	return *(x.current_ + y.offset_);
    }

	    /** read pixel with offset (dx, dy) from current pixel. \\
	        usage: # SomePixelType value = iterator.current(dx, dy) #
		@memo
	    */
    inline PixelType const & current(int const & dx, int const & dy) const
    {
	return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }

	    /** read pixel with offset d from current pixel. \\
	        usage: # SomePixelType value = iterator.current(Diff2D(1,1)) #
		@memo
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

/** Standard 2D random access iterator for images that store the
    data in a plain C array.
    Adds #operator*()#, #operator()(dx, dy)#, #operator[](dy)# and #operator[](Diff2D)#
    to access the current
    pixel and the pixel at distance (dx, dy).
    The usage examples assume that you constructed two iterators like
    this:
    
    \begin{verbatim}
    ImageIteratorBase<SomePixelType> iterator(base, width);
    ImageIteratorBase<SomePixelType> iterator1(base, width);
    \end{verbatim}
    
    See the paper: U. Koethe: 
    \URL[Reusable Algorithms in Image Processing]{documents/GenericProg2D.ps} 
    for a discussion of the concepts behind ImageIterators.

    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
*/    
template <class PIXELTYPE>
class ImageIterator: public ImageIteratorBase<PIXELTYPE>
{
  public:

	/** The underlying image's pixel type.
	    @memo
	*/
    typedef PIXELTYPE PixelType;
    
    /** @name Construction and Assignment */
    //@{
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
    //@}

    /** @name Computed Assignment */
    //@{
        /** Add offset via Diff2D
	    @memo
	*/
    inline ImageIterator & operator+=(Diff2D const & s)
    {
	ImageIteratorBase<PIXELTYPE>::operator+=(s);
	return *this;
    }
        /** Subtract offset via Diff2D
	    @memo
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
	    @memo
	*/
    inline ImageIterator operator+(Diff2D const & s) const
    {
	ImageIterator ret(*this);
	
	ret += s;
	
	return ret;
    }
        /** Subtract a distance
	    @memo
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
	    /** Access current pixel. \\
	        usage: # SomePixelType value = *iterator #
		@memo
	    */
    inline PixelType & operator*()
    {
	return this->current();
    }
    
	    /** Read current pixel. \\
	        usage: # SomePixelType value = *iterator #
		@memo
	    */
    inline PixelType const & operator*() const
    {
	return this->current();
    }
    
	    /** Call member of current pixel. \\
	        usage: # iterator->pixelMemberFunction() #
		@memo
	    */
    inline PixelType * operator->()
    {
	return &(this->current)();
    }
    
	    /** Call const member of current pixel. \\
	        usage: # iterator->constPixelMemberFunction() #
	        Use this operation in inner loop.
		@memo
	    */
    inline PixelType const * operator->() const
    {
	return &(this->current)();
    }
    
	    /** Access pixel at offset from current location. \\
	        usage: # SomePixelType value = iterator[Diff2D(1,1)] #
		@memo
	    */
    inline PixelType & operator[](Diff2D const & d)
    {
	return this->current(d);
    }
    
	    /** Read pixel at offset from current location. \\
	        usage: # SomePixelType value = iterator[Diff2D(1,1)] #
		@memo
	    */
    inline PixelType const & operator[](Diff2D const & d) const
    {
	return this->current(d);
    }
    
	    /** Access pixel at offset (dx, dy) from current location. \\
	        usage: # SomePixelType value = iterator(dx, dy) #
		@memo
	    */
    inline PixelType & operator()(int const & dx, int const & dy)
    {
	return this->current(dx, dy);
    }
    
	    /** Read pixel at offset (dx, dy) from current location. \\
	        usage: # SomePixelType value = iterator(dx, dy) #
		@memo
	    */
    inline PixelType const & operator()(int const & dx, int const & dy) const
    {
	return this->current(dx, dy);
    }
    
	    /** Access pixel with offset [dy][dx] from current pixel.
	        Note that the 'x' index is the trailing index. \\
	        usage: # SomePixelType value = iterator[dy][dx] #
		@memo
	    */
    inline PixelType * operator[](int const & dy)
    {
	return &(this->current())+dy*width();
    }
    
	    /** Read pixel with offset [dy][dx] from current pixel.
	        Note that the 'x' index is the trailing index. \\
	        usage: # SomePixelType value = iterator[dy][dx] #
		@memo
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

/** Standard 2D random access const iterator for images that
    store the data as a plain C array.
    Adds #operator*()#, #operator()(dx, dy)#, and #operator[](dy)#
    to access the current
    pixel and the pixel at distance (dx, dy).
    The usage examples assume that you constructed two iterators like
    this:
    
    #      ConstImageIterator<SomePixelType> iterator(base, width)  #
    #      ConstImageIterator<SomePixelType> iterator1(base, width) #

    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
*/
template <class PIXELTYPE>
class ConstImageIterator: public ImageIteratorBase<PIXELTYPE>
{
  public:

	/** The underlying image's pixel type.
	    @memo
	*/
    typedef PIXELTYPE PixelType;
    
    /** @name Construction and Assignment */
    //@{
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
    //@}

    /** @name Computed Assignment */
    //@{
        /** Add offset via Diff2D
	    @memo
	*/
    inline ConstImageIterator & operator+=(Diff2D const & s)
    {
	ImageIteratorBase<PIXELTYPE>::operator+=(s);
	return *this;
    }
        /** Subtract offset via Diff2D
	    @memo
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
	    @memo
	*/
    inline ConstImageIterator operator+(Diff2D const & s) const
    {
	ConstImageIterator ret(*this);
	
	ret += s;
	
	return ret;
    }
        /** Subtract a distance
	    @memo
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
	    /** Read current pixel. \\
	        usage: # SomePixelType value = *iterator #
		@memo
	    */
    inline PixelType const & operator*() const
    {
	return this->current();
    }
    
	    /** Call member of current pixel. \\
	        usage: # iterator->pixelMemberFunction() #
		@memo
	    */
    inline PixelType const * operator->() const
    {
	return &(this->current)();
    }
    
	    /** Read pixel at offset (dx, dy) from current location. \\
	        usage: # SomePixelType value = iterator(dx, dy) #
		@memo
	    */
    inline PixelType const & operator()(int const & dx, int const & dy) const
    {
	return this->current(dx, dy);
    }

	    /** Read pixel at offset from current location. \\
	        usage: # SomePixelType value = iterator[Diff2D(1,1)] #
		@memo
	    */
    inline PixelType const & operator[](Diff2D const & d) const
    {
	return this->current(d);
    }
    
	    /** Read pixel with offset [dy][dx] from current pixel.
	        Note that the 'x' index is the trailing index. \\
	        usage: # SomePixelType value = iterator[dy][dx] #
		@memo
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

/** Iterator that always returns the constant specified in the 
    constructor. This iterator can be used to simulate an image that
    does not actually exist.
    
    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
*/
template <class PIXELTYPE>
class ConstValueIterator
{
  public:
        /** the iterator's PixelType
	    @memo
	*/
    typedef PIXELTYPE PixelType;
    
	/** Let operations act in X direction
	*/
    typedef int MoveX;

	/** Let operations act in Y direction
	*/
    typedef int MoveY;
    
    /** @name Construction and Assignment */
    //@{
        /** Default Constructor. (the constant is set to
	    #NumericTraits<PIXELTYPE>::zero()# )
	    @memo
	*/
    ConstValueIterator()
    : value_(NumericTraits<PIXELTYPE>::zero()), x(0), y(0)
    {}
    
        /** Construct with given constant.
	    @memo
	*/
    ConstValueIterator(PixelType const & v)
    : value_(v), x(0), y(0)
    {}
    
        /** Copy Constructor.
	    @memo
	*/
    ConstValueIterator(ConstValueIterator const & v)
    : value_(v.value_), x(v.x), y(v.y)
    {}
    
        /** Copy Assigment.
	    @memo
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
	    @memo
	*/
    ConstValueIterator & operator+=(Diff2D const & d)
    {
        x += d.x;
	y += d.y;
        return *this;
    }
    
        /** Move iterator by specified distance.
	    @memo
	*/
    ConstValueIterator & operator-=(Diff2D const & d)
    {
        x -= d.x;
	y -= d.y;
        return *this;
    }
    //@}
    
    /** @name Methods */
    //@{

        /** Create iterator at specified distance.
	    @memo
	*/
    ConstValueIterator operator+(Diff2D const & d) const
    {
        ConstValueIterator ret(*this);
	ret += d;
	return ret;
    }
    
        /** Create iterator at specified distance.
	    @memo
	*/
    ConstValueIterator operator-(Diff2D const & d) const
    {
        ConstValueIterator ret(*this);
	ret -= d;
	return ret;
    }
    
        /** Compute distance between two iterators
	    @memo
	*/
    Diff2D operator-(ConstValueIterator const & r) const
    {
        return Diff2D(x - r.x, y - r.y);
    }
    
        /** Equality.
	    @memo
	*/
    bool operator==(ConstValueIterator const & r) const
    {
        return (x == r.x) && (y == r.y);
    }
    
        /** Inequality.
	    @memo
	*/
    bool operator!=(ConstValueIterator const & r) const
    {
        return (x != r.x) || (y != r.y);
    }
    
        /** Read current pixel (return specified constant).
	    @memo
	*/
    PixelType const & operator*() const
    {
        return value_;
    }
    
        /** Read pixel at a distance (return specified constant).
	    @memo
	*/
    PixelType const & operator()(int const &, int const &) const
    {
        return value_;
    }

        /** Read pixel at a distance (return specified constant).
	    @memo
	*/
    PixelType const & operator[](Diff2D const &) const
    {
        return value_;
    }
    //@}
  
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


#if 0
/********************************************************/
/*                                                      */
/*                 CoordinateIterator                   */
/*                                                      */
/********************************************************/

/* This iterator returns the its current coordinate instead of a 
    value. It behaves like an image where each pixel contains this
    coordinate. This can, for example, be usd to initialize an image with
    a function that depends on the coordinates. Suupose we want to 
    init a chessboard image with squares of 1 pixel by 1 pixel.
    This is realized by the following functor:
    
    \begin{verbatim}
    struct ChessBoard
    {
        int operator()(Diff2D const & d) const
	{
	    return (d.x + d.y) & 1;
	}
    };
    \end{verbatim}
    
    Using this functor and the #CoordinateIterator#, we can apply the 
    \Ref{transformImage}() algorithm to do the initialization:
    
    \begin{verbatim}
    BImage img(w,h);
    
    transformImage(CoordinateIterator(), CoordinateIterator() + Diff2D(w,h),
                   StandardConstAccessor<Diff2D>(),
		   destImage(img), ChessBoard());
    \end{verbatim}
    
    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
*/
class CoordinateIterator
{
  public:
        /* the iterator's PixelType
	    @memo
	*/
    typedef Diff2D PixelType;
    
	/* Let operations act in X direction
	*/
    typedef int MoveX;

	/* Let operations act in Y direction
	*/
    typedef int MoveY;
    
        /* Default Constructor. Init iterator at position (0,0)
	    @memo
	*/
    CoordinateIterator()
    : x(0), y(0)
    {}
    
        /* Construct at given position.
	    @memo
	*/
    CoordinateIterator(int ax, int ay)
    : x(ax), y(ay)
    {}
    
        /* Copy Constructor.
	    @memo
	*/
    CoordinateIterator(CoordinateIterator const & v)
    : x(v.x), y(v.y)
    {}
    
        /* Copy Assigment.
	    @memo
	*/
    CoordinateIterator & operator=(CoordinateIterator const & v)
    {
        if(this != &v)
	{
	    x = v.x;
	    y = v.y;
	}
	return *this;
    }
    
        /* Move iterator by specified distance.
	    @memo
	*/
    CoordinateIterator & operator+=(Diff2D const & d)
    {
        x += d.x;
	y += d.y;
        return *this;
    }
    
        /* Move iterator by specified distance.
	    @memo
	*/
    CoordinateIterator & operator-=(Diff2D const & d)
    {
        x -= d.x;
	y -= d.y;
        return *this;
    }

        /* Create iterator at specified distance.
	    @memo
	*/
    CoordinateIterator operator+(Diff2D const & d) const
    {
        CoordinateIterator ret(*this);
	ret += d;
	return ret;
    }
    
        /* Create iterator at specified distance.
	    @memo
	*/
    CoordinateIterator operator-(Diff2D const & d) const
    {
        CoordinateIterator ret(*this);
	ret -= d;
	return ret;
    }
    
        /* Compute distance between two iterators
	    @memo
	*/
    Diff2D operator-(CoordinateIterator const & r) const
    {
        return Diff2D(x - r.x, y - r.y);
    }
    
        /* Equality.
	    @memo
	*/
    bool operator==(CoordinateIterator const & r) const
    {
        return (x == r.x) && (y == r.y);
    }
    
        /* Inequality.
	    @memo
	*/
    bool operator!=(CoordinateIterator const & r) const
    {
        return (x != r.x) || (y != r.y);
    }
    
        /* Read current coordinate.
	    @memo
	*/
    PixelType operator*() const
    {
        return PixelType(x,y);
    }
    
        /* Read coordinate at a distance.
	    @memo
	*/
    PixelType operator()(int const & dx, int const & dy) const
    {
        return PixelType(x + dx, y + dy);
    }

        /* Read coordinate at a distance.
	    @memo
	*/
    PixelType operator[](Diff2D const & d) const
    {
        return PixelType(x + d.x, y + d.y);
    }
    
    int x;
    int y;
};


template <> 
struct IteratorTraits<CoordinateIterator > 
{
    typedef StandardConstValueAccessor<Diff2D> DefaultAccessor;
};
#endif /* #if 0 */

typedef Diff2D CoordinateIterator;

template <> 
struct IteratorTraits<Diff2D > 
{
    typedef StandardConstValueAccessor<Diff2D> DefaultAccessor;
};


/** @name CoordinateIterator 

    This used to be a separate class, 
    but has now become an alias for \Ref{Diff2D}. This is possible because
    Diff2D now provides all the necessary functionality. 
    
    CoordinateIterator behaves like a read-only \Ref{ImageIterator} for 
    an image in which each pixel contains its coordinate. This is useful for
    algorithms that need access to the current pixel's location.
    For example, you can use CoordinateIterator/Diff2D to 
    find the center of mass of an image region. To implement this,
    we first need a functor for center-of-mass calculations:

    \begin{verbatim}
    
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
    \end{verbatim}
    
    Using this functor, we find the center of mass like so:
    
    \begin{verbatim}
    BImage img(w,h);
    ... // mark a region in the image with '1', background with '0'
    
    CenterOfMassFunctor center;
    
    inspectImageIf(
        srcIterRange(Diff2D(), Diff2D() + img.size()),
        srcImage(img),
        center);
        
    std::cout << "Center of mass: " << center.xCenter() << 
                                ", " << center.yCenter() << std::endl;
    \end{verbatim}    
    
    Include-File:
    \URL[vigra/imageiterator.hxx]{../include/vigra/imageiterator.hxx}
    
    @memo Simulate an image where each pixel contains its coordinate
*/
//@}



#endif // VIGRA_IMAGEITERATOR_HXX
