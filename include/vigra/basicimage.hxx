/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.1.4a, Nov 20 2001 )                                   */
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
 
 
#ifndef VIGRA_BASICIMAGE_HXX
#define VIGRA_BASICIMAGE_HXX


#include <new>
#include <memory>
#include "vigra/utilities.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                  BasicImageIteratorBase              */
/*                                                      */
/********************************************************/

template <class PIXELTYPE, class ITERATOR>
class BasicImageIteratorBase
{
  protected:
    typedef ITERATOR Iterator;

  public:
    typedef PIXELTYPE PixelType;
    typedef PIXELTYPE value_type;
    
    typedef int MoveX;
    
    typedef Iterator MoveY;
    
    bool operator==(BasicImageIteratorBase const & rhs) const
    {
        return (x == rhs.x) && (y == rhs.y);
    }

    bool operator!=(BasicImageIteratorBase const & rhs) const
    {
        return (x != rhs.x) || (y != rhs.y);
    }
    
    Diff2D operator-(BasicImageIteratorBase const & rhs) const
    {
        return Diff2D(x - rhs.x, y - rhs.y);
    }
    
    MoveX x;
    MoveY y;

  protected:
    BasicImageIteratorBase(Iterator line)
    : x(0),
      y(line)
    {}
    
    BasicImageIteratorBase(BasicImageIteratorBase const & rhs)
    : x(rhs.x),
      y(rhs.y)
    {}
    
    BasicImageIteratorBase()
    : x(0),
      y(0)
    {}
    
    BasicImageIteratorBase & operator=(BasicImageIteratorBase const & rhs)
    {
        if(this != &rhs)
        {
            x = rhs.x;
            y = rhs.y;
        }
        return *this;
    }
    
    inline BasicImageIteratorBase & operator+=(Diff2D const & s)
    {
        x += s.x;
        y += s.y;
        return *this;
    }

    inline BasicImageIteratorBase & operator-=(Diff2D const & s)
    {
        x -= s.x;
        y -= s.y;
        return *this;
    }
    
    inline PixelType & current()
    {
        return *(*y + x );
    }
    
    inline PixelType & current(Diff2D const & dist)
    {
        return *(*(y + dist.y) + x + dist.x);
    }
    
    inline PixelType & current(int const & dx, int const & dy)
    {
        return *(*(y + dy) + x + dx);
    }
    
    inline PixelType const & current() const
    {
        return *(*y + x );
    }

    inline PixelType const & current(Diff2D const & dist) const
    {
        return *(*(y + dist.y) + x + dist.x);
    }

    inline PixelType const & current(int const & dx, int const & dy) const
    {
        return *(*(y + dy) + x + dx);
    }
    
    inline Iterator line() const
    {
        return y;
    }
}; 

/********************************************************/
/*                                                      */
/*                    BasicImageIterator                */
/*                                                      */
/********************************************************/

/** Implementation of the standard image iterator for \ref vigra::BasicImage.
    See \ref vigra::ImageIterator for documentation.

    <b>\#include</b> "<a href="basicimage_8hxx-source.html">vigra/basicimage.hxx</a>"    
    Namespace: vigra
*/    
template <class PIXELTYPE, class ITERATOR>
class BasicImageIterator
: public BasicImageIteratorBase<PIXELTYPE, ITERATOR>
{
  public:

    typedef PIXELTYPE PixelType;
    
    inline BasicImageIterator(ITERATOR line)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(line)
    {}
    
    inline BasicImageIterator()
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>()
    {}  
    
    inline BasicImageIterator(const BasicImageIterator & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    inline BasicImageIterator & operator=(const BasicImageIterator & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    inline BasicImageIterator & operator+=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator+=(s);
        return *this;
    }

    inline BasicImageIterator & operator-=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator-=(s);
        return *this;
    }
    
    Diff2D operator-(BasicImageIterator const & rhs) const
    {
        return Diff2D(this->x - rhs.x, this->y - rhs.y);
    }

    inline BasicImageIterator operator+(Diff2D const & s) const
    {
        BasicImageIterator ret(*this);
        
        ret += s;
        
        return ret;
    }

    inline BasicImageIterator operator-(Diff2D const & s) const
    {
        BasicImageIterator ret(*this);
        
        ret -= s;
        
        return ret;
    }

    inline PixelType & operator*()
    {
        return this->current();
    }
    
    inline PixelType const & operator*() const
    {
        return this->current();
    }
    
    inline PixelType * operator->()
    {
        return &(this->current)();
    }
    
    inline PixelType const * operator->() const
    {
        return &(this->current)();
    }
    
    inline PixelType & operator[](Diff2D const & d)
    {
        return this->current(d);
    }
    
    inline PixelType const & operator[](Diff2D const & d) const
    {
        return this->current(d);
    }
    
    inline PixelType * operator[](int const & dy) 
    {
        return this->y[dy] + this->x;
    }
    
    inline PixelType const * operator[](int const & dy) const
    {
        return this->y[dy] + this->x;
    }
    
    inline PixelType & operator()(int const & dx, int const & dy)
    {
        return this->current(dx, dy);
    }
    
    inline PixelType const & operator()(int const & dx, int const & dy) const
    {
        return this->current(dx, dy);
    }
};


/********************************************************/
/*                                                      */
/*                ConstBasicImageIterator               */
/*                                                      */
/********************************************************/

/** Implementation of the standard const image iterator for \ref vigra::BasicImage.
    See \ref vigra::ConstImageIterator for documentation.

    <b>\#include</b> "<a href="basicimage_8hxx-source.html">vigra/basicimage.hxx</a>"    
    Namespace: vigra
*/    
template <class PIXELTYPE, class ITERATOR>
class ConstBasicImageIterator: public BasicImageIteratorBase<PIXELTYPE, ITERATOR>
{
  public:

    typedef PIXELTYPE const PixelType;
    
    inline ConstBasicImageIterator(ITERATOR line)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(line)
    {}
    
    inline ConstBasicImageIterator()
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>()
    {}  
    
    inline ConstBasicImageIterator(const ConstBasicImageIterator & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    inline ConstBasicImageIterator(const BasicImageIteratorBase<PIXELTYPE, ITERATOR> & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    inline ConstBasicImageIterator & operator=(const ConstBasicImageIterator & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    inline ConstBasicImageIterator & operator=(const BasicImageIteratorBase<PIXELTYPE, ITERATOR> & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    inline ConstBasicImageIterator & operator+=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator+=(s);
        return *this;
    }

    inline ConstBasicImageIterator & operator-=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator-=(s);
        return *this;
        }
    
    Diff2D operator-(ConstBasicImageIterator const & rhs) const
    {
        return Diff2D(this->x - rhs.x, this->y - rhs.y);
    }

    inline ConstBasicImageIterator operator+(Diff2D const & s) const
    {
        ConstBasicImageIterator ret(*this);
        
        ret += s;
        
        return ret;
    }

    inline ConstBasicImageIterator operator-(Diff2D const & s) const
    {
        ConstBasicImageIterator ret(*this);
        
        ret -= s;
        
        return ret;
    }

    inline PixelType & operator*() const
    {
        return this->current();
    }
    
    inline PixelType * operator->() const
    {
        return &(this->current)();
    }
    
    inline PixelType & operator[](Diff2D const & d) const
    {
        return this->current(d.x, d.y);
    }
    
    inline PixelType * operator[](int const & dy) const
    {
        return this->y[dy] + this->x;
    }
    
    inline PixelType & operator()(int const & dx, int const & dy) const
    {
        return this->current(dx, dy);
    }
};

template <class T> struct IteratorTraits;

/********************************************************/
/*                                                      */
/*                      RowIterator                     */
/*                                                      */
/********************************************************/

/* Specialization for higher performance.
*/
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION
template <class T>
class RowIterator;

template <class ValueType>
class RowIterator<BasicImageIterator<ValueType, ValueType **> >
{
    typedef RowIterator Self;
    typedef BasicImageIterator<ValueType, ValueType **> Adaptee;

    ValueType * data_;
    
  public:

    typedef ValueType PixelType;
    
    RowIterator(Adaptee const & i)
    : data_(&const_cast<ValueType &>(*i))
    {}
    
    RowIterator(ValueType * d)
    : data_(d)
    {}
    
#else

template <class IMAGE_ITERATOR>
class BasicImageRowIterator
{
    typedef BasicImageRowIterator Self;
    typedef typename IMAGE_ITERATOR::PixelType ValueType;
    typedef IMAGE_ITERATOR Adaptee;

    ValueType * data_;
    
  public:

    typedef ValueType PixelType;
    
    BasicImageRowIterator(IMAGE_ITERATOR const & i)
    : data_(&const_cast<ValueType &>(*i))
    {}
    
    BasicImageRowIterator(ValueType * d)
    : data_(d)
    {}
    
#endif

    Self & operator=(Adaptee const & i)
    {
        data_ = &const_cast<ValueType &>(*i);
        return *this;
    }
    
    Self & operator++()
    {
        ++data_;
        return *this;
    }

    Self operator++(int)
    {
        Self ret(*this);
        data_++;
        return ret;
    }
    
    Self & operator--()
    {
        --data_;
        return *this;
    }
    
    Self operator--(int)
    {
        Self ret(*this);
        data_--;
        return ret;
    }
    
    Self & operator+=(int d)
    {
        data_ += d;
        return *this;
    }
    
    Self & operator-=(int d)
    {
        data_ -= d;
        return *this;
    }

    Self operator+(int d) const
    {
        return Self(data_ + d);
    }

    Self operator-(int d) const
    {
        return Self(data_ - d);
    }

    int operator-(Self const & c) const 
    {
        return data_ - c.data_;
    }
    
    bool operator==(Self const & c) const
    {
        return data_ == c.data_;
    }
    
    bool operator!=(Self const & c) const
    {
        return data_ != c.data_;
    }
    
    bool operator<(Self const & c) const
    {
        return data_ < c.data_;
    }
    
    PixelType & operator*()
    {
        return *data_; 
    }
    
    PixelType const & operator*() const
    {
        return *data_; 
    }
    
    PixelType * operator->()
    {
        return data_; 
    }
    
    PixelType const * operator->() const
    {
        return data_; 
    }
    
    PixelType & operator[](int d)
    {
        return data_[d];
    }

    PixelType const & operator[](int d) const
    {
        return data_[d];
    }
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION
template <class T>
class ConstRowIterator;

template <class ValueType>
class ConstRowIterator<ConstBasicImageIterator<ValueType, ValueType **> >
{
    typedef ConstRowIterator Self;
    typedef ConstBasicImageIterator<ValueType, ValueType **> Adaptee;

    ValueType const * data_;
    
  public:

    typedef ValueType const PixelType;
    
    ConstRowIterator(Adaptee const & i)
    : data_(&(*i))
    {}
    
    ConstRowIterator(ValueType const * d)
    : data_(d)
    {}
    
#else

template <class IMAGE_ITERATOR>
class ConstBasicImageRowIterator
{
    typedef ConstBasicImageRowIterator Self;
    typedef typename IMAGE_ITERATOR::PixelType ValueType;
    typedef IMAGE_ITERATOR Adaptee;

    ValueType const * data_;
    
  public:

    typedef ValueType const PixelType;
    
    ConstBasicImageRowIterator(IMAGE_ITERATOR  const & i)
    : data_(&(*i))
    {}
    
    ConstBasicImageRowIterator(ValueType * d)
    : data_(d)
    {}
    
#endif

    Self & operator=(Adaptee const & i)
    {
        data_ = &(*i);
        return *this;
    }
    
    Self & operator++()
    {
        ++data_;
        return *this;
    }

    Self operator++(int)
    {
        Self ret(*this);
        data_++;
        return ret;
    }
    
    Self & operator--()
    {
        --data_;
        return *this;
    }
    
    Self operator--(int)
    {
        Self ret(*this);
        data_--;
        return ret;
    }
    
    Self & operator+=(int d)
    {
        data_ += d;
        return *this;
    }
    
    Self & operator-=(int d)
    {
        data_ -= d;
        return *this;
    }

    Self operator+(int d) const
    {
        return Self(data_ + d);
    }

    Self operator-(int d) const
    {
        return Self(data_ - d);
    }

    int operator-(Self const & c) const 
    {
        return data_ - c.data_;
    }
    
    bool operator==(Self const & c) const
    {
        return data_ == c.data_;
    }
    
    bool operator!=(Self const & c) const
    {
        return data_ != c.data_;
    }
    
    bool operator<(Self const & c) const
    {
        return data_ < c.data_;
    }
    
    PixelType & operator*() const
    {
        return *data_; 
    }
    
    PixelType * operator->() const
    {
        return data_; 
    }
    
    PixelType & operator[](int d) const
    {
        return data_[d];
    }
};

/********************************************************/
/*                                                      */
/*                     BasicImage                       */
/*                                                      */
/********************************************************/

/** \brief Fundamental class template for images .

    <b>\#include</b> "<a href="basicimage_8hxx-source.html">vigra/basicimage.hxx</a>"
    
    Namespace: vigra
*/    
template <class PIXELTYPE>
class BasicImage
{
  public:
    
        /** the BasicImage's pixel type
        */
    typedef PIXELTYPE PixelType;
    
        /** the BasicImage's pixel type
        */
    typedef PIXELTYPE value_type;
    
        /** the BasicImage's 1D random access iterator
        */
    typedef PIXELTYPE * ScanOrderIterator;
    
        /** the BasicImage's 1D random access const iterator
        */
    typedef PIXELTYPE const * ConstScanOrderIterator;
    
        /** the BasicImage's 2D random access iterator
        */
    typedef BasicImageIterator<PixelType, PIXELTYPE **> Iterator;
    
        /** the BasicImage's 2D random access const iterator
        */
    typedef ConstBasicImageIterator<PixelType, PIXELTYPE **> ConstIterator; 

        /** the BasicImage's default accessor
        */
    typedef typename
          IteratorTraits<Iterator>::DefaultAccessor Accessor;

        /** the BasicImage's default const accessor
        */
    typedef typename 
          IteratorTraits<ConstIterator>::DefaultAccessor ConstAccessor;
    
    struct Allocator
    {
        static PixelType * allocate(int n) { 
                  return (PixelType *)::operator new(n*sizeof(PixelType)); }
        static void deallocate(PixelType * p) {
                 ::operator delete(p); }
    };
    
        /** construct image of size 0x0 
        */
    BasicImage()
    : data_(0),
      width_(0),
      height_(0)      
    {}
    
        /** construct image of size width x height 
        */
    BasicImage(int width, int height)
    : data_(0),
      width_(0),
      height_(0)      
    {
        vigra_precondition((width >= 0) && (height >= 0),
             "BasicImage::BasicImage(int width, int height): "
             "width and height must be >= 0.\n");
        
        resize(width, height, PixelType());
    }
    
        /** construct image of size size.x x size.y 
        */
    BasicImage(Diff2D size)
    : data_(0),
      width_(0),
      height_(0)      
    {
        vigra_precondition((size.x >= 0) && (size.y >= 0),
             "BasicImage::BasicImage(Diff2D size): "
             "size.x and size.y must be >= 0.\n");
        
        resize(size.x, size.y, PixelType());
    }
    
        /** construct image of size width*height and initialize every
        pixel with given data (use this constructor, if 
        PixelType doesn't have a default constructor)
        */
    BasicImage(int width, int height, PixelType d)
    : data_(0),
      width_(0),
      height_(0)      
    {
        vigra_precondition((width >= 0) && (height >= 0),
             "BasicImage::BasicImage(int width, int height, PixelType): "
             "width and height must be >= 0.\n");
        
        resize(width, height, d);
    }
    
        /** copy rhs image 
        */
    BasicImage(const BasicImage & rhs)
    : data_(0),
      width_(0),
      height_(0)      
    {
        resizeCopy(rhs);
    }
    
        /** destructor 
        */
    ~BasicImage() 
    { 
        deallocate(); 
    }

        /** copy rhs image (image is resized if necessary)
        */
    BasicImage & operator=(const BasicImage & rhs)
    {
        if(this != &rhs)
        {
            if((width() != rhs.width()) || 
               (height() != rhs.height()))
            {
                resizeCopy(rhs);
            }
            else
            {
                ConstScanOrderIterator is = rhs.begin();
                ConstScanOrderIterator iend = rhs.end();
                ScanOrderIterator id = begin();

                for(; is != iend; ++is, ++id) *id = *is;
            }
        }
        return *this;
    }
     
        /** set Image with const value 
        */
    BasicImage & operator=(PixelType pixel)
    {
        ScanOrderIterator i = begin();
        ScanOrderIterator iend = end();
        
        for(; i != iend; ++i) *i = pixel;
        
        return *this;
    }
    
        /** reset image to specified size (dimensions must not be negative)
        (old data are destroyed) 
        */
    void resize(int width, int height)
    {
        resize(width, height, PixelType());
    }

        /** reset image to specified size (dimensions must not be negative)
        (old data are destroyed) 
        */
    void resize(Diff2D size)
    {
        resize(size.x, size.y, PixelType());
    }

        /** reset image to specified size and initialize it with
            given data (use this if PixelType doesn't have a default
            constructor, dimensions must not be negative, old data are destroyed) 
        */
    void resize(int width, int height, PixelType d)
    {
        vigra_precondition((width >= 0) && (height >= 0),
             "BasicImage::resize(int width, int height, PixelType): "
             "width and height must be >= 0.\n");
        
        PixelType * newdata = 0;
        PixelType ** newlines = 0;
        if(width*height > 0)
        {
            newdata = Allocator::allocate(width*height);

            std::uninitialized_fill_n(newdata, width*height, d);

            newlines = initLineStartArray(newdata, width, height);
        }
        
        deallocate();
        data_ = newdata;
        lines_ = newlines;
        width_ = width;
        height_ = height;
    }


        /** resize image to size of other image and copy it's data 
        */
    void resizeCopy(const BasicImage & rhs)
    {
        PixelType * newdata = 0;  
        PixelType ** newlines = 0;
        if(rhs.width()*rhs.height() > 0)
        {
            newdata = Allocator::allocate(rhs.width()*rhs.height());

            std::uninitialized_copy(rhs.begin(), rhs.end(), newdata);

            newlines = 
               initLineStartArray(newdata, rhs.width(), rhs.height());
        }

        deallocate();
        data_ = newdata;
        lines_ = newlines;
        width_ = rhs.width();
        height_ = rhs.height(); 
    }

        /** width of Image 
        */
    int width() const
    { 
        return width_; 
    }
    
        /** height of Image 
        */
    int height() const
    { 
        return height_; 
    }
        
        /** size of Image 
        */
    inline Diff2D size() const
    { 
        return Diff2D(width(), height()); 
    }
    
        /** test whether a given coordinate is inside the image
        */
    inline bool isInside(Diff2D const & d) const
    {
        return d.x >= 0 && d.y >= 0 &&
               d.x < width() && d.y < height();
    }
        
        /** access pixel at given location. <br>
	    usage: <TT> PixelType value = image[Diff2D(1,2)] </TT>
        */
    inline PixelType & operator[](Diff2D const & d)
    { 
        return lines_[d.y][d.x]; 
    }
        
        /** read pixel at given location. <br>
	    usage: <TT> PixelType value = image[Diff2D(1,2)] </TT>
        */
    inline PixelType const & operator[](Diff2D const & d) const
    { 
        return lines_[d.y][d.x]; 
    }
        
        /** access pixel at given location. <br>
	    usage: <TT> PixelType value = image(1,2) </TT>
        */
    inline PixelType & operator()(int const & dx, int const & dy)
    { 
        return lines_[dy][dx]; 
    }
        
        /** read pixel at given location. <br>
	    usage: <TT> PixelType value = image(1,2) </TT>
        */
    inline PixelType const & operator()(int const & dx, int const & dy) const
    { 
        return lines_[dy][dx]; 
    }
        
        /** access pixel at given location. 
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> PixelType value = image[2][1] </TT>
        */
    inline PixelType * operator[](int const & dy)
    { 
        return lines_[dy]; 
    }
        
        /** read pixel at given location. 
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> PixelType value = image[2][1] </TT>
        */
    inline PixelType const * operator[](int const & dy) const
    { 
        return lines_[dy]; 
    }
        
        /** init 2D random access iterator poining to upper left pixel
        */
    Iterator upperLeft()
    {
        return Iterator(lines_);
    }
    
        /** init 2D random access iterator poining to 
         pixel(width, height), i.e. one pixel right and below lower right 
         corner of the image as is common in C/C++.
        */
    Iterator lowerRight()
    {
        return upperLeft() + size();
    }   
    
        /** init 2D random access const iterator poining to upper left pixel
        */
    ConstIterator upperLeft() const
    {
        return ConstIterator(lines_);
    }
    
        /** init 2D random access const iterator poining to 
         pixel(width, height), i.e. one pixel right and below lower right 
         corner of the image as is common in C/C++.
        */
    ConstIterator lowerRight() const
    {
        return upperLeft() + size();
    }   

        /** init 1D random access iterator pointing to first pixel
        */
    ScanOrderIterator begin()
    {
        return data_;
    }
    
        /** init 1D random access iterator pointing past the end
        */
    ScanOrderIterator end()
    {
        return data_ + width() * height();
    }
    
        /** init 1D random access const iterator pointing to first pixel
        */
    ConstScanOrderIterator begin() const
    {
        return data_;
    }
    
        /** init 1D random access const iterator pointing past the end
        */
    ConstScanOrderIterator end() const
    {
        return data_ + width() * height();
    }
    
        /** return default accessor
        */
    Accessor accessor() 
    {
        return Accessor();
    }

        /** return default const accessor
        */
    ConstAccessor accessor() const
    {
        return ConstAccessor();
    }
    
  private:
  
    void deallocate()
    {
        if(data_) 
        {
            ScanOrderIterator i = begin();
            ScanOrderIterator iend = end();

            for(; i != iend; ++i)   (*i).~PIXELTYPE();

            Allocator::deallocate(data_);
            delete[] lines_;
        }
    }
    
    static PixelType ** initLineStartArray(PixelType * data, int width, int height)
    {
        PixelType ** lines = new PIXELTYPE*[height];
        for(int y=0; y<height; ++y) 
             lines[y] = data + y*width;
        return lines;
    }
      
    PIXELTYPE * data_;
    PIXELTYPE ** lines_;
    int width_, height_;
};

} // namespace vigra

#endif // VIGRA_BASICIMAGE_HXX
