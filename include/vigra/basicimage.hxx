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

template <class IMAGEITERATOR>
class LineBasedColumnIteratorPolicy
{
  public:
    typedef IMAGEITERATOR                             ImageIterator;
    typedef typename IMAGEITERATOR::LineStartIterator LineStartIterator;
    typedef typename ImageIterator::value_type        value_type;
    typedef typename ImageIterator::difference_type::MoveY            
                                                      difference_type;
    typedef typename ImageIterator::reference         reference;
    typedef typename ImageIterator::index_reference   index_reference;
    typedef typename ImageIterator::pointer           pointer;
    typedef std::random_access_iterator_tag           iterator_category;
    
    
    struct BaseType
    {
        explicit BaseType(LineStartIterator c = LineStartIterator(), 
                          difference_type o = 0)
        : line_start_(c), offset_(o)
        {}
        
        LineStartIterator line_start_;
        difference_type offset_;
    };
    
    static void initialize(BaseType & d) {}
    
    static reference dereference(BaseType const & d)
        { return const_cast<reference>(*(*d.line_start_ + d.offset_)); }
    
    static index_reference dereference(BaseType const & d, difference_type n)
    { 
        return const_cast<index_reference>(*(d.line_start_[n] + d.offset_));
    }
    
    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.line_start_ == d2.line_start_; }
    
    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.line_start_ < d2.line_start_; }
    
    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return d1.line_start_ - d2.line_start_; }
    
    static void increment(BaseType & d)
        { ++d.line_start_; }
    
    static void decrement(BaseType & d)
        { --d.line_start_; }
    
    static void advance(BaseType & d, difference_type n)
        { d.line_start_ += n; }
};

/********************************************************/
/*                                                      */
/*                  BasicImageIteratorBase              */
/*                                                      */
/********************************************************/

template <class PIXELTYPE, class ITERATOR>
class BasicImageIteratorBase
{
  public:
    typedef ITERATOR  LineStartIterator;
    typedef PIXELTYPE PixelType;
    typedef PIXELTYPE value_type;
    
    typedef int MoveX;
    
    typedef LineStartIterator MoveY;
    
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
    BasicImageIteratorBase(LineStartIterator line)
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
    
    BasicImageIteratorBase & operator+=(Diff2D const & s)
    {
        x += s.x;
        y += s.y;
        return *this;
    }

    BasicImageIteratorBase & operator-=(Diff2D const & s)
    {
        x -= s.x;
        y -= s.y;
        return *this;
    }
    
    PixelType & current()
    {
        return *(*y + x );
    }
    
    PixelType & current(Diff2D const & dist)
    {
        return *(*(y + dist.y) + x + dist.x);
    }
    
    PixelType & current(int const & dx, int const & dy)
    {
        return *(*(y + dy) + x + dx);
    }
    
    PixelType const & current() const
    {
        return *(*y + x );
    }

    PixelType const & current(Diff2D const & dist) const
    {
        return *(*(y + dist.y) + x + dist.x);
    }

    PixelType const & current(int const & dx, int const & dy) const
    {
        return *(*(y + dy) + x + dx);
    }
    
    LineStartIterator line() const
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
    typedef PIXELTYPE value_type;
    
    typedef Diff2D               difference_type;
    typedef PIXELTYPE &          reference;
    typedef PIXELTYPE &          index_reference;
    typedef PIXELTYPE *          pointer;
    typedef image_traverser_tag  iterator_category;
    typedef PIXELTYPE *          row_iterator;
    typedef IteratorAdaptor<LineBasedColumnIteratorPolicy<BasicImageIterator> > 
                                 column_iterator;
    
    BasicImageIterator(ITERATOR line)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(line)
    {}
    
    BasicImageIterator()
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>()
    {}  
    
    BasicImageIterator(const BasicImageIterator & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    BasicImageIterator & operator=(const BasicImageIterator & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    BasicImageIterator & operator+=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator+=(s);
        return *this;
    }

    BasicImageIterator & operator-=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator-=(s);
        return *this;
    }
    
    Diff2D operator-(BasicImageIterator const & rhs) const
    {
        return Diff2D(this->x - rhs.x, this->y - rhs.y);
    }

    BasicImageIterator operator+(Diff2D const & s) const
    {
        BasicImageIterator ret(*this);
        
        ret += s;
        
        return ret;
    }

    BasicImageIterator operator-(Diff2D const & s) const
    {
        BasicImageIterator ret(*this);
        
        ret -= s;
        
        return ret;
    }

    reference operator*() const
    {
        return const_cast<reference>(current());
    }
    
    pointer operator->() const
    {
        return const_cast<pointer>(&current());
    }
    
    index_reference operator[](Diff2D const & d) const
    {
        return const_cast<index_reference>(current(d));
    }
    
    index_reference operator()(int dx, int dy) const
    {
        return const_cast<index_reference>(current(dx, dy));
    }
    
    pointer operator[](int dy) const
    {
        return const_cast<pointer>(y[dy] + x);
    }
    
    row_iterator rowIterator() const
        { return const_cast<row_iterator>(&current()); }
    
    column_iterator columnIterator() const
    { 
        return 
           column_iterator(typename column_iterator::BaseType(line(), x));
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

    typedef PIXELTYPE            PixelType;
    typedef PIXELTYPE            value_type;
    typedef Diff2D               difference_type;
    typedef PIXELTYPE const &    reference;
    typedef PIXELTYPE const &    index_reference;
    typedef PIXELTYPE const *    pointer;
    typedef image_traverser_tag  iterator_category;
    typedef PIXELTYPE const *    row_iterator;
    typedef IteratorAdaptor<LineBasedColumnIteratorPolicy<ConstBasicImageIterator> > 
                                 column_iterator;
    
    ConstBasicImageIterator(ITERATOR line)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(line)
    {}
    
    ConstBasicImageIterator()
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>()
    {}  
    
    ConstBasicImageIterator(const ConstBasicImageIterator & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    ConstBasicImageIterator(const BasicImageIteratorBase<PIXELTYPE, ITERATOR> & rhs)
    : BasicImageIteratorBase<PIXELTYPE, ITERATOR>(rhs)
    {}  
    
    ConstBasicImageIterator & operator=(const ConstBasicImageIterator & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    ConstBasicImageIterator & operator=(const BasicImageIteratorBase<PIXELTYPE, ITERATOR> & rhs)
    {
        if(this != &rhs)
        {
            BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator=(rhs);
        }
        return *this;
    }

    ConstBasicImageIterator & operator+=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator+=(s);
        return *this;
    }

    ConstBasicImageIterator & operator-=(Diff2D const & s)
    {
        BasicImageIteratorBase<PIXELTYPE, ITERATOR>::operator-=(s);
        return *this;
        }
    
    Diff2D operator-(ConstBasicImageIterator const & rhs) const
    {
        return Diff2D(this->x - rhs.x, this->y - rhs.y);
    }

    ConstBasicImageIterator operator+(Diff2D const & s) const
    {
        ConstBasicImageIterator ret(*this);
        
        ret += s;
        
        return ret;
    }

    ConstBasicImageIterator operator-(Diff2D const & s) const
    {
        ConstBasicImageIterator ret(*this);
        
        ret -= s;
        
        return ret;
    }

    reference operator*() const
    {
        return current();
    }
    
    pointer operator->() const
    {
        return &current();
    }
    
    index_reference operator[](Diff2D const & d) const
    {
        return current(d.x, d.y);
    }
    
    index_reference operator()(int dx, int dy) const
    {
        return current(dx, dy);
    }
    
    pointer operator[](int dy) const
    {
        return y[dy] + x;
    }
    
    row_iterator rowIterator() const
        { return &current(); }
    
    column_iterator columnIterator() const
    { 
        return 
           column_iterator(typename column_iterator::BaseType(line(), x));
    }
};

template <class T> struct IteratorTraits;

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
    typedef PIXELTYPE value_type;
    
        /** the BasicImage's pixel type
        */
    typedef PIXELTYPE PixelType;

        /** the BasicImage's reference type (i.e. the
            return type of image[diff] and image(dx,dy))
        */
    typedef PIXELTYPE &       reference;

        /** the BasicImage's const reference type (i.e. the
            return type of image[diff] and image(dx,dy) when image is const)
        */
    typedef PIXELTYPE const & const_reference;

        /** the BasicImage's pointer type 
        */
    typedef PIXELTYPE *       pointer;

        /** the BasicImage's const pointer type 
        */
    typedef PIXELTYPE const * const_pointer;

        /** the BasicImage's 1D random access iterator
            (note: lower case 'iterator' is a STL compatible 1D random 
             access iterator, don't confuse with capitalized Iterator)
        */
    typedef PIXELTYPE * iterator;
 
        /** the BasicImage's 1D random access iterator
        */
   typedef PIXELTYPE * ScanOrderIterator;

        /** the BasicImage's 1D random access const iterator
            (note: lower case 'const_iterator' is a STL compatible 1D 
            random access const iterator)
        */
    typedef PIXELTYPE const * const_iterator;

        /** the BasicImage's 1D random access const iterator
        */
    typedef PIXELTYPE const * ConstScanOrderIterator;

        /** the BasicImage's 2D random access iterator ('traverser')
        */
    typedef BasicImageIterator<PixelType, PIXELTYPE **> traverser;

        /** the BasicImage's 2D random access iterator
            (note: capitalized 'Iterator' is a 2D image iterator,
             don't confuse with lower case iterator)
        */
    typedef BasicImageIterator<PixelType, PIXELTYPE **> Iterator;
 
        /** the BasicImage's 2D random access const iterator ('conts traverser')
        */
    typedef ConstBasicImageIterator<PixelType, PIXELTYPE **> const_traverser;
 
        /** the BasicImage's 2D random access const iterator
        */
    typedef ConstBasicImageIterator<PixelType, PIXELTYPE **> ConstIterator; 

        /** the BasicImage's difference type (argument type of image[diff])
        */
    typedef Diff2D difference_type;

         /** the BasicImage's size type (result type of image.size())
        */
    typedef Diff2D size_type;

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
     
        /** \deprecated set Image with const value 
        */
    BasicImage & operator=(PixelType pixel)
    {
        ScanOrderIterator i = begin();
        ScanOrderIterator iend = end();
        
        for(; i != iend; ++i) *i = pixel;
        
        return *this;
    }
    
        /** set Image with const value 
        */
    BasicImage & init(PixelType const & pixel)
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
    Diff2D size() const
    { 
        return Diff2D(width(), height()); 
    }
    
        /** test whether a given coordinate is inside the image
        */
    bool isInside(Diff2D const & d) const
    {
        return d.x >= 0 && d.y >= 0 &&
               d.x < width() && d.y < height();
    }
        
        /** access pixel at given location. <br>
	    usage: <TT> PixelType value = image[Diff2D(1,2)] </TT>
        */
    PixelType & operator[](Diff2D const & d)
    { 
        return lines_[d.y][d.x]; 
    }
        
        /** read pixel at given location. <br>
	    usage: <TT> PixelType value = image[Diff2D(1,2)] </TT>
        */
    PixelType const & operator[](Diff2D const & d) const
    { 
        return lines_[d.y][d.x]; 
    }
        
        /** access pixel at given location. <br>
	    usage: <TT> PixelType value = image(1,2) </TT>
        */
    PixelType & operator()(int dx, int dy)
    { 
        return lines_[dy][dx]; 
    }
        
        /** read pixel at given location. <br>
	    usage: <TT> PixelType value = image(1,2) </TT>
        */
    PixelType const & operator()(int dx, int dy) const
    { 
        return lines_[dy][dx]; 
    }
        
        /** access pixel at given location. 
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> PixelType value = image[2][1] </TT>
        */
    PixelType * operator[](int dy)
    { 
        return lines_[dy]; 
    }
        
        /** read pixel at given location. 
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> PixelType value = image[2][1] </TT>
        */
    PixelType const * operator[](int dy) const
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
