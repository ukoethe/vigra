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
 
 
#ifndef VIGRA_IMAGEITERATOR_HXX
#define VIGRA_IMAGEITERATOR_HXX

#include "vigra/utilities.hxx"
#include "vigra/accessor.hxx"
#include "vigra/iteratortraits.hxx"

namespace vigra { 

template <class IMAGEITERATOR>
class ContigousMemoryColumnIteratorPolicy
{
  public:
    typedef IMAGEITERATOR                            ImageIterator;
    typedef typename IMAGEITERATOR::value_type       value_type;
    typedef typename IMAGEITERATOR::difference_type::MoveY            
                                                     difference_type;
    typedef typename IMAGEITERATOR::reference        reference;
    typedef typename IMAGEITERATOR::index_reference  index_reference;
    typedef typename IMAGEITERATOR::pointer          pointer;
    typedef std::random_access_iterator_tag iterator_category;
    
    
    struct BaseType
    {
        explicit BaseType(pointer c = 0, difference_type o = 0)
        : current_(c), offset_(o)
        {}
        
        pointer current_;
        difference_type offset_;
    };
    
    static void initialize(BaseType & d) {}
    
    static reference dereference(BaseType const & d)
        { return const_cast<reference>(*d.current_); }
    
    static index_reference dereference(BaseType const & d, difference_type n)
    { 
        return const_cast<index_reference>(d.current_[n*d.offset_]);
    }
    
    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.current_ == d2.current_; }
    
    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.current_ < d2.current_; }
    
    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return (d1.current_ - d2.current_) / d1.offset_; }
    
    static void increment(BaseType & d)
        { d.current_ += d.offset_; }
    
    static void decrement(BaseType & d)
        { d.current_ -= d.offset_; }
    
    static void advance(BaseType & d, difference_type n)
        { d.current_ += d.offset_*n; }
};

/** \addtogroup ImageIterators  Image Iterators

    \brief General image iterator definition and implementations.
    
<p>
    The following tables describe the general requirements for image iterators
    and their iterator traits. The iterator implementations provided here
    may be used for any image data type that stores its
    data as a linear array of pixels. The array will be interpreted as a
    row-major matrix with a particular width. 
</p>
<h3>Requirements for Image Iterators</h3>
<p>

<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr>
    <th colspan=2>Local Types</th><th>Meaning</th>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::value_type</tt></td><td>the underlying image's pixel type</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::PixelType</tt></td><td>the underlying image's pixel type</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::reference</tt></td>
    <td>the iterator's reference type (return type of <TT>*iter</TT>). Will be 
    <tt>value_type &</tt> for a mutable iterator, and convertible to 
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::index_reference</tt></td>
    <td>the iterator's index reference type (return type of <TT>iter[diff]</TT>). Will be 
    <tt>value_type &</tt> for a mutable iterator, and convertible to 
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::pointer</tt></td>
    <td>the iterator's pointer type (return type of <TT>iter.operator->()</TT>). Will be 
    <tt>value_type *</tt> for a mutable iterator, and convertible to 
    <tt>value_type const *</tt> for a const iterator.</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::difference_type</tt></td>
    <td>the iterator's difference type (<TT>vigra::Diff2D</TT>)</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::iterator_category</tt></td>
    <td>the iterator tag (<tt>vigra::image_traverser_tag</tt>)</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::row_iterator</tt></td><td>the associated row iterator</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::column_iterator</tt></td><td>the associated column iterator</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::MoveX</td><td>type of the horizontal navigator</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::MoveY</td><td>type of the vertical navigator</td>
</tr>
<tr>
  <th>Operation</th><th>Result</th><th>Semantics</th>
</tr>
<tr>
    <td><tt>++i.x<br>i.x--</tt></td><td><tt>void</tt></td><td>increment x-coordinate</td>
</tr>
<tr>
    <td><tt>--i.x<br>i.x--</tt></td><td><tt>void</tt></td><td>decrement x-coordinate</td>
</tr>
<tr>
    <td><tt>i.x += dx</tt></td><td><tt>ImageIterator::MoveX &</tt></td>
    <td>add <tt>dx</tt> to x-coordinate</td>
</tr>
<tr>
    <td><tt>i.x -= dx</tt></td><td><tt>ImageIterator::MoveX &</tt></td>
    <td>subtract <tt>dx</tt> from x-coordinate</td>
</tr>
<tr>
    <td><tt>i.x - j.x</tt></td><td><tt>int</tt></td>
    <td>difference of the x-coordinates of <tt>i</tt> and <tt>j</tt></td>
</tr>
<tr>
    <td><tt>i.x = j.x</tt></td><td><tt>ImageIterator::MoveX &</tt></td><td><tt>i.x += j.x - i.x</tt></td>
</tr>
<tr>
    <td><tt>i.x == i.y</tt></td><td><tt>bool</tt></td><td><tt>j.x - i.x == 0</tt></td>

</tr>
<tr>
    <td><tt>i.x < j.x</tt></td><td><tt>bool</tt></td><td><tt>j.x - i.x > 0</tt></td>

</tr>
<tr>
    <td><tt>++i.y<br>i.y++</tt></td><td><tt>void</tt></td><td>increment y-coordinate</td>
</tr>
<tr>
    <td><tt>--i.y<br>i.y--</tt></td><td><tt>void</tt></td><td>decrement y-coordinate</td>
</tr>
<tr>
    <td><tt>i.y += dy</tt></td><td><tt>ImageIterator::MoveY &</tt></td>
    <td>add <tt>dy</tt> to y-coordinate</td>
</tr>
<tr>
    <td><tt>i.y -= dy</tt></td><td><tt>ImageIterator::MoveY &</tt></td>
    <td>subtract <tt>dy</tt> from y-coordinate</td>
</tr>
<tr>
    <td><tt>i.y - j.y</tt></td><td><tt>int</tt></td>
    <td>difference of the y-coordinates of <tt>i</tt> and <tt>j</tt></td>
</tr>
<tr>
    <td><tt>i.y = j.y</tt></td><td><tt>ImageIterator::MoveY &</tt></td><td><tt>i.y += j.y - i.y</tt></td>
</tr>
<tr>
    <td><tt>i.y == j.y</tt></td><td><tt>bool</tt></td><td><tt>j.y - i.y == 0</tt></td>

</tr>
<tr>
    <td><tt>i.y < j.y</tt></td><td><tt>bool</tt></td><td><tt>j.y - i.y > 0</tt></td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator k(i)</tt></td><td>copy constructor</td>
</tr>
<tr>
    <td><tt>k = i</tt></td><td><tt>ImageIterator &</tt></td><td>assignment</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator k</tt></td><td>default constructor</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::row_iterator r(i)</tt></td><td>construction of row iterator</td>
</tr>
<tr>
    <td colspan=2><tt>ImageIterator::column_iterator c(i)</tt></td><td>construction of column iterator</td>
</tr>
<tr>
    <td><tt>i += diff</tt></td><td><tt>ImageIterator &</tt></td>
    <td><tt>{ i.x += diff.x<br>i.y += diff.y; }</tt></td>
</tr>
<tr>
    <td><tt>i -= diff</tt></td><td><tt>ImageIterator &</tt></td>
    <td><tt>{ i.x -= diff.x<br>i.y -= diff.y; }</tt></td>
</tr>
<tr>
    <td><tt>i + diff</tt></td><td><tt>ImageIterator</tt></td>
    <td><tt>{ ImageIterator tmp(i);<br>tmp += diff;<br>return tmp; }</tt></td>
</tr>
<tr>
    <td><tt>i - diff</tt></td><td><tt>ImageIterator</tt></td>
    <td><tt>{ ImageIterator tmp(i);<br>tmp -= diff;<br>return tmp; }</tt></td>
</tr>
<tr>
    <td><tt>i - j</tt></td><td><tt>ImageIterator::difference_type</tt></td>
    <td><tt>{ ImageIterator::difference_type tmp(i.x - j.x, i.y - j.y);<br>return tmp; }</tt></td>
</tr>
<tr>
    <td><tt>i == j</tt></td><td><tt>bool</tt></td>
    <td><tt>i.x == j.x && i.y == j.y</tt></td>
</tr>
<tr>
    <td><tt>*i</tt></td><td><tt>ImageIterator::reference</tt></td>
    <td>access the current pixel</td>
</tr>
<tr>
    <td><tt>i[diff]</tt></td><td><tt>ImageIterator::index_reference</tt></td>
    <td>access pixel at offset <tt>diff</tt></td>
</tr>
<tr>
    <td><tt>i(dx, dy)</tt></td><td><tt>ImageIterator::index_reference</tt></td>
    <td>access pixel at offset <tt>(dx, dy)</tt></td>
</tr>
<tr>
    <td><tt>i->member()</tt></td><td>depends on operation</td>
    <td>call member function of underlying pixel type via <tt>operator-></tt> of iterator</td>
</tr>
<tr>
    <td colspan=3>
       <tt>i, j, k</tt> are of type <tt>ImageIterator</tt><br>
       <tt>diff</tt> is of type <tt>ImageIterator::difference_type</tt><br>
       <tt>dx, dy</tt> are of type <tt>int</tt><br>       
    </td>
</tr>
</table>
</p>
<h3>Requirements for Image Iterator Traits</h3>
<p>
The following iterator traits must be defined for an image iterator:
</p>
<p>
<table border=2 cellspacing=0 cellpadding=2 width="100%">
<tr>
    <th>Types</th><th>Meaning</th>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::Iterator</tt></td><td>the iterator type the traits are referring to</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::iterator</tt></td><td>the iterator type the traits are referring to</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::value_type</tt></td><td>the underlying image's pixel type</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::reference</tt></td>
    <td>the iterator's reference type (return type of <TT>*iter</TT>)</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::index_reference</tt></td>
    <td>the iterator's index reference type (return type of <TT>iter[diff]</TT>)</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::pointer</tt></td>
    <td>the iterator's pointer type (return type of <TT>iter.operator->()</TT>)</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::difference_type</tt></td>
    <td>the iterator's difference type</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::iterator_category</tt></td>
    <td>the iterator tag (<tt>vigra::image_traverser_tag</tt>)</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::row_iterator</tt></td><td>the associated row iterator</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::column_iterator</tt></td><td>the associated column iterator</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::DefaultAccessor</tt></td>
    <td>the default accessor to be used with the iterator</td>
</tr>
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::default_accessor</tt></td>
    <td>the default accessor to be used with the iterator</td>
</tr>
</table>
</p>
*/
//@{
/********************************************************/
/*                                                      */
/*                      ImageIteratorBase               */
/*                                                      */
/********************************************************/

/** \brief Base class for 2D random access iterators.

    This class contains the navigational part of the iterator.
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
    ImageIteratorBase(PixelType const * base, int width)
    : x(const_cast<PixelType *>(base)),
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
    ImageIteratorBase & operator+=(Diff2D const & s)
    {
        x += s.x;
        y += s.y;
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    ImageIteratorBase & operator-=(Diff2D const & s)
    {
        x -= s.x;
        y -= s.y;
        return *this;
    }
    
        /** access current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current() </TT>
        */
    PixelType & current()
    {
        return *(x.current_ + y.offset_);
    }
    
        /** access pixel at offset (dx, dy) from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(dx, dy) </TT>
        */
    PixelType & current(int const & dx, int const & dy)
    {
        return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }
    
        /** access pixel at offset d from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(Diff2D(1,1)) </TT>
        */
    PixelType & current(Diff2D const & d)
    {
        return *(x.current_ + d.x + y.offset_ + y.width_ * d.y);
    }
    
        /** read current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current() </TT>
        */
    PixelType const & current() const
    {
        return *(x.current_ + y.offset_);
    }

        /** read pixel with offset (dx, dy) from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(dx, dy) </TT>
        */
    PixelType const & current(int dx, int dy) const
    {
        return *(x.current_ + dx + y.offset_ + y.width_ * dy);
    }

        /** read pixel with offset d from current pixel. <br>
            usage: <TT> SomePixelType value = iterator.current(Diff2D(1,1)) </TT>
        */
    PixelType const & current(Diff2D const & d) const
    {
        return *(x.current_ + d.x + y.offset_ + y.width_ * d.y);
    }
    
    int width() const { return y.width_; }
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
    typedef PIXELTYPE value_type;

        /** The underlying image's pixel type.
        */
    typedef PIXELTYPE PixelType;
            
        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef PIXELTYPE &          reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef PIXELTYPE &          index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef PIXELTYPE *          pointer;
    
        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               difference_type;

        /** the iterator tag (image traverser)
        */
    typedef image_traverser_tag  iterator_category;
    
        /** The associated row iterator.
        */
    typedef PIXELTYPE * row_iterator;
    
        /** The associated column iterator.
        */
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<ImageIterator> > 
        column_iterator;
    
        /** Construct from raw memory with 'offset' between lines.
        If the raw memory is encapsulated in an image object this
        object should have a fatory function that constructs the
        iterator.
        */
    ImageIterator(pointer base, int offset)
    : ImageIteratorBase<PIXELTYPE>(base, offset)
    {}
    
        /** Copy constructor */
    ImageIterator(const ImageIterator & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Default constructor */
    ImageIterator()
    : ImageIteratorBase<PIXELTYPE>()
    {}    
    
        /** Copy assignment */
    ImageIterator & operator=(const ImageIterator & rhs)
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
    ImageIterator & operator+=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator+=(s);
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    ImageIterator & operator-=(Diff2D const & s)
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
    ImageIterator operator+(Diff2D const & s) const
    {
        ImageIterator ret(*this);

        ret += s;

        return ret;
    }
        /** Subtract a distance
        */
    ImageIterator operator-(Diff2D const & s) const
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
    reference operator*() const
    {
        return const_cast<reference>(current());
    }
    
        /** Call member of current pixel. <br>
            usage: <TT> iterator->pixelMemberFunction() </TT>
        */
    pointer operator->() const
    {
        return const_cast<pointer>(&current());
    }
    
        /** Access pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    index_reference operator[](Diff2D const & d) const
    {
        return const_cast<index_reference>(current(d));
    }
    
        /** Access pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    index_reference operator()(int dx, int dy) const
    {
        return const_cast<index_reference>(current(dx, dy));
    }
    
        /** Read pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    pointer operator[](int dy) const
    {
        return const_cast<pointer>(&current()+dy*width());
    }
    //@}
    
    row_iterator rowIterator() const
        { return const_cast<row_iterator>(&current()); }
        
    column_iterator columnIterator() const
    { 
        typedef typename column_iterator::BaseType Iter;
        return column_iterator(Iter(const_cast<pointer>(&current()), width()));
    }
    
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
    typedef PIXELTYPE value_type;

        /** The underlying image's pixel type.
        */
    typedef PIXELTYPE PixelType;
    
        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef PIXELTYPE const &    reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef PIXELTYPE const &    index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef PIXELTYPE const *    pointer;
    
        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               difference_type;

        /** the iterator tag (image traverser)
        */
    typedef image_traverser_tag  iterator_category;
    
        /** The associated row iterator.
        */
    typedef PIXELTYPE const * row_iterator;
    
        /** The associated column iterator.
        */
    typedef IteratorAdaptor<ContigousMemoryColumnIteratorPolicy<ConstImageIterator> > 
        column_iterator;
    
        /** Construct from raw memory with 'offset' between lines.
        If the raw memory is encapsulated in an image object this
        object should have a fatory function that constructs the
        iterator.
        */
    ConstImageIterator(pointer base, int offset)
    : ImageIteratorBase<PIXELTYPE>(base, offset)
    {}
    
        /** Copy constructor */
    ConstImageIterator(const ConstImageIterator & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Defult constructor */
    ConstImageIterator()
    : ImageIteratorBase<PIXELTYPE>()
    {}    
    
        /** Constructor from mutable iterator */
    ConstImageIterator(const ImageIterator<PIXELTYPE> & rhs)
    : ImageIteratorBase<PIXELTYPE>(rhs)
    {}    
    
        /** Copy assignment */
    ConstImageIterator & operator=(const ConstImageIterator & rhs)
    {
        if(this != &rhs)
        {
            ImageIteratorBase<PIXELTYPE>::operator=(rhs);
        }
        return *this;
    }
    
        /** Assignment from mutable iterator */
    ConstImageIterator & 
    operator=(const ImageIterator<PIXELTYPE> & rhs)
    {
        ImageIteratorBase<PIXELTYPE>::operator=(rhs);
        return *this;
    }

    /** @name Computed Assignment */
    //@{
        /** Add offset via Diff2D
        */
    ConstImageIterator & operator+=(Diff2D const & s)
    {
        ImageIteratorBase<PIXELTYPE>::operator+=(s);
        return *this;
    }
        /** Subtract offset via Diff2D
        */
    ConstImageIterator & operator-=(Diff2D const & s)
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
    ConstImageIterator operator+(Diff2D const & s) const
    {
        ConstImageIterator ret(*this);
    
        ret += s;
    
        return ret;
    }
        /** Subtract a distance
        */
    ConstImageIterator operator-(Diff2D const & s) const
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
    reference operator*() const
    {
       return current();
    }
    
        /** Call member of current pixel. <br>
            usage: <TT> iterator->pixelMemberFunction() </TT>
        */
    pointer operator->() const
    {
        return &current();
    }
    
        /** Read pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    index_reference operator()(int dx, int dy) const
    {
       return current(dx, dy);
    }

        /** Read pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    index_reference operator[](Diff2D const & d) const
    {
        return current(d);
    }
    
        /** Read pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    pointer operator[](int dy) const
    {
        return &current()+dy*width();
    }
    //@}
    
    row_iterator rowIterator() const
        { return &(this->current()); }
        
    column_iterator columnIterator() const
    { 
        typedef typename column_iterator::BaseType Iter;
        return column_iterator(Iter(&current(), width()));
    }
    
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T>
struct IteratorTraits<ImageIterator<T> >
{
    typedef ImageIterator<T>                    Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardAccessor<T>                  DefaultAccessor; 
    typedef StandardAccessor<T>                  default_accessor; 
};  

template <class T>
struct IteratorTraits<ConstImageIterator<T> >
{
    typedef ConstImageIterator<T>               Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardConstAccessor<T>             DefaultAccessor; 
    typedef StandardConstAccessor<T>             default_accessor; 
};  

#endif

template <class PIXELTYPE>
class ConstValueIteratorPolicy
{
  public:
    
    typedef PIXELTYPE                       value_type;
    typedef int                             difference_type;
    typedef PIXELTYPE const &               reference;
    typedef PIXELTYPE const &               index_reference;
    typedef PIXELTYPE const *               pointer;
    typedef std::random_access_iterator_tag iterator_category;
    
    struct BaseType
    {
        BaseType(PIXELTYPE const & v = PIXELTYPE(), int p = 0)
        : value(v), pos(p)
        {}
        
        PIXELTYPE value;
        int pos;
    };
    
    static void initialize(BaseType & d) {}
    
    static reference dereference(BaseType const & d)
        { return d.value; }
    
    static index_reference dereference(BaseType d, difference_type)
    { 
        return d.value; 
    }
    
    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.pos == d2.pos; }
    
    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.pos < d2.pos; }
    
    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return d1.pos - d2.pos; }
    
    static void increment(BaseType & d)
        { ++d.pos; }
    
    static void decrement(BaseType & d)
        { --d.pos; }
    
    static void advance(BaseType & d, difference_type n)
        { d.pos += n; }
};

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
        /** The type of the constant the iterator holds.
        */
   typedef PIXELTYPE value_type;

        /** The type of the constant the iterator holds.
        */
    typedef PIXELTYPE PixelType;
    
        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef PIXELTYPE const &    reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef PIXELTYPE const &    index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef PIXELTYPE const *    pointer;
    
        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               difference_type;

        /** the iterator tag (image traverser)
        */
    typedef image_traverser_tag  iterator_category;
    
        /** The associated row iterator.
        */
    typedef IteratorAdaptor<ConstValueIteratorPolicy<PIXELTYPE> > row_iterator;
    
        /** The associated column iterator.
        */
    typedef IteratorAdaptor<ConstValueIteratorPolicy<PIXELTYPE> > column_iterator;
    
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
    reference operator*() const
    {
        return value_;
    }
    
        /** Call member function for stored constant.
        */
    pointer operator->() const
    {
        return &value_;
    }
    
        /** Read pixel at a distance (return specified constant).
        */
    index_reference operator()(int const &, int const &) const
    {
        return value_;
    }

        /** Read pixel at a distance (return specified constant).
        */
    index_reference operator[](Diff2D const &) const
    {
        return value_;
    }
    
        /** Get row iterator at current position (which will also hold the constant).
        */
    row_iterator rowIterator() const
        { return row_iterator(typename row_iterator::BaseType(value_, x)); }
  
        /** Get column iterator at current position (which will also hold the constant).
        */
    column_iterator columnIterator() const
        { return column_iterator(typename column_iterator::BaseType(value_, y)); }
  
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
    typedef ConstValueIterator<T>                Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardConstAccessor<T>             DefaultAccessor; 
    typedef StandardConstAccessor<T>             default_accessor; 
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
