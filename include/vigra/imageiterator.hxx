/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_IMAGEITERATOR_HXX
#define VIGRA_IMAGEITERATOR_HXX

#include "utilities.hxx"
#include "accessor.hxx"
#include "iteratortraits.hxx"
#include "metaprogramming.hxx"

namespace vigra {

template <class IMAGEITERATOR>
class StridedIteratorPolicy
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
        explicit BaseType(pointer c = 0, difference_type stride = 0)
        : current_(c), stride_(stride)
        {}

        pointer current_;
        difference_type stride_;
    };

    static void initialize(BaseType & /* d */) {}

    static reference dereference(BaseType const & d)
        { return const_cast<reference>(*d.current_); }

    static index_reference dereference(BaseType const & d, difference_type n)
    {
        return const_cast<index_reference>(d.current_[n*d.stride_]);
    }

    static bool equal(BaseType const & d1, BaseType const & d2)
        { return d1.current_ == d2.current_; }

    static bool less(BaseType const & d1, BaseType const & d2)
        { return d1.current_ < d2.current_; }

    static difference_type difference(BaseType const & d1, BaseType const & d2)
        { return (d1.current_ - d2.current_) / d1.stride_; }

    static void increment(BaseType & d)
        { d.current_ += d.stride_; }

    static void decrement(BaseType & d)
        { d.current_ -= d.stride_; }

    static void advance(BaseType & d, difference_type n)
        { d.current_ += d.stride_*n; }
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
<tr><td>
    \htmlonly
    <th colspan=2>
    \endhtmlonly
    Local Types
    \htmlonly
    </th><th>
    \endhtmlonly
    Meaning
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
<tr><td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::value_type</tt></td><td>the underlying image's pixel type</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::PixelType</tt></td><td>the underlying image's pixel type</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::reference</tt></td>
    <td>the iterator's reference type (return type of <TT>*iter</TT>). Will be
    <tt>value_type &</tt> for a mutable iterator, and convertible to
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::index_reference</tt></td>
    <td>the iterator's index reference type (return type of <TT>iter[diff]</TT>). Will be
    <tt>value_type &</tt> for a mutable iterator, and convertible to
    <tt>value_type const &</tt> for a const iterator.</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::pointer</tt></td>
    <td>the iterator's pointer type (return type of <TT>iter.operator->()</TT>). Will be
    <tt>value_type *</tt> for a mutable iterator, and convertible to
    <tt>value_type const *</tt> for a const iterator.</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::difference_type</tt></td>
    <td>the iterator's difference type (<TT>vigra::Diff2D</TT>)</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::iterator_category</tt></td>
    <td>the iterator tag (<tt>vigra::image_traverser_tag</tt>)</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::row_iterator</tt></td><td>the associated row iterator</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::column_iterator</tt></td><td>the associated column iterator</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::MoveX</tt></td><td>type of the horizontal navigator</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::MoveY</tt></td><td>type of the vertical navigator</td>
</tr>
<tr><td>
    \htmlonly
    <th>
    \endhtmlonly
    Operation
    \htmlonly
    </th><th>
    \endhtmlonly
    Result
    \htmlonly
    </th><th>
    \endhtmlonly
    Semantics
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
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
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator k(i)</tt></td><td>copy constructor</td>
</tr>
<tr>
    <td><tt>k = i</tt></td><td><tt>ImageIterator &</tt></td><td>assignment</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator k</tt></td><td>default constructor</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::row_iterator r(i)</tt></td><td>construction of row iterator</td>
</tr>
<tr>
    <td>
    \htmlonly
    <td colspan=2>
    \endhtmlonly
    <tt>ImageIterator::column_iterator c(i)</tt></td><td>construction of column iterator</td>
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
    <td>
    \htmlonly
    <td colspan=3>
    \endhtmlonly
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
<tr><td>
    \htmlonly
    <th>
    \endhtmlonly
    Types
    \htmlonly
    </th><th>
    \endhtmlonly
    Meaning
    \htmlonly
    </th>
    \endhtmlonly
</td></tr>
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
<tr>
    <td><tt>IteratorTraits&lt;ImageIterator&gt;::hasConstantStrides</tt></td>
    <td>whether the iterator uses constant strides on the underlying memory
        (always <tt>VigraTrueType</tt> for <tt>ImageIterator</tt>s).</td>
</tr>
</table>
</p>
*/
//@{

namespace detail {

template <class StridedOrUnstrided>
class DirectionSelector;

template <>
class DirectionSelector<UnstridedArrayTag>
{
  public:

    template <class T>
    class type
    {
      public:
        type(T base)
        : current_(base)
        {}

        type(type const & rhs)
        : current_(rhs.current_)
        {}

        type & operator=(type const & rhs)
        {
            current_ = rhs.current_;
            return *this;
        }

        void operator++() {++current_;}
        void operator++(int) {++current_;}
        void operator--() {--current_;}
        void operator--(int) {--current_;}
        void operator+=(int dx) {current_ += dx; }
        void operator-=(int dx) {current_ -= dx; }

        bool operator==(type const & rhs) const
         { return current_ == rhs.current_; }

        bool operator!=(type const & rhs) const
         { return current_ != rhs.current_; }

        bool operator<(type const & rhs) const
         { return current_ < rhs.current_; }

        bool operator<=(type const & rhs) const
         { return current_ <= rhs.current_; }

        bool operator>(type const & rhs) const
         { return current_ > rhs.current_; }

        bool operator>=(type const & rhs) const
         { return current_ >= rhs.current_; }

        int operator-(type const & rhs) const
         { return current_ - rhs.current_; }

        T operator()() const
        { return current_; }

        T operator()(int d) const
        { return current_ + d; }

        T current_;
    };
};

template <>
class DirectionSelector<StridedArrayTag>
{
  public:

    template <class T>
    class type
    {
      public:
        type(int stride, T base = 0)
        : stride_(stride),
          current_(base)
        {}

        type(type const & rhs)
        : stride_(rhs.stride_),
          current_(rhs.current_)
        {}

        type & operator=(type const & rhs)
        {
            stride_ = rhs.stride_;
            current_ = rhs.current_;
            return *this;
        }

        void operator++() {current_ += stride_; }
        void operator++(int) {current_ += stride_; }
        void operator--() {current_ -= stride_; }
        void operator--(int) {current_ -= stride_; }
        void operator+=(int dy) {current_ += dy*stride_; }
        void operator-=(int dy) {current_ -= dy*stride_; }

        bool operator==(type const & rhs) const
         { return (current_ == rhs.current_); }

        bool operator!=(type const & rhs) const
         { return (current_ != rhs.current_); }

        bool operator<(type const & rhs) const
         { return (current_ < rhs.current_); }

        bool operator<=(type const & rhs) const
         { return (current_ <= rhs.current_); }

        bool operator>(type const & rhs) const
         { return (current_ > rhs.current_); }

        bool operator>=(type const & rhs) const
         { return (current_ >= rhs.current_); }

        int operator-(type const & rhs) const
         { return (current_ - rhs.current_) / stride_; }

        T operator()() const
        { return current_; }

        T operator()(int d) const
        { return current_ + d*stride_; }

        int stride_;
        T current_;
    };
};

template <class StridedOrUnstrided>
class LinearIteratorSelector;

template <>
class LinearIteratorSelector<UnstridedArrayTag>
{
  public:
    template <class IMAGEITERATOR>
    class type
    {
      public:
        typedef typename IMAGEITERATOR::pointer res;

        template <class DirSelect>
        static res construct(typename IMAGEITERATOR::pointer data, DirSelect const &)
        {
            return data;
        }
    };
};

template <>
class LinearIteratorSelector<StridedArrayTag>
{
  public:
    template <class IMAGEITERATOR>
    class type
    {
      public:
        typedef IteratorAdaptor<StridedIteratorPolicy<IMAGEITERATOR> > res;

        template <class DirSelect>
        static res construct(typename IMAGEITERATOR::pointer data, DirSelect const & d)
        {
            typedef typename res::BaseType Base;
            return res(Base(data, d.stride_));
        }
    };
};


} // namespace detail

/********************************************************/
/*                                                      */
/*                      ImageIteratorBase               */
/*                                                      */
/********************************************************/

/** \brief Base class for 2D random access iterators.

    This class contains the navigational part of the iterator.
    It is usually not constructed directly, but via some derived class such as
    \ref ImageIterator or \ref StridedImageIterator.

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"

    Namespace: vigra

    The usage examples assume that you constructed two iterators like
    this:

    \code
    vigra::ImageIterator<SomePixelType> iterator(base, width);
    vigra::ImageIterator<SomePixelType> iterator1(base, width);
    \endcode

    See the paper: U. Koethe:
    <a href="documents/GenericProg2D.ps">Reusable Algorithms in Image Processing</a>
    for a discussion of the concepts behind ImageIterators.

*/
template <class IMAGEITERATOR,
          class PIXELTYPE, class REFERENCE, class POINTER,
          class StridedOrUnstrided = UnstridedArrayTag>
class ImageIteratorBase
{
    typedef typename
        detail::LinearIteratorSelector<StridedOrUnstrided>::template type<ImageIteratorBase>
        RowIteratorSelector;
    typedef typename
        detail::LinearIteratorSelector<StridedArrayTag>::template type<ImageIteratorBase>
        ColumnIteratorSelector;
  public:
    typedef ImageIteratorBase<IMAGEITERATOR,
                 PIXELTYPE, REFERENCE, POINTER, StridedOrUnstrided> self_type;

        /** The underlying image's pixel type.
        */
    typedef PIXELTYPE value_type;

        /** deprecated, use <TT>value_type</TT> instead.
        */
    typedef PIXELTYPE PixelType;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef REFERENCE            reference;

        /** the iterator's index reference type (return type of <TT>iter[diff]</TT>)
        */
    typedef REFERENCE            index_reference;

        /** the iterator's pointer type (return type of <TT>iter.operator->()</TT>)
        */
    typedef POINTER              pointer;

        /** the iterator's difference type (argument type of <TT>iter[diff]</TT>)
        */
    typedef Diff2D               difference_type;

        /** the iterator tag (image traverser)
        */
    typedef image_traverser_tag  iterator_category;

        /** The associated row iterator.
        */
    typedef typename RowIteratorSelector::res row_iterator;

        /** The associated column iterator.
        */
    typedef typename ColumnIteratorSelector::res column_iterator;

        /** Let operations act in X direction
        */
    typedef typename
        detail::DirectionSelector<StridedOrUnstrided>::template type<pointer> MoveX;

        /** Let operations act in Y direction
        */
    typedef typename
        detail::DirectionSelector<StridedArrayTag>::template type<int> MoveY;

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
    difference_type operator-(ImageIteratorBase const & rhs) const
    {
        return difference_type(x - rhs.x, y - rhs.y);
    }

    //@}

    /** @name Specify coordinate to operate on */
    //@{
        /** Refer to iterator's x coordinate.
            Usage examples:<br>
            \code
            ++iterator.x;        // move one step to the right
            --iterator.x;        // move one step to the left
            iterator.x += dx;    // move dx steps to the right
            iterator.x -= dx;    // move dx steps to the left
            bool notAtEndOfRow = iterator.x < lowerRight.x;   // compare x coordinates of two iterators
            int width = lowerRight.x - upperLeft.x;           // calculate difference of x coordinates
                                                              // between two iterators
            \endcode
        */
    MoveX x;
        /** Refer to iterator's y coordinate.
            Usage examples:<br>
            \code
            ++iterator.y;        // move one step down
            --iterator.y;        // move one step up
            iterator.y += dy;    // move dy steps down
            iterator.y -= dy;    // move dy steps up
            bool notAtEndOfColumn = iterator.y < lowerRight.y; // compare y coordinates of two iterators
            int height = lowerRight.y - upperLeft.y;           // calculate difference of y coordinates
                                                               // between two iterators
            \endcode
        */
    MoveY y;
    //@}

  protected:
        /** Construct from raw memory with a vertical stride of <TT>ystride</TT>.
        <TT>ystride</TT> must equal the physical image width (row length),
        even if the iterator will only be used for a sub image. This constructor
        must only be called for unstrided iterators
        (<tt>StridedOrUnstrided == UnstridedArrayTag</tt>)
        */
    ImageIteratorBase(pointer base, int ystride)
    : x(base),
      y(ystride)
    {}

        /** Construct from raw memory with a horizontal stride of <TT>xstride</TT>
        and a vertical stride of <TT>ystride</TT>. This constructor
        may be used for iterators that shall skip pixels. Thus, it
        must only be called for strided iterators
        (<tt>StridedOrUnstrided == StridedArrayTag</tt>)
        */
    ImageIteratorBase(pointer base, int xstride, int ystride)
    : x(xstride, base),
      y(ystride)
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

  public:
    /** @name Random navigation */
    //@{
        /** Add offset via Diff2D
        */
    IMAGEITERATOR & operator+=(difference_type const & s)
    {
        x += s.x;
        y += s.y;
        return static_cast<IMAGEITERATOR &>(*this);
    }
        /** Subtract offset via Diff2D
        */
    IMAGEITERATOR & operator-=(difference_type const & s)
    {
        x -= s.x;
        y -= s.y;
        return static_cast<IMAGEITERATOR &>(*this);
    }

        /** Add a distance
        */
    IMAGEITERATOR operator+(difference_type const & s) const
    {
        IMAGEITERATOR ret(static_cast<IMAGEITERATOR const &>(*this));

        ret += s;

        return ret;
    }

        /** Subtract a distance
        */
    IMAGEITERATOR operator-(difference_type const & s) const
    {
        IMAGEITERATOR ret(static_cast<IMAGEITERATOR const &>(*this));

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
        return *current();
    }

        /** Call member of current pixel. <br>
            usage: <TT> iterator->pixelMemberFunction() </TT>
        */
    pointer operator->() const
    {
        return current();
    }

        /** Access pixel at offset from current location. <br>
            usage: <TT> SomePixelType value = iterator[Diff2D(1,1)] </TT>
        */
    index_reference operator[](Diff2D const & d) const
    {
        return *current(d.x, d.y);
    }

        /** Access pixel at offset (dx, dy) from current location. <br>
            usage: <TT> SomePixelType value = iterator(dx, dy) </TT>
        */
    index_reference operator()(int dx, int dy) const
    {
        return *current(dx, dy);
    }

        /** Read pixel with offset [dy][dx] from current pixel.
            Note that the 'x' index is the trailing index. <br>
            usage: <TT> SomePixelType value = iterator[dy][dx] </TT>
        */
    pointer operator[](int dy) const
    {
        return x() + y(dy);
    }
    //@}

    row_iterator rowIterator() const
    {
        return RowIteratorSelector::construct(current(), x);
    }

    column_iterator columnIterator() const
    {
        return ColumnIteratorSelector::construct(current(), y);
    }

  private:

    pointer current() const
        { return x() + y(); }

    pointer current(int dx, int dy) const
        { return x(dx) + y(dy); }
};

/********************************************************/
/*                                                      */
/*                      ImageIterator                   */
/*                                                      */
/********************************************************/

/** \brief Standard 2D random access iterator for images that store the
    data in a linear array.

    Most functions and local types are inherited from ImageIteratorBase.

    See the paper: U. Koethe:
    <a href="documents/GenericProg2D.ps">Reusable Algorithms in Image Processing</a>
    for a discussion of the concepts behind ImageIterators.

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"

    Namespace: vigra

*/
template <class PIXELTYPE>
class ImageIterator
: public ImageIteratorBase<ImageIterator<PIXELTYPE>,
                           PIXELTYPE, PIXELTYPE &, PIXELTYPE *>
{
  public:
    typedef ImageIteratorBase<ImageIterator<PIXELTYPE>,
                              PIXELTYPE, PIXELTYPE &, PIXELTYPE *> Base;

    typedef typename Base::pointer         pointer;
    typedef typename Base::difference_type difference_type;

        /** Construct from raw memory with a vertical stride of <TT>ystride</TT>.
        <TT>ystride</TT> must equal the physical image width (row length),
        even if the iterator will only be used for a sub image.
        If the raw memory is encapsulated in an image object this
        object should have a factory function that constructs the
        iterator.
        */
    ImageIterator(pointer base, int ystride)
    : Base(base, ystride)
    {}

        /** Default constructor */
    ImageIterator()
    : Base()
    {}

};

/********************************************************/
/*                                                      */
/*                   ConstImageIterator                 */
/*                                                      */
/********************************************************/

/** \brief Standard 2D random access const iterator for images that
    store the data as a linear array.

    Most functions are inherited from ImageIteratorBase.

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"

    Namespace: vigra

*/
template <class PIXELTYPE>
class ConstImageIterator
: public ImageIteratorBase<ConstImageIterator<PIXELTYPE>,
                           PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *>
{
  public:
    typedef ImageIteratorBase<ConstImageIterator<PIXELTYPE>,
                        PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *> Base;

    typedef typename Base::pointer         pointer;
    typedef typename Base::difference_type difference_type;

        /** Construct from raw memory with a vertical stride of <TT>ystride</TT>.
        <TT>ystride</TT> must equal the physical image width (row length),
        even if the iterator will only be used for a sub image.
        If the raw memory is encapsulated in an image object this
        object should have a factory function that constructs the
        iterator.
        */
    ConstImageIterator(pointer base, int ystride)
    : Base(base, ystride)
    {}

    ConstImageIterator(ImageIterator<PIXELTYPE> const & o)
    : Base(o.x, o.y)
    {}

        /** Default constructor */
    ConstImageIterator()
    : Base()
    {}

    ConstImageIterator & operator=(ImageIterator<PIXELTYPE> const & o)
    {
        Base::x = o.x;
        Base::y = o.y;
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*                 StridedImageIterator                 */
/*                                                      */
/********************************************************/

/** \brief Iterator to be used when pixels are to be skipped.

    This iterator can be used when some pixels shall be automatically skipped, for example
    if an image is to be sub-sampled: instead of advancing to the next pixel,
    <tt>++iterator.x</tt> jumps to the pixel at a horizontal offset of <tt>xskip</tt>.
    Likewise with <tt>yskip</tt> in vertical direction. Most functions and local types
    are inherited from ImageIteratorBase.

    <b> Usage:</b>

    \code
    BImage img(w,h);
    ...
    int xskip = 2, yskip = 2;
    int wskip = w / xskip + 1, hskip = h / yskip + 1;

    StridedImageIterator<BImage::value_type> upperLeft(&img(0,0), w, xskip, yskip);
    StridedImageIterator<BImage::value_type> lowerRight = upperLeft + Diff2D(wskip, hskip);

    // now navigation with upperLeft and lowerRight lets the image appear to have half
    // the original resolution in either dimension
    \endcode

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"

    Namespace: vigra

*/
template <class PIXELTYPE>
class StridedImageIterator
: public ImageIteratorBase<StridedImageIterator<PIXELTYPE>,
                           PIXELTYPE, PIXELTYPE &, PIXELTYPE *, StridedArrayTag>
{
  public:
    typedef ImageIteratorBase<StridedImageIterator<PIXELTYPE>,
                              PIXELTYPE, PIXELTYPE &, PIXELTYPE *, StridedArrayTag> Base;

    typedef typename Base::pointer         pointer;
    typedef typename Base::difference_type difference_type;

        /** Construct from raw memory with a vertical stride of <TT>ystride</TT>,
        jumping by <tt>xskip</tt> horizontally and <tt>yskip</tt> vertically.
        <tt>ystride</tt> must be the physical width (row length) of the image.
        */
    StridedImageIterator(pointer base, int ystride, int xskip, int yskip)
    : Base(base, xskip, ystride*yskip)
    {}

        /** Default constructor */
    StridedImageIterator()
    : Base()
    {}

};

/********************************************************/
/*                                                      */
/*               ConstStridedImageIterator              */
/*                                                      */
/********************************************************/

/** \brief Const iterator to be used when pixels are to be skipped.

    This iterator can be used when some pixels shall be automatically skipped, for example
    if an image is to be sub-sampled: instead of advancing to the next pixel,
    <tt>++iterator.x</tt> jumps to the pixel at a horizontal offset of <tt>xskip</tt>.
    Likewise with <tt>yskip</tt> in vertical direction. Most functions and local types
    are inherited from ImageIteratorBase.

    <b> Usage:</b>

    \code
    BImage img(w,h);
    ...
    int xskip = 2, yskip = 2;
    int wskip = w / xskip + 1, hskip = h / yskip + 1;

    ConstStridedImageIterator<BImage::value_type> upperLeft(&img(0,0), w, xskip, yskip);
    ConstStridedImageIterator<BImage::value_type> lowerRight = upperLeft + Diff2D(wskip, hskip);

    // now navigation with upperLeft and lowerRight lets the image appear to have half
    // the original resolution in either dimension
    \endcode

    <b>\#include</b> "<a href="imageiterator_8hxx-source.html">vigra/imageiterator.hxx</a>"

    Namespace: vigra

*/
template <class PIXELTYPE>
class ConstStridedImageIterator
: public ImageIteratorBase<ConstStridedImageIterator<PIXELTYPE>,
                           PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *,
                           StridedArrayTag>
{
  public:
    typedef ImageIteratorBase<ConstStridedImageIterator<PIXELTYPE>,
                        PIXELTYPE, PIXELTYPE const &, PIXELTYPE const *,
                        StridedArrayTag> Base;

    typedef typename Base::pointer         pointer;
    typedef typename Base::difference_type difference_type;

        /** Construct from raw memory with a vertical stride of <TT>ystride</TT>,
        jumping by <tt>xskip</tt> horizontally and <tt>yskip</tt> vertically.
        <tt>ystride</tt> must be the physical width (row length) of the image.
        */
    ConstStridedImageIterator(pointer base, int ystride, int xskip, int yskip)
    : Base(base, xskip, ystride*yskip)
    {}

        /** Copy-construct from mutable iterator */
    ConstStridedImageIterator(StridedImageIterator<PIXELTYPE> const & o)
    : Base(o.x, o.y)
    {}

        /** Default constructor */
    ConstStridedImageIterator()
    : Base()
    {}

        /** Assign mutable iterator */
    ConstStridedImageIterator & operator=(StridedImageIterator<PIXELTYPE> const & o)
    {
        Base::x = o.x;
        Base::y = o.y;
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*             definition of iterator traits            */
/*                                                      */
/********************************************************/


#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T>
struct IteratorTraits<ImageIterator<T> >
: public IteratorTraitsBase<ImageIterator<T> >
{
    typedef ImageIterator<T>                              mutable_iterator;
    typedef ConstImageIterator<T>                         const_iterator;
    typedef typename AccessorTraits<T>::default_accessor  DefaultAccessor;
    typedef DefaultAccessor                               default_accessor;
    typedef VigraTrueType                                 hasConstantStrides;
};

template <class T>
struct IteratorTraits<ConstImageIterator<T> >
: public IteratorTraitsBase<ConstImageIterator<T> >
{
    typedef ImageIterator<T>                              mutable_iterator;
    typedef ConstImageIterator<T>                         const_iterator;
    typedef typename AccessorTraits<T>::default_const_accessor  DefaultAccessor;
    typedef DefaultAccessor                               default_accessor;
    typedef VigraTrueType                                 hasConstantStrides;
};

template <class T>
struct IteratorTraits<StridedImageIterator<T> >
: public IteratorTraitsBase<StridedImageIterator<T> >
{
    typedef StridedImageIterator<T>                       mutable_iterator;
    typedef ConstStridedImageIterator<T>                  const_iterator;
    typedef typename AccessorTraits<T>::default_accessor  DefaultAccessor;
    typedef DefaultAccessor                               default_accessor;
    typedef VigraTrueType                                 hasConstantStrides;
};

template <class T>
struct IteratorTraits<ConstStridedImageIterator<T> >
: public IteratorTraitsBase<ConstStridedImageIterator<T> >
{
    typedef StridedImageIterator<T>                       mutable_iterator;
    typedef ConstStridedImageIterator<T>                  const_iterator;
    typedef typename AccessorTraits<T>::default_const_accessor  DefaultAccessor;
    typedef DefaultAccessor                               default_accessor;
    typedef VigraTrueType                                 hasConstantStrides;
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

#define VIGRA_DEFINE_ITERATORTRAITS(VALUETYPE) \
    template <>  \
    struct IteratorTraits<ImageIterator<VALUETYPE > > \
    : public IteratorTraitsBase<ImageIterator<VALUETYPE > > \
    { \
        typedef ImageIterator<VALUETYPE>                         mutable_iterator; \
        typedef ConstImageIterator<VALUETYPE>                    const_iterator; \
        typedef typename AccessorTraits<VALUETYPE >::default_accessor  DefaultAccessor; \
        typedef DefaultAccessor                               default_accessor; \
        typedef VigraTrueType                                 hasConstantStrides; \
    }; \
    \
    template <>  \
    struct IteratorTraits<ConstImageIterator<VALUETYPE > > \
    : public IteratorTraitsBase<ConstImageIterator<VALUETYPE > > \
    { \
        typedef ImageIterator<VALUETYPE>                         mutable_iterator; \
        typedef ConstImageIterator<VALUETYPE>                    const_iterator; \
        typedef typename AccessorTraits<VALUETYPE >::default_const_accessor  DefaultAccessor; \
        typedef DefaultAccessor                               default_accessor; \
        typedef VigraTrueType                                 hasConstantStrides; \
    }; \
    template <>  \
    struct IteratorTraits<StridedImageIterator<VALUETYPE > > \
    : public IteratorTraitsBase<StridedImageIterator<VALUETYPE > > \
    { \
        typedef StridedImageIterator<VALUETYPE>                         mutable_iterator; \
        typedef ConstStridedImageIterator<VALUETYPE>                    const_iterator; \
        typedef typename AccessorTraits<VALUETYPE >::default_accessor  DefaultAccessor; \
        typedef DefaultAccessor                               default_accessor; \
        typedef VigraTrueType                                 hasConstantStrides; \
    }; \
    \
    template <>  \
    struct IteratorTraits<ConstStridedImageIterator<VALUETYPE > > \
    : public IteratorTraitsBase<ConstStridedImageIterator<VALUETYPE > > \
    { \
        typedef StridedImageIterator<VALUETYPE>                         mutable_iterator; \
        typedef ConstStridedImageIterator<VALUETYPE>                    const_iterator; \
        typedef typename AccessorTraits<VALUETYPE >::default_const_accessor  DefaultAccessor; \
        typedef DefaultAccessor                               default_accessor; \
        typedef VigraTrueType                                 hasConstantStrides; \
    };

VIGRA_DEFINE_ITERATORTRAITS(RGBValue<unsigned char>)
VIGRA_DEFINE_ITERATORTRAITS(RGBValue<short>)
VIGRA_DEFINE_ITERATORTRAITS(RGBValue<int>)
VIGRA_DEFINE_ITERATORTRAITS(RGBValue<float>)
VIGRA_DEFINE_ITERATORTRAITS(RGBValue<double>)

#define VIGRA_PIXELTYPE TinyVector<unsigned char, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<unsigned char, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<unsigned char, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<short, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<int, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<float, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 2>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 3>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE
#define VIGRA_PIXELTYPE TinyVector<double, 4>
VIGRA_DEFINE_ITERATORTRAITS(VIGRA_PIXELTYPE)
#undef VIGRA_PIXELTYPE

#undef VIGRA_DEFINE_ITERATORTRAITS

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

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
    typedef ConstValueIterator<T>                  Iterator;
    typedef Iterator                               iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
    typedef StandardConstAccessor<T>               DefaultAccessor;
    typedef StandardConstAccessor<T>               default_accessor;
    typedef VigraTrueType                                 hasConstantStrides;
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
