/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
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

#ifndef VIGRA_BASICIMAGEVIEW_HXX
#define VIGRA_BASICIMAGEVIEW_HXX

#include "imageiterator.hxx"
#include "initimage.hxx"

// Bounds checking Macro used if VIGRA_CHECK_BOUNDS is defined.
#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(this->isInside(diff), "Index out of bounds")
#else
#define VIGRA_ASSERT_INSIDE(diff)
#endif


namespace vigra {


/********************************************************/
/*                                                      */
/*                     BasicImageView                   */
/*                                                      */
/********************************************************/

/** \brief BasicImage using foreign memory.

    This class provides the same interface as \ref vigra::BasicImage
    (with the exception of <tt>resize()</tt>) but the image's
    memory is provided from the outside instead of allocated internally.

    A <tt>BasicImageView</tt> can also be created from a
    \ref vigra::MultiArrayView with the appropriate shape -- see
    \ref MultiArrayToImage.

    <b>\#include</b> \<vigra/basicimageview.hxx\> <br/>
    Namespace: vigra
*/
template <class PIXELTYPE>
class BasicImageView
{
  public:

        /** the BasicImageView's pixel type
        */
    typedef PIXELTYPE value_type;

        /** the BasicImageView's pixel type
        */
    typedef PIXELTYPE PixelType;

        /** the BasicImageView's reference type (i.e. the
            return type of image[diff] and image(dx,dy))
        */
    typedef PIXELTYPE &       reference;

        /** the BasicImageView's const reference type (i.e. the
            return type of image[diff] and image(dx,dy) when image is const)
        */
    typedef PIXELTYPE const & const_reference;

        /** the BasicImageView's pointer type
        */
    typedef PIXELTYPE *       pointer;

        /** the BasicImageView's const pointer type
        */
    typedef PIXELTYPE const * const_pointer;

        /** the BasicImageView's 1D random access iterator
            (note: lower case 'iterator' is a STL compatible 1D random
             access iterator, don't confuse with capitalized Iterator)
        */
    typedef PIXELTYPE * iterator;

        /** deprecated, use <TT>iterator</TT> instead
        */
    typedef PIXELTYPE * ScanOrderIterator;

        /** the BasicImageView's 1D random access const iterator
            (note: lower case 'const_iterator' is a STL compatible 1D
            random access const iterator)
        */
    typedef PIXELTYPE const * const_iterator;

        /** deprecated, use <TT>const_iterator</TT> instead
        */
    typedef PIXELTYPE const * ConstScanOrderIterator;

        /** the BasicImageView's 2D random access iterator ('traverser')
        */
    typedef ImageIterator<value_type> traverser;

        /** deprecated, use <TT>traverser</TT> instead
        */
    typedef ImageIterator<value_type> Iterator;

        /** the BasicImageView's 2D random access const iterator ('const traverser')
        */
    typedef ConstImageIterator<value_type> const_traverser;

        /** deprecated, use <TT>const_traverser</TT> instead
        */
    typedef ConstImageIterator<value_type> ConstIterator;

        /** the row iterator associated with the traverser
        */
    typedef typename traverser::row_iterator row_iterator;

        /** the const row iterator associated with the const_traverser
        */
    typedef typename const_traverser::row_iterator const_row_iterator;

        /** the column iterator associated with the traverser
        */
    typedef typename traverser::column_iterator column_iterator;

        /** the const column iterator associated with the const_traverser
        */
    typedef typename const_traverser::column_iterator const_column_iterator;

        /** the BasicImageView's difference type (argument type of image[diff])
        */
    typedef Diff2D difference_type;

         /** the BasicImageView's size type (result type of image.size())
        */
    typedef Size2D size_type;

       /** the BasicImageView's default accessor
        */
    typedef typename
          IteratorTraits<traverser>::DefaultAccessor Accessor;

        /** the BasicImageView's default const accessor
        */
    typedef typename
          IteratorTraits<const_traverser>::DefaultAccessor ConstAccessor;

        /** construct image of size 0x0
        */
    BasicImageView()
    : data_(0),
      width_(0),
      height_(0),
      stride_(0)
    {}

        /** construct view of size w x h
        */
    BasicImageView(const_pointer data, std::ptrdiff_t w, std::ptrdiff_t h, std::ptrdiff_t stride = 0)
    : data_(const_cast<pointer>(data)),
      width_(w),
      height_(h),
      stride_(stride == 0 ? w : stride)
    {}

        /** construct view of size size.x x size.y
        */
    BasicImageView(const_pointer data, difference_type const & size, std::ptrdiff_t stride = 0)
    : data_(const_cast<pointer>(data)),
      width_(size.x),
      height_(size.y),
      stride_(stride == 0 ? size.x : stride)
    {}

        /** set Image with const value
        */
    BasicImageView & init(value_type const & pixel)
    {
        initImage(upperLeft(), lowerRight(), accessor(), pixel);

        return *this;
    }

        /** width of Image
        */
    std::ptrdiff_t width() const
    {
        return width_;
    }

        /** height of Image
        */
    std::ptrdiff_t height() const
    {
        return height_;
    }

        /** stride of Image.
            Memory offset between the start of two successive rows.
        */
    std::ptrdiff_t stride() const
    {
        return stride_;
    }

        /** size of Image
        */
    size_type size() const
    {
        return size_type(width(), height());
    }

        /** test whether a given coordinate is inside the image
        */
    bool isInside(difference_type const & d) const
    {
        return d.x >= 0 && d.y >= 0 &&
               d.x < width() && d.y < height();
    }

        /** access pixel at given location. <br>
        usage: <TT> value_type value = image[Diff2D(1,2)] </TT>
        */
    reference operator[](difference_type const & d)
    {
        VIGRA_ASSERT_INSIDE(d);
        return data_[d.y*stride_ + d.x];
    }

        /** read pixel at given location. <br>
        usage: <TT> value_type value = image[Diff2D(1,2)] </TT>
        */
    const_reference operator[](difference_type const & d) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(d));
        return data_[d.y*stride_ + d.x];
    }

        /** access pixel at given location. <br>
        usage: <TT> value_type value = image(1,2) </TT>
        */
    reference operator()(std::ptrdiff_t dx, std::ptrdiff_t dy)
    {
        VIGRA_ASSERT_INSIDE(difference_type(dx,dy));
        return data_[dy*stride_ + dx];
    }

        /** read pixel at given location. <br>
        usage: <TT> value_type value = image(1,2) </TT>
        */
    const_reference operator()(std::ptrdiff_t dx, std::ptrdiff_t dy) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(dx, dy));
        return data_[dy*stride_ + dx];
    }

        /** access pixel at given location.
            Note that the 'x' index is the trailing index. <br>
        usage: <TT> value_type value = image[2][1] </TT>
        */
    pointer operator[](std::ptrdiff_t dy)
    {
        VIGRA_ASSERT_INSIDE(difference_type(0, dy));
        return data_ + dy*stride_;
    }

        /** read pixel at given location.
            Note that the 'x' index is the trailing index. <br>
        usage: <TT> value_type value = image[2][1] </TT>
        */
    const_pointer operator[](std::ptrdiff_t dy) const
    {
        VIGRA_ASSERT_INSIDE(difference_type(0,dy));
        return data_ + dy*stride_;
    }

        /** init 2D random access iterator pointing to upper left pixel
        */
    traverser upperLeft()
    {
        return traverser(data_, stride_);
    }

        /** init 2D random access iterator pointing to
         pixel(width, height), i.e. one pixel right and below lower right
         corner of the image as is common in C/C++.
        */
    traverser lowerRight()
    {
        return upperLeft() + size();
    }

        /** init 2D random access const iterator pointing to upper left pixel
        */
    const_traverser upperLeft() const
    {
        return const_traverser(data_, stride_);
    }

        /** init 2D random access const iterator pointing to
         pixel(width, height), i.e. one pixel right and below lower right
         corner of the image as is common in C/C++.
        */
    const_traverser lowerRight() const
    {
        return upperLeft() + size();
    }

        /** init 1D random access iterator pointing to first pixel.
            Note: Only works if stride equals width.
        */
    iterator begin()
    {
        vigra_precondition(stride_ == width_,
            "BasicImageView::begin(): "
            "can only create scan order iterator if width() == stride().");
        return data_;
    }

        /** init 1D random access iterator pointing past the end.
            Note: Only works if stride equals width.
        */
    iterator end()
    {
        vigra_precondition(stride_ == width_,
            "BasicImageView::end(): "
            "can only create scan order iterator if width() == stride().");
        return data_ + width() * height();
    }

        /** init 1D random access const iterator pointing to first pixel.
            Note: Only works if stride equals width.
        */
    const_iterator begin() const
    {
        vigra_precondition(stride_ == width_,
            "BasicImageView::begin(): "
            "can only create scan order iterator if width() == stride().");
        return data_;
    }

        /** init 1D random access const iterator pointing past the end.
            Note: Only works if stride equals width.
        */
    const_iterator end() const
    {
        vigra_precondition(stride_ == width_,
            "BasicImageView::end(): "
            "can only create scan order iterator if width() == stride().");
        return data_ + width() * height();
    }

        /** init 1D random access iterator pointing to first pixel of row \a y
        */
    row_iterator rowBegin(std::ptrdiff_t y)
    {
        return data_ + stride_ * y;
    }

        /** init 1D random access iterator pointing past the end of row \a y
        */
    row_iterator rowEnd(std::ptrdiff_t y)
    {
        return rowBegin(y) + width();
    }

        /** init 1D random access const iterator pointing to first pixel of row \a y
        */
    const_row_iterator rowBegin(std::ptrdiff_t y) const
    {
        return data_ + stride_ * y;
    }

        /** init 1D random access const iterator pointing past the end of row \a y
        */
    const_row_iterator rowEnd(std::ptrdiff_t y) const
    {
        return rowBegin(y) + width();
    }

        /** init 1D random access iterator pointing to first pixel of column \a x
        */
    column_iterator columnBegin(std::ptrdiff_t x)
    {
        typedef typename column_iterator::BaseType Iter;
        return column_iterator(Iter(data_ + x, stride_));
    }

        /** init 1D random access iterator pointing past the end of column \a x
        */
    column_iterator columnEnd(std::ptrdiff_t x)
    {
        return columnBegin(x) + height();
    }

        /** init 1D random access const iterator pointing to first pixel of column \a x
        */
    const_column_iterator columnBegin(std::ptrdiff_t x) const
    {
        typedef typename const_column_iterator::BaseType Iter;
        return const_column_iterator(Iter(data_ + x, stride_));
    }

        /** init 1D random access const iterator pointing past the end of column \a x
        */
    const_column_iterator columnEnd(std::ptrdiff_t x) const
    {
        return columnBegin(x) + height();
    }

        /** get a pointer to the internal data
        */
    const_pointer data() const
    {
        return data_;
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

    pointer data_;
    std::ptrdiff_t width_, height_, stride_;
};


/********************************************************/
/*                                                      */
/*              argument object factories               */
/*                                                      */
/********************************************************/

template <class PixelType, class Accessor>
inline triple<typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::const_traverser, Accessor>
srcImageRange(BasicImageView<PixelType> const & img, Accessor a)
{
    return triple<typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::const_traverser,
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline triple<typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::const_traverser, Accessor>
srcImageRange(BasicImageView<PixelType> const & img, Rect2D const & roi, Accessor a)
{
    vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "srcImageRange(): ROI rectangle outside image.");
    return triple<typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::const_traverser,
          Accessor>(img.upperLeft() + roi.upperLeft(),
                    img.upperLeft() + roi.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
srcImage(BasicImageView<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::const_traverser,
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
srcImage(BasicImageView<PixelType> const & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "srcImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::const_traverser,
                Accessor>(img.upperLeft() + ul, a);
}

template <class PixelType, class Accessor>
inline triple<typename BasicImageView<PixelType>::traverser,
              typename BasicImageView<PixelType>::traverser, Accessor>
destImageRange(BasicImageView<PixelType> & img, Accessor a)
{
    return triple<typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::traverser,
          Accessor>(img.upperLeft(),
                    img.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline triple<typename BasicImageView<PixelType>::traverser,
              typename BasicImageView<PixelType>::traverser, Accessor>
destImageRange(BasicImageView<PixelType> & img, Rect2D const & roi, Accessor a)
{
    vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "destImageRange(): ROI rectangle outside image.");
    return triple<typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::traverser,
          Accessor>(img.upperLeft() + roi.upperLeft(),
                    img.upperLeft() + roi.lowerRight(),
                    a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::traverser, Accessor>
destImage(BasicImageView<PixelType> & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::traverser,
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::traverser, Accessor>
destImage(BasicImageView<PixelType> & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "destImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::traverser,
                Accessor>(img.upperLeft() + ul, a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
maskImage(BasicImageView<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::const_traverser,
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
maskImage(BasicImageView<PixelType> const & img, Point2D const & ul, Accessor a)
{
    vigra_precondition(img.isInside(ul),
                       "maskImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::const_traverser,
                Accessor>(img.upperLeft() + ul, a);
}

/****************************************************************/

template <class PixelType>
inline triple<typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::ConstAccessor>
srcImageRange(BasicImageView<PixelType> const & img)
{
    return triple<typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft(),
                                                                        img.lowerRight(),
                                                                        img.accessor());
}

template <class PixelType>
inline triple<typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::const_traverser,
              typename BasicImageView<PixelType>::ConstAccessor>
srcImageRange(BasicImageView<PixelType> const & img, Rect2D const & roi)
{
    vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "srcImageRange(): ROI rectangle outside image.");
    return triple<typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::const_traverser,
                  typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft() + roi.upperLeft(),
                                                                        img.upperLeft() + roi.lowerRight(),
                                                                        img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::const_traverser,
             typename BasicImageView<PixelType>::ConstAccessor>
srcImage(BasicImageView<PixelType> const & img)
{
    return pair<typename BasicImageView<PixelType>::const_traverser,
                typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft(),
                                                                      img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::const_traverser,
             typename BasicImageView<PixelType>::ConstAccessor>
srcImage(BasicImageView<PixelType> const & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "srcImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::const_traverser,
                typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft() + ul,
                                                                      img.accessor());
}

template <class PixelType>
inline triple< typename BasicImageView<PixelType>::traverser,
               typename BasicImageView<PixelType>::traverser,
               typename BasicImageView<PixelType>::Accessor>
destImageRange(BasicImageView<PixelType> & img)
{
    return triple<typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::Accessor>(img.upperLeft(),
                                                                   img.lowerRight(),
                                                                   img.accessor());
}

template <class PixelType>
inline triple< typename BasicImageView<PixelType>::traverser,
               typename BasicImageView<PixelType>::traverser,
               typename BasicImageView<PixelType>::Accessor>
destImageRange(BasicImageView<PixelType> & img, Rect2D const & roi)
{
    vigra_precondition(roi.left() >= 0 && roi.top() >= 0 &&
                       roi.right() <= img.width() && roi.bottom() <= img.height(),
                       "destImageRange(): ROI rectangle outside image.");
    return triple<typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::traverser,
                  typename BasicImageView<PixelType>::Accessor>(img.upperLeft() + roi.upperLeft(),
                                                                   img.upperLeft() + roi.lowerRight(),
                                                                   img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::traverser,
             typename BasicImageView<PixelType>::Accessor>
destImage(BasicImageView<PixelType> & img)
{
    return pair<typename BasicImageView<PixelType>::traverser,
                typename BasicImageView<PixelType>::Accessor>(img.upperLeft(),
                                                          img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::traverser,
             typename BasicImageView<PixelType>::Accessor>
destImage(BasicImageView<PixelType> & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "destImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::traverser,
                typename BasicImageView<PixelType>::Accessor>(img.upperLeft() + ul,
                                                                 img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::const_traverser,
             typename BasicImageView<PixelType>::ConstAccessor>
maskImage(BasicImageView<PixelType> const & img)
{
    return pair<typename BasicImageView<PixelType>::const_traverser,
                typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft(),
                                                                      img.accessor());
}

template <class PixelType>
inline pair< typename BasicImageView<PixelType>::const_traverser,
             typename BasicImageView<PixelType>::ConstAccessor>
maskImage(BasicImageView<PixelType> const & img, Point2D const & ul)
{
    vigra_precondition(img.isInside(ul),
                       "maskImage(): ROI rectangle outside image.");
    return pair<typename BasicImageView<PixelType>::const_traverser,
                typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft() + ul,
                                                                      img.accessor());
}

} // namespace vigra
#undef VIGRA_ASSERT_INSIDE
#endif /* VIGRA_BASICIMAGEVIEW_HXX */
