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

#ifndef VIGRA_BASICIMAGEVIEW_HXX
#define VIGRA_BASICIMAGEVIEW_HXX

#include "vigra/imageiterator.hxx"
#include "vigra/initimage.hxx"

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
    
    <b>\#include</b> "<a href="basicimageview_8hxx-source.html">vigra/basicimageview.hxx</a>"

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
    BasicImageView(const_pointer data, int w, int h, int stride = 0)
    : data_(const_cast<pointer>(data)),
      width_(w),
      height_(h),
      stride_(stride == 0 ? w : stride)
    {}

        /** construct view of size size.x x size.y
        */
    BasicImageView(const_pointer data, difference_type const & size, int stride = 0)
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

        /** stride of Image. 
            Memory offset between the start of two successive rows.
        */
    int stride() const
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
        return data_[d.y*stride_ + d.x];
    }

        /** read pixel at given location. <br>
	    usage: <TT> value_type value = image[Diff2D(1,2)] </TT>
        */
    const_reference operator[](difference_type const & d) const
    {
        return data_[d.y*stride_ + d.x];
    }

        /** access pixel at given location. <br>
	    usage: <TT> value_type value = image(1,2) </TT>
        */
    reference operator()(int dx, int dy)
    {
        return data_[dy*stride_ + dx];
    }

        /** read pixel at given location. <br>
	    usage: <TT> value_type value = image(1,2) </TT>
        */
    const_reference operator()(int dx, int dy) const
    {
        return data_[dy*stride_ + dx];
    }

        /** access pixel at given location.
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> value_type value = image[2][1] </TT>
        */
    pointer operator[](int dy)
    {
        return data_ + dy*stride_;
    }

        /** read pixel at given location.
	        Note that the 'x' index is the trailing index. <br>
	    usage: <TT> value_type value = image[2][1] </TT>
        */
    const_pointer operator[](int dy) const
    {
        return data_ + dy*stride_;
    }

        /** init 2D random access iterator poining to upper left pixel
        */
    traverser upperLeft()
    {
        return traverser(data_, stride_);
    }

        /** init 2D random access iterator poining to
         pixel(width, height), i.e. one pixel right and below lower right
         corner of the image as is common in C/C++.
        */
    traverser lowerRight()
    {
        return upperLeft() + size();
    }

        /** init 2D random access const iterator poining to upper left pixel
        */
    const_traverser upperLeft() const
    {
        return const_traverser(data_, stride_);
    }

        /** init 2D random access const iterator poining to
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
    int width_, height_, stride_;
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
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
srcImage(BasicImageView<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
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
inline pair<typename BasicImageView<PixelType>::traverser, Accessor>
destImage(BasicImageView<PixelType> & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::traverser, 
                Accessor>(img.upperLeft(), a);
}

template <class PixelType, class Accessor>
inline pair<typename BasicImageView<PixelType>::const_traverser, Accessor>
maskImage(BasicImageView<PixelType> const & img, Accessor a)
{
    return pair<typename BasicImageView<PixelType>::const_traverser, 
                Accessor>(img.upperLeft(), a);
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
inline pair< typename BasicImageView<PixelType>::const_traverser, 
             typename BasicImageView<PixelType>::ConstAccessor>
srcImage(BasicImageView<PixelType> const & img)
{
    return pair<typename BasicImageView<PixelType>::const_traverser, 
                typename BasicImageView<PixelType>::ConstAccessor>(img.upperLeft(), 
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
inline pair< typename BasicImageView<PixelType>::traverser, 
             typename BasicImageView<PixelType>::Accessor>
destImage(BasicImageView<PixelType> & img)
{
    return pair<typename BasicImageView<PixelType>::traverser, 
                typename BasicImageView<PixelType>::Accessor>(img.upperLeft(), 
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

} // namespace vigra

#endif /* VIGRA_BASICIMAGEVIEW_HXX */
