/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#ifndef VIGRA_IMAGECONTAINER_HXX
#define VIGRA_IMAGECONTAINER_HXX

#include "utilities.hxx"
#include "array_vector.hxx"
#include "copyimage.hxx"

namespace vigra {

/** \addtogroup ImageContainers Image Containers
    Classes to manage multiple images (ImageArray..)
*/
//@{

/********************************************************/
/*                                                      */
/*                      ImageArray                      */
/*                                                      */
/********************************************************/

/** \brief Fundamental class template for arrays of equal-sized images.

    An ImageArray manages an array of images of the type given as
    template parameter. Use it like a ArrayVector<ImageType>, it has
    the same interface, only operator< is missing from ImageArray. It
    offers additional functions for resizing the images and querying
    their common size. See \ref imageSize() for additional notes.

    A custimized allocator can be passed as a template argument and via the constructor.
    By default, the allocator of the <tt>ImageType</tt> is reused.

    <b>\#include</b> \<<a href="imagecontainer_8hxx-source.html">vigra/imagecontainer.hxx</a>\>

    Namespace: vigra
*/
template <class ImageType,
      class Alloc = typename ImageType::allocator_type::template rebind<ImageType>::other >
class ImageArray
{
    Size2D imageSize_;

protected:
    typedef ArrayVector<ImageType, Alloc> ImageVector;
    ImageVector images_;

public:
        /** the type of the contained values/images
         */
    typedef ImageType    value_type;

    typedef typename ImageVector::iterator iterator;
    typedef typename ImageVector::const_iterator const_iterator;
    typedef typename ImageVector::reverse_iterator reverse_iterator;
    typedef typename ImageVector::const_reverse_iterator const_reverse_iterator;
    typedef typename ImageVector::reference reference;
    typedef typename ImageVector::const_reference const_reference;
#if !defined(_MSC_VER) || _MSC_VER >= 1300
    typedef typename ImageVector::pointer pointer;
#endif
    typedef typename ImageVector::difference_type difference_type;
    typedef typename ImageVector::size_type size_type;

        /** init an array of numImages equal-sized images; use the specified allocator.
         */
    ImageArray(unsigned int numImages, const Diff2D &imageSize,
               Alloc const & alloc = Alloc())
        : imageSize_(imageSize),
          images_(numImages, ImageType(), alloc)
    {
        for(unsigned int i=0; i<numImages; i++)
            images_[i].resize(Size2D(imageSize));
    }

        /** Init an array of numImages equal-sized images. The size
            depends on ImageType's default constructor (so it will
            usually be 0x0); use the specified allocator.
         */
    ImageArray(unsigned int numImages= 0, Alloc const & alloc = Alloc())
        : images_(numImages, alloc)
    {
        imageSize_= empty()? Size2D(0, 0) : front().size();
    }

        /** fill constructor: Init an array with numImages copies of
            the given image. (STL-Sequence interface); use the specified allocator.
         */
    ImageArray(unsigned int numImages, const ImageType &image, Alloc const & alloc = Alloc())
        : imageSize_(image.size()),
          images_(numImages, image, alloc)
    {
    }

        /** range constructor: Construct an array containing copies of
            the images in [begin, end). Those images must all have the
            same size, see \ref imageSize(). (STL-Sequence interface);
            use the specified allocator.
         */
    template<class InputIterator>
    ImageArray(InputIterator begin, InputIterator end, Alloc const & alloc = Alloc())
        : imageSize_(begin!=end? (*begin).size() : Size2D(0,0)),
          images_(begin, end, alloc)
    {
    }

    virtual ~ImageArray() {}

        /** Operator for a vector-like access to the contained images
            (STL-Vector interface)
         */
    reference operator [](size_type index)
    {
        return images_[index];
    }

        /** Operator for a vector-like access to the contained images
            (STL-Vector interface)
         */
    const_reference operator [](size_type index) const
    {
        return images_[index];
    }

        /** Returns an iterator pointing to the first image
            (STL-Container interface)
         */
    iterator begin()
    {
        return images_.begin();
    }

        /** Returns an iterator pointing to the first image
            (STL-Container interface)
         */
    const_iterator begin() const
    {
        return images_.begin();
    }

        /** Returns an iterator pointing behind the last image
            (STL-Container interface)
         */
    iterator end()
    {
        return images_.end();
    }

        /** Returns an iterator pointing behind the last image
            (STL-Container interface)
         */
    const_iterator end() const
    {
        return images_.end();
    }

        /** Returns a reverse_iterator pointing to the first image of
            the reversed view of this array (STL-Reversable Container
            interface)
         */
    reverse_iterator rbegin()
    {
        return images_.rbegin();
    }

        /** Returns a reverse_iterator pointing to the first image of
            the reversed view of this array (STL-Reversable Container
            interface)
         */
    const_reverse_iterator rbegin() const
    {
        return images_.rbegin();
    }

        /** Returns a reverse_iterator pointing behind the last image
            of the reversed view of this array (STL-Reversable
            Container interface)
         */
    reverse_iterator rend()
    {
        return images_.rend();
    }

        /** Returns a reverse_iterator pointing behind the last image
            of the reversed view of this array (STL-Reversable
            Container interface)
         */
    const_reverse_iterator rend() const
    {
        return images_.rend();
    }

        /** Query size of this ImageArray, that is: the number of
            images. (STL-Container interface)
        */
    size_type size() const
    {
        return images_.size();
    }

        /** Query maximum size of this ImageArray, that is: the
            max. parameter you may pass to resize(). (STL-Container
            interface)
        */
    size_type max_size() const
    {
        return images_.max_size();
    }

        /** Returns true if and only if there are no contained
            images. (STL-Container interface)
        */
    bool empty()
    {
        return images_.empty();
    }

        /** Returns true if and only if both ImageArrays have exactly
            the same contents and all images did compare equal with the
            corresponding image in the other ImageArray. (STL-Forward
            Container interface)
         */
    bool operator ==(const ImageArray<ImageType> &other)
    {
        return (imageSize() == other.imageSize())
                && (images_ == other.images_);
    }

        /** Insert image at/before pos. (STL-Sequence interface)
         */
    iterator insert(iterator pos, const_reference image)
    {
        return images_.insert(pos, image);
    }

        /** Insert count copies of image at/before pos. (STL-Sequence
            interface)
         */
    void insert (iterator pos, size_type count, const_reference image);

        /** Insert copies of images from [begin, end) at/before
            pos. (STL-Sequence interface)
         */
    template<class InputIterator>
    void insert(iterator pos, InputIterator begin, InputIterator end)
    {
        images_.insert(pos, begin, end);
    }

        /** Removes the image at pos from this array. (STL-Sequence
            interface)
         */
    iterator erase(iterator pos)
    {
        return images_.erase(pos);
    }

        /** Removes the images from [begin, end) from this
            array. (STL-Sequence interface)
         */
    iterator erase(iterator begin, iterator end)
    {
        return images_.erase(begin, end);
    }

        /** Empty this array. (STL-Sequence interface)
         */
    void clear()
    {
        images_.clear();
    }

        /** Resize this ImageArray, throwing the last images away if
            you make the array smaller or appending new images of the
            right size at the end of the array if you make it
            larger. (STL-Sequence interface)
        */
    void resize(size_type newSize)
    {
        if (newSize != size())
        {
            size_type oldSize= size();
            images_.resize(newSize);
            for (size_type i= oldSize; i<newSize; i++)
                images_[i].resize(imageSize());
        }
    }

        /** Resize this ImageArray, throwing the last images away if
            you make the array smaller or appending new copies of image
            at the end of the array if you make it larger.
            precondition: <tt>image.size() == imageSize()</tt>
            (STL-Sequence interface)
        */
    void resize(size_type newSize, ImageType &image)
    {
        if (newSize != size())
        {
            vigra_precondition(image.size() == imageSize(),
                               "trying to append images of wrong size to ImageArray with resize()");
            images_.resize(newSize, image);
        }
    }

        /** return the first image. (STL-Sequence interface)
         */
    reference front()
    {
        return images_.front();
    }

        /** return the first image. (STL-Sequence interface)
         */
    const_reference front() const
    {
        return images_.front();
    }

        /** return the last image. (STL-Vector interface)
         */
    reference back()
    {
        return images_.back();
    }

        /** return the last image. (STL-Vector interface)
         */
    const_reference back() const
    {
        return images_.back();
    }

        /** append image to array (STL-Back Insertion Sequence interface)
         */
    void push_back(const_reference image)
    {
        images_.push_back(image);
    }

        /** remove last image from array (STL-Back Insertion Sequence interface)
         */
    void pop_back()
    {
        images_.pop_back();
    }

        /** swap contents of this array with the contents of other
            (STL-Container interface)
         */
    void swap(const_reference other)
    {
        Size2D oldImageSize = imageSize_;
        images_.swap(other.images_);
        imageSize_ = other.imageSize_;
        other.imageSize_ = oldImageSize;
    }

        /** number of image objects for which memory has been allocated
            (STL-Vector interface)
        */
    size_type capacity() const
    {
        return images_.capacity();
    }

        /** increase capacity(). (STL-Vector interface)
         */
    void reserve(size_type n)
    {
        images_.reserve(n);
    }

        /** Query the size of the contained images. ImageArray will
            maintain an array of equal-sized images of this
            size. However, <em>do not resize the contained images
            manually</em>. ImageArray currently has no way to detect or
            prevent this.
         */
    Size2D imageSize() const
        { return imageSize_; }

        /** Resize all images to a common new size (No-op if
            <tt>newSize == imageSize()</tt>). See \ref imageSize() for
            an important note about resizing the images.
        */
    virtual void resizeImages(const Diff2D &newSize)
    {
        if (newSize!=imageSize())
        {
            for(unsigned int i=0; i<size(); i++)
                images_[i].resize(Size2D(newSize));
            imageSize_= newSize;
        }
    }

        /** Resize all images to a common new size (No-op if
            <tt>newSize == imageSize()</tt>). See \ref imageSize() for
            an important note about resizing the images.

            (Convenience function, same as calling
            <tt>resizeImages(Diff2D(width, height));</tt>.)
        */
    void resizeImages(int width, int height)
    {
        resizeImages(Size2D(width, height));
    }
};

/********************************************************/
/*                                                      */
/*                      ImagePyramid                    */
/*                                                      */
/********************************************************/

/** \brief Class template for logarithmically tapering image pyramids.

    An ImagePyramid manages an array of images of the type given as
    template parameter, where each level has half the width and height
    of its predecessor.  It actually represents a sequence of pyramid
    levels whose start and end index are configurable.  For Burt-style
    pyramids, see also \ref pyramidReduceBurtFilter and \ref
    pyramidExpandBurtFilter.

    A custimized allocator can be passed as a template argument and
    via the constructor.  By default, the allocator of the
    <tt>ImageType</tt> is reused.

    <b>\#include</b> \<<a href="imagecontainer_8hxx-source.html">vigra/imagecontainer.hxx</a>\>

    Namespace: vigra
*/
template <class ImageType,
      class Alloc = typename ImageType::allocator_type::template rebind<ImageType>::other >
class ImagePyramid
{
    int lowestLevel_, highestLevel_;

protected:
    typedef ArrayVector<ImageType, Alloc> ImageVector;
    ImageVector images_;

public:
        /** the type of the contained values/images
         */
    typedef ImageType    value_type;

    typedef typename ImageVector::iterator iterator;
    typedef typename ImageVector::const_iterator const_iterator;
    typedef typename ImageVector::reverse_iterator reverse_iterator;
    typedef typename ImageVector::const_reverse_iterator const_reverse_iterator;
    typedef typename ImageVector::reference reference;
    typedef typename ImageVector::const_reference const_reference;
#if !defined(_MSC_VER) || _MSC_VER >= 1300
    typedef typename ImageVector::pointer pointer;
#endif
    typedef typename ImageVector::difference_type difference_type;
    typedef int size_type;

        /** Init a pyramid between the given levels (inclusive).
         *
         * Allocate the given \a imageSize at the pyramid level given
         * in \a sizeAppliesToLevel (default: level 0 / bottom) and
         * size the other levels using recursive reduction/expansion
         * by factors of 2.  Use the specified allocator for image
         * creation.  The image type must be default constructible and
         * resizable.  sizeAppliesToLevel must be the in range
         * lowestLevel..highestLevel (inclusive).
         */
    ImagePyramid(int lowestLevel, int highestLevel,
                 const Diff2D &imageSize, int sizeAppliesToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, imageSize, sizeAppliesToLevel);
    }

        /**
         * Init a pyramid between the given levels (inclusive).
         *
         * Copy the given \a image into the pyramid level given in \a
         * copyImageToLevel (default: level 0 / bottom) and size the
         * other levels using recursive reduction/expansion by factors
         * of 2 (their image data is not initialized).  Use the
         * specified allocator for image creation.  The image type
         * must be default constructible and resizable.
         * sizeAppliesToLevel must be the in range
         * lowestLevel..highestLevel (inclusive).
         */
    ImagePyramid(int lowestLevel, int highestLevel,
                 const ImageType &image, int copyImageToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, image.size(), copyImageToLevel);
        copyImage(srcImageRange(image), destImage((*this)[copyImageToLevel]));
    }

        /**
         * Init a pyramid between the given levels (inclusive).
         *
         * Copy the image given by the range \a ul to \a lr into the
         * pyramid level given in \a copyImageToLevel (default: level
         * 0 / bottom) and size the other levels using recursive
         * reduction/expansion by factors of 2 (their image data is
         * not initialized).  Use the specified allocator for image
         * creation.  The image type must be default constructible and
         * resizable.  sizeAppliesToLevel must be the in range
         * lowestLevel..highestLevel (inclusive).
         */
    template <class SrcIterator, class SrcAccessor>
    ImagePyramid(int lowestLevel, int highestLevel,
                 SrcIterator ul, SrcIterator lr, SrcAccessor src,
                 int copyImageToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, lr - ul, copyImageToLevel);
        copyImage(srcIterRange(ul, lr, src), destImage((*this)[copyImageToLevel]));
    }

        /** Init an empty pyramid.  Use the specified allocator.
         */
    ImagePyramid(Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {}

    virtual ~ImagePyramid() {}

        /** Get the index of the lowest allocated level of the pyramid.
        */
    int lowestLevel() const
    {
        return lowestLevel_;
    }

        /** Get the index of the highest allocated level of the pyramid.
        */
    int highestLevel() const
    {
        return highestLevel_;
    }

        /** Operator for a vector-like access to the contained images
            (STL-Vector interface)
         */
    reference operator [](size_type index)
    {
        return images_[index - lowestLevel_];
    }

        /** Operator for a vector-like access to the contained images
            (STL-Vector interface)
         */
    const_reference operator [](size_type index) const
    {
        return images_[index - lowestLevel_];
    }

        /** Returns an iterator pointing to the first image
            (STL-Container interface)
         */
    iterator begin()
    {
        return images_.begin();
    }

        /** Returns an iterator pointing to the first image
            (STL-Container interface)
         */
    const_iterator begin() const
    {
        return images_.begin();
    }

        /** Returns an iterator pointing behind the last image
            (STL-Container interface)
         */
    iterator end()
    {
        return images_.end();
    }

        /** Returns an iterator pointing behind the last image
            (STL-Container interface)
         */
    const_iterator end() const
    {
        return images_.end();
    }

        /** Returns a reverse_iterator pointing to the first image of
            the reversed view of this array (STL-Reversable Container
            interface)
         */
    reverse_iterator rbegin()
    {
        return images_.rbegin();
    }

        /** Returns a reverse_iterator pointing to the first image of
            the reversed view of this array (STL-Reversable Container
            interface)
         */
    const_reverse_iterator rbegin() const
    {
        return images_.rbegin();
    }

        /** Returns a reverse_iterator pointing behind the last image
            of the reversed view of this array (STL-Reversable
            Container interface)
         */
    reverse_iterator rend()
    {
        return images_.rend();
    }

        /** Returns a reverse_iterator pointing behind the last image
            of the reversed view of this array (STL-Reversable
            Container interface)
         */
    const_reverse_iterator rend() const
    {
        return images_.rend();
    }

        /** Query size of this ImageArray, that is: the number of
            images. (STL-Container interface)
        */
    size_type size() const
    {
        return images_.size();
    }

        /** Returns true if and only if there are no contained
            images. (STL-Container interface)
        */
    bool empty()
    {
        return images_.empty();
    }

        /** Returns true if and only if both ImageArrays have exactly
            the same contents and all images did compare equal with the
            corresponding image in the other ImageArray. (STL-Forward
            Container interface)
         */
    bool operator ==(const ImagePyramid<ImageType, Alloc> &other) const
    {
        return (lowestLevel_ == other.lowestLevel_) && (highestLevel_ == other.highestLevel_) &&
                (images_ == other.images_);
    }

        /** Empty this array. (STL-Sequence interface)
         */
    void clear()
    {
        images_.clear();
        lowestLevel_ = 0;
        highestLevel_ = -1;
    }

        /** Resize this ImageArray, throwing the last images away if
            you make the array smaller or appending new images of the
            right size at the end of the array if you make it
            larger. (STL-Sequence interface)
        */
    void resize(int lowestLevel, int highestLevel,
                const Diff2D &imageSize, int sizeAppliesToLevel = 0)
    {
        vigra_precondition(lowestLevel <= highestLevel,
           "ImagePyramid::resize(): lowestLevel <= highestLevel required.");
        vigra_precondition(lowestLevel <= sizeAppliesToLevel && sizeAppliesToLevel <= highestLevel,
           "ImagePyramid::resize(): sizeAppliesToLevel must be between lowest and highest level (inclusive).");

        ImageVector images(highestLevel - lowestLevel + 1, ImageType());

        images[sizeAppliesToLevel - lowestLevel].resize(imageSize);
        for(int i=sizeAppliesToLevel + 1; i<=highestLevel; ++i)
        {
            unsigned int w = (images[i - 1 - lowestLevel].width() + 1) / 2;
            unsigned int h = (images[i - 1 - lowestLevel].height() + 1) / 2;
            images[i - lowestLevel].resize(w, h);
        }
        for(int i=sizeAppliesToLevel - 1; i>=lowestLevel; --i)
        {
            unsigned int w = 2*images[i + 1 - lowestLevel].width() - 1;
            unsigned int h = 2*images[i + 1 - lowestLevel].height() - 1;
            images[i - lowestLevel].resize(w, h);
        }

        images_.swap(images);
        lowestLevel_ = lowestLevel;
        highestLevel_ = highestLevel;
    }

        /** return the first image (lowestLevel()). (STL-Sequence interface)
         */
    reference front()
    {
        return images_.front();
    }

        /** return the first image (lowestLevel()). (STL-Sequence interface)
         */
    const_reference front() const
    {
        return images_.front();
    }

        /** return the last image (highestLevel()). (STL-Vector interface)
         */
    reference back()
    {
        return images_.back();
    }

        /** return the last image (highestLevel()). (STL-Vector interface)
         */
    const_reference back() const
    {
        return images_.back();
    }

        /** swap contents of this array with the contents of other
            (STL-Container interface)
         */
    void swap(const ImagePyramid<ImageType, Alloc> &other)
    {
        images_.swap(other.images_);
        std::swap(lowestLevel_, other.lowestLevel_);
        std::swap(highestLevel_, other.highestLevel_);
    }
};

//@}

} // namespace vigra

#endif // VIGRA_IMAGECONTAINER_HXX
