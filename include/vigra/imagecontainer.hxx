#ifndef IMAGECONTAINER_HXX // -*- indent-tabs-mode: nil -*-
#define IMAGECONTAINER_HXX

#include "vigra/utilities.hxx"
#include <vector>

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
    template parameter. Use it like a std::vector<ImageType>, it has
    the same interface, only operator< is missing from ImageArray. It
    offers additional functions for resizing the images and querying
    their common size. See \ref imageSize() for additional notes.

    <b>\#include</b> "<a href="imagecontainer_8hxx-source.html">vigra/imagecontainer.hxx</a>"

    Namespace: vigra
*/
template <class ImageType>
class ImageArray
{
    Diff2D imageSize_;

protected:
    std::vector<ImageType> images_;

public:
        /** the type of the contained values/images
         */
    typedef ImageType    value_type;

    typedef typename std::vector<ImageType>::iterator iterator;
    typedef typename std::vector<ImageType>::const_iterator const_iterator;
    typedef typename std::vector<ImageType>::reverse_iterator reverse_iterator;
    typedef typename std::vector<ImageType>::const_reverse_iterator const_reverse_iterator;
    typedef typename std::vector<ImageType>::reference reference;
    typedef typename std::vector<ImageType>::const_reference const_reference;
    typedef typename std::vector<ImageType>::pointer pointer;
    typedef typename std::vector<ImageType>::difference_type difference_type;
    typedef typename std::vector<ImageType>::size_type size_type;

        /** init an array of numImages equal-sized images
         */
    ImageArray(unsigned int numImages, const Diff2D &imageSize)
        : imageSize_(imageSize),
          images_(numImages)
    {
        for(unsigned int i=0; i<numImages; i++)
            images_[i].resize(imageSize);
    }

        /** Init an array of numImages equal-sized images, the size
            depends on ImageType's default constructor.
         */
    ImageArray(unsigned int numImages= 0)
        : images_(numImages)
    {
        imageSize_= empty()? Diff2D(0, 0) : front().size();
    }

        /** fill constructor: Init an array with numImages copies of
            the given image. (STL-Sequence interface)
         */
    ImageArray(unsigned int numImages, const ImageType &image)
        : imageSize_(image.size()),
          images_(numImages, image)
    {
    }

        /** range constructor: Construct an array containing copies of
            the images in [begin, end). Those images must all have the
            same size, see \ref imageSize(). (STL-Sequence interface)
         */
    template<class InputIterator>
    ImageArray(InputIterator begin, InputIterator end)
        : imageSize_(begin!=end? (*begin).size() : Diff2D(0,0)),
          images_(begin, end)
    {
    }

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
        Diff2D oldImageSize = imageSize_;
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
    Diff2D imageSize() const
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
                images_[i].resize(newSize);
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
        resizeImages(Diff2D(width, height));
    }
};

//@}

} // namespace vigra

#endif // IMAGECONTAINER_HXX
