/************************************************************************/
/*                                                                      */
/*         Copyright 1998-2016 by Sven Peter and Ullrich Koethe         */
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

#ifndef VIGRA_MULTIPYRAMID_HXX
#define VIGRA_MULTIPYRAMID_HXX

#include "utilities.hxx"
#include "array_vector.hxx"
#include "copyimage.hxx"
#include "multi_array.hxx"
#include "multi_resampling_convolution.hxx"

namespace vigra {

/** \brief Class template for logarithmically tapering image pyramids.

    An MultiArrayPyramid manages an array of multi arrays of the type given as
    template parameter, where each level has half the width and height
    of its predecessor.  It actually represents a sequence of pyramid
    levels whose start and end index are configurable.  
    
    To initialize all pyramid levels in the sense of a Gaussian pyramid,
    use \ref multiPyramidReduceBurtFilter() and \ref multiPyramidExpandBurtFilter().
    To create and reconstruct a Laplcaian pyramid, use 
    \ref multiPyramidReduceBurtLaplacian() and \ref multiPyramidExpandBurtLaplacian().

    A customized allocator can be passed as a template argument and
    via the constructor.  By default, the allocator of the
    <tt>MultiArrayType</tt> is reused.

    <b>\#include</b> \<vigra/imagecontainer.hxx\> <br/>
    Namespace: vigra
*/
template <class MultiArrayType,
      class Alloc = typename MultiArrayType::allocator_type::template rebind<MultiArrayType>::other >
class MultiArrayPyramid
{
    int lowestLevel_, highestLevel_;

protected:
    typedef ArrayVector<MultiArrayType, Alloc> MultiArrayVector;
    MultiArrayVector images_;

public:
        /** the type of the contained values/images
         */
    typedef MultiArrayType    value_type;

    typedef typename MultiArrayType::difference_type difference_type;

    typedef typename MultiArrayVector::iterator iterator;
    typedef typename MultiArrayVector::const_iterator const_iterator;
    typedef typename MultiArrayVector::reverse_iterator reverse_iterator;
    typedef typename MultiArrayVector::const_reverse_iterator const_reverse_iterator;
    typedef typename MultiArrayVector::reference reference;
    typedef typename MultiArrayVector::const_reference const_reference;
#if !defined(_MSC_VER) || _MSC_VER >= 1300
    typedef typename MultiArrayVector::pointer pointer;
#endif
#if 0
    //typedef typename MultiArrayVector::difference_type difference_type;
#endif
    typedef int size_type;

        /** Init a pyramid between the given levels (inclusive).
         *
         * Allocate the given \a multiArrayShape at the pyramid level given
         * in \a sizeAppliesToLevel (default: level 0 / bottom) and
         * size the other levels using recursive reduction/expansion
         * by factors of 2.  Use the specified allocator for image
         * creation.  The image type must be default constructible and
         * resizable.  sizeAppliesToLevel must be the in range
         * lowestLevel..highestLevel (inclusive).
         */
    MultiArrayPyramid(int lowestLevel, int highestLevel,
                 difference_type multiArrayShape, int sizeAppliesToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, multiArrayShape, sizeAppliesToLevel);
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
    MultiArrayPyramid(int lowestLevel, int highestLevel,
                 const MultiArrayType &image, int copyImageToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, image.shape(), copyImageToLevel);
        copyMultiArray(srcMultiArrayRange(image), destMultiArray((*this)[copyImageToLevel]));
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
    MultiArrayPyramid(int lowestLevel, int highestLevel,
                 SrcIterator ul, SrcIterator lr, SrcAccessor src,
                 int copyImageToLevel = 0,
                 Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {
        resize(lowestLevel, highestLevel, lr - ul, copyImageToLevel);
        copyMultiArray(srcIterRange(ul, lr, src), destMultiArray((*this)[copyImageToLevel]));
    }

        /** Init an empty pyramid.  Use the specified allocator.
         */
    MultiArrayPyramid(Alloc const & alloc = Alloc())
        : lowestLevel_(0), highestLevel_(-1),
          images_(alloc)
    {}

    virtual ~MultiArrayPyramid() {}

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

        /** Query size of this MultiArrayPyramid, that is: the number of
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

        /** Returns true if and only if both MultiArrayPyramids have exactly
            the same contents and all images did compare equal with the
            corresponding image in the other MultiArrayPyramid. (STL-Forward
            Container interface)
         */
    bool operator ==(const MultiArrayPyramid<MultiArrayType, Alloc> &other) const
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

        /** Resize this MultiArrayPyramid, throwing the last images away if
            you make the array smaller or appending new images of the
            right size at the end of the array if you make it
            larger. (STL-Sequence interface)
        */
    void resize(int lowestLevel, int highestLevel,
                const difference_type multiArrayShape, int sizeAppliesToLevel = 0)
    {
        vigra_precondition(lowestLevel <= highestLevel,
           "MultiArrayPyramid::resize(): lowestLevel <= highestLevel required.");
        vigra_precondition(lowestLevel <= sizeAppliesToLevel && sizeAppliesToLevel <= highestLevel,
           "MultiArrayPyramid::resize(): sizeAppliesToLevel must be between lowest and highest level (inclusive).");


        MultiArrayVector images(highestLevel - lowestLevel + 1, MultiArrayType());

        images[sizeAppliesToLevel - lowestLevel].reshape(multiArrayShape);

        for(int i=sizeAppliesToLevel + 1; i<=highestLevel; ++i)
            images[i - lowestLevel].reshape(floor((images[i - 1 - lowestLevel].shape() + 1) / 2));
        for(int i=sizeAppliesToLevel - 1; i>=lowestLevel; --i)
            images[i - lowestLevel].reshape(floor(2*images[i + 1 - lowestLevel].shape() - 1));

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
    void swap(const MultiArrayPyramid<MultiArrayType, Alloc> &other)
    {
        images_.swap(other.images_);
        std::swap(lowestLevel_, other.lowestLevel_);
        std::swap(highestLevel_, other.highestLevel_);
    }
};

/** \brief Two-fold down-sampling for image pyramid construction.

    This function implements the reduction by one resolution level (first signature)
    or across several pyramid levels (last signature) of a Gaussian pyramid as described in

    P. Burt, E. Adelson: <i>"The Laplacian Pyramid as a Compact Image Code"</i>, IEEE Trans. Communications, 9(4):532–540, 1983

    That is, it applies the smoothing filter
    \code
    [0.25 - centerValue / 2.0, 0.25, centerValue, 0.25, 0.25 - centerValue / 2.0]
    \endcode
    to the source image and then copies all pixels with even coordinates to the destination
    image. The destination image shape must be <tt>dest_shape = ceil(src_shape / 2.0)</tt>.
    <tt>centerValue</tt> must be between 0.25 and 0.5 and determines the strength of smoothing
    (bigger values correspond to less smoothing). If <tt>toLevel - fromLevel > 1</tt> in the
    pyramid variant of the function, this process is repeated until <tt>toLevel</tt> is
    reached.

    Typically, this functions is used in connection with a \ref vigra::MultiArrayPyramid
    (last signature below) to perform several levels of reduction in one go.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/multi_pyramid.hxx\><br>
    Namespace: vigra

    pass multi array iterators, shapes and accessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestShape>
        void multiPyramidReduceBurtFilter(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                     DestIterator d, DestShape const & dshape, DestAccessor dest,
                                     double centerValue = 0.4)
    }
    \endcode

    pass multi array views directly:
    \code
    namespace vigra {
        template <unsigned int srcN, class srcT, class srcStrideTag,
                  unsigned int destN, class destT, class destStrideTag>
        void multiPyramidReduceBurtFilter(MultiArrayView<srcN, srcT, srcStrideTag> src,
                                     MultiArrayView<destN, destT, destStrideTag> dest,
                                     double centerValue = 0.4)
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor>
        void multiPyramidReduceBurtFilter(triple<SrcIterator, SrcShape, SrcAccessor> a,
                                     triple<DestIterator, DestShape, DestAccessor> b,
                                     double centerValue = 0.4)
    }
    \endcode

    use a \ref vigra::MultiArrayPyramid :
    \code
    namespace vigra {
        template <class MultiArrayType, class Alloc>
        inline void multiPyramidReduceBurtFilter(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                                     int fromLevel, int toLevel,
                                     double centerValue = 0.4)
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void multiPyramidReduceBurtFilter)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestShape>
void multiPyramidReduceBurtFilter(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestShape const & dshape, DestAccessor dest,
                             double centerValue = 0.4)
{
    vigra_precondition(0.25 <= centerValue && centerValue <= 0.5,
             "multiPyramidReduceBurtFilter(): centerValue must be between 0.25 and 0.5.");

    Kernel1D<double> kern;
    kern.initExplicitly(-2, 2) = 0.25 - centerValue / 2.0, 0.25, centerValue, 0.25, 0.25 - centerValue / 2.0;

    resamplingSeparableConvolveMultiArray(s, shape, src, d, dshape, dest, kern, Rational<int>(1, 2));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void multiPyramidReduceBurtFilter(triple<SrcIterator, SrcShape, SrcAccessor> a,
                             triple<DestIterator, DestShape, DestAccessor> b,
                             double centerValue = 0.4)
{
    multiPyramidReduceBurtFilter(a.first, a.second, a.third, b.first, b.second, b.third, centerValue);
}

template <unsigned int srcN, class srcT, class srcStrideTag,
          unsigned int destN, class destT, class destStrideTag>
void multiPyramidReduceBurtFilter(MultiArrayView<srcN, srcT, srcStrideTag> src,
                             MultiArrayView<destN, destT, destStrideTag> dest,
                             double centerValue = 0.4)
{
    multiPyramidReduceBurtFilter(srcMultiArrayRange(src), destMultiArrayRange(dest), centerValue);
}

template <class MultiArrayType, class Alloc>
inline void multiPyramidReduceBurtFilter(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                             int fromLevel, int toLevel,
                             double centerValue = 0.4)
{
    vigra_precondition(fromLevel  < toLevel,
       "multiPyramidReduceBurtFilter(): fromLevel must be smaller than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= fromLevel && toLevel <= pyramid.highestLevel(),
       "multiPyramidReduceBurtFilter(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel+1; i <= toLevel; ++i)
        multiPyramidReduceBurtFilter(pyramid[i-1], pyramid[i], centerValue);
}


/** \brief Two-fold up-sampling for image pyramid reconstruction.

    This function implements the expansion by one resolution level (first signature)
    or across several pyramid levels (last signature) of a Gaussian pyramid as described in

    P. Burt, E. Adelson: <i>"The Laplacian Pyramid as a Compact Image Code"</i>, IEEE Trans. Communications, 9(4):532–540, 1983

    That is, the function first places the pixel values of the low-resolution
    image at the even pixel coordinates of the high-resolution image (pixels with
    at least one odd coordinate are zero-initialized) and then applies the
    interpolation filter
    \code
    [0.5 - centerValue, 0.5, 2*centerValue, 0.5, 0.5 - centerValue]
    \endcode
    to the high-resolution image. The source image shape must be
    <tt>src_shape = ceil(dest_shape / 2.0)</tt>.
    <tt>centerValue</tt> must be between 0.25 and 0.5 and determines the sharpness
    of the interpolation (bigger values correspond to sharper images).
    If <tt>fromLevel - toLevel > 1</tt> in the pyramid variant of the function,
    this process is repeated until <tt>toLevel</tt> is reached.

    Typically, this functions is used in connection with a \ref vigra::MultiArrayPyramid
    (last signature below) to perform several levels of expansion in one go.

    <b> Declarations:</b>

    <b>\#include</b> \<vigra/multi_pyramid.hxx\><br>
    Namespace: vigra

    pass multi array iterators, shapes and accessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class DestShape>
        void multiPyramidExpandBurtFilter(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                     DestIterator d, DestShape const & dshape, DestAccessor dest,
                                     double centerValue = 0.4)
    }
    \endcode

    pass multi array views directly:
    \code
    namespace vigra {
        template <unsigned int srcN, class srcT, class srcStrideTag,
                  unsigned int destN, class destT, class destStrideTag>
        void multiPyramidExpandBurtFilter(MultiArrayView<srcN, srcT, srcStrideTag> src,
                                     MultiArrayView<destN, destT, destStrideTag> dest,
                                     double centerValue = 0.4)
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestShape, class DestAccessor>
        void multiPyramidExpandBurtFilter(triple<SrcIterator, SrcShape, SrcAccessor> a,
                                     triple<DestIterator, DestShape, DestAccessor> b,
                                     double centerValue = 0.4)
    }
    \endcode

    use a \ref vigra::MultiArrayPyramid :
    \code
    namespace vigra {
        template <class MultiArrayType, class Alloc>
        inline void multiPyramidExpandBurtFilter(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                                     int fromLevel, int toLevel,
                                     double centerValue = 0.4)
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void multiPyramidExpandBurtFilter)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestShape>
void multiPyramidExpandBurtFilter(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestShape const & dshape, DestAccessor dest,
                             double centerValue = 0.4)
{
    vigra_precondition(0.25 <= centerValue && centerValue <= 0.5,
             "multiPyramidExpandBurtFilter(): centerValue must be between 0.25 and 0.5.");

    ArrayVector<Kernel1D<double> > kernels(2);
    kernels[0].initExplicitly(-1, 1) = 0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue;
    kernels[1].initExplicitly(-1, 0) = 0.5, 0.5;

    resamplingSeparableConvolveMultiArray(s, shape, src, d, dshape, dest, kernels, Rational<int>(2, 1));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void multiPyramidExpandBurtFilter(triple<SrcIterator, SrcShape, SrcAccessor> a,
                             triple<DestIterator, DestShape, DestAccessor> b,
                             double centerValue = 0.4)
{
    multiPyramidExpandBurtFilter(a.first, a.second, a.third, b.first, b.second, b.third, centerValue);
}

template <unsigned int srcN, class srcT, class srcStrideTag,
          unsigned int destN, class destT, class destStrideTag>
void multiPyramidExpandBurtFilter(MultiArrayView<srcN, srcT, srcStrideTag> src,
                             MultiArrayView<destN, destT, destStrideTag> dest,
                             double centerValue = 0.4)
{
    multiPyramidExpandBurtFilter(srcMultiArrayRange(src), destMultiArrayRange(dest), centerValue);
}

template <class MultiArrayType, class Alloc>
inline void multiPyramidExpandBurtFilter(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                             int fromLevel, int toLevel,
                             double centerValue = 0.4)
{
    vigra_precondition(fromLevel  > toLevel,
       "multiPyramidExpandBurtFilter(): fromLevel must be larger than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= toLevel && fromLevel <= pyramid.highestLevel(),
       "multiPyramidExpandBurtFilter(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel-1; i >= toLevel; --i)
        multiPyramidExpandBurtFilter(pyramid[i+1], pyramid[i], centerValue);
}

/** \brief Create a Laplacian pyramid.

    This function implements the reduction across several resolution levels of
    a Laplacian pyramid as described in

    P. Burt, E. Adelson: <i>"The Laplacian Pyramid as a Compact Image Code"</i>, IEEE Trans. Communications, 9(4):532–540, 1983

    It first creates a Gaussian pyramid using \ref multiPyramidReduceBurtFilter(), then
    upsamples each level once using \ref multiPyramidExpandBurtFilter(), and finally
    stores the difference between the upsampled and original versions of
    each level (i.e. the Laplacian of Gaussian is approximated by a difference
    of Gaussian).

    <b>\#include</b> \<vigra/multi_pyramid.hxx\><br>
    Namespace: vigra
*/
template <class MultiArrayType, class Alloc>
inline void
multiPyramidReduceBurtLaplacian(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                           int fromLevel, int toLevel,
                           double centerValue = 0.4)
{
    using namespace functor;

    vigra_precondition(toLevel  > fromLevel,
       "multiPyramidReduceBurtLaplacian(): toLevel must be larger than fromLevel.");
    vigra_precondition(pyramid.lowestLevel() <= toLevel && fromLevel <= pyramid.highestLevel(),
       "multiPyramidReduceBurtLaplacian(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    multiPyramidReduceBurtFilter(pyramid, fromLevel, toLevel, centerValue);
    for(int i=fromLevel; i < toLevel; ++i)
    {
        typename MultiArrayPyramid<MultiArrayType, Alloc>::value_type tmpImage(pyramid[i].shape());
        multiPyramidExpandBurtFilter(pyramid[i+1], tmpImage, centerValue);
        combineTwoMultiArrays(tmpImage, pyramid[i], pyramid[i],
                       Arg1() - Arg2());
    }
}

/** \brief Reconstruct a Laplacian pyramid.

    This function implements the reconstruction of a Gaussian pyramid
    across several resolution levels of a Laplacian pyramid as described in

    P. Burt, E. Adelson: <i>"The Laplacian Pyramid as a Compact Image Code"</i>, IEEE Trans. Communications, 9(4):532–540, 1983

    At each level starting from <tt>fromLevel</tt>, this function calls
    \ref multiPyramidExpandBurtFilter() to interpolate the image to the next highest
    resolution, and then adds the interpolated image to the image stored at the
    next level.

    <b>\#include</b> \<vigra/multi_pyramid.hxx\><br>
    Namespace: vigra
*/
template <class MultiArrayType, class Alloc>
inline void
multiPyramidExpandBurtLaplacian(MultiArrayPyramid<MultiArrayType, Alloc> & pyramid,
                           int fromLevel, int toLevel,
                           double centerValue = 0.4)
{
    using namespace functor;

    vigra_precondition(fromLevel  > toLevel,
       "multiPyramidExpandBurtLaplacian(): fromLevel must be larger than toLevel.");
    vigra_precondition(pyramid.lowestLevel() <= toLevel && fromLevel <= pyramid.highestLevel(),
       "multiPyramidExpandBurtLaplacian(): fromLevel and toLevel must be between the lowest and highest pyramid levels (inclusive).");

    for(int i=fromLevel-1; i >= toLevel; --i)
    {
        typename MultiArrayPyramid<MultiArrayType, Alloc>::value_type tmpImage(pyramid[i].shape());
        multiPyramidExpandBurtFilter(pyramid[i+1], tmpImage, centerValue);
        combineTwoMultiArrays(tmpImage, pyramid[i], pyramid[i],
                       Arg1() - Arg2());
    }
}

}

#endif