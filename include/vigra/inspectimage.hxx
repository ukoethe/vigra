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
 
 
#ifndef VIGRA_INSPECTIMAGE_HXX
#define VIGRA_INSPECTIMAGE_HXX

#include <vector>
#include <algorithm>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/rgbvalue.hxx"

namespace vigra {

/** \addtogroup InspectAlgo Algorithms to Inspect Images

    Apply read-only functor to every pixel
*/
//@{

/********************************************************/
/*                                                      */
/*                        inspectImage                  */
/*                                                      */
/********************************************************/

/** \brief Apply read-only functor to every pixel in the image.

    This function can be used to collect statistics of the image etc.
    The results must be stored in the functor, which serves as a return
    value.
    The function uses an accessor to access the pixel data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class Functor>
        void
        inspectImage(ImageIterator upperleft, ImageIterator lowerright,
                     Accessor a, Functor & f)
    }
    \endcode
        
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class Functor>
        void
        inspectImage(triple<ImageIterator, ImageIterator, Accessor> img,
             Functor & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    // init functor
    vigra::BImage img;
    
    vigra::FindMinMax<vigra::BImage::PixelType> minmax();   
    
    vigra::inspectImage(srcImageRange(img), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \endcode

    <b> Required Interface:</b>
    
    \code
    ConstImageIterator upperleft, lowerright;
    IteratorTraits<ConstImageIterator>::RowIterator ix(upperleft);
    
    Accessor accessor;
    Functor functor;
    
    functor(accessor(ix));         // return not used
    \endcode
    
*/
template <class ImageIterator, class Accessor, class Functor>
void
inspectImage(ImageIterator upperleft, ImageIterator lowerright,
         Accessor a, Functor & f)
{
    int w = lowerright.x - upperleft.x;
    
    for(; upperleft.y<lowerright.y; ++upperleft.y)
    {
        typename IteratorTraits<ImageIterator>::RowIterator ix(upperleft);
        typename IteratorTraits<ImageIterator>::RowIterator ixend = ix + w;
    
        for(; ix < ixend; ++ix)
        {
            f(a(ix));
        }
    }
}
    
template <class ImageIterator, class Accessor, class Functor>
inline
void
inspectImage(triple<ImageIterator, ImageIterator, Accessor> img,
         Functor & f)
{
    inspectImage(img.first, img.second, img.third, f);
}

/********************************************************/
/*                                                      */
/*                      inspectImageIf                  */
/*                                                      */
/********************************************************/

/** \brief Apply read-only functor to every pixel in the ROI.

    This function can be used to collect statistics of the roi etc.
    The functor is called whenever the return value of the mask's 
    accessor is not zero.
    The results must be stored in the functor, which serves as a return
    value.
    Accessors are used to access the pixel and mask data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor,
                  class MaskImageIterator, class MaskAccessor, class Functor>
        void
        inspectImageIf(ImageIterator upperleft, ImageIterator lowerright, 
               MaskImageIterator mask_upperleft, MaskAccessor ma,
               Functor & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor,
              class MaskImageIterator, class MaskAccessor, class Functor>
        void
        inspectImageIf(triple<ImageIterator, ImageIterator, Accessor> img,
               pair<MaskImageIterator, MaskAccessor> mask,
               Functor & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::BImage img(100, 100);
    vigra::BImage mask(100, 100);
    
    // init functor
    vigra::FindMinMax<vigra::BImage::PixelType> minmax();    
    
    vigra::inspectImageIf(srcImageRange(img), 
                          maskImage(mask), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \endcode

    <b> Required Interface:</b>
    
    \code
    ConstImageIterator upperleft, lowerright;
    MaskImageIterator mask_upperleft;
    IteratorTraits<ConstImageIterator>::RowIterator ix(upperleft);
    IteratorTraits<MaskImageIterator>::RowIterator mx(mask_upperleft);
    
    Accessor accessor;
    MaskAccessor mask_accessor;
    
    Functor functor;

    if(mask_accessor(mx)) functor(accessor(ix)); 
    \endcode
    
*/
template <class ImageIterator, class Accessor,
      class MaskImageIterator, class MaskAccessor, class Functor>
void
inspectImageIf(ImageIterator upperleft, 
               ImageIterator lowerright, Accessor a,
           MaskImageIterator mask_upperleft, MaskAccessor ma,
           Functor & f)
{
    int w = lowerright.x - upperleft.x;
            
    for(; upperleft.y<lowerright.y; ++upperleft.y, ++mask_upperleft.y)
    {
        typename IteratorTraits<ImageIterator>::RowIterator ix(upperleft);
        typename IteratorTraits<ImageIterator>::RowIterator ixend = ix + w;
        typename IteratorTraits<MaskImageIterator>::RowIterator mx(mask_upperleft);
    
        for(; ix < ixend; ++ix, ++mx)
        {
            if(ma(mx)) 
                f(a(ix));
        }
    }
}

template <class ImageIterator, class Accessor,
      class MaskImageIterator, class MaskAccessor, class Functor>
inline
void
inspectImageIf(triple<ImageIterator, ImageIterator, Accessor> img,
               pair<MaskImageIterator, MaskAccessor> mask,
               Functor & f)
{
    inspectImageIf(img.first, img.second, img.third, 
                   mask.first, mask.second, f);
}

/********************************************************/
/*                                                      */
/*                  inspectTwoImages                    */
/*                                                      */
/********************************************************/

/** \brief Apply read-only functor to every pixel of both images.

    This function can be used to collect statistics for each region of a 
    labeled image, especially in conjunction with 
    the \ref ArrayOfRegionStatistics functor. The results must be 
    stored in the functor which serves as a return value. 
    Accessors are used to access the pixel data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator1, class Accessor1,
              class ImageIterator2, class Accessor2, 
              class Functor>
        void
        inspectTwoImages(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                 ImageIterator2 upperleft2, Accessor2 a2,
                 Functor & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator1, class Accessor1,
              class ImageIterator2, class Accessor2, 
              class Functor>
        void
        inspectTwoImages(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
                         pair<ImageIterator2, Accessor2> img2,
                 Functor & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::BImage image1;
    vigra::BImage image2;
    
    SomeStatisticsFunctor stats(...);     // init functor
    
    vigra::inspectTwoImages(srcImageRange(image1), srcImage(image2),
                            region_stats);
    
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    ImageIterator1 upperleft1, lowerright1;
    ImageIterator2 upperleft2;
    IteratorTraits<ImageIterator1>::RowIterator ix1(upperleft1);
    IteratorTraits<ImageIterator2>::RowIterator ix2(upperleft2);
    
    Accessor1 accessor1;
    Accessor2 accessor2;
    
    Functor functor;
    functor(accessor1(ix1), accessor2(ix2));  // return not used
    \endcode
    
*/
template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
      class Functor>
void
inspectTwoImages(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                 ImageIterator2 upperleft2, Accessor2 a2,
         Functor & f)
{
    int w = lowerright1.x - upperleft1.x;
    
    for(; upperleft1.y<lowerright1.y; ++upperleft1.y, ++upperleft2.y)
    {
        typename IteratorTraits<ImageIterator1>::RowIterator ix1(upperleft1);
        typename IteratorTraits<ImageIterator1>::RowIterator ix1end = ix1 + w;
        typename IteratorTraits<ImageIterator2>::RowIterator ix2(upperleft2);
    
        for(; ix1 < ix1end; ++ix1, ++ix2)
        {
            f(a1(ix1), a2(ix2));
        }
    }
}

template <class ImageIterator1, class Accessor1,
      class ImageIterator2, class Accessor2, 
      class Functor>
inline
void
inspectTwoImages(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
         pair<ImageIterator2, Accessor2> img2,
         Functor & f)
{
    inspectTwoImages(img1.first, img1.second, img1.third, 
                     img2.first, img2.second, f);
}

/********************************************************/
/*                                                      */
/*                inspectTwoImagesIf                    */
/*                                                      */
/********************************************************/

/** \brief Apply read-only functor to those pixels of both images where
    the mask image is non-zero.
    
    This function can be used to collect statistics for selected regions of a 
    labeled image, especially in conjunction with 
    the \ref ArrayOfRegionStatistics functor. The results must be 
    stored in the functor which serves as a return value. 
    Accessors are used to access the pixel data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator1, class Accessor1,
                  class ImageIterator2, class Accessor2, 
                  class MaskImageIterator, class MaskAccessor, 
                  class Functor>
        void
        inspectTwoImagesIf(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                         ImageIterator2 upperleft2, Accessor2 a2,
                         MaskImageIterator mupperleft, MaskAccessor mask,
                         Functor & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator1, class Accessor1,
                  class ImageIterator2, class Accessor2, 
                  class MaskImageIterator, class MaskAccessor, 
                  class Functor>
        void
        inspectTwoImagesIf(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
                 pair<ImageIterator2, Accessor2> img2,
                 pair<MaskImageIterator, MaskAccessor> mimg,
                 Functor & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::BImage image1;
    vigra::BImage image2;
    vigra::BImage maskimage;
    
    SomeStatisticsFunctor stats(...);     // init functor
    
    vigra::inspectTwoImagesIf(srcImageRange(image1), srcImage(image2),
                              srcImage(maskimage), region_stats);
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    ImageIterator1 upperleft1, lowerright1;
    ImageIterator2 upperleft2;
    MaskImageIterator upperleftm;
    IteratorTraits<ImageIterator1>::RowIterator ix1(upperleft1);
    IteratorTraits<ImageIterator2>::RowIterator ix2(upperleft2);
    IteratorTraits<MaskImageIterator>::RowIterator mx(mupperleft);
    
    Accessor1 accessor1;
    Accessor2 accessor2;
    MaskAccessor mask;
    
    Functor functor;
    if(mask(mx))
        functor(accessor1(ix1), accessor2(ix2)); 
    \endcode
    
*/
template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class MaskImageIterator, class MaskAccessor, 
      class Functor>
void
inspectTwoImagesIf(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                 ImageIterator2 upperleft2, Accessor2 a2,
                 MaskImageIterator mupperleft, MaskAccessor mask,
                 Functor & f)
{
    int w = lowerright1.x - upperleft1.x;
    
    for(; upperleft1.y<lowerright1.y; ++upperleft1.y, ++upperleft2.y, ++mupperleft.y)
    {
        typename IteratorTraits<ImageIterator1>::RowIterator ix1(upperleft1);
        typename IteratorTraits<ImageIterator1>::RowIterator ix1end = ix1 + w;
        typename IteratorTraits<ImageIterator2>::RowIterator ix2(upperleft2);
        typename IteratorTraits<MaskImageIterator>::RowIterator ixm(mupperleft);
    
        for(; ix1 < ix1end; ++ix1, ++ix2, ++ixm)
        {
            if(mask(ixm)) 
                f(a1(ix1), a2(ix2));
        }
    }
}

template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class MaskImageIterator, class MaskAccessor, 
          class Functor>
inline
void
inspectTwoImagesIf(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
         pair<ImageIterator2, Accessor2> img2,
         pair<MaskImageIterator, MaskAccessor> m,
         Functor & f)
{
    inspectTwoImagesIf(img1.first, img1.second, img1.third, 
                     img2.first, img2.second, 
                     m.first, m.second, 
                     f);
}

//@}

/** \addtogroup InspectFunctor Functors To Inspect Images
    Functors which report image statistics
*/
//@{

/********************************************************/
/*                                                      */
/*                     FindMinMax                       */
/*                                                      */
/********************************************************/

/** \brief Find the minimum and maximum pixel value in an image or ROI.

    In addition the size of the ROI is calculated.
    This Functor can also be used in conjunction with
    \ref ArrayOfRegionStatistics to find the extremes of all regions in
    a labeled image.
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::BImage img;
    
    vigra::FindMinMax<BImage::PixelType> minmax;   // init functor
    
    vigra::inspectImage(srcImageRange(img), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    VALUETYPE v1, v2(v1);
    
    v1 < v2;
    v1 = v2;
    \endcode
    
*/
template <class VALUETYPE>
class FindMinMax
{
   public:
   
        /** the functor's value type
    */
    typedef VALUETYPE value_type;
   
    /** init min and max
    */
    FindMinMax()
    : count(0)
    {}
    
    /** update min and max
    */
    void operator()(VALUETYPE const & v) 
    {
        if(count)
        {
            if(v < min) min = v;
            if(max < v) max = v;
        }
        else
        {
            min = v;
            max = v;
        }
        ++count;
    }
    
    /** update min and max with components of RGBValue<VALUETYPE>
    */
    void operator()(RGBValue<VALUETYPE> const & v) 
    {
        operator()(v.red());
        operator()(v.green());
        operator()(v.blue());
    }
    
    /** merge two statistics
    */
    void operator()(FindMinMax const & v) 
    {
        if(v.count)
        {
            if(count)
                {
                    if(v.min < min) min = v.min;
                if((this->max) < v.max) max = v.max;
                }
                else
                {
                min = v.min;
                max = v.max;
                }
        }
        count += v.count;
    }
    
    /** the current min
    */
    VALUETYPE min;
    
    /** the current max
    */
    VALUETYPE max;
    
    /** the number of values processed so far
    */
    unsigned int count;
    
}; 

/********************************************************/
/*                                                      */
/*                    FindAverage                       */
/*                                                      */
/********************************************************/

/** \brief  Find the average pixel value in an image or ROI.

    In addition the size of the ROI is calculated.
    This Functor can also be used in conjunction with
    \ref ArrayOfRegionStatistics to find the average of all regions in
    a labeled image.
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::BImage img;
    
    vigra::FindAverage<vigra::BImage::PixelType> average;   // init functor
    
    vigra::inspectImage(srcImageRange(img), average);
    
    cout << "Average: " << average();
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    VALUETYPE v1, v2(v1);
    
    v1 < v2;
    v1 = v2;
    \endcode
    
*/
template <class VALUETYPE>
class FindAverage
{
   public:
   
        /** the functor's value type
    */
    typedef typename NumericTraits<VALUETYPE>::RealPromote value_type;
   
    /** init average
    */
    FindAverage()
    : count(0), sum(NumericTraits<value_type>::zero())
    {}
    
    /** update average
    */
    void operator()(VALUETYPE const & v) 
    {
        sum += v;
        ++count;
    }
    
    /** merge two statistics
    */
    void operator()(FindAverage const & v) 
    {
        sum += v.sum;
        count += v.count;
    }
    
    /** return current average
    */
    value_type average() const
    {
        return sum / (double)count;
    }
    
    /** return current average
    */
    value_type operator()() const
    {
        return sum / (double)count;
    }
    
    unsigned int count;
    value_type sum;
    
}; 

/********************************************************/
/*                                                      */
/*                    FindROISize                       */
/*                                                      */
/********************************************************/

/** \brief Calculate the size of an ROI in an image. 

    This Functor is often used in conjunction with
    \ref ArrayOfRegionStatistics to find the sizes of all regions in
    a labeled image.
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img, mask;
    
    vigra::FindROISize<vigra::BImage::PixelType> roisize;   // init functor
    
    vigra::inspectImageIf(srcImageRange(img), srcImage(mask), roisize);
    
    cout << "Size of ROI: " << roisize.count;
    
    \endcode

*/
template <class VALUETYPE>
class FindROISize
{
   public:
   
        /** the functor's value type
    */
    typedef VALUETYPE value_type;
   
    /** init counter to 0
    */
    FindROISize()
    : count(0)
    {}
    
    /** update counter
    */
    void operator()(VALUETYPE const &) 
    {
        ++count;
    }
    
    /** return current size
    */
    int operator()() const 
    {
        return count;
    }
    
    /** return current size
    */
    int size() const 
    {
        return count;
    }
    
    /** merge two statistics
    */
    void operator()(FindROISize const & o) 
    {
        count += o.count;
    }
    
    /** the current counter
    */
    unsigned int count;
    
}; 

/********************************************************/
/*                                                      */
/*                FindBoundingRectangle                 */
/*                                                      */
/********************************************************/

/** \brief Calculate the bounding rectangle of an ROI in an image. 

    As always in VIGRA, <TT>roiRect.lowerRight</TT> is <em> just outside the rectangle</em>.
    That is, the last pixel actually in the rectangle is <TT>roiRect.lowerRight - Diff2D(1,1)</TT>.
    This Functor is often used in conjunction with
    \ref ArrayOfRegionStatistics to find the bounding rectangles 
    of all regions in a labeled image.
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img, mask;
    ...
    
    vigra::FindBoundingRectangle roiRect;   // init functor
    
    // Diff2D is used as the iterator for the source image. This
    // simulates an image where each pixel value equals that pixel's 
    // coordinates. Tha image 'mask' determines the ROI.
    vigra::inspectImageIf(srcIterRange(Diff2D(0,0), img.size()), 
                          srcImage(mask), roiRect);
    
    cout << "Upper left of ROI: " << 
        roiRect.upperLeft.x << ", " << roiRect.upperLeft.y << endl;
    cout << "Lower right of ROI: " << 
        roiRect.lowerRight.x << ", " << roiRect.lowerRight.y << endl;
    
    \endcode

*/
class FindBoundingRectangle
{
  public:
        /** the functor's value type
        */
    typedef Diff2D value_type;
    
        /** Upper left of the region as seen so far
        */
    Diff2D upperLeft;
        /** Lower right of the region as seen so far
        */
    Diff2D lowerRight;
    
        /** are the functors contents valid ?
        */
    bool valid;
    
        /** init rectangle to invalid values
        */
    FindBoundingRectangle()
    : valid(false)
    {}
    
        /** update rectangle by including the coordinate coord
        */
    void operator()(Diff2D const & coord) 
    {
        if(!valid)
        {
            upperLeft = coord;
            lowerRight = coord + Diff2D(1,1);
            valid = true;
        }
        else
        {
            upperLeft.x = std::min(upperLeft.x, coord.x);
            upperLeft.y = std::min(upperLeft.y, coord.y);
            lowerRight.x = std::max(lowerRight.x, coord.x + 1);
            lowerRight.y = std::max(lowerRight.y, coord.y + 1);
        }
    }
    
        /** update rectangle by merging it with another rectangle
        */
    void operator()(FindBoundingRectangle const & otherRegion) 
    {
        if(!valid)
        {
            upperLeft = otherRegion.upperLeft;
            lowerRight = otherRegion.lowerRight;
            valid = otherRegion.valid;
        }
        else if(otherRegion.valid)
        {
            upperLeft.x = std::min(upperLeft.x, otherRegion.upperLeft.x);
            upperLeft.y = std::min(upperLeft.y, otherRegion.upperLeft.y);
            lowerRight.x = std::max(lowerRight.x, otherRegion.lowerRight.x);
            lowerRight.y = std::max(lowerRight.y, otherRegion.lowerRight.y);
        }
    }
};

/********************************************************/
/*                                                      */
/*                 LastValueFunctor                     */
/*                                                      */
/********************************************************/

/** \brief Stores and returns the last value it has seen. 

    This Functor is best used in conjunction with
    \ref ArrayOfRegionStatistics to realize a look-up table.
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img;
    
    vigra::ArrayOfRegionStatistics<LastValueFunctor<unsigned char> > lut(255);
    
    for(int i=0; i<256; ++i)
    {
        lut[i] = ...; // init look-up table
    }
    
    vigra::transformImage(srcImageRange(img), destImage(img), lut);
    
    \endcode

*/
template <class VALUETYPE>
class LastValueFunctor
{
   public:
   
        /** the functor's value type
    */
    typedef VALUETYPE value_type;
   
    /** default initialization of value
    */
    LastValueFunctor()
    {}

        /** replace value
    */
    void operator=(VALUETYPE const & v) { value = v; }
        
    /** replace value
    */
    void operator()(VALUETYPE const & v) { value = v; }
    
    /** return current value
    */
    VALUETYPE operator()() const { return value; }

    /** the current value
    */
    VALUETYPE value;
    
}; 

/********************************************************/
/*                                                      */
/*              ArrayOfRegionStatistics                 */
/*                                                      */
/********************************************************/

/** \brief Calculate statistics for all regions of a labeled image. 

    This Functor encapsulates an array of statistics functors, one
    for each label, and selects the one to be updated according to the 
    pixel's label. 
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img;
    vigra::IImage labels;
    int max_label;
    ...
    
    // init functor as an array of 'max_label' FindMinMax-Functors
    vigra::ArrayOfRegionStatistics<vigra::FindMinMax<vigra::BImage::PixelType> > 
                                                         minmax(max_label);   
        
    vigra::inspectTwoImages(srcImageRange(img), srcImage(labels), minmax);
    
    for(int i=0; i<= max_label; ++i)
    {
        cout << "Max gray lavel of region " << i << ": " 
             << minmax.region[i].max << endl;
    }
    
    // init functor as an array of 'max_label' FindAverage-Functors
    vigra::ArrayOfRegionStatistics<vigra::FindAverage<vigra::BImage::PixelType> > 
                                                         average(max_label);   
        
    vigra::inspectTwoImages(srcImageRange(img), srcImage(labels), average);
    
    // write back the average of each region into the original image
    vigra::transformImage(srcImageRange(labels), destImage(img), average);
    
    \endcode
    
    <b> Required Interface:</b>

    \code
    RegionStatistics region;
    RegionStatistics::value_type v;
    
    region(v);     // update statistics
    v = region();  // return statistics
    
    \endcode
*/
template <class RegionStatistics>
class ArrayOfRegionStatistics
{
    typedef std::vector<RegionStatistics> RegionArray;

  public:
         ///
    typedef typename RegionStatistics::value_type value_type;
         ///
    typedef typename RegionArray::iterator iterator;

    /** init array of RegionStatistics with default size 0.
    */
    ArrayOfRegionStatistics()
    {}
    
    /** init array of RegionStatistics with index domain
        0...max_region_label.
    */
    ArrayOfRegionStatistics(int max_region_label)
    : regions(max_region_label+1)
    {}
    
    /** resize array to new index domain 0...max_region_label.
        All bin are re-initialized.
    */
    void resize(int max_region_label)
    {
        RegionArray newRegions(max_region_label+1);
        regions.swap(newRegions);
    }
    
    /** reset the contained functors to their initial state.
    */
    void reset()
    {
        RegionArray newRegions(regions.size());
        regions.swap(newRegions);
    }
    
    /** update regions statistics for region 'label'
    */
    void operator()(value_type const & v, int label) {
        regions[label](v);
    }
    
    /** merge second region into first
    */
    void merge(int label1, int label2) {
        regions[label1](regions[label2]);
    }
    
    /** ask for maximal index (label) allowed
    */
    int maxRegionLabel() const { return size() - 1; } 

    /** ask for array size (i.e. maxRegionLabel() + 1)
    */
    int size() const { return regions.size(); } 

    /** access the statistics for a region via its label
    */
    value_type operator()(int label) const { return regions[label](); }
    
    /** read the statistics functor for a region via its label
    */
    RegionStatistics const & operator[](int label) const { return regions[label]; }
    
    /** access the statistics functor for a region via its label
    */
    RegionStatistics & operator[](int label) { return regions[label]; }
    
    /** iterator to the begin of the region array
    */
    iterator begin() { return regions.begin(); } 
    
    /** iterator to the end of the region array
    */
    iterator end() { return regions.end(); } 
    
  private:
    std::vector<RegionStatistics> regions;
};


//@}


} // namespace vigra

#endif // VIGRA_INSPECTIMAGE_HXX
