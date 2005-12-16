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


#ifndef VIGRA_INSPECTIMAGE_HXX
#define VIGRA_INSPECTIMAGE_HXX

#include <vector>
#include <algorithm>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/iteratortraits.hxx"
#include "vigra/functortraits.hxx"
#include "vigra/rgbvalue.hxx"

namespace vigra {

/** \addtogroup InspectAlgo Algorithms to Inspect Images

    Apply read-only functor to every pixel
*/
//@{

/********************************************************/
/*                                                      */
/*                      inspectLine                     */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor, class Functor>
void
inspectLine(SrcIterator s,
            SrcIterator send, SrcAccessor src,
            Functor & f)
{
    for(; s != send; ++s)
        f(src(s));
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class Functor>
void
inspectLineIf(SrcIterator s,
              SrcIterator send, SrcAccessor src,
              MaskIterator m, MaskAccessor mask,
              Functor & f)
{
    for(; s != send; ++s, ++m)
        if(mask(m))
            f(src(s));
}

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class Functor>
void
inspectTwoLines(SrcIterator1 s1,
                SrcIterator1 s1end, SrcAccessor1 src1,
                SrcIterator2 s2, SrcAccessor2 src2,
                Functor & f)
{
    for(; s1 != s1end; ++s1, ++s2)
        f(src1(s1), src2(s2));
}

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class MaskIterator, class MaskAccessor,
          class Functor>
void
inspectTwoLinesIf(SrcIterator1 s1,
                  SrcIterator1 s1end, SrcAccessor1 src1,
                  SrcIterator2 s2, SrcAccessor2 src2,
                  MaskIterator m, MaskAccessor mask,
                  Functor & f)
{
    for(; s1 != s1end; ++s1, ++s2, ++m)
        if(mask(m))
            f(src1(s1), src2(s2));
}

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

    use argument objects in conjunction with \ref ArgumentObjectFactories:
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

    vigra::FindMinMax<vigra::BImage::PixelType> minmax;

    vigra::inspectImage(srcImageRange(img), minmax);

    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \endcode

    <b> Required Interface:</b>

    \code
    ConstImageIterator upperleft, lowerright;
    ConstImageIterator::row_iterator ix = upperleft.rowIterator();

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
        inspectLine(upperleft.rowIterator(),
                    upperleft.rowIterator() + w, a, f);
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

namespace functor
{
    template <class T> class UnaryAnalyser;
}

template <class ImageIterator, class Accessor, class Functor>
inline
void
inspectImage(ImageIterator upperleft, ImageIterator lowerright,
         Accessor a, functor::UnaryAnalyser<Functor> const & f)
{
    inspectImage(upperleft, lowerright, a,
                 const_cast<functor::UnaryAnalyser<Functor> &>(f));
}

template <class ImageIterator, class Accessor, class Functor>
inline
void
inspectImage(triple<ImageIterator, ImageIterator, Accessor> img,
         functor::UnaryAnalyser<Functor> const & f)
{
    inspectImage(img.first, img.second, img.third,
                 const_cast<functor::UnaryAnalyser<Functor> &>(f));
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


    use argument objects in conjunction with \ref ArgumentObjectFactories:
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
    ConstImageIterator::row_iterator ix = upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();

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
        inspectLineIf(upperleft.rowIterator(),
                      upperleft.rowIterator() + w, a,
                      mask_upperleft.rowIterator(), ma, f);
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


    use argument objects in conjunction with \ref ArgumentObjectFactories:
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
                            stats);


    \endcode

    <b> Required Interface:</b>

    \code
    ImageIterator1 upperleft1, lowerright1;
    ImageIterator2 upperleft2;
    ImageIterator1::row_iterator ix1 = upperleft1.rowIterator();
    ImageIterator2::row_iterator ix2 = upperleft2.rowIterator();

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
        inspectTwoLines(upperleft1.rowIterator(),
                        upperleft1.rowIterator() + w, a1,
                        upperleft2.rowIterator(), a2, f);
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


    use argument objects in conjunction with \ref ArgumentObjectFactories:
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
    ImageIterator1::row_iterator ix1 = upperleft1.rowIterator();
    ImageIterator2::row_iterator ix2 = upperleft2.rowIterator();
    MaskImageIterator::row_iterator mx = mupperleft.rowIterator();

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
        inspectTwoLinesIf(upperleft1.rowIterator(),
                          upperleft1.rowIterator() + w, a1,
                          upperleft2.rowIterator(), a2,
                          mupperleft.rowIterator(), mask, f);
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
    These functors can also be used in conjunction with
    \ref ArrayOfRegionStatistics to find the extremes of all regions in
    a labeled image.
    
    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt> is true (<tt>VigraTrueType<tt>)
    
    <b> Usage:</b>

        <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img;

    vigra::FindMinMax<vigra::BImage::PixelType> minmax;   // init functor

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

        /** the functor's argument type
        */
    typedef VALUETYPE argument_type;

        /** the functor's result type
        */
    typedef VALUETYPE result_type;

        /** \deprecated use argument_type
        */
    typedef VALUETYPE value_type;

        /** init min and max
        */
    FindMinMax()
    : count(0)
    {}

        /** (re-)init functor (clear min, max)
        */
    void reset()
    {
        count = 0;
    }

        /** update min and max
        */
    void operator()(argument_type const & v)
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

template <class VALUETYPE>
class FunctorTraits<FindMinMax<VALUETYPE> >
: public FunctorTraitsBase<FindMinMax<VALUETYPE> >
{
  public:
    typedef VigraTrueType isUnaryAnalyser;
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

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt> and <tt>FunctorTraits::isInitializer</tt>
    are true (<tt>VigraTrueType<tt>)
    
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

        /** the functor's argument type
        */
    typedef VALUETYPE argument_type;

        /** the functor's result type
        */
    typedef typename NumericTraits<VALUETYPE>::RealPromote result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef typename NumericTraits<VALUETYPE>::RealPromote value_type;

        /** init average
        */
    FindAverage()
    : count(0), sum(NumericTraits<result_type>::zero())
    {}

        /** (re-)init average
        */
    void reset()
    {
        count = 0;
        sum = NumericTraits<result_type>::zero();
    }

        /** update average
        */
    void operator()(argument_type const & v)
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
    result_type average() const
    {
        return sum / (double)count;
    }

        /** return current average
        */
    result_type operator()() const
    {
        return sum / (double)count;
    }

    unsigned int count;
    result_type sum;

};

template <class VALUETYPE>
class FunctorTraits<FindAverage<VALUETYPE> >
: public FunctorTraitsBase<FindAverage<VALUETYPE> >
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryAnalyser;
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

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt> and <tt>FunctorTraits::isInitializer</tt>
    are true (<tt>VigraTrueType<tt>)
    
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

        /** the functor's argument type
        */
    typedef VALUETYPE argument_type;

        /** the functor's result type
        */
    typedef unsigned int result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef VALUETYPE value_type;

        /** init counter to 0
        */
    FindROISize()
    : count(0)
    {}

        /** (re-)init ROI size with 0
        */
    void reset()
    {
        count = 0;
    }

        /** update counter
        */
    void operator()(argument_type const &)
    {
        ++count;
    }

        /** return current size
        */
    result_type operator()() const
    {
        return count;
    }

        /** return current size
        */
    result_type size() const
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
    result_type count;

};

template <class VALUETYPE>
class FunctorTraits<FindROISize<VALUETYPE> >
: public FunctorTraitsBase<FindROISize<VALUETYPE> >
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryAnalyser;
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

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt> and <tt>FunctorTraits::isInitializer</tt>
    are true (<tt>VigraTrueType<tt>)
    
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

        /** the functor's argument type
        */
    typedef Diff2D argument_type;

        /** the functors result type
        */
    typedef Rect2D result_type;

        /** \deprecated use argument_type
        */
    typedef Diff2D value_type;

        /** Upper left of the region as seen so far
        */
    Point2D upperLeft;

        /** Lower right of the region as seen so far
        */
    Point2D lowerRight;

        /** are the functors contents valid ?
        */
    bool valid;

        /** init rectangle to invalid values
        */
    FindBoundingRectangle()
    : valid(false)
    {}

        /** (re-)init functor to find other bounds
        */
    void reset()
    {
        valid = false;
    }

        /** update rectangle by including the coordinate coord
        */
    void operator()(argument_type const & coord)
    {
        if(!valid)
        {
            upperLeft = Point2D(coord);
            lowerRight = Point2D(coord + Diff2D(1,1));
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

        /** Get size of current rectangle.
        */
    Size2D size() const
    {
        return lowerRight - upperLeft;
    }

        /** Get current rectangle. <TT>result_type::first</TT> is the upper
            left corner of the rectangle, <TT>result_type::second</TT>
            the lower right.
        */
    result_type operator()() const
    {
        return result_type(upperLeft, lowerRight);
    }
};

template <>
class FunctorTraits<FindBoundingRectangle>
: public FunctorTraitsBase<FindBoundingRectangle>
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryAnalyser;
};

/********************************************************/
/*                                                      */
/*                 LastValueFunctor                     */
/*                                                      */
/********************************************************/

/** \brief Stores and returns the last value it has seen.

    This Functor is best used in conjunction with
    \ref ArrayOfRegionStatistics to realize a look-up table.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt> and <tt>FunctorTraits::isInitializer</tt>
    are true (<tt>VigraTrueType<tt>)
    
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

        /** the functor's argument type
        */
    typedef VALUETYPE argument_type;

        /** the functor's result type
        */
    typedef VALUETYPE result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef VALUETYPE value_type;

        /** default initialization of value
        */
    LastValueFunctor()
    {}

        /** replace value
        */
    void operator=(argument_type const & v) { value = v; }

        /** reset to initia÷ value
        */
    void reset() { value = VALUETYPE(); }

        /** replace value
        */
    void operator()(argument_type const & v) { value = v; }

        /** return current value
        */
    result_type const & operator()() const { return value; }

        /** the current value
        */
    VALUETYPE value;

};

template <class VALUETYPE>
class FunctorTraits<LastValueFunctor<VALUETYPE> >
: public FunctorTraitsBase<LastValueFunctor<VALUETYPE> >
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryAnalyser;
};

/********************************************************/
/*                                                      */
/*                     ReduceFunctor                    */
/*                                                      */
/********************************************************/

/** \brief Apply a functor to reduce the dimensionality of an array.

    This functor can be used to emulate the <tt>reduce</tt> standard function of
    functional programming using <tt>std::for_each()</tt> or <tt>inspectImage()</tt>
    and similar functions. This functor is initialized with a functor encoding
    the expression to be applied, and an accumulator storing the current state
    of the reduction. For each element of the array, the embedded functor is called
    with the accumulator and the current element(s) of the array. The result
    of the reduction is available by calling <tt>reduceFunctor()</tt>. 

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryAnalyser</tt>, <tt>FunctorTraits::isBinaryAnalyser</tt> 
    and <tt>FunctorTraits::isInitializer</tt>
    are true (<tt>VigraTrueType<tt>)
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="inspectimage_8hxx-source.html">vigra/inspectimage.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::BImage img;
    ... // fill the image

    // create a functor to sum the elements of the image
    vigra::ReduceFunctor<std::plus<int>, int> sumElements(std::plus<int>, 0);
    
    vigra::inspectImage(srcImageRange(img), sumElements);

    cout << "The sum of the elements " << sumElements() << endl;

    \endcode

    <b> Required Interface:</b>

    \code
    FUNCTOR f;
    VALUETYPE accumulator, current1, current2;
    
    f(accumulator, current1); // for inspectImage()
    f(accumulator, current1, current2); // for inspectTwoImages()
    \endcode
*/
template <class FUNCTOR, class VALUETYPE>
class ReduceFunctor
{
    FUNCTOR f_;
    VALUETYPE start_, accumulator_;
   public:

        /** the functor's argument type
            when used as a unary inspector.
            (This is not strictly correct since the argument type
            is actuall a template parameter.)
        */
    typedef VALUETYPE argument_type;

        /** the functor's first argument type
            when used as a binary inspector.
            (This is not strictly correct since the argument type
            is actuall a template parameter.)
        */
    typedef VALUETYPE first_argument_type;

        /** the functor's second argument type
            when used as a binary inspector.
            (This is not strictly correct since the argument type
            is actuall a template parameter.)
        */
    typedef VALUETYPE second_argument_type;

        /** the functor's result type
        */
    typedef VALUETYPE result_type;

        /** create with the given functor and initial value \a initial
            for the accumulator.
        */
    ReduceFunctor(FUNCTOR const & f, VALUETYPE const & initial)
    : f_(f),
      start_(initial),
      accumulator_(initial)
    {}
    
        /** Reset accumulator to the initial value.
        */
    void reset()
      { accumulator_ = start_; }

        /** Use binary functor to connect given value with the accumulator.
            The accumulator is used as the first argument, the value \a v
            as the second.
        */
    template <class T>
    void operator()(T const & v) 
    { 
        accumulator_ = f_(accumulator_, v); 
    }

        /** Use ternary functor to connect given values with accumulator.
            The accumulator is used as the first argument, the values \a v1
            ans \a v2 as the second and third.
        */
    template <class T1, class T2>
    void operator()(T1 const & v1, T2 const & v2) 
    { 
        accumulator_ = f_(accumulator_, v1, v2); 
    }

        /** return current value
        */
    result_type const & operator()() const 
      { return accumulator_; }
};

template <class FUNCTOR, class VALUETYPE>
ReduceFunctor<FUNCTOR, VALUETYPE>
reduceFunctor(FUNCTOR const & f, VALUETYPE const & initial)
{
    return ReduceFunctor<FUNCTOR, VALUETYPE>(f, initial);
}

template <class FUNCTOR, class VALUETYPE>
class FunctorTraits<ReduceFunctor<FUNCTOR, VALUETYPE> >
: public FunctorTraitsBase<ReduceFunctor<FUNCTOR, VALUETYPE> >
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryAnalyser;
    typedef VigraTrueType isBinaryAnalyser;
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

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isBinaryAnalyser</tt> and <tt>FunctorTraits::isUnaryFunctor</tt>
    are true (<tt>VigraTrueType<tt>)
    
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
    RegionStatistics::argument_type a;
    RegionStatistics::result_type r;

    region(a);     // update statistics
    r = region();  // return statistics

    \endcode
*/
template <class RegionStatistics, class LabelType = int>
class ArrayOfRegionStatistics
{
    typedef std::vector<RegionStatistics> RegionArray;

  public:
         /** argument type of the contained statistics object
             becomes first argument of the analyser
         */
    typedef typename RegionStatistics::argument_type first_argument_type;

         /** label type is used to determine the region to be updated
         */
    typedef LabelType second_argument_type;

         /** label type is also used to determine the region to be
             returned by the 1 argument operator()
         */
    typedef LabelType argument_type;

         /** result type of the contained statistics object
             becomes result type of the analyser
         */
    typedef typename RegionStatistics::result_type result_type;

         /** the value type of the array: the contained statistics object.
             <b>Note:</b> this definition was different in older
             VIGRA versions. The old definition was wrong.
         */
    typedef RegionStatistics value_type;

         /** the array's reference type
         */
    typedef RegionStatistics & reference;

         /** the array's const reference type
         */
    typedef RegionStatistics const & const_reference;

         /** type to iterate over the statistics array
         */
    typedef typename RegionArray::iterator iterator;

         /** type to iterate over a const statistics array
         */
    typedef typename RegionArray::const_iterator const_iterator;

        /** init array of RegionStatistics with default size 0.
        */
    ArrayOfRegionStatistics()
    {}

        /** init array of RegionStatistics with index domain
            0...max_region_label.
        */
    ArrayOfRegionStatistics(unsigned int max_region_label)
    : regions(max_region_label+1)
    {}

        /** resize array to new index domain 0...max_region_label.
            All bin are re-initialized.
        */
    void resize(unsigned int max_region_label)
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

        /** update regions statistics for region <TT>label</TT>. The label type
            is converted to <TT>unsigned int</TT>.
        */
    void operator()(first_argument_type const & v, second_argument_type label) {
        regions[static_cast<unsigned int>(label)](v);
    }

        /** merge second region into first
        */
    void merge(argument_type label1, argument_type label2) {
        regions[static_cast<unsigned int>(label1)](regions[static_cast<unsigned int>(label2)]);
    }

        /** ask for maximal index (label) allowed
        */
    unsigned int maxRegionLabel() const
        { return size() - 1; }

        /** ask for array size (i.e. maxRegionLabel() + 1)
        */
    unsigned int size() const
        { return regions.size(); }

        /** access the statistics for a region via its label. The label type
            is converted to <TT>unsigned int</TT>.
        */
    result_type operator()(argument_type label) const
        { return regions[static_cast<unsigned int>(label)](); }

        /** read the statistics functor for a region via its label
        */
    const_reference operator[](argument_type label) const
        { return regions[static_cast<unsigned int>(label)]; }

        /** access the statistics functor for a region via its label
        */
    reference operator[](argument_type label)
        { return regions[static_cast<unsigned int>(label)]; }

        /** iterator to the begin of the region array
        */
    iterator begin()
        { return regions.begin(); }

        /** const iterator to the begin of the region array
        */
    const_iterator begin() const
        { return regions.begin(); }

        /** iterator to the end of the region array
        */
    iterator end()
        { return regions.end(); }

        /** const iterator to the end of the region array
        */
    const_iterator end() const
        { return regions.end(); }

  private:
    std::vector<RegionStatistics> regions;
};

template <class RegionStatistics, class LabelType>
class FunctorTraits<ArrayOfRegionStatistics<RegionStatistics, LabelType> >
: public FunctorTraitsBase<ArrayOfRegionStatistics<RegionStatistics, LabelType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
    typedef VigraTrueType isBinaryAnalyser;
};

//@}

} // namespace vigra

#endif // VIGRA_INSPECTIMAGE_HXX
