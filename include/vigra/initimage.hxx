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
 
 
#ifndef VIGRA_INITIMAGE_HXX
#define VIGRA_INITIMAGE_HXX

#include "utilities.hxx"
#include "iteratortraits.hxx"
#include "functortraits.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup InitAlgo Algorithms to Initialize Images
    
    Init images or image borders
*/
//@{

/********************************************************/
/*                                                      */
/*                       initLine                       */
/*                                                      */
/********************************************************/

template <class DestIterator, class DestAccessor, class VALUETYPE>
inline void
initLineImpl(DestIterator d, DestIterator dend, DestAccessor dest,
             VALUETYPE const & v, VigraFalseType)
{
    for(; d != dend; ++d)
        dest.set(v, d);
}

template <class DestIterator, class DestAccessor, class FUNCTOR>
inline void
initLineImpl(DestIterator d, DestIterator dend, DestAccessor dest,
             FUNCTOR const & f, VigraTrueType)
{
    for(; d != dend; ++d)
        dest.set(f(), d);
}

template <class DestIterator, class DestAccessor, class VALUETYPE>
inline void
initLine(DestIterator d, DestIterator dend, DestAccessor dest,
         VALUETYPE const & v)
{
    initLineImpl(d, dend, dest, v, typename FunctorTraits<VALUETYPE>::isInitializer());
}

template <class DestIterator, class DestAccessor, class FUNCTOR>
inline void
initLineFunctor(DestIterator d, DestIterator dend, DestAccessor dest,
         FUNCTOR & f)
{
    for(; d != dend; ++d)
        dest.set(f(), d);
}

template <class DestIterator, class DestAccessor, 
          class MaskIterator, class MaskAccessor, 
          class VALUETYPE>
inline void
initLineIfImpl(DestIterator d, DestIterator dend, DestAccessor dest,
               MaskIterator m, MaskAccessor mask,
               VALUETYPE const & v, VigraFalseType)
{
    for(; d != dend; ++d, ++m)
        if(mask(m))
            dest.set(v, d);
}

template <class DestIterator, class DestAccessor, 
          class MaskIterator, class MaskAccessor, 
          class FUNCTOR>
inline void
initLineIfImpl(DestIterator d, DestIterator dend, DestAccessor dest,
               MaskIterator m, MaskAccessor mask,
               FUNCTOR const & f, VigraTrueType)
{
    for(; d != dend; ++d, ++m)
        if(mask(m))
            dest.set(f(), d);
}

template <class DestIterator, class DestAccessor, 
          class MaskIterator, class MaskAccessor, 
          class VALUETYPE>
inline void
initLineIf(DestIterator d, DestIterator dend, DestAccessor dest,
           MaskIterator m, MaskAccessor mask,
           VALUETYPE const & v)
{
    initLineIfImpl(d, dend, dest, m, mask, v, typename FunctorTraits<VALUETYPE>::isInitializer());
}

template <class DestIterator, class DestAccessor, 
          class MaskIterator, class MaskAccessor, 
          class FUNCTOR>
inline void
initLineFunctorIf(DestIterator d, DestIterator dend, DestAccessor dest,
                  MaskIterator m, MaskAccessor mask,
                  FUNCTOR & f)
{
    for(; d != dend; ++d, ++m)
        if(mask(m))
            dest.set(f(), d);
}

/********************************************************/
/*                                                      */
/*                        initImage                     */
/*                                                      */
/********************************************************/

/** \brief Write a value to every pixel in an image or rectangular ROI.

    This function can be used to init the image.
    
    The initial value can either be a constant of appropriate type (compatible with 
    the destination's value_type), or a functor with compatible result_type. These two 
    cases are automatically distinguished when <tt>FunctorTraits<FUNCTOR>::isInitializer</tt>
    yields <tt>VigraTrueType</tt>. Since the functor is passed by <tt>const</tt> reference, its 
    <tt>operator()</tt> must be const, and its internal state may need to be <tt>mutable</tt>.
    
    Function \ref initMultiArray() implements the same functionality for arbitrary dimensional
    arrays. In many situations, the assignment functions of \ref vigra::MultiArrayView offer
    a simpler and more readable alternative to the init functions.

    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class VALUETYPE>
        void
        initImage(MultiArrayView<2, T, S> img, VALUETYPE const & v);
        
        template <class T, class S, class FUNCTOR>
        void
        initImage(MultiArrayView<2, T, S> img, FUNCTOR const & v);
    }
    \endcode
    
    \deprecatedAPI{initImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void initImage(ImageIterator upperleft, ImageIterator lowerright, 
                       Accessor a, VALUETYPE const & v);

        template <class ImageIterator, class Accessor, class FUNCTOR>
        void initImage(ImageIterator upperleft, ImageIterator lowerright, 
                       Accessor a, FUNCTOR const & v);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void initImage(triple<ImageIterator, ImageIterator, Accessor> img, VALUETYPE const & v);

        template <class ImageIterator, class Accessor, class FUNCTOR>
        void initImage(triple<ImageIterator, ImageIterator, Accessor> img, FUNCTOR const & v);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/initimage.hxx\><br>
    Namespace: vigra
    
    Initialize with a constant:
    \code
    MultiArray<2, unsigned char> img(100, 100);
    
    // init the image with the value 128
    initImage(img, 128);
    
    // init the interior with the value 1
    initImage(img.subarray(Shape2(10), Shape2(-10)), 1);
    
    // equivalent to
    img = 128;
    img.init(128);
    img.subarray(Shape2(10), Shape2(-10)) = 1;
    \endcode

    Initialize with a functor:
    \code
    struct Counter {
        Counter() : count(0) {}
        
        int operator()() const { return count++; }
    
        mutable int count;
    };
    
    MultiArray<2, int> img(100, 100);
        
    // write the current count in every pixel
    initImage(img, Counter());
    
    // equivalent to
    #include <vigra/algorithm.hxx>
    
    linearSequence(img.begin(), img.end());
    \endcode

    \deprecatedUsage{initImage}
    \code
    vigra::BImage img(100, 100);
    
    // init the image with the value 128
    vigra::initImage(destImageRange(img), 128);

    // Initialize with a functor:
    struct Counter {
        Counter() : count(0) {}
        
        int operator()() const { return count++; }
    
        mutable int count;
    };
    
    // write the current count in every pixel
    vigra::initImage(destImageRange(img), Counter());
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator upperleft, lowerright;
    ImageIterator::row_iterator ix = upperleft.rowIterator();
    
    Accessor accessor;
    VALUETYPE v;
    
    accessor.set(v, ix); 
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void initImage)

template <class ImageIterator, class Accessor, class VALUETYPE>
void
initImage(ImageIterator upperleft, ImageIterator lowerright, 
          Accessor a,  VALUETYPE const & v)
{
    int w = lowerright.x - upperleft.x;
    
    for(; upperleft.y < lowerright.y; ++upperleft.y)
    {
        initLineImpl(upperleft.rowIterator(), upperleft.rowIterator() + w, a, 
                     v, typename FunctorTraits<VALUETYPE>::isInitializer());
    }
}
    
template <class ImageIterator, class Accessor, class VALUETYPE>
inline void
initImage(triple<ImageIterator, ImageIterator, Accessor> img, VALUETYPE const & v)
{
    initImage(img.first, img.second, img.third, v);
}
    
template <class T, class S, class VALUETYPE>
inline void
initImage(MultiArrayView<2, T, S> img, VALUETYPE const & v)
{
    initImage(destImageRange(img), v);
}
    
/********************************************************/
/*                                                      */
/*                 initImageWithFunctor                 */
/*                                                      */
/********************************************************/

/** \brief Write the result of a functor call to every pixel in an image or rectangular ROI.

    This function can be used to init the image by calling the given 
    functor for each pixel. The functor is 
    passed by reference, so that its internal state can be updated in each call.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class FUNCTOR>
        void
        initImageWithFunctor(MultiArrayView<2, T, S> img, FUNCTOR & f);
    }
    \endcode
    
    \deprecatedAPI{initImageWithFunctor}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class FUNCTOR>
        void
        initImageWithFunctor(ImageIterator upperleft, ImageIterator lowerright, 
                  Accessor a,  FUNCTOR & f);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class FUNCTOR>
        void
        initImageWithFunctor(triple<ImageIterator, ImageIterator, Accessor> img, FUNCTOR & f);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/initimage.hxx\><br>
    Namespace: vigra

    \code
    struct Counter {
        Counter() : count(0) {}
        
        int operator()() const { return count++; }
    
        int count;
    };
    
    MultiArray<2, int> img(100, 100);
        
    // write the current count in every pixel
    Counter counter;
    initImageWithFunctor(img, counter);
    
    // equivalent to
    #include <vigra/algorithm.hxx>
    
    linearSequence(img.begin(), img.end());
    \endcode

    \deprecatedUsage{initImageWithFunctor}
    \code
    struct Counter {
        Counter() : count(0) {}
        
        int operator()() const { return count++; }
    
        mutable int count;
    };
    
    vigra::IImage img(100, 100);
        
    // write the current count in every pixel
    Counter counter;
    vigra::initImageWithFunctor(destImageRange(img), counter);
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator upperleft, lowerright;
    ImageIterator::row_iterator ix = upperleft.rowIterator();
    
    Accessor accessor;
    Functor f;
    
    accessor.set(f(), ix); 
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void initImageWithFunctor)

template <class ImageIterator, class Accessor, class FUNCTOR>
void
initImageWithFunctor(ImageIterator upperleft, ImageIterator lowerright, 
          Accessor a,  FUNCTOR & f)
{
    int w = lowerright.x - upperleft.x;
    
    for(; upperleft.y < lowerright.y; ++upperleft.y)
    {
        initLineFunctor(upperleft.rowIterator(), upperleft.rowIterator() + w, a, f);
    }
}
    
template <class ImageIterator, class Accessor, class FUNCTOR>
inline void
initImageWithFunctor(triple<ImageIterator, ImageIterator, Accessor> img, FUNCTOR & f)
{
    initImageWithFunctor(img.first, img.second, img.third, f);
}
    
template <class T, class S, class FUNCTOR>
inline void
initImageWithFunctor(MultiArrayView<2, T, S> img, FUNCTOR & f)
{
    initImageWithFunctor(destImageRange(img), f);
}
    
/********************************************************/
/*                                                      */
/*                      initImageIf                     */
/*                                                      */
/********************************************************/

/** \brief Write value to pixel in the image if mask is true.

    This function can be used to init a region-of-interest of the image.
    
    The initial value can either be a constant of appropriate type (compatible with 
    the destination's value_type), or a functor with compatible result_type. These two 
    cases are automatically distinguished when <tt>FunctorTraits<FUNCTOR>::isInitializer</tt>
    yields <tt>VigraTrueType</tt>. Since the functor is passed by <tt>const</tt> reference, its 
    <tt>operator()</tt> must be const, and its internal state may need to be <tt>mutable</tt>.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, 
                  class TM, class SM,
                  class VALUETYPE>
        void
        initImageIf(MultiArrayView<2, T, S> img, 
                    MultiArrayView<2, TM, SM> const & mask,
                    VALUETYPE const & v);
                    
        template <class T, class S, 
                  class TM, class SM,
                  class FUNCTOR>
        void
        initImageIf(MultiArrayView<2, T, S> img, 
                    MultiArrayView<2, TM, SM> const & mask,
                    FUNCTOR const & v);
    }
    \endcode     
    \deprecatedAPI{initImageIf}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, 
                  class MaskImageIterator, class MaskAccessor,
                  class VALUETYPE>
        void initImageIf(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
                         MaskImageIterator mask_upperleft, MaskAccessor ma,
                         VALUETYPE const & v);

        template <class ImageIterator, class Accessor, 
                  class MaskImageIterator, class MaskAccessor,
                  class FUNCTOR>
        void initImageIf(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
                         MaskImageIterator mask_upperleft, MaskAccessor ma,
                         FUNCTOR const & v);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, 
                  class MaskImageIterator, class MaskAccessor,
                  class VALUETYPE>
        void initImageIf(triple<ImageIterator, ImageIterator, Accessor> img, 
                         pair<MaskImageIterator, MaskAccessor> mask,
                         VALUETYPE const & v);

        template <class ImageIterator, class Accessor, 
                  class MaskImageIterator, class MaskAccessor,
                  class FUNCTOR>
        void initImageIf(triple<ImageIterator, ImageIterator, Accessor> img, 
                         pair<MaskImageIterator, MaskAccessor> mask,
                         FUNCTOR const & v);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/initimage.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, RGBValue<unsigned char> >  img(100, 100),
    MultiArray<2, unsigned char>             mask(100, 100);
    ... // init the ROI mask
    
    // set the ROI to one
    initImageIf(img, mask,
                NumericTraits<RGBValue<unsigned char> >::one());
    \endcode

    \deprecatedUsage{initImageIf}
    \code
    vigra::BImage img(100, 100);
    vigra::BImage mask(100, 100);
    
    // zero the ROI
    vigra::initImageIf(destImageRange(img), 
                maskImage(mask),
                vigra::NumericTraits<vigra::BImage::PixelType>::zero());
    \endcode
    <b> Required Interface:</b>
    \code
    ImageIterator upperleft, lowerright;
    MaskImageIterator mask_upperleft;
    ImageIterator::row_iterator ix = upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();
    
    Accessor accessor;
    MaskAccessor mask_accessor;
    VALUETYPE v;
    
    if(mask_accessor(mx)) accessor.set(v, ix); 
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void initImageIf)

template <class ImageIterator, class Accessor, 
          class MaskImageIterator, class MaskAccessor,
          class VALUETYPE>
void
initImageIf(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
          MaskImageIterator mask_upperleft, MaskAccessor ma,
          VALUETYPE const & v)
{
    int w = lowerright.x - upperleft.x;
        
    for(; upperleft.y < lowerright.y; ++upperleft.y, ++mask_upperleft.y)
    {
        initLineIfImpl(upperleft.rowIterator(), 
                   upperleft.rowIterator() + w, a, 
                   mask_upperleft.rowIterator(), ma, 
                   v, typename FunctorTraits<VALUETYPE>::isInitializer());
    }
}
    
template <class ImageIterator, class Accessor, 
          class MaskImageIterator, class MaskAccessor,
          class VALUETYPE>
inline void
initImageIf(triple<ImageIterator, ImageIterator, Accessor> img, 
            pair<MaskImageIterator, MaskAccessor> mask,
            VALUETYPE const & v)
{
    initImageIf(img.first, img.second, img.third, mask.first, mask.second, v);
}
    
template <class T, class S, 
          class TM, class SM,
          class VALUETYPE>
inline void
initImageIf(MultiArrayView<2, T, S> img, 
            MultiArrayView<2, TM, SM> const & mask,
            VALUETYPE const & v)
{
    vigra_precondition(img.shape() == mask.shape(),
        "initImageIf(): shape mismatch between input and mask.");
    initImageIf(destImageRange(img), maskImage(mask), v);
}
    
/********************************************************/
/*                                                      */
/*                    initImageBorder                   */
/*                                                      */
/********************************************************/

/** \brief Write value to the specified border pixels in the image.

    A pixel is initialized if its distance to the border 
    is at most 'borderwidth'. 
    
    The initial value can either be a constant of appropriate type (compatible with 
    the destination's value_type), or a functor with compatible result_type. These two 
    cases are automatically distinguished when <tt>FunctorTraits<FUNCTOR>::isInitializer</tt>
    yields <tt>VigraTrueType</tt>. Since the functor is passed by <tt>const</tt> reference, its 
    <tt>operator()</tt> must be const, and its internal state may need to be <tt>mutable</tt>.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class VALUETYPE>
        void
        initImageBorder(MultiArrayView<2, T, S> img, 
                        int border_width, VALUETYPE const & v);
                        
        template <class T, class S, class FUNCTOR>
        void
        initImageBorder(MultiArrayView<2, T, S> img, 
                        int border_width, FUNCTOR const & v);
    }
    \endcode
    
    \deprecatedAPI{initImageBorder}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void initImageBorder(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
                             int border_width, VALUETYPE const & v);

        template <class ImageIterator, class Accessor, class FUNCTOR>
        void initImageBorder(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
                             int border_width, FUNCTOR const & v);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void initImageBorder(triple<ImageIterator, ImageIterator, Accessor> img, 
                             int border_width, VALUETYPE const & v);

        template <class ImageIterator, class Accessor, class FUNCTOR>
        void initImageBorder(triple<ImageIterator, ImageIterator, Accessor> img, 
                             int border_width, FUNCTOR const & v);
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/initimage.hxx\><br>
    Namespace: vigra
    
    \code
    #include <vigra/random.hxx>
    
    MultiArray<2, int> img(100, 100);
    
    // fill a border of 5 pixels with random numbers
    initImageBorder(img, 5, MersenneTwister());
    \endcode

    \deprecatedUsage{initImageBorder}
    \code
    vigra::BImage img(100, 100);
    
    // zero a border of 5 pixel
    vigra::initImageBorder(destImageRange(img),
                    5, vigra::NumericTraits<vigra::BImage::PixelType>::zero());
    \endcode
    <b> Required Interface:</b>
    <br/>see \ref initImage()
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void initImageBorder)

template <class ImageIterator, class Accessor, class VALUETYPE>
inline 
void
initImageBorder(ImageIterator upperleft, ImageIterator lowerright, 
                Accessor a,  int border_width, VALUETYPE const & v)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    int hb = (border_width > h) ? h : border_width;
    int wb = (border_width > w) ? w : border_width;
    
    initImage(upperleft, upperleft+Diff2D(w,hb), a, v);
    initImage(upperleft, upperleft+Diff2D(wb,h), a, v);
    initImage(upperleft+Diff2D(0,h-hb), lowerright, a, v);
    initImage(upperleft+Diff2D(w-wb,0), lowerright, a, v);
}
    
template <class ImageIterator, class Accessor, class VALUETYPE>
inline void
initImageBorder(triple<ImageIterator, ImageIterator, Accessor> img, 
                int border_width, VALUETYPE const & v)
{
    initImageBorder(img.first, img.second, img.third, border_width, v);
}
    
template <class T, class S, class VALUETYPE>
inline void
initImageBorder(MultiArrayView<2, T, S> img, 
                int border_width, VALUETYPE const & v)
{
    initImageBorder(destImageRange(img), border_width, v);
}
    
//@}


} // namespace vigra

#endif // VIGRA_INITIMAGE_HXX
