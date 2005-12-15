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
 
 
#ifndef VIGRA_COPYIMAGE_HXX
#define VIGRA_COPYIMAGE_HXX

#include "vigra/utilities.hxx"

namespace vigra {

/** \addtogroup CopyAlgo Algorithms to Copy Images
    Copy images or regions
*/
//@{

/********************************************************/
/*                                                      */
/*                       copyLine                       */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
copyLine(SrcIterator s, 
         SrcIterator send, SrcAccessor src,
         DestIterator d, DestAccessor dest)
{
    for(; s != send; ++s, ++d)
        dest.set(src(s), d);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor, 
          class DestIterator, class DestAccessor>
void
copyLineIf(SrcIterator s, 
           SrcIterator send, SrcAccessor src,
           MaskIterator m, MaskAccessor mask,
           DestIterator d, DestAccessor dest)
{
    for(; s != send; ++s, ++d, ++m)
        if(mask(m))
            dest.set(src(s), d);
}

/********************************************************/
/*                                                      */
/*                        copyImage                     */
/*                                                      */
/********************************************************/

/** \brief Copy source image into destination image.

    If necessary, type conversion takes place.
    The function uses accessors to access the pixel data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void
        copyImage(SrcImageIterator src_upperleft, 
              SrcImageIterator src_lowerright, SrcAccessor sa,
              DestImageIterator dest_upperleft, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        void
        copyImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              pair<DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="copyimage_8hxx-source.html">vigra/copyimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::copyImage(srcImageRange(src), destImage(dest));
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator      dest_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    dest_accessor.set(src_accessor(sx), dx);

    \endcode
    
*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void
copyImage(SrcImageIterator src_upperleft, 
          SrcImageIterator src_lowerright, SrcAccessor sa,
          DestImageIterator dest_upperleft, DestAccessor da)
{
    int w = src_lowerright.x - src_upperleft.x;
    
    for(; src_upperleft.y<src_lowerright.y; ++src_upperleft.y, ++dest_upperleft.y)
    {
        copyLine(src_upperleft.rowIterator(), 
                 src_upperleft.rowIterator() + w, sa, 
                 dest_upperleft.rowIterator(), da);
    }
}
    
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline
void
copyImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
      pair<DestImageIterator, DestAccessor> dest)
{
    copyImage(src.first, src.second, src.third, 
                   dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                       copyImageIf                    */
/*                                                      */
/********************************************************/

/** \brief Copy source ROI into destination image.

    Pixel values are copied whenever the return value of the mask's
    accessor is not zero.
    If necessary, type conversion takes place.
    The function uses accessors to access the pixel data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor>
        void
        copyImageIf(SrcImageIterator src_upperleft, 
            SrcImageIterator src_lowerright, SrcAccessor sa,
            MaskImageIterator mask_upperleft, MaskAccessor ma,
            DestImageIterator dest_upperleft, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor>
        void
        copyImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<MaskImageIterator, MaskAccessor> mask,
            pair<DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="copyimage_8hxx-source.html">vigra/copyimage.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::copyImageIf(srcImageRange(src), maskImage(mask), destImage(dest));

    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    MaskImageIterator mask_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    MaskAccessor mask_accessor;
    Functor functor;
    
    if(mask_accessor(mx))
        dest_accessor.set(src_accessor(sx), dx);

    \endcode
    
*/
template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor>
void
copyImageIf(SrcImageIterator src_upperleft, 
            SrcImageIterator src_lowerright, SrcAccessor sa,
            MaskImageIterator mask_upperleft, MaskAccessor ma,
        DestImageIterator dest_upperleft, DestAccessor da)
{
    int w = src_lowerright.x - src_upperleft.x;
    
    for(; src_upperleft.y<src_lowerright.y; 
             ++src_upperleft.y, ++mask_upperleft.y, ++dest_upperleft.y)
    {
        copyLineIf(src_upperleft.rowIterator(), 
                   src_upperleft.rowIterator() + w, sa, 
                   mask_upperleft.rowIterator(), ma, 
                   dest_upperleft.rowIterator(), da);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor>
inline
void
copyImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        pair<MaskImageIterator, MaskAccessor> mask,
        pair<DestImageIterator, DestAccessor> dest)
{
    copyImageIf(src.first, src.second, src.third, 
                mask.first, mask.second, 
        dest.first, dest.second);
}

//@}

} // namespace vigra

#endif // VIGRA_COPYIMAGE_HXX
