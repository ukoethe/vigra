/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_INITIMAGE_HXX
#define VIGRA_INITIMAGE_HXX

#include "vigra/utilities.hxx"

namespace vigra {

/** @name Functions to Initialize Images
    @memo init images or image borders
*/
//@{

/********************************************************/
/*                                                      */
/*                        initImage                     */
/*                                                      */
/********************************************************/

/** Write a value to every pixel in an image or rectangular ROI.
    This function can be used to init the image.
    It uses an accessor to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void
        initImage(ImageIterator upperleft, ImageIterator lowerright, 
              Accessor a, VALUETYPE v)
    }
    \end{verbatim}

    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void
        initImage(triple<ImageIterator, ImageIterator, Accessor> img, VALUETYPE v)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/initimage.hxx]{../include/vigra/initimage.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage img(100, 100);
    
    // zero the image
    vigra::initImage(destImageRange(img),
                     vigra::NumericTraits<vigra::BImage::PixelType>::zero());
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator upperleft, lowerright;
    Accessor accessor;
    VALUETYPE v;
    
    accessor.set(v, upperleft); 
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator, class Accessor, class VALUETYPE>
void
initImage(ImageIterator upperleft, ImageIterator lowerright, 
          Accessor a,  VALUETYPE v)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    for(int y=0; y<h; ++y, ++upperleft.y)
    {
        ImageIterator ix(upperleft);
    
        for(int x=0; x<w; ++x, ++ix.x)
    {
        a.set(v, ix);
    }
    }
}
    
template <class ImageIterator, class Accessor, class VALUETYPE>
inline 
void
initImage(triple<ImageIterator, ImageIterator, Accessor> img, VALUETYPE v)
{
    initImage(img.first, img.second, img.third, v);
}
    
/********************************************************/
/*                                                      */
/*                      initImageIf                     */
/*                                                      */
/********************************************************/

/** Write value to pixel in the image if mask is true.
    This function can be used to init a region-of-interest of the image.
    It uses an accessor to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, 
              class MaskImageIterator, class MaskAccessor,
              class VALUETYPE>
        void
        initImageIf(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
              MaskImageIterator mask_upperleft, MaskAccessor ma,
              VALUETYPE v)
    }
    \end{verbatim}    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, 
              class MaskImageIterator, class MaskAccessor,
              class VALUETYPE>
        void
        initImageIf(triple<ImageIterator, ImageIterator, Accessor> img, 
            pair<MaskImageIterator, MaskAccessor> mask,
            VALUETYPE v)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/initimage.hxx]{../include/vigra/initimage.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage img(100, 100);
    vigra::BImage mask(100, 100);
    
    // zero the ROI
    vigra::initImageIf(destImageRange(img), 
                maskImage(mask),
                vigra::NumericTraits<vigra::BImage::PixelType>::zero());
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator upperleft, lowerright;
    MaskImageIterator mask_upperleft;
    Accessor accessor;
    MaskAccessor mask_accessor;
    VALUETYPE v;
    
    if(mask_accessor(mask_upperleft)) accessor.set(v, upperleft); 
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator, class Accessor, 
          class MaskImageIterator, class MaskAccessor,
          class VALUETYPE>
void
initImageIf(ImageIterator upperleft, ImageIterator lowerright, Accessor a,
          MaskImageIterator mask_upperleft, MaskAccessor ma,
          VALUETYPE v)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
        
    for(int y=0; y<h; ++y, ++upperleft.y, ++mask_upperleft.y)
    {
        ImageIterator ix(upperleft);
        MaskImageIterator mx(mask_upperleft);
    
        for(int x=0; x<w; ++x, ++ix.x, ++mx.x)
        {
            if(ma(mx)) a.set(v, ix);
        }
    }
}
    
template <class ImageIterator, class Accessor, 
          class MaskImageIterator, class MaskAccessor,
          class VALUETYPE>
inline 
void
initImageIf(triple<ImageIterator, ImageIterator, Accessor> img, 
            pair<MaskImageIterator, MaskAccessor> mask,
            VALUETYPE v)
{
    initImageIf(img.first, img.second, img.third, mask.first, mask.second, v);
}
    
/********************************************************/
/*                                                      */
/*                    initImageBorder                   */
/*                                                      */
/********************************************************/

/** Write value to the specified border pixels in the image.
    A pixel is initialized if its distance to the border 
    is at most 'borderwidth'.
    It uses an accessor to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void
        initImageBorder(ImageIterator upperleft, ImageIterator lowerright, 
                Accessor a,  int border_width, VALUETYPE v)
    }
    \end{verbatim}

    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class ImageIterator, class Accessor, class VALUETYPE>
        void
        initImageBorder(triple<ImageIterator, ImageIterator, Accessor> img, 
                int border_width, VALUETYPE v)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/initimage.hxx]{../include/vigra/initimage.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage img(100, 100);
    
    // zero a border of 5 pixel
    vigra::initImageBorder(destImageRange(img),
                    5, vigra::NumericTraits<vigra::BImage::PixelType>::zero());
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator upperleft, lowerright;
    Accessor accessor;
    VALUETYPE v;
    
    accessor.set(v, upperleft); 
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator, class Accessor, class VALUETYPE>
inline 
void
initImageBorder(ImageIterator upperleft, ImageIterator lowerright, 
                Accessor a,  int border_width, VALUETYPE v)
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
inline 
void
initImageBorder(triple<ImageIterator, ImageIterator, Accessor> img, 
                int border_width, VALUETYPE v)
{
    initImageBorder(img.first, img.second, img.third, border_width, v);
}
    
//@}


} // namespace vigra

#endif // VIGRA_INITIMAGE_HXX
