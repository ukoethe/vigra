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
 
 
#ifndef VIGRA_TRANSFORMIMAGE_HXX
#define VIGRA_TRANSFORMIMAGE_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"

/** @name Functions to Transform Images
    @memo apply functor to calculate a pixelwise transformation of one image
*/
//@{

/********************************************************/
/*                                                      */
/*                      transformImage                  */
/*                                                      */
/********************************************************/

/** Transform source image into destination image.
    The transformation given by the functor is applied to every source 
    pixel and the result written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
    void
    transformImage(SrcImageIterator src_upperleft, 
           SrcImageIterator src_lowerright, SrcAccessor sa,
           DestImageIterator dest_upperleft, DestAccessor da,
           Functor f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
    void
    transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
           pair<DestImageIterator, DestAccessor> dest,
           Functor f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/transformimage.hxx]{../include/vigra/transformimage.hxx}
    
    \begin{verbatim}
    
    #include <math.h>         // for sqrt()
    #include <function.h>     // for ptr_fun()
    
    transformImage(srcImageRange(src), 
                   destImage(dest), 
           ptr_fun(sqrt) );
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator      dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(functor(src_accessor(src_upperleft)), dest_upperleft);

    \end{verbatim}
    
    @memo
*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
void
transformImage(SrcImageIterator src_upperleft, 
               SrcImageIterator src_lowerright, SrcAccessor sa,
               DestImageIterator dest_upperleft, DestAccessor da,
           Functor f)
{
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;
    
    for(int y=0; y<h; ++y, ++src_upperleft.y, ++dest_upperleft.y)
    {
        SrcImageIterator six(src_upperleft);
        DestImageIterator dix(dest_upperleft);

        for(int x=0; x<w; ++x, ++six.x, ++dix.x)
        {
            da.set(static_cast<typename 
                 DestAccessor::value_type>(f(sa(six))), dix);
        }
    }
}
    
template <class SrcImageIterator, class SrcAccessor,
      class DestImageIterator, class DestAccessor, class Functor>
inline
void
transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
           pair<DestImageIterator, DestAccessor> dest,
           Functor f)
{
    transformImage(src.first, src.second, src.third, 
                   dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                   transformImageIf                   */
/*                                                      */
/********************************************************/

/** Transform source ROI into destination ROI.
    The transformation given by the functor is applied to every source 
    pixel in the ROI (i.e. when the return vlaue of the mask's accessor
    is not zero)
    and the result is written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, clas DestAccessor,
          class Functor>
    void
    transformImageIf(SrcImageIterator src_upperleft, 
        SrcImageIterator src_lowerright, SrcAccessor sa,
        MaskImageIterator mask_upperleft, MaskAccessor ma,
        DestImageIterator dest_upperleft, DestAccessor da,
        Functor f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, clas DestAccessor,
          class Functor>
    void
    transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
             pair<MaskImageIterator, MaskAccessor> mask,
             pair<DestImageIterator, DestAccessor> dest,
             Functor f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/transformimage.hxx]{../include/vigra/transformimage.hxx}
    
    \begin{verbatim}
    #include <math.h>         // for sqrt()
    #include <function.h>     // for ptr_fun()
    
    transformImageIf(srcImageRange(src), 
                     maskImage(mask), 
                     destImage(dest), 
             ptr_fun(sqrt) );

    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator  dest_upperleft;
    MaskImageIterator mask_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    MaskAccessor mask_accessor;
    Functor functor;
    
    if(mask_accessor(mask_uppeleft))
       dest_accessor.set(functor(src_accessor(src_upperleft)), dest_upperleft);

    \end{verbatim}
    
    @memo
*/
template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
      class Functor>
void
transformImageIf(SrcImageIterator src_upperleft, 
            SrcImageIterator src_lowerright, SrcAccessor sa,
            MaskImageIterator mask_upperleft, MaskAccessor ma,
        DestImageIterator dest_upperleft, DestAccessor da,
        Functor f)
{
    int w = src_lowerright.x - src_upperleft.x;
    int h = src_lowerright.y - src_upperleft.y;
    
    for(int y=0; y<h; ++y, 
             ++src_upperleft.y, ++mask_upperleft.y, ++dest_upperleft.y)
    {
        SrcImageIterator six(src_upperleft);
        MaskImageIterator mx(mask_upperleft);
        DestImageIterator       dix(dest_upperleft);

        for(int x=0; x<w; ++x, ++six.x, ++mx.x, ++dix.x)
        {
            if(ma(mx)) da.set(static_cast<typename 
                DestAccessor::value_type>(f(sa(six))), dix);
        }
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
      class Functor>
inline
void
transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
             pair<MaskImageIterator, MaskAccessor> mask,
             pair<DestImageIterator, DestAccessor> dest,
             Functor f)
{
    transformImageIf(src.first, src.second, src.third, 
                     mask.first, mask.second, 
             dest.first, dest.second, f);
}

//@}

/** @name Functors to Transform Images

    Note that the unary functors of the STL can also be used.
    
    @memo frequently used pixel transformations
*/
//@{

template <class SrcValueType>
class LinearIntensityTransform
{
   public:
    /** a real type convertible from SrcValueType
        @memo
    */
    typedef typename 
            NumericTraits<SrcValueType>::RealPromote SrcPromote;
    
    /** init scale and offset
        @memo
    */
    LinearIntensityTransform(double scale, SrcPromote offset)
    : scale_(scale), offset_(offset)
    {}
    
    /** calculate transform
        @memo
    */
    SrcPromote operator()(SrcValueType s) const
    {
        return (scale_ * (s + offset_));
    }
    
  private:
  
    double scale_;
    SrcPromote offset_;
};

/********************************************************/
/*                                                      */
/*              linearIntensityTransform                */
/*                                                      */
/********************************************************/

/** Factory function for a functor that applies a linear transform to the 
    source pixel values. The functor applies the transform 
    '#destvalue = scale * (srcvalue + offset)#' to every pixel.
    This can, for example, be used to transform images into the visible 
    range 0...255 or to invert an image.
    
    {\bf Declaration:}
    
    template <class SrcValueType>
    LinearIntensityTransform<SrcValueType>
    linearIntensityTransform(double scale, SrcValueType offset)
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/transformimage.hxx]{../include/vigra/transformimage.hxx}
    
    \begin{verbatim}
    IImage src(width, height);
    BImage dest(width, height);
    ...
    FindMinMax<IImage::PixelType> minmax;   // functor to find range
    
    inspectImage(srcImageRange(src), minmax); // find original range
    
    // transform to range 0...255
    transformImage(srcImageRange(src), destImage(dest),
                   linearIntensityTransform(
                       255.0 / (minmax.max - minmax.min), // scaling
                       - minmax.min));                    // offset
     
    \end{verbatim}

    {\bf Required Interface:}
    
    The source value type must be a model of \Ref{LinearSpace}.
    
    @memo apply a linear transform to the source pixel values 
    
*/
template <class SrcValueType>
LinearIntensityTransform<SrcValueType>
linearIntensityTransform(double scale, SrcValueType offset)
{
    return LinearIntensityTransform<SrcValueType>(scale, offset);
}


/********************************************************/
/*                                                      */
/*                      Threshold                       */
/*                                                      */
/********************************************************/

/** Threshold an image.
    If a source pixel is above or equal the lower and below 
    or equal the higher threshold (i.e. within the closed interval
    [lower, heigher]) the destination pixel is set to 'yesresult',
    otherwise to 'noresult'.

    {\bf Usage:}
    
        Include-File:
        \URL[vigra/transformimage.hxx]{../include/vigra/transformimage.hxx}
    
    \begin{verbatim}
    BImage src(width, height), dest(width, height);
    ...
    transformImage(src.upperLeft(), src.lowerRight(), src.accessor(),
                   dest.upperLeft(), dest.accessor(),
     Threshold<BImage::PixelType, BImage::PixelType>(10, 100, 0, 255));
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    
    SrcValueType   src;
    DestValueType  dest, yesresult, noresult;
    
    dest = ((src < lower) || (higher < src)) ? noresult : yesresult;
    
    \end{verbatim}
    
*/
template <class SrcValueType, class DestValueType>
class Threshold
{
   public:
    /** init thresholds and resturn values
        @memo
    */
    Threshold(SrcValueType lower, SrcValueType higher,
              DestValueType noresult, DestValueType yesresult)
    : lower_(lower), higher_(higher),
      yesresult_(yesresult), noresult_(noresult)
    {}
    
    /** calculate transform
        @memo
    */
    DestValueType operator()(SrcValueType s) const
    {
        return ((s < lower_) || (higher_ < s)) ? noresult_ : yesresult_;
    }
    
  private:
  
    SrcValueType lower_, higher_;
    DestValueType yesresult_, noresult_;
};

//@}

#endif // VIGRA_TRANSFORMIMAGE_HXX
