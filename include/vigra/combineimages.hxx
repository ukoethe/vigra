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
 
 
#ifndef VIGRA_COMBINEIMAGES_HXX
#define VIGRA_COMBINEIMAGES_HXX

#include "vigra/utilities.hxx"

namespace vigra {

/** @heading Algorithms to Combine Images

    Note that the binary functors of the STL may be used with these functions.

    @memo apply functor to calculate a pixelwise transformation depending on multiple images
*/
//@{

/********************************************************/
/*                                                      */
/*                    combineTwoImages                  */
/*                                                      */
/********************************************************/

/** Combine two source images into destination image.
    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    The function uses accessors to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineTwoImages(SrcImageIterator1 src1_upperleft, 
                 SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                 SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                 DestImageIterator dest_upperleft, DestAccessor da,
                 Functor f)
    }
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineTwoImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor f)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/combineimages.hxx]{../include/vigra/combineimages.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    #include <functional>     // for plus

    vigra::combineTwoImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                destIter(dest.upperLeft()),  
                std::plus<SrcValueType>());
    
    \end{verbatim}
    
    Note that #SrcValueType# must be replaced with the appropriate type (e.g. 
    the promote type of the input images' pixel type, see also 
    \URef{NumericandPromotionTraits})
    
    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    DestImageIterator dest_upperleft;
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(src1_upperleft), src2_accessor(src2_upperleft)), 
      dest_upperleft);

    \end{verbatim}
    
    @memo
    
*/
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
      class DestImageIterator, class DestAccessor,
          class Functor>
void
combineTwoImages(SrcImageIterator1 src1_upperleft, 
                 SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                 SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                 DestImageIterator dest_upperleft, DestAccessor da,
         Functor f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    int h = src1_lowerright.y - src1_upperleft.y;
    
    for(int y=0; y<h; 
            ++y, ++src1_upperleft.y, ++src2_upperleft.y, ++dest_upperleft.y)
    {
        SrcImageIterator1 ix1(src1_upperleft);
        SrcImageIterator2 ix2(src2_upperleft);
        DestImageIterator dix(dest_upperleft);

        for(int x=0; x<w; ++x, ++ix1.x, ++ix2.x, ++dix.x)
        {
            da.set(static_cast<typename 
           DestAccessor::value_type>(f(sa1(ix1), sa2(ix2))), dix);
        }
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
      class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
combineTwoImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
             pair<SrcImageIterator2, SrcAccessor2> src2,
             pair<DestImageIterator, DestAccessor> dest,
         Functor f)
{
    combineTwoImages(src1.first, src1.second, src1.third, 
                     src2.first, src2.second, 
             dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                  combineTwoImagesIf                  */
/*                                                      */
/********************************************************/

/** Combine ROI of two source images into destination image.
    The transformation given by the functor is applied to all source 
    pixels in the ROI (i.e. whenever the return value of the mask's accessor
    is not zero)
    and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    The function uses accessors to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        combineTwoImagesIf(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   MaskImageIterator mask_upperleft, MaskAccessor ma,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor f)
    }
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        combineTwoImagesIf(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                   pair<SrcImageIterator2, SrcAccessor2> src2,
                   pair<MaskImageIterator, MaskAccessor> mask,
                   pair<DestImageIterator, DestAccessor> dest,
                   Functor f)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/combineimages.hxx]{../include/vigra/combineimages.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    #include <functional>     // for plus

    vigra::combineTwoImagesIf(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                maskIter(mask.upperLeft()), 
                destIter(dest.upperLeft()),
                std::plus<SrcValueType>());
    
    \end{verbatim}

    Note that #SrcValueType# must be replaced with the appropriate type (e.g. 
    the promote type of the input images' pixel type, see also 
    \URef{NumericandPromotionTraits})
    
    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    MaskImageIterator mask_upperleft;
    DestImageIterator dest_upperleft;
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    MaskAccessor mask_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;
    
    if(mask_accessor(mask_upperleft))
       dest_accessor.set(
          functor(src1_accessor(src1_upperleft), src2_accessor(src2_upperleft)), 
      dest_upperleft);

    \end{verbatim}
    
    @memo
*/
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
      class Functor>
void
combineTwoImagesIf(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   MaskImageIterator mask_upperleft, MaskAccessor ma,
               DestImageIterator dest_upperleft, DestAccessor da,
               Functor f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    int h = src1_lowerright.y - src1_upperleft.y;
    
    for(int y=0; y<h; ++y, ++src1_upperleft.y, ++src2_upperleft.y, 
                           ++dest_upperleft.y, ++mask_upperleft.y)
    {
        SrcImageIterator1 ix1(src1_upperleft);
        SrcImageIterator2 ix2(src2_upperleft);
        MaskImageIterator mix(mask_upperleft);
        DestImageIterator dix(dest_upperleft);

        for(int x=0; x<w; ++x, ++ix1.x, ++ix2.x, ++dix.x, ++mix.x)
        {
            if(ma(mix)) da.set(static_cast<typename 
             DestAccessor::value_type>(f(sa1(ix1), sa2(ix2))), dix);
        }
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
      class Functor>
inline
void
combineTwoImagesIf(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
               pair<SrcImageIterator2, SrcAccessor2> src2,
               pair<MaskImageIterator, MaskAccessor> mask,
               pair<DestImageIterator, DestAccessor> dest,
               Functor f)
{
    combineTwoImagesIf(src1.first, src1.second, src1.third, 
                       src2.first, src2.second, 
                       mask.first, mask.second, 
                       dest.first, dest.second, f);
}

/********************************************************/
/*                                                      */
/*                  combineThreeImages                  */
/*                                                      */
/********************************************************/

/** Combine three source images into destination image.
    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class SrcImageIterator3, class SrcAccessor3,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineThreeImages(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   SrcImageIterator3 src2_upperleft, SrcAccessor3 sa3,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor f)
    }
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class SrcImageIterator3, class SrcAccessor3,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineThreeImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<SrcImageIterator3, SrcAccessor3> src3,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor f)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/combineimages.hxx]{../include/vigra/combineimages.hxx}\\
        Namespace: vigra
    
    \begin{verbatim}
    vigra::combineThreeImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                srcIter(src3.upperLeft()), 
                destIter(dest.upperLeft()),  
                SomeThreeArgumentFunctor());
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    SrcImageIterator3 src3_upperleft;
    DestImageIterator dest_upperleft;
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    SrcAccessor3 src3_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(src1_upperleft), 
              src2_accessor(src2_upperleft), 
              src3_accessor(src3_upperleft)), 
      dest_upperleft);

    \end{verbatim}
    
    @memo
    
*/
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
combineThreeImages(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   SrcImageIterator3 src3_upperleft, SrcAccessor3 sa3,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    int h = src1_lowerright.y - src1_upperleft.y;
    
    for(int y=0; y<h; 
            ++y, ++src1_upperleft.y, ++src2_upperleft.y, ++src3_upperleft.y, 
        ++dest_upperleft.y)
    {
        SrcImageIterator1 ix1(src1_upperleft);
        SrcImageIterator2 ix2(src2_upperleft);
        SrcImageIterator3 ix3(src3_upperleft);
        DestImageIterator dix(dest_upperleft);

        for(int x=0; x<w; ++x, ++ix1.x, ++ix2.x, ++ix3.x, ++dix.x)
        {
            da.set(static_cast<typename 
           DestAccessor::value_type>(f(sa1(ix1), sa2(ix2), sa3(ix3))), 
           dix);
        }
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
combineThreeImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
             pair<SrcImageIterator2, SrcAccessor2> src2,
             pair<SrcImageIterator3, SrcAccessor3> src3,
             pair<DestImageIterator, DestAccessor> dest,
             Functor f)
{
    combineThreeImages(src1.first, src1.second, src1.third, 
                     src2.first, src2.second, 
                     src3.first, src3.second, 
                     dest.first, dest.second, f);
}

    
//@}

/** @heading Functors to Combine Images

    @memo common functors with several arguments
*/
//@{

/********************************************************/
/*                                                      */
/*                    MagnitudeFunctor                  */
/*                                                      */
/********************************************************/

/** Calculate the magnitude from two arguments.
*/
template <class ValueType>
class MagnitudeFunctor
{
  public:
        /* the functor's value type
            @memo
        */
    typedef ValueType value_type;
    
        /** calculate transform '#sqrt(v1*v1 + v2*v2)#'. Note that
            the appropriate #sqrt()# function must be included
            
            @memo
        */
    value_type operator()(value_type const & v1, value_type const & v2) const
    {
        return sqrt(v1*v1 + v2*v2);
    }
};
//@}

} // namespace vigra

#endif // VIGRA_COMBINEIMAGES_HXX
