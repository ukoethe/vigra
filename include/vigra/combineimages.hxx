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
 
 
#ifndef VIGRA_COMBINEIMAGES_HXX
#define VIGRA_COMBINEIMAGES_HXX

#include "vigra/utilities.hxx"
#include "vigra/iteratortraits.hxx"
#include <cmath>

namespace vigra {

/** \addtogroup CombineAlgo Algorithms to Combine Images

    Apply functor to calculate a pixelwise transformation depending on multiple images.
    Note that the binary functors of the STL may be used with these functions.
*/
//@{

/********************************************************/
/*                                                      */
/*                    combineTwoImages                  */
/*                                                      */
/********************************************************/

/** \brief Combine two source images into destination image.

    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    The function uses accessors to access the pixel data.
    Note that the binary functors of the STL can be used in addition to
    the functors specifically defined in \ref CombineFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
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
                 Functor const & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineTwoImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="combineimages_8hxx-source.html">vigra/combineimages.hxx</a>"<br>
        Namespace: vigra
    
    \code
    #include <functional>     // for plus

    vigra::combineTwoImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                destIter(dest.upperLeft()),  
                std::plus<SrcValueType>());
    
    \endcode
    
    Note that <TT>SrcValueType</TT> must be replaced with the appropriate type (e.g. 
    the promote type of the input images' pixel type, see also 
    \URef{NumericandPromotionTraits})
    
    <b> Required Interface:</b>
    
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    DestImageIterator dest_upperleft;
    IteratorTraits<SrcImageIterator1>::RowIterator sx1(src1_upperleft);
    IteratorTraits<SrcImageIterator2>::RowIterator sx2(src2_upperleft);
    IteratorTraits<DestImageIterator>::RowIterator dx(dest_upperleft);
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(sx1), src2_accessor(sx2)), 
          dx);

    \endcode
    
    
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
         Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y; 
            ++src1_upperleft.y, ++src2_upperleft.y, ++dest_upperleft.y)
    {
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1(src1_upperleft);
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1end = ix1 + w;
        typename IteratorTraits<SrcImageIterator2>::RowIterator ix2(src2_upperleft);
        typename IteratorTraits<DestImageIterator>::RowIterator dix(dest_upperleft);

        for(; ix1 != ix1end; ++ix1, ++ix2, ++dix)
        {
            da.set(detail::RequiresExplicitCast<typename DestAccessor::value_type>::cast(f(sa1(ix1), sa2(ix2))), dix);
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
         Functor const & f)
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

/** \brief Combine ROI of two source images into destination image.

    The transformation given by the functor is applied to all source 
    pixels in the ROI (i.e. whenever the return value of the mask's accessor
    is not zero)
    and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    The function uses accessors to access the pixel data.
    Note that the binary functors of the STL can be used in addition to
    the functors specifically defined in \ref CombineFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
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
                   Functor const & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
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
                   Functor const & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="combineimages_8hxx-source.html">vigra/combineimages.hxx</a>"<br>
        Namespace: vigra
    
    \code
    #include <functional>     // for plus

    vigra::combineTwoImagesIf(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                maskIter(mask.upperLeft()), 
                destIter(dest.upperLeft()),
                std::plus<SrcValueType>());
    
    \endcode

    Note that <TT>SrcValueType</TT> must be replaced with the appropriate type (e.g. 
    the promote type of the input images' pixel type, see also 
    \URef{NumericandPromotionTraits})
    
    <b> Required Interface:</b>
    
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    MaskImageIterator mask_upperleft;
    DestImageIterator dest_upperleft;
    IteratorTraits<SrcImageIterator1>::RowIterator sx1(src1_upperleft);
    IteratorTraits<SrcImageIterator2>::RowIterator sx2(src2_upperleft);
    IteratorTraits<MaskImageIterator>::RowIterator mx(mask_upperleft);
    IteratorTraits<DestImageIterator>::RowIterator dx(dest_upperleft);
    
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    MaskAccessor mask_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;
    
    if(mask_accessor(mx))
       dest_accessor.set(
          functor(src1_accessor(sx1), src2_accessor(sx2)), 
          dx);

    \endcode
    
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
               Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y;
          ++src1_upperleft.y, ++src2_upperleft.y, 
          ++dest_upperleft.y, ++mask_upperleft.y)
    {
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1(src1_upperleft);
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1end = ix1 + w;
        typename IteratorTraits<SrcImageIterator2>::RowIterator ix2(src2_upperleft);
        typename IteratorTraits<MaskImageIterator>::RowIterator mix(mask_upperleft);
        typename IteratorTraits<DestImageIterator>::RowIterator dix(dest_upperleft);

        for(; ix1 != ix1end; ++ix1, ++ix2, ++dix, ++mix)
        {
            if(ma(mix)) da.set(detail::RequiresExplicitCast<typename DestAccessor::value_type>::cast(f(sa1(ix1), sa2(ix2))), dix);
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
               Functor const & f)
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

/** \brief Combine three source images into destination image.

    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
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
                   Functor const & f)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
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
                 Functor const & f)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="combineimages_8hxx-source.html">vigra/combineimages.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::combineThreeImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                srcIter(src3.upperLeft()), 
                destIter(dest.upperLeft()),  
                SomeThreeArgumentFunctor());
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    SrcImageIterator3 src3_upperleft;
    DestImageIterator dest_upperleft;
    IteratorTraits<SrcImageIterator1>::RowIterator sx1(src1_upperleft);
    IteratorTraits<SrcImageIterator2>::RowIterator sx2(src2_upperleft);
    IteratorTraits<SrcImageIterator3>::RowIterator sx3(src3_upperleft);
    IteratorTraits<DestImageIterator>::RowIterator dx(dest_upperleft);
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    SrcAccessor3 src3_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(sx1), 
                  src2_accessor(sx2), 
                  src3_accessor(sx3)), 
          dx);

    \endcode
    
    
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
                   Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y;
        ++src1_upperleft.y, ++src2_upperleft.y, ++src3_upperleft.y, 
        ++dest_upperleft.y)
    {
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1(src1_upperleft);
        typename IteratorTraits<SrcImageIterator1>::RowIterator ix1end = ix1 + w;
        typename IteratorTraits<SrcImageIterator2>::RowIterator ix2(src2_upperleft);
        typename IteratorTraits<SrcImageIterator3>::RowIterator ix3(src3_upperleft);
        typename IteratorTraits<DestImageIterator>::RowIterator dix(dest_upperleft);

        for(; ix1 != ix1end; ++ix1, ++ix2, ++ix3, ++dix)
        {
            da.set(detail::RequiresExplicitCast<typename DestAccessor::value_type>::cast(f(sa1(ix1), sa2(ix2), sa3(ix3))), 
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
             Functor const & f)
{
    combineThreeImages(src1.first, src1.second, src1.third, 
                     src2.first, src2.second, 
                     src3.first, src3.second, 
                     dest.first, dest.second, f);
}

    
//@}

/** \addtogroup CombineFunctor Functors to Combine Images

    Common functors with several arguments
*/
//@{

/********************************************************/
/*                                                      */
/*                    MagnitudeFunctor                  */
/*                                                      */
/********************************************************/

/** Calculate the magnitude from two arguments.
    Can be used in conjunction with \ref gradientBasedTransform().
*/
template <class ValueType>
class MagnitudeFunctor
{
  public:
        /* the functor's value type
        */
    typedef ValueType value_type;
    
        /** calculate transform '<TT>sqrt(v1*v1 + v2*v2)</TT>'. 
            
        */
    value_type operator()(value_type const & v1, value_type const & v2) const
    {
#ifndef CMATH_NOT_IN_STD
        return std::sqrt(v1*v1 + v2*v2);
#else
        return sqrt(v1*v1 + v2*v2);
#endif
    }
};

/********************************************************/
/*                                                      */
/*             RGBGradientMagnitudeFunctor              */
/*                                                      */
/********************************************************/


/** Calculate the gradient magnitude from RGB arguments.
    Can be used in conjunction with \ref gradientBasedTransform().
*/
template <class ValueType>
class RGBGradientMagnitudeFunctor
{
  public:
        /* the functor's value type
        */
    typedef ValueType value_type;
    
        /** Calculate the gradient magnitude form given RGB components.
            The function returns 
            
            \f[ \sqrt{|\nabla red|^2 + |\nabla green|^2 + |\nabla blue|^2}
            \f]
            
            where \f$|\nabla red|^2\f$ is defined by <TT>gx.red()*gx.red() + gy.red()*gy.red()</TT>.
            
            <TT>ValueType</TT> (the RGB's component type) must support addition, multiplication, 
            abd <TT>sqrt()</TT>.
        */
    ValueType 
    operator()(RGBValue<ValueType> const & gx, RGBValue<ValueType> const & gy) const
    {
#ifndef CMATH_NOT_IN_STD
        return std::sqrt(gx.red()*gx.red() + gx.green()*gx.green() +
                    gx.blue()*gx.blue() + gy.red()*gy.red() + 
                    gy.green()*gy.green() + gy.blue()*gy.blue());
#else
        return sqrt(gx.red()*gx.red() + gx.green()*gx.green() +
                    gx.blue()*gx.blue() + gy.red()*gy.red() + 
                    gy.green()*gy.green() + gy.blue()*gy.blue());
#endif
    }
};

//@}

} // namespace vigra

#endif // VIGRA_COMBINEIMAGES_HXX
