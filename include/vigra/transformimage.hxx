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


#ifndef VIGRA_TRANSFORMIMAGE_HXX
#define VIGRA_TRANSFORMIMAGE_HXX

#include "utilities.hxx"
#include "numerictraits.hxx"
#include "iteratortraits.hxx"
#include "rgbvalue.hxx"
#include "functortraits.hxx"
#include "inspectimage.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup TransformAlgo Algorithms to Transform Images
    Apply functor to calculate a pixelwise transformation of one image

    @{
*/

/********************************************************/
/*                                                      */
/*                      transformLine                   */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class Functor>
void
transformLine(SrcIterator s,
              SrcIterator send, SrcAccessor src,
              DestIterator d, DestAccessor dest,
              Functor const & f)
{
    for(; s != send; ++s, ++d)
        dest.set(f(src(s)), d);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor,
          class Functor>
void
transformLineIf(SrcIterator s,
                SrcIterator send, SrcAccessor src,
                MaskIterator m, MaskAccessor mask,
                DestIterator d, DestAccessor dest,
                Functor const & f)
{
    for(; s != send; ++s, ++d, ++m)
        if(mask(m))
            dest.set(f(src(s)), d);
}

/********************************************************/
/*                                                      */
/*                      transformImage                  */
/*                                                      */
/********************************************************/

/** \brief Apply unary point transformation to each pixel.

    After the introduction of arithmetic and algebraic \ref MultiMathModule "array expressions",
    this function is rarely needed. Moreover, \ref transformMultiArray() provides the 
    same functionality for arbitrary dimensional arrays.

    The transformation given by the functor is applied to every source
    pixel and the result written into the corresponding destination pixel.
    Note that the unary functors of the STL can be used in addition to
    the functors specifically defined in \ref TransformFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
              class T2, class S2, class Functor>
        void
        transformImage(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                       Functor const & f);
    }
    \endcode

    \deprecatedAPI{transformImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        transformImage(SrcImageIterator src_upperleft,
               SrcImageIterator src_lowerright, SrcAccessor sa,
               DestImageIterator dest_upperleft, DestAccessor da,
               Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
               pair<DestImageIterator, DestAccessor> dest,
               Functor const & f)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    #include <cmath>         // for sqrt()
    MultiArray<2, float>  src(100, 200),
                          dest(100, 200);
    ...
    
    transformImage(src, dest, &std::sqrt );
    \endcode

    \deprecatedUsage{transformImage}
    \code
    #include <cmath>         // for sqrt()
    FImage  src(100, 200),
            dest(100, 200);

    vigra::transformImage(srcImageRange(src),
                          destImage(dest),
                          (double(*)(double))&std::sqrt );

    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator      dest_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    Functor functor;

    dest_accessor.set(functor(src_accessor(sx)), dx);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void transformImage)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
void
transformImage(SrcImageIterator src_upperleft,
               SrcImageIterator src_lowerright, SrcAccessor sa,
               DestImageIterator dest_upperleft, DestAccessor da,
               Functor const & f)
{
    int w = src_lowerright.x - src_upperleft.x;

    for(; src_upperleft.y < src_lowerright.y; ++src_upperleft.y, ++dest_upperleft.y)
    {
        transformLine(src_upperleft.rowIterator(),
                      src_upperleft.rowIterator() + w, sa,
                      dest_upperleft.rowIterator(), da, f);
    }
}

template <class SrcImageIterator, class SrcAccessor,
      class DestImageIterator, class DestAccessor, class Functor>
inline void
transformImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
               pair<DestImageIterator, DestAccessor> dest,
               Functor const & f)
{
    transformImage(src.first, src.second, src.third,
                   dest.first, dest.second, f);
}

template <class T1, class S1,
      class T2, class S2, class Functor>
inline void
transformImage(MultiArrayView<2, T1, S1> const & src,
               MultiArrayView<2, T2, S2> dest,
               Functor const & f)
{
    vigra_precondition(src.shape() == dest.shape(),
        "transformImage(): shape mismatch between input and output.");
    transformImage(srcImageRange(src),
                   destImage(dest), f);
}

/********************************************************/
/*                                                      */
/*                   transformImageIf                   */
/*                                                      */
/********************************************************/

/** \brief Apply unary point transformation to each pixel within the ROI
    (i.e., where the mask is non-zero).

    After the introduction of arithmetic and algebraic \ref MultiMathModule "array expressions",
    this function is rarely needed. Moreover, \ref combineTwoMultiArrays() provides the 
    same functionality for arbitrary dimensional arrays.

    The transformation given by the functor is applied to every source
    pixel in the ROI (i.e. when the return value of the mask's accessor
    is not zero)
    and the result is written into the corresponding destination pixel.
    The function uses accessors to access the pixel data.
    Note that the unary functors of the STL can be used in addition to
    the functors specifically defined in \ref TransformFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class TM, class SM,
                  class T2, class S2,
                  class Functor>
        void
        transformImageIf(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, TM, SM> const & mask,
                         MultiArrayView<2, T2, S2> dest,
                         Functor const & f);
    }
    \endcode

    \deprecatedAPI{transformImageIf}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, clas DestAccessor,
                  class Functor>
        void
        transformImageIf(SrcImageIterator src_upperleft,
                         SrcImageIterator src_lowerright, SrcAccessor sa,
                         MaskImageIterator mask_upperleft, MaskAccessor ma,
                         DestImageIterator dest_upperleft, DestAccessor da,
                         Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, clas DestAccessor,
                  class Functor>
        void
        transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                         pair<MaskImageIterator, MaskAccessor> mask,
                         pair<DestImageIterator, DestAccessor> dest,
                         Functor const & f)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    #include <cmath>         // for sqrt()
    
    MultiArray<2, unsigned char>  mask(100, 200),
    MultiArray<2, float>          src(100, 200),
                                  dest(100, 200);
    ... // fill src and mask
    
    transformImageIf(src, mask, dest, &std::sqrt );
    \endcode

    \deprecatedUsage{transformImageIf}
    \code
    #include <cmath>         // for sqrt()

    vigra::transformImageIf(srcImageRange(src),
                            maskImage(mask),
                            destImage(dest),
                            (double(*)(double))&std::sqrt );

    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator  dest_upperleft;
    MaskImageIterator mask_upperleft;
    SrcImageIterator::row_iterator sx = src_upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    MaskAccessor mask_accessor;
    Functor functor;

    if(mask_accessor(mx))
       dest_accessor.set(functor(src_accessor(sx)), dx);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void transformImageIf)

template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
transformImageIf(SrcImageIterator src_upperleft,
                 SrcImageIterator src_lowerright, SrcAccessor sa,
                 MaskImageIterator mask_upperleft, MaskAccessor ma,
                 DestImageIterator dest_upperleft, DestAccessor da,
                 Functor const & f)
{
    int w = src_lowerright.x - src_upperleft.x;

    for(; src_upperleft.y < src_lowerright.y;
             ++src_upperleft.y, ++mask_upperleft.y, ++dest_upperleft.y)
    {
        transformLineIf(src_upperleft.rowIterator(),
                        src_upperleft.rowIterator() + w, sa,
                        mask_upperleft.rowIterator(), ma,
                        dest_upperleft.rowIterator(), da, f);
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline void
transformImageIf(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                 pair<MaskImageIterator, MaskAccessor> mask,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
{
    transformImageIf(src.first, src.second, src.third,
                     mask.first, mask.second,
                     dest.first, dest.second, f);
}

template <class T1, class S1,
          class TM, class SM,
          class T2, class S2,
          class Functor>
inline void
transformImageIf(MultiArrayView<2, T1, S1> const & src,
                 MultiArrayView<2, TM, SM> const & mask,
                 MultiArrayView<2, T2, S2> dest,
                 Functor const & f)
{
    vigra_precondition(src.shape() == mask.shape() && src.shape() == dest.shape(),
        "transformImageIf(): shape mismatch between input and output.");
    transformImageIf(srcImageRange(src),
                     maskImage(mask),
                     destImage(dest), f);
}

/********************************************************/
/*                                                      */
/*               gradientBasedTransform                 */
/*                                                      */
/********************************************************/

/** \brief Calculate a function of the image gradient.

    The gradient and the function represented by <TT>Functor f</TT>
    are calculated in one go: for each location, the symmetric
    difference in x- and y-directions (asymmetric difference at the
    image borders) are passed to the given functor, and the result is
    written to the destination image. Functors to be used with this
    function include \ref MagnitudeFunctor and \ref
    RGBGradientMagnitudeFunctor.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2, 
                  class Functor>
        void
        gradientBasedTransform(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest, 
                               Functor const & grad);
    }
    \endcode

    \deprecatedAPI{gradientBasedTransform}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        gradientBasedTransform(SrcImageIterator srcul, SrcImageIterator srclr, SrcAccessor sa,
                               DestImageIterator destul, DestAccessor da, Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor, class Functor>
        void
        gradientBasedTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                               pair<DestImageIterator, DestAccessor> dest, Functor const & const & f)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), magnitude(w,h);
    ...

    gradientBasedTransform(src, magnitude,
                           MagnitudeFunctor<float>());
    \endcode

    \deprecatedUsage{gradientBasedTransform}
    \code
    vigra::FImage src(w,h), magnitude(w,h);
    ...

    gradientBasedTransform(srcImageRange(src), destImage(magnitude),
                                vigra::MagnitudeFunctor<float>());
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator is, isend;
    DestImageIterator id;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    typename NumericTraits<typename SrcAccessor::value_type>::RealPromote
        diffx, diffy;

    diffx = src_accessor(is, Diff2D(-1,0)) - src_accessor(is, Diff2D(1,0));
    diffy = src_accessor(is, Diff2D(0,-1)) - src_accessor(is, Diff2D(0,1));

    Functor f;

    dest_accessor.set(f(diffx, diffy), id);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void gradientBasedTransform)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
void
gradientBasedTransform(SrcImageIterator srcul, SrcImageIterator srclr, SrcAccessor sa,
                       DestImageIterator destul, DestAccessor da, Functor const & grad)
{
    int w = srclr.x - srcul.x;
    int h = srclr.y - srcul.y;
    int x,y;

    SrcImageIterator sy = srcul;
    DestImageIterator dy = destul;

    const Diff2D left(-1,0);
    const Diff2D right(1,0);
    const Diff2D top(0,-1);
    const Diff2D bottom(0,1);

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    TmpType diffx, diffy;

    SrcImageIterator sx = sy;
    DestImageIterator dx = dy;

    diffx = sa(sx) - sa(sx, right);
    diffy = sa(sx) - sa(sx, bottom);
    da.set(grad(diffx, diffy), dx);

    for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
    {
        diffx = (sa(sx, left) - sa(sx, right)) / TmpType(2.0);
        diffy = sa(sx) - sa(sx, bottom);
        da.set(grad(diffx, diffy), dx);
    }

    diffx = sa(sx, left) - sa(sx);
    diffy = sa(sx) - sa(sx, bottom);
    da.set(grad(diffx, diffy), dx);

    ++sy.y;
    ++dy.y;

    for(y=2; y<h; ++y, ++sy.y, ++dy.y)
    {
        sx = sy;
        dx = dy;

        diffx = sa(sx) - sa(sx, right);
        diffy = (sa(sx, top) - sa(sx, bottom)) / TmpType(2.0);
        da.set(grad(diffx, diffy), dx);

        for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
        {
            diffx = (sa(sx, left) - sa(sx, right)) / TmpType(2.0);
            diffy = (sa(sx, top) - sa(sx, bottom)) / TmpType(2.0);
            da.set(grad(diffx, diffy), dx);
        }

        diffx = sa(sx, left) - sa(sx);
        diffy = (sa(sx, top) - sa(sx, bottom)) / TmpType(2.0);
        da.set(grad(diffx, diffy), dx);
    }

    sx = sy;
    dx = dy;

    diffx = sa(sx) - sa(sx, right);
    diffy = sa(sx, top) - sa(sx);
    da.set(grad(diffx, diffy), dx);

    for(x=2, ++sx.x, ++dx.x; x<w; ++x, ++sx.x, ++dx.x)
    {
        diffx = (sa(sx, left) - sa(sx, right)) / TmpType(2.0);
        diffy = sa(sx, top) - sa(sx);
        da.set(grad(diffx, diffy), dx);
    }

    diffx = sa(sx, left) - sa(sx);
    diffy = sa(sx, top) - sa(sx);
    da.set(grad(diffx, diffy), dx);
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor, class Functor>
inline void
gradientBasedTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                       pair<DestImageIterator, DestAccessor> dest, Functor const & grad)
{
    gradientBasedTransform(src.first, src.second, src.third,
                           dest.first, dest.second, grad);
}

template <class T1, class S1,
          class T2, class S2, class Functor>
inline void
gradientBasedTransform(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest, Functor const & grad)
{
    vigra_precondition(src.shape() == dest.shape(),
        "gradientBasedTransform(): shape mismatch between input and output.");
    gradientBasedTransform(srcImageRange(src),
                           destImage(dest), grad);
}

/** @} */
/** \addtogroup TransformFunctor Functors to Transform Images

    Note that the unary functors of the STL can also be used in
    connection with \ref transformImage().
*/
//@{

template <class DestValueType, class Multiplier = double>
class LinearIntensityTransform
{
  public:
        /* the functors argument type (actually, since
           <tt>operator()</tt> is a template, much more types are possible)
        */
    typedef DestValueType argument_type;

        /* the functors result type
        */
    typedef DestValueType result_type;

        /* \deprecated use argument_type and result_type
        */
    typedef DestValueType value_type;

        /* type of the offset (used in internal calculations to prevent
            overflows and minimize round-off errors).
        */
    typedef typename
            NumericTraits<DestValueType>::RealPromote argument_promote;

        /* type of the scale factor
        */
    typedef Multiplier scalar_multiplier_type;

        /* init scale and offset
        */
    LinearIntensityTransform(scalar_multiplier_type scale, argument_promote offset)
    : scale_(scale), offset_(offset)
    {}

        /* calculate transform
        */
    template <class SrcValueType>
    result_type operator()(SrcValueType const & s) const
    {
        return NumericTraits<result_type>::fromRealPromote(scale_ * (s + offset_));
    }

  private:

    scalar_multiplier_type scale_;
    argument_promote offset_;
};

template <class DestValueType, class Multiplier>
class FunctorTraits<LinearIntensityTransform<DestValueType, Multiplier> >
: public FunctorTraitsBase<LinearIntensityTransform<DestValueType, Multiplier> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

template <class DestValueType, class Multiplier = double>
class ScalarIntensityTransform
{
  public:
        /* the functors argument type (actually, since
           <tt>operator()</tt> is a template, much more types are possible)
        */
    typedef DestValueType argument_type;

        /* the functors result type
        */
    typedef DestValueType result_type;

        /* \deprecated use argument_type and result_type
        */
    typedef DestValueType value_type;

        /* type of the scale factor
        */
    typedef Multiplier scalar_multiplier_type;

        /* init scale
        */
    ScalarIntensityTransform(scalar_multiplier_type scale)
    : scale_(scale)
    {}

        /* calculate transform
        */
    template <class SrcValueType>
    result_type operator()(SrcValueType const & s) const
    {
        return NumericTraits<result_type>::fromRealPromote(scale_ * s);
    }

  private:
    scalar_multiplier_type scale_;
};

template <class DestValueType, class Multiplier>
class FunctorTraits<ScalarIntensityTransform<DestValueType, Multiplier> >
: public FunctorTraitsBase<ScalarIntensityTransform<DestValueType, Multiplier> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/********************************************************/
/*                                                      */
/*              linearIntensityTransform                */
/*                                                      */
/********************************************************/

/** \brief Apply a linear transform to the source pixel values

    Factory function for a functor that linearly transforms the
    source pixel values. The functor applies the transform
    '<TT>destvalue = scale * (srcvalue + offset)</TT>' to every pixel.
    This can, for example, be used to transform images into the visible
    range 0...255 or to invert an image.

    If you leave out the second parameter / offset, you will get an
    optimized version of the functor which only scales by the given
    factor, however you have to make the template parameter (pixel
    type) explicit then.

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <class Multiplier, class DestValueType>
        LinearIntensityTransform<DestValueType, Multiplier>
        linearIntensityTransform(Multiplier scale, DestValueType offset);

        template <class DestValueType, class Multiplier>
        ScalarIntensityTransform<DestValueType, Multiplier>
        linearIntensityTransform(Multiplier scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    vigra::IImage src(width, height);
    vigra::BImage dest(width, height);
    ...
    vigra::FindMinMax<IImage::PixelType> minmax;   // functor to find range

    vigra::inspectImage(srcImageRange(src), minmax); // find original range

    // transform to range 0...255
    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearIntensityTransform(
                            255.0 / (minmax.max - minmax.min), // scaling
                          - minmax.min));                    // offset
    \endcode

    The one-parameter version can be used like this:

    \code
    // scale from 0..255 to 0..1.0
    FImage dest(src.size());

    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearIntensityTransform<float>(1.0 / 255));
    \endcode

    <b> Required Interface:</b>

    The source and destination value types must be models of \ref LinearSpace in both cases.

*/
template <class Multiplier, class DestValueType>
LinearIntensityTransform<DestValueType, Multiplier>
linearIntensityTransform(Multiplier scale, DestValueType offset)
{
    return LinearIntensityTransform<DestValueType, Multiplier>(scale, offset);
}

template <class DestValueType, class Multiplier>
ScalarIntensityTransform<DestValueType, Multiplier>
linearIntensityTransform(Multiplier scale)
{
    return ScalarIntensityTransform<DestValueType, Multiplier>(scale);
}

/********************************************************/
/*                                                      */
/*                   linearRangeMapping                 */
/*                                                      */
/********************************************************/

/** \brief Map a source intensity range linearly to a destination range.

    Factory function for a functor that linearly transforms the
    source pixel values. The functor applies the transform
    '<TT>destvalue = scale * (srcvalue + offset)</TT>' to every pixel,
    where <tt>scale = (dest_max - dest_min) / (src_max - src_min)</tt>
    and <tt>offset = dest_min / scale - src_min</tt>. As a result,
    the pixel values <tt>src_max</tt>, <tt>src_min</tt> in the source image
    are mapped onto <tt>dest_max</tt>, <tt>dest_min</tt> respectively.
    This works for scalar as well as vector pixel types. Instead of 
    <tt>src_min</tt> and <tt>src_max</tt>, you may also pass a functor
    \ref FindMinMax. 

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <class SrcValueType, class DestValueType>
        LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
        linearRangeMapping(SrcValueType src_min, SrcValueType src_max,
                           DestValueType dest_min, DestValueType dest_max );

                           template <class SrcValueType, class DestValueType>
        LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
        linearRangeMapping(SrcValueType src_min, SrcValueType src_max,
                           DestValueType dest_min, DestValueType dest_max );
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    vigra::IImage src(width, height);
    vigra::BImage dest(width, height);
    ...
    vigra::FindMinMax<IImage::PixelType> minmax;   // functor to find range

    vigra::inspectImage(srcImageRange(src), minmax); // find original range

    // transform to range 0...255
    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearRangeMapping(
                            minmax.min, minmax.max,  // src range
                            0, 255)                  // dest range
                          );

    // equivalent, but shorter
    vigra::transformImage(srcImageRange(src), destImage(dest),
                          linearRangeMapping(
                            minmax,                 // src range
                            0, 255)                 // dest range
                          );
    \endcode

    <b> Required Interface:</b>

    The source and destination value types must be models of \ref LinearSpace in both cases.

*/
template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(SrcValueType src_min, SrcValueType src_max,
                   DestValueType dest_min, DestValueType dest_max )
{
    return linearRangeMapping(src_min, src_max, dest_min, dest_max,
            typename NumericTraits<DestValueType>::isScalar());
}

template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(FindMinMax<SrcValueType> const & src,
                   DestValueType dest_min, DestValueType dest_max )
{
    return linearRangeMapping(src.min, src.max, dest_min, dest_max,
            typename NumericTraits<DestValueType>::isScalar());
}

template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(
    SrcValueType src_min, SrcValueType src_max,
    DestValueType dest_min, DestValueType dest_max,
    VigraTrueType /* isScalar */ )
{
    typedef typename NumericTraits<DestValueType>::RealPromote Multiplier;
    Multiplier diff = src_max - src_min;
    Multiplier scale = diff == NumericTraits<Multiplier>::zero()
                     ? NumericTraits<Multiplier>::one()
                     : (dest_max - dest_min) / diff;
    return LinearIntensityTransform<DestValueType, Multiplier>(
                                   scale, dest_min / scale - src_min );
}

template <class SrcValueType, class DestValueType>
LinearIntensityTransform<DestValueType, typename NumericTraits<DestValueType>::RealPromote>
linearRangeMapping(
    SrcValueType src_min, SrcValueType src_max,
    DestValueType dest_min, DestValueType dest_max,
    VigraFalseType /* isScalar */ )
{
    typedef typename NumericTraits<DestValueType>::RealPromote Multiplier;
    typedef typename Multiplier::value_type MComponent;
    Multiplier scale(dest_max), offset(dest_max);
    for(unsigned int i=0; i<src_min.size(); ++i)
    {
        MComponent diff = src_max[i] - src_min[i];
        scale[i] = diff == NumericTraits<MComponent>::zero()
                     ? NumericTraits<MComponent>::one()
                     : (dest_max[i] - dest_min[i]) / diff;
        offset[i] = dest_min[i] / scale[i] - src_min[i];
    }
    return LinearIntensityTransform<DestValueType, Multiplier>(scale, offset);
}

/********************************************************/
/*                                                      */
/*                      Threshold                       */
/*                                                      */
/********************************************************/

/** \brief Threshold an image.

    If a source pixel is above or equal the lower and below
    or equal the higher threshold (i.e. within the closed interval
    [lower, higher]) the destination pixel is set to 'yesresult',
    otherwise to 'noresult'.

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(width, height), dest(width, height);
    ...
    vigra::transformImage(src.upperLeft(), src.lowerRight(), src.accessor(),
       dest.upperLeft(), dest.accessor(),
       vigra::Threshold<
          vigra::BImage::PixelType, vigra::BImage::PixelType>(10, 100, 0, 255));

    \endcode

    <b> Required Interface:</b>

    \code

    SrcValueType   src;
    DestValueType  dest, yesresult, noresult;

    dest = ((src < lower) || (higher < src)) ? noresult : yesresult;

    \endcode

*/
template <class SrcValueType, class DestValueType>
class Threshold
{
   public:

        /** the functor's argument type
        */
    typedef SrcValueType argument_type;

        /** the functor's result type
        */
    typedef DestValueType result_type;

        /** init thresholds and return values
        */
    Threshold(argument_type lower, argument_type higher,
              result_type noresult, result_type yesresult)
    : lower_(lower), higher_(higher),
      yesresult_(yesresult), noresult_(noresult)
    {}

        /** calculate transform
        */
    result_type operator()(argument_type s) const
    {
        return ((s < lower_) || (higher_ < s)) ? noresult_ : yesresult_;
    }

  private:

    argument_type lower_, higher_;
    result_type yesresult_, noresult_;
};

template <class SrcValueType, class DestValueType>
class FunctorTraits<Threshold<SrcValueType, DestValueType> >
: public FunctorTraitsBase<Threshold<SrcValueType, DestValueType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/********************************************************/
/*                                                      */
/*                BrightnessContrastFunctor             */
/*                                                      */
/********************************************************/

/** \brief Adjust brightness and contrast of an image.

    This functor applies a gamma correction to each pixel in order to
    modify the brightness of the image. To the result of the gamma
    correction, another transform is applied that modifies the
    contrast. The brightness and contrast parameters must be
    positive. Values greater than 1 will increase image brightness or
    contrast respectively, values smaller than 1 decrease them.  A
    value of exactly 1 will have no effect.  If contrast is set to 1,
    the result is equivalent to that of the GammaFunctor with gamma =
    1./brightness.

    For \ref RGBValue "RGBValue's", the transforms are applied
    component-wise. The pixel values are assumed to lie between the
    given minimum and maximum values (in case of RGB, this is again
    understood component-wise). In case of <TT>unsigned char</TT>, min
    and max default to 0 and 255 respectively.  Precisely, the
    following transform is applied to each <em> PixelValue</em>:

    \f[
    \begin{array}{rcl}
    V_1 & = & \frac{PixelValue - min}{max - min} \\
    V_2 & = & V_1^\frac{1}{brightness} \\
    V_3 & = & 2 V_2 - 1 \\
    V_4 & = & \left\lbrace
        \begin{array}{l}
         V_3^\frac{1}{contrast} \mbox{\rm \quad if  } V_3 \ge 0 \\
         - (-V_3)^\frac{1}{contrast} \mbox{\rm \quad otherwise}
        \end{array} \right. \\
    Result & = & \frac{V_4 + 1}{2} (max - min) + min
    \end{array}
    \f]

    If the <TT>PixelType</TT> is <TT>unsigned char</TT>, a look-up-table is used
    for faster computation.

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage bimage(width, height);
    double brightness, contrast;
    ...
    vigra::transformImage(srcImageRange(bimage), destImage(bimage),
       vigra::BrightnessContrastFunctor<unsigned char>(brightness, contrast));



    vigra::FImage fimage(width, height);
    ...

    vigra::FindMinMax<float> minmax;
    vigra::inspectImage(srcImageRange(fimage), minmax);

    vigra::transformImage(srcImageRange(fimage), destImage(fimage),
       vigra::BrightnessContrastFunctor<float>(brightness, contrast, minmax.min, minmax.max));


    \endcode

    <b> Required Interface:</b>

    Scalar types: must be a linear algebra (+, - *, NumericTraits),
    strict weakly ordered (<), and <TT>pow()</TT> must be defined.

    RGB values: the component type must meet the above requirements.
*/
template <class PixelType>
class BrightnessContrastFunctor
{
    typedef typename
        NumericTraits<PixelType>::RealPromote promote_type;

 public:

        /** the functor's argument type
        */
    typedef PixelType argument_type;

        /** the functor's result type
        */
    typedef PixelType result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef PixelType value_type;

        /** Init functor for argument range <TT>[min, max]</TT>.
            <TT>brightness</TT> and <TT>contrast</TT> values > 1 will
            increase brightness and contrast, < 1 will decrease them, and == 1 means
            no change.
        */
    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              argument_type const & min, argument_type const & max)
    : b_(1.0/brightness),
      c_(1.0/contrast),
      min_(min),
      diff_(max - min),
      zero_(NumericTraits<promote_type>::zero()),
      one_(NumericTraits<promote_type>::one())
    {}

        /** Calculate modified gray or color value
        */
    result_type operator()(argument_type const & v) const
    {
        promote_type v1 = (v - min_) / diff_;
        promote_type brighter = VIGRA_CSTD::pow(v1, b_);
        promote_type v2 = 2.0 * brighter - one_;
        promote_type contrasted = (v2 < zero_) ?
                                     -VIGRA_CSTD::pow(-v2, c_) :
                                      VIGRA_CSTD::pow(v2, c_);
        return result_type(0.5 * diff_ * (contrasted + one_) + min_);
    }

  private:
    promote_type b_, c_;
    argument_type min_;
    promote_type diff_, zero_, one_;
};

template <>
class BrightnessContrastFunctor<unsigned char>
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
     unsigned char lut[256];

 public:

    typedef unsigned char value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min = 0, value_type const & max = 255)
    {
        BrightnessContrastFunctor<promote_type> f(brightness, contrast, min, max);

        for(int i = min; i <= max; ++i)
        {
            lut[i] = static_cast<unsigned char>(f(i)+0.5);
        }
    }

    value_type operator()(value_type const & v) const
    {

        return lut[v];
    }
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class ComponentType>
class BrightnessContrastFunctor<RGBValue<ComponentType> >
{
    typedef typename
        NumericTraits<ComponentType>::RealPromote promote_type;
    BrightnessContrastFunctor<ComponentType> red, green, blue;

 public:

    typedef RGBValue<ComponentType> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class BrightnessContrastFunctor<RGBValue<int> >
{
    typedef NumericTraits<int>::RealPromote promote_type;
    BrightnessContrastFunctor<int> red, green, blue;

 public:

    typedef RGBValue<int> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

template <>
class BrightnessContrastFunctor<RGBValue<float> >
{
    typedef NumericTraits<float>::RealPromote promote_type;
    BrightnessContrastFunctor<float> red, green, blue;

 public:

    typedef RGBValue<float> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
                              value_type const & min, value_type const & max)
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

template <class PixelType>
class FunctorTraits<BrightnessContrastFunctor<PixelType> >
: public FunctorTraitsBase<BrightnessContrastFunctor<PixelType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class BrightnessContrastFunctor<RGBValue<unsigned char> >
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
    BrightnessContrastFunctor<unsigned char> red, green, blue;

 public:

    typedef RGBValue<unsigned char> value_type;

    BrightnessContrastFunctor(promote_type brightness, promote_type contrast,
       value_type const & min = value_type(0,0,0),
       value_type const & max = value_type(255, 255, 255))
    : red(brightness, contrast, min.red(), max.red()),
      green(brightness, contrast, min.green(), max.green()),
      blue(brightness, contrast, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {

        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};



/********************************************************/
/*                                                      */
/*                     GammaFunctor                     */
/*                                                      */
/********************************************************/

/** \brief Perform gamma correction of an image.

    This functor applies a gamma correction to each pixel in order to
    modify the brightness of the image.  Gamma values smaller than 1
    will increase image brightness, whereas values greater than 1
    decrease it. A value of gamma = 1 will have no effect.  (See also
    BrightnessContrastFunctor, which additionally changes the
    contrast.)

    For \ref RGBValue "RGBValue's", the transforms are applied
    component-wise.  For ease of use, the pixel values are assumed to
    lie between the given minimum and maximum values (in case of RGB,
    this is again understood component-wise). In case of <TT>unsigned
    char</TT>, min and max default to 0 and 255 respectively.
    Precisely, the following transform is applied to each <em>
    PixelValue</em>:

    \f[
    \begin{array}{rcl}
    V_1 & = & \frac{PixelValue - min}{max - min} \\
    V_2 & = & V_1^{gamma} \\
    Result & = & V_2 (max - min) + min
    \end{array}
    \f]

    If the <TT>PixelType</TT> is <TT>unsigned char</TT>, a
    look-up-table is used for faster computation.

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage bimage(width, height);
    double gamma;
    ...
    vigra::transformImage(srcImageRange(bimage), destImage(bimage),
       vigra::GammaFunctor<unsigned char>(gamma));



    vigra::FImage fimage(width, height);
    ...

    vigra::FindMinMax<float> minmax;
    vigra::inspectImage(srcImageRange(fimage), minmax);

    vigra::transformImage(srcImageRange(fimage), destImage(fimage),
       vigra::GammaFunctor<float>(gamma, minmax.min, minmax.max));

    \endcode

    <b> Required Interface:</b>

    Scalar types: must be a linear algebra (+, - *, NumericTraits),
    strict weakly ordered (<), and <TT>pow()</TT> must be defined.

    RGB values: the component type must meet the above requirements.
*/
template <class PixelType>
class GammaFunctor
{
    typedef typename
        NumericTraits<PixelType>::RealPromote promote_type;

 public:

        /** the functor's argument type
        */
    typedef PixelType argument_type;

        /** the functor's result type
        */
    typedef PixelType result_type;

        /** \deprecated use argument_type and result_type
        */
    typedef PixelType value_type;

        /** Init functor for argument range <TT>[min, max]</TT>.
            <TT>gamma</TT> values < 1 will increase brightness, > 1
            will decrease it (gamma == 1 means no change).
        */
    GammaFunctor(double gamma,
                 argument_type const & min, argument_type const & max)
    : gamma_((promote_type)gamma),
      min_(min),
      diff_(max - min),
      zero_(NumericTraits<promote_type>::zero()),
      one_(NumericTraits<promote_type>::one())
    {}

        /** Calculate modified gray or color value
        */
    result_type operator()(argument_type const & v) const
    {
        promote_type v1 = (v - min_) / diff_;
        promote_type brighter = VIGRA_CSTD::pow(v1, gamma_);
        return result_type(diff_ * brighter + min_);
    }

  private:
    promote_type gamma_;
    argument_type min_;
    promote_type diff_, zero_, one_;
};

template <>
class GammaFunctor<unsigned char>
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
     unsigned char lut[256];

 public:

    typedef unsigned char value_type;

    GammaFunctor(promote_type gamma,
                 value_type const & min = 0, value_type const & max = 255)
    {
        GammaFunctor<promote_type> f(gamma, min, max);

        for(int i = min; i <= max; ++i)
        {
            lut[i] = static_cast<unsigned char>(f(i)+0.5);
        }
    }

    value_type operator()(value_type const & v) const
    {

        return lut[v];
    }
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class ComponentType>
class GammaFunctor<RGBValue<ComponentType> >
{
    typedef typename
        NumericTraits<ComponentType>::RealPromote promote_type;
    GammaFunctor<ComponentType> red, green, blue;

 public:

    typedef RGBValue<ComponentType> value_type;

    GammaFunctor(promote_type gamma,
                 value_type const & min, value_type const & max)
    : red(gamma, min.red(), max.red()),
      green(gamma, min.green(), max.green()),
      blue(gamma, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {
        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class GammaFunctor<RGBValue<int> >
{
    typedef NumericTraits<int>::RealPromote promote_type;
    GammaFunctor<int> red, green, blue;

 public:

    typedef RGBValue<int> value_type;

    GammaFunctor(promote_type gamma,
                 value_type const & min, value_type const & max)
    : red(gamma, min.red(), max.red()),
      green(gamma, min.green(), max.green()),
      blue(gamma, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {
        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

template <>
class GammaFunctor<RGBValue<float> >
{
    typedef NumericTraits<float>::RealPromote promote_type;
    GammaFunctor<float> red, green, blue;

 public:

    typedef RGBValue<float> value_type;

    GammaFunctor(promote_type gamma,
                 value_type const & min, value_type const & max)
    : red(gamma, min.red(), max.red()),
      green(gamma, min.green(), max.green()),
      blue(gamma, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {
        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};

template <class PixelType>
class FunctorTraits<GammaFunctor<PixelType> >
: public FunctorTraitsBase<GammaFunctor<PixelType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <>
class GammaFunctor<RGBValue<unsigned char> >
{
    typedef NumericTraits<unsigned char>::RealPromote promote_type;
    GammaFunctor<unsigned char> red, green, blue;

 public:
    typedef RGBValue<unsigned char> value_type;

    GammaFunctor(promote_type gamma,
                 value_type const & min = value_type(0,0,0),
                 value_type const & max = value_type(255, 255, 255))
    : red(gamma, min.red(), max.red()),
      green(gamma, min.green(), max.green()),
      blue(gamma, min.blue(), max.blue())
    {}

    value_type operator()(value_type const & v) const
    {
        return value_type(red(v.red()), green(v.green()), blue(v.blue()));
    }
};


/********************************************************/
/*                                                      */
/*                     VectorNormFunctor                */
/*                                                      */
/********************************************************/

/** \brief A functor for computing the vector norm

    Calculate the magnitude or norm from a given vector-valued
    entity. The vector type will typically be some sort of
    ref vigra::TinyVector. If the vector is represented by a pair of
    scalar-valued images, use \ref vigra::MagnitudeFunctor instead.

    At least, the vector type is required to have a function
    '<em>result</em><TT> = dot(v,v)</TT>'.

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/transformimage.hxx\><br>
    Namespace: vigra

    \code
    typedef vigra::TinyVector<float, 2> Vector;
    vigra::BasicImage<Vector> grad(width, height);
    vigra::FImage magn(width,height);
    ...
    vigra::transformImage(srcImageRange(grad), destImage(magn),
                          VectorNormFunctor<float>()
                          );
    \endcode

    \see vigra::TinyVector, dot(), vigra::MagnitudeFunctor
*/
template <class ValueType>
class VectorNormFunctor
{
public:
  /** the functor's argument type
   */
  typedef ValueType argument_type;

  /** the functor's result type
   */
  typedef typename NumericTraits<typename ValueType::value_type>::RealPromote result_type;

  /** calculate transform '<TT>sqrt(v1*v1 + v2*v2 + ...)</TT>'.
   */
  result_type operator()( const argument_type &a ) const
  {
    return VIGRA_CSTD::sqrt( dot(a,a) );
  }
};    //-- class VectorNormFunctor

template <class ValueType>
class FunctorTraits<VectorNormFunctor<ValueType> >
: public FunctorTraitsBase<VectorNormFunctor<ValueType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief A functor for computing the squared vector norm

    Calculate the squared magnitude or norm from a given
    vector-valued entity. The vector type will typically be some
    sort of TinyVector.

    At least, the vector type is required to have a function
    '<em>result</em><TT> = dot(v,v)</TT>'.

    For an example of its usage see VectorNormFunctor

    <b> Traits defined:</b>

    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType</tt>)

    \see TinyVector, dot()
*/
template <class ValueType>
class VectorNormSqFunctor
{
public:
  /** the functor's argument type
   */
  typedef ValueType argument_type;

  /** the functor's result type
   */
  typedef typename NumericTraits<typename ValueType::value_type>::RealPromote result_type;

  /** calculate transform '<TT>v1*v1 + v2*v2 + ...</TT>'.
   */
  result_type operator()( const argument_type &a ) const
  {
    return dot(a,a);
  }
};    //-- class VectorNormSqFunctor

template <class ValueType>
class FunctorTraits<VectorNormSqFunctor<ValueType> >
: public FunctorTraitsBase<VectorNormSqFunctor<ValueType> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

//@}

} // namespace vigra

#endif // VIGRA_TRANSFORMIMAGE_HXX
