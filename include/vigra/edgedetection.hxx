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


#ifndef VIGRA_EDGEDETECTION_HXX
#define VIGRA_EDGEDETECTION_HXX

#include <vector>
#include <queue>
#include <cmath>     // sqrt(), abs()
#include "utilities.hxx"
#include "numerictraits.hxx"
#include "stdimage.hxx"
#include "stdimagefunctions.hxx"
#include "recursiveconvolution.hxx"
#include "separableconvolution.hxx"
#include "convolution.hxx"
#include "labelimage.hxx"
#include "mathutil.hxx"
#include "pixelneighborhood.hxx"
#include "linear_solve.hxx"
#include "functorexpression.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup EdgeDetection Edge Detection
    Edge detectors based on first and second derivatives,
          and related post-processing.
*/
//@{

/********************************************************/
/*                                                      */
/*           differenceOfExponentialEdgeImage           */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in an edge image using the Shen/Castan zero-crossing detector.

    This operator applies an exponential filter to the source image
    at the given <TT>scale</TT> and subtracts the result from the original image.
    Zero crossings are detected in the resulting difference image. Whenever the
    gradient at a zero crossing is greater than the given <TT>gradient_threshold</TT>,
    an edge point is marked (using <TT>edge_marker</TT>) in the destination image on
    the darker side of the zero crossing (note that zero crossings occur
    <i>between</i> pixels). For example:

    \code
    sign of difference image     resulting edge points (*)

        + - -                          * * .
        + + -               =>         . * *
        + + +                          . . .
    \endcode

    Non-edge pixels (<TT>.</TT>) remain untouched in the destination image.
    The result can be improved by the post-processing operation \ref removeShortEdges().
    A more accurate edge placement can be achieved with the function
    \ref differenceOfExponentialCrackEdgeImage().

    The source value type
    (<TT>SrcAccessor::value_type</TT>) must be a linear algebra, i.e. addition,
    subtraction and multiplication of the type with itself, and multiplication
    with double and
    \ref NumericTraits "NumericTraits" must
    be defined. In addition, this type must be less-comparable.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class GradValue, class DestValue>
        void 
        differenceOfExponentialEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                         MultiArrayView<2, T2, S2> dest,
                                         double scale,
                                         GradValue gradient_threshold,
                                         DestValue edge_marker = NumericTraits<DestValue>::one());
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor,
              class GradValue,
              class DestValue = DestAccessor::value_type>
        void differenceOfExponentialEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               DestIterator dul, DestAccessor da,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor,
              class GradValue,
              class DestValue = DestAccessor::value_type>
        void differenceOfExponentialEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(w,h);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::differenceOfExponentialEdgeImage(srcImageRange(src), destImage(edges),
                                     0.8, 4.0, 1);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);
    double d;
    GradValue gradient_threshold;

    u = u + u
    u = u - u
    u = u * u
    u = d * u
    u < gradient_threshold

    DestValue edge_marker;
    dest_accessor.set(edge_marker, dest_upperleft);
    \endcode

    <b> Preconditions:</b>

    \code
    scale > 0
    gradient_threshold > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void differenceOfExponentialEdgeImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void differenceOfExponentialEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    vigra_precondition(scale > 0,
                 "differenceOfExponentialEdgeImage(): scale > 0 required.");

    vigra_precondition(gradient_threshold > 0,
                 "differenceOfExponentialEdgeImage(): "
         "gradient_threshold > 0 required.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int x,y;

    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
    TMPTYPE;
    typedef BasicImage<TMPTYPE> TMPIMG;

    TMPIMG tmp(w,h);
    TMPIMG smooth(w,h);

    recursiveSmoothX(srcIterRange(sul, slr, sa), destImage(tmp), scale / 2.0);
    recursiveSmoothY(srcImageRange(tmp), destImage(tmp), scale / 2.0);

    recursiveSmoothX(srcImageRange(tmp), destImage(smooth), scale);
    recursiveSmoothY(srcImageRange(smooth), destImage(smooth), scale);

    typename TMPIMG::Iterator iy = smooth.upperLeft();
    typename TMPIMG::Iterator ty = tmp.upperLeft();
    DestIterator              dy = dul;

    const Diff2D right(1, 0);
    const Diff2D bottom(0, 1);


    TMPTYPE thresh = detail::RequiresExplicitCast<TMPTYPE>::cast((gradient_threshold * gradient_threshold) *
                     NumericTraits<TMPTYPE>::one());
    TMPTYPE zero = NumericTraits<TMPTYPE>::zero();

    for(y=0; y<h-1; ++y, ++iy.y, ++ty.y, ++dy.y)
    {
        typename TMPIMG::Iterator ix = iy;
        typename TMPIMG::Iterator tx = ty;
        DestIterator              dx = dy;

        for(x=0; x<w-1; ++x, ++ix.x, ++tx.x, ++dx.x)
        {
            TMPTYPE diff = *tx - *ix;
            TMPTYPE gx = tx[right] - *tx;
            TMPTYPE gy = tx[bottom] - *tx;

            if((gx * gx > thresh) &&
                (diff * (tx[right] - ix[right]) < zero))
            {
                if(gx < zero)
                {
                    da.set(edge_marker, dx, right);
                }
                else
                {
                    da.set(edge_marker, dx);
                }
            }
            if(((gy * gy > thresh) &&
                (diff * (tx[bottom] - ix[bottom]) < zero)))
            {
                if(gy < zero)
                {
                    da.set(edge_marker, dx, bottom);
                }
                else
                {
                    da.set(edge_marker, dx);
                }
            }
        }
    }

    typename TMPIMG::Iterator ix = iy;
    typename TMPIMG::Iterator tx = ty;
    DestIterator              dx = dy;

    for(x=0; x<w-1; ++x, ++ix.x, ++tx.x, ++dx.x)
    {
        TMPTYPE diff = *tx - *ix;
        TMPTYPE gx = tx[right] - *tx;

        if((gx * gx > thresh) &&
           (diff * (tx[right] - ix[right]) < zero))
        {
            if(gx < zero)
            {
                da.set(edge_marker, dx, right);
            }
            else
            {
                da.set(edge_marker, dx);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue>
inline
void differenceOfExponentialEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold)
{
    differenceOfExponentialEdgeImage(sul, slr, sa, dul, da,
                                        scale, gradient_threshold, 1);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
inline void 
differenceOfExponentialEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 double scale, GradValue gradient_threshold,
                                 DestValue edge_marker)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                     dest.first, dest.second,
                                     scale, gradient_threshold, edge_marker);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue>
inline void
differenceOfExponentialEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 double scale, GradValue gradient_threshold)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                     dest.first, dest.second,
                                     scale, gradient_threshold, 1);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue, class DestValue>
inline void 
differenceOfExponentialEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                 MultiArrayView<2, T2, S2> dest,
                                 double scale,
                                 GradValue gradient_threshold,
                                 DestValue edge_marker)
{
    vigra_precondition(src.shape() == dest.shape(),
        "differenceOfExponentialEdgeImage(): shape mismatch between input and output.");
    differenceOfExponentialEdgeImage(srcImageRange(src),
                                     destImage(dest),
                                     scale, gradient_threshold, edge_marker);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue>
inline void
differenceOfExponentialEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                 MultiArrayView<2, T2, S2> dest,
                                 double scale, GradValue gradient_threshold)
{
    vigra_precondition(src.shape() == dest.shape(),
        "differenceOfExponentialEdgeImage(): shape mismatch between input and output.");
    differenceOfExponentialEdgeImage(srcImageRange(src),
                                     destImage(dest),
                                     scale, gradient_threshold, T2(1));
}

/********************************************************/
/*                                                      */
/*        differenceOfExponentialCrackEdgeImage         */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in a crack edge image using the Shen/Castan zero-crossing detector.

    This operator applies an exponential filter to the source image
    at the given <TT>scale</TT> and subtracts the result from the original image.
    Zero crossings are detected in the resulting difference image. Whenever the
    gradient at a zero crossing is greater than the given <TT>gradient_threshold</TT>,
    an edge point is marked (using <TT>edge_marker</TT>) in the destination image
    <i>between</i> the corresponding original pixels. Topologically, this means we
    must insert additional pixels between the original ones to represent the
    boundaries between the pixels (the so called zero- and one-cells, with the original
    pixels being two-cells). Within VIGRA, such an image is called \ref CrackEdgeImage.
    To allow insertion of the zero- and one-cells, the destination image must have twice the
    size of the original (precisely, <TT>(2*w-1)</TT> by <TT>(2*h-1)</TT> pixels). Then the algorithm
    proceeds as follows:

    \code
    sign of difference image     insert zero- and one-cells     resulting edge points (*)
    
                                         + . - . -                   . * . . .
          + - -                          . . . . .                   . * * * .
          + + -               =>         + . + . -           =>      . . . * .
          + + +                          . . . . .                   . . . * *
                                         + . + . +                   . . . . .
    \endcode

    Thus the edge points are marked where they actually are - in between the pixels.
    An important property of the resulting edge image is that it conforms to the notion
    of well-composedness as defined by Latecki et al., i.e. connected regions and edges
    obtained by a subsequent \ref Labeling do not depend on
    whether 4- or 8-connectivity is used.
    The non-edge pixels (<TT>.</TT>) in the destination image remain unchanged.
    The result conforms to the requirements of a \ref CrackEdgeImage. It can be further
    improved by the post-processing operations \ref removeShortEdges() and
    \ref closeGapsInCrackEdgeImage().

    The source value type (<TT>SrcAccessor::value_type</TT>) must be a linear algebra, i.e. addition,
    subtraction and multiplication of the type with itself, and multiplication
    with double and
    \ref NumericTraits "NumericTraits" must
    be defined. In addition, this type must be less-comparable.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class GradValue, class DestValue>
        void
        differenceOfExponentialCrackEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                              MultiArrayView<2, T2, S2> dest,
                                              double scale,
                                              GradValue gradient_threshold,
                                              DestValue edge_marker = NumericTraits<DestValue>::one());
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor,
              class GradValue,
              class DestValue = DestAccessor::value_type>
        void differenceOfExponentialCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               DestIterator dul, DestAccessor da,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor,
              class GradValue,
              class DestValue = DestAccessor::value_type>
        void differenceOfExponentialCrackEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(2*w-1,2*h-1);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::differenceOfExponentialCrackEdgeImage(srcImageRange(src), destImage(edges),
                                     0.8, 4.0, 1);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);
    double d;
    GradValue gradient_threshold;

    u = u + u
    u = u - u
    u = u * u
    u = d * u
    u < gradient_threshold

    DestValue edge_marker;
    dest_accessor.set(edge_marker, dest_upperleft);
    \endcode

    <b> Preconditions:</b>

    \code
    scale > 0
    gradient_threshold > 0
    \endcode

    The destination image must have twice the size of the source:
    \code
    w_dest = 2 * w_src - 1
    h_dest = 2 * h_src - 1
    \endcode
*/
doxygen_overloaded_function(template <...> void differenceOfExponentialCrackEdgeImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void differenceOfExponentialCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold,
           DestValue edge_marker)
{
    vigra_precondition(scale > 0,
                 "differenceOfExponentialCrackEdgeImage(): scale > 0 required.");

    vigra_precondition(gradient_threshold > 0,
                 "differenceOfExponentialCrackEdgeImage(): "
         "gradient_threshold > 0 required.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int x, y;

    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
    TMPTYPE;
    typedef BasicImage<TMPTYPE> TMPIMG;

    TMPIMG tmp(w,h);
    TMPIMG smooth(w,h);

    TMPTYPE zero = NumericTraits<TMPTYPE>::zero();

    const Diff2D right(1,0);
    const Diff2D bottom(0,1);
    const Diff2D left(-1,0);
    const Diff2D top(0,-1);

    recursiveSmoothX(srcIterRange(sul, slr, sa), destImage(tmp), scale / 2.0);
    recursiveSmoothY(srcImageRange(tmp), destImage(tmp), scale / 2.0);

    recursiveSmoothX(srcImageRange(tmp), destImage(smooth), scale);
    recursiveSmoothY(srcImageRange(smooth), destImage(smooth), scale);

    typename TMPIMG::Iterator iy = smooth.upperLeft();
    typename TMPIMG::Iterator ty = tmp.upperLeft();
    DestIterator              dy = dul;

    TMPTYPE thresh = detail::RequiresExplicitCast<TMPTYPE>::cast((gradient_threshold * gradient_threshold) *
                     NumericTraits<TMPTYPE>::one());

    // find zero crossings above threshold
    for(y=0; y<h-1; ++y, ++iy.y, ++ty.y, dy.y+=2)
    {
        typename TMPIMG::Iterator ix = iy;
        typename TMPIMG::Iterator tx = ty;
        DestIterator              dx = dy;

        for(int x=0; x<w-1; ++x, ++ix.x, ++tx.x, dx.x+=2)
        {
            TMPTYPE diff = *tx - *ix;
            TMPTYPE gx = tx[right] - *tx;
            TMPTYPE gy = tx[bottom] - *tx;

            if((gx * gx > thresh) &&
               (diff * (tx[right] - ix[right]) < zero))
            {
                da.set(edge_marker, dx, right);
            }
            if((gy * gy > thresh) &&
               (diff * (tx[bottom] - ix[bottom]) < zero))
            {
                da.set(edge_marker, dx, bottom);
            }
        }

        TMPTYPE diff = *tx - *ix;
        TMPTYPE gy = tx[bottom] - *tx;

        if((gy * gy > thresh) &&
           (diff * (tx[bottom] - ix[bottom]) < zero))
        {
            da.set(edge_marker, dx, bottom);
        }
    }

    typename TMPIMG::Iterator ix = iy;
    typename TMPIMG::Iterator tx = ty;
    DestIterator              dx = dy;

    for(x=0; x<w-1; ++x, ++ix.x, ++tx.x, dx.x+=2)
    {
        TMPTYPE diff = *tx - *ix;
        TMPTYPE gx = tx[right] - *tx;

        if((gx * gx > thresh) &&
           (diff * (tx[right] - ix[right]) < zero))
        {
            da.set(edge_marker, dx, right);
        }
    }

    iy = smooth.upperLeft() + Diff2D(0,1);
    ty = tmp.upperLeft() + Diff2D(0,1);
    dy = dul + Diff2D(1,2);

    const Diff2D topleft(-1,-1);
    const Diff2D topright(1,-1);
    const Diff2D bottomleft(-1,1);
    const Diff2D bottomright(1,1);

    // find missing 1-cells below threshold (x-direction)
    for(y=0; y<h-2; ++y, ++iy.y, ++ty.y, dy.y+=2)
    {
        typename TMPIMG::Iterator ix = iy;
        typename TMPIMG::Iterator tx = ty;
        DestIterator              dx = dy;

        for(int x=0; x<w-2; ++x, ++ix.x, ++tx.x, dx.x+=2)
        {
            if(da(dx) == edge_marker) continue;

            TMPTYPE diff = *tx - *ix;

            if((diff * (tx[right] - ix[right]) < zero) &&
               (((da(dx, bottomright) == edge_marker) &&
                 (da(dx, topleft) == edge_marker)) ||
                ((da(dx, bottomleft) == edge_marker) &&
                 (da(dx, topright) == edge_marker))))

            {
                da.set(edge_marker, dx);
            }
        }
    }

    iy = smooth.upperLeft() + Diff2D(1,0);
    ty = tmp.upperLeft() + Diff2D(1,0);
    dy = dul + Diff2D(2,1);

    // find missing 1-cells below threshold (y-direction)
    for(y=0; y<h-2; ++y, ++iy.y, ++ty.y, dy.y+=2)
    {
        typename TMPIMG::Iterator ix = iy;
        typename TMPIMG::Iterator tx = ty;
        DestIterator              dx = dy;

        for(int x=0; x<w-2; ++x, ++ix.x, ++tx.x, dx.x+=2)
        {
            if(da(dx) == edge_marker) continue;

            TMPTYPE diff = *tx - *ix;

            if((diff * (tx[bottom] - ix[bottom]) < zero) &&
               (((da(dx, bottomright) == edge_marker) &&
                 (da(dx, topleft) == edge_marker)) ||
                ((da(dx, bottomleft) == edge_marker) &&
                 (da(dx, topright) == edge_marker))))

            {
                da.set(edge_marker, dx);
            }
        }
    }

    dy = dul + Diff2D(1,1);

    // find missing 0-cells
    for(y=0; y<h-1; ++y, dy.y+=2)
    {
        DestIterator              dx = dy;

        for(int x=0; x<w-1; ++x, dx.x+=2)
        {
            const Diff2D dist[] = {right, top, left, bottom };

            int i;
            for(i=0; i<4; ++i)
            {
            if(da(dx, dist[i]) == edge_marker) break;
            }

            if(i < 4) da.set(edge_marker, dx);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
inline void
differenceOfExponentialCrackEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                      pair<DestIterator, DestAccessor> dest,
                                      double scale, GradValue gradient_threshold,
                                      DestValue edge_marker)
{
    differenceOfExponentialCrackEdgeImage(src.first, src.second, src.third,
                                          dest.first, dest.second,
                                          scale, gradient_threshold, edge_marker);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue, class DestValue>
inline void
differenceOfExponentialCrackEdgeImage(MultiArrayView<2, T1, S1> const & src,
                                      MultiArrayView<2, T2, S2> dest,
                                      double scale,
                                      GradValue gradient_threshold,
                                      DestValue edge_marker)
{
    vigra_precondition(2*src.shape() - Shape2(1) == dest.shape(),
        "differenceOfExponentialCrackEdgeImage(): shape mismatch between input and output.");
    differenceOfExponentialCrackEdgeImage(srcImageRange(src),
                                          destImage(dest),
                                          scale, gradient_threshold, edge_marker);
}

/********************************************************/
/*                                                      */
/*                  removeShortEdges                    */
/*                                                      */
/********************************************************/

/** \brief Remove short edges from an edge image.

    This algorithm can be applied as a post-processing operation of
    \ref differenceOfExponentialEdgeImage() and \ref differenceOfExponentialCrackEdgeImage().
    It removes all edges that are shorter than <TT>min_edge_length</TT>. The corresponding
    pixels are set to the <TT>non_edge_marker</TT>. The idea behind this algorithms is
    that very short edges are probably caused by noise and don't represent interesting
    image structure. Technically, the algorithms executes a connected components labeling,
    so the image's value type must be equality comparable.

    If the source image fulfills the requirements of a \ref CrackEdgeImage,
    it will still do so after application of this algorithm.

    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place,
    i.e. on only one image. Also, the algorithm assumes that all non-edges pixels are already
    marked with the given <TT>non_edge_marker</TT> value.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class Value>
        void
        removeShortEdges(MultiArrayView<2, T, S> image,
                         unsigned int min_edge_length, Value non_edge_marker);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        void removeShortEdges(
               Iterator sul, Iterator slr, Accessor sa,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        void removeShortEdges(
               triple<Iterator, Iterator, Accessor> src,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(w,h);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::differenceOfExponentialEdgeImage(srcImageRange(src), destImage(edges),
                                     0.8, 4.0, 1);

    // zero edges shorter than 10 pixels
    vigra::removeShortEdges(srcImageRange(edges), 10, 0);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u

    SrcValue non_edge_marker;
    src_accessor.set(non_edge_marker, src_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> void removeShortEdges)

template <class Iterator, class Accessor, class Value>
void removeShortEdges(
               Iterator sul, Iterator slr, Accessor sa,
           unsigned int min_edge_length, Value non_edge_marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    int x,y;

    IImage labels(w, h);
    labels = 0;

    int number_of_regions =
                labelImageWithBackground(srcIterRange(sul,slr,sa),
                                     destImage(labels), true, non_edge_marker);

    ArrayOfRegionStatistics<FindROISize<int> >
                                         region_stats(number_of_regions);

    inspectTwoImages(srcImageRange(labels), srcImage(labels), region_stats);

    IImage::Iterator ly = labels.upperLeft();
    Iterator oy = sul;

    for(y=0; y<h; ++y, ++oy.y, ++ly.y)
    {
        Iterator ox(oy);
        IImage::Iterator lx(ly);

        for(x=0; x<w; ++x, ++ox.x, ++lx.x)
        {
            if(sa(ox) == non_edge_marker) continue;
            if((region_stats[*lx].count) < min_edge_length)
            {
                 sa.set(non_edge_marker, ox);
            }
        }
    }
}

template <class Iterator, class Accessor, class Value>
inline void
removeShortEdges(triple<Iterator, Iterator, Accessor> src,
                 unsigned int min_edge_length, Value non_edge_marker)
{
    removeShortEdges(src.first, src.second, src.third,
                     min_edge_length, non_edge_marker);
}

template <class T, class S, class Value>
inline void
removeShortEdges(MultiArrayView<2, T, S> image,
                 unsigned int min_edge_length, Value non_edge_marker)
{
    removeShortEdges(destImageRange(image),
                     min_edge_length, non_edge_marker);
}

/********************************************************/
/*                                                      */
/*             closeGapsInCrackEdgeImage                */
/*                                                      */
/********************************************************/

/** \brief Close one-pixel wide gaps in a cell grid edge image.

    This algorithm is typically applied as a post-processing operation of
    \ref differenceOfExponentialCrackEdgeImage(). The source image must fulfill
    the requirements of a \ref CrackEdgeImage, and will still do so after
    application of this algorithm.

    It closes one pixel wide gaps in the edges resulting from this algorithm.
    Since these gaps are usually caused by zero crossing slightly below the gradient
    threshold used in edge detection, this algorithms acts like a weak hysteresis
    thresholding. The newly found edge pixels are marked with the given <TT>edge_marker</TT>.
    The image's value type must be equality comparable.

    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place,
    i.e. on only one image.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class Value>
        void
        closeGapsInCrackEdgeImage(MultiArrayView<2, T, S> image, Value edge_marker);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void closeGapsInCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               SrcValue edge_marker)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void closeGapsInCrackEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(2*w-1, 2*h-1);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::differenceOfExponentialCrackEdgeImage(srcImageRange(src), destImage(edges),
                                         0.8, 4.0, 1);

    // close gaps, mark with 1
    vigra::closeGapsInCrackEdgeImage(srcImageRange(edges), 1);

    // zero edges shorter than 20 pixels
    vigra::removeShortEdges(srcImageRange(edges), 10, 0);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u
    u != u

    SrcValue edge_marker;
    src_accessor.set(edge_marker, src_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> void closeGapsInCrackEdgeImage)

template <class SrcIterator, class SrcAccessor, class SrcValue>
void closeGapsInCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           SrcValue edge_marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    vigra_precondition(w % 2 == 1 && h % 2 == 1,
        "closeGapsInCrackEdgeImage(): Input is not a crack edge image (must have odd-numbered shape).");
        
    int w2 = w / 2, h2 = h / 2, x, y;

    int count1, count2, count3;

    const Diff2D right(1,0);
    const Diff2D bottom(0,1);
    const Diff2D left(-1,0);
    const Diff2D top(0,-1);

    const Diff2D leftdist[] = { Diff2D(0, 0), Diff2D(-1, 1), Diff2D(-2, 0), Diff2D(-1, -1)};
    const Diff2D rightdist[] = { Diff2D(2, 0), Diff2D(1, 1), Diff2D(0, 0), Diff2D(1, -1)};
    const Diff2D topdist[] = { Diff2D(1, -1), Diff2D(0, 0), Diff2D(-1, -1), Diff2D(0, -2)};
    const Diff2D bottomdist[] = { Diff2D(1, 1), Diff2D(0, 2), Diff2D(-1, 1), Diff2D(0, 0)};

    int i;

    SrcIterator sy = sul + Diff2D(0,1);
    SrcIterator sx;

    // close 1-pixel wide gaps (x-direction)
    for(y=0; y<h2; ++y, sy.y+=2)
    {
        sx = sy + Diff2D(2,0);

        for(x=2; x<w2; ++x, sx.x+=2)
        {
            if(sa(sx) == edge_marker) continue;

            if(sa(sx, left) != edge_marker) continue;
            if(sa(sx, right) != edge_marker) continue;

            count1 = 0;
            count2 = 0;
            count3 = 0;

            for(i=0; i<4; ++i)
            {
                if(sa(sx, leftdist[i]) == edge_marker)
                {
                    ++count1;
                    count3 ^= 1 << i;
                }
                if(sa(sx, rightdist[i]) == edge_marker)
                {
                    ++count2;
                    count3 ^= 1 << i;
                }
            }

            if(count1 <= 1 || count2 <= 1 || count3 == 15)
            {
                sa.set(edge_marker, sx);
            }
        }
   }

    sy = sul + Diff2D(1,2);

    // close 1-pixel wide gaps (y-direction)
    for(y=2; y<h2; ++y, sy.y+=2)
    {
        sx = sy;

        for(x=0; x<w2; ++x, sx.x+=2)
        {
            if(sa(sx) == edge_marker) continue;

            if(sa(sx, top) != edge_marker) continue;
            if(sa(sx, bottom) != edge_marker) continue;

            count1 = 0;
            count2 = 0;
            count3 = 0;

            for(i=0; i<4; ++i)
            {
                if(sa(sx, topdist[i]) == edge_marker)
                {
                    ++count1;
                    count3 ^= 1 << i;
                }
                if(sa(sx, bottomdist[i]) == edge_marker)
                {
                    ++count2;
                    count3 ^= 1 << i;
                }
            }

            if(count1 <= 1 || count2 <= 1 || count3 == 15)
            {
                sa.set(edge_marker, sx);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor, class SrcValue>
inline void
closeGapsInCrackEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          SrcValue edge_marker)
{
    closeGapsInCrackEdgeImage(src.first, src.second, src.third,
                              edge_marker);
}

template <class T, class S, class Value>
inline void
closeGapsInCrackEdgeImage(MultiArrayView<2, T, S> image, Value edge_marker)
{
    closeGapsInCrackEdgeImage(destImageRange(image), edge_marker);
}

/********************************************************/
/*                                                      */
/*              beautifyCrackEdgeImage                  */
/*                                                      */
/********************************************************/

/** \brief Beautify crack edge image for visualization.

    This algorithm is applied as a post-processing operation of
    \ref differenceOfExponentialCrackEdgeImage(). The source image must fulfill
    the requirements of a \ref CrackEdgeImage, but will <b> not</b> do so after
    application of this algorithm. In particular, the algorithm removes zero-cells
    marked as edges to avoid staircase effects on diagonal lines like this:

    \code
    original edge points (*)     resulting edge points

          . * . . .                   . * . . .
          . * * * .                   . . * . .
          . . . * .           =>      . . . * .
          . . . * *                   . . . . *
          . . . . .                   . . . . .
    \endcode

    Therefore, this algorithm should only be applied as a visualization aid, i.e.
    for human inspection. The algorithm assumes that edges are marked with <TT>edge_marker</TT>,
    and background pixels with <TT>background_marker</TT>. The image's value type must be
    equality comparable.

    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place,
    i.e. on only one image.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S, class Value>
        void
        beautifyCrackEdgeImage(MultiArrayView<2, T, S> image,
                               Value edge_marker, Value background_marker);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void beautifyCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               SrcValue edge_marker, SrcValue background_marker)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void beautifyCrackEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker, SrcValue background_marker)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(2*w-1, 2*h-1);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::differenceOfExponentialCrackEdgeImage(srcImageRange(src), destImage(edges),
                                         0.8, 4.0, 1);

    // beautify edge image for visualization
    vigra::beautifyCrackEdgeImage(destImageRange(edges), 1, 0);

    // show to the user
    window.open(edges);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u
    u != u

    SrcValue background_marker;
    src_accessor.set(background_marker, src_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> void beautifyCrackEdgeImage)

template <class SrcIterator, class SrcAccessor, class SrcValue>
void beautifyCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           SrcValue edge_marker, SrcValue background_marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    
    vigra_precondition(w % 2 == 1 && h % 2 == 1,
        "beautifyCrackEdgeImage(): Input is not a crack edge image (must have odd-numbered shape).");
        
    int w2 = w / 2, h2 = h / 2, x, y;

    SrcIterator sy = sul + Diff2D(1,1);
    SrcIterator sx;

    const Diff2D right(1,0);
    const Diff2D bottom(0,1);
    const Diff2D left(-1,0);
    const Diff2D top(0,-1);

    //  delete 0-cells at corners
    for(y=0; y<h2; ++y, sy.y+=2)
    {
        sx = sy;

        for(x=0; x<w2; ++x, sx.x+=2)
        {
            if(sa(sx) != edge_marker) continue;

            if(sa(sx, right) == edge_marker && sa(sx, left) == edge_marker) continue;
            if(sa(sx, bottom) == edge_marker && sa(sx, top) == edge_marker) continue;

            sa.set(background_marker, sx);
        }
    }
}

template <class SrcIterator, class SrcAccessor, class SrcValue>
inline void
beautifyCrackEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       SrcValue edge_marker, SrcValue background_marker)
{
    beautifyCrackEdgeImage(src.first, src.second, src.third,
                           edge_marker, background_marker);
}

template <class T, class S, class Value>
inline void
beautifyCrackEdgeImage(MultiArrayView<2, T, S> image,
                       Value edge_marker, Value background_marker)
{
    beautifyCrackEdgeImage(destImageRange(image),
                           edge_marker, background_marker);
}


/** Helper class that stores edgel attributes.
*/
class Edgel
{
  public:
  
        /** The type of an Edgel's members.
        */
    typedef float value_type;

        /** The edgel's sub-pixel x coordinate.
        */
    value_type x;

        /** The edgel's sub-pixel y coordinate.
        */
    value_type y;

        /** The edgel's strength (magnitude of the gradient vector).
        */
    value_type strength;

        /** \brief The edgel's orientation. 
        
        This is the clockwise angle in radians
        between the x-axis and the edge, so that the bright side of the
        edge is on the left when one looks along the orientation vector. 
        The angle is measured clockwise because the y-axis increases 
        downwards (left-handed coordinate system):

        \code

        edgel axis
             \  
        (dark \  (bright side)
        side)  \ 
                \ 
                 +------------> x-axis
                 |\    |
                 | \ /_/  orientation angle
                 |  \\
                 |   \
                 |   
          y-axis V  
        \endcode

        So, for example a vertical edge with its dark side on the left
        has orientation PI/2, and a horizontal edge with dark side on top
        has orientation PI. Obviously, the edge's orientation changes
        by PI if the contrast is reversed.
        
        Note that this convention changed as of VIGRA version 1.7.0.

        */
    value_type orientation;

    Edgel()
    : x(0.0), y(0.0), strength(0.0), orientation(0.0)
    {}

    Edgel(value_type ix, value_type iy, value_type is, value_type io)
    : x(ix), y(iy), strength(is), orientation(io)
    {}
};

template <class SrcIterator, class SrcAccessor, 
          class MagnitudeImage, class BackInsertable, class GradValue>
void internalCannyFindEdgels(SrcIterator ul, SrcAccessor grad,
                             MagnitudeImage const & magnitude,
                             BackInsertable & edgels, GradValue grad_thresh)
{
    typedef typename SrcAccessor::value_type PixelType;
    typedef typename PixelType::value_type ValueType;

    vigra_precondition(grad_thresh >= NumericTraits<GradValue>::zero(),
         "cannyFindEdgels(): gradient threshold must not be negative.");
    
    double t = 0.5 / VIGRA_CSTD::sin(M_PI/8.0);

    ul += Diff2D(1,1);
    for(int y=1; y<magnitude.height()-1; ++y, ++ul.y)
    {
        SrcIterator ix = ul;
        for(int x=1; x<magnitude.width()-1; ++x, ++ix.x)
        {
            double mag = magnitude(x, y);
            if(mag <= grad_thresh)
                   continue;
            ValueType gradx = grad.getComponent(ix, 0);
            ValueType grady = grad.getComponent(ix, 1);

            int dx = (int)VIGRA_CSTD::floor(gradx*t/mag + 0.5);
            int dy = (int)VIGRA_CSTD::floor(grady*t/mag + 0.5);

            int x1 = x - dx,
                x2 = x + dx,
                y1 = y - dy,
                y2 = y + dy;

            double m1 = magnitude(x1, y1);
            double m3 = magnitude(x2, y2);

            if(m1 < mag && m3 <= mag)
            {
                Edgel edgel;

                // local maximum => quadratic interpolation of sub-pixel location
                double del = 0.5 * (m1 - m3) / (m1 + m3 - 2.0*mag);
                edgel.x = Edgel::value_type(x + dx*del);
                edgel.y = Edgel::value_type(y + dy*del);
                edgel.strength = Edgel::value_type(mag);
                double orientation = VIGRA_CSTD::atan2(grady, gradx) + 0.5*M_PI;
                if(orientation < 0.0)
                    orientation += 2.0*M_PI;
                edgel.orientation = Edgel::value_type(orientation);
                edgels.push_back(edgel);
            }
        }
    }
}

/********************************************************/
/*                                                      */
/*                      cannyEdgelList                  */
/*                                                      */
/********************************************************/

/** \brief Simple implementation of Canny's edge detector.
    
    The function can be called in two modes: If you pass a 'scale', it is assumed that the 
    original image is scalar, and the Gaussian gradient is internally computed at the
    given 'scale'. If the function is called without scale parameter, it is assumed that
    the given image already contains the gradient (i.e. its value_type must be 
    a vector of length 2).

    On the basis of the gradient image, a simple non-maxima suppression is performed:
    for each 3x3 neighborhood, it is determined whether the center pixel has
    larger gradient magnitude than its two neighbors in gradient direction
    (where the direction is rounded into octants). If this is the case,
    a new \ref Edgel is appended to the given vector of <TT>edgels</TT>. The subpixel
    edgel position is determined by fitting a parabola to the three gradient 
    magnitude values mentioned above. The sub-pixel location of the parabola's tip
    and the gradient magnitude and direction (from the pixel center)
    are written in the newly created edgel.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class T, class S, class BackInsertable>
        void
        cannyEdgelList(MultiArrayView<2, T, S> const & src,
                       BackInsertable & edgels,
                       double scale);

        // compute edgels from a pre-computed gradient image
        template <class T, class S, class BackInsertable>
        void
        cannyEdgelList(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                       BackInsertable & edgels);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void 
        cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                       BackInsertable & edgels);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                       BackInsertable & edgels, double scale);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       BackInsertable & edgels);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       BackInsertable & edgels, double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h);

    // empty edgel list
    std::vector<vigra::Edgel> edgels;
    ...

    // find edgels at scale 0.8
    vigra::cannyEdgelList(srcImageRange(src), edgels, 0.8);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft;
    SrcAccessor src_accessor;

    src_accessor(src_upperleft);

    BackInsertable edgels;
    edgels.push_back(Edgel());
    \endcode

    SrcAccessor::value_type must be a type convertible to float

    <b> Preconditions:</b>

    \code
    scale > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgelList)

template <class SrcIterator, class SrcAccessor, class BackInsertable>
void 
cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
               BackInsertable & edgels, double scale)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TinyVector<TmpType, 2> > grad(lr-ul);
    gaussianGradient(srcIterRange(ul, lr, src), destImage(grad), scale);

    cannyEdgelList(srcImageRange(grad), edgels);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
void 
cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
               BackInsertable & edgels)
{
    using namespace functor;
    
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename NumericTraits<typename SrcType::value_type>::RealPromote TmpType;
    BasicImage<TmpType> magnitude(lr-ul);
    transformImage(srcIterRange(ul, lr, src), destImage(magnitude), norm(Arg1()));

    // find edgels
    internalCannyFindEdgels(ul, src, magnitude, edgels, NumericTraits<TmpType>::zero());
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
inline void
cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
               BackInsertable & edgels, double scale)
{
    cannyEdgelList(src.first, src.second, src.third, edgels, scale);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
inline void
cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
               BackInsertable & edgels)
{
    cannyEdgelList(src.first, src.second, src.third, edgels);
}

template <class T, class S, class BackInsertable>
inline void
cannyEdgelList(MultiArrayView<2, T, S> const & src,
               BackInsertable & edgels, double scale)
{
    cannyEdgelList(srcImageRange(src), edgels, scale);
}

template <class T, class S, class BackInsertable>
inline void
cannyEdgelList(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
               BackInsertable & edgels)
{
    cannyEdgelList(srcImageRange(src), edgels);
}

/********************************************************/
/*                                                      */
/*              cannyEdgelListThreshold                 */
/*                                                      */
/********************************************************/

/** \brief Canny's edge detector with thresholding.
    
    This function works exactly like \ref cannyEdgelList(), but 
    you also pass a threshold for the minimal gradient magnitude, 
    so that edgels whose strength is below the threshold are not 
    inserted into the edgel list.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class T, class S, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelListThreshold(MultiArrayView<2, T, S> const & src,
                                BackInsertable & edgels,
                                double scale,
                                GradValue grad_threshold);

        // compute edgels from a pre-computed gradient image
        template <class T, class S, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelListThreshold(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                                BackInsertable & edgels,
                                GradValue grad_threshold);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelListThreshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                BackInsertable & edgels, GradValue grad_threshold);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelListThreshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                BackInsertable & edgels, double scale, GradValue grad_threshold);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelListThreshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                BackInsertable & edgels, GradValue grad_threshold);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelListThreshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                BackInsertable & edgels, double scale, GradValue grad_threshold);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h);

    // empty edgel list
    std::vector<vigra::Edgel> edgels;
    ...

    // find edgels at scale 0.8, only considering gradient above 2.0
    vigra::cannyEdgelListThreshold(srcImageRange(src), edgels, 0.8, 2.0);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft;
    SrcAccessor src_accessor;

    src_accessor(src_upperleft);

    BackInsertable edgels;
    edgels.push_back(Edgel());
    \endcode

    SrcAccessor::value_type must be a type convertible to float

    <b> Preconditions:</b>

    \code
    scale > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgelListThreshold)

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
void 
cannyEdgelListThreshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                        BackInsertable & edgels, double scale, GradValue grad_threshold)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TinyVector<TmpType, 2> > grad(lr-ul);
    gaussianGradient(srcIterRange(ul, lr, src), destImage(grad), scale);

    cannyEdgelListThreshold(srcImageRange(grad), edgels, grad_threshold);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
void 
cannyEdgelListThreshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                        BackInsertable & edgels, GradValue grad_threshold)
{
    using namespace functor;
    
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename NumericTraits<typename SrcType::value_type>::RealPromote TmpType;
    BasicImage<TmpType> magnitude(lr-ul);
    transformImage(srcIterRange(ul, lr, src), destImage(magnitude), norm(Arg1()));

    // find edgels
    internalCannyFindEdgels(ul, src, magnitude, edgels, grad_threshold);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelListThreshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        BackInsertable & edgels, double scale, GradValue grad_threshold)
{
    cannyEdgelListThreshold(src.first, src.second, src.third, edgels, scale, grad_threshold);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelListThreshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        BackInsertable & edgels, GradValue grad_threshold)
{
    cannyEdgelListThreshold(src.first, src.second, src.third, edgels, grad_threshold);
}

template <class T, class S, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelListThreshold(MultiArrayView<2, T, S> const & src,
                        BackInsertable & edgels,
                        double scale,
                        GradValue grad_threshold)
{
    cannyEdgelListThreshold(srcImageRange(src), edgels, scale, grad_threshold);
}

template <class T, class S, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelListThreshold(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                        BackInsertable & edgels,
                        GradValue grad_threshold)
{
    cannyEdgelListThreshold(srcImageRange(src), edgels, grad_threshold);
}


/********************************************************/
/*                                                      */
/*                       cannyEdgeImage                 */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in an edge image using Canny's algorithm.

    This operator first calls \ref cannyEdgelList() with the given scale to generate an
    edgel list for the given image. Then it scans this list and selects edgels
    whose strength is above the given <TT>gradient_threshold</TT>. For each of these
    edgels, the edgel's location is rounded to the nearest pixel, and that
    pixel marked with the given <TT>edge_marker</TT>.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class GradValue, class DestValue>
        void 
        cannyEdgeImage(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                       double scale,
                       GradValue gradient_threshold,
                       DestValue edge_marker);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImage(
                   SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                   DestIterator dul, DestAccessor da,
                   double scale, GradValue gradient_threshold, DestValue edge_marker);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   double scale, GradValue gradient_threshold, DestValue edge_marker);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(w,h);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1
    vigra::cannyEdgeImage(srcImageRange(src), destImage(edges),
                                     0.8, 4.0, 1);
    \endcode

    <b> Required Interface:</b>

    see also: \ref cannyEdgelList().

    \code
    DestImageIterator dest_upperleft;
    DestAccessor dest_accessor;
    DestValue edge_marker;

    dest_accessor.set(edge_marker, dest_upperleft, vigra::Diff2D(1,1));
    \endcode

    <b> Preconditions:</b>

    \code
    scale > 0
    gradient_threshold > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgeImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void cannyEdgeImage(
           SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    std::vector<Edgel> edgels;

    cannyEdgelListThreshold(sul, slr, sa, edgels, scale, gradient_threshold);
    
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    for(unsigned int i=0; i<edgels.size(); ++i)
    {
        Diff2D pix((int)(edgels[i].x + 0.5), (int)(edgels[i].y + 0.5));
        
        if(pix.x < 0 || pix.x >= w || pix.y < 0 || pix.y >= h)
            continue;

        da.set(edge_marker, dul, pix);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
inline void 
cannyEdgeImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    cannyEdgeImage(src.first, src.second, src.third,
                   dest.first, dest.second,
                   scale, gradient_threshold, edge_marker);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue, class DestValue>
inline void 
cannyEdgeImage(MultiArrayView<2, T1, S1> const & src,
               MultiArrayView<2, T2, S2> dest,
               double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    vigra_precondition(src.shape() == dest.shape(),
        "cannyEdgeImage(): shape mismatch between input and output.");
    cannyEdgeImage(srcImageRange(src),
                   destImage(dest),
                   scale, gradient_threshold, edge_marker);
}

/********************************************************/

namespace detail {

template <class DestIterator>
int neighborhoodConfiguration(DestIterator dul)
{
    int v = 0;
    NeighborhoodCirculator<DestIterator, EightNeighborCode> c(dul, EightNeighborCode::SouthEast);
    for(int i=0; i<8; ++i, --c)
    {
        v = (v << 1) | ((*c != 0) ? 1 : 0);
    }

    return v;
}

template <class GradValue>
struct SimplePoint
{
    Diff2D point;
    GradValue grad;

    SimplePoint(Diff2D const & p, GradValue g)
    : point(p), grad(g)
    {}

    bool operator<(SimplePoint const & o) const
    {
        return grad < o.grad;
    }

    bool operator>(SimplePoint const & o) const
    {
        return grad > o.grad;
    }
};

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void cannyEdgeImageFromGrad(
           SrcIterator sul, SrcIterator slr, SrcAccessor grad,
           DestIterator dul, DestAccessor da,
           GradValue gradient_threshold, DestValue edge_marker)
{
    typedef typename SrcAccessor::value_type PixelType;
    typedef typename NormTraits<PixelType>::SquaredNormType NormType;

    NormType zero = NumericTraits<NormType>::zero();
    double tan22_5 = M_SQRT2 - 1.0;
    typename NormTraits<GradValue>::SquaredNormType g2thresh = squaredNorm(gradient_threshold);

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    sul += Diff2D(1,1);
    dul += Diff2D(1,1);
    Diff2D p(0,0);

    for(int y = 1; y < h-1; ++y, ++sul.y, ++dul.y)
    {
        SrcIterator sx = sul;
        DestIterator dx = dul;
        for(int x = 1; x < w-1; ++x, ++sx.x, ++dx.x)
        {
            PixelType g = grad(sx);
            NormType g2n = squaredNorm(g);
            if(g2n < g2thresh)
                continue;

            NormType g2n1, g2n3;
            // find out quadrant
            if(abs(g[1]) < tan22_5*abs(g[0]))
            {
                // north-south edge
                g2n1 = squaredNorm(grad(sx, Diff2D(-1, 0)));
                g2n3 = squaredNorm(grad(sx, Diff2D(1, 0)));
            }
            else if(abs(g[0]) < tan22_5*abs(g[1]))
            {
                // west-east edge
                g2n1 = squaredNorm(grad(sx, Diff2D(0, -1)));
                g2n3 = squaredNorm(grad(sx, Diff2D(0, 1)));
            }
            else if(g[0]*g[1] < zero)
            {
                // north-west-south-east edge
                g2n1 = squaredNorm(grad(sx, Diff2D(1, -1)));
                g2n3 = squaredNorm(grad(sx, Diff2D(-1, 1)));
            }
            else
            {
                // north-east-south-west edge
                g2n1 = squaredNorm(grad(sx, Diff2D(-1, -1)));
                g2n3 = squaredNorm(grad(sx, Diff2D(1, 1)));
            }

            if(g2n1 < g2n && g2n3 <= g2n)
            {
                da.set(edge_marker, dx);
            }
        }
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*          cannyEdgeImageFromGradWithThinning          */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in an edge image using Canny's algorithm.

    The input pixels of this algorithm must be vectors of length 2 (see Required Interface below).
    It first searches for all pixels whose gradient magnitude is larger
    than the given <tt>gradient_threshold</tt> and larger than the magnitude of its two neighbors
    in gradient direction (where these neighbors are determined by nearest neighbor
    interpolation, i.e. according to the octant where the gradient points into).
    The resulting edge pixel candidates are then subjected to topological thinning
    so that the remaining edge pixels can be linked into edgel chains with a provable,
    non-heuristic algorithm. Thinning is performed so that the pixels with highest gradient
    magnitude survive. Optionally, the outermost pixels are marked as edge pixels
    as well when <tt>addBorder</tt> is true. The remaining pixels will be marked in the destination
    image with the value of <tt>edge_marker</tt> (all non-edge pixels remain untouched).

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class GradValue, class DestValue>
        void 
        cannyEdgeImageFromGradWithThinning(MultiArrayView<2, TinyVector<T1, 2>, S1> const & src,
                                           MultiArrayView<2, T2, S2> dest,
                                           GradValue gradient_threshold,
                                           DestValue edge_marker,
                                           bool addBorder = true);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImageFromGradWithThinning(
                   SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                   DestIterator dul, DestAccessor da,
                   GradValue gradient_threshold,
                   DestValue edge_marker, bool addBorder = true);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImageFromGradWithThinning(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   GradValue gradient_threshold,
                   DestValue edge_marker, bool addBorder = true);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(w,h);

    vigra::FVector2Image grad(w,h);
    // compute the image gradient at scale 0.8
    vigra::gaussianGradient(srcImageRange(src), destImage(grad), 0.8);

    // empty edge image
    edges = 0;
    // find edges gradient larger than 4.0, mark with 1, and add border
    vigra::cannyEdgeImageFromGradWithThinning(srcImageRange(grad), destImage(edges),
                                              4.0, 1, true);
    \endcode

    <b> Required Interface:</b>

    \code
    // the input pixel type must be a vector with two elements
    SrcImageIterator src_upperleft;
    SrcAccessor src_accessor;
    typedef SrcAccessor::value_type SrcPixel;
    typedef NormTraits<SrcPixel>::SquaredNormType SrcSquaredNormType;

    SrcPixel g = src_accessor(src_upperleft);
    SrcPixel::value_type g0 = g[0];
    SrcSquaredNormType gn = squaredNorm(g);

    DestImageIterator dest_upperleft;
    DestAccessor dest_accessor;
    DestValue edge_marker;

    dest_accessor.set(edge_marker, dest_upperleft, vigra::Diff2D(1,1));
    \endcode

    <b> Preconditions:</b>

    \code
    gradient_threshold > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgeImageFromGradWithThinning)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void cannyEdgeImageFromGradWithThinning(
           SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           GradValue gradient_threshold,
           DestValue edge_marker, bool addBorder = true)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    BImage edgeImage(w, h, BImage::value_type(0));
    BImage::traverser eul = edgeImage.upperLeft();
    BImage::Accessor ea = edgeImage.accessor();
    if(addBorder)
        initImageBorder(destImageRange(edgeImage), 1, 1);
    detail::cannyEdgeImageFromGrad(sul, slr, sa, eul, ea, gradient_threshold, 1);

    bool isSimplePoint[256] = {
        0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
        1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1,
        0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0,
        0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
        1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
        0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0,
        1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
        0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1,
        1, 0, 1, 0 };

    eul += Diff2D(1,1);
    sul += Diff2D(1,1);
    int w2 = w-2;
    int h2 = h-2;

    typedef detail::SimplePoint<GradValue> SP;
    // use std::greater because we need the smallest gradients at the top of the queue
    std::priority_queue<SP, std::vector<SP>, std::greater<SP> >  pqueue;

    Diff2D p(0,0);
    for(; p.y < h2; ++p.y)
    {
        for(p.x = 0; p.x < w2; ++p.x)
        {
            BImage::traverser e = eul + p;
            if(*e == 0)
                continue;
            int v = detail::neighborhoodConfiguration(e);
            if(isSimplePoint[v])
            {
                pqueue.push(SP(p, norm(sa(sul+p))));
                *e = 2; // remember that it is already in queue
            }
        }
    }

    const Diff2D dist[] = { Diff2D(-1,0), Diff2D(0,-1), Diff2D(1,0),  Diff2D(0,1) };

    while(pqueue.size())
    {
        p = pqueue.top().point;
        pqueue.pop();

        BImage::traverser e = eul + p;
        int v = detail::neighborhoodConfiguration(e);
        if(!isSimplePoint[v])
            continue; // point may no longer be simple because its neighbors changed

        *e = 0; // delete simple point

        for(int i=0; i<4; ++i)
        {
            Diff2D pneu = p + dist[i];
            if(pneu.x == -1 || pneu.y == -1 || pneu.x == w2 || pneu.y == h2)
                continue; // do not remove points at the border

            BImage::traverser eneu = eul + pneu;
            if(*eneu == 1) // point is boundary and not yet in the queue
            {
                int v = detail::neighborhoodConfiguration(eneu);
                if(isSimplePoint[v])
                {
                    pqueue.push(SP(pneu, norm(sa(sul+pneu))));
                    *eneu = 2; // remember that it is already in queue
                }
            }
        }
    }

    initImageIf(destIterRange(dul, dul+Diff2D(w,h), da),
                maskImage(edgeImage), edge_marker);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
inline void cannyEdgeImageFromGradWithThinning(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest,
           GradValue gradient_threshold,
           DestValue edge_marker, bool addBorder = true)
{
    cannyEdgeImageFromGradWithThinning(src.first, src.second, src.third,
                               dest.first, dest.second,
                               gradient_threshold, edge_marker, addBorder);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue, class DestValue>
inline void 
cannyEdgeImageFromGradWithThinning(MultiArrayView<2, TinyVector<T1, 2>, S1> const & src,
                                   MultiArrayView<2, T2, S2> dest,
                                   GradValue gradient_threshold,
                                   DestValue edge_marker, bool addBorder = true)
{
    vigra_precondition(src.shape() == dest.shape(),
        "cannyEdgeImageFromGradWithThinning(): shape mismatch between input and output.");
    cannyEdgeImageFromGradWithThinning(srcImageRange(src),
                                       destImage(dest),
                                       gradient_threshold, edge_marker, addBorder);
}

/********************************************************/
/*                                                      */
/*              cannyEdgeImageWithThinning              */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in an edge image using Canny's algorithm.

    This operator first calls \ref gaussianGradient() to compute the gradient of the input
    image, ad then \ref cannyEdgeImageFromGradWithThinning() to generate an
    edge image. See there for more detailed documentation.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class GradValue, class DestValue>
        void 
        cannyEdgeImageWithThinning(MultiArrayView<2, T1, S1> const & src,
                                   MultiArrayView<2, T2, S2> dest,
                                   double scale, 
                                   GradValue gradient_threshold,
                                   DestValue edge_marker, 
                                   bool addBorder = true);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImageWithThinning(
                   SrcIterator sul, SrcIterator slr, SrcAccessor sa,
                   DestIterator dul, DestAccessor da,
                   double scale, GradValue gradient_threshold,
                   DestValue edge_marker, bool addBorder = true);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class GradValue, class DestValue>
        void cannyEdgeImageWithThinning(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   double scale, GradValue gradient_threshold,
                   DestValue edge_marker, bool addBorder = true);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h), edges(w,h);

    // empty edge image
    edges = 0;
    ...

    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1, annd add border
    vigra::cannyEdgeImageWithThinning(srcImageRange(src), destImage(edges),
                                     0.8, 4.0, 1, true);
    \endcode

    <b> Required Interface:</b>

    see also: \ref cannyEdgelList().

    \code
    DestImageIterator dest_upperleft;
    DestAccessor dest_accessor;
    DestValue edge_marker;

    dest_accessor.set(edge_marker, dest_upperleft, vigra::Diff2D(1,1));
    \endcode

    <b> Preconditions:</b>

    \code
    scale > 0
    gradient_threshold > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgeImageWithThinning)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
void cannyEdgeImageWithThinning(
           SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold,
           DestValue edge_marker, bool addBorder = true)
{
    // mark pixels that are higher than their neighbors in gradient direction
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TinyVector<TmpType, 2> > grad(slr-sul);
    gaussianGradient(srcIterRange(sul, slr, sa), destImage(grad), scale);
    cannyEdgeImageFromGradWithThinning(srcImageRange(grad), destIter(dul, da),
                               gradient_threshold, edge_marker, addBorder);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue, class DestValue>
inline void 
cannyEdgeImageWithThinning(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           pair<DestIterator, DestAccessor> dest,
                           double scale, GradValue gradient_threshold,
                           DestValue edge_marker, bool addBorder = true)
{
    cannyEdgeImageWithThinning(src.first, src.second, src.third,
                               dest.first, dest.second,
                               scale, gradient_threshold, edge_marker, addBorder);
}

template <class T1, class S1,
          class T2, class S2,
          class GradValue, class DestValue>
inline void 
cannyEdgeImageWithThinning(MultiArrayView<2, T1, S1> const & src,
                           MultiArrayView<2, T2, S2> dest,
                           double scale, GradValue gradient_threshold,
                           DestValue edge_marker, bool addBorder = true)
{
    vigra_precondition(src.shape() == dest.shape(),
        "cannyEdgeImageWithThinning(): shape mismatch between input and output.");
    cannyEdgeImageWithThinning(srcImageRange(src),
                               destImage(dest),
                               scale, gradient_threshold, edge_marker, addBorder);
}

/********************************************************/

template <class SrcIterator, class SrcAccessor, 
          class MaskImage, class BackInsertable, class GradValue>
void internalCannyFindEdgels3x3(SrcIterator ul, SrcAccessor grad,
                                MaskImage const & mask,
                                BackInsertable & edgels,
                                GradValue grad_thresh)
{
    typedef typename SrcAccessor::value_type PixelType;
    typedef typename PixelType::value_type ValueType;

    vigra_precondition(grad_thresh >= NumericTraits<GradValue>::zero(),
         "cannyFindEdgels3x3(): gradient threshold must not be negative.");
    
    ul += Diff2D(1,1);
    for(int y=1; y<mask.height()-1; ++y, ++ul.y)
    {
        SrcIterator ix = ul;
        for(int x=1; x<mask.width()-1; ++x, ++ix.x)
        {
            if(!mask(x,y))
                continue;

            ValueType gradx = grad.getComponent(ix, 0);
            ValueType grady = grad.getComponent(ix, 1);
            double mag = hypot(gradx, grady);
            if(mag <= grad_thresh)
                   continue;
            double c = gradx / mag,
                   s = grady / mag;

            Matrix<double> ml(3,3), mr(3,1), l(3,1), r(3,1);
            l(0,0) = 1.0;

            for(int yy = -1; yy <= 1; ++yy)
            {
                for(int xx = -1; xx <= 1; ++xx)
                {
                    double u = c*xx + s*yy;
                    double v = norm(grad(ix, Diff2D(xx, yy)));
                    l(1,0) = u;
                    l(2,0) = u*u;
                    ml += outer(l);
                    mr += v*l;
                }
            }

            linearSolve(ml, mr, r);

            Edgel edgel;

            // local maximum => quadratic interpolation of sub-pixel location
            double del = -r(1,0) / 2.0 / r(2,0);
            if(std::fabs(del) > 1.5)  // don't move by more than about a pixel diameter
                del = 0.0;
            edgel.x = Edgel::value_type(x + c*del);
            edgel.y = Edgel::value_type(y + s*del);
            edgel.strength = Edgel::value_type(mag);
            double orientation = VIGRA_CSTD::atan2(grady, gradx) + 0.5*M_PI;
            if(orientation < 0.0)
                orientation += 2.0*M_PI;
            edgel.orientation = Edgel::value_type(orientation);
            edgels.push_back(edgel);
        }
    }
}


/********************************************************/
/*                                                      */
/*                   cannyEdgelList3x3                  */
/*                                                      */
/********************************************************/

/** \brief Improved implementation of Canny's edge detector.

    This operator first computes pixels which are crossed by the edge using
    cannyEdgeImageWithThinning(). The gradient magnitudes in the 3x3 neighborhood of these
    pixels are then projected onto the normal of the edge (as determined
    by the gradient direction). The edgel's subpixel location is found by fitting a
    parabola through the 9 gradient values and determining the parabola's tip.
    A new \ref Edgel is appended to the given vector of <TT>edgels</TT>. Since the parabola
    is fitted to 9 points rather than 3 points as in cannyEdgelList(), the accuracy is higher.
    
    The function can be called in two modes: If you pass a 'scale', it is assumed that the 
    original image is scalar, and the Gaussian gradient is internally computed at the
    given 'scale'. If the function is called without scale parameter, it is assumed that
    the given image already contains the gradient (i.e. its value_type must be 
    a vector of length 2).

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class T, class S, class BackInsertable>
        void
        cannyEdgelList3x3(MultiArrayView<2, T, S> const & src,
                          BackInsertable & edgels, 
                          double scale);

        // compute edgels from a pre-computed gradient image
        template <class T, class S, class BackInsertable>
        void
        cannyEdgelList3x3(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                          BackInsertable & edgels);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void cannyEdgelList3x3(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                               BackInsertable & edgels);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void cannyEdgelList3x3(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                               BackInsertable & edgels, double scale);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        cannyEdgelList3x3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          BackInsertable & edgels);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, class BackInsertable>
        void
        cannyEdgelList3x3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          BackInsertable & edgels, double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h);

    // empty edgel list
    std::vector<vigra::Edgel> edgels;
    ...

    // find edgels at scale 0.8
    vigra::cannyEdgelList3x3(srcImageRange(src), edgels, 0.8);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft;
    SrcAccessor src_accessor;

    src_accessor(src_upperleft);

    BackInsertable edgels;
    edgels.push_back(Edgel());
    \endcode

    SrcAccessor::value_type must be a type convertible to float

    <b> Preconditions:</b>

    \code
    scale > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgelList3x3)

template <class SrcIterator, class SrcAccessor, class BackInsertable>
void 
cannyEdgelList3x3(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                  BackInsertable & edgels, double scale)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TinyVector<TmpType, 2> > grad(lr-ul);
    gaussianGradient(srcIterRange(ul, lr, src), destImage(grad), scale);

    cannyEdgelList3x3(srcImageRange(grad), edgels);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
void 
cannyEdgelList3x3(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                  BackInsertable & edgels)
{
    typedef typename NormTraits<typename SrcAccessor::value_type>::NormType NormType;
    
    UInt8Image edges(lr-ul);
    cannyEdgeImageFromGradWithThinning(srcIterRange(ul, lr, src), destImage(edges),
                                       0.0, 1, false);

    // find edgels
    internalCannyFindEdgels3x3(ul, src, edges, edgels, NumericTraits<NormType>::zero());
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
inline void
cannyEdgelList3x3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  BackInsertable & edgels, double scale)
{
    cannyEdgelList3x3(src.first, src.second, src.third, edgels, scale);
}

template <class SrcIterator, class SrcAccessor, class BackInsertable>
inline void
cannyEdgelList3x3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  BackInsertable & edgels)
{
    cannyEdgelList3x3(src.first, src.second, src.third, edgels);
}

template <class T, class S, class BackInsertable>
inline void
cannyEdgelList3x3(MultiArrayView<2, T, S> const & src,
                  BackInsertable & edgels, double scale)
{
    cannyEdgelList3x3(srcImageRange(src), edgels, scale);
}

template <class T, class S, class BackInsertable>
inline void
cannyEdgelList3x3(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                  BackInsertable & edgels)
{
    cannyEdgelList3x3(srcImageRange(src), edgels);
}

/********************************************************/
/*                                                      */
/*             cannyEdgelList3x3Threshold               */
/*                                                      */
/********************************************************/

/** \brief Improved implementation of Canny's edge detector with thresholding.

    This function works exactly like \ref cannyEdgelList3x3(), but 
    you also pass a threshold for the minimal gradient magnitude, 
    so that edgels whose strength is below the threshold are not 
    inserted into the edgel list.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                   BackInsertable & edgels, GradValue grad_thresh);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                   BackInsertable & edgels, double scale, GradValue grad_thresh);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                   BackInsertable & edgels, GradValue grad_thresh);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void 
        cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                   BackInsertable & edgels, double scale, GradValue grad_thresh);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // compute edgels from a gradient image
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelList3x3Threshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   BackInsertable & edgels, GradValue grad_thresh);

        // compute edgels from a scalar image (determine gradient internally at 'scale')
        template <class SrcIterator, class SrcAccessor, 
                  class BackInsertable, class GradValue>
        void
        cannyEdgelList3x3Threshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   BackInsertable & edgels, double scale, GradValue grad_thresh);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/edgedetection.hxx\><br>
    Namespace: vigra

    \code
    vigra::BImage src(w,h);

    // empty edgel list
    std::vector<vigra::Edgel> edgels;
    ...

    // find edgels at scale 0.8
    vigra::cannyEdgelList3x3(srcImageRange(src), edgels, 0.8);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft;
    SrcAccessor src_accessor;

    src_accessor(src_upperleft);

    BackInsertable edgels;
    edgels.push_back(Edgel());
    \endcode

    SrcAccessor::value_type must be a type convertible to float

    <b> Preconditions:</b>

    \code
    scale > 0
    \endcode
*/
doxygen_overloaded_function(template <...> void cannyEdgelList3x3Threshold)

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
void 
cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                           BackInsertable & edgels, double scale, GradValue grad_thresh)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TinyVector<TmpType, 2> > grad(lr-ul);
    gaussianGradient(srcIterRange(ul, lr, src), destImage(grad), scale);

    cannyEdgelList3x3Threshold(srcImageRange(grad), edgels, grad_thresh);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
void 
cannyEdgelList3x3Threshold(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                           BackInsertable & edgels, GradValue grad_thresh)
{
    UInt8Image edges(lr-ul);
    cannyEdgeImageFromGradWithThinning(srcIterRange(ul, lr, src), destImage(edges),
                                       0.0, 1, false);

    // find edgels
    internalCannyFindEdgels3x3(ul, src, edges, edgels, grad_thresh);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelList3x3Threshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           BackInsertable & edgels, double scale, GradValue grad_thresh)
{
    cannyEdgelList3x3Threshold(src.first, src.second, src.third, edgels, scale, grad_thresh);
}

template <class SrcIterator, class SrcAccessor, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelList3x3Threshold(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           BackInsertable & edgels, GradValue grad_thresh)
{
    cannyEdgelList3x3Threshold(src.first, src.second, src.third, edgels, grad_thresh);
}

template <class T, class S, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelList3x3Threshold(MultiArrayView<2, T, S> const & src,
                           BackInsertable & edgels, double scale, GradValue grad_thresh)
{
    cannyEdgelList3x3Threshold(srcImageRange(src), edgels, scale, grad_thresh);
}

template <class T, class S, 
          class BackInsertable, class GradValue>
inline void
cannyEdgelList3x3Threshold(MultiArrayView<2, TinyVector<T, 2>, S> const & src,
                           BackInsertable & edgels,
                           GradValue grad_thresh)
{
    cannyEdgelList3x3Threshold(srcImageRange(src), edgels, grad_thresh);
}

//@}

/** \page CrackEdgeImage Crack Edge Image

Crack edges are marked <i>between</i> the pixels of an image.
A Crack Edge Image is an image that represents these edges. In order
to accommodate the cracks, the Crack Edge Image must be twice as large
as the original image (precisely (2*w - 1) by (2*h - 1)). A Crack Edge Image
can easily be derived from a binary image or from the signs of the
response of a Laplacian filter. Consider the following sketch, where
<TT>+</TT> encodes the foreground, <TT>-</TT> the background, and
<TT>*</TT> the resulting crack edges.

    \code
sign of difference image         insert cracks         resulting CrackEdgeImage

                                   + . - . -              . * . . .
      + - -                        . . . . .              . * * * .
      + + -               =>       + . + . -      =>      . . . * .
      + + +                        . . . . .              . . . * *
                                   + . + . +              . . . . .
    \endcode

Starting from the original binary image (left), we insert crack pixels
to get to the double-sized image (center). Finally, we mark all
crack pixels whose non-crack neighbors have different signs as
crack edge points, while all other pixels (crack and non-crack) become
region pixels.

<b>Requirements on a Crack Edge Image:</b>

<ul>
    <li>Crack Edge Images have odd width and height.
    <li>Crack pixels have at least one odd coordinate.
    <li>Only crack pixels may be marked as edge points.
    <li>Crack pixels with two odd coordinates must be marked as edge points
        whenever any of their neighboring crack pixels was marked.
</ul>

The last two requirements ensure that both edges and regions are 4-connected.
Thus, 4-connectivity and 8-connectivity yield identical connected
components in a Crack Edge Image (so called <i>well-composedness</i>).
This ensures that Crack Edge Images have nice topological properties
(cf. L. J. Latecki: "Well-Composed Sets", Academic Press, 2000).
*/


} // namespace vigra

#endif // VIGRA_EDGEDETECTION_HXX
