/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_EDGEDETECTION_HXX
#define VIGRA_EDGEDETECTION_HXX

#include <vector>
#include <cmath>     // sqrt(), abs()
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/separableconvolution.hxx"
#include "vigra/labelimage.hxx"


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
    
    pass arguments explicitly:
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
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor, 
              class GradValue,
              class DestValue = DestAccessor::value_type>
        inline 
        void differenceOfExponentialEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
    
    static const Diff2D right(1, 0);
    static const Diff2D bottom(0, 1);
    
    
    TMPTYPE thresh = (gradient_threshold * gradient_threshold) * 
                     NumericTraits<TMPTYPE>::one();
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
inline 
void differenceOfExponentialEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, GradValue gradient_threshold,
       DestValue edge_marker)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold,
                    edge_marker);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class GradValue>
inline 
void differenceOfExponentialEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, GradValue gradient_threshold)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold, 1);
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
    <i>between</i>} the corresponding original pixels. Topologically, this means we 
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
    The result conformes to the requirements of a \ref CrackEdgeImage. It can be further
    improved by the post-processing operations \ref removeShortEdges() and
    \ref closeGapsInCrackEdgeImage().
    
    The source value type (<TT>SrcAccessor::value_type</TT>) must be a linear algebra, i.e. addition, 
    subtraction and multiplication of the type with itself, and multiplication 
    with double and 
    \ref NumericTraits "NumericTraits" must 
    be defined. In addition, this type must be less-comparable.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
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
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor, 
              class GradValue,
              class DestValue = DestAccessor::value_type>
        inline 
        void differenceOfExponentialCrackEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
    
    static const Diff2D right(1,0);
    static const Diff2D bottom(0,1);
    static const Diff2D left(-1,0);
    static const Diff2D top(0,-1);
    
    recursiveSmoothX(srcIterRange(sul, slr, sa), destImage(tmp), scale / 2.0);
    recursiveSmoothY(srcImageRange(tmp), destImage(tmp), scale / 2.0);

    recursiveSmoothX(srcImageRange(tmp), destImage(smooth), scale);
    recursiveSmoothY(srcImageRange(smooth), destImage(smooth), scale);

    typename TMPIMG::Iterator iy = smooth.upperLeft();
    typename TMPIMG::Iterator ty = tmp.upperLeft();
    DestIterator              dy = dul;
    
    TMPTYPE thresh = (gradient_threshold * gradient_threshold) * 
                     NumericTraits<TMPTYPE>::one();

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
    
    static const Diff2D topleft(-1,-1);
    static const Diff2D topright(1,-1);
    static const Diff2D bottomleft(-1,1);
    static const Diff2D bottomright(1,1);

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
            static const Diff2D dist[] = {right, top, left, bottom };

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
inline 
void differenceOfExponentialCrackEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, GradValue gradient_threshold,
       DestValue edge_marker)
{
    differenceOfExponentialCrackEdgeImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold,
                    edge_marker);
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
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        void removeShortEdges(
               Iterator sul, Iterator slr, Accessor sa,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        inline 
        void removeShortEdges(
               triple<Iterator, Iterator, Accessor> src,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
inline 
void removeShortEdges(
           triple<Iterator, Iterator, Accessor> src,
       unsigned int min_edge_length, Value non_edge_marker)
{
    removeShortEdges(src.first, src.second, src.third,
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
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void closeGapsInCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
               SrcValue edge_marker)
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        inline
        void closeGapsInCrackEdgeImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
template <class SrcIterator, class SrcAccessor, class SrcValue>
void closeGapsInCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
           SrcValue edge_marker)
{
    int w = (slr.x - sul.x) / 2;
    int h = (slr.y - sul.y) / 2;
    int x, y;

    int count1, count2, count3;

    static const Diff2D right(1,0);
    static const Diff2D bottom(0,1);
    static const Diff2D left(-1,0);
    static const Diff2D top(0,-1);
    
    static const Diff2D leftdist[] = { 
        Diff2D(0, 0), Diff2D(-1, 1), Diff2D(-2, 0), Diff2D(-1, -1)};
    static const Diff2D rightdist[] = { 
        Diff2D(2, 0), Diff2D(1, 1), Diff2D(0, 0), Diff2D(1, -1)};
    static const Diff2D topdist[] = { 
        Diff2D(1, -1), Diff2D(0, 0), Diff2D(-1, -1), Diff2D(0, -2)};
    static const Diff2D bottomdist[] = { 
        Diff2D(1, 1), Diff2D(0, 2), Diff2D(-1, 1), Diff2D(0, 0)};

    int i;

    SrcIterator sy = sul + Diff2D(0,1);
    SrcIterator sx;
    
    // close 1-pixel wide gaps (x-direction)
    for(y=0; y<h; ++y, sy.y+=2)
    {
        sx = sy + Diff2D(2,0);

        for(x=2; x<w; ++x, sx.x+=2)
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
    for(y=2; y<h; ++y, sy.y+=2)
    {
        sx = sy;

        for(x=0; x<w; ++x, sx.x+=2)
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
inline
void closeGapsInCrackEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       SrcValue edge_marker)
{
    closeGapsInCrackEdgeImage(src.first, src.second, src.third,
                    edge_marker);
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
    
    Therfore, this algorithm should only be applied as a vizualization aid, i.e. 
    for human inspection. The algorithm assumes that edges are marked with <TT>edge_marker</TT>, 
    and background pixels with <TT>background_marker</TT>. The image's value type must be 
    equality comparable. 
    
    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place, 
    i.e. on only one image.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void beautifyCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
               SrcValue edge_marker, SrcValue background_marker)
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        inline
        void beautifyCrackEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker, SrcValue background_marker)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
    vigra::beautifyCrackEdgeImage(srcImageRange(edges), 1, 0);
    
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
template <class SrcIterator, class SrcAccessor, class SrcValue>
void beautifyCrackEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
           SrcValue edge_marker, SrcValue background_marker)
{
    int w = (slr.x - sul.x) / 2;
    int h = (slr.y - sul.y) / 2;
    int x, y;

    SrcIterator sy = sul + Diff2D(1,1);
    SrcIterator sx;
    
    static const Diff2D right(1,0);
    static const Diff2D bottom(0,1);
    static const Diff2D left(-1,0);
    static const Diff2D top(0,-1);
    
    //  delete 0-cells at corners
    for(y=0; y<h; ++y, sy.y+=2)
    {
        sx = sy;

        for(x=0; x<w; ++x, sx.x+=2)
        {
            if(sa(sx) != edge_marker) continue;

            if(sa(sx, right) == edge_marker && sa(sx, left) == edge_marker) continue;
            if(sa(sx, bottom) == edge_marker && sa(sx, top) == edge_marker) continue;

            sa.set(background_marker, sx);
        }
    }
}

template <class SrcIterator, class SrcAccessor, class SrcValue>
inline
void beautifyCrackEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
           SrcValue edge_marker, SrcValue background_marker)
{
    beautifyCrackEdgeImage(src.first, src.second, src.third,
                    edge_marker, background_marker);
}


/** Helper class that stores edgel attributes.
*/
class Edgel
{
  public:
        /** The edgel's sub-pixel x coordinate.
        */    
    float x; 

        /** The edgel's sub-pixel y coordinate.
        */    
    float y;

        /** The edgel's strength (magnitude of the gradient vector).
        */    
    float strength;
    
        /**
        The edgel's orientation. This is the angle 
        between the x-axis and the edge, so that the bright side of the 
        edge is on the right. The angle is measured
        counter-clockwise in radians like this:


        \code

  edgel axis
      \  (bright side)
 (dark \
 side)  \ /__ 
         \\  \ orientation angle
          \  |
           +------------> x-axis
           |
           |
           |
           |
    y-axis V
        \endcode

        So, for example a vertical edge with its dark side on the left
        has orientation PI/2, and a horizontal edge with dark side on top
        has orientation 0. Obviously, the edge's orientation changes
        by PI if the contrast is reversed.

        */
    float orientation;
};

template <class Image>
void internalCannyFindEdgels(Image const & dx,
                             Image const & dy,
                             Image const & magnitude,
                             std::vector<Edgel> & edgels)
{
    typedef typename Image::PixelType PixelType;
    
    PixelType zero = NumericTraits<PixelType>::zero();
    double tan22_5 = M_SQRT2 - 1.0;
    
    for(int y=1; y<dx.height()-1; ++y)
    {
        for(int x=1; x<dx.width()-1; ++x)
        {
            bool maximum_found = false;
            Edgel edgel;
            
            PixelType gradx = dx(x,y);
            PixelType grady = dy(x,y);
            
            // find out quadrant
            if(abs(grady) < tan22_5*abs(gradx))
            {
                // north-south edge
                PixelType m1 = magnitude(x-1, y);
                PixelType m2 = magnitude(x, y);
                PixelType m3 = magnitude(x+1, y);
                
                if(m1 < m2 && m3 <= m2)
                {
                    edgel.y = y;
                
                    // local maximum => quadratic interpolation of sub-pixel location
                    PixelType del = (m1 - m3) / 2.0;
                    del /= (m1 + m3 - 2.0*m2);
                    edgel.x = x + del;
                    edgel.strength = m2;
                    
                    maximum_found = true;                    
                }
            }
            else if(abs(gradx) < tan22_5*abs(grady))
            {
                // west-east edge
                PixelType m1 = magnitude(x, y-1);
                PixelType m2 = magnitude(x, y);
                PixelType m3 = magnitude(x, y+1);
                
                if(m1 < m2 && m3 <= m2)
                {
                    edgel.x = x;
                
                    // local maximum => quadratic interpolation of sub-pixel location
                    PixelType del = (m1 - m3) / 2.0;
                    del /= (m1 + m3 - 2.0*m2);
                    edgel.y = y + del;
                    edgel.strength = m2;
                    
                    maximum_found = true;                    
                }
            }
            else if(gradx*grady < zero)
            {
                // north-west-south-east edge
                PixelType m1 = magnitude(x+1, y-1);
                PixelType m2 = magnitude(x, y);
                PixelType m3 = magnitude(x-1, y+1);
                
                if(m1 < m2 && m3 <= m2)
                {
                    // local maximum => quadratic interpolation of sub-pixel location
                    PixelType del = (m1 - m3) / 2.0;
                    del /= (m1 + m3 - 2.0*m2);
                    edgel.x = x - del;
                    edgel.y = y + del;
                    edgel.strength = m2;
                    
                    maximum_found = true;                    
                }
            }
            else
            {
                // north-east-south-west edge
                PixelType m1 = magnitude(x-1, y-1);
                PixelType m2 = magnitude(x, y);
                PixelType m3 = magnitude(x+1, y+1);
                
                if(m1 < m2 && m3 <= m2)
                {
                    // local maximum => quadratic interpolation of sub-pixel location
                    PixelType del = (m1 - m3) / 2.0;
                    del /= (m1 + m3 - 2.0*m2);
                    edgel.x = x + del;
                    edgel.y = y + del;
                    edgel.strength = m2;
                    
                    maximum_found = true;                    
                }
            }
            
            if(maximum_found)
            {
                double orientation = atan2(-grady, gradx) - M_PI * 1.5;
                if(orientation < 0.0)
                    orientation += 2.0*M_PI;
                edgel.orientation = orientation;
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

    This operator first calculates the gradient vector for each
    pixel of the image using first derivatives of a Gaussian at the 
    given scale. Then a very simple non-maxima supression is performed: 
    for each 3x3 neighborhood, it is determined whether the center pixel has 
    larger gradient magnitude than its two neighbors in gradient direction
    (where the direction is rounded into octands). If this is the case,
    a new \ref Edgel is appended to the given vector of <TT>edgels</TT>. The subpixel
    edgel position is determined by fitting a parabola 
    to the three gradient magnitude values 
    mentioned above. The sub-pixel location of the parabola's tip 
    and the gradient magnitude and direction are written in the newly created edgel.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                                std::vector<Edgel> & edgels, double scale);
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void 
        cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       std::vector<Edgel> & edgels, double scale);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
    \endcode
    
    SrcAccessor::value_type must be a type convertible to float
    
    <b> Preconditions:</b>
    
    \code
    scale > 0
    \endcode
*/
template <class SrcIterator, class SrcAccessor>
void cannyEdgelList(SrcIterator ul, SrcIterator lr, SrcAccessor src,
                        std::vector<Edgel> & edgels, double scale)
{
    int w = lr.x - ul.x;
    int h = lr.y - ul.y;
    
    // calculate image gradients
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
        
    BasicImage<TmpType> tmp(w,h), dx(w,h), dy(w,h);
    
    Kernel1D<double> smooth, grad;
    
    smooth.initGaussian(scale);
    grad.initGaussianDerivative(scale, 1);
    
    separableConvolveX(srcIterRange(ul, lr, src), destImage(tmp), kernel1d(grad));
    separableConvolveY(srcImageRange(tmp), destImage(dx), kernel1d(smooth));
    
    separableConvolveY(srcIterRange(ul, lr, src), destImage(tmp), kernel1d(grad));
    separableConvolveX(srcImageRange(tmp), destImage(dy), kernel1d(smooth));
    
    combineTwoImages(srcImageRange(dx), srcImage(dy), destImage(tmp),
                     MagnitudeFunctor<TmpType>());
    
    
    // find edgels
    internalCannyFindEdgels(dx, dy, tmp, edgels);
}

template <class SrcIterator, class SrcAccessor>
inline void 
cannyEdgelList(triple<SrcIterator, SrcIterator, SrcAccessor> src,
               std::vector<Edgel> & edgels, double scale)
{
    cannyEdgelList(src.first, src.second, src.third, edgels, scale);
}

/********************************************************/
/*                                                      */
/*                       cannyEdgeImage                 */
/*                                                      */
/********************************************************/

/** \brief Detect and mark edges in an edge image using Canny's algorithm.

    This operator first calls \ref cannyEdgelList() to generate an 
    edgel list for the given image. Than it scans this list and selects edgels
    whose strength is above the given <TT>gradient_threshold</TT>. For each of these 
    edgels, the edgel's location is rounded to the nearest pixel, and that
    pixel marked with the given <TT>edge_marker</TT>.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
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
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor, 
                  class GradValue, class DestValue>
        inline void cannyEdgeImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   double scale, GradValue gradient_threshold, DestValue edge_marker);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="edgedetection_8hxx-source.html">vigra/edgedetection.hxx</a>"<br>
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
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class GradValue, class DestValue>
void cannyEdgeImage(
           SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    std::vector<Edgel> edgels;
    
    cannyEdgelList(sul, slr, sa, edgels, scale);
    
    for(unsigned int i=0; i<edgels.size(); ++i)
    {
        if(gradient_threshold < edgels[i].strength)
        {
            Diff2D pix((int)(edgels[i].x + 0.5), (int)(edgels[i].y + 0.5));
            
            da.set(edge_marker, dul, pix);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class GradValue, class DestValue>
inline void cannyEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest,
           double scale, GradValue gradient_threshold, DestValue edge_marker)
{
    cannyEdgeImage(src.first, src.second, src.third,
                   dest.first, dest.second,
                   scale, gradient_threshold, edge_marker);
}

//@}

/** \page CrackEdgeImage Crack Edge Image

Crack edges are marked <i>between</i> the pixels of an image. 
A Crack Edge Image is an image that represents these edges. In order
to accomodate the cracks, the Crack Edge Image must be twice as large
as the original image (precisely (2*w - 1) by (2*h - 1)). A Crack Edge Image
can easily be derived from a binary image or from the signs of the 
response of a Laplacean filter. Consider the following sketch, where
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
