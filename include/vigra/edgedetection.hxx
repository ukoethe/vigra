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
 
 
#ifndef VIGRA_EDGEDETECTION_HXX
#define VIGRA_EDGEDETECTION_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/labelimage.hxx"

namespace vigra {

/** @heading Edge Detection
    @memo zero-crossing based edge detectors and related post-processing
*/
//@{ 
                                    
/********************************************************/
/*                                                      */
/*           differenceOfExponentialEdgeImage           */
/*                                                      */
/********************************************************/

/** Detect and mark edges in an image.
    This operator applies an exponential filter to the source image 
    at the given #scale# and subtracts the result from the original image. 
    Zero crossings are detected in the reulting difference image. Whenever the
    gradient at a zero crossing is greater than the given #gradient_threshold#,
    an edge point is marked (using #edge_marker#) in the destination image on
    the darker side of the zero crossing (note that zero crossings occur 
    {\it between} pixels). For example:
    
    \begin{verbatim}
    sign of difference image     resulting edge points (*)
    
        + - -                          * * .
        + + -               =>         . * *
        + + +                          . . .
    \end{verbatim}
    
    Non-edge pixels (#.#) remain untouched in the destination image. 
    The result can be improved by the post-processing operation \Ref{removeShortEdges}.
    A more accurate edge placement can be achieved with the function 
    \Ref{differenceOfExponentialCellGridImage}. 

    The source value type 
    (#SrcAccessor::value_type#) must be a linear algabra, i.e. addition, 
    subtraction and multiplication of the type with itself, and multiplication 
    with double and 
    \URL[NumericTraits]{templatestructNumericTraitsArithmeticType.html} must 
    be defined. In addition, this type must be less-comparable.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
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
    \end{verbatim}
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
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
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/edgedetection.hxx]{../include/vigra/edgedetection.hxx}\\
    Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage src(w,h), edges(w,h);
    
    // empty edge image
    edges = 0;
    ...
    
    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1 
    vigra::differenceOfExponentialEdgeImage(srcImageRange(src), destImage(edges), 
                                     0.8, 4.0, 1);
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
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
    \end{verbatim}
    
    {\bf Preconditions:}
    
    \begin{verbatim}
    scale > 0
    gradient_threshold > 0
    \end{verbatim}
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class SrcValue, class DestValue>
void differenceOfExponentialEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, SrcValue gradient_threshold, DestValue edge_marker)
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
          class SrcValue>
inline 
void differenceOfExponentialEdgeImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, SrcValue gradient_threshold)
{
    differenceOfExponentialEdgeImage(sul, slr, sa, dul, da, 
                                        scale, gradient_threshold, 1);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, 
      class SrcValue, class DestValue>
inline 
void differenceOfExponentialEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, SrcValue gradient_threshold,
       DestValue edge_marker)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold,
                    edge_marker);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class SrcValue>
inline 
void differenceOfExponentialEdgeImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, SrcValue gradient_threshold)
{
    differenceOfExponentialEdgeImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold, 1);
}

/********************************************************/
/*                                                      */
/*         differenceOfExponentialCellGridImage         */
/*                                                      */
/********************************************************/

/** Detect and mark edges in an cell grid image.
    This operator applies an exponential filter to the source image 
    at the given #scale# and subtracts the result from the original image. 
    Zero crossings are detected in the reulting difference image. Whenever the
    gradient at a zero crossing is greater than the given #gradient_threshold#,
    an edge point is marked (using #edge_marker#) in the destination image 
    {\it between} the corresponding original pixels. Topologically, this means we 
    must insert additional pixels between the original ones to represent the
    boundaries between the pixels (the so called zero- and one-cells, with the original
    pixels being two-cells). Within VIGRA, such an image is called \Ref{Cell Grid Image}.
    To allow insertion of the zero- and one-cells, the destination image must have twice the 
    size of the original (precisely, #(2*w-1)# by #(2*h-1)# pixels). Then the algorithm 
    proceeds as follows:
    
    \begin{verbatim}
sign of difference image     insert zero- and one-cells     resulting edge points (*)

                                     + . - . -                   . * . . .
      + - -                          . . . . .                   . * * * .
      + + -               =>         + . + . -           =>      . . . * .
      + + +                          . . . . .                   . . . * *
                                     + . + . +                   . . . . .
    \end{verbatim}
    
    Thus the edge points are marked where they actually are - in between the pixels. 
    An important property of the resulting edge image is that it conforms to the notion 
    of well-composedness as defined by Latecki et al., i.e. connected regions and edges 
    obtained by a subsequent \Ref{Connected Components Labeling} do not depend on 
    whether 4- or 8-connectivity is used.
    The non-edge pixels (#.#) in the destination image remain unchanged. 
    The result conformes to the requirements of a \Ref{Cell Grid Image}. It can be further
    improved by the post-processing operations \Ref{removeShortEdges} and
    \Ref{closeGapsInCellGridImage}.
    
    The source value type (#SrcAccessor::value_type#) must be a linear algabra, i.e. addition, 
    subtraction and multiplication of the type with itself, and multiplication 
    with double and 
    \URL[NumericTraits]{templatestructNumericTraitsArithmeticType.html} must 
    be defined. In addition, this type must be less-comparable.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor,
              class GradValue,
              class DestValue = DestAccessor::value_type>
        void differenceOfExponentialCellGridImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
               DestIterator dul, DestAccessor da,
               double scale, GradValue gradient_threshold, 
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \end{verbatim}
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor, 
              class GradValue,
              class DestValue = DestAccessor::value_type>
        inline 
        void differenceOfExponentialCellGridImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               pair<DestIterator, DestAccessor> dest,
               double scale, GradValue gradient_threshold,
               DestValue edge_marker = NumericTraits<DestValue>::one())
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/edgedetection.hxx]{../include/vigra/edgedetection.hxx}\\
    Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage src(w,h), edges(2*w-1,2*h-1);
    
    // empty edge image
    edges = 0;
    ...
    
    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1 
    vigra::differenceOfExponentialCellGridImage(srcImageRange(src), destImage(edges), 
                                     0.8, 4.0, 1);
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
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
    \end{verbatim}
    
    {\bf Preconditions:}
    
    \begin{verbatim}
    scale > 0
    gradient_threshold > 0
    \end{verbatim}
    
    The destination image must have twice the size of the source:
    \begin{verbatim}
    w_dest = 2 * w_src - 1
    h_dest = 2 * h_src - 1
    \end{verbatim}
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor, class DestValue>
void differenceOfExponentialCellGridImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa,
           DestIterator dul, DestAccessor da,
           double scale, double gradient_threshold, 
           DestValue edge_marker)
{
    vigra_precondition(scale > 0,
                 "differenceOfExponentialCellGridImage(): scale > 0 required.");
         
    vigra_precondition(gradient_threshold > 0,
                 "differenceOfExponentialCellGridImage(): "
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
          class DestIterator, class DestAccessor, class DestValue>
inline 
void differenceOfExponentialCellGridImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       double scale, double gradient_threshold,
       DestValue edge_marker)
{
    differenceOfExponentialCellGridImage(src.first, src.second, src.third,
                                        dest.first, dest.second,
                    scale, gradient_threshold,
                    edge_marker);
}

/********************************************************/
/*                                                      */
/*                  removeShortEdges                    */
/*                                                      */
/********************************************************/

/** Remove short edges from an edge image.
    This algorithm is applied as a post-processing operation of 
    \Ref{differenceOfExponentialEdgeImage} and \Ref{differenceOfExponentialCellGridImage}. 
    It removes all edges that are shorter than #min_edge_length#. The corresponding
    pixels are set to the #non_edge_marker#. The idea behind this algorithms is
    that very short edges are probably caused by noise and don't represent interesting
    image structure. Technically, the algorithms executes a connected components labeling,
    so the image's value type must be equality comparable. 
    
    If the source image fulfilled the requirements of a \Ref{Cell Grid Image},
    it will still do so after application of this algorithm.
    
    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place, 
    i.e. on only one image. Also, the algorithm assumes that all non-edges pixels are already
    marked with the given #non_edge_marker# value.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        void removeShortEdges(
               Iterator sul, Iterator slr, Accessor sa,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \end{verbatim}
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class Iterator, class Accessor, class SrcValue>
        inline 
        void removeShortEdges(
               triple<Iterator, Iterator, Accessor> src,
               int min_edge_length, SrcValue non_edge_marker)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/edgedetection.hxx]{../include/vigra/edgedetection.hxx}\\
    Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage src(w,h), edges(w,h);
    
    // empty edge image
    edges = 0;
    ...
    
    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1 
    vigra::differenceOfExponentialEdgeImage(srcImageRange(src), destImage(edges), 
                                     0.8, 4.0, 1);
                    
    // zero edges shorter than 10 pixels
    vigra::removeShortEdges(srcImageRange(edges), 10, 0);
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u
    
    SrcValue non_edge_marker;
    src_accessor.set(non_edge_marker, src_upperleft);
    \end{verbatim}
*/
template <class Iterator, class Accessor, class Value>
void removeShortEdges(
               Iterator sul, Iterator slr, Accessor sa,
           int min_edge_length, Value non_edge_marker)
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
       int min_edge_length, Value non_edge_marker)
{
    removeShortEdges(src.first, src.second, src.third,
                     min_edge_length, non_edge_marker);
}

/********************************************************/
/*                                                      */
/*             closeGapsInCellGridImage                 */
/*                                                      */
/********************************************************/

/** Close one-pixel wide gaps in a cell grid edge image.
    This algorithm is typically applied as a post-processing operation of 
    \Ref{differenceOfExponentialCellGridImage}. The source image must fulfill
    the requirements of a \Ref{Cell Grid Image}, and will still do so after 
    application of this algorithm.

    It closes one pixel wide gaps in the edges resulting from this algorithm. 
    Since these gaps are usually caused by zero crossing slightly below the gradient 
    threshold used in edge detection, this algorithms acts like a weak hysteresis 
    thresholding. The newly found edge pixels are marked with the givem #edge_marker#.
    The image's value type must be equality comparable. 
    
    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place, 
    i.e. on only one image.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void closeGapsInCellGridImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
               SrcValue edge_marker)
    }
    \end{verbatim}
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        inline
        void closeGapsInCellGridImage(
               triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/edgedetection.hxx]{../include/vigra/edgedetection.hxx}\\
    Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage src(w,h), edges(2*w-1, 2*h-1);
    
    // empty edge image
    edges = 0;
    ...
    
    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1 
    vigra::differenceOfExponentialCellGridImage(srcImageRange(src), destImage(edges), 
                                         0.8, 4.0, 1);
                    
    // close gaps, mark with 1
    vigra::closeGapsInCellGridImage(srcImageRange(edges), 1);
                    
    // zero edges shorter than 20 pixels
    vigra::removeShortEdges(srcImageRange(edges), 10, 0);
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator src_upperleft, src_lowerright;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u
    u != u
    
    SrcValue edge_marker;
    src_accessor.set(edge_marker, src_upperleft);
    \end{verbatim}
*/
template <class SrcIterator, class SrcAccessor, class SrcValue>
void closeGapsInCellGridImage(
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
void closeGapsInCellGridImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
       SrcValue edge_marker)
{
    closeGapsInCellGridImage(src.first, src.second, src.third,
                    edge_marker);
}

/********************************************************/
/*                                                      */
/*               beautifyCellGridImage                  */
/*                                                      */
/********************************************************/

/** Beautify Cell Grid Edge Image for visualization.
    This algorithm is applied as a post-processing operation of 
    \Ref{differenceOfExponentialCellGridImage}. The source image must fulfill
    the requirements of a \Ref{Cell Grid Image}, but will {\bf not} do so after 
    application of this algorithm. In particular, the algorithm removes zero-cells 
    marked as edges to avoid staircase effects on diagonal lines like this:
    
    \begin{verbatim}
    original edge points (*)     resulting edge points
    
          . * . . .                   . * . . .
          . * * * .                   . . * . .
          . . . * .           =>      . . . * .
          . . . * *                   . . . . *
          . . . . .                   . . . . .
    \end{verbatim}
    
    Therfore, this algorithm should only be applied as a vizualization aid, i.e. 
    for human inspection. The algorithm assumes that edges are marked with #edge_marker#, 
    and background pixels with #background_marker#. The image's value type must be 
    equality comparable. 
    
    Note that this algorithm, unlike most other algorithms in VIGRA, operates in-place, 
    i.e. on only one image.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        void beautifyCellGridImage(
               SrcIterator sul, SrcIterator slr, SrcAccessor sa, 
               SrcValue edge_marker, SrcValue background_marker)
    }
    \end{verbatim}
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, class SrcValue>
        inline
        void beautifyCellGridImage(
                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
               SrcValue edge_marker, SrcValue background_marker)
    }
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/edgedetection.hxx]{../include/vigra/edgedetection.hxx}\\
    Namespace: vigra
    
    \begin{verbatim}
    vigra::BImage src(w,h), edges(2*w-1, 2*h-1);
    
    // empty edge image
    edges = 0;
    ...
    
    // find edges at scale 0.8 with gradient larger than 4.0, mark with 1 
    vigra::differenceOfExponentialCellGridImage(srcImageRange(src), destImage(edges), 
                                         0.8, 4.0, 1);
                    
    // beautify edge image for visualization
    vigra::beautifyCellGridImage(srcImageRange(edges), 1, 0);
    
    // show to the user
    window.open(edges);
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    SrcImageIterator src_upperleft, src_lowerright;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);

    u == u
    u != u
    
    SrcValue background_marker;
    src_accessor.set(background_marker, src_upperleft);
    \end{verbatim}
*/
template <class SrcIterator, class SrcAccessor, class SrcValue>
void beautifyCellGridImage(
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
void beautifyCellGridImage(
           triple<SrcIterator, SrcIterator, SrcAccessor> src,
           SrcValue edge_marker, SrcValue background_marker)
{
    beautifyCellGridImage(src.first, src.second, src.third,
                    edge_marker, background_marker);
}

//@}

} // namespace vigra

#endif // VIGRA_EDGEDETECTION_HXX
