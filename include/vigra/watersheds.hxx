/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2005 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#ifndef VIGRA_WATERSHEDS_HXX
#define VIGRA_WATERSHEDS_HXX

#include <functional>
#include "mathutil.hxx"
#include "stdimage.hxx"
#include "pixelneighborhood.hxx"

namespace vigra {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int watershedLabeling(SrcIterator upperlefts,
                        SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da,
                        Neighborhood neighborhood)
{
    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y,i;

    SrcIterator ys(upperlefts);
    SrcIterator xs(ys);

    // temporary image to store region labels
    typedef IImage LabelImage;
    typedef LabelImage::traverser LabelTraverser;
    
    LabelImage labelimage(w, h);

    LabelTraverser yt = labelimage.upperLeft();
    LabelTraverser xt(yt);

    // Kovalevsky's clever idea to use
    // image iterator and scan order iterator simultaneously
    LabelImage::ScanOrderIterator label = labelimage.begin();
    
    // initialize the neighborhood circulators
    NeighborOffsetCirculator<Neighborhood> ncstart(Neighborhood::CausalFirst);
    NeighborOffsetCirculator<Neighborhood> ncstartBorder(Neighborhood::North);
    NeighborOffsetCirculator<Neighborhood> ncend(Neighborhood::CausalLast);
    ++ncend;
    NeighborOffsetCirculator<Neighborhood> ncendBorder(Neighborhood::North);
    ++ncendBorder;

    // pass 1: scan image from upper left to lower right
    // to find connected components

    // Each component will be represented by a tree of pixels. Each
    // pixel contains the scan order address of its parent in the
    // tree.  In order for pass 2 to work correctly, the parent must
    // always have a smaller scan order address than the child.
    // Therefore, we can merge trees only at their roots, because the
    // root of the combined tree must have the smallest scan order
    // address among all the tree's pixels/ nodes.  The root of each
    // tree is distinguished by pointing to itself (it contains its
    // own scan order address). This condition is enforced whenever a
    // new region is found or two regions are merged
    xs = ys;
    xt = yt;
    
    *xt = 0;
    
    ++xs.x;
    ++xt.x;
    for(x = 1; x != w; ++x, ++xs.x, ++xt.x)
    {
        if((*xs & Neighborhood::directionBit(Neighborhood::West)) ||
           (xs[Neighborhood::west()] & Neighborhood::directionBit(Neighborhood::East)))
        {
            *xt = xt[Neighborhood::west()];
        }
        else
        {
            *xt = x;
        }
    }
    
    ++ys.y;
    ++yt.y;
    for(y = 1; y != h; ++y, ++ys.y, ++yt.y)
    {
        xs = ys;
        xt = yt;

        for(x = 0; x != w; ++x, ++xs.x, ++xt.x)
        {
            NeighborOffsetCirculator<Neighborhood> nc(x == w-1
                                                        ? ncstartBorder
                                                        : ncstart);
            NeighborOffsetCirculator<Neighborhood> nce(x == 0 
                                                         ? ncendBorder 
                                                         : ncend);
            *xt = x + w*y; // default: new region            
            for(; nc != nce; ++nc)
            {
                if((*xs & nc.directionBit()) || (xs[*nc] & nc.oppositeDirectionBit()))
                {
                    int neighborLabel = xt[*nc];
                    // find the root label of a label tree
                    while(neighborLabel != label[neighborLabel])
                    {
                        neighborLabel = label[neighborLabel];
                    }
                    if(neighborLabel < *xt) // always keep the smallest among the possible labels
                    {
                        label[*xt] = neighborLabel;
                        *xt = neighborLabel;
                    }
                    else
                    {
                        label[neighborLabel] = *xt;
                    }
                }
            }
        }
    }

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    DestIterator yd(upperleftd);

    unsigned int count = 0;
    i = 0;
    for(y=0; y != h; ++y, ++yd.y)
    {
        DestIterator xd(yd);
        for(x = 0; x != w; ++x, ++xd.x, ++i)
        {
            if(label[i] == i)
            {
                label[i] = ++count;
            }
            else
            {
                label[i] = label[label[i]];
            }
            da.set(label[i], xd);
        }
    }
    return count;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void prepareWatersheds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                      DestIterator upperleftd, DestAccessor da, 
                      FourNeighborCode)
{
    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y;

    SrcIterator ys(upperlefts);
    SrcIterator xs(ys);

    DestIterator yd = upperleftd;

    for(y = 0; y != h; ++y, ++ys.y, ++yd.y)
    {
        xs = ys;
        DestIterator xd = yd;

        for(x = 0; x != w; ++x, ++xs.x, ++xd.x)
        {
            AtImageBorder atBorder = isAtImageBorder(x,y,w,h);
            typename SrcAccessor::value_type v = sa(xs);
            // the following choice causes minima to point
            // to their lowest neighbor -- would this be better???
            // typename SrcAccessor::value_type v = NumericTraits<typename SrcAccessor::value_type>::max();
            int o = 0; // means center is minimum
            if(atBorder == NotAtBorder)
            {
                NeighborhoodCirculator<SrcIterator, FourNeighborCode>  c(xs), cend(c);
                do {
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
                while(++c != cend);
            }
            else
            {
                RestrictedNeighborhoodCirculator<SrcIterator, FourNeighborCode>  c(xs, atBorder), cend(c);
                do {
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
                while(++c != cend);
            }
            da.set(o, xd);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void prepareWatersheds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                      DestIterator upperleftd, DestAccessor da, 
                      EightNeighborCode)
{
    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y;

    SrcIterator ys(upperlefts);
    SrcIterator xs(ys);

    DestIterator yd = upperleftd;

    for(y = 0; y != h; ++y, ++ys.y, ++yd.y)
    {
        xs = ys;
        DestIterator xd = yd;

        for(x = 0; x != w; ++x, ++xs.x, ++xd.x)
        {
            AtImageBorder atBorder = isAtImageBorder(x,y,w,h);
            typename SrcAccessor::value_type v = sa(xs);
            // the following choice causes minima to point
            // to their lowest neighbor -- would this be better???
            // typename SrcAccessor::value_type v = NumericTraits<typename SrcAccessor::value_type>::max();
            int o = 0; // means center is minimum
            if(atBorder == NotAtBorder)
            {
                // handle diagonal and principal neighbors separately
                // so that principal neighbors are preferred when there are 
                // candidates with equal strength
                NeighborhoodCirculator<SrcIterator, EightNeighborCode>  
                                      c(xs, EightNeighborCode::NorthEast);
                for(int i = 0; i < 4; ++i, c += 2)
                {
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
                --c;
                for(int i = 0; i < 4; ++i, c += 2)
                {
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
            }
            else
            {
                RestrictedNeighborhoodCirculator<SrcIterator, EightNeighborCode>  
                             c(xs, atBorder), cend(c);
                do 
                {
                    if(!c.isDiagonal())
                        continue;
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
                while(++c != cend);
                do 
                {
                    if(c.isDiagonal())
                        continue;
                    if(sa(c) <= v)
                    {
                        v = sa(c);
                        o = c.directionBit();
                    }
                }
                while(++c != cend);
            }
            da.set(o, xd);
        }
    }
}

/** \addtogroup SeededRegionGrowing Region Segmentation Algorithms
    Region growing, watersheds, and voronoi tesselation
*/
//@{

/********************************************************/
/*                                                      */
/*                      watersheds                      */
/*                                                      */
/********************************************************/

/** \brief Region Segmentation by means of the watershed algorithm.

    This function implements the union-find version of the watershed algorithms
    as described in

    J. Roerdink, R. Meijster: "<em>The watershed transform: definitions, algorithms,
    and parallelization stretegies</em>", Fundamenta Informaticae, 41:187-228, 2000

    The source image is a boundary indicator such as the gradient magnitude
    of the trace of the \ref boundaryTensor(). Local minima of the boundary indicator
    are used as region seeds, and all other pixels are recursively assigned to the same 
    region as their lowest neighbor. Pass \ref vigra::EightNeighborCode or 
    \ref vigra::FourNeighborCode to determine the neighborhood where pixel values 
    are compared. The pixel type of the input image must be <tt>LessThanComparable</tt>.
    The function uses accessors. 
    
    Note that VIGRA provides an alternative implementaion of the watershed transform via
    \ref seededRegionGrowing(). It is slower, but handles plateaus better 
    and allows to keep a one pixel wide boundary between regions.
    
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int 
        watersheds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                   DestIterator upperleftd, DestAccessor da, 
                   Neighborhood neighborhood = EightNeighborCode())
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int 
        watersheds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest, 
                   Neighborhood neighborhood = EightNeighborCode())
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<<a href="watersheds_8hxx-source.html">vigra/watersheds.hxx</a>\><br>
    Namespace: vigra

    Example: watersheds of the gradient magnitude.

    \code
    vigra::BImage in(w,h);
    ... // read input data
    
    vigra::FImage gradx(x,y), grady(x,y), gradMag(x,y);
    gaussianGradient(srcImageRange(src), destImage(gradx), destImage(grady), 3.0);
    combineTwoImages(srcImageRange(gradx), srcImage(grady), destImage(gradMag),
                     vigra::MagnitudeFunctor<float>());
    
    // the pixel type of the destination image must be large enough to hold
    // numbers up to 'max_region_label' to prevent overflow
    vigra::IImage labeling(x,y);
    int max_region_label = watersheds(srcImageRange(gradMag), destImage(labeling));

    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    // compare src values
    src_accessor(src_upperleft) <= src_accessor(src_upperleft)

    // set result
    int label;
    dest_accessor.set(label, dest_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> unsigned int watersheds)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int 
watersheds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
           DestIterator upperleftd, DestAccessor da, Neighborhood neighborhood)
{
    SImage orientationImage(lowerrights - upperlefts);
    SImage::traverser yo = orientationImage.upperLeft();

    prepareWatersheds(upperlefts, lowerrights, sa, 
                     orientationImage.upperLeft(), orientationImage.accessor(), neighborhood);
    return watershedLabeling(orientationImage.upperLeft(), orientationImage.lowerRight(), orientationImage.accessor(),
                             upperleftd, da, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int 
watersheds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
           DestIterator upperleftd, DestAccessor da)
{
    return watersheds(upperlefts, lowerrights, sa, upperleftd, da, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline unsigned int 
watersheds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest, Neighborhood neighborhood)
{
    return watersheds(src.first, src.second, src.third, dest.first, dest.second, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int 
watersheds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest)
{
    return watersheds(src.first, src.second, src.third, dest.first, dest.second);
}

//@}

} // namespace vigra

#endif // VIGRA_WATERSHEDS_HXX
