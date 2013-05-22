/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2005 by Ullrich Koethe                  */
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

#ifndef VIGRA_WATERSHEDS_HXX
#define VIGRA_WATERSHEDS_HXX

#include <functional>
#include "mathutil.hxx"
#include "stdimage.hxx"
#include "pixelneighborhood.hxx"
#include "localminmax.hxx"
#include "labelimage.hxx"
#include "seededregiongrowing.hxx"
#include "functorexpression.hxx"
#include "union_find.hxx"

namespace vigra {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int watershedLabeling(SrcIterator upperlefts,
                        SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da,
                        Neighborhood)
{
    typedef typename DestAccessor::value_type LabelType;

    int w = lowerrights.x - upperlefts.x;
    int h = lowerrights.y - upperlefts.y;
    int x,y;

    SrcIterator ys(upperlefts);
    SrcIterator xs(ys);
    DestIterator yd(upperleftd);
    DestIterator xd(yd);

    // temporary image to store region labels
    detail::UnionFindArray<LabelType> labels;

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
    da.set(labels.finalizeLabel(labels.nextFreeLabel()), xd);

    ++xs.x;
    ++xd.x;
    for(x = 1; x != w; ++x, ++xs.x, ++xd.x)
    {
        if((sa(xs) & Neighborhood::directionBit(Neighborhood::West)) ||
           (sa(xs, Neighborhood::west()) & Neighborhood::directionBit(Neighborhood::East)))
        {
            da.set(da(xd, Neighborhood::west()), xd);
        }
        else
        {
            da.set(labels.finalizeLabel(labels.nextFreeLabel()), xd);
        }
    }

    ++ys.y;
    ++yd.y;
    for(y = 1; y != h; ++y, ++ys.y, ++yd.y)
    {
        xs = ys;
        xd = yd;

        for(x = 0; x != w; ++x, ++xs.x, ++xd.x)
        {
            NeighborOffsetCirculator<Neighborhood> nc(x == w-1
                                                        ? ncstartBorder
                                                        : ncstart);
            NeighborOffsetCirculator<Neighborhood> nce(x == 0
                                                         ? ncendBorder
                                                         : ncend);
            LabelType currentLabel = labels.nextFreeLabel();
            for(; nc != nce; ++nc)
            {
                if((sa(xs) & nc.directionBit()) || (sa(xs, *nc) & nc.oppositeDirectionBit()))
                {
                    currentLabel = labels.makeUnion(da(xd,*nc), currentLabel);
                }
            }
            da.set(labels.finalizeLabel(currentLabel), xd);
        }
    }

    unsigned int count = labels.makeContiguous();

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    yd = upperleftd;
    for(y=0; y != h; ++y, ++yd.y)
    {
        DestIterator xd(yd);
        for(x = 0; x != w; ++x, ++xd.x)
        {
            da.set(labels[da(xd)], xd);
        }
    }
    return count;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int watershedLabeling(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               pair<DestIterator, DestAccessor> dest,
                               Neighborhood neighborhood)
{
    return watershedLabeling(src.first, src.second, src.third,
                             dest.first, dest.second, neighborhood);
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

    /**\brief Options object for generateWatershedSeeds().
     *
        <b> Usage:</b>

        <b>\#include</b> \<vigra/watersheds.hxx\><br>
        Namespace: vigra
        
        \code
        IImage seeds(boundary_indicator.size());
        
        // detect all minima in 'boundary_indicator' that are below gray level 22
        generateWatershedSeeds(srcImageRange(boundary_indicator),
                               destImage(seeds),
                               SeedOptions().minima().threshold(22.0));
        \endcode
     */
class SeedOptions
{
public:
    enum DetectMinima { LevelSets, Minima, ExtendedMinima, Unspecified };
    
    double thresh;
    DetectMinima mini;
    
        /**\brief Construct default options object.
         *
            Defaults are: detect minima without thresholding (i.e. all minima).
         */
    SeedOptions()
    : thresh(NumericTraits<double>::max()),
      mini(Minima)
    {}
    
        /** Generate seeds at minima.
        
            Default: true
         */
    SeedOptions & minima()
    {
        mini = Minima;
        return *this;
    }
    
        /** Generate seeds at minima and minimal plateaus.
        
            Default: false
         */
    SeedOptions & extendedMinima()
    {
        mini = ExtendedMinima;
        return *this;
    }
    
        /** Generate seeds as level sets.
        
            Note that you must also set a threshold to define which level set is to be used.<br>
            Default: false
         */
    SeedOptions & levelSets()
    {
        mini = LevelSets;
        return *this;
    }
    
        /** Generate seeds as level sets at given threshold.
        
            Equivalent to <tt>SeedOptions().levelSet().threshold(threshold)</tt><br>
            Default: false
         */
    SeedOptions & levelSets(double threshold)
    {
        mini = LevelSets;
        thresh = threshold;
        return *this;
    }
    
        /** Set threshold.
        
            The threshold will be used by both the minima and level set variants
            of seed generation.<br>
            Default: no thresholding
         */
    SeedOptions & threshold(double threshold)
    {
        thresh = threshold;
        return *this;
    }
    
        // check whether the threshold has been set for the target type T
    template <class T>
    bool thresholdIsValid() const
    {
        return thresh < double(NumericTraits<T>::max());
    }
    
        // indicate that this option object is invalid (for internal use in watersheds)
    SeedOptions & unspecified()
    {
        mini = Unspecified;
        return *this;
    }
};

/** \brief Generate seeds for watershed computation and seeded region growing.

    The source image is a boundary indicator such as the gradient magnitude
    or the trace of the \ref boundaryTensor(). Seeds are generally generated
    at locations where the boundaryness (i.e. the likelihood of the point being on the
    boundary) is very small. In particular, seeds can be placed by either
    looking for local minima (possibly including minimal plateaus) of the boundaryness,
    of by looking at level sets (i.e. regions where the boundaryness is below a threshold).
    Both methods can also be combined, so that only minima below a threshold are returned.
    The particular seeding strategy is specified by the <tt>options</tt> object 
    (see \ref SeedOptions).
    
    The pixel type of the input image must be <tt>LessThanComparable</tt>.
    The pixel type of the output image must be large enough to hold the labels for all seeds.
    (typically, you will use <tt>UInt32</tt>). The function will label seeds by consecutive integers
    (starting from 1) and returns the largest label it used.
    
    Pass \ref vigra::EightNeighborCode or \ref vigra::FourNeighborCode to determine the 
    neighborhood where pixel values are compared. 
    
    The function uses accessors.

    <b> Declarations:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        template <unsigned int N, class T, class S1,
                                  class Label, class S2>
        Label
        generateWatershedSeeds(MultiArrayView<N, T, S1> const & data,
                               MultiArrayView<N, Label, S2> seeds,
                               NeighborhoodType neighborhood = IndirectNeighborhood,
                               SeedOptions const & options = SeedOptions());
    }
    \endcode

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int
        generateWatershedSeeds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                               DestIterator upperleftd, DestAccessor da, 
                               Neighborhood neighborhood = EightNeighborCode(),
                               SeedOptions const & options = SeedOptions());
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int
        generateWatershedSeeds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               pair<DestIterator, DestAccessor> dest, 
                               Neighborhood neighborhood = EightNeighborCode(),
                               SeedOptions const & options = SeedOptions());
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/_multi_watersheds.hxx\> (first form)<br>
    <b>\#include</b> \<vigra/watersheds.hxx\> (second and third form)<br>
    Namespace: vigra

    For detailed examples see watershedsRegionGrowing() and watershedsMultiArray().
*/
doxygen_overloaded_function(template <...> unsigned int generateWatershedSeeds)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int
generateWatershedSeeds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                       DestIterator upperleftd, DestAccessor da, 
                       Neighborhood neighborhood,
                       SeedOptions const & options = SeedOptions())
{
    using namespace functor;
    typedef typename SrcAccessor::value_type SrcType;
    
    vigra_precondition(options.mini != SeedOptions::LevelSets || 
                       options.thresholdIsValid<SrcType>(),
        "generateWatershedSeeds(): SeedOptions.levelSets() must be specified with threshold.");
    
    Diff2D shape = lowerrights - upperlefts;
    BImage seeds(shape);
    
    if(options.mini == SeedOptions::LevelSets)
    {
        transformImage(srcIterRange(upperlefts, lowerrights, sa),
                       destImage(seeds),
                       ifThenElse(Arg1() <= Param(options.thresh), Param(1), Param(0)));
    }
    else
    {
        LocalMinmaxOptions lm_options;
        lm_options.neighborhood(Neighborhood::DirectionCount)
                  .markWith(1.0)
                  .allowAtBorder()
                  .allowPlateaus(options.mini == SeedOptions::ExtendedMinima);
        if(options.thresholdIsValid<SrcType>())
            lm_options.threshold(options.thresh);
            
        localMinima(srcIterRange(upperlefts, lowerrights, sa), destImage(seeds),
                    lm_options);
    }
    
    return labelImageWithBackground(srcImageRange(seeds), destIter(upperleftd, da), 
                                    Neighborhood::DirectionCount == 8, 0);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
generateWatershedSeeds(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                       DestIterator upperleftd, DestAccessor da, 
                       SeedOptions const & options = SeedOptions())
{
    return generateWatershedSeeds(upperlefts, lowerrights, sa, upperleftd, da, 
                                   EightNeighborCode(), options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline unsigned int
generateWatershedSeeds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       pair<DestIterator, DestAccessor> dest, 
                       Neighborhood neighborhood,
                       SeedOptions const & options = SeedOptions())
{
    return generateWatershedSeeds(src.first, src.second, src.third,
                                   dest.first, dest.second,    
                                   neighborhood, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
generateWatershedSeeds(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                       pair<DestIterator, DestAccessor> dest, 
                       SeedOptions const & options = SeedOptions())
{
    return generateWatershedSeeds(src.first, src.second, src.third,
                                   dest.first, dest.second,    
                                   EightNeighborCode(), options);
}

/********************************************************/
/*                                                      */
/*                      watersheds                      */
/*                                                      */
/********************************************************/

/** \brief Region segmentation by means of the union-find watershed algorithm.

    This function implements the union-find version of the watershed algorithms
    described as algorithm 4.7 in

    J. Roerdink, R. Meijster: "<em>The watershed transform: definitions, algorithms,
    and parallelization strategies</em>", Fundamenta Informaticae, 41:187-228, 2000

    The source image is a boundary indicator such as the gaussianGradientMagnitude()
    or the trace of the \ref boundaryTensor(). Local minima of the boundary indicator
    are used as region seeds, and all other pixels are recursively assigned to the same
    region as their lowest neighbor. Pass \ref vigra::EightNeighborCode or
    \ref vigra::FourNeighborCode to determine the neighborhood where pixel values
    are compared. The pixel type of the input image must be <tt>LessThanComparable</tt>.
    The function uses accessors.

    Note that VIGRA provides an alternative implementation of the watershed transform via
    \ref watershedsRegionGrowing(). It is slower, but offers many more configuration options.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int
        watershedsUnionFind(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
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
        watershedsUnionFind(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            Neighborhood neighborhood = EightNeighborCode())
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/watersheds.hxx\><br>
    Namespace: vigra

    Example: watersheds of the gradient magnitude.

    \code
    vigra::BImage in(w,h);
    ... // read input data

    // compute gradient magnitude as boundary indicator
    vigra::FImage gradMag(w, h);
    gaussianGradientMagnitude(srcImageRange(src), destImage(gradMag), 3.0);

    // the pixel type of the destination image must be large enough to hold
    // numbers up to 'max_region_label' to prevent overflow
    vigra::IImage labeling(w,h);
    int max_region_label = watershedsUnionFind(srcImageRange(gradMag), destImage(labeling));

    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_upperleft, src_lowerright;
    DestIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    // compare src values
    src_accessor(src_upperleft) <= src_accessor(src_upperleft)

    // set result
    int label;
    dest_accessor.set(label, dest_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> unsigned int watershedsUnionFind)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int
watershedsUnionFind(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                    DestIterator upperleftd, DestAccessor da, 
                    Neighborhood neighborhood)
{
    SImage orientationImage(lowerrights - upperlefts);

    prepareWatersheds(upperlefts, lowerrights, sa,
                     orientationImage.upperLeft(), orientationImage.accessor(), neighborhood);
    return watershedLabeling(orientationImage.upperLeft(), orientationImage.lowerRight(), orientationImage.accessor(),
                             upperleftd, da, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
watershedsUnionFind(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
           DestIterator upperleftd, DestAccessor da)
{
    return watershedsUnionFind(upperlefts, lowerrights, sa, upperleftd, da, EightNeighborCode());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline unsigned int
watershedsUnionFind(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest, Neighborhood neighborhood)
{
    return watershedsUnionFind(src.first, src.second, src.third, 
                                dest.first, dest.second, neighborhood);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
watershedsUnionFind(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest)
{
    return watershedsUnionFind(src.first, src.second, src.third, 
                                dest.first, dest.second);
}

/** \brief Options object for watershedsRegionGrowing().

    <b> Usage:</b>

    see watershedsRegionGrowing() for detailed examples.
*/
class WatershedOptions
{
  public:
    enum Method { RegionGrowing, UnionFind };
  
    double max_cost, bias;
    SRGType terminate;
    Method method;
    unsigned int biased_label, bucket_count;
    SeedOptions seed_options;
    
    
    
        /** \brief Create options object with default settings.

            Defaults are: perform complete grow (all pixels are assigned to regions),
            use standard algorithm, assume that the destination image already contains 
            region seeds.
        */
    WatershedOptions()
    : max_cost(0.0),
      bias(1.0),
      terminate(CompleteGrow),
      method(RegionGrowing),
      biased_label(0),
      bucket_count(0),
      seed_options(SeedOptions().unspecified())
    {}    
    
        /** \brief Perform complete grow.

            That is, all pixels are assigned to regions, without explicit contours
            in between.
            
            Default: true
        */
    WatershedOptions & completeGrow()
    {
        terminate = SRGType(CompleteGrow | (terminate & StopAtThreshold));
        return *this;
    }
    
        /** \brief Keep one-pixel wide contour between regions.
        
            Note that this option is unsupported by the turbo algorithm.

            Default: false
        */
    WatershedOptions & keepContours()
    {
        terminate = SRGType(KeepContours | (terminate & StopAtThreshold));
        return *this;
    }
    
        /** \brief Set \ref SRGType explicitly.
        
            Default: CompleteGrow
        */
    WatershedOptions & srgType(SRGType type)
    {
        terminate = type;
        return *this;
    }
    
        /** \brief Stop region growing when the boundaryness exceeds the threshold.
        
            This option may be combined with completeGrow() and keepContours().
        
            Default: no early stopping
        */
    WatershedOptions & stopAtThreshold(double threshold)
    {
        terminate = SRGType(terminate | StopAtThreshold);
        max_cost = threshold;
        return *this;
    }
    
        /** \brief Use a simpler, but faster region growing algorithm.
        
            The algorithm internally uses a \ref BucketQueue to determine
            the processing order of the pixels. This is only useful,
            when the input boundary indicator image contains integers
            in the range <tt>[0, ..., bucket_count-1]</tt>. Since
            these boundary indicators are typically represented as
            UInt8 images, the default <tt>bucket_count</tt> is 256.
        
            Default: don't use the turbo algorithm
        */
    WatershedOptions & turboAlgorithm(unsigned int bucket_count = 256)
    {
        this->bucket_count = bucket_count;
        method = RegionGrowing;
        return *this;
    }
    
        /** \brief Specify seed options.
        
            In this case, watershedsRegionGrowing() assumes that the destination
            image does not yet contain seeds. It will therefore call 
            generateWatershedSeeds() and pass on the seed options.
        
            Default: don't compute seeds (i.e. assume that destination image already
            contains seeds).
        */
    WatershedOptions & seedOptions(SeedOptions const & s)
    {
        seed_options = s;
        return *this;
    }
    
        /** \brief Bias the cost of the specified region by the given factor.
        
            In certain applications, one region (typically the background) should
            be preferred in region growing. This is most easily achieved
            by adjusting the assignment cost for that region as <tt>factor*cost</tt>,
            with a factor slightly below 1.
        
            Default: don't bias any region.
        */
    WatershedOptions & biasLabel(unsigned int label, double factor)
    {
        biased_label = label;
        bias = factor;
        return *this;
    }
    
        /** \brief Specify the algorithm to be used.
        
            Possible values are <tt>WatershedOptions::RegionGrowing</tt> and
            <tt>WatershedOptions::UnionFind</tt>. The latter algorithms is fastest
            but doesn't support seeds and any other option.
        
            Default: RegionGrowing.
        */
    WatershedOptions & useMethod(Method method)
    {
        this->method = method;
        return *this;
    }
    
        /** \brief Use region-growing watershed.
        
            Use this method when you want to specify seeds explicitly (seeded watersheds) 
            or use any of the other options.
        
            Default: true.
        */
    WatershedOptions & regionGrowing()
    {
        method = RegionGrowing;
        return *this;
    }
    
        /** \brief Use union-find watershed.
        
            This is the fasted method, but it doesn't support seeds and any of the other 
            options (they will be silently ignored).
        
            Default: false.
        */
    WatershedOptions & unionFind()
    {
        method = UnionFind;
        return *this;
    }
};

namespace detail {

template <class CostType, class LabelType>
class WatershedStatistics
{
  public:
  
    typedef SeedRgDirectValueFunctor<CostType> value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    
    typedef CostType  first_argument_type;
    typedef LabelType second_argument_type;
    typedef LabelType argument_type;
    
    WatershedStatistics()
    {}

    void resize(unsigned int)
    {}

    void reset()
    {}

        /** update regions statistics (do nothing in the watershed algorithm)
        */
    template <class T1, class T2>
    void operator()(first_argument_type const &, second_argument_type const &) 
    {}

        /** ask for maximal index (label) allowed
        */
    LabelType maxRegionLabel() const
        { return size() - 1; }

        /** ask for array size (i.e. maxRegionLabel() + 1)
        */
    LabelType size() const
        { return NumericTraits<LabelType>::max(); }

        /** read the statistics functor for a region via its label
        */
    const_reference operator[](argument_type label) const
        { return stats; }

        /** access the statistics functor for a region via its label
        */
    reference operator[](argument_type label)
        { return stats; }

    value_type stats;
};

template <class Value>
class SeedRgBiasedValueFunctor
{
  public:
    double bias;

        /* the functor's argument type
        */
    typedef Value argument_type;

        /* the functor's result type (unused, only necessary for
            use of SeedRgDirectValueFunctor in \ref vigra::ArrayOfRegionStatistics
        */
    typedef Value result_type;

        /* the return type of the cost() function
        */
    typedef Value cost_type;
    
    SeedRgBiasedValueFunctor(double b = 1.0)
    : bias(b)
    {}

        /* Do nothing (since we need not update region statistics).
        */
    void operator()(argument_type const &) const {}

        /* Return scaled argument
        */
    cost_type cost(argument_type const & v) const
    {
        return cost_type(bias*v);
    }
};

template <class CostType, class LabelType>
class BiasedWatershedStatistics
{
  public:
  
    typedef SeedRgBiasedValueFunctor<CostType> value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    
    typedef CostType  first_argument_type;
    typedef LabelType second_argument_type;
    typedef LabelType argument_type;
    
    BiasedWatershedStatistics(LabelType biasedLabel, double bias)
    : biased_label(biasedLabel),
      biased_stats(bias)
    {}

    void resize(unsigned int)
    {}

    void reset()
    {}

        /** update regions statistics (do nothing in the watershed algorithm)
        */
    template <class T1, class T2>
    void operator()(first_argument_type const &, second_argument_type const &) 
    {}

        /** ask for maximal index (label) allowed
        */
    LabelType maxRegionLabel() const
        { return size() - 1; }

        /** ask for array size (i.e. maxRegionLabel() + 1)
        */
    LabelType size() const
        { return NumericTraits<LabelType>::max(); }

        /** read the statistics functor for a region via its label
        */
    const_reference operator[](argument_type label) const
    { 
        return (label == biased_label)
                    ? biased_stats
                    : stats; 
    }

        /** access the statistics functor for a region via its label
        */
    reference operator[](argument_type label)
    { 
        return (label == biased_label)
                    ? biased_stats
                    : stats; 
    }

    LabelType biased_label;
    value_type stats, biased_stats;
};

} // namespace detail

/** \brief Region segmentation by means of a flooding-based watershed algorithm.

    This function implements variants of the watershed algorithm
    described in

    L. Vincent and P. Soille: "<em>Watersheds in digital spaces: An efficient algorithm
    based on immersion simulations</em>", IEEE Trans. Patt. Analysis Mach. Intell. 13(6):583-598, 1991

    The source image is a boundary indicator such as the gaussianGradientMagnitude()
    or the trace of the \ref boundaryTensor(), and the destination is a label image
    designating membership of each pixel in one of the regions. Plateaus in the boundary
    indicator (i.e. regions of constant gray value) are handled via a Euclidean distance
    transform by default.
    
    By default, the destination image is assumed to hold seeds for a seeded watershed 
    transform. Seeds may, for example, be created by means of generateWatershedSeeds(). 
    Note that the seeds will be overridden with the final watershed segmentation.
    
    Alternatively, you may provide \ref SeedOptions in order to instruct 
    watershedsRegionGrowing() to generate its own seeds (it will call generateWatershedSeeds()
    internally). In that case, the destination image should be zero-initialized.
    
    You can specify the neighborhood system to be used by passing \ref FourNeighborCode 
    or \ref EightNeighborCode (default).
    
    Further options to be specified via \ref WatershedOptions are:
    
    <ul>
    <li> Whether to keep a 1-pixel-wide contour (with label 0) between regions or 
         perform complete grow (i.e. all pixels are assigned to a region).
    <li> Whether to stop growing when the boundaryness exceeds a threshold (remaining
         pixels keep label 0).
    <li> Whether to use a faster, but less powerful algorithm ("turbo algorithm"). It
         is faster because it orders pixels by means of a \ref BucketQueue (therefore,
         the boundary indicator must contain integers in the range 
         <tt>[0, ..., bucket_count-1]</tt>, where <tt>bucket_count</tt> is specified in
         the options object), it only supports complete growing (no contour between regions
         is possible), and it handles plateaus in a simplistic way. It also saves some
         memory because it allocates less temporary storage.
    <li> Whether one region (label) is to be preferred or discouraged by biasing its cost 
         with a given factor (smaller than 1 for preference, larger than 1 for discouragement).
    </ul>

    Note that VIGRA provides an alternative implementation of the watershed transform via
    \ref watershedsUnionFind(). 

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int
        watershedsRegionGrowing(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                                DestIterator upperleftd, DestAccessor da, 
                                Neighborhood neighborhood = EightNeighborCode(),
                                WatershedOptions const & options = WatershedOptions());

        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        unsigned int
        watershedsRegionGrowing(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                                DestIterator upperleftd, DestAccessor da, 
                                WatershedOptions const & options = WatershedOptions());
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class Neighborhood = EightNeighborCode>
        unsigned int
        watershedsRegionGrowing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                pair<DestIterator, DestAccessor> dest, 
                                Neighborhood neighborhood = EightNeighborCode(),
                                WatershedOptions const & options = WatershedOptions());
                                
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        unsigned int
        watershedsRegionGrowing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                pair<DestIterator, DestAccessor> dest, 
                                WatershedOptions const & options = WatershedOptions());
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/watersheds.hxx\><br>
    Namespace: vigra

    Example: watersheds of the gradient magnitude.

    \code
    vigra::BImage src(w, h);
    ... // read input data
    
    // compute gradient magnitude at scale 1.0 as a boundary indicator
    vigra::FImage gradMag(w, h);
    gaussianGradientMagnitude(srcImageRange(src), destImage(gradMag), 1.0);

    // example 1
    {
        // the pixel type of the destination image must be large enough to hold
        // numbers up to 'max_region_label' to prevent overflow
        vigra::IImage labeling(w, h);
        
        // call watershed algorithm for 4-neighborhood, leave a 1-pixel boundary between regions,
        // and autogenerate seeds from all gradient minima where the magnitude is below 2.0
        unsigned int max_region_label = 
              watershedsRegionGrowing(srcImageRange(gradMag), destImage(labeling),
                                      FourNeighborCode(),
                                      WatershedOptions().keepContours()
                                           .seedOptions(SeedOptions().minima().threshold(2.0)));
    }
    
    // example 2
    {
        vigra::IImage labeling(w, h);
        
        // compute seeds beforehand (use connected components of all pixels 
        // where the gradient  is below 4.0)
        unsigned int max_region_label = 
              generateWatershedSeeds(srcImageRange(gradMag), destImage(labeling),
                                     SeedOptions().levelSets(4.0));
        
        // quantize the gradient image to 256 gray levels
        vigra::BImage gradMag256(w, h);
        vigra::FindMinMax<float> minmax; 
        inspectImage(srcImageRange(gradMag), minmax); // find original range
        transformImage(srcImageRange(gradMag), destImage(gradMag256),
                       linearRangeMapping(minmax, 0, 255));
        
        // call the turbo algorithm with 256 bins, using 8-neighborhood
        watershedsRegionGrowing(srcImageRange(gradMag256), destImage(labeling),
                                WatershedOptions().turboAlgorithm(256));
    }
    
    // example 3
    {
        vigra::IImage labeling(w, h);
        
        .. // get seeds from somewhere, e.g. an interactive labeling program,
           // make sure that label 1 corresponds to the background
        
        // bias the watershed algorithm so that the background is preferred
        // by reducing the cost for label 1 to 90%
        watershedsRegionGrowing(srcImageRange(gradMag), destImage(labeling),
                                WatershedOptions().biasLabel(1, 0.9));
    }
    \endcode

    <b> Required Interface:</b>

    \code
    SrcIterator src_upperleft, src_lowerright;
    DestIterator dest_upperleft;

    SrcAccessor src_accessor;
    DestAccessor dest_accessor;

    // compare src values
    src_accessor(src_upperleft) <= src_accessor(src_upperleft)

    // set result
    int label;
    dest_accessor.set(label, dest_upperleft);
    \endcode
*/
doxygen_overloaded_function(template <...> unsigned int watershedsRegionGrowing)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
unsigned int
watershedsRegionGrowing(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da, 
                        Neighborhood neighborhood,
                        WatershedOptions const & options = WatershedOptions())
{
    typedef typename SrcAccessor::value_type ValueType; 
    typedef typename DestAccessor::value_type LabelType; 
    
    unsigned int max_region_label = 0;
    
    if(options.seed_options.mini != SeedOptions::Unspecified)
    {
        // we are supposed to compute seeds
        max_region_label = 
            generateWatershedSeeds(srcIterRange(upperlefts, lowerrights, sa), 
                                   destIter(upperleftd, da),
                                   neighborhood, options.seed_options);
    }
    
    if(options.biased_label != 0)
    {
        // create a statistics functor for biased region growing
        detail::BiasedWatershedStatistics<ValueType, LabelType> 
                                 regionstats(options.biased_label, options.bias);

        // perform region growing, starting from the seeds computed above
        if(options.bucket_count == 0)
        {
            max_region_label = 
            seededRegionGrowing(srcIterRange(upperlefts, lowerrights, sa),
                                srcIter(upperleftd, da),
                                destIter(upperleftd, da), 
                                regionstats, options.terminate, neighborhood, options.max_cost);
        }
        else
        {
            max_region_label = 
            fastSeededRegionGrowing(srcIterRange(upperlefts, lowerrights, sa),
                                    destIter(upperleftd, da), 
                                    regionstats, options.terminate, 
                                    neighborhood, options.max_cost, options.bucket_count);
        }
    }
    else
    {
        // create a statistics functor for region growing
        detail::WatershedStatistics<ValueType, LabelType> regionstats;

        // perform region growing, starting from the seeds computed above
        if(options.bucket_count == 0)
        {
            max_region_label = 
            seededRegionGrowing(srcIterRange(upperlefts, lowerrights, sa),
                                srcIter(upperleftd, da),
                                destIter(upperleftd, da), 
                                regionstats, options.terminate, neighborhood, options.max_cost);
        }
        else
        {
            max_region_label = 
            fastSeededRegionGrowing(srcIterRange(upperlefts, lowerrights, sa),
                                    destIter(upperleftd, da), 
                                    regionstats, options.terminate, 
                                    neighborhood, options.max_cost, options.bucket_count);
        }
    }
    
    return max_region_label;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
watershedsRegionGrowing(SrcIterator upperlefts, SrcIterator lowerrights, SrcAccessor sa,
                        DestIterator upperleftd, DestAccessor da, 
                        WatershedOptions const & options = WatershedOptions())
{
    return watershedsRegionGrowing(upperlefts, lowerrights, sa, upperleftd,  da,
                                   EightNeighborCode(), options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Neighborhood>
inline unsigned int
watershedsRegionGrowing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        pair<DestIterator, DestAccessor> dest, 
                        Neighborhood neighborhood,
                        WatershedOptions const & options = WatershedOptions())
{
    return watershedsRegionGrowing(src.first, src.second, src.third,
                                   dest.first, dest.second,    
                                   neighborhood, options);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline unsigned int
watershedsRegionGrowing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        pair<DestIterator, DestAccessor> dest, 
                        WatershedOptions const & options = WatershedOptions())
{
    return watershedsRegionGrowing(src.first, src.second, src.third,
                                    dest.first, dest.second,    
                                    EightNeighborCode(), options);
}


//@}

} // namespace vigra

#endif // VIGRA_WATERSHEDS_HXX
