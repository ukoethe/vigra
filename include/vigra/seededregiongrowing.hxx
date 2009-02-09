/************************************************************************/
/*                                                                      */
/*         Copyright 1998-2003 by Ullrich Koethe, Hans Meine            */
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

#ifndef VIGRA_SEEDEDREGIONGROWING_HXX
#define VIGRA_SEEDEDREGIONGROWING_HXX

#include <vector>
#include <stack>
#include <queue>
#include "utilities.hxx"
#include "stdimage.hxx"
#include "stdimagefunctions.hxx"

namespace vigra {

namespace detail {

template <class COST>
class SeedRgPixel
{
public:
    Point2D location_, nearest_;
    COST cost_;
    int count_;
    int label_;
    int dist_;

    SeedRgPixel()
    : location_(0,0), nearest_(0,0), cost_(0), count_(0), label_(0)
    {}

    SeedRgPixel(Point2D const & location, Point2D const & nearest,
                COST const & cost, int const & count, int const & label)
    : location_(location), nearest_(nearest),
      cost_(cost), count_(count), label_(label)
    {
        int dx = location_.x - nearest_.x;
        int dy = location_.y - nearest_.y;
        dist_ = dx * dx + dy * dy;
    }

    void set(Point2D const & location, Point2D const & nearest,
             COST const & cost, int const & count, int const & label)
    {
        location_ = location;
        nearest_ = nearest;
        cost_ = cost;
        count_ = count;
        label_ = label;

        int dx = location_.x - nearest_.x;
        int dy = location_.y - nearest_.y;
        dist_ = dx * dx + dy * dy;
    }

    struct Compare
    {
        // must implement > since priority_queue looks for largest element
        bool operator()(SeedRgPixel const & l,
                        SeedRgPixel const & r) const
        {
            if(r.cost_ == l.cost_)
            {
                if(r.dist_ == l.dist_) return r.count_ < l.count_;

                return r.dist_ < l.dist_;
            }

            return r.cost_ < l.cost_;
        }
        bool operator()(SeedRgPixel const * l,
                        SeedRgPixel const * r) const
        {
            if(r->cost_ == l->cost_)
            {
                if(r->dist_ == l->dist_) return r->count_ < l->count_;

                return r->dist_ < l->dist_;
            }

            return r->cost_ < l->cost_;
        }
    };

    struct Allocator
    {
        ~Allocator()
        {
            while(!freelist_.empty())
            {
                delete freelist_.top();
                freelist_.pop();
            }
        }

        SeedRgPixel *
        create(Point2D const & location, Point2D const & nearest,
               COST const & cost, int const & count, int const & label)
        {
            if(!freelist_.empty())
            {
                SeedRgPixel * res = freelist_.top();
                freelist_.pop();
                res->set(location, nearest, cost, count, label);
                return res;
            }

            return new SeedRgPixel(location, nearest, cost, count, label);
        }

        void dismiss(SeedRgPixel * p)
        {
            freelist_.push(p);
        }

        std::stack<SeedRgPixel<COST> *> freelist_;
    };
};

struct UnlabelWatersheds
{
    int operator()(int label) const
    {
        return label < 0 ? 0 : label;
    }
};

} // namespace detail

enum SRGType { KeepContours, CompleteGrow, SRGWatershedLabel = -1 };

/** \addtogroup SeededRegionGrowing Region Segmentation Algorithms
    Region growing, watersheds, and voronoi tesselation
*/
//@{

/********************************************************/
/*                                                      */
/*                    seededRegionGrowing               */
/*                                                      */
/********************************************************/

/** \brief Region Segmentation by means of Seeded Region Growing.

    This algorithm implements seeded region growing as described in

    R. Adams, L. Bischof: "<em> Seeded Region Growing</em>", IEEE Trans. on Pattern
    Analysis and Maschine Intelligence, vol 16, no 6, 1994, and

    Ullrich K&ouml;the:
    <em><a href="http://kogs-www.informatik.uni-hamburg.de/~koethe/papers/#primary">Primary Image Segmentation</a></em>,
    in: G. Sagerer, S.
    Posch, F. Kummert (eds.): Mustererkennung 1995, Proc. 17. DAGM-Symposium,
    Springer 1995

    The seed image is a partly segmented image which contains uniquely
    labeled regions (the seeds) and unlabeled pixels (the candidates, label 0).
    Seed regions can be as large as you wish and as small as one pixel. If
    there are no candidates, the algorithm will simply copy the seed image
    into the output image. Otherwise it will aggregate the candidates into
    the existing regions so that a cost function is minimized. This
    works as follows:

    <ol>

    <li> Find all candidate pixels that are 4-adjacent to a seed region.
    Calculate the cost for aggregating each candidate into its adajacent region
    and put the candidates into a priority queue.

    <li> While( priority queue is not empty)

        <ol>

        <li> Take the candidate with least cost from the queue. If it has not
        already been merged, merge it with it's adjacent region.

        <li> Put all candidates that are 4-adjacent to the pixel just processed
        into the priority queue.

        </ol>

    </ol>

    If <tt>SRGType == CompleteGrow</tt> (the default), this algorithm
    will produce a complete 4-connected tesselation of the image.
    If <tt>SRGType == KeepContours</tt>, a one-pixel-wide border will be left
    between the regions. The border pixels get label 0 (zero).

    The cost is determined jointly by the source image and the
    region statistics functor. The source image contains feature values for each
    pixel which will be used by the region statistics functor to calculate and
    update statistics for each region and to calculate the cost for each
    candidate. The <TT>RegionStatisticsArray</TT> must be compatible to the
    \ref ArrayOfRegionStatistics functor and contains an <em> array</em> of
    statistics objects for each region. The indices must correspond to the
    labels of the seed regions. The statistics for the initial regions must have
    been calculated prior to calling <TT>seededRegionGrowing()</TT> (for example by
    means of \ref inspectTwoImagesIf()).

    For each candidate
    <TT>x</TT> that is adjacent to region <TT>i</TT>, the algorithm will call
    <TT>stats[i].cost(as(x))</TT> to get the cost (where <TT>x</TT> is a <TT>SrcImageIterator</TT>
    and <TT>as</TT> is
    the SrcAccessor). When a candidate has been merged with a region, the
    statistics are updated by calling <TT>stats[i].operator()(as(x))</TT>. Since
    the <TT>RegionStatisticsArray</TT> is passed by reference, this will overwrite
    the original statistics.

    If a candidate could be merged into more than one regions with identical
    cost, the algorithm will favour the nearest region.

    In some cases, the cost only depends on the feature value of the current
    pixel. Then the update operation will simply be a no-op, and the <TT>cost()</TT>
    function returns its argument. This behavior is implemented by the
    \ref SeedRgDirectValueFunctor. With <tt>SRGType == KeepContours</tt>,
    this is equivalent to the watershed algorithm.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray>
        void seededRegionGrowing(SrcImageIterator srcul,
                                 SrcImageIterator srclr, SrcAccessor as,
                                 SeedImageIterator seedsul, SeedAccessor aseeds,
                                 DestImageIterator destul, DestAccessor ad,
                                 RegionStatisticsArray & stats,
                                 SRGType srgType = CompleteGrow);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray>
        void
        seededRegionGrowing(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> img1,
                            pair<SeedImageIterator, SeedAccessor> img3,
                            pair<DestImageIterator, DestAccessor> img4,
                            RegionStatisticsArray & stats,
                            SRGType srgType = CompleteGrow);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<<a href="seededregiongrowing_8hxx-source.html">vigra/seededregiongrowing.hxx</a>\><br>
    Namespace: vigra

    Example: implementation of the voronoi tesselation

    \code
    vigra::BImage points(w,h);
    vigra::FImage dist(x,y);

    // empty edge image
    points = 0;
    dist = 0;

    int max_region_label = 100;

    // throw in some random points:
    for(int i = 1; i <= max_region_label; ++i)
           points(w * rand() / RAND_MAX , h * rand() / RAND_MAX) = i;

    // calculate Euclidean distance transform
    vigra::distanceTransform(srcImageRange(points), destImage(dist), 2);

    // init statistics functor
    vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> >
                                              stats(max_region_label);

    // find voronoi region of each point
   vigra:: seededRegionGrowing(srcImageRange(dist), srcImage(points),
                               destImage(points), stats);
    \endcode

    <b> Required Interface:</b>

    \code
    SrcImageIterator src_upperleft, src_lowerright;
    SeedImageIterator seed_upperleft;
    DestImageIterator dest_upperleft;

    SrcAccessor src_accessor;
    SeedAccessor seed_accessor;
    DestAccessor dest_accessor;

    RegionStatisticsArray stats;

    // calculate costs
    RegionStatisticsArray::value_type::cost_type cost =
        stats[seed_accessor(seed_upperleft)].cost(src_accessor(src_upperleft));

    // compare costs
    cost < cost;

    // update statistics
    stats[seed_accessor(seed_upperleft)](src_accessor(src_upperleft));

    // set result
    dest_accessor.set(seed_accessor(seed_upperleft), dest_upperleft);
    \endcode

    Further requirements are determined by the <TT>RegionStatisticsArray</TT>.
*/
doxygen_overloaded_function(template <...> void seededRegionGrowing)

template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
void seededRegionGrowing(SrcImageIterator srcul,
                         SrcImageIterator srclr, SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats,
                         const SRGType srgType)
{
    int w = srclr.x - srcul.x;
    int h = srclr.y - srcul.y;
    int count = 0;

    SrcImageIterator isy = srcul, isx = srcul;  // iterators for the src image

    typedef typename RegionStatisticsArray::value_type RegionStatistics;
    typedef typename RegionStatistics::cost_type CostType;
    typedef detail::SeedRgPixel<CostType> Pixel;

    typename Pixel::Allocator allocator;

    typedef std::priority_queue<Pixel *, std::vector<Pixel *>,
                                typename Pixel::Compare>  SeedRgPixelHeap;

    // copy seed image in an image with border
    IImage regions(w+2, h+2);
    IImage::Iterator ir = regions.upperLeft() + Diff2D(1,1);
    IImage::Iterator iry, irx;

    initImageBorder(destImageRange(regions), 1, SRGWatershedLabel);
    copyImage(seedsul, seedsul+Diff2D(w,h), aseeds, ir, regions.accessor());

    // allocate and init memory for the results

    SeedRgPixelHeap pheap;
    int cneighbor;

    static const Diff2D dist[] = { Diff2D(-1,0), Diff2D(0,-1),
                                   Diff2D(1,0),  Diff2D(0,1) };

    Point2D pos(0,0);
    for(isy=srcul, iry=ir, pos.y=0; pos.y<h;
        ++pos.y, ++isy.y, ++iry.y)
    {
        for(isx=isy, irx=iry, pos.x=0; pos.x<w;
            ++pos.x, ++isx.x, ++irx.x)
        {
            if(*irx == 0)
            {
                // find candidate pixels for growing and fill heap
                for(int i=0; i<4; i++)
                {
                    cneighbor = irx[dist[i]];
                    if(cneighbor > 0)
                    {
                        CostType cost = stats[cneighbor].cost(as(isx));

                        Pixel * pixel =
                            allocator.create(pos, pos+dist[i], cost, count++, cneighbor);
                        pheap.push(pixel);
                    }
                }
            }
        }
    }

    // perform region growing
    while(pheap.size() != 0)
    {
        Pixel * pixel = pheap.top();
        pheap.pop();

        Point2D pos = pixel->location_;
        Point2D nearest = pixel->nearest_;
        int lab = pixel->label_;

        allocator.dismiss(pixel);

        irx = ir + pos;
        isx = srcul + pos;

        if(*irx) // already labelled region / watershed?
            continue;

        if(srgType == KeepContours)
        {
            for(int i=0; i<4; i++)
            {
                cneighbor = irx[dist[i]];
                if((cneighbor>0) && (cneighbor != lab))
                {
                    lab = SRGWatershedLabel;
                    break;
                }
            }
        }

        *irx = lab;

        if((srgType != KeepContours) || (lab > 0))
        {
            // update statistics
            stats[*irx](as(isx));

            // search neighborhood
            // second pass: find new candidate pixels
            for(int i=0; i<4; i++)
            {
                if(irx[dist[i]] == 0)
                {
                    CostType cost = stats[lab].cost(as(isx, dist[i]));

                    Pixel * new_pixel =
                        allocator.create(pos+dist[i], nearest, cost, count++, lab);
                    pheap.push(new_pixel);
                }
            }
        }
    }

    // write result
    transformImage(ir, ir+Point2D(w,h), regions.accessor(), destul, ad,
                   detail::UnlabelWatersheds());
}

template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing(SrcImageIterator srcul,
                    SrcImageIterator srclr, SrcAccessor as,
                    SeedImageIterator seedsul, SeedAccessor aseeds,
                    DestImageIterator destul, DestAccessor ad,
                    RegionStatisticsArray & stats)
{
    seededRegionGrowing(srcul, srclr, as,
                        seedsul, aseeds,
                        destul, ad,
                        stats, CompleteGrow);
}

template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats,
                    SRGType srgType)
{
    seededRegionGrowing(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, srgType);
}

template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats)
{
    seededRegionGrowing(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, CompleteGrow);
}

/********************************************************/
/*                                                      */
/*               SeedRgDirectValueFunctor               */
/*                                                      */
/********************************************************/

/** \brief Statistics functor to be used for seeded region growing.

    This functor can be used if the cost of a candidate during
    \ref seededRegionGrowing() is equal to the feature value of that
    candidate and does not depend on properties of the region it is going to
    be merged with.

    <b>\#include</b> \<<a href="seededregiongrowing_8hxx-source.html">vigra/seededregiongrowing.hxx</a>\><br>
    Namespace: vigra


     <b> Required Interface:</b>

     no requirements
*/
template <class Value>
class SeedRgDirectValueFunctor
{
  public:
        /** the functor's argument type
        */
    typedef Value argument_type;

        /** the functor's result type (unused, only necessary for
            use of SeedRgDirectValueFunctor in \ref vigra::ArrayOfRegionStatistics
        */
    typedef Value result_type;

        /** \deprecated use argument_type
        */
    typedef Value value_type;

        /** the return type of the cost() function
        */
    typedef Value cost_type;

        /** Do nothing (since we need not update region statistics).
        */
    void operator()(argument_type const &) const {}

        /** Return argument (since cost is identical to feature value)
        */
    cost_type const & cost(argument_type const & v) const
    {
        return v;
    }
};

//@}

} // namespace vigra

#endif // VIGRA_SEEDEDREGIONGROWING_HXX

