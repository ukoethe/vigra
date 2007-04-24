/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2007 by Kasim Terzic, Christian-Dennis Rahn       */
/*                        and Ullrich Koethe                            */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_SEEDEDREGIONGROWING_3D_HXX
#define VIGRA_SEEDEDREGIONGROWING_3D_HXX

#include <vector>
#include <stack>
#include <queue>
#include "utilities.hxx"
#include "stdimage.hxx"
#include "stdimagefunctions.hxx"
#include "seededregiongrowing.hxx"
#include "multi_pointoperators.hxx"

namespace vigra {

namespace detail {

template <class COST, class Diff_type>
class SeedRgVoxel
{
public:
    Diff_type location_, nearest_;
    COST cost_;
    int count_;
    int label_;
    int dist_;

    SeedRgVoxel()
    //: location_(0,0,0), nearest_(0,0,0), cost_(0), count_(0), label_(0)
    {
        location_ = Diff_type(0,0,0);
        nearest_ = Diff_type(0,0,0);
        cost_ = 0;
        count_ = 0;
        label_ = 0;
    }

    SeedRgVoxel(Diff_type const & location, Diff_type const & nearest,
                COST const & cost, int const & count, int const & label)
    : location_(location), nearest_(nearest),
      cost_(cost), count_(count), label_(label)
    {
        int dx = location_[0] - nearest_[0];
        int dy = location_[1] - nearest_[1];
        int dz = location_[2] - nearest_[2];
        dist_ = dx * dx + dy * dy + dz * dz;
    }

    void set(Diff_type const & location, Diff_type const & nearest,
             COST const & cost, int const & count, int const & label)
    {
        location_ = location;
        nearest_ = nearest;
        cost_ = cost;
        count_ = count;
        label_ = label;

        int dx = location_[0] - nearest_[0];
        int dy = location_[1] - nearest_[1];
        int dz = location_[2] - nearest_[2];
        dist_ = dx * dx + dy * dy + dz * dz;
    }

    struct Compare
    {
        // must implement > since priority_queue looks for largest element
        bool operator()(SeedRgVoxel const & l,
                        SeedRgVoxel const & r) const
        {
            if(r.cost_ == l.cost_)
            {
                if(r.dist_ == l.dist_) return r.count_ < l.count_;

                return r.dist_ < l.dist_;
            }

            return r.cost_ < l.cost_;
        }
        bool operator()(SeedRgVoxel const * l,
                        SeedRgVoxel const * r) const
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

        SeedRgVoxel * create(Diff_type const & location, Diff_type const & nearest,
                             COST const & cost, int const & count, int const & label)
        {
            if(!freelist_.empty())
            {
                SeedRgVoxel * res = freelist_.top();
                freelist_.pop();
                res->set(location, nearest, cost, count, label);
                return res;
            }

            return new SeedRgVoxel(location, nearest, cost, count, label);
        }

        void dismiss(SeedRgVoxel * p)
        {
            freelist_.push(p);
        }

        std::stack<SeedRgVoxel<COST,Diff_type> *> freelist_;
    };
};

} // namespace detail

/** \addtogroup SeededRegionGrowing Seeded Region Growing
    Region segmentation and voronoi tesselation
*/
//@{

/********************************************************/
/*                                                      */
/*                    seededRegionGrowing3d             */
/*                                                      */
/********************************************************/

/** \brief Three-dimensional Region Segmentation by means of Seeded Region Growing.

    This algorithm implements seeded region growing as described in

    The seed image is a partly segmented multi-dimensional array which contains uniquely
    labeled regions (the seeds) and unlabeled voxels (the candidates, label 0).
    Seed regions can be as large as you wish and as small as one voxel. If
    there are no candidates, the algorithm will simply copy the seed array
    into the output array. Otherwise it will aggregate the candidates into
    the existing regions so that a cost function is minimized. This
    works as follows:

    <ol>

    <li> Find all candidate pixels that are 6-adjacent to a seed region.
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

    If <tt>SRGType == CompleteGrow</tt> (the default), this algorithm will 
    produce a complete 6-connected tesselation of the array.
    Other grow types (such as keeping contours for watersheds) are currently not
    supported 

    The cost is determined jointly by the source array and the
    region statistics functor. The source array contains feature values for each
    pixel which will be used by the region statistics functor to calculate and
    update statistics for each region and to calculate the cost for each
    candidate. The <TT>RegionStatisticsArray</TT> must be compatible to the
    \ref ArrayOfRegionStatistics functor and contains an <em> array</em> of
    statistics objects for each region. The indices must correspond to the
    labels of the seed regions. The statistics for the initial regions must have
    been calculated prior to calling <TT>seededRegionGrowing3d()</TT>

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
    voxel. Then the update operation will simply be a no-op, and the <TT>cost()</TT>
    function returns its argument. This behavior is implemented by the
    \ref SeedRgDirectValueFunctor.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class Diff_type, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray >
        void seededRegionGrowing3d(SrcImageIterator srcul, Diff_type shape,
                                 SrcAccessor as,
                                 SeedImageIterator seedsul, SeedAccessor aseeds,
                                 DestImageIterator destul, DestAccessor ad,
                                 RegionStatisticsArray & stats, 
                                 const SRGType srgType == CompleteGrow)

       }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
            template <class SrcImageIterator, class Shape, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray>
            inline void
            seededRegionGrowing3d(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                            pair<SeedImageIterator, SeedAccessor> img3,
                            pair<DestImageIterator, DestAccessor> img4,
                            RegionStatisticsArray & stats, 
                            const SRGType srgType == CompleteGrow)
    }
    \endcode

*/

/*
template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
void seededRegionGrowing3d(SrcImageIterator srcul,
                         SrcImageIterator srclr, SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats,
                         const SRGType srgType)
*/
template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
void seededRegionGrowing3d(SrcImageIterator srcul, Diff_type shape,
                         SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats, int max_cost, 
                         const SRGType srgType)
{
    SrcImageIterator srclr = srcul + shape;
    //int w = srclr.x - srcul.x;
    int w = shape[0];
    //int h = srclr.y - srcul.y;
    int h = shape[1];
    //int d = srclr.z - srcul.z;
    int d = shape[2];
    int count = 0;

    SrcImageIterator isy = srcul, isx = srcul, isz = srcul;  // iterators for the src image

    typedef typename RegionStatisticsArray::value_type RegionStatistics;
    typedef typename RegionStatistics::cost_type CostType;
    typedef detail::SeedRgVoxel<CostType, Diff_type> Voxel;

    typename Voxel::Allocator allocator;

    typedef std::priority_queue< Voxel *,
                                 std::vector<Voxel *>,
                                 typename Voxel::Compare >  SeedRgVoxelHeap;
    typedef MultiArray<3, int> IVolume;

    // copy seed image in an image with border
    Diff_type regionshape = shape + Diff_type(2,2,2);
    IVolume regions(regionshape);
    MultiIterator<3,int> ir = regions.traverser_begin();
    ir = ir + Diff_type(1,1,1);
    
    //IVolume::Iterator iry, irx, irz;
    MultiIterator<3,int> iry, irx, irz;

    //initImageBorder(destImageRange(regions), 1, SRGWatershedLabel);
    initMultiArrayBorder(destMultiArrayRange(regions), 1, SRGWatershedLabel); 
    
    copyMultiArray(seedsul, Diff_type(w,h,d), aseeds, ir, AccessorTraits<int>::default_accessor()/* vigra::StandardValueAccessor<int>*/);

    // allocate and init memory for the results

    SeedRgVoxelHeap pheap;
    int cneighbor;

    static const Diff_type dist[] = { Diff_type(-1, 0, 0), Diff_type( 0,-1, 0),
                                      Diff_type( 1, 0, 0), Diff_type( 0, 1, 0),
                                      Diff_type( 0, 0,-1), Diff_type( 0, 0, 1) };
    
    Diff_type pos(0,0,0);

    for(isz=srcul, irz=ir, pos[2]=0; pos[2]<d;
            pos[2]++, isz.dim2()++, irz.dim2()++)
    {
        //std::cerr << "Z = " << pos[2] << std::endl;

        for(isy=isz, iry=irz, pos[1]=0; pos[1]<h;
            pos[1]++, isy.dim1()++, iry.dim1()++)
        {
            //std::cerr << "Y = " << pos[1] << std::endl;
            
            for(isx=isy, irx=iry, pos[0]=0; pos[0]<w;
                pos[0]++, isx.dim0()++, irx.dim0()++)
            {
                //std::cerr << "X = " << pos[0] << std::endl;
                
                if(*irx == 0)
                {
                    // find candidate pixels for growing and fill heap
                    for(int i=0; i<6; i++)
                    {
                        cneighbor = *(irx + dist[i]);
                        if(cneighbor > 0)
                        {
                            CostType cost = stats[cneighbor].cost(as(isx));

                            Voxel * voxel =
                                allocator.create(pos, pos+dist[i], cost, count++, cneighbor);
                            pheap.push(voxel);
                        }
                    }
                }
            }
        }
    }
    
    // perform region growing
    while(pheap.size() != 0)
    {
        Voxel * voxel = pheap.top();
        pheap.pop();

        if(max_cost > 0 && voxel->cost_ > max_cost) break;

        Diff_type pos = voxel->location_;
        Diff_type nearest = voxel->nearest_;
        int lab = voxel->label_;

        allocator.dismiss(voxel);

        irx = ir + pos;
        isx = srcul + pos;

        if(*irx) // already labelled region / watershed?
            continue;

        if(srgType == KeepContours)
        {
            for(int i=0; i<6; i++)
            {
                cneighbor = * (irx + dist[i]);
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
            for(int i=0; i<6; i++)
            {
                if(*(irx + dist[i]) == 0)
                {
                    CostType cost = stats[lab].cost(as(isx, dist[i]));

                    Voxel * new_voxel =
                        allocator.create(pos+dist[i], nearest, cost, count++, lab);
                    pheap.push(new_voxel);
                }
            }
        }
    }

    // write result
    transformMultiArray(ir, Diff_type(w,h,d), AccessorTraits<int>::default_accessor()/* regions.accessor()*/, destul, ad,
                        detail::UnlabelWatersheds());
}
/*
template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(SrcImageIterator srcul,
                    SrcImageIterator srclr, SrcAccessor as,
                    SeedImageIterator seedsul, SeedAccessor aseeds,
                    DestImageIterator destul, DestAccessor ad,
                    RegionStatisticsArray & stats)
{
    seededRegionGrowing3d(srcul, srclr, as,
                        seedsul, aseeds,
                        destul, ad,
                        stats, CompleteGrow);
}*/
/*
template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats,
                    SRGType srgType)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, srgType);
}

template <class SrcImageIterator, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, CompleteGrow);
}
*/

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
void seededRegionGrowing3d(SrcImageIterator srcul, Diff_type shape,
                         SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats, int max_cost)
{
    seededRegionGrowing3d( srcul, shape, as, seedsul, aseeds, destul, ad, stats, max_cost, CompleteGrow);
}

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
void seededRegionGrowing3d(SrcImageIterator srcul, Diff_type shape,
                         SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats)
{
    seededRegionGrowing3d( srcul, shape, as, seedsul, aseeds, destul, ad, stats, -1, CompleteGrow);
}

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
void seededRegionGrowing3d(SrcImageIterator srcul, Diff_type shape,
                         SrcAccessor as,
                         SeedImageIterator seedsul, SeedAccessor aseeds,
                         DestImageIterator destul, DestAccessor ad,
                         RegionStatisticsArray & stats, SRGType srgType)
{
    seededRegionGrowing3d( srcul, shape, as, seedsul, aseeds, destul, ad, stats, -1, srgType);
}



template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats, int max_cost)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, max_cost, CompleteGrow);
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats, int max_cost, const SRGType srgType)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, max_cost, srgType);
}


template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats, const SRGType srgType)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, -1, srgType);
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3d(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                    pair<SeedImageIterator, SeedAccessor> img3,
                    pair<DestImageIterator, DestAccessor> img4,
                    RegionStatisticsArray & stats)
{
    seededRegionGrowing3d(img1.first, img1.second, img1.third,
                        img3.first, img3.second,
                        img4.first, img4.second,
                        stats, -1, CompleteGrow);
}



} // namespace vigra

#endif // VIGRA_SEEDEDREGIONGROWING_HXX

