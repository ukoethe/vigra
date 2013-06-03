/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2007 by Kasim Terzic, Christian-Dennis Rahn       */
/*                        and Ullrich Koethe                            */
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

#ifndef VIGRA_SEEDEDREGIONGROWING_3D_HXX
#define VIGRA_SEEDEDREGIONGROWING_3D_HXX

#include <vector>
#include <stack>
#include <queue>
#include "utilities.hxx"
#include "stdimage.hxx"
#include "stdimagefunctions.hxx"
#include "seededregiongrowing.hxx"
#include "multi_shape.hxx"
#include "multi_pointoperators.hxx"
#include "voxelneighborhood.hxx"

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

/** \addtogroup SeededRegionGrowing
*/
//@{

/********************************************************/
/*                                                      */
/*                    seededRegionGrowing3D             */
/*                                                      */
/********************************************************/

/** \brief Three-dimensional Region Segmentation by means of Seeded Region Growing.

    This algorithm implements seeded region growing as described in

    The seed image is a partly segmented multi-dimensional array which contains uniquely
    labeled regions (the seeds) and unlabeled voxels (the candidates, label 0).
    Seed regions can be as large as you wish and as small as one voxel. If
    there are no candidates, the algorithm will simply copy the seed array
    into the output array. Otherwise it will aggregate the candidates into
    the existing regions so that a cost function is minimized.
    Candidates are taken from the neighborhood of the already assigned pixels, 
    where the type of neighborhood is determined by parameter <tt>neighborhood</tt>
    which can take the values <tt>NeighborCode3DSix()</tt> (the default) 
    or <tt>NeighborCode3DTwentySix()</tt>. The algorithm basically works as follows 
    (illustrated for 6-neighborhood, but 26-neighborhood works in the same way):

    <ol>

    <li> Find all candidate pixels that are 6-adjacent to a seed region.
    Calculate the cost for aggregating each candidate into its adjacent region
    and put the candidates into a priority queue.

    <li> While( priority queue is not empty)

        <ol>

        <li> Take the candidate with least cost from the queue. If it has not
        already been merged, merge it with it's adjacent region.

        <li> Put all candidates that are 4-adjacent to the pixel just processed
        into the priority queue.

        </ol>

    </ol>

    <tt>SRGType</tt> can take the following values:
    
    <DL>
    <DT><tt>CompleteGrow</tt> <DD> produce a complete tesselation of the volume (default).
    <DT><tt>KeepContours</tt> <DD> keep a 1-voxel wide unlabeled contour between all regions.
    <DT><tt>StopAtThreshold</tt> <DD> stop when the boundary indicator values exceed the 
                             threshold given by parameter <tt>max_cost</tt>.
    <DT><tt>KeepContours | StopAtThreshold</tt> <DD> keep 1-voxel wide contour and stop at given <tt>max_cost</tt>.
    </DL>

    The cost is determined jointly by the source array and the
    region statistics functor. The source array contains feature values for each
    pixel which will be used by the region statistics functor to calculate and
    update statistics for each region and to calculate the cost for each
    candidate. The <TT>RegionStatisticsArray</TT> must be compatible to the
    \ref ArrayOfRegionStatistics functor and contains an <em> array</em> of
    statistics objects for each region. The indices must correspond to the
    labels of the seed regions. The statistics for the initial regions must have
    been calculated prior to calling <TT>seededRegionGrowing3D()</TT>

    For each candidate
    <TT>x</TT> that is adjacent to region <TT>i</TT>, the algorithm will call
    <TT>stats[i].cost(as(x))</TT> to get the cost (where <TT>x</TT> is a <TT>SrcImageIterator</TT>
    and <TT>as</TT> is
    the SrcAccessor). When a candidate has been merged with a region, the
    statistics are updated by calling <TT>stats[i].operator()(as(x))</TT>. Since
    the <TT>RegionStatisticsArray</TT> is passed by reference, this will overwrite
    the original statistics.

    If a candidate could be merged into more than one regions with identical
    cost, the algorithm will favour the nearest region. If <tt>StopAtThreshold</tt> is active, 
    and the cost of the current candidate at any point in the algorithm exceeds the optional 
    <tt>max_cost</tt> value (which defaults to <tt>NumericTraits<double>::max()</tt>), 
    region growing is aborted, and all voxels not yet assigned to a region remain unlabeled.

    In some cases, the cost only depends on the feature value of the current
    voxel. Then the update operation will simply be a no-op, and the <TT>cost()</TT>
    function returns its argument. This behavior is implemented by the
    \ref SeedRgDirectValueFunctor.

    <b> Declarations:</b>

    pass 3D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class TS, class AS,
                  class T2, class S2,
                  class RegionStatisticsArray, class Neighborhood>
        void
        seededRegionGrowing3D(MultiArrayView<3, T1, S1> const & src,
                              MultiArrayView<3, TS, AS> const & seeds,
                              MultiArrayView<3, T2, S2>         labels,
                              RegionStatisticsArray &           stats, 
                              SRGType                           srgType = CompleteGrow,
                              Neighborhood                      neighborhood = NeighborCode3DSix(),
                              double                            max_cost = NumericTraits<double>::max());
    }
    \endcode

    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class Shape, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray, class Neighborhood>
        void 
        seededRegionGrowing3D(SrcImageIterator srcul, Shape shape, SrcAccessor as,
                              SeedImageIterator seedsul, SeedAccessor aseeds,
                              DestImageIterator destul, DestAccessor ad,
                              RegionStatisticsArray & stats, 
                              SRGType srgType = CompleteGrow,
                              Neighborhood neighborhood = NeighborCode3DSix(),
                              double max_cost = NumericTraits<double>::max());
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class Shape, class SrcAccessor,
                  class SeedImageIterator, class SeedAccessor,
                  class DestImageIterator, class DestAccessor,
                  class RegionStatisticsArray, class Neighborhood>
        void
        seededRegionGrowing3D(triple<SrcImageIterator, Shape, SrcAccessor> src,
                              pair<SeedImageIterator, SeedAccessor> seeds,
                              pair<DestImageIterator, DestAccessor> dest,
                              RegionStatisticsArray & stats, 
                              SRGType srgType = CompleteGrow,
                              Neighborhood neighborhood = NeighborCode3DSix(), 
                              double max_cost = NumericTraits<double>::max());
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void seededRegionGrowing3D)

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray, class Neighborhood>
void 
seededRegionGrowing3D(SrcImageIterator srcul, Diff_type shape, SrcAccessor as,
                      SeedImageIterator seedsul, SeedAccessor aseeds,
                      DestImageIterator destul, DestAccessor ad,
                      RegionStatisticsArray & stats, 
                      SRGType srgType,
                      Neighborhood,
                      double max_cost)
{
    //SrcImageIterator srclr = srcul + shape;
    //int w = srclr.x - srcul.x;
    int w = shape[0];
    //int h = srclr.y - srcul.y;
    int h = shape[1];
    //int d = srclr.z - srcul.z;
    int d = shape[2];
    int count = 0;

    SrcImageIterator isy = srcul, isx = srcul, isz = srcul;  // iterators for the src image

    typedef typename RegionStatisticsArray::value_type RegionStatistics;
    typedef typename PromoteTraits<typename RegionStatistics::cost_type, double>::Promote CostType;
    typedef detail::SeedRgVoxel<CostType, Diff_type> Voxel;

    typename Voxel::Allocator allocator;

    typedef std::priority_queue< Voxel *,
                                 std::vector<Voxel *>,
                                 typename Voxel::Compare >  SeedRgVoxelHeap;
    typedef MultiArray<3, int> IVolume;
    typedef IVolume::traverser Traverser;

    // copy seed image in an image with border
    Diff_type regionshape = shape + Diff_type(2,2,2);
    IVolume regions(regionshape);
    Traverser ir = regions.traverser_begin();
    ir = ir + Diff_type(1,1,1);
    
    //IVolume::Iterator iry, irx, irz;
   Traverser iry, irx, irz;

    //initImageBorder(destImageRange(regions), 1, SRGWatershedLabel);
    initMultiArrayBorder(destMultiArrayRange(regions), 1, SRGWatershedLabel); 
    
    copyMultiArray(seedsul, Diff_type(w,h,d), aseeds, ir, AccessorTraits<int>::default_accessor());

    // allocate and init memory for the results

    SeedRgVoxelHeap pheap;
    int cneighbor;

    typedef typename Neighborhood::Direction Direction;
    int directionCount = Neighborhood::DirectionCount;

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
                    for(int i=0; i<directionCount; i++)
                    {
                        cneighbor = *(irx + Neighborhood::diff((Direction)i));
                        if(cneighbor > 0)
                        {
                            CostType cost = stats[cneighbor].cost(as(isx));

                            Voxel * voxel =
                                allocator.create(pos, pos+Neighborhood::diff((Direction)i), cost, count++, cneighbor);
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

        Diff_type pos = voxel->location_;
        Diff_type nearest = voxel->nearest_;
        int lab = voxel->label_;
        CostType cost = voxel->cost_;

        allocator.dismiss(voxel);

        if((srgType & StopAtThreshold) != 0 && cost > max_cost)
            break;

        irx = ir + pos;
        isx = srcul + pos;

        if(*irx) // already labelled region / watershed?
            continue;

        if((srgType & KeepContours) != 0)
        {
            for(int i=0; i<directionCount; i++)
            {
                cneighbor = * (irx + Neighborhood::diff((Direction)i));
                if((cneighbor>0) && (cneighbor != lab))
                {
                    lab = SRGWatershedLabel;
                    break;
                }
            }
        }

        *irx = lab;

        if((srgType & KeepContours) == 0 || lab > 0)
        {
            // update statistics
            stats[*irx](as(isx));

            // search neighborhood
            // second pass: find new candidate pixels
            for(int i=0; i<directionCount; i++)
            {
                if(*(irx + Neighborhood::diff((Direction)i)) == 0)
                {
                    CostType cost = stats[lab].cost(as(isx, Neighborhood::diff((Direction)i)));

                    Voxel * new_voxel =
                        allocator.create(pos+Neighborhood::diff((Direction)i), nearest, cost, count++, lab);
                    pheap.push(new_voxel);
                }
            }
        }
    }
    
    // free temporary memory
    while(pheap.size() != 0)
    {
        allocator.dismiss(pheap.top());
        pheap.pop();
    }

    // write result
    transformMultiArray(ir, Diff_type(w,h,d), AccessorTraits<int>::default_accessor(), 
                        destul, ad, detail::UnlabelWatersheds());
}

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray, class Neighborhood >
inline void
seededRegionGrowing3D(SrcImageIterator srcul, Diff_type shape, SrcAccessor as,
                      SeedImageIterator seedsul, SeedAccessor aseeds,
                      DestImageIterator destul, DestAccessor ad,
                      RegionStatisticsArray & stats, SRGType srgType, Neighborhood n)
{
    seededRegionGrowing3D( srcul, shape, as, seedsul, aseeds, 
                           destul, ad, stats, srgType, n, NumericTraits<double>::max());
}

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
inline void
seededRegionGrowing3D(SrcImageIterator srcul, Diff_type shape, SrcAccessor as,
                      SeedImageIterator seedsul, SeedAccessor aseeds,
                      DestImageIterator destul, DestAccessor ad,
                      RegionStatisticsArray & stats, SRGType srgType)
{
    seededRegionGrowing3D( srcul, shape, as, seedsul, aseeds, 
                           destul, ad, stats, srgType, NeighborCode3DSix());
}

template <class SrcImageIterator, class Diff_type, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray >
inline void
seededRegionGrowing3D(SrcImageIterator srcul, Diff_type shape, SrcAccessor as,
                      SeedImageIterator seedsul, SeedAccessor aseeds,
                      DestImageIterator destul, DestAccessor ad,
                      RegionStatisticsArray & stats)
{
    seededRegionGrowing3D( srcul, shape, as, seedsul, aseeds, destul, ad, 
                           stats, CompleteGrow);
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray, class Neighborhood>
inline void
seededRegionGrowing3D(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                      pair<SeedImageIterator, SeedAccessor> img3,
                      pair<DestImageIterator, DestAccessor> img4,
                      RegionStatisticsArray & stats, 
                      SRGType srgType, Neighborhood n, double max_cost)
{
    seededRegionGrowing3D(img1.first, img1.second, img1.third,
                          img3.first, img3.second,
                          img4.first, img4.second,
                          stats, srgType, n, max_cost);
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray, class Neighborhood>
inline void
seededRegionGrowing3D(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                      pair<SeedImageIterator, SeedAccessor> img3,
                      pair<DestImageIterator, DestAccessor> img4,
                      RegionStatisticsArray & stats, 
                      SRGType srgType, Neighborhood n)
{
    seededRegionGrowing3D(img1.first, img1.second, img1.third,
                          img3.first, img3.second,
                          img4.first, img4.second,
                          stats, srgType, n, NumericTraits<double>::max());
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3D(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                      pair<SeedImageIterator, SeedAccessor> img3,
                      pair<DestImageIterator, DestAccessor> img4,
                      RegionStatisticsArray & stats, SRGType srgType)
{
    seededRegionGrowing3D(img1.first, img1.second, img1.third,
                          img3.first, img3.second,
                          img4.first, img4.second,
                          stats, srgType, NeighborCode3DSix());
}

template <class SrcImageIterator, class Shape, class SrcAccessor,
          class SeedImageIterator, class SeedAccessor,
          class DestImageIterator, class DestAccessor,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3D(triple<SrcImageIterator, Shape, SrcAccessor> img1,
                      pair<SeedImageIterator, SeedAccessor> img3,
                      pair<DestImageIterator, DestAccessor> img4,
                      RegionStatisticsArray & stats)
{
    seededRegionGrowing3D(img1.first, img1.second, img1.third,
                          img3.first, img3.second,
                          img4.first, img4.second,
                          stats);
}

template <class T1, class S1,
          class TS, class AS,
          class T2, class S2,
          class RegionStatisticsArray, class Neighborhood>
inline void
seededRegionGrowing3D(MultiArrayView<3, T1, S1> const & img1,
                      MultiArrayView<3, TS, AS> const & img3,
                      MultiArrayView<3, T2, S2> img4,
                      RegionStatisticsArray & stats, 
                      SRGType srgType, Neighborhood n, double max_cost)
{
    vigra_precondition(img1.shape() == img3.shape(),
        "seededRegionGrowing3D(): shape mismatch between input and output.");
    seededRegionGrowing3D(srcMultiArrayRange(img1),
                          srcMultiArray(img3),
                          destMultiArray(img4),
                          stats, srgType, n, max_cost);
}

template <class T1, class S1,
          class TS, class AS,
          class T2, class S2,
          class RegionStatisticsArray, class Neighborhood>
inline void
seededRegionGrowing3D(MultiArrayView<3, T1, S1> const & img1,
                      MultiArrayView<3, TS, AS> const & img3,
                      MultiArrayView<3, T2, S2> img4,
                      RegionStatisticsArray & stats, 
                      SRGType srgType, Neighborhood n)
{
    vigra_precondition(img1.shape() == img3.shape(),
        "seededRegionGrowing3D(): shape mismatch between input and output.");
    seededRegionGrowing3D(srcMultiArrayRange(img1),
                          srcMultiArray(img3),
                          destMultiArray(img4),
                          stats, srgType, n, NumericTraits<double>::max());
}

template <class T1, class S1,
          class TS, class AS,
          class T2, class S2,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3D(MultiArrayView<3, T1, S1> const & img1,
                      MultiArrayView<3, TS, AS> const & img3,
                      MultiArrayView<3, T2, S2> img4,
                      RegionStatisticsArray & stats, SRGType srgType)
{
    vigra_precondition(img1.shape() == img3.shape(),
        "seededRegionGrowing3D(): shape mismatch between input and output.");
    seededRegionGrowing3D(srcMultiArrayRange(img1),
                          srcMultiArray(img3),
                          destMultiArray(img4),
                          stats, srgType, NeighborCode3DSix());
}

template <class T1, class S1,
          class TS, class AS,
          class T2, class S2,
          class RegionStatisticsArray>
inline void
seededRegionGrowing3D(MultiArrayView<3, T1, S1> const & img1,
                      MultiArrayView<3, TS, AS> const & img3,
                      MultiArrayView<3, T2, S2> img4,
                      RegionStatisticsArray & stats)
{
    vigra_precondition(img1.shape() == img3.shape(),
        "seededRegionGrowing3D(): shape mismatch between input and output.");
    seededRegionGrowing3D(srcMultiArrayRange(img1),
                          srcMultiArray(img3),
                          destMultiArray(img4),
                          stats);
}

} // namespace vigra

#endif // VIGRA_SEEDEDREGIONGROWING_HXX

