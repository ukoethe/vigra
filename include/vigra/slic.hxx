/************************************************************************/
/*                                                                      */
/*    Copyright 2012-2013 by Ullrich Koethe and Thorsten Beier          */
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

#ifndef VIGRA_SLIC_HXX
#define VIGRA_SLIC_HXX

#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <numeric>
#include <iostream>
#include <assert.h>

#include "multi_array.hxx"
#include "multi_iterator_coupled.hxx"
#include "tinyvector.hxx"
#include "labelimage.hxx"
#include "imageiteratoradapter.hxx"
#include "copyimage.hxx"
#include "numerictraits.hxx"
#include "accumulator.hxx"

namespace vigra {

// move slic seeds to a minima of the boundary indicator image (gradient-magnitude)
template <unsigned int N, class T, class Label>
int generateSlicSeedsImpl(
    MultiArrayView<N, T> const &  boundaryIndicatorImage,
    MultiArrayView<N, Label>      seeds,
    int                           seedDist,
    int                           searchRadius)
{
    typedef typename MultiArrayShape<N>::type   Shape;

    seeds.init(0);
    Shape shape(boundaryIndicatorImage.shape());   
    Shape seedShape(floor(shape / double(seedDist)));
    
    int label = 0;
    MultiCoordinateIterator<N> iter(seedShape),
                               end = iter.getEndIterator();
    for(; iter != end; ++iter)
    {
        // define search window around current seed center
        Shape center = (*iter)*seedDist + Shape(seedDist/2);
        Shape startCoord = max(Shape(0), center-Shape(searchRadius));
        Shape endCoord   = min(center+Shape(searchRadius+1), shape);
        
        // find the coordinate of minimum boundary indicator in window
        using namespace acc;
        AccumulatorChain<CoupledArrays<N, T>,
                         Select<WeightArg<1>, Coord<ArgMinWeight> > > a;
        extractFeatures(boundaryIndicatorImage.subarray(startCoord, endCoord), a);

        // add seed at minimum position, if not already set
        Shape minCoord = get<Coord<ArgMinWeight> >(a) + startCoord;
        if(seeds[minCoord] == 0)
            seeds[minCoord] = ++label;
    }
    return label;
}

struct SlicOptions
{
    SlicOptions()
    : iterations(40),
      sizeLimit(4)
    {}
    
    SlicOptions & maxIterations(unsigned int i)
    {
        iterations = i;
        return *this;
    }
    
    SlicOptions & minSize(unsigned int s)
    {
        sizeLimit = s;
        return *this;
    }
    
    unsigned int iterations;
    unsigned int sizeLimit;
};

template <unsigned int N, class T, class Label>
class Slic
{
  public: 
    // 
    typedef MultiArrayView<N, T>                    DataImageType;
    typedef MultiArrayView<N, Label>                LabelImageType; 
    typedef typename DataImageType::difference_type ShapeType;
    typedef typename PromoteTraits<
                   typename NormTraits<T>::NormType,
                   typename NormTraits<MultiArrayIndex>::NormType
             >::Promote                             DistanceType;

    Slic(DataImageType dataImage, 
         LabelImageType labelImage, 
         DistanceType intensityScaling, 
         int maxRadius, 
         SlicOptions const & options = SlicOptions());
    
    int execute();

  private:
    void getWindowLimits(ShapeType const & centerCoord,
                         ShapeType & startCoord, ShapeType & endCoord, const int radius) const;

    DistanceType updateAssigments();
    void updateMeans();
    size_t postProcessing();
    
    typedef MultiArray<N,DistanceType>  DistanceImageType;

    DataImageType                   dataImage_;
    LabelImageType                  labelImage_;
    SlicOptions                     options_;
    MultiArray<N,Label>             tmpLabelImage_;
    DistanceImageType               distance_;
    const ShapeType                 shape_;
    int                             max_radius_;
    DistanceType                    normalization_;
    
    typedef acc::Select<acc::DataArg<1>, acc::LabelArg<2>, acc::Mean, acc::RegionCenter> Statistics;
    typedef acc::AccumulatorChainArray<CoupledArrays<N, T, Label>, Statistics> RegionFeatures;
    RegionFeatures centers_;
};



template <unsigned int N, class T, class Label>
Slic<N, T, Label>::Slic(
    DataImageType         dataImage, 
    LabelImageType        labelImage,
    DistanceType          intensityScaling,
    int                   maxRadius,
    SlicOptions const &   options)
:   dataImage_(dataImage),
    labelImage_(labelImage),
    options_(options),
    tmpLabelImage_(dataImage.shape()),
    distance_(dataImage.shape()),
    shape_(dataImage.shape()),
    max_radius_(maxRadius),
    normalization_(sq(intensityScaling) / sq(max_radius_))
{
    centers_.ignoreLabel(0);
}

template <unsigned int N, class T, class Label>
int Slic<N, T, Label>::execute()
{
    // Do SLIC
    DistanceType err = NumericTraits< DistanceType >::max();
    for(size_t i=0; i<options_.iterations; ++i)
    {
        updateMeans();
        // update which pixels gets assigned to which cluster
        const DistanceType err2 = updateAssigments();

        // convergence?
        if(err2+std::numeric_limits<DistanceType>::epsilon()>=err)
        {
            break;
        }
        err=err2;
    }
    // update assignments bevore postprocessing
    // FIXME: is this really necessary?
    updateMeans();
    err = updateAssigments();
    size_t nBlocked=1;
    while(nBlocked!=0)
    {
        // remove all regions which are smaller than a sizeLimit
        nBlocked = postProcessing();
    }
    // FIXME: use dimension-independent labeling algorithm
    tmpLabelImage_ = labelImage_;
    return labelImage(srcImageRange(tmpLabelImage_), destImage(labelImage_), false);
}

template <unsigned int N, class T, class Label>
typename Slic<N, T, Label>::DistanceType 
Slic<N, T, Label>::updateAssigments()
{
    using namespace acc;
    distance_.init(NumericTraits<DistanceType>::max());
    for(unsigned int c=1; c<=centers_.maxRegionLabel(); ++c)
    {
        if(get<Count>(centers_, c) == 0) // label doesn't exist
            continue;
            
        typedef typename LookupTag<RegionCenter, RegionFeatures>::value_type CenterType;
        CenterType center = get<RegionCenter>(centers_, c);

        // get ROI limits around region center
        ShapeType pixelCoord, startCoord, endCoord;
        getWindowLimits(round(center), startCoord, endCoord, max_radius_);
        center -= startCoord; // need center relative to ROI
        
        // setup iterators for ROI
        typedef typename CoupledArrays<N, T, Label, DistanceType>::IteratorType Iterator;
        Iterator iter = createCoupledIterator(dataImage_, labelImage_, distance_).
                            restrictToSubarray(startCoord, endCoord),
                 end = iter.getEndIterator();
        
        // only pixels within the ROI can be assigned to a cluster
        for(; iter != end; ++iter)
        {
            // compute distance between cluster center and pixel
            DistanceType spatialDist   = squaredNorm(center-iter.point());
            DistanceType colorDist     = squaredNorm(get<Mean>(centers_, c)-iter.get<1>());
            DistanceType dist =  colorDist + normalization_*spatialDist;
            // update label?
            if(dist < iter.get<3>())
            {
                iter.get<2>() = static_cast<Label>(c);
                iter.get<3>() = dist;
            }
        }
    }
    // return total distance
    return distance_.sum<DistanceType>();
}

template <unsigned int N, class T, class Label>
inline void 
Slic<N, T, Label>::getWindowLimits
(
    const ShapeType &  centerCoord,
    ShapeType &        startCoord,
    ShapeType &        endCoord,
    const int radius) const
{
    startCoord = max(ShapeType(0), centerCoord - ShapeType(radius));
    endCoord   = min(shape_, centerCoord + ShapeType(radius+1));
}


template <unsigned int N, class T, class Label>
void 
Slic<N, T, Label>::updateMeans()
{
    // get mean for each cluster
    centers_.reset();
    extractFeatures(dataImage_, labelImage_, centers_);
}

template <unsigned int N, class T, class Label>
size_t 
Slic<N, T, Label>::postProcessing()
{
    // get rid of disjoint regions
    // The original SLIC algorithm is much simpler: whenever it finds a new region, it remembers
    // any of the adjacent existing labels and relabels the new regions with the existing label 
    // in case the new region is too small (i.e. less than 1/4 of the average region size).
    // Labeling in the original algorithm is performed by flood-fill, which makes relabeling
    // easy because we still know the coordinates of the new region and the replacement label
    // when we find out its size.
    // They use 4-neighborhood in 2D, but 10-neighborhood in 3D (indirect nh within slices and 
    // direct nh between slices)
    
    // FIXME: use dimension-independent labeling algorithm
    tmpLabelImage_ = labelImage_;
    size_t numLabels = 1 + labelImage(srcImageRange(tmpLabelImage_), destImage(labelImage_), false);
    
    // FIXME: refactor this in terms of a graph-based algorithm
    std::vector< std::vector< ShapeType> >  regionsPixels(numLabels);
    std::vector< std::vector<Label> >   regionAdjacency(numLabels);
    
    typename LabelImageType::iterator iter = labelImage_.begin(),
                                      end  = iter.getEndIterator();
    for(; iter != end; ++iter)
    {
        Label l1 = *iter;
        ShapeType point = iter.point();
        
        for(int d=0; d<N; ++d)
        {
            if(point[d]+1 < shape_[d])
            {
                point[d] += 1;
                Label l2 = labelImage_[point];
                if(l1!=l2)
                {
                    regionAdjacency[l1].push_back(l2);
                    regionAdjacency[l2].push_back(l1);
                }
                point[d] -= 1;
            }
        }
        regionsPixels[l1].push_back(point);
    }
    // fill region size
    const size_t sizeLimit_ = options_.sizeLimit == 0
                                 ? (size_t)std::pow((double)labelImage_.size() / numLabels, 1.0 / N)
                                 : options_.sizeLimit;
    std::vector<size_t> cSize(numLabels);
    for(size_t c=1;c<numLabels;++c)
    {
        cSize[c]=regionsPixels[c].size();
    }
    size_t numChanges=0,blocked=0;
    std::vector<bool> merged(numLabels,false);
    //search for regions which are smaller than size limit
    for(size_t c=1;c<numLabels;++c)
    {
        const size_t regionSize=regionsPixels[c].size();
        // a region is to small?
        if(regionSize<sizeLimit_ )
        {
            // the region is not a merged one (get ride of this check?)
            if(merged[c]==false)
            {
                size_t mergeWith=0,maxSize=0;
                bool found=false; // search for a region to merge the smaller region with
                for(size_t i=0;i<regionAdjacency[c].size();++i)
                {
                    const size_t c2=regionAdjacency[c][i];
                    assert(c2!=c);
                    const size_t size2=regionsPixels[c2].size();
                    if(size2 >=maxSize && merged[c2]==false)
                    {
                        found=true;
                        mergeWith=c2;
                        maxSize=size2;
                    }
                }
                // found a region to merge with?
                if(found)
                {
                    cSize[mergeWith]+=regionSize;
                    // is there size > sizeLimit
                    if(cSize[mergeWith]<sizeLimit_)
                    {
                        merged[mergeWith]=true;
                        ++blocked;
                    }
                    merged[c]=true;
                    assert(c!=mergeWith);
                    ++numChanges;
                    for(size_t p=0;p<regionsPixels[c].size();++p)
                    {
                        ////std::cout<<"change stuff\n";
                        assert(labelImage_(regionsPixels[c][p][0],regionsPixels[c][p][1])!=mergeWith);
                        labelImage_(regionsPixels[c][p][0],regionsPixels[c][p][1])=mergeWith;
                    }
                }
                // did not find one.
                else
                {
                    ++blocked;
                }
            }
            // did not find one.
            else
            {
                ++blocked;
            }
        }
    }
    if(numChanges==0)
        return 0;
    return blocked;
}

template <unsigned int N, class T, class Label>
inline int 
generateSlicSeeds(
    MultiArrayView<N, T> const &  boundaryIndicatorImage,
    MultiArrayView<N, Label>      seeds,
    int                           seedDistance,
    int                           searchRadius=1)
{
    return generateSlicSeedsImpl( boundaryIndicatorImage, seeds, seedDistance, searchRadius);
}


template <unsigned int N, class T, class Label>
inline int 
slicSuperpixels(
    const MultiArrayView<N, T> &    dataImage,
    MultiArrayView<N, Label>        labelImage,
    typename Slic<N, T, Label>::DistanceType intensityScaling,
    int                             seedDistance, 
    const SlicOptions &             parameter = SlicOptions())
{
    return Slic<N, T, Label>(dataImage, labelImage, intensityScaling, seedDistance, parameter).execute();
}

} // namespace vigra

#endif // VIGRA_SLIC_HXX
