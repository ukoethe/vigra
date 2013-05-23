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
#include "multi_labeling.hxx"
#include "numerictraits.hxx"
#include "accumulator.hxx"

namespace vigra {

// move slic seeds to the smallest boundary indicator within search resius
template <unsigned int N, class T, class Label>
int generateSlicSeedsImpl(
    MultiArrayView<N, T> const &  boundaryIndicatorImage,
    MultiArrayView<N, Label>      seeds,
    int                           seedDist,
    int                           searchRadius)
{
    typedef typename MultiArrayShape<N>::type   Shape;

    seeds.init(0);
    Shape shape(boundaryIndicatorImage.shape()),
          seedShape(floor(shape / double(seedDist))),
          offset((shape - (seedShape - Shape(1))*seedDist) / 2);
    
    int label = 0;
    MultiCoordinateIterator<N> iter(seedShape),
                               end = iter.getEndIterator();
    for(; iter != end; ++iter)
    {
        // define search window around current seed center
        Shape center = (*iter)*seedDist + offset;
        Shape startCoord = max(Shape(0), center-Shape(searchRadius));
        Shape endCoord   = min(center+Shape(searchRadius+1), shape);
        
        // find the coordinate of minimum boundary indicator in window
        using namespace acc;
        AccumulatorChain<CoupledArrays<N, T>,
                         Select<WeightArg<1>, Coord<ArgMinWeight> > > a;
        extractFeatures(boundaryIndicatorImage.subarray(startCoord, endCoord), a);

        // add seed at minimum position, if not already occupied
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
    
    unsigned int execute();

  private:
    DistanceType updateAssigments();
    unsigned int postProcessing();
    
    typedef MultiArray<N,DistanceType>  DistanceImageType;

    ShapeType                       shape_;
    DataImageType                   dataImage_;
    LabelImageType                  labelImage_;
    DistanceImageType               distance_;
    int                             max_radius_;
    DistanceType                    normalization_;
    SlicOptions                     options_;
    
    typedef acc::Select<acc::DataArg<1>, acc::LabelArg<2>, acc::Mean, acc::RegionCenter> Statistics;
    typedef acc::AccumulatorChainArray<CoupledArrays<N, T, Label>, Statistics> RegionFeatures;
    RegionFeatures clusters_;
};



template <unsigned int N, class T, class Label>
Slic<N, T, Label>::Slic(
    DataImageType         dataImage, 
    LabelImageType        labelImage,
    DistanceType          intensityScaling,
    int                   maxRadius,
    SlicOptions const &   options)
:   shape_(dataImage.shape()),
    dataImage_(dataImage),
    labelImage_(labelImage),
    distance_(shape_),
    max_radius_(maxRadius),
    normalization_(sq(intensityScaling) / sq(max_radius_)),
    options_(options)
{
    clusters_.ignoreLabel(0);
}

template <unsigned int N, class T, class Label>
unsigned int Slic<N, T, Label>::execute()
{
    // Do SLIC
    DistanceType err = NumericTraits< DistanceType >::max();
    for(size_t i=0; i<options_.iterations; ++i)
    {
        // update mean for each cluster
        clusters_.reset();
        extractFeatures(dataImage_, labelImage_, clusters_);
        
        // update which pixels get assigned to which cluster
        const DistanceType err2 = updateAssigments();

        // convergence?
        if(err2+std::numeric_limits<DistanceType>::epsilon()>=err)
        {
            break;
        }
        err=err2;
    }

    return postProcessing();
}

template <unsigned int N, class T, class Label>
typename Slic<N, T, Label>::DistanceType 
Slic<N, T, Label>::updateAssigments()
{
    using namespace acc;
    distance_.init(NumericTraits<DistanceType>::max());
    for(unsigned int c=1; c<=clusters_.maxRegionLabel(); ++c)
    {
        if(get<Count>(clusters_, c) == 0) // label doesn't exist
            continue;
            
        typedef typename LookupTag<RegionCenter, RegionFeatures>::value_type CenterType;
        CenterType center = get<RegionCenter>(clusters_, c);

        // get ROI limits around region center
        ShapeType pixelCenter(round(center)), 
                  startCoord(max(ShapeType(0), pixelCenter - ShapeType(max_radius_))), 
                  endCoord(min(shape_, pixelCenter + ShapeType(max_radius_+1)));
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
            DistanceType colorDist     = squaredNorm(get<Mean>(clusters_, c)-iter.get<1>());
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
unsigned int 
Slic<N, T, Label>::postProcessing()
{
    // get rid of regions below a size limit
    MultiArray<N,Label> tmpLabelImage(labelImage_);
    unsigned int maxLabel = labelMultiArray(tmpLabelImage, labelImage_, DirectNeighborhood);
    
    unsigned int sizeLimit = options_.sizeLimit == 0
                                 ? (unsigned int)(0.25 * labelImage_.size() / maxLabel)
                                 : options_.sizeLimit;
    if(sizeLimit == 1)
        return maxLabel;
        
    // determine region size
    using namespace acc;
    AccumulatorChainArray<CoupledArrays<N, Label>, Select<LabelArg<1>, Count> > sizes;
    extractFeatures(labelImage_, sizes);
        
    typedef GridGraph<N, undirected_tag> Graph;
    Graph graph(labelImage_.shape(), DirectNeighborhood);
    
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;

    ArrayVector<Label> regions(maxLabel+1);

    for (graph_scanner node(graph); node != lemon::INVALID; ++node) 
    {
        Label label = labelImage_[*node];
        
        if(regions[label] > 0)
            continue;   // already processed
            
        regions[label] = label;
        
        // merge region into neighborng one if too small
        if(get<Count>(sizes, label) < sizeLimit)
        {
            for (neighbor_iterator arc(graph, node); arc != lemon::INVALID; ++arc)
            {
                regions[label] = regions[labelImage_[graph.target(*arc)]];
                break;
            }
        }
    }
    
    // make labels contiguous after possible merging
    maxLabel = 0; 
    for(unsigned int i=1; i<=maxLabel; ++i)
    {
        if(regions[i] == i)
        {
                regions[i] = (Label)++maxLabel;
        }
        else
        {
                regions[i] = regions[regions[i]]; 
        }
    }

    // update labels
    for (graph_scanner node(graph); node != lemon::INVALID; ++node) 
    {
        labelImage_[*node] = regions[labelImage_[*node]];
    }
    
    return maxLabel;
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


template <unsigned int N, class T, class Label, class DistanceType>
inline int 
slicSuperpixels(
    MultiArrayView<N, T> const &  dataImage,
    MultiArrayView<N, Label>      labelImage,
    DistanceType                  intensityScaling,
    int                           seedDistance, 
    SlicOptions const &           options = SlicOptions())
{
    return Slic<N, T, Label>(dataImage, labelImage, intensityScaling, seedDistance, options).execute();
}

} // namespace vigra

#endif // VIGRA_SLIC_HXX
