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

#include "multi_array.hxx"
#include "multi_convolution.hxx"
#include "multi_labeling.hxx"
#include "numerictraits.hxx"
#include "accumulator.hxx"
#include "array_vector.hxx"

namespace vigra {

/** \addtogroup SeededRegionGrowing
*/
//@{

/********************************************************/
/*                                                      */
/*                  generateSlicSeeds                   */
/*                                                      */
/********************************************************/

/** \brief Generate seeds for SLIC superpixel computation in arbitrary dimensions.

    The source array \a src must be a scalar boundary indicator such as the gradient 
    magnitude. Seeds are initially placed on a regular Cartesian grid with spacing
    \a seedDist und then moved to the point with smallest boundary indicator within
    a search region of radius \a searchRadius around the initial position. The resulting
    points are then marked in the output array \a seeds by consecutive labels.
    
    The function returns the number of selected seeds, which equals the largest seed label 
    because labeling starts at 1.

    <b> Declaration:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        template <unsigned int N, class T, class S1,
                                  class Label, class S2>
        unsigned int 
        generateSlicSeeds(MultiArrayView<N, T, S1> const & src,
                          MultiArrayView<N, Label, S2>     seeds,
                          unsigned int                     seedDist,
                          unsigned int                     searchRadius = 1);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/slic.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, RGBValue<float> > src(Shape2(w, h));
    ... // fill src image
    
    // transform image to Lab color space
    transformImage(srcImageRange(src), destImage(src), RGBPrime2LabFunctor<float>());
    
    // compute image gradient magnitude at scale 1.0 as a boundary indicator
    MultiArray<2, float> grad(src.shape());
    gaussianGradientMagnitude(srcImageRange(src), destImage(grad), 1.0);
    
    MultiArray<2, unsigned int>  seeds(src.shape());
    int seedDistance = 15;
    
    // place seeds on a grid with distance 15, but then move it to the lowest gradient
    // poistion in a 3x3 window
    generateSlicSeeds(grad, seeds, seedDistance);
    \endcode

    For more details and examples see slicSuperpixels().
*/
template <unsigned int N, class T, class S1,
                          class Label, class S2>
unsigned int 
generateSlicSeeds(MultiArrayView<N, T, S1> const & boundaryIndicatorImage,
                  MultiArrayView<N, Label, S2>     seeds,
                  unsigned int                     seedDist,
                  unsigned int                     searchRadius = 1)
{
    typedef typename MultiArrayShape<N>::type   Shape;

    seeds.init(0);
    Shape shape(boundaryIndicatorImage.shape()),
          seedShape(floor(shape / double(seedDist))),
          offset((shape - (seedShape - Shape(1))*seedDist) / 2);
    
    unsigned int label = 0;
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

/** \brief Options object for slicSuperpixels().

    <b> Usage:</b>

    see slicSuperpixels() for detailed examples.
*/
struct SlicOptions
{
        /** \brief Create options object with default settings.

            Defaults are: perform 10 iterations, determine a size limit for superpixels automatically.
        */
    SlicOptions()
    : iter(10),
      sizeLimit(0)
    {}
    
        /** \brief Number of iterations.

            Default: 10
        */
    SlicOptions & iterations(unsigned int i)
    {
        iter = i;
        return *this;
    }
    
        /** \brief Minimum superpixel size.
        
            If you set this to 1, no size filtering will be performed.

            Default: 0 (determine size limit automatically as <tt>average size / 4</tt>)
        */
    SlicOptions & minSize(unsigned int s)
    {
        sizeLimit = s;
        return *this;
    }
    
    unsigned int iter;
    unsigned int sizeLimit;
};

namespace detail {

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
    void updateAssigments();
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
    for(size_t i=0; i<options_.iter; ++i)
    {
        // update mean for each cluster
        clusters_.reset();
        extractFeatures(dataImage_, labelImage_, clusters_);
        
        // update which pixels get assigned to which cluster
        updateAssigments();
    }

    return postProcessing();
}

template <unsigned int N, class T, class Label>
void
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

    // make sure that all regions exceed the sizeLimit
    for (graph_scanner node(graph); node != lemon::INVALID; ++node) 
    {
        Label label = labelImage_[*node];
        
        if(regions[label] > 0)
            continue;   // already processed
            
        regions[label] = label;
        
        if(get<Count>(sizes, label) < sizeLimit)
        {
            // region is too small => merge into an existing neighbor
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

} // namespace detail


/** \brief Compute SLIC superpixels in arbitrary dimensions.

    This function implements the algorithm described in 
    
    R. Achanta et al.: "<i>SLIC Superpixels Compared to State-of-the-Art 
    Superpixel Methods</i>", IEEE Trans. Patt. Analysis Mach. Intell. 34(11):2274-2281, 2012
    
    The value type <tt>T</tt> of the source array \a src must provide the necessary functionality 
    to compute averages and squared distances (i.e. it must fulfill the requirements of a 
    \ref LinearSpace and support squaredNorm()). This is true for all scalar types as well as
    \ref vigra::TinyVector and \ref vigra::RGBValue. The output array \a labels will be filled
    with labels designating membership of each point in one of the superpixel regions.
    
    The output array can optionally contain seeds (which will be overwritten by the output) 
    to give you full control over seed placement. If \a labels is empty, seeds will be created
    automatically by an internal call to generateSlicSeeds(). 
    
    The parameter \a seedDistance specifies the radius of the window around each seed (or, more
    precisely, around the present regions centers) where the algorithm looks for potential members 
    of the corresponding superpixel. It thus places an upper limit on the superpixel size. When seeds
    are computed automatically, this parameter also determines the grid spacing for seed placement.
    
    The parameter \a intensityScaling is used to normalize (i.e. divide) the color/intensity difference 
    before it is compared with the spatial distance. This corresponds to parameter <i>m</i> in equation
    (2) of the paper.
    
    The options object can be used to specify the number of iterations (<tt>SlicOptions::iterations()</tt>)
    and an explicit minimal superpixel size (<tt>SlicOptions::minSize()</tt>). By default, the algorithm 
    merges all regions that are smaller than 1/4 the average superpixel size.
    
    The function returns the number of superpixels, which equals the largest label 
    because labeling starts at 1.

    <b> Declaration:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class DistanceType>
        unsigned int 
        slicSuperpixels(MultiArrayView<N, T, S1> const &  src,
                        MultiArrayView<N, Label, S2>      labels,
                        DistanceType                      intensityScaling,
                        unsigned int                      seedDistance, 
                        SlicOptions const &               options = SlicOptions());
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/slic.hxx\><br>
    Namespace: vigra

    \code
    MultiArray<2, RGBValue<float> > src(Shape2(w, h));
    ... // fill src image
    
    // transform image to Lab color space
    transformMultiArray(srcMultiArrayRange(src), destMultiArray(src), RGBPrime2LabFunctor<float>());    
    
    MultiArray<2, unsigned int>  labels(src.shape());
    int seedDistance = 15;
    double intensityScaling = 20.0;
    
    // compute seeds automatically, perform 40 iterations, scaling intensity differences
    // down to 1/20 before comparing with spatial distances
    slicSuperpixels(src, labels, intensityScaling, seedDistance, SlicOptions().iterations(40));
    \endcode
    
    This works for arbitrary-dimensional arrays.
*/
doxygen_overloaded_function(template <...> unsigned int slicSuperpixels)

template <unsigned int N, class T, class S1,
                          class Label, class S2,
          class DistanceType>
unsigned int 
slicSuperpixels(MultiArrayView<N, T, S1> const &  src,
                MultiArrayView<N, Label, S2>      labels,
                DistanceType                      intensityScaling,
                unsigned int                      seedDistance, 
                SlicOptions const &               options = SlicOptions())
{
    if(!labels.any())
    {
        typedef typename NormTraits<T>::NormType TmpType;
        MultiArray<N, TmpType> grad(src.shape());
        gaussianGradientMagnitude(src, grad, 1.0);
        generateSlicSeeds(grad, labels, seedDistance);
    }
    return detail::Slic<N, T, Label>(src, labels, intensityScaling, seedDistance, options).execute();
}

//@}

} // namespace vigra

#endif // VIGRA_SLIC_HXX
