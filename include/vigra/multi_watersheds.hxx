/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2013 by Ullrich Koethe                  */
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

#ifndef VIGRA_MULTI_WATERSHEDS_HXX
#define VIGRA_MULTI_WATERSHEDS_HXX

#include <functional>
#include "mathutil.hxx"
#include "multi_array.hxx"
#include "multi_math.hxx"
#include "multi_gridgraph.hxx"
#include "multi_localminmax.hxx"
#include "multi_labeling.hxx"
#include "watersheds.hxx"
#include "bucket_queue.hxx"
#include "union_find.hxx"

namespace vigra {

namespace lemon_graph {

namespace graph_detail {

template <class Graph, class T1Map, class T2Map>
void
prepareWatersheds(Graph const & g, 
                  T1Map const & data,
                  T2Map & lowestNeighborIndex)
{
    typedef typename Graph::NodeIt    graph_scanner;
    typedef typename Graph::OutArcIt  neighbor_iterator;

    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        typename T1Map::value_type lowestValue  = data[*node];
        typename T2Map::value_type lowestIndex  = -1;

        for(neighbor_iterator arc(g, node); arc != INVALID; ++arc) 
        {
            if(data[g.target(*arc)] <= lowestValue)
            {
                lowestValue = data[g.target(*arc)];
                lowestIndex = arc.neighborIndex();
            }
        }
        lowestNeighborIndex[*node] = lowestIndex;
    }
}

// template <class Graph, class T1Map, class T2Map>
// typename T2Map::value_type
// unionFindWatersheds(Graph const & g, 
                    // T1Map const & lowestNeighborIndex,
                    // T2Map & labels)
// {
    // typedef typename Graph::NodeIt        graph_scanner;
    // typedef typename Graph::OutBackArcIt  neighbor_iterator;
    // typedef typename T2Map::value_type    LabelType;

    // vigra::detail::UnionFindArray<LabelType>  regions;

    // // pass 1: find connected components
    // for (graph_scanner node(g); node != INVALID; ++node) 
    // {
        // // define tentative label for current node
        // LabelType currentLabel = regions.nextFreeLabel();
        
        // for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        // {
            // // merge regions if current target is center's lowest neighbor or vice versa
            // if(lowestNeighborIndex[*node] == arc.neighborIndex() || 
               // lowestNeighborIndex[g.target(*arc)] == g.oppositeIndex(arc.neighborIndex()))
            // {
                // LabelType neighborLabel = regions[labels[g.target(*arc)]];
                // currentLabel = regions.makeUnion(neighborLabel, currentLabel);
            // }
        // }
        // // set label of current node
        // labels[*node] = regions.finalizeLabel(currentLabel);
    // }
    
    // LabelType count = regions.makeContiguous();

    // // pass 2: make component labels contiguous
    // for (graph_scanner node(g); node != INVALID; ++node) 
    // {
        // labels[*node] = regions[labels[*node]];
    // }
    // return count;
// }

template <class Graph, class T1Map, class T2Map, class T3Map>
typename T2Map::value_type
unionFindWatersheds(Graph const & g,
                    T1Map const & data, 
                    T2Map const & lowestNeighborIndex,
                    T3Map & labels)
{
    typedef typename Graph::Node          Node;
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T1Map::value_type    DataType;
    typedef typename T3Map::value_type    LabelType;

    vigra::detail::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        // define tentative label for current node
        LabelType currentLabel = regions.nextFreeLabel();
        Node plateauNeighbor(INVALID);
        
        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            // merge regions if current target is center's lowest neighbor or vice versa
            if(lowestNeighborIndex[*node] == arc.neighborIndex() || 
               lowestNeighborIndex[g.target(*arc)] == g.oppositeIndex(arc.neighborIndex()))
            {
                if(data[*node] == data[g.target(*arc)])
                    plateauNeighbor = g.target(*arc);
                LabelType neighborLabel = regions[labels[g.target(*arc)]];
                currentLabel = regions.makeUnion(neighborLabel, currentLabel);
            }
        }
        
        if(plateauNeighbor != INVALID)
        {
            // we are on a plateau => link all plateau points
            for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
            {
                if(data[*node] == data[g.target(*arc)])
                {
                    LabelType neighborLabel = regions[labels[g.target(*arc)]];
                    currentLabel = regions.makeUnion(neighborLabel, currentLabel);
                }
            }
        }
        
        // set label of current node
        labels[*node] = regions.finalizeLabel(currentLabel);
    }
    
    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        labels[*node] = regions[labels[*node]];
    }
    return count;
}

template <class Graph, class T1Map, class T2Map>
typename T2Map::value_type
generateWatershedSeeds(Graph const & g, 
                       T1Map const & data,
                       T2Map & seeds,
                       SeedOptions const & options = SeedOptions())
{
    typedef typename T1Map::value_type DataType;
    typedef unsigned char MarkerType;
    
    typename Graph::template NodeMap<MarkerType>  minima(g);
    
    if(options.mini == SeedOptions::LevelSets)
    {
        vigra_precondition(options.thresholdIsValid<DataType>(),
            "generateWatershedSeeds(): SeedOptions.levelSets() must be specified with threshold.");
    
        using namespace multi_math;
        minima = data <= DataType(options.thresh);
    }
    else
    {
        DataType threshold = options.thresholdIsValid<DataType>()
                                ? options.thresh
                                : NumericTraits<DataType>::max();
        
        if(options.mini == SeedOptions::ExtendedMinima)
            extendedLocalMinMaxGraph(g, data, minima, MarkerType(1), threshold, 
                                     std::less<DataType>(), std::equal_to<DataType>(), true);
        else
            localMinMaxGraph(g, data, minima, MarkerType(1), threshold, 
                             std::less<DataType>(), true);
    }
    return labelGraphWithBackground(g, minima, seeds, MarkerType(0), std::equal_to<MarkerType>());
}


template <class Graph, class T1Map, class T2Map>
typename T2Map::value_type 
seededWatersheds(Graph const & g, 
                 T1Map const & data,
                 T2Map & labels,
                 WatershedOptions const & options)
{
    typedef typename Graph::Node        Node;
    typedef typename Graph::NodeIt      graph_scanner;
    typedef typename Graph::OutArcIt    neighbor_iterator;
    typedef typename T1Map::value_type  CostType;
    typedef typename T2Map::value_type  LabelType;

    PriorityQueue<Node, CostType, true> pqueue;
    
    bool keepContours = ((options.terminate & KeepContours) != 0);
    LabelType maxRegionLabel = 0;
    
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        LabelType label = labels[*node];
        if(label != 0)
        {
            if(maxRegionLabel < label)
                maxRegionLabel = label;
                
            for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
            {
                if(labels[g.target(*arc)] == 0)
                {
                    // register all seeds that have an unlabeled neighbor
                    if(label == options.biased_label)
                        pqueue.push(*node, data[*node] * options.bias);
                    else
                        pqueue.push(*node, data[*node]);
                    break;
                }
            }
        }
    }
    
    LabelType contourLabel = maxRegionLabel + 1;  // temporary contour label
    
    // perform region growing
    while(!pqueue.empty())
    {
        Node node = pqueue.top();
        CostType cost = pqueue.topPriority();
        pqueue.pop();
        
        if((options.terminate & StopAtThreshold) && (cost > options.max_cost))
            break;

        LabelType label = labels[node];
        
        if(label == contourLabel)
            continue;

        // Put the unlabeled neighbors in the priority queue.
        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            LabelType neighborLabel = labels[g.target(*arc)];
            if(neighborLabel == 0)
            {
                labels[g.target(*arc)] = label;
                CostType priority = (label == options.biased_label)
                                       ? data[g.target(*arc)] * options.bias
                                       : data[g.target(*arc)];
                if(priority < cost)
                    priority = cost;
                pqueue.push(g.target(*arc), priority);
            }
            else if(keepContours && (label != neighborLabel) && (neighborLabel != contourLabel))
            {
                // Check if the present point is adjacent to more than one region
                // and if so, mark it as contour.
                CostType priority = (neighborLabel == options.biased_label)
                                       ? data[g.target(*arc)] * options.bias
                                       : data[g.target(*arc)];
                if(cost < priority) // neighbor not yet processed
                    labels[g.target(*arc)] = contourLabel;
            }
        }
    }
    
    if(keepContours)
    {
        // Replace the temporary contour label with label 0.
        typename T2Map::iterator k   = labels.begin(),
                                 end = labels.end();
        for(; k != end; ++k)
            if(*k == contourLabel)
                *k = 0;
    }
    
    return maxRegionLabel;
}

} // namespace graph_detail

template <class Graph, class T1Map, class T2Map>
typename T2Map::value_type 
watershedsGraph(Graph const & g, 
                T1Map const & data,
                T2Map & labels,
                WatershedOptions const & options)
{
    typedef typename T2Map::value_type LabelType;
    
    if(options.method == WatershedOptions::UnionFind)
    {
        vigra_precondition(g.maxDegree() <= NumericTraits<unsigned short>::max(),
            "watershedsGraph(): cannot handle nodes with degree > 65535.");
            
        typename Graph::template NodeMap<unsigned short>  lowestNeighborIndex(g);
        
        graph_detail::prepareWatersheds(g, data, lowestNeighborIndex);
        return graph_detail::unionFindWatersheds(g, data, lowestNeighborIndex, labels);
    }
    else if(options.method == WatershedOptions::RegionGrowing)
    {
        SeedOptions seed_options;
        
        // check if the user has explicitly requested seed computation
        if(options.seed_options.mini != SeedOptions::Unspecified)
        {
            seed_options = options.seed_options;
        }
        else
        {
            // otherwise, don't compute seeds if 'labels' already contains them 
            using namespace multi_math;
            if(any(labels > LabelType(0)))
                seed_options.mini = SeedOptions::Unspecified;
        }

        if(seed_options.mini != SeedOptions::Unspecified)
        {
            graph_detail::generateWatershedSeeds(g, data, labels, seed_options);
        }
        
        return graph_detail::seededWatersheds(g, data, labels, options);
    }
    else
    {
        vigra_precondition(false,
           "watershedsGraph(): invalid method in watershed options.");
        return 0;
    }
}


} // namespace lemon_graph

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

    <b>\#include</b> \<vigra/watersheds.hxx\><br>
    Namespace: vigra

    For detailed examples see watershedsRegionGrowing().
*/
doxygen_overloaded_function(template <...> unsigned int generateWatershedSeeds)

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
generateWatershedSeeds(MultiArrayView<N, T, S1> const & data,
                       MultiArrayView<N, Label, S2> seeds,
                       NeighborhoodType neighborhood = IndirectNeighborhood,
                       SeedOptions const & options = SeedOptions())
{
    vigra_precondition(data.shape() == seeds.shape(),
        "generateWatershedSeeds(): Shape mismatch between input and output.");
    
    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    return lemon_graph::graph_detail::generateWatershedSeeds(graph, data, seeds, options);
}


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
doxygen_overloaded_function(template <...> unsigned int watershedsMultiArray)

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
watershedsMultiArray(MultiArrayView<N, T, S1> const & data,
                     MultiArrayView<N, Label, S2> labels,  // may also hold input seeds
                     NeighborhoodType neighborhood = IndirectNeighborhood,
                     WatershedOptions const & options = WatershedOptions())
{
    vigra_precondition(data.shape() == labels.shape(),
        "watershedsMultiArray(): Shape mismatch between input and output.");
    
    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    return lemon_graph::watershedsGraph(graph, data, labels, options);
}

//@}

} // namespace vigra

#endif // VIGRA_MULTI_WATERSHEDS_HXX
