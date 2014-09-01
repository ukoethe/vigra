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

/** \addtogroup SeededRegionGrowing
*/
//@{
namespace lemon_graph {

namespace graph_detail {

    // select the neighbor ID for union-find watersheds
    // standard behavior: use global node ID
template <class Graph>
struct NeighborIndexFunctor
{
    typedef typename Graph::index_type index_type;
    
    template <class NodeIter, class ArcIter>
    static index_type get(Graph const & g, NodeIter const &, ArcIter const & a)
    {
        return g.id(g.target(*a));
    }
    
    template <class NodeIter, class ArcIter>
    static index_type getOpposite(Graph const & g, NodeIter const & n, ArcIter const &)
    {
        return g.id(*n);
    }
};

    // select the neighbor ID for union-find watersheds
    // GridGraph optimization: use local neighbor index (needs only 1/4 of the memory)
template<unsigned int N, class DirectedTag>
struct NeighborIndexFunctor<GridGraph<N, DirectedTag> >
{
    typedef GridGraph<N, DirectedTag> Graph;
    typedef UInt16 index_type;
    
    template <class NodeIter, class ArcIter>
    static index_type get(Graph const & g, NodeIter const &, ArcIter const & a)
    {
        return a.neighborIndex();
    }
    
    template <class NodeIter, class ArcIter>
    static index_type getOpposite(Graph const & g, NodeIter const &, ArcIter const & a)
    {
        return g.oppositeIndex(a.neighborIndex());
    }
};

template <class Graph, class T1Map, class T2Map>
void
prepareWatersheds(Graph const & g, 
                  T1Map const & data,
                  T2Map & lowestNeighborIndex)
{
    typedef typename Graph::NodeIt    graph_scanner;
    typedef typename Graph::OutArcIt  neighbor_iterator;
    typedef NeighborIndexFunctor<Graph> IndexFunctor;

    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        typename T1Map::value_type lowestValue  = data[*node];
        typename T2Map::value_type lowestIndex  = NumericTraits<typename T2Map::value_type>::max();

        for(neighbor_iterator arc(g, *node); arc != INVALID; ++arc) 
        {
            if(data[g.target(*arc)] <= lowestValue)
            {
                lowestValue = data[g.target(*arc)];
                lowestIndex = IndexFunctor::get(g, node, arc);
            }
        }
        lowestNeighborIndex[*node] = lowestIndex;
    }
}


template <class Graph, class T1Map, class T2Map, class T3Map>
typename T2Map::value_type
unionFindWatersheds(Graph const & g,
                    T1Map const & data, 
                    T2Map const & lowestNeighborIndex,
                    T3Map & labels)
{
    typedef typename Graph::NodeIt        graph_scanner;
    typedef typename Graph::OutBackArcIt  neighbor_iterator;
    typedef typename T3Map::value_type    LabelType;
    typedef NeighborIndexFunctor<Graph> IndexFunctor;

    vigra::UnionFindArray<LabelType>  regions;

    // pass 1: find connected components
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        // define tentative label for current node
        LabelType currentIndex = regions.nextFreeIndex();
        bool hasPlateauNeighbor = false;
        
        for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
        {
            // merge regions if current target is center's lowest neighbor or vice versa
            if(lowestNeighborIndex[*node] == IndexFunctor::get(g, node, arc) || 
               lowestNeighborIndex[g.target(*arc)] == IndexFunctor::getOpposite(g, node, arc))
            {
                if(data[*node] == data[g.target(*arc)])
                    hasPlateauNeighbor = true;
                currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
            }
        }
        
        if(hasPlateauNeighbor)
        {
            // we are on a plateau => link all plateau points
            for (neighbor_iterator arc(g, node); arc != INVALID; ++arc)
            {
                if(data[*node] == data[g.target(*arc)])
                {
                    currentIndex = regions.makeUnion(labels[g.target(*arc)], currentIndex);
                }
            }
        }
        
        // set label of current node
        labels[*node] = regions.finalizeIndex(currentIndex);
    }
    
    LabelType count = regions.makeContiguous();

    // pass 2: make component labels contiguous
    for (graph_scanner node(g); node != INVALID; ++node) 
    {
        labels[*node] = regions.findLabel(labels[*node]);
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
        for(typename Graph::NodeIt iter(g);iter!=lemon::INVALID;++iter){
            minima[*iter]= data[*iter] <= DataType(options.thresh);
        }
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
                
            for (neighbor_iterator arc(g, *node); arc != INVALID; ++arc)
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
                // The present neighbor is adjacent to more than one region
                // => mark it as contour.
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
        ///typename T2Map::iterator k   = labels.begin(),
        ///                         end = labels.end();
        ///for(; k != end; ++k)
        ///    if(*k == contourLabel)
        ///        *k = 0;

        for(typename Graph::NodeIt iter(g);iter!=lemon::INVALID;++iter){
            if(labels[*iter]==contourLabel)
                labels[*iter]=0;
        }
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
    if(options.method == WatershedOptions::UnionFind)
    {
        typedef typename graph_detail::NeighborIndexFunctor<Graph>::index_type index_type;
        
        vigra_precondition(g.maxDegree() <= NumericTraits<index_type>::max(),
            "watershedsGraph(): cannot handle nodes with degree > 65535.");
            
        typename Graph::template NodeMap<index_type>  lowestNeighborIndex(g);
        
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
            if(labels.any())
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

    // documentation is in watersheds.hxx
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


/** \brief Watershed segmentation of an arbitrary-dimensional array.

    This function implements variants of the watershed algorithms
    described in

    [1] L. Vincent and P. Soille: <em>"Watersheds in digital spaces: An efficient algorithm
    based on immersion simulations"</em>, IEEE Trans. Patt. Analysis Mach. Intell. 13(6):583-598, 1991

    [2] J. Roerdink, R. Meijster: <em>"The watershed transform: definitions, algorithms,
    and parallelization strategies"</em>, Fundamenta Informaticae, 41:187-228, 2000

    The source array \a data is a boundary indicator such as the gaussianGradientMagnitude()
    or the trace of the \ref boundaryTensor(), and the destination \a labels is a label array
    designating membership of each point in one of the regions found. Plateaus in the boundary
    indicator are handled via simple tie breaking strategies. Argument \a neighborhood 
    specifies the connectivity between points and can be <tt>DirectNeighborhood</tt> (meaning 
    4-neighborhood in 2D and 6-neighborhood in 3D, default) or <tt>IndirectNeighborhood</tt> 
    (meaning 8-neighborhood in 2D and 26-neighborhood in 3D).
    
    The watershed variant to be applied can be selected in the \ref WatershedOptions
    object: When you call <tt>WatershedOptions::regionGrowing()</tt> (default), the flooding
    algorithm from [1] is used. Alternatively, <tt>WatershedOptions::unionFind()</tt> uses
    the scan-line algorithm 4.7 from [2]. The latter is faster, but does not support any options 
    (if you pass options nonetheless, they are silently ignored).
    
    The region growing algorithm needs a seed for each region. Seeds can either be provided in
    the destination array \a labels (which will then be overwritten with the result) or computed
    automatically by an internal call to generateWatershedSeeds(). In the former case you have 
    full control over seed placement, while the latter is more convenient. Automatic seed 
    computation is performed when you provide seeding options via <tt>WatershedOptions::seedOptions()</tt> 
    or when the array \a labels is empty (all zeros), in which case default seeding options 
    are chosen. The destination image should be zero-initialized for automatic seed computation.
    
    Further options to be specified via \ref WatershedOptions are:
    
    <ul>
    <li> <tt>keepContours()</tt>: Whether to keep a 1-pixel-wide contour (with label 0) between 
         regions (otherwise, a complete grow is performed, i.e. all pixels are assigned to a region).
    <li> <tt>stopAtThreshold()</tt>: Whether to stop growing when the boundaryness exceeds a threshold 
         (remaining pixels keep label 0).
    <li> <tt>biasLabel()</tt>: Whether one region (label) is to be preferred or discouraged by biasing its cost 
         with a given factor (smaller than 1 for preference, larger than 1 for discouragement).
    </ul>
    
    The option <tt>turboAlgorithm()</tt> is implied by method <tt>regionGrowing()</tt> (this is
    in contrast to watershedsRegionGrowing(), which supports an additional algorithm in 2D only).

    watershedsMultiArray() returns the number of regions found (= the highest region label, because 
    labels start at 1). 

    <b> Declaration:</b>

    \code
    namespace vigra {
        template <unsigned int N, class T, class S1,
                                  class Label, class S2>
        Label
        watershedsMultiArray(MultiArrayView<N, T, S1> const & data,
                             MultiArrayView<N, Label, S2> labels,  // may also hold input seeds
                             NeighborhoodType neighborhood = DirectNeighborhood,
                             WatershedOptions const & options = WatershedOptions());
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_watersheds.hxx\><br>
    Namespace: vigra

    Example: watersheds of the gradient magnitude (the example works likewise for higher dimensions).

    \code
    MultiArray<2, unsigned char> src(Shape2(w, h));
    ... // read input data
    
    // compute gradient magnitude at scale 1.0 as a boundary indicator
    MultiArray<2, float> gradMag(src.shape());
    gaussianGradientMagnitude(srcImageRange(src), destImage(gradMag), 1.0);

    // example 1
    {
        // the pixel type of the destination image must be large enough to hold
        // numbers up to 'max_region_label' to prevent overflow
        MultiArray<2, unsigned int> labeling(src.shape());
        
        // call region-growing algorithm for 4-neighborhood, leave a 1-pixel boundary between 
        // regions, and autogenerate seeds from all gradient minima where the magnitude is 
        // less than 2.0.
        unsigned int max_region_label = 
              watershedsMultiArray(gradMag, labeling, DirectNeighborhood,
                                   WatershedOptions().keepContours()
                                      .seedOptions(SeedOptions().minima().threshold(2.0)));
    }
    
    // example 2
    {
        MultiArray<2, unsigned int> labeling(src.shape());
        
        // compute seeds beforehand (use connected components of all pixels 
        // where the gradient is below 4.0)
        unsigned int max_region_label = generateWatershedSeeds(gradMag, labeling,
                                                       SeedOptions().levelSets(4.0));
        
        // quantize the gradient image to 256 gray levels
        float m, M;
        gradMag.minmax(&m, &M);
        
        using namespace multi_math;
        MultiArray<2, unsigned char> gradMag256(255.0 / (M - m) * (gradMag - m));
        
        // call region-growing algorithm with 8-neighborhood,
        // since the data are 8-bit, a faster priority queue will be used
        watershedsMultiArray(gradMag256, labeling, IndirectNeighborhood);
    }
    
    // example 3
    {
        MultiArray<2, unsigned int> labeling(src.shape());
        
        .. // put seeds in 'labeling', e.g. from an interactive labeling program,
           // make sure that label 1 corresponds to the background
        
        // bias the watershed algorithm so that the background is preferred
        // by reducing the cost for label 1 to 90%
        watershedsMultiArray(gradMag, labeling, 
                             WatershedOptions().biasLabel(1, 0.9));
    }
    
    // example 4
    {
        MultiArray<2, unsigned int> labeling(src.shape());
        
        // use the fast union-find algorithm with 4-neighborhood
        watershedsMultiArray(gradMag, labeling, WatershedOptions().unionFind());
    }
    \endcode
*/
doxygen_overloaded_function(template <...> Label watershedsMultiArray)

template <unsigned int N, class T, class S1,
                          class Label, class S2>
inline Label
watershedsMultiArray(MultiArrayView<N, T, S1> const & data,
                     MultiArrayView<N, Label, S2> labels,  // may also hold input seeds
                     NeighborhoodType neighborhood = DirectNeighborhood,
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
