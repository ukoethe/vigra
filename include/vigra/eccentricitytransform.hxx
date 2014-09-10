/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Philip Schill and Ullrich Koethe               */
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


#ifndef VIGRA_ECCENTRICITYTRANSFORM_HXX
#define VIGRA_ECCENTRICITYTRANSFORM_HXX

/*std*/
#include <algorithm>
#include <set>

/*vigra*/
#include "accumulator.hxx"
#include "multi_labeling.hxx"
#include "multi_distance.hxx"
#include "multi_resize.hxx"
#include "graph_algorithms.hxx"


namespace vigra
{
#if 0
    /// \brief Creates edge weights using a slightly modified inverse distance transform.
    ///
    /// \param img : labeled image
    /// \param g : graph of the image
    /// \param edgeWeights[out] : the computed edge weights
    /// \param background : background value of img
    /// \param INFVAL : optional edge weights of background
    template<unsigned int N, class T, class S, class GRAPH>
    void createEdgeWeightsFromInverseDistanceTransform(
            const MultiArrayView<N, T> & img,
            const GRAPH & g,
            typename GRAPH:: template EdgeMap<S> & edgeWeights,
            const T & background,
            const S & INFVAL = 1000000.0
    ){
        typedef typename MultiArrayShape<N>::type Shape;

        // Get the distance transform and invert the values.
        MultiArray<N, S> distImage(img.shape());
        separableMultiDistance(img, distImage, true);
        S m = *std::max_element(distImage.begin(), distImage.end());
        distImage *= -1;
        distImage += m+10;

        // Set the outside (where the background was) to a very high value.
        for (int i=0; i<distImage.size(); ++i)
        {
            if (img[i] == background)
            {
                distImage[i] = INFVAL;
            }
        }

        // Create the edge weights.
        Shape graphShape = img.shape()*2-1;
        MultiArray<N, S> interpolatedImage(graphShape);
        BSpline< 0, S > spline;
        resizeMultiArraySplineInterpolation(distImage, interpolatedImage, spline);
        edgeWeightsFromInterpolatedImage(g, interpolatedImage, edgeWeights, true);
    }

    /// \brief Compute shortest path between each two vertices, find the longest of those paths and return all vertices on that path.
    ///
    /// \param pathFinder : dijkstra path finder
    /// \param edgeWeights : edge weights for dijkstra shortest path algorithm
    /// \param edges : edges of the given object
    /// \param shortestPath [out] : shortest path between each two vertices
    /// \param maxIterations : optional number of iterations for the approximation
    template<class GRAPH, class T, class EDGEMAP>
    void findLongestShortestPath(
            ShortestPathDijkstra<GRAPH, T> & pathFinder,
            const EDGEMAP & edgeWeights,
            const MultiArray<2, MultiArrayIndex> & edges,
            MultiArray<2, MultiArrayIndex> & shortestPath,
            unsigned int maxIterations = 4
    ){
        typedef typename GRAPH::Node Node;
        typedef typename GRAPH:: template NodeMap<Node> PredecessorsMap;

        const int nEdges = edges.shape()[0];
        const int N = edges.shape()[1];

        // For only one edge, no computation must be done.
        // This makes sure, that our path later on consists of at least 2 points.
        if (nEdges == 1)
        {
            shortestPath = MultiArray<2, MultiArrayIndex>(edges.shape());
            std::copy(edges.begin(), edges.end(), shortestPath.begin());
            return;
        }

        PredecessorsMap predecessors;
        Node srcOfMaxPath;
        Node tarOfMaxPath;

        bool visited[nEdges];
        std::fill_n(visited, nEdges, false);
        int srcIndex = 0;
        T maxPathLength = -1;
        for (int i=0; i<maxIterations; ++i)
        {
            visited[srcIndex] = true;
            Node srcNode;
            for (int d=0; d<N; ++d)
            {
                srcNode[d] = edges(srcIndex, d);
            }
            pathFinder.run(edgeWeights, srcNode);
            int maxIndex;
            bool maxChanged = false;
            for (int j=0; j<nEdges; ++j)
            {
                Node targetNode;
                for (int d=0; d<N; ++d)
                {
                    targetNode[d] = edges(j, d);
                }
                double pathLength = pathFinder.distance(targetNode);
                if (pathLength > maxPathLength)
                {
                    maxChanged = true;
                    maxPathLength = pathLength;
                    maxIndex = j;
                }
            }
            if (maxChanged)
            {
                predecessors = PredecessorsMap(pathFinder.predecessors());
                for (int d=0; d<N; ++d)
                {
                    srcOfMaxPath[d] = edges(srcIndex, d);
                    tarOfMaxPath[d] = edges(maxIndex, d);
                }
                srcIndex = maxIndex;
            }
            if ( (!maxChanged) || visited[srcIndex])
            {
                break;
            }
        }

        // Write the predecessors map into the MultiArray.
        int len = pathLength(srcOfMaxPath, tarOfMaxPath, predecessors);
        shortestPath = MultiArray<2, MultiArrayIndex>(Shape2(len, N));
        int k=0;
        Node currentNode = tarOfMaxPath;
        while(currentNode != srcOfMaxPath)
        {
            for (int d=0; d<N; ++d)
            {
                shortestPath(k, d) = currentNode[d];
            }
            currentNode=predecessors[currentNode];
            k++;
        }
        for (int d=0; d<N; ++d)
        {
            shortestPath(k, d) = srcOfMaxPath[d];
        }
    }

    /// \brief Finds the the eccentricity center for a single region.
    ///
    /// \param src : labeled image
    /// \param label : label of region
    template <unsigned int N, class T, class S>
    TinyVector<MultiArrayIndex, N> findEccentricityCenterSingleRegion(
            const MultiArrayView<N, T> & src,
            const T & label
    ){
        typedef typename MultiArrayShape<N>::type Shape;
        typedef GridGraph<N, undirected_tag> Graph;
        //typedef typename Graph::Node Node;
        typedef typename Graph:: template NodeMap<UInt32> UInt32NodeMap;
        typedef typename Graph:: template EdgeMap<S> EdgeMap;
        typedef typename AdjacencyListGraph:: template EdgeMap< std::vector<typename Graph::Edge> > RagAffiliatedEdges;
        typedef typename AdjacencyListGraph::Node RagNode;

        const UInt32 BACKGROUNDVALUE = 1;

        // For small images, no computation has to be done.
        if (src.size() <= 1) {
            return TinyVector<T, N>(0);
        }

        // Enlarge the image by 1 pixel on each side and set it 0 on the object.
        Shape bigShape = src.shape()+2;
        MultiArray<N, UInt32> bigImg(bigShape);
        bigImg.init(BACKGROUNDVALUE);
        Shape upperLeftCut(1);
        Shape lowerRightCut = bigShape-1;
        MultiArrayView<N, UInt32> sliceBigImg = bigImg.subarray(upperLeftCut, lowerRightCut);
        for (int i=0; i<sliceBigImg.size(); ++i)
        {
            if (src[i] == label)
            {
                sliceBigImg[i] = 0;
            }
        }

        // Relabel the image and find the label for the inside.
        MultiArray<N, UInt32> bigLblImg(bigShape);
        labelMultiArray(bigImg, bigLblImg, IndirectNeighborhood);
        UInt32 insideLabel=0;
        for (int i=0; i<bigImg.size(); ++i)
        {
            if (bigImg[i] == 0)
            {
                insideLabel = bigLblImg[i];
                break;
            }
        }
        vigra_precondition(insideLabel > 0, "The given label does not occure in the source image.");

        // Create edge weights and path finder.
        Graph gGraph(bigImg.shape(), IndirectNeighborhood);
        EdgeMap edgeWeights(gGraph);
        createEdgeWeightsFromInverseDistanceTransform< N, T, S, Graph >(bigImg, gGraph, edgeWeights, BACKGROUNDVALUE);
        ShortestPathDijkstra< Graph , S > pathFinder(gGraph);

        // Find the border in the image.
        AdjacencyListGraph rag;
        RagAffiliatedEdges affEdges;
        UInt32NodeMap bigLblNodes(gGraph);
        for (int i=0; i<bigLblNodes.size(); ++i)
        {
            bigLblNodes[i] = bigLblImg[i];
        }
        makeRegionAdjacencyGraph(gGraph, bigLblNodes, rag, affEdges);
        RagNode node = rag.nodeFromId(insideLabel);
        MultiArray<2, MultiArrayIndex> edges = ragFindEdges(rag, gGraph, affEdges, bigLblImg, node);

        // Get the diameter (path) of the object.
        MultiArray<2, MultiArrayIndex> diameter;
        findLongestShortestPath(pathFinder, edgeWeights, edges, diameter);

        // Find euclidean length of the diameter (unweighted edges).
        S diameterLength = 0;
        for (int k=1; k<diameter.shape()[0]; ++k)
        {
            int diffCount = 0;
            for (int d=0; d<N; ++d)
            {
                if (diameter(k, d) != diameter(k-1, d))
                {
                    diffCount++;
                }
            }
            diameterLength += sqrt(diffCount);
        }

        // Find index of the midpoint of the diameter.
        S midLength = diameterLength/2;
        S currentLength = 0;
        int midPointIndex = 0;
        for (int k=1; k<diameter.shape()[0]; ++k)
        {
            S addLength = 0;
            for (int d=0; d<N; ++d)
            {
                addLength += (diameter(k, d) - diameter(k-1, d)) * (diameter(k, d) - diameter(k-1, d));
            }
            currentLength += sqrt(addLength);
            if (currentLength > midLength)
            {
                if ( currentLength-midLength < midLength - currentLength + addLength )
                {
                    midPointIndex = k;
                }
                else
                {
                    midPointIndex = k-1;
                }
                break;
            }
        }

        // The image was enlarged by 1 one all sides, so subtract that from the midpoint coordinates.
        TinyVector<MultiArrayIndex, N> midpoint;
        for (int d=0; d<N; ++d)
        {
            midpoint[d] = diameter(midPointIndex, d)-1;
        }

        // Sometimes, the midpoint lies outside the object. In this case, take the next point from the diameter.
        if (bigImg[midpoint+1] == BACKGROUNDVALUE)
        {
            for (int d=0; d<N; ++d)
            {
                midpoint[d] = diameter(midPointIndex+1, d)-1;
            }
        }

        return midpoint;
    }

    /// \brief Computes the eccentricity transform (approximation). In dest, only pixels of the region with the given label are touched.
    ///
    /// \param src : labeled image
    /// \param dest[out] : eccentricity transform of the region "label" in src
    /// \param label : label of the region
    template <unsigned int N, class T, class S>
    void eccentricityTransform(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest,
            const T & label
    ) {
        TinyVector<MultiArrayIndex, N> center;
        eccentricityTransform<N, T, S>(src, dest, label, center);
    }

    /// \brief Computes the eccentricity transform (approximation). In dest, only pixels of the region with the given label are touched.
    ///
    /// \param src : labeled image
    /// \param dest[out] : eccentricity transform of the region "label" in src
    /// \param label : label of the region
    /// \param center[out] : eccentricity center
    template <unsigned int N, class T, class S>
    void eccentricityTransform(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest,
            const T & label,
            TinyVector<MultiArrayIndex, N> & center
    ){
        typedef typename MultiArrayShape<N>::type Shape;
        typedef GridGraph<N, undirected_tag> Graph;
        //typedef typename Graph::Node Node;
        typedef typename Graph:: template EdgeMap<S> EdgeMap;

        const S INFVAL = 1000000;

        // Check if the image size is correct.
        for (int d=0; d<N; ++d)
        {
            vigra_precondition(src.shape()[d] == dest.shape()[d], "Image sizes must match.");
        }

        // Create edge weights for the dijkstra algorithm.
        Graph gGraph(src.shape(), IndirectNeighborhood);
        EdgeMap edgeWeights(gGraph);
        Shape graphShape = src.shape()*2-1;
        MultiArray<N, S> interpolatedImage(graphShape);
        interpolatedImage.init(1);
        for (int i=0; i<src.size(); ++i)
        {
            if (src[i] != label)
            {
                interpolatedImage[2*i+1] = INFVAL;
            }
        }
        edgeWeightsFromInterpolatedImage(gGraph, interpolatedImage, edgeWeights, true);

        // Run dijkstra from the eccentricity center.
        TinyVector<MultiArrayIndex, N> eccCenter = findEccentricityCenterSingleRegion< N, T, S >(src, label);
        for (int i=0; i<N; ++i)
        {
            center[i] = eccCenter[i];
        }
        ShortestPathDijkstra< Graph , S > pathFinder(gGraph);
        pathFinder.run(edgeWeights, center);

        // Write the distances into the dest image.
        for (int i=0; i<src.size(); ++i)
        {
            if (src[i] == label)
            {
                dest[i] = pathFinder.distances()[i];
            }
        }
    }

    /// \brief Computes the eccentricity transform (approximation) on each region of a labeled image.
    ///
    /// \param src : labeled array
    /// \param dest[out] : eccentricity transform of src
    template <unsigned int N, class T, class S>
    void eccentricityTransformOnLabels(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest
    ){
        ArrayVector<TinyVector<float, N> > centers;
        eccentricityTransformOnLabels<N, T, S>(src, dest, centers);
    }

    /// \brief Computes the eccentricity transform (approximation) on each region of a labeled image.
    ///
    /// \param src : labeled array
    /// \param dest[out] : eccentricity transform of src
    /// \param centers[out] : list of eccentricity centers (<tt>centers[k] = TinyVector<float, N>()</tt>
    ///                        must be supported)    
    template <unsigned int N, class T, class S, class Array>
    void 
    eccentricityTransformOnLabels(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest,
            Array & centers)
    {
        using namespace acc;
        typedef typename MultiArrayShape<N>::type Shape;
        typedef GridGraph<N> Graph;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::EdgeIt EdgeIt;
        typedef float WeightType;
        
        vigra_precondition(src.shape() == dest.shape(), 
            "eccentricityTransformOnLabels(): Shape mismatch between src and dest.");
            
        MultiArray<N, WeightType> distances(src.shape());
        boundaryMultiDistance(src, distances, true);

        AccumulatorChainArray<CoupledArrays<N, WeightType, T>,
                Select< DataArg<1>, LabelArg<2>,
                        Maximum, Count, BoundingBox, RegionAnchor> > a;
        extractFeatures(distances, src, a);
        
        Graph g(src.shape(), IndirectNeighborhood);
        typename Graph::template EdgeMap<WeightType> weights(g);
        
        for (EdgeIt edge(g); edge != lemon::INVALID; ++edge)
        {
            const Node u(g.u(*edge)), v(g.v(*edge));
            const T label = src[u];
            if(label != src[v])
                weights[*edge] = NumericTraits<WeightType>::max();
            else
                weights[*edge] = norm(u - v) * 
                                  (get<Maximum>(a, label) + 10.0 - 0.5*(distances[u] + distances[v]));
        }
        
        T maxLabel = a.maxRegionLabel();
        centers.resize(maxLabel+1);

        for (T i=0; i <= maxLabel; ++i)
        {
            if(get<Count>(a, i) == 0)
                continue;
            eccentricityTransformGraph(g, weights, get<RegionAnchor>(a, i), get<BoundingBox>(a, i), 
                                       dest, centers[i]);
        }
    }


    /// \brief Finds the eccentricity center for each region in the source image (first column will be empty, since labels start with 1).
    ///
    /// \param src : labeled image
    /// \param centers [out] : eccentricity centers
    /// \param maxLabel : largest value in src
    template <unsigned int N, class T>
    void findEccentricityCenters(
            const MultiArrayView<N, T> & src,
            MultiArrayView<2, MultiArrayIndex> & centers,
            T maxLabel = 0
    ){
        using namespace acc;

        typedef typename MultiArrayShape<N>::type Shape;

        if (maxLabel == 0)
        {
            maxLabel = *std::max_element(src.begin(), src.end());
        }
        vigra_precondition(centers.shape() == Shape2(maxLabel, N), "Input array must have the shape (maxLabel, N).");

        // Extract bounding boxes.
        MultiArray<N, T> data(src);
        AccumulatorChainArray<CoupledArrays<N, T, T>,
                Select< DataArg<1>, LabelArg<2>,
                    Coord< Minimum >,
                    Coord< Maximum > > > a;
        extractFeatures(data, src, a);

        // Fill the array.
        for (T i=0; i<maxLabel; ++i)
        {
            Shape topLeftBB = get< Coord< Minimum > >(a, i+1);
            Shape bottomRightBB = get< Coord< Maximum > >(a, i+1)+1;
            MultiArrayView<N, T> region = src.subarray(topLeftBB, bottomRightBB);

            // Call the function for a single region.
            TinyVector<MultiArrayIndex, N> center = findEccentricityCenterSingleRegion< N, T, float >(region, i+1);
            for (int d=0; d<N; ++d)
            {
                centers(i,d) = center[d] + topLeftBB[d];
            }
        }
    }
#endif

    template <class Graph, class WeightType, 
              class EdgeMap, class Shape>
    TinyVector<float, Shape::static_size> 
    eccentricityCentersOneRegionImpl(ShortestPathDijkstra<Graph, WeightType> & pathFinder,
                            const EdgeMap & weights, WeightType maxWeight,
                            Shape anchor, Shape const & start, Shape const & stop)
    {
        int maxIterations = 4;
        for(int k=0; k < maxIterations; ++k)
        {
            pathFinder.run(weights, anchor, lemon::INVALID, maxWeight, start, stop);
            anchor = pathFinder.target();
            // FIXME: implement early stopping when source and target don't change anymore
        }
        
        Polygon<TinyVector<float, Shape::static_size> > path;
        path.push_back_unsafe(anchor);
        while(pathFinder.predecessors()[path.back()] != path.back())
            path.push_back_unsafe(pathFinder.predecessors()[path.back()]);
        if(path.size() == 1)
            return TinyVector<float, Shape::static_size>(anchor);
        ArrayVector<double> arcLength;
        path.arcLengthList(arcLength);
        double halfLength = 0.5*arcLength.back();
        unsigned int k=0;
        for(; k<arcLength.size(); ++k)
            if(arcLength[k] >= halfLength)
                break;
        double fraction = (halfLength - arcLength[k-1]) / (arcLength[k] - arcLength[k-1]);
        return path.interpolate(k-1, fraction);
    }

    template <unsigned int N, class T, class S, class Graph,
              class ACCUMULATOR, class DIJKSTRA, class Array>
    void 
    eccentricityCentersImpl(const MultiArrayView<N, T, S> & src,
                            Graph const & g,
                            ACCUMULATOR const & r,
                            DIJKSTRA & pathFinder,
                            Array & centers)
    {
        using namespace acc;
        typedef typename MultiArrayShape<N>::type Shape;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::EdgeIt EdgeIt;
        typedef float WeightType;
        
        typename Graph::template EdgeMap<WeightType> weights(g);
        WeightType maxWeight = 0.0;
        {
            AccumulatorChainArray<CoupledArrays<N, WeightType, T>,
                                  Select< DataArg<1>, LabelArg<2>, Maximum> > a;
        
            MultiArray<N, WeightType> distances(src.shape());
            boundaryMultiDistance(src, distances, true);
            extractFeatures(distances, src, a);
            for (EdgeIt edge(g); edge != lemon::INVALID; ++edge)
            {
                const Node u(g.u(*edge)), v(g.v(*edge));
                const T label = src[u];
                if(label != src[v])
                {
                    weights[*edge] = NumericTraits<WeightType>::max();
                }
                else
                {
                    WeightType weight = norm(u - v) * 
                                      (get<Maximum>(a, label) + 10.0 - 0.5*(distances[u] + distances[v]));
                    weights[*edge] = weight;
                    maxWeight = std::max(weight, maxWeight);
                }
            }
        }
        maxWeight *= src.size();
        
        T maxLabel = r.maxRegionLabel();
        centers.resize(maxLabel+1);

        for (T i=0; i <= maxLabel; ++i)
        {
            if(get<Count>(r, i) == 0)
                continue;
            centers[i] = eccentricityCentersOneRegionImpl(pathFinder, weights, maxWeight, 
                                                 get<RegionAnchor>(r, i), 
                                                 get<Coord<Minimum> >(r, i), 
                                                 get<Coord<Maximum> >(r, i) + Shape(1));
        }
    }

    /// \brief Find the (approximate) eccentricitycenter in each region of a labeled image.
    ///
    /// \param src : labeled array
    /// \param centers[out] : list of eccentricity centers (<tt>centers[k] = TinyVector<float, N>()</tt>
    ///                        must be supported)    
    template <unsigned int N, class T, class S, class Array>
    void 
    eccentricityCenters(const MultiArrayView<N, T, S> & src,
                        Array & centers)
    {
        using namespace acc;
        typedef GridGraph<N> Graph;
        typedef float WeightType;
        
        Graph g(src.shape(), IndirectNeighborhood);
        ShortestPathDijkstra<Graph, WeightType> pathFinder(g);

        AccumulatorChainArray<CoupledArrays<N, T>,
                              Select< DataArg<1>, LabelArg<1>,
                                      Count, BoundingBox, RegionAnchor> > a;
        extractFeatures(src, a);
        
        eccentricityCentersImpl(src, g, a, pathFinder, centers);
    }

    /// \brief Computes the eccentricity transform (approximation) on each region of a labeled image.
    ///
    /// \param src : labeled array
    /// \param dest[out] : eccentricity transform of src
    template <unsigned int N, class T, class S>
    void eccentricityTransformOnLabels(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest
    ){
        ArrayVector<TinyVector<float, N> > centers;
        eccentricityTransformOnLabels(src, dest, centers);
    }

    /// \brief Computes the eccentricity transform (approximation) on each region of a labeled image.
    ///
    /// \param src : labeled array
    /// \param dest[out] : eccentricity transform of src
    /// \param centers[out] : list of eccentricity centers (<tt>centers[k] = TinyVector<float, N>()</tt>
    ///                        must be supported)    
    template <unsigned int N, class T, class S, class Array>
    void 
    eccentricityTransformOnLabels(
            const MultiArrayView<N, T> & src,
            MultiArrayView<N, S> & dest,
            Array & centers)
    {
        using namespace acc;
        typedef typename MultiArrayShape<N>::type Shape;
        typedef GridGraph<N> Graph;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::EdgeIt EdgeIt;
        typedef float WeightType;
        
        vigra_precondition(src.shape() == dest.shape(), 
            "eccentricityTransformOnLabels(): Shape mismatch between src and dest.");
            
        using namespace acc;
        typedef GridGraph<N> Graph;
        typedef float WeightType;
        
        Graph g(src.shape(), IndirectNeighborhood);
        ShortestPathDijkstra<Graph, WeightType> pathFinder(g);

        AccumulatorChainArray<CoupledArrays<N, T>,
                              Select< DataArg<1>, LabelArg<1>,
                                      Count, BoundingBox, RegionAnchor> > a;
        extractFeatures(src, a);
        
        eccentricityCentersImpl(src, g, a, pathFinder, centers);

        typename Graph::template EdgeMap<WeightType> weights(g);
        for (EdgeIt edge(g); edge != lemon::INVALID; ++edge)
        {
            const Node u(g.u(*edge)), v(g.v(*edge));
            const T label = src[u];
            if(label != src[v])
            {
                weights[*edge] = NumericTraits<WeightType>::max();
            }
            else
            {
                weights[*edge] = norm(u - v);
            }
        }
        WeightType maxWeight = M_SQRT2*src.size();
        T maxLabel = a.maxRegionLabel();

        for (T i=0; i <= maxLabel; ++i)
        {
            if(get<Count>(a, i) == 0)
                continue;
            Shape center = roundi(centers[i]),
                  start = get<Coord<Minimum> >(a, i),
                  stop  = get<Coord<Maximum> >(a, i) + Shape(1);
            vigra_invariant(src[center] == i,
                "eccentricityTransformOnLabels(): internal error -- center rounding failed");
            pathFinder.run(weights, center, lemon::INVALID, maxWeight, start, stop);
            
            typename CoupledIteratorType<N, T, WeightType, S>::type j = 
                createCoupledIterator(src, pathFinder.distances(), dest);
            j.restrictToSubarray(start, stop);
            for( ; j.isValid(); ++j)
                if(get<1>(*j) == i)
                    get<3>(*j) = get<2>(*j);
        }
    }

} // namespace vigra


#endif // VIGRA_ECCENTRICITYTRANSFORM_HXX
