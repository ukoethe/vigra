#include <iostream>

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/colorconversions.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/hierarchical_clustering.hxx>
#include <vigra/metrics.hxx>

using namespace vigra;

int main (int argc, char ** argv)
{
    // parameters
    float sigmaGradMag = 3.0f;      // sigma Gaussian gradient
    float beta = 0.5f;              // node vs edge weight
    float wardness = 0.8f;          // ???
    float gamma = 10000000.0f;      // ???
    int nodeNumStop = 50;          // desired num. nodes in result

    if(argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile outfile2" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        std::cout << "(only color images)" << std::endl;

        return 1;
    }
    try
    {
        // read image given as first argument
        ImageImportInfo info(argv[1]);

        // instantiate arrays for image data and for smoothed image of appropriate size
        MultiArray<2, TinyVector<float, 3> > imageArray(info.shape()),
                                             imageArrayBig(info.shape()*2-Shape2(1));
        MultiArray<2, float>                 gradMag(info.shape());

        // read image data
        importImage(info, imageArray);

        // upsample to twice the resolution for visualization of the results
        resizeMultiArraySplineInterpolation(imageArray, imageArrayBig);

        // convert to Lab color space to have better color similarity estimates
        transformMultiArray(imageArray, imageArray, RGB2LabFunctor<float>());

        // compute gradient magnitude as an edge indicator
        gaussianGradientMagnitude(imageArray, gradMag, sigmaGradMag);

        // create grid-graph of appropriate size
        typedef GridGraph<2, undirected_tag > Graph;
        Graph graph(info.shape());

        // create watershed superpixels with the fast union-find algorithm;
        // we use a NodeMap (a subclass of MultiArray) to store the labels so
        // that they can be passed to makeRegionAdjacencyGraph() ??? directly
        Graph::NodeMap<unsigned int> labelArray(graph);
        unsigned int max_label =
            watershedsMultiArray(gradMag, labelArray, DirectNeighborhood,
                                 WatershedOptions().unionFind());

        // visualize the watershed superpixels
        regionImageToCrackEdgeImage(labelArray, imageArrayBig, RGBValue<float>( 255, 0, 0 ), edgeOverlayOnly);
        exportImage(imageArrayBig, argv[2]);

        // determine size and average color of each superpixel
        using namespace acc;
        AccumulatorChainArray<CoupledArrays<2, TinyVector<float, 3>, unsigned int>,
                              Select<DataArg<1>, LabelArg<2>, // where to look for data and region labels
                                     Count, Mean> >           // what statistics to compute
        stats;
        extractFeatures(imageArray, labelArray, stats);

        // create region adjacency graph (RAG) for the superpixels
        // 'affiliatedEdges' stores the grid graph edges belonging to each RAG edge
        typedef AdjacencyListGraph RAG;
        RAG rag;
        RAG::EdgeMap<std::vector<Graph::Edge>> affiliatedEdges(rag);
        makeRegionAdjacencyGraph(graph, labelArray, rag, affiliatedEdges);

        // copy superpixel features into NodeMaps to be passed to agglomerativeClustering???
        RAG::NodeMap<TinyVector<float, 3>> means(rag);
        RAG::NodeMap<unsigned int>         sizes(rag);
        for(unsigned int k=0; k<=max_label; ++k)
        {
            means[k] = get<Mean>(stats, k);
            sizes[k] = get<Count>(stats, k);
        }

        // compute edge length and strength and put into EdgeMaps
        RAG::EdgeMap<float> edgeWeights(rag),
                            edgeLengths(rag);
        for(RAG::EdgeIt rag_edge(rag); rag_edge != lemon::INVALID; ++rag_edge)
        {
            edgeLengths[*rag_edge] = affiliatedEdges[*rag_edge].size();
            // determine RAG edge strength by avaraging the gradients
            // along the RAG edge, as represented by the edges affiliated grid edges
            for(unsigned int k = 0; k < affiliatedEdges[*rag_edge].size(); ++k)
            {
                auto const & grid_edge = affiliatedEdges[*rag_edge][k];
                edgeWeights[*rag_edge] += gradMag[graph.u(grid_edge)] + gradMag[graph.v(grid_edge)];
            }
            edgeWeights[*rag_edge] /= (2.0 * edgeLengths[*rag_edge]);
        }

        // create a merge graph that performs hierarchical clustering
        typedef MergeGraphAdaptor<RAG> MergeGraph;
        MergeGraph mergeGraph(rag);

        RAG::EdgeMap<float> edgeIndicators(rag);
        RAG::EdgeMap<float> edgeMinWeightMap(rag);
        RAG::NodeMap<unsigned int>  nodeLabelMap(rag);

        // create an operator that stores all property maps needed for
        // hierarchical clustering and updates them after every merge step
        typedef cluster_operators::EdgeWeightNodeFeatures<
            MergeGraph,
            RAG::EdgeMap<float>,
            RAG::EdgeMap<float>,
            RAG::NodeMap<TinyVector<float, 3>>,
            RAG::NodeMap<unsigned int>,
            RAG::EdgeMap<float>,
            RAG::NodeMap<unsigned int>>
        RAGFeatures;
        typedef HierarchicalClustering<RAGFeatures> Clustering;

        // are these parameters all needed? What's their meaning?
        RAGFeatures ragFeatures(mergeGraph,
                                edgeIndicators, edgeLengths,
                                means, sizes,
                                edgeMinWeightMap, nodeLabelMap,
                                beta, metrics::ManhattanMetric,
                                wardness, gamma);

        Clustering::Parameter param;
        param.nodeNumStopCond_=nodeNumStop;
        param.buildMergeTreeEncoding_=false;
        param.verbose_=true;
        Clustering clustering(ragFeatures, param);

        clustering.cluster();

        transformMultiArray(labelArray, labelArray,
            [&mergeGraph](unsigned int oldlabel)
            {
                return mergeGraph.reprNodeId(oldlabel);
            });

        regionImageToCrackEdgeImage(labelArray, imageArrayBig, RGBValue<float>( 0, 255, 0), edgeOverlayOnly);
        exportImage(imageArrayBig, argv[3]);

    }
    catch (std::exception & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}