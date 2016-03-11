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
    int nodeNumStop = 30;          // desired num. nodes in result

    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        std::cout << "(only color images)" << std::endl;

        return 1;
    }
    try
    {
        // read metadata of image file given in argv[1]
        ImageImportInfo info(argv[1]);

        // instantiate image arrays of appropriate size
        MultiArray<2, TinyVector<float, 3> > imageArray(info.shape()),
                                             imageArrayBig(info.shape()*2-Shape2(1));
        MultiArray<2, float>                 gradMag(info.shape());

        // read image data
        importImage(info, imageArray);

        // double the image resolution for better visualization of the results
        resizeMultiArraySplineInterpolation(imageArray, imageArrayBig);

        // convert to Lab color space for better color similarity estimates
        transformMultiArray(imageArray, imageArray, RGB2LabFunctor<float>());

        // compute gradient magnitude as an indicator of edge strength
        gaussianGradientMagnitude(imageArray, gradMag, sigmaGradMag);

        // create grid-graph of appropriate size
        typedef GridGraph<2, undirected_tag > Graph;
        Graph graph(info.shape());

        // create watershed superpixels with the fast union-find algorithm;
        // we use a NodeMap (a subclass of MultiArray) to store the labels so
        // that they can be passed to hierarchicalClustering() directly
        Graph::NodeMap<unsigned int> labelArray(graph);
        unsigned int max_label =
            watershedsMultiArray(gradMag, labelArray, DirectNeighborhood,
                                 WatershedOptions().unionFind());

        // visualize the watersheds as a red overlay over the enlarged image
        regionImageToCrackEdgeImage(labelArray, imageArrayBig,
                                    RGBValue<float>( 255, 0, 0 ), EdgeOverlayOnly);

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

        // copy superpixel features into NodeMaps to be passed to hierarchicalClustering()
        RAG::NodeMap<TinyVector<float, 3>> means(rag);
        RAG::NodeMap<unsigned int>         sizes(rag);
        for(unsigned int k=0; k<=max_label; ++k)
        {
            means[k] = get<Mean>(stats, k);
            sizes[k] = get<Count>(stats, k);
        }

        // create edge maps for weights and lengths of the RAG edges
        RAG::EdgeMap<float> edgeWeights(rag),
                            edgeLengths(rag);
        for(RAG::EdgeIt rag_edge(rag); rag_edge != lemon::INVALID; ++rag_edge)
        {
            // the RAG edge length equala the number of grid edges
            // belonging to current RAG edge
            edgeLengths[*rag_edge] = affiliatedEdges[*rag_edge].size();

            // determine RAG edge weight by avaraging the gradients
            // along the current RAG edge, i.e. over the end pixels of the
            // corresponding grid graph edges
            for(unsigned int k = 0; k < affiliatedEdges[*rag_edge].size(); ++k)
            {
                auto const & grid_edge = affiliatedEdges[*rag_edge][k];
                edgeWeights[*rag_edge] += gradMag[graph.u(grid_edge)] + gradMag[graph.v(grid_edge)];
            }
            edgeWeights[*rag_edge] /= (2.0 * edgeLengths[*rag_edge]);
        }

        // create a node map for the new (clustered) region labels and perform
        // clustering to remove unimportant watershed edges
        RAG::NodeMap<unsigned int>  nodeLabelMap(rag);
        hierarchicalClustering(rag, edgeWeights, edgeLengths, means, sizes, nodeLabelMap,
                               ClusteringOptions().minRegionCount(nodeNumStop)
                                                  .nodeFeatureImportance(beta)
                                                  .sizeImportance(wardness)
                                                  .nodeFeatureMetric(metrics::L2Norm)
                               );

        // create label image with the new labels
        transformMultiArray(labelArray, labelArray,
            [&nodeLabelMap](unsigned int oldlabel)
            {
                return nodeLabelMap[oldlabel];
            });

        // visualize the salient edges as a green overlay
        regionImageToCrackEdgeImage(labelArray, imageArrayBig,
                                    RGBValue<float>( 0, 255, 0), EdgeOverlayOnly);

        // write result into image file given by argv[2]
        exportImage(imageArrayBig, argv[2]);

    }
    catch (std::exception & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}