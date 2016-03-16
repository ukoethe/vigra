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
    // parameters of the hierarchical clustering algorithm
    float sigmaGradMag = 3.0f;   // scale of the Gaussian gradient
    float beta = 0.5f;           // importance of node features relative to edge weights
    float wardness = 0.8f;       // importance of cluster size
    int numClusters = 30;        // desired number of resulting regions (clusters)

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

        vigra_precondition(info.numBands() == 3, "an RGB image is required.");

        // instantiate image arrays of appropriate size
        MultiArray<2, TinyVector<float, 3> > imageArrayRGB(info.shape()),
                                             imageArrayLab(info.shape());

        // read image data
        importImage(info, imageArrayRGB);

        // convert to Lab color space for better color similarity estimates
        transformMultiArray(imageArrayRGB, imageArrayLab, RGB2LabFunctor<float>());

        // compute gradient magnitude as an indicator of edge strength
        MultiArray<2, float>   gradMag(imageArrayLab.shape());
        gaussianGradientMagnitude(imageArrayLab, gradMag, sigmaGradMag);

        // create watershed superpixels with the fast union-find algorithm;
        // we use a NodeMap (a subclass of MultiArray) to store the labels so
        // that they can be passed to hierarchicalClustering() directly
        MultiArray<2, unsigned int> labelArray(gradMag.shape());
        unsigned int max_label =
            watershedsMultiArray(gradMag, labelArray, DirectNeighborhood,
                                 WatershedOptions().unionFind());

        // double the image resolution for better visualization of the results
        MultiArray<2, TinyVector<float, 3> > imageArrayBig(info.shape()*2-Shape2(1));
        resizeMultiArraySplineInterpolation(imageArrayRGB, imageArrayBig);

        // visualize the watersheds as a red overlay over the enlarged image
        regionImageToCrackEdgeImage(labelArray, imageArrayBig,
                                    RGBValue<float>( 255, 0, 0 ), EdgeOverlayOnly);

        // create grid-graph of appropriate size
        typedef GridGraph<2, undirected_tag > ImageGraph;
        ImageGraph imageGraph(labelArray.shape());

        // construct empty  region adjacency graph (RAG) for the superpixels
        typedef AdjacencyListGraph RAG;
        RAG rag;

        // create mapping 'affiliatedEdges' from edges in the RAG to
        // corresponding edges in imageGraph and build the RAG
        RAG::EdgeMap<std::vector<ImageGraph::Edge>> affiliatedEdges(rag);
        makeRegionAdjacencyGraph(imageGraph, labelArray, rag, affiliatedEdges);

        // create edge maps for weights and lengths of the RAG edges (zero initialized)
        RAG::EdgeMap<float> edgeWeights(rag),
                            edgeLengths(rag);

        // iterate over all RAG edges (this loop follows a standard LEMON idiom)
        for(RAG::EdgeIt rag_edge(rag); rag_edge != lemon::INVALID; ++rag_edge)
        {
            // iterate over all grid edges that constitute the present RAG edge
            for(unsigned int k = 0; k < affiliatedEdges[*rag_edge].size(); ++k)
            {
                // look up the current grid edge and its end points
                auto const & grid_edge = affiliatedEdges[*rag_edge][k];
                auto start = imageGraph.u(grid_edge),
                     end   = imageGraph.v(grid_edge);

                // compute gradient by linear interpolation between end points
                double grid_edge_gradient = 0.5 * (gradMag[start] + gradMag[end]);
                // aggregate the total
                edgeWeights[*rag_edge] += grid_edge_gradient;
            }

            // the length of the RAG edge equals the number of constituent grid edges
            edgeLengths[*rag_edge] = affiliatedEdges[*rag_edge].size();
            // define edge weight by the average gradient
            edgeWeights[*rag_edge] /= edgeLengths[*rag_edge];
        }

        // determine size and average color of each superpixel
        using namespace acc;
        AccumulatorChainArray<CoupledArrays<2, TinyVector<float, 3>, unsigned int>,
                              Select<DataArg<1>, LabelArg<2>, // where to look for data and region labels
                                     Count, Mean> >           // what statistics to compute
            features;
        extractFeatures(imageArrayLab, labelArray, features);

        // copy superpixel features into NodeMaps to be passed to hierarchicalClustering()
        RAG::NodeMap<TinyVector<float, 3>> meanColor(rag);
        RAG::NodeMap<unsigned int>         regionSize(rag);
        for(unsigned int k=0; k<=max_label; ++k)
        {
            meanColor[k] = get<Mean>(features, k);
            regionSize[k] = get<Count>(features, k);
        }

        // create a node map for the new (clustered) region labels and perform
        // clustering to remove unimportant watershed edges
        RAG::NodeMap<unsigned int>  nodeLabels(rag);
        hierarchicalClustering(rag,          // input: the superpixel adjacency graph
                               edgeWeights, edgeLengths, meanColor, regionSize, // features
                               nodeLabels,   // output: a cluster labeling of the RAG
                               ClusteringOptions().minRegionCount(numClusters)
                                                  .nodeFeatureImportance(beta)
                                                  .sizeImportance(wardness)
                                                  .nodeFeatureMetric(metrics::L2Norm)
                               );

        // create label image with the new labels
        transformMultiArray(labelArray, labelArray,
            [&nodeLabels](unsigned int oldlabel)
            {
                return nodeLabels[oldlabel];
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