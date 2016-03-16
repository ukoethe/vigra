#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/convolution.hxx>
#include <vigra/multi_blockwise.hxx>
#include <iostream>
#include <vigra/adjacency_list_graph.hxx>

#include <vigra/watersheds.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/rgbvalue.hxx>


using namespace vigra;

int main (int argc, char ** argv) 
{
    // parameters
    float sigmaGradMag = 5.0;      // sigma Gaussian gradient
    float watershedScale = 1.0;    // SLIC color - spatial weight
    float beta = 0.5;              // node vs edge weight
    int superpixelDiameter = 10;   // super-pixel size
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
        MultiArray<2, RGBValue<float> > imageArray(info.shape()),
                                        imageArrayBig(info.shape()*2-Shape2(1)),
                                        exportArray(info.shape());
        MultiArray<2, float> gradMag(info.shape());
        MultiArray<2, float> gradMagBig(info.shape()*2-Shape2(1));
        MultiArray<2, unsigned int> labelArray(info.shape());
    
        // copy image data into array
        importImage(info, imageArray);
        resizeMultiArraySplineInterpolation(imageArray, imageArrayBig);
        
        // compute gradient magnitude as a boundary indicator
        gaussianGradientMagnitude(imageArray, gradMag, sigmaGradMag);
        resizeMultiArraySplineInterpolation(gradMag, gradMagBig);
        
        // watershed ... use the fast union-find algorithm with 4-neighborhood
        watershedsMultiArray(gradMag, labelArray, DirectNeighborhood, WatershedOptions().unionFind());
        
        
        // regionImageToCrackEdgeImage(imageArray, imageArrayBig, RGBValue<int>( 256, 0, 0 ));
        regionImageToCrackEdgeImage(labelArray, imageArrayBig, RGBValue<float>( 255, 0, 0 ), edgeOverlayOnly);
        
        /*
        // get 2D grid graph and edgeMap for grid graph
        // from gradMag of interpolated image
        GridGraph<2, undirected_tag > gidGraph (info.shape());
        // somehow fill graph...
        EdgeMap affiliatedEdges (gridGraph);
        
        // get region adjacency graph
        AdjacencyListGraph rag;
        makeRegionAdjacencyGraph(gridGraph, labels, rag, affiliatedEdges)
        */
        /*
        # accumulate edge weights from gradient magnitude
        edgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)

        # accumulate node features from grid graph node map
        # which is just a plain image (with channels)
        nodeFeatures = rag.accumulateNodeFeatures(imgLab)

        # do agglomerativeClustering
        labels = graphs.agglomerativeClustering(graph=rag, edgeWeights=edgeWeights,
                                                beta=beta, nodeFeatures=nodeFeatures,
                                                nodeNumStop=nodeNumStop,wardness=0.8)
             */                                   
        
        
        // write image data to the file given as second argument
        exportImage(labelArray, ImageExportInfo(argv[2]));
        exportImage(imageArray, ImageExportInfo(argv[3]));
        
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}