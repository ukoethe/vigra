#include <valgrind/callgrind.h>

#include <vigra/hdf5impex.hxx>
#include <vigra/slic.hxx>

#include <vigra/convolution.hxx>
#include <vigra/labelvolume.hxx>

#include <vigra/hiracical_clustering.hxx>
#include <vigra/rag/rag.hxx>

int main() {   
    using namespace vigra;

    // import volume fromn hdf5
    const std::string path = "data.h5";
    const std::string dset = "data";


    std::cout<<"read in data\n";
    HDF5ImportInfo info(path.c_str(), dset.c_str());
    vigra_precondition(info.numDimensions() == 3, "Dataset must be 3-dimensional.");
    MultiArrayShape<3>::type shape(info.shape().begin());
    MultiArray<3, float> rawDataBig(shape);
    readHDF5(info, rawDataBig);

    // get a subarray set is smaller by one element at all sides
    typedef MultiArray<3, float>::difference_type Shape;
    MultiArrayView <3, float> rawData = rawDataBig.subarray(Shape(0,0,0), Shape(100, 100, 100));

    std::cout<<"oversegment with slic\n";
    // oversegment with slic
    MultiArray<3, unsigned int> labelsRaw(rawData.shape());
    int seedDistance = 3;
    double intensityScaling = 20.0;
    // compute seeds automatically, perform 40 iterations, and scale intensity differences
    // down to 1/20 before comparing with spatial distances
    slicSuperpixels(rawData, labelsRaw, intensityScaling, seedDistance, SlicOptions().iterations(3));

    std::cout<<"get label volumen \n";
    MultiArray<3, unsigned int> labels(labelsRaw.shape());
    // find 6-connected regions
    int max_region_label = labelVolumeSix(labelsRaw, labels);

    std::cout<<"max_region_label "<<max_region_label<<"\n";

    std::cout<<"gradmag\n";
    // use a 3-dimensional float array
    MultiArray<3, float>  gradmag(rawData.shape());
    // calculate gradient magnitude at scale = 3.0
    gaussianGradientMagnitude(rawData, gradmag, 2.0);





    std::cout<<"get rag\n";
    typedef Rag<3,unsigned int> RagType;
    typedef typename GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
    typedef DenseEdgeReferenceMap<RagType,float>                        RagEdgeFloatMap;
    typedef MultiArray<RagType::Dimension , float  >                    SingleBandFloatImage;



    RagType rag(labels);
    std::cout<<" alloc maps\n";
    RagEdgeCoordinatesMap coordMapmap(rag);
    RagEdgeFloatMap       edgeIndicatorMap(rag);
    RagEdgeFloatMap       edgeSizeMap(rag);

    std::cout<<"extract smaps\n";
    extractEdgeCoordinates(rag,coordMapmap);
    edgeSizeMapFromCoordMap(rag,coordMapmap,edgeSizeMap);
    extractEdgeFeaturesFromImage(rag,coordMapmap,gradmag,edgeIndicatorMap);


    CALLGRIND_START_INSTRUMENTATION;


    std::cout<<"Do clustering\n";
    typedef HiracicalClustering<RagType,RagEdgeFloatMap,RagEdgeFloatMap> HiracicalClusteringType;
    HiracicalClusteringType hcluster(rag,edgeIndicatorMap,edgeSizeMap);

    std::cout<<"Do clustering\n";
    hcluster.cluster();


    CALLGRIND_STOP_INSTRUMENTATION;
   
}
