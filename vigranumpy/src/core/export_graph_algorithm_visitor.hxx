#ifndef VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX
#define VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX
//#define NO_IMPORT_ARRAY

/*std*/
#include <sstream>
#include <string>

/*vigra*/
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <boost/python.hpp>
#include <vigra/graphs.hxx>
#include <vigra/graph_helper/dense_map.hxx>
#include <vigra/graph_helper/on_the_fly_edge_map.hxx>
#include <vigra/python_graph_generalization.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/metrics.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/error.hxx>
namespace python = boost::python;

namespace vigra{


template<class GRAPH>
class LemonGraphAlgorithmVisitor 
:   public boost::python::def_visitor<LemonGraphAlgorithmVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonGraphAlgorithmVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;

    typedef typename Graph::NodeIt              NodeIt;
    typedef typename Graph::EdgeIt              EdgeIt;
    typedef typename Graph::ArcIt               ArcIt;


    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;


    // predefined array (for map usage)
    const static unsigned int EdgeMapDim = IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension;
    const static unsigned int NodeMapDim = IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension;

    typedef NumpyArray<EdgeMapDim,   Singleband<float > > FloatEdgeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<float > > FloatNodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<UInt32> > UInt32NodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<Int32 > > Int32NodeArray;
    typedef NumpyArray<NodeMapDim +1,Multiband <float > > MultiFloatNodeArray;

    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>         FloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<Graph,FloatNodeArray>         FloatNodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,UInt32NodeArray>        UInt32NodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,Int32NodeArray>         Int32NodeArrayMap;
    typedef NumpyMultibandNodeMap<Graph,MultiFloatNodeArray> MultiFloatNodeArrayMap;


    typedef ShortestPathDijkstra<Graph,float> ShortestPathDijkstraType;


    typedef typename GraphDescriptorToMultiArrayIndex<Graph>::IntrinsicNodeMapShape NodeCoordinate;
    typedef NumpyArray<1,NodeCoordinate>  NodeCoorinateArray;

    LemonGraphAlgorithmVisitor(const std::string clsName)
    :clsName_(clsName){

    }


    void exportSegmentationAlgorithms()const{
        python::def("watershedsSegmentation",registerConverters(&pyWatershedSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );

        python::def("carvingSegmentation",registerConverters(&pyCarvingSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("backgroundLabel"),
                python::arg("backgroundBias"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );

        python::def("felzenszwalbSegmentation",registerConverters(&pyFelzenszwalbSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("nodeSizes")=python::object(),
                python::arg("k")=300.0f,
                python::arg("out")=python::object()
            ),
            "Felzenwalb graph based segmentation"
        );
    }

    void exportMiscAlgorithms()const{
        python::def("_nodeIdsLabels",registerConverters(&pyNodeIdsLabels),
            (
                python::arg("graph"),
                python::arg("nodeIds"),
                python::arg("labels"),
                python::arg("out")=python::object()
            )
        );

        python::def("_nodeIdsFeatures",registerConverters(&pyNodeIdsFeatures),
            (
                python::arg("graph"),
                python::arg("nodeIds"),
                python::arg("features"),
                python::arg("out")=python::object()
            )
        );

        python::def("nodeFeatureDistToEdgeWeight",registerConverters(&pyNodeFeatureDistToEdgeWeight),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("metric"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights with the given metric"
        );
        python::def("nodeFeatureSumToEdgeWeight",registerConverters(&pyNodeFeatureSumToEdgeWeight),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights with the given metric"
        );
    }

    void exportSmoothingAlgorithms()const{

        python::def("recursiveGraphSmoothing",registerConverters(&pyRecursiveGraphSmoothing),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("edgeIndicator"),
                python::arg("gamma"),
                python::arg("edgeThreshold"),
                python::arg("scale"),
                python::arg("iterations")=1,
                python::arg("buffer")=python::object(),
                python::arg("out")=python::object()
            ),
            "recursive edge weighted guided graph smoothing"
        );

        python::def("dynamicRecursiveGraphSmoothing",registerConverters(&pyDynamicRecursiveGraphSmoothing),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("metric"),
                python::arg("gamma"),
                python::arg("edgeThreshold"),
                python::arg("scale"),
                python::arg("iterations")=1,
                python::arg("buffer")=python::object(),
                python::arg("out")=python::object()
            ),
            "recursive edge weighted guided graph smoothing"
        );
    }

    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {   
        // - watersheds-segmentation
        // - carving-segmentation
        // - felzenwalb-segmentation
        // - labeling
        exportSegmentationAlgorithms();

        // - node Labels (usefull to make grid graph labels from rag labels)
        // - node feature distance to edge weights
        exportMiscAlgorithms();

        // - recursiveGraphSmoothing
        // - dynamicRecursiveGraphSmoothing
        exportSmoothingAlgorithms();
    }

    static NumpyAnyArray pyNodeIdsLabels(
        const GRAPH & g,
        NumpyArray<1,Singleband<UInt32> >  nodeIds,
        UInt32NodeArray                    nodeLabelArray,
        NumpyArray<1,Singleband<UInt32> >  out
    ){
        // reshape out
        out.reshapeIfEmpty(nodeIds.shape());

        // numpy arrays => lemon maps
        UInt32NodeArrayMap nodeLabelArrayMap(g,nodeLabelArray);

        for(size_t i=0;i<nodeIds.shape(0);++i)
            out(i)=nodeLabelArrayMap[g.nodeFromId(nodeIds(i))];
        return out;
    }
    
    static NumpyAnyArray pyNodeIdsFeatures(
        const GRAPH & g,
        NumpyArray<1,Singleband<UInt32> >  nodeIds,
        MultiFloatNodeArray                nodeFeaturesArray,
        NumpyArray<2,Multiband<float >  >  out
    ){
        //  reshape out ?
        typename NumpyArray<2,Multiband<float> >::difference_type outShape(nodeIds.shape(0),nodeFeaturesArray.shape(NodeMapDim));
        out.reshapeIfEmpty(  NumpyArray<2,Multiband<float >  >::ArrayTraits::taggedShape(outShape,"xc"));

        // numpy arrays => lemon maps
        MultiFloatNodeArrayMap nodeFeaturesArrayMap(g,nodeFeaturesArray);

        typedef typename  NumpyArray<1,int>::difference_type Coord1;
        for(size_t i=0;i<nodeIds.shape(0);++i)
            out[Coord1(i)]=nodeFeaturesArrayMap[g.nodeFromId(nodeIds(i))];
        return out;
    }

    static NumpyAnyArray pyNodeFeatureDistToEdgeWeight(
        const GRAPH & g,
        const MultiFloatNodeArray & nodeFeaturesArray,
        const std::string & functor,
        FloatEdgeArray edgeWeightsArray
    ){
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        if(functor==std::string("eucledian") || functor==std::string("norm") || functor==std::string("l2")){
            typedef  metrics::Norm<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        if(functor==std::string("squaredNorm")){
            typedef  metrics::SquaredNorm<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else if (functor==std::string("manhattan") || functor==std::string("l1")){
            typedef  metrics::Manhattan<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else if (functor==std::string("chiSquared")){
            typedef  metrics::ChiSquared<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else{
            throw std::runtime_error(
                "distance not supported\n"
                "supported distance types:\n"
                "- eucledian/norm/l2\n"
                "- squaredNorm\n"
                "- manhattan/l1\n"
                "- chiSquared\n"
            );
        }
    }



    static NumpyAnyArray pyNodeFeatureSumToEdgeWeight(
        const GRAPH & g,
        const FloatNodeArray & nodeFeaturesArray,
        FloatEdgeArray edgeWeightsArray
    ){
        // reshape out?
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        // numpy arrays => lemon maps
        FloatNodeArrayMap  nodeFeatureArrayMap(g,nodeFeaturesArray);
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        
        for(EdgeIt e(g);e!=lemon::INVALID;++e){
            const Edge edge(*e);
            const Node u=g.u(edge);
            const Node v=g.v(edge);
            edgeWeightsArrayMap[edge]=nodeFeatureArrayMap[u]+nodeFeatureArrayMap[v];
        }
        return edgeWeightsArray;
    }

    template<class FUNCTOR>
    static NumpyAnyArray pyNodeFeatureDistToEdgeWeightT(
        const GRAPH & g,
        const MultiFloatNodeArray & nodeFeaturesArray,
        FUNCTOR & functor,
        FloatEdgeArray edgeWeightsArray
    ){
        // reshape out?
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        // numpy arrays => lemon maps
        MultiFloatNodeArrayMap nodeFeatureArrayMap(g,nodeFeaturesArray);
        FloatEdgeArrayMap      edgeWeightsArrayMap(g,edgeWeightsArray);
        
        for(EdgeIt e(g);e!=lemon::INVALID;++e){
            const Edge edge(*e);
            const Node u=g.u(edge);
            const Node v=g.v(edge);
            edgeWeightsArrayMap[edge]=functor(nodeFeatureArrayMap[u],nodeFeatureArrayMap[v]);
        }
        return edgeWeightsArray;
    }

    static NumpyAnyArray pyWatershedSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        UInt32NodeArray seedsArray,
        UInt32NodeArray labelsArray
    ){
        // resize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        UInt32NodeArrayMap seedsArrayMap(g,seedsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        watershedsSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyCarvingSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        UInt32NodeArray seedsArray,
        const UInt32    backgroundLabel,
        const float     backgroundBias,
        UInt32NodeArray labelsArray
    ){
        // resize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        UInt32NodeArrayMap seedsArrayMap(g,seedsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        carvingSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,backgroundLabel,backgroundBias,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyFelzenszwalbSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        FloatNodeArray nodeSizesArray,
        const float k,
        UInt32NodeArray labelsArray
    ){
        // resize output ? 
        labelsArray.reshapeIfEmpty(  IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // is size array empty?
        // if empty fill with ones
        if(nodeSizesArray.shape(0)==0){
            nodeSizesArray.reshapeIfEmpty(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g));
            std::fill(nodeSizesArray.begin(),nodeSizesArray.end(),1.0);
        }

        // numpy arrays => lemon maps
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        FloatNodeArrayMap  nodeSizesArrayMap(g,nodeSizesArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        felzenszwalbSegmentation(g,edgeWeightsArrayMap,nodeSizesArrayMap,k,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyRecursiveGraphSmoothing(
        const GRAPH & g,
        MultiFloatNodeArray nodeFeaturesArray,
        FloatEdgeArray      edgeIndicatorArray,
        const float         lambda,
        const float         edgeThreshold,
        const float         scale,
        const size_t        iterations,
        MultiFloatNodeArray nodeFeaturesBufferArray,
        MultiFloatNodeArray nodeFeaturesOutArray
    ){
        // resize output ? 
        nodeFeaturesOutArray.reshapeIfEmpty( nodeFeaturesArray.taggedShape());

        // resize buffer ? 
        nodeFeaturesBufferArray.reshapeIfEmpty( nodeFeaturesArray.taggedShape());

        // numpy arrays => lemon maps
        MultiFloatNodeArrayMap nodeFeaturesArrayMap(g,nodeFeaturesArray);
        FloatEdgeArrayMap edgeIndicatorArrayMap(g,edgeIndicatorArray);
        MultiFloatNodeArrayMap nodeFeaturesBufferArrayMap(g,nodeFeaturesBufferArray);
        MultiFloatNodeArrayMap nodeFeaturesOutArrayMap(g,nodeFeaturesOutArray);

        // call algorithm itself
        recursiveGraphSmoothing(g,nodeFeaturesArrayMap,edgeIndicatorArrayMap,lambda,edgeThreshold,scale,iterations,nodeFeaturesBufferArrayMap,nodeFeaturesOutArrayMap);

        // retun smoothed features
        return nodeFeaturesOutArray;
    }


    static NumpyAnyArray pyDynamicRecursiveGraphSmoothing(
        const GRAPH & g,
        MultiFloatNodeArray nodeFeaturesArray,
        const std::string & functor,
        const float         lambda,
        const float         edgeThreshold,
        const float         scale,
        const size_t        iterations,
        MultiFloatNodeArray nodeFeaturesBufferArray,
        MultiFloatNodeArray nodeFeaturesOutArray
    ){
        // resize output ? 
        nodeFeaturesOutArray.reshapeIfEmpty( nodeFeaturesArray.taggedShape());

        // resize buffer ? 
        nodeFeaturesBufferArray.reshapeIfEmpty( nodeFeaturesArray.taggedShape());

        // numpy arrays => lemon maps
        MultiFloatNodeArrayMap nodeFeaturesArrayMap(g,nodeFeaturesArray);
        MultiFloatNodeArrayMap nodeFeaturesBufferArrayMap(g,nodeFeaturesBufferArray);
        MultiFloatNodeArrayMap nodeFeaturesOutArrayMap(g,nodeFeaturesOutArray);

        // call algorithm itself
        if(functor==std::string("eucledian") || functor==std::string("norm") || functor==std::string("l2")){
            typedef  metrics::Norm<float> DistFunctor;
            DistFunctor f;
            dynamicRecursiveGraphSmoothing(g,nodeFeaturesArrayMap,lambda,edgeThreshold,scale,f,iterations,nodeFeaturesBufferArrayMap,nodeFeaturesOutArrayMap);
        }
        else if(functor==std::string("squaredNorm")){
            typedef  metrics::SquaredNorm<float> DistFunctor;
            DistFunctor f;
            dynamicRecursiveGraphSmoothing(g,nodeFeaturesArrayMap,lambda,edgeThreshold,scale,f,iterations,nodeFeaturesBufferArrayMap,nodeFeaturesOutArrayMap);
        }
        else if (functor==std::string("manhattan") || functor==std::string("l1")){
            typedef  metrics::Manhattan<float> DistFunctor;
            DistFunctor f;
            dynamicRecursiveGraphSmoothing(g,nodeFeaturesArrayMap,lambda,edgeThreshold,scale,f,iterations,nodeFeaturesBufferArrayMap,nodeFeaturesOutArrayMap);
        }
        else if (functor==std::string("chiSquared")){
            typedef  metrics::ChiSquared<float> DistFunctor;
            DistFunctor f;
            dynamicRecursiveGraphSmoothing(g,nodeFeaturesArrayMap,lambda,edgeThreshold,scale,f,iterations,nodeFeaturesBufferArrayMap,nodeFeaturesOutArrayMap);
        }
        else{
            std::cout<<"selected distance "<<functor<<"\n";
            throw std::runtime_error(
                "distance not supported\n"
                "supported distance types:\n"
                "- eucledian/norm/l2\n"
                "- squaredNorm\n"
                "- manhattan/l1\n"
                "- chiSquared\n"
            );
        }


        // retun smoothed features
        return nodeFeaturesOutArray;
    }

};


template<class GRAPH>
class LemonGridGraphAlgorithmAddonVisitor 
:   public boost::python::def_visitor<LemonGridGraphAlgorithmAddonVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonGraphAlgorithmVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;

    typedef typename Graph::NodeIt              NodeIt;
    typedef typename Graph::EdgeIt              EdgeIt;
    typedef typename Graph::ArcIt               ArcIt;


    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;


    // predefined array (for map usage)
    const static unsigned int EdgeMapDim = IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension;
    const static unsigned int NodeMapDim = IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension;

    typedef NumpyArray<EdgeMapDim,   Singleband<float > > FloatEdgeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<float > > FloatNodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<UInt32> > UInt32NodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<Int32 > > Int32NodeArray;
    typedef NumpyArray<NodeMapDim +1,Multiband <float > > MultiFloatNodeArray;

    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>         FloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<Graph,FloatNodeArray>         FloatNodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,UInt32NodeArray>        UInt32NodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,Int32NodeArray>         Int32NodeArrayMap;
    typedef NumpyMultibandNodeMap<Graph,MultiFloatNodeArray> MultiFloatNodeArrayMap;


    typedef ShortestPathDijkstra<Graph,float> ShortestPathDijkstraType;


    typedef typename GraphDescriptorToMultiArrayIndex<Graph>::IntrinsicNodeMapShape NodeCoordinate;
    typedef NumpyArray<1,NodeCoordinate>  NodeCoorinateArray;

    LemonGridGraphAlgorithmAddonVisitor(const std::string & clsName){}


    template <class classT>
    void visit(classT& c) const
    {   

        // - edge weights from interpolated image
        exportMiscAlgorithms();

    }


    void exportMiscAlgorithms()const{
        


        python::def("edgeFeaturesFromInterpolatedImage",registerConverters(&pyEdgeWeightsFromInterpolatedImage),
            (
                python::arg("graph"),
                python::arg("image"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights with the given metric"
        );
    }



    static NumpyAnyArray pyEdgeWeightsFromInterpolatedImage(
        const GRAPH & g,
        const FloatNodeArray & interpolatedImage,
        FloatEdgeArray edgeWeightsArray
    ){

        for(size_t d=0;d<NodeMapDim;++d){
            //std::cout<<"is "<<interpolatedImage.shape(d)<<"gs "<<2*g.shape()[d]-1<<"\n";
            vigra_precondition(interpolatedImage.shape(d)==2*g.shape()[d]-1, "interpolated shape must be shape*2 -1");
        }


        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        typedef typename FloatNodeArray::difference_type CoordType;
        for(EdgeIt iter(g); iter!=lemon::INVALID; ++ iter){

            const Edge edge(*iter);
            const CoordType uCoord(g.u(edge));
            const CoordType vCoord(g.u(edge));
            const CoordType tCoord = uCoord+vCoord;
            edgeWeightsArrayMap[edge]=interpolatedImage[tCoord];
        }
        return edgeWeightsArray;
    }
};



} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX