#ifndef VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX
#define VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX
//#define NO_IMPORT_ARRAY

/*boost python before anything else*/
#include <boost/python.hpp>

/*std*/
#include <sstream>
#include <string>

/*vigra*/
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/graphs.hxx>
#include <vigra/graph_maps.hxx>
#include <vigra/python_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/metrics.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/error.hxx>
#include <vigra/multi_watersheds.hxx>
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
    typedef NumpyArray<EdgeMapDim,   Singleband<UInt32> > UInt32EdgeArray;
    typedef NumpyArray<EdgeMapDim,   Singleband<Int32 > > Int32EdgeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<float > > FloatNodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<UInt32> > UInt32NodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<Int32 > > Int32NodeArray;
    typedef NumpyArray<NodeMapDim +1,Multiband <float > > MultiFloatNodeArray;
    typedef NumpyArray<EdgeMapDim +1,Multiband <float > > MultiFloatEdgeArray;


    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>         FloatEdgeArrayMap;
    typedef NumpyScalarEdgeMap<Graph,UInt32EdgeArray>        UInt32EdgeArrayMap;
    typedef NumpyScalarEdgeMap<Graph,Int32EdgeArray>         Int32EdgeArrayMap;
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
        python::def("_edgeWeightedWatershedsSegmentation",registerConverters(&pyEdgeWeightedWatershedsSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );

        python::def("_nodeWeightedWatershedsSegmentation",registerConverters(&pyNodeWeightedWatershedsSegmentation),
            (
                python::arg("graph"),
                python::arg("nodeWeights"),
                python::arg("seeds"),
                python::arg("method")=std::string("regionGrowing"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a node weighted graph"
        );
        python::def("_nodeWeightedWatershedsSeeds",registerConverters(&pyNodeWeightedWatershedsSeeds),
            (
                python::arg("graph"),
                python::arg("nodeWeights"),
                python::arg("out")=python::object()
            ),
            "Generate seeds for node weighted watersheds"
        );

        python::def("_carvingSegmentation",registerConverters(&pyCarvingSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("backgroundLabel"),
                python::arg("backgroundBias"),
                python::arg("noBiasBelow") = 0.0,
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );


        python::def("_shortestPathSegmentation",registerConverters(&pyShortestPathSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("nodeWeights"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            ),
            "Seeded shorted path segmentation on a edge and node weighted graph"
        );


        python::def("_felzenszwalbSegmentation",registerConverters(&pyFelzenszwalbSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("nodeSizes"),
                python::arg("k")=300.0f,
                python::arg("nodeNumStop")=-1,
                python::arg("out")=python::object()
            ),
            "Felzenwalb graph based segmentation"
        );
    }

    void exportMiscAlgorithms()const{

        python::def("_nodeFeatureDistToEdgeWeight",registerConverters(&pyNodeFeatureDistToEdgeWeight),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("metric"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights with the given metric"
        );
        python::def("_nodeFeatureSumToEdgeWeight",registerConverters(&pyNodeFeatureSumToEdgeWeight),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights"
        );

        python::def("_opengmMulticutDataStructure",registerConverters(&pyMulticutDataStructure),
            (
                python::arg("graph"),
                python::arg("edgeWeights")
            )
        );

        

        python::def("nodeGtToEdgeGt",registerConverters(&pyNodeGtToEdgeGt),
            (
                python::arg("graph"),
                python::arg("nodeGt"),
                python::arg("ignoreLabel"),
                python::arg("out")=python::object()
            )
        );

        python::def("_opengmArgToLabeling",registerConverters(&pyMulticutArgToLabeling),
            (
                python::arg("graph"),
                python::arg("arg"),
                python::arg("out")=python::object()
            )
        );
        python::def("_wardCorrection",registerConverters(&pyWardCorrection),
            (
                python::arg("graph"),
                python::arg("edgeIndicator"),
                python::arg("nodeSize"),
                python::arg("out")=python::object()
            ),
            "apply wards method to an edgeIndicator"
        );

        python::def("find3Cycles", registerConverters(&pyFind3Cycles));
        python::def("find3CyclesEdges", registerConverters(&pyFind3CyclesEdges));
        python::def("cyclesEdges", registerConverters(&pyCyclesEdges),
            (
                python::arg("graph"),
                python::arg("graph"),
                python::arg("out") = python::object()
            )
        );
    }

    void exportSmoothingAlgorithms()const{

        python::def("_recursiveGraphSmoothing",registerConverters(&pyRecursiveGraphSmoothing),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("edgeIndicator"),
                python::arg("gamma"),
                python::arg("edgeThreshold"),
                python::arg("scale"),
                python::arg("iterations")=1,
                python::arg("outBuffer")=python::object(),
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
        exportSmoothingAlgorithms();
    }


    static NumpyAnyArray pyFind3Cycles(
        const GRAPH & graph
    ){
        NumpyArray<1, vigra::TinyVector<Int32, 3> > cycles;

        MultiArray<1, vigra::TinyVector<Int32, 3> > cyclesArray;
        find3Cycles(graph, cyclesArray);
        cycles.reshapeIfEmpty(cyclesArray.shape());
        cycles = cyclesArray;
        return cycles;
    }


    static NumpyAnyArray pyFind3CyclesEdges(
        const GRAPH & graph
    ){
        NumpyArray<1, vigra::TinyVector<Int32, 3> > cyclesEdges;

        MultiArray<1, vigra::TinyVector<Int32, 3> > cyclesNodes;

        find3Cycles(graph, cyclesNodes);
        cyclesEdges.reshapeIfEmpty(cyclesNodes.shape());
        
        Node nodes[3];
        Edge edges[3];

        for(size_t i=0; i<cyclesNodes.size(); ++i){
            for(size_t j=0; j<3; ++j){
                nodes[j] = graph.nodeFromId(cyclesNodes(i)[j]);
            }
            edges[0] = graph.findEdge(nodes[0],nodes[1]);
            edges[1] = graph.findEdge(nodes[0],nodes[2]);
            edges[2] = graph.findEdge(nodes[1],nodes[2]);
            for(size_t j=0; j<3; ++j){
                cyclesEdges(i)[j] = graph.id(edges[j]);
            }
        }

        return cyclesEdges;
    }

    static NumpyAnyArray pyCyclesEdges(
        const GRAPH & graph,
        NumpyArray<1, vigra::TinyVector<Int32, 3> > cycles,
        NumpyArray<1, vigra::TinyVector<Int32, 3> > edgesOut       
    ){

        Node nodes[3];
        Edge edges[3];
        Int32 ids[3];

        edgesOut.reshapeIfEmpty(cycles.shape());
        for(size_t i=0; i<cycles.size(); ++i){
            for(size_t j=0; j<3; ++j){
                nodes[j] = graph.nodeFromId(cycles(i)[j]);
            }
            edges[0] = graph.findEdge(nodes[0],nodes[1]);
            edges[1] = graph.findEdge(nodes[0],nodes[2]);
            edges[2] = graph.findEdge(nodes[1],nodes[2]);
            for(size_t j=0; j<3; ++j){
                edgesOut(i)[j] = graph.id(edges[j]);
            }
        }   
        return edgesOut;
    }

    

    static NumpyAnyArray pyWardCorrection(
        const Graph &           g,
        const FloatEdgeArray    edgeWeightsArray,
        const FloatNodeArray    nodeSizeArray,
        const float             wardness,
        FloatEdgeArray    outArray
    ){
        outArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g));

        // numpy arrays => lemon maps
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        FloatNodeArrayMap  nodeSizeArrayMap(g,nodeSizeArray);
        FloatEdgeArrayMap  outArrayMap(g,outArray);

        for(EdgeIt iter(g);iter!=lemon::INVALID;++iter){
            const float uSize=nodeSizeArrayMap[g.u(*iter)];
            const float vSize=nodeSizeArrayMap[g.v(*iter)];
            const float w = edgeWeightsArrayMap[*iter];
            const float ward  = 1.0f/(1.0f/std::log(uSize) + 1.0f/std::log(vSize)  );
            const float wardF = wardness*ward + (1.0-wardness);
            outArrayMap[*iter]=w*wardF;
        }
        return outArray;

    }




    static python::tuple pyMulticutDataStructure(
        const Graph &           g,
        const FloatEdgeArray    edgeWeightsArray
    ){
        UInt32NodeArray toDenseArray( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g));

        // numpy arrays => lemon maps
        UInt32NodeArrayMap toDenseArrayMap(g,toDenseArray);
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);

        NumpyArray<2,UInt32> vis      ((    typename NumpyArray<2,UInt64>::difference_type(g.edgeNum(),2)));
        NumpyArray<1,float > weights  ((    typename NumpyArray<1,double>::difference_type(g.edgeNum()  )));
        
        size_t denseIndex = 0 ;
        for(NodeIt iter(g);iter!=lemon::INVALID;++iter){
            toDenseArrayMap[*iter]=denseIndex;
            ++denseIndex;
        }
        denseIndex=0;
        for(EdgeIt iter(g);iter!=lemon::INVALID;++iter){
            const size_t dU=toDenseArrayMap[g.u(*iter)];
            const size_t dV=toDenseArrayMap[g.v(*iter)];
            vis(denseIndex,0)=std::min(dU,dV);
            vis(denseIndex,1)=std::max(dU,dV);
            weights(denseIndex)=edgeWeightsArrayMap[*iter];
            ++denseIndex;
        }
        return python::make_tuple(vis,weights);

    }


    static NumpyAnyArray pyNodeGtToEdgeGt(
        const Graph &           g,
        const UInt32NodeArray & nodeGt,
        const Int64 ignoreLabel,
        UInt32EdgeArray edgeGt
    ){
        edgeGt.reshapeIfEmpty(IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g));

        // numpy arrays => lemon maps
        UInt32NodeArrayMap nodeGtMap(g,nodeGt);
        UInt32EdgeArrayMap  edgeGtMap(g,edgeGt);
        nodeGtToEdgeGt(g, nodeGtMap, ignoreLabel, edgeGtMap);
        return edgeGt;
    }



    static NumpyAnyArray pyMulticutArgToLabeling(
        const Graph &              g,
        const NumpyArray<1,UInt32> arg,
        UInt32NodeArray            labelsArray
    ){
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g));
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);
        size_t denseIndex = 0 ;
        for(NodeIt iter(g);iter!=lemon::INVALID;++iter){
            labelsArrayMap[*iter]=arg(denseIndex);
            ++denseIndex;
        }
        return labelsArray;
    }
    /*
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
    */
    static NumpyAnyArray pyNodeFeatureDistToEdgeWeight(
        const GRAPH & g,
        const MultiFloatNodeArray & nodeFeaturesArray,
        const std::string & functor,
        FloatEdgeArray edgeWeightsArray
    ){
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        if(functor==std::string("euclidean") || functor==std::string("norm") || functor==std::string("l2")){
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
                "- euclidean/norm/l2\n"
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

    static NumpyAnyArray pyEdgeWeightedWatershedsSegmentation(
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
        edgeWeightedWatershedsSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyNodeWeightedWatershedsSegmentation(
        const Graph &       g,
        FloatNodeArray      nodeWeightsArray,
        UInt32NodeArray     seedsArray,
        const std::string & method,
        UInt32NodeArray     labelsArray
    ){

        // resize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );


        WatershedOptions watershedsOption;
        if(method==std::string("regionGrowing"))
            watershedsOption.regionGrowing();
        else
            watershedsOption.unionFind();

        // numpy arrays => lemon maps
        FloatNodeArrayMap  nodeWeightsArrayMap(g,nodeWeightsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        std::copy(seedsArray.begin(),seedsArray.end(),labelsArray.begin());

        //lemon_graph::graph_detail::generateWatershedSeeds(g, nodeWeightsArrayMap, labelsArrayMap, watershedsOption.seed_options);
        lemon_graph::watershedsGraph(g, nodeWeightsArrayMap, labelsArrayMap, watershedsOption);
        //lemon_graph::graph_detail::seededWatersheds(g, nodeWeightsArrayMap, seedsArrayMap, watershedsOption);
        
        return labelsArray;
    }





    static NumpyAnyArray pyNodeWeightedWatershedsSeeds(
        const Graph &       g,
        FloatNodeArray      nodeWeightsArray,
        UInt32NodeArray     seedsArray
    ){
        const std::string method="regionGrowing";
        // resize output ? 
        seedsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        WatershedOptions watershedsOption;
        if(method==std::string("regionGrowing"))
            watershedsOption.regionGrowing();

        // numpy arrays => lemon maps
        FloatNodeArrayMap  nodeWeightsArrayMap(g,nodeWeightsArray);
        UInt32NodeArrayMap seedsArrayMap(g,seedsArray);

        lemon_graph::graph_detail::generateWatershedSeeds(g, nodeWeightsArrayMap, seedsArrayMap, watershedsOption.seed_options);

        return seedsArray;
    }

    static NumpyAnyArray pyCarvingSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        UInt32NodeArray seedsArray,
        const UInt32    backgroundLabel,
        const float     backgroundBias,
        const float     noBiasBelow,
        UInt32NodeArray labelsArray
    ){
        // resize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        UInt32NodeArrayMap seedsArrayMap(g,seedsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        carvingSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,backgroundLabel,backgroundBias,noBiasBelow,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyShortestPathSegmentation(
        const Graph &       g,
        FloatEdgeArray      edgeWeightsArray,
        FloatNodeArray      nodeWeightsArray,
        UInt32NodeArray     seedsArray,
        UInt32NodeArray     labelsArray
    ){

        // resize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        FloatNodeArrayMap  nodeWeightsArrayMap(g,nodeWeightsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);



        std::copy(seedsArray.begin(),seedsArray.end(),labelsArray.begin());

        shortestPathSegmentation<
            Graph,FloatEdgeArrayMap, FloatNodeArrayMap, UInt32NodeArrayMap, float
        >(g, edgeWeightsArrayMap, nodeWeightsArrayMap, labelsArrayMap);
     
        
        return labelsArray;
    }


    static NumpyAnyArray pyFelzenszwalbSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        FloatNodeArray nodeSizesArray,
        const float k,
        const int nodeNumStop,
        UInt32NodeArray labelsArray
    ){
        // resize output ? 
        labelsArray.reshapeIfEmpty(  IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        FloatNodeArrayMap  nodeSizesArrayMap(g,nodeSizesArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        felzenszwalbSegmentation(g,edgeWeightsArrayMap,nodeSizesArrayMap,k,labelsArrayMap,nodeNumStop);

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
        TaggedShape inShape  = nodeFeaturesArray.taggedShape();
        TaggedShape outShape = TaggedGraphShape<Graph>::taggedNodeMapShape(g);
        if(inShape.hasChannelAxis()){
            outShape.setChannelCount(inShape.channelCount());
        }
        nodeFeaturesBufferArray.reshapeIfEmpty(outShape);
        nodeFeaturesOutArray.reshapeIfEmpty(outShape);

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
    typedef NumpyArray<EdgeMapDim +1,Multiband <float > > MultiFloatEdgeArray;


    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>         FloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<Graph,FloatNodeArray>         FloatNodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,UInt32NodeArray>        UInt32NodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,Int32NodeArray>         Int32NodeArrayMap;
    typedef NumpyMultibandNodeMap<Graph,MultiFloatNodeArray> MultiFloatNodeArrayMap;
    typedef NumpyMultibandEdgeMap<Graph,MultiFloatEdgeArray> MultiFloatEdgeArrayMap;

    typedef ShortestPathDijkstra<Graph,float> ShortestPathDijkstraType;


    typedef typename GraphDescriptorToMultiArrayIndex<Graph>::IntrinsicNodeMapShape NodeCoordinate;
    typedef NumpyArray<1,NodeCoordinate>  NodeCoorinateArray;

    LemonGridGraphAlgorithmAddonVisitor(const std::string & clsName){}


    template <class classT>
    void visit(classT& c) const
    {   

        // - edge weights from interpolated image
        exportMiscAlgorithms(c);

    }

    template <class classT>
    void exportMiscAlgorithms(classT & c)const{
        


        python::def("edgeFeaturesFromInterpolatedImage",registerConverters(&pyEdgeWeightsFromInterpolatedImage),
            (
                python::arg("graph"),
                python::arg("image"),
                python::arg("out")=python::object()
            ),
            "convert an image with with ``shape = graph.shape*2 - 1`` to an edge weight array"
        );

        python::def("edgeFeaturesFromImage",registerConverters(&pyEdgeWeightsFromImage),
            (
                python::arg("graph"),
                python::arg("image"),
                python::arg("out")=python::object()
            ),
            "convert an image with with shape = graph.shape OR shape = graph.shape *2 -1 to an edge weight array"
        );

        python::def("edgeFeaturesFromImage",registerConverters(&pyEdgeWeightsFromImageMb),
            (
                python::arg("graph"),
                python::arg("image"),
                python::arg("out")=python::object()
            ),
            "convert an image with with shape = graph.shape OR shape = graph.shape *2 -1 to an edge weight array"
        );


        c
        .def("affiliatedEdgesSerializationSize",&pyAffiliatedEdgesSerializationSize,
            (
                python::arg("rag"),
                python::arg("affiliatedEdges")
            )
        );


        //'python::def("edgeFeaturesFromInterpolatedImageCorrected",registerConverters(&pyEdgeWeightsFromInterpolatedImageCorrected),
        //'    (
        //'        python::arg("graph"),
        //'        python::arg("image"),
        //'        python::arg("out")=python::object()
        //'    ),
        //'    "convert an image with with shape = graph.shape *2 -1 to an edge weight array"
        //'    ""
        //');


    }


    


    static size_t pyAffiliatedEdgesSerializationSize(
        const GRAPH & gridGraph,
        const AdjacencyListGraph & rag,
        const typename AdjacencyListGraph:: template EdgeMap< std::vector<Edge> > & affiliatedEdges
    ){
        return affiliatedEdgesSerializationSize(gridGraph, rag, affiliatedEdges);
    }


    static NumpyAnyArray pyEdgeWeightsFromImage(
        const GRAPH & g,
        const FloatNodeArray & image,
        FloatEdgeArray edgeWeightsArray
    ){

        bool regularShape=true;
        bool topologicalShape=true;

        for(size_t d=0;d<NodeMapDim;++d){
            if(image.shape(d)!=g.shape()[d]){
                regularShape=false;
            }
            if(image.shape(d)!=2*g.shape()[d]-1){
                topologicalShape=false;
            }
        }
       
        if(regularShape)
            return pyEdgeWeightsFromOrginalSizeImage(g,image,edgeWeightsArray);
        else if(topologicalShape)
            return pyEdgeWeightsFromInterpolatedImage(g,image,edgeWeightsArray);
        else{
            vigra_precondition(false, "shape of edge image does not match graph shape");
            // to avid no return warnings
            return pyEdgeWeightsFromOrginalSizeImage(g,image,edgeWeightsArray);
        }
    }


    static NumpyAnyArray pyEdgeWeightsFromImageMb(
        const GRAPH & g,
        const MultiFloatNodeArray & image,
        MultiFloatEdgeArray edgeWeightsArray
    ){

        bool regularShape=true;
        bool topologicalShape=true;

        for(size_t d=0;d<NodeMapDim;++d){
            if(image.shape(d)!=g.shape()[d]){
                regularShape=false;
            }
            if(image.shape(d)!=2*g.shape()[d]-1){
                topologicalShape=false;
            }
        }
       
        if(regularShape)
            return pyEdgeWeightsFromOrginalSizeImageMb(g,image,edgeWeightsArray);
        else if(topologicalShape)
            return pyEdgeWeightsFromInterpolatedImageMb(g,image,edgeWeightsArray);
        else{
            vigra_precondition(false, "shape of edge image does not match graph shape");
            // to avid no return warnings
            return pyEdgeWeightsFromOrginalSizeImageMb(g,image,edgeWeightsArray);
        }
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
            const CoordType vCoord(g.v(edge));
            const CoordType tCoord = uCoord+vCoord;
            edgeWeightsArrayMap[edge]=interpolatedImage[tCoord];
        }
        return edgeWeightsArray;
    }

    static NumpyAnyArray pyEdgeWeightsFromOrginalSizeImage(
        const GRAPH & g,
        const FloatNodeArray & image,
        FloatEdgeArray edgeWeightsArray
    ){

        for(size_t d=0;d<NodeMapDim;++d){
            //std::cout<<"is "<<image.shape(d)<<"gs "<<2*g.shape()[d]-1<<"\n";
            vigra_precondition(image.shape(d)==g.shape()[d], "interpolated shape must be shape*2 -1");
        }
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        typedef typename FloatNodeArray::difference_type CoordType;
        for(EdgeIt iter(g); iter!=lemon::INVALID; ++ iter){

            const Edge edge(*iter);
            const CoordType uCoord(g.u(edge));
            const CoordType vCoord(g.v(edge));
            edgeWeightsArrayMap[edge]=(image[uCoord]+image[vCoord])/2.0;
        }
        return edgeWeightsArray;
    }



    static NumpyAnyArray pyEdgeWeightsFromInterpolatedImageMb(
        const GRAPH & g,
        const MultiFloatNodeArray & interpolatedImage,
        MultiFloatEdgeArray edgeWeightsArray
    ){

        for(size_t d=0;d<NodeMapDim;++d){
            //std::cout<<"is "<<interpolatedImage.shape(d)<<"gs "<<2*g.shape()[d]-1<<"\n";
            vigra_precondition(interpolatedImage.shape(d)==2*g.shape()[d]-1, "interpolated shape must be shape*2 -1");
        }

        // resize out
        typename MultiArray<EdgeMapDim+1,int>::difference_type outShape;
        for(size_t d=0;d<EdgeMapDim;++d){
            outShape[d]=IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g)[d];
        }
        outShape[EdgeMapDim] = interpolatedImage.shape(NodeMapDim);

        edgeWeightsArray.reshapeIfEmpty(   MultiFloatEdgeArray::ArrayTraits::taggedShape(outShape,"nc") );

        
        // numpy arrays => lemon maps
        MultiFloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        typedef typename FloatNodeArray::difference_type CoordType;
        for(EdgeIt iter(g); iter!=lemon::INVALID; ++ iter){

            const Edge edge(*iter);
            const CoordType uCoord(g.u(edge));
            const CoordType vCoord(g.v(edge));
            const CoordType tCoord = uCoord+vCoord;
            edgeWeightsArrayMap[edge]=interpolatedImage[tCoord];
        }
        return edgeWeightsArray;
    }


    static NumpyAnyArray pyEdgeWeightsFromOrginalSizeImageMb(
        const GRAPH & g,
        const MultiFloatNodeArray & image,
        MultiFloatEdgeArray edgeWeightsArray
    ){

        for(size_t d=0;d<NodeMapDim;++d){
            //std::cout<<"is "<<image.shape(d)<<"gs "<<2*g.shape()[d]-1<<"\n";
            vigra_precondition(image.shape(d)==g.shape()[d], "interpolated shape must be shape*2 -1");
        }

        // resize out
        typename MultiArray<EdgeMapDim+1,int>::difference_type outShape;
        for(size_t d=0;d<EdgeMapDim;++d){
            outShape[d]=IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g)[d];
        }
        outShape[EdgeMapDim] = image.shape(NodeMapDim);


        //edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(outShape),"ec" );


        edgeWeightsArray.reshapeIfEmpty(   MultiFloatEdgeArray::ArrayTraits::taggedShape(outShape,"nc") );


        // numpy arrays => lemon maps
        MultiFloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        typedef typename FloatNodeArray::difference_type CoordType;
        for(EdgeIt iter(g); iter!=lemon::INVALID; ++ iter){

            const Edge edge(*iter);
            const CoordType uCoord(g.u(edge));
            const CoordType vCoord(g.v(edge));
            MultiArray<1, float>  val = image[uCoord];
            val+=image[vCoord];
            val/=2.0;
            edgeWeightsArrayMap[edge]=val;
        }
        return edgeWeightsArray;
    }


    /*
    static NumpyAnyArray pyEdgeWeightsFromInterpolatedImageCorrected(
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
            const CoordType vCoord(g.v(edge));
            const CoordType tCoord = uCoord+vCoord;
            int diffCounter = 0;
            for(int i=0; i<NodeMapDim; ++i) {
                if (uCoord[i] != vCoord[i]) {
                    diffCounter++;
                }
            }
            edgeWeightsArrayMap[edge]=sqrt(diffCounter)*interpolatedImage[tCoord];
        }
        return edgeWeightsArray;
    }
    */
};



} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX
