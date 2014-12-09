#ifndef VIGRA_EXPORT_GRAPH_RAG_VISITOR_HXX
#define VIGRA_EXPORT_GRAPH_RAG_VISITOR_HXX
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
#include <vigra/graph_rag_project_back.hxx>
namespace python = boost::python;

namespace vigra{



template<class GRAPH>
class LemonGraphRagVisitor 
:   public boost::python::def_visitor<LemonGraphRagVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;
    typedef AdjacencyListGraph RagGraph;

    typedef LemonGraphRagVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;
    typedef typename Graph::NodeIt           NodeIt;
    typedef typename Graph::EdgeIt           EdgeIt;
    typedef typename Graph::ArcIt            ArcIt;

    typedef typename RagGraph::Edge             RagEdge;
    typedef typename RagGraph::Node             RagNode;
    typedef typename RagGraph::Arc              RagArc;
    typedef typename RagGraph::OutArcIt         RagOutArcIt;
    typedef typename RagGraph::NodeIt           RagNodeIt;
    typedef typename RagGraph::EdgeIt           RagEdgeIt;
    typedef typename RagGraph::ArcIt            RagArcIt;


    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;

    typedef NodeHolder<RagGraph> PyRagNode;

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


   
    const static unsigned int RagEdgeMapDim = IntrinsicGraphShape<RagGraph>::IntrinsicEdgeMapDimension;
    const static unsigned int RagNodeMapDim = IntrinsicGraphShape<RagGraph>::IntrinsicNodeMapDimension;

    typedef NumpyArray<RagEdgeMapDim,   Singleband<float > > RagFloatEdgeArray;
    typedef NumpyArray<RagNodeMapDim,   Singleband<float > > RagFloatNodeArray;
    typedef NumpyArray<RagNodeMapDim,   Singleband<UInt32> > RagUInt32NodeArray;
    typedef NumpyArray<RagNodeMapDim,   Singleband<Int32 > > RagInt32NodeArray;
    typedef NumpyArray<RagNodeMapDim +1,Multiband <float > > RagMultiFloatNodeArray;

    typedef NumpyScalarEdgeMap<RagGraph,RagFloatEdgeArray>         RagFloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagFloatNodeArray>         RagFloatNodeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagUInt32NodeArray>        RagUInt32NodeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagInt32NodeArray>         RagInt32NodeArrayMap;
    typedef NumpyMultibandNodeMap<RagGraph,RagMultiFloatNodeArray> RagMultiFloatNodeArrayMap;



    typedef typename RagGraph:: template EdgeMap< std::vector<Edge> > RagAffiliatedEdges;


    typedef typename GraphDescriptorToMultiArrayIndex<Graph>::IntrinsicNodeMapShape NodeCoordinate;
    typedef NumpyArray<1,NodeCoordinate>  NodeCoorinateArray;

    LemonGraphRagVisitor(const std::string clsName)
    :clsName_(clsName){

    }

    void exportRagAffiliatedEdges()const{

        const std::string hyperEdgeMapNamClsName = clsName_ + std::string("RagAffiliatedEdges");
        python::class_<RagAffiliatedEdges>(hyperEdgeMapNamClsName.c_str(),python::init<const RagGraph &>())
        ;

    }

    template <class classT>
    void visit(classT& c) const
    {   

        // something like RagEdgeMap< std::vector< Edge > >
        exportRagAffiliatedEdges();

        // make the region adjacency graph
        python::def("_regionAdjacencyGraph",registerConverters(&pyMakeRegionAdjacencyGraph),
            python::return_value_policy<  python::manage_new_object >()
        );



        // on the fly rag edge mean 
        {


            typedef OnTheFlyEdgeMap2<
                Graph, typename PyNodeMapTraits<Graph,float>::Map,
                MeanFunctor<float>, float
            > ImplicitEdgeMap;

            
            python::def("_ragEdgeFeatures",
                registerConverters(
                    &pyRagEdgeMeanFromImplicit< float, float, ImplicitEdgeMap >
                ),
                (
                    python::arg("rag"),
                    python::arg("graph"),
                    python::arg("affiliatedEdges"),
                    python::arg("edgeFeatures"),
                    python::arg("out")=python::object()
                )
            );

        }

        



        // explicit rag features
        python::def("_ragEdgeFeatures",registerConverters(&pyRagEdgeFeatures<Singleband<float> >),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("affiliatedEdges"),
                python::arg("edgeFeatures"),
                python::arg("edgeSizes"),
                python::arg("acc"),
                python::arg("out")=python::object()
            )
        );

        python::def("_ragFindEdges",registerConverters(&pyRagFindEdges<Singleband<float> >),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("affiliatedEdges"),
                python::arg("labels"),
                python::arg("node")
            )
        );

        python::def("_ragNodeFeatures",registerConverters(&pyRagNodeFeaturesMultiband),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("labels"),
                python::arg("nodeFeatures"),
                python::arg("nodeSizes"),
                python::arg("acc"),
                python::arg("ignoreLabel")=-1,
                python::arg("out")=python::object()
            )
        );
        python::def("_ragNodeFeatures",registerConverters(&pyRagNodeFeaturesSingleband),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("labels"),
                python::arg("nodeFeatures"),
                python::arg("nodeSizes"),
                python::arg("acc"),
                python::arg("ignoreLabel")=-1,
                python::arg("out")=python::object()
            )
        );

        python::def("_ragNodeSize",registerConverters(&pyRagNodeSize),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("labels"),
                python::arg("ignoreLabel")=-1,
                python::arg("out")=python::object()
            )
        );
        python::def("_ragEdgeSize",registerConverters(&pyRagEdgeSize),
            (
                python::arg("rag"),
                python::arg("affiliatedEdges"),
                python::arg("out")=python::object()
            )
        );


        python::def("_ragProjectGroundTruth",registerConverters(& pyProjectGroundTruth),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("labels"),
                python::arg("gt"),
                python::arg("ragGt")=python::object(),
                python::arg("ragGtQuality")=python::object()
            )
        );

        python::def("_pyAccNodeSeeds",registerConverters(&pyAccNodeSeeds),
            (
                python::arg("rag"),
                python::arg("graph"),
                python::arg("labels"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            )
        );


        exportPyRagProjectNodeFeaturesToBaseGraph< Singleband<float > >();
        exportPyRagProjectNodeFeaturesToBaseGraph< Singleband<UInt32> >();
        exportPyRagProjectNodeFeaturesToBaseGraph< Multiband< float > >();
        exportPyRagProjectNodeFeaturesToBaseGraph< Multiband< UInt32> >();

    }


    static NumpyAnyArray pyAccNodeSeeds(
        const RagGraph &           rag,
        const Graph &              graph,
        UInt32NodeArray            labelsArray,
        UInt32NodeArray            seedsArray,
        typename PyNodeMapTraits<RagGraph, UInt32>::Array  ragSeedsArray=RagUInt32NodeArray()
    ){
        ragSeedsArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedNodeMapShape(rag));
        std::fill(ragSeedsArray.begin(),ragSeedsArray.end(),0);

        UInt32NodeArrayMap labelsArrayMap(graph,labelsArray);
        UInt32NodeArrayMap seedsArrayMap(graph,seedsArray);

        typename PyNodeMapTraits<RagGraph, UInt32>::Map ragSeedsArrayMap(rag, ragSeedsArray);


        for(NodeIt iter(graph); iter!=lemon::INVALID; ++iter){
            const UInt32 label = labelsArrayMap[*iter];
            const UInt32 seed  = seedsArrayMap[*iter];
            if(seed!=0){
                RagNode node = rag.nodeFromId(label);
                ragSeedsArrayMap[node] = seed;
            } 
        }

        return ragSeedsArray;
    }




    static python::tuple 
    pyProjectGroundTruth(
        const RagGraph &    rag,
        const Graph &       baseGraph,
        UInt32NodeArray     baseGraphRagLabels,
        UInt32NodeArray     baseGraphGt,
        RagUInt32NodeArray  ragGt,
        RagFloatNodeArray   ragGtQt
    ){

        // reshape both output arrays
        ragGt.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedNodeMapShape(rag));
        ragGtQt.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedNodeMapShape(rag));

        // make lemon maps
        UInt32NodeArrayMap baseGraphRagLabelsMap(baseGraph, baseGraphRagLabels);
        UInt32NodeArrayMap baseGraphGtMap(baseGraph, baseGraphGt);
        RagUInt32NodeArrayMap ragGtMap(rag, ragGt);
        RagFloatNodeArrayMap ragGtQtMap(rag, ragGtQt);

        // call algorithm
        projectGroundTruth(rag, baseGraph, baseGraphRagLabelsMap,
                           baseGraphGtMap, ragGtMap, ragGtQtMap);


        return python::make_tuple(ragGt, ragGtQt);
    }


    static RagAffiliatedEdges * pyMakeRegionAdjacencyGraph(
        const Graph &   graph,
        UInt32NodeArray labelsArray,
        RagGraph &      rag,
        const Int32 ignoreLabel=-1
    ){
        // numpy arrays => lemon maps
        UInt32NodeArrayMap labelsArrayMap(graph,labelsArray);

        // allocate a new RagAffiliatedEdges
        RagAffiliatedEdges * affiliatedEdges = new RagAffiliatedEdges(rag);

        // call algorithm itself
        makeRegionAdjacencyGraph(graph,labelsArrayMap,rag,*affiliatedEdges,ignoreLabel);

        return affiliatedEdges;
    }


    static RagAffiliatedEdges * pyMakeRegionAdjacencyGraphFast(
        const Graph &   graph,
        UInt32NodeArray labelsArray,
        RagGraph &      rag,
        const UInt32 maxLabel,
        const UInt32 reserveEdges
    ){
        // numpy arrays => lemon maps
        UInt32NodeArrayMap labelsArrayMap(graph,labelsArray);

        // allocate a new RagAffiliatedEdges
        RagAffiliatedEdges * affiliatedEdges = new RagAffiliatedEdges(rag);

        // call algorithm itself
        makeRegionAdjacencyGraphFast(graph,labelsArrayMap,rag,*affiliatedEdges,maxLabel,reserveEdges);

        return affiliatedEdges;
    }


    template<class T>
    static NumpyAnyArray  pyRagEdgeFeatures(
        const RagGraph &           rag,
        const Graph &              graph,
        const RagAffiliatedEdges & affiliatedEdges,
        typename PyEdgeMapTraits<Graph,T >::Array edgeFeaturesArray ,
        typename PyEdgeMapTraits<Graph,T >::Array edgeSizesArray,
        const std::string &        accumulator,
        typename PyEdgeMapTraits<RagGraph,T >::Array ragEdgeFeaturesArray
    ){

        vigra_precondition(rag.edgeNum()>=1,"rag.edgeNum()>=1 is violated");

        vigra_precondition(accumulator==std::string("mean") || accumulator==std::string("sum") || 
                           accumulator==std::string("min")  || accumulator==std::string("max"),
            "currently the accumulators are limited to mean and sum and min and max"
        );

        // resize out
        ragEdgeFeaturesArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedEdgeMapShape(rag));
        std::fill(ragEdgeFeaturesArray.begin(),ragEdgeFeaturesArray.end(),0.0f);
        // numpy arrays => lemon maps
        typename PyEdgeMapTraits<Graph   ,T >::Map edgeFeaturesArrayMap(graph,edgeFeaturesArray);
        typename PyEdgeMapTraits<Graph   ,T >::Map edgeSizesArrayMap(graph,edgeSizesArray);
        typename PyEdgeMapTraits<RagGraph,T >::Map ragEdgeFeaturesArrayMap(rag,ragEdgeFeaturesArray);


        if(accumulator == std::string("mean") ){
            for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagEdge ragEdge = *iter;
                const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
                float weightSum=0.0;
                for(size_t i=0;i<affEdges.size();++i){
                    const float weight = edgeSizesArrayMap[affEdges[i]];
                    ragEdgeFeaturesArrayMap[ragEdge]+=weight*edgeFeaturesArrayMap[affEdges[i]];
                    weightSum+=weight;
                }

                ragEdgeFeaturesArrayMap[ragEdge]/=weightSum;
            }
        }
        else if( accumulator == std::string("sum")){
            for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagEdge ragEdge = *iter;
                const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
                for(size_t i=0;i<affEdges.size();++i){
                    ragEdgeFeaturesArrayMap[ragEdge]+=edgeFeaturesArrayMap[affEdges[i]];
                }
            }
        }
        else if(accumulator == std::string("min")){
            for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagEdge ragEdge = *iter;
                const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
                float minVal=std::numeric_limits<float>::infinity();
                for(size_t i=0;i<affEdges.size();++i){
                    minVal  = std::min(minVal,edgeFeaturesArrayMap[affEdges[i]]);
                }
                ragEdgeFeaturesArrayMap[ragEdge]=minVal;
            }
        }
        else if(accumulator == std::string("max")){
            for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagEdge ragEdge = *iter;
                const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
                float maxVal=-1.0*std::numeric_limits<float>::infinity();
                for(size_t i=0;i<affEdges.size();++i){
                    maxVal  = std::max(maxVal,edgeFeaturesArrayMap[affEdges[i]]);
                }
                ragEdgeFeaturesArrayMap[ragEdge]=maxVal;
            }
        }
        else{
            throw std::runtime_error("not supported accumulator");
        }

        return ragEdgeFeaturesArray;
    }



    template<class T_PIXEL, class T, class OTF_EDGES>
    static NumpyAnyArray pyRagEdgeMeanFromImplicit(
        const RagGraph &           rag,
        const Graph &              graph,
        const RagAffiliatedEdges & affiliatedEdges,
        const OTF_EDGES & otfEdgeMap,
        typename PyEdgeMapTraits<RagGraph,T >::Array ragEdgeFeaturesArray
    ){

        // preconditions
        vigra_precondition(rag.edgeNum()>=1,"rag.edgeNum()>=1 is violated");

        // resize out
        ragEdgeFeaturesArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedEdgeMapShape(rag));
        std::fill(ragEdgeFeaturesArray.begin(),ragEdgeFeaturesArray.end(),0.0f);


        // numpy arrays => lemon maps
        typename PyEdgeMapTraits<RagGraph,T >::Map ragEdgeFeaturesArrayMap(rag,ragEdgeFeaturesArray);



        for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
            const RagEdge ragEdge = *iter;
            const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
            for(size_t i=0;i<affEdges.size();++i){
                ragEdgeFeaturesArrayMap[ragEdge]+=otfEdgeMap[affEdges[i]];
            }
            ragEdgeFeaturesArrayMap[ragEdge]/=affEdges.size();
        }
        



        // return 
        return ragEdgeFeaturesArray;

    }




    template<class T>
    static NumpyAnyArray  pyRagFindEdges(
        const RagGraph &           rag,
        const Graph &              graph,
        const RagAffiliatedEdges & affiliatedEdges,
        UInt32NodeArray            labelsArray,
        const PyRagNode &          ragNode
    ){
        UInt32NodeArrayMap   labelsArrayMap(graph,labelsArray);
        RagNode node = ragNode;
        UInt32 nodeLabel = rag.id(node);

        // Get number of points
        UInt32 nPoints = 0;
        for (RagOutArcIt iter(rag, node); iter != lemon::INVALID; ++iter) {
            const RagEdge ragEdge(*iter);
            const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
            nPoints += affEdges.size();
        }
        NumpyArray<2, UInt32> edgePoints(NumpyArray<2, UInt32>::difference_type(nPoints, NodeMapDim));

        // Find edges
        size_t nNext = 0;
        for(RagOutArcIt iter(rag, node); iter != lemon::INVALID; ++iter) {
            const RagEdge ragEdge(*iter);
            const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
            for (size_t i=0; i<affEdges.size(); ++i) {
                Node u = graph.u(affEdges[i]);
                Node v = graph.v(affEdges[i]);
                UInt32 uLabel = labelsArrayMap[u];
                UInt32 vLabel = labelsArrayMap[v];

                NodeCoordinate coords;
                if (uLabel == nodeLabel) {
                    coords = GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(graph, u);
                } else if (vLabel == nodeLabel) {
                    coords = GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(graph, v);
                } else {
                    // If you get here, then there's an error. Maybe print a message?
                }
                for(size_t k=0; k<coords.size(); ++k) {
                    edgePoints(nNext, k) = coords[k];
                }
                nNext++;
            }
        }
        return edgePoints;
    }



    static NumpyAnyArray  pyRagNodeFeaturesSingleband(
        const RagGraph &           rag,
        const Graph &              graph,
        UInt32NodeArray            labelsArray,
        FloatNodeArray             nodeFeaturesArray,
        FloatNodeArray             nodeSizesArray,
        const std::string &        accumulator,
        const Int32                ignoreLabel=-1,
        RagFloatNodeArray          ragNodeFeaturesArray=RagFloatNodeArray()
    ){

        vigra_precondition(accumulator==std::string("mean") || accumulator==std::string("sum") || 
                           accumulator==std::string("min")  || accumulator==std::string("max"),
            "currently the accumulators are limited to mean and sum and min and max "
        );

        // resize out

        ragNodeFeaturesArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedNodeMapShape(rag));
        std::fill(ragNodeFeaturesArray.begin(),ragNodeFeaturesArray.end(),0.0f);

        // numpy arrays => lemon maps
        UInt32NodeArrayMap   labelsArrayMap(graph,labelsArray);
        FloatNodeArrayMap    nodeFeaturesArrayMap(graph,nodeFeaturesArray);
        FloatNodeArrayMap    nodeSizesArrayMap(graph,nodeSizesArray);
        RagFloatNodeArrayMap ragNodeFeaturesArrayMap(rag,ragNodeFeaturesArray);

        if(accumulator == std::string("mean")){
            typename RagGraph:: template NodeMap<float> counting(rag,0.0f);
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const float  weight = nodeSizesArrayMap[*iter];
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]+= weight*nodeFeaturesArrayMap[*iter];
                    counting[ragNode]+=weight;
                }
            }
            for(RagNodeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagNode ragNode   = *iter;
                ragNodeFeaturesArrayMap[ragNode]/=counting[ragNode];
            }
        }
        else if(accumulator == std::string("sum")){
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]+=nodeFeaturesArrayMap[*iter];
                }
            }
        }
        else if(accumulator == std::string("min")){
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]=std::numeric_limits<float>::infinity();
                }
            }
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]=std::min(nodeFeaturesArrayMap[*iter],ragNodeFeaturesArrayMap[ragNode]);
                }
            }
        }
        else if(accumulator == std::string("max")){
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]= -1.0*std::numeric_limits<float>::infinity();
                }
            }
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]=std::max(nodeFeaturesArrayMap[*iter],ragNodeFeaturesArrayMap[ragNode]);
                }
            }
        }
        else{
           
        }
        return ragNodeFeaturesArray;
    }


    static NumpyAnyArray  pyRagNodeFeaturesMultiband(
        const RagGraph &           rag,
        const Graph &              graph,
        UInt32NodeArray            labelsArray,
        MultiFloatNodeArray        nodeFeaturesArray,
        FloatNodeArray             nodeSizesArray,
        const std::string &        accumulator,
        const Int32                ignoreLabel=-1,
        RagMultiFloatNodeArray     ragNodeFeaturesArray=RagMultiFloatNodeArray()
    ){
        vigra_precondition(accumulator==std::string("mean") || accumulator==std::string("sum"),
            "currently the accumulators are limited to mean and sum"
        );

        // resize out
        typename MultiArray<RagNodeMapDim+1,int>::difference_type outShape;
        for(size_t d=0;d<RagNodeMapDim;++d){
            outShape[d]=IntrinsicGraphShape<RagGraph>::intrinsicNodeMapShape(rag)[d];
        }
        outShape[RagNodeMapDim]=nodeFeaturesArray.shape(NodeMapDim);

        ragNodeFeaturesArray.reshapeIfEmpty(   RagMultiFloatNodeArray::ArrayTraits::taggedShape(outShape,"nc") );
        std::fill(ragNodeFeaturesArray.begin(),ragNodeFeaturesArray.end(),0.0f);

        // numpy arrays => lemon maps
        UInt32NodeArrayMap        labelsArrayMap(graph,labelsArray);
        MultiFloatNodeArrayMap    nodeFeaturesArrayMap(graph,nodeFeaturesArray);
        FloatNodeArrayMap         nodeSizesArrayMap(graph,nodeSizesArray);
        RagMultiFloatNodeArrayMap ragNodeFeaturesArrayMap(rag,ragNodeFeaturesArray);

        if(accumulator == std::string("mean")){
            typename RagGraph:: template NodeMap<float> counting(rag,0.0f);
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const float weight = nodeSizesArrayMap[*iter];
                    const RagNode ragNode   = rag.nodeFromId(l);
                    typename MultiFloatNodeArrayMap::Value feat = nodeFeaturesArrayMap[*iter];
                    feat*=weight;
                    ragNodeFeaturesArrayMap[ragNode]+=feat;
                    counting[ragNode]+=weight;
                }
            }
            for(RagNodeIt iter(rag);iter!=lemon::INVALID;++iter){
                const RagNode ragNode   = *iter;
                ragNodeFeaturesArrayMap[ragNode]/=counting[ragNode];
            }
        }
        else if(accumulator == std::string("sum")){
            for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
                UInt32 l = labelsArrayMap[*iter];
                if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                    const RagNode ragNode   = rag.nodeFromId(l);
                    ragNodeFeaturesArrayMap[ragNode]+=nodeFeaturesArrayMap[*iter];
                }
            }
        }
        else{
            throw std::runtime_error("for multiband only mean and sum is implemented");
        }
        return ragNodeFeaturesArray;
    }

    static NumpyAnyArray  pyRagNodeSize(
        const RagGraph &           rag,
        const Graph &              graph,
        UInt32NodeArray            labelsArray,
        const Int32                ignoreLabel=-1,
        RagFloatNodeArray          ragNodeSizeArray=RagFloatNodeArray()
    ){
        // resize out
        ragNodeSizeArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedNodeMapShape(rag));
        std::fill(ragNodeSizeArray.begin(),ragNodeSizeArray.end(),0.0f);

        // numpy arrays => lemon maps
        UInt32NodeArrayMap labelsArrayMap(graph,labelsArray);
        RagFloatNodeArrayMap ragNodeSizeArrayMap(rag,ragNodeSizeArray);
        for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
            UInt32 l = labelsArrayMap[*iter];
            if(ignoreLabel==-1 || static_cast<Int32>(l)!=ignoreLabel){
                const RagNode ragNode   = rag.nodeFromId(l);
                ragNodeSizeArrayMap[ragNode]+=1.0f;
            }
        }

        return ragNodeSizeArray;
    }

    static NumpyAnyArray  pyRagEdgeSize(
        const RagGraph &           rag,
        const RagAffiliatedEdges & affiliatedEdges,
        RagFloatEdgeArray          ragEdgeFeaturesArray
    ){
        // reshape out
        ragEdgeFeaturesArray.reshapeIfEmpty(TaggedGraphShape<RagGraph>::taggedEdgeMapShape(rag));
        // numpy arrays => lemon maps
        RagFloatEdgeArrayMap ragEdgeFeaturesArrayMap(rag,ragEdgeFeaturesArray);

        for(RagEdgeIt iter(rag);iter!=lemon::INVALID;++iter){
            const RagEdge ragEdge = *iter;
            const std::vector<Edge> & affEdges = affiliatedEdges[ragEdge];
            ragEdgeFeaturesArrayMap[ragEdge]=static_cast<float>(affEdges.size());
        }
        return ragEdgeFeaturesArray;
    }



    template<class T>
    static void exportPyRagProjectNodeFeaturesToBaseGraph(){
        python::def("_ragProjectNodeFeaturesToBaseGraph",
            registerConverters(  &pyRagProjectNodeFeaturesToBaseGraph< T > ),
            (
                python::arg("rag"),
                python::arg("baseGraph"),
                python::arg("baseGraphLabels"),
                python::arg("ragNodeFeatures"),
                python::arg("ignoreLabel")=-1,
                python::arg("out")=python::object()
            )
        );
    }

    template<class T>
    static NumpyAnyArray pyRagProjectNodeFeaturesToBaseGraph(
        const RagGraph &                                         rag,
        const Graph    &                                         graph,
        const typename PyNodeMapTraits<Graph,   UInt32>::Array & labelsWhichGeneratedRagArray,
        const typename PyNodeMapTraits<RagGraph,T     >::Array & ragNodeFeaturesArray,
        const Int32                                              ignoreLabel=-1,
        typename PyNodeMapTraits<Graph,T>::Array                 graphNodeFeaturesArray=(typename PyNodeMapTraits<Graph,T>::Array()) // out
    ){

        TaggedShape ragNodeFeaturesArrayShape = ragNodeFeaturesArray.taggedShape();
        TaggedShape graphNodeFeaturesArrayShape = TaggedGraphShape<Graph>::taggedNodeMapShape(graph);
        if(ragNodeFeaturesArrayShape.hasChannelAxis()){
            graphNodeFeaturesArrayShape.setChannelCount(ragNodeFeaturesArrayShape.channelCount());
        }
        graphNodeFeaturesArray.reshapeIfEmpty(graphNodeFeaturesArrayShape);


        // reshape out  ( last argument (out) will be reshaped if empty, and #channels is taken from second argument)
        //reshapeNodeMapIfEmpty(graph,ragNodeFeaturesArray,graphNodeFeaturesArray);
        // numpy arrays => lemon maps 
        typename PyNodeMapTraits<Graph,   UInt32>::Map labelsWhichGeneratedRagArrayMap(graph, labelsWhichGeneratedRagArray);
        typename PyNodeMapTraits<RagGraph,T     >::Map ragNodeFeaturesArrayMap(rag,ragNodeFeaturesArray);
        typename PyNodeMapTraits<Graph,   T     >::Map graphNodeFeaturesArrayMap(graph,graphNodeFeaturesArray);
        

        projectBack(rag, graph, ignoreLabel, labelsWhichGeneratedRagArrayMap, 
                    ragNodeFeaturesArrayMap, graphNodeFeaturesArrayMap);


        /*

        // run algorithm
        for(typename Graph::NodeIt iter(graph);iter!=lemon::INVALID;++iter){
            if(ignoreLabel==-1 || static_cast<Int32>(labelsWhichGeneratedRagArrayMap[*iter])!=ignoreLabel)
                graphNodeFeaturesArrayMap[*iter]=ragNodeFeaturesArrayMap[rag.nodeFromId(labelsWhichGeneratedRagArrayMap[*iter])];
            else{
                // do nothing
            }
        }
        */
        return graphNodeFeaturesArray; // out
        
    }

private:
    std::string clsName_;
};






} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_RAG_VISITOR_HXX
