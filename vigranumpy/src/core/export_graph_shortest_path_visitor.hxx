#ifndef VIGRA_EXPORT_GRAPH_SHORTEST_PATH_VISITOR_HXX
#define VIGRA_EXPORT_GRAPH_SHORTEST_PATH_VISITOR_HXX
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
namespace python = boost::python;

namespace vigra{


template<class GRAPH,class PREDECESSORS,class IDS_ARRAY>
void pathIds(
    const GRAPH & g,
    const typename  GRAPH::Node source,
    const typename  GRAPH::Node target,
    const PREDECESSORS & predecessors,
    IDS_ARRAY & ids
){
    if(predecessors[target]!=lemon::INVALID){

        ids(0)=g.id(target);
        typename GRAPH::Node currentNode = target;
        size_t length=1;
        while(currentNode!=source){
            currentNode=predecessors[currentNode];
            ids(length)=g.id(currentNode);
            length+=1;
        }
        // reverse such ids[0]==g.id(source)
        std::reverse(ids.begin(),ids.begin()+length);
    }
}

template<class GRAPH,class PREDECESSORS,class COORDINATE_ARRAY>
void pathCoordinates(
    const GRAPH & g,
    const typename  GRAPH::Node source,
    const typename  GRAPH::Node target,
    const PREDECESSORS & predecessors,
    COORDINATE_ARRAY & coords
){
    typedef GraphDescriptorToMultiArrayIndex<GRAPH> DescToCoord;
    if(predecessors[target]!=lemon::INVALID){
        coords(0)=DescToCoord::intrinsicNodeCoordinate(g,target);
        typename GRAPH::Node currentNode = target;
        size_t length=1;
        while(currentNode!=source){
            currentNode=predecessors[currentNode];
            coords(length)=DescToCoord::intrinsicNodeCoordinate(g,currentNode);
            length+=1;
        }
        // reverse such ids[0]==DescToCoord::intrinsicNodeCoordinate(g,source);
        std::reverse(coords.begin(),coords.begin()+length);
    }
}






template<class GRAPH>
class LemonGraphShortestPathVisitor 
:   public boost::python::def_visitor<LemonGraphShortestPathVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonGraphShortestPathVisitor<GRAPH> VisitorType;
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

    LemonGraphShortestPathVisitor(const std::string clsName)
    :clsName_(clsName){

    }

    void exportShortestPathAlgorithms()const{
        // class based algorithms (with factories):
        // - ShortestPathDijkstra
        const std::string dijkstraClsName = std::string("ShortestPathDijkstra")+clsName_;
        python::class_<ShortestPathDijkstraType, boost::noncopyable >(dijkstraClsName.c_str(),python::init<const Graph &>())
        .def("run",registerConverters(&runShortestPathNoTarget),
            (
                python::arg("edgeWeights"),
                python::arg("source")
            )
        )
        .def("run",registerConverters(&runShortestPath),
            (
                python::arg("edgeWeights"),
                python::arg("source"),
                python::arg("target")
            )
        )
        .def("nodeIdPath",registerConverters(&makeNodeIdPath),
            (
                python::arg("target"),
                python::arg("out")=python::object()
            ) 
        )
        .def("nodeCoordinatePath",registerConverters(&makeNodeCoordinatePath),
            (
                python::arg("target"),
                python::arg("out")=python::object()
            ) 
        )
        .def("distance",&pyShortestPathSingleDist,
            (
                python::arg("target")
            )
        )
        .def("distances",registerConverters(&pyShortestPathDistance),
            (
                python::arg("out")=python::object()
            )
        )
        .def("predecessors",registerConverters(&pyShortestPathPredecessors),
            (
                python::arg("out")=python::object()
            )
        )
        ;

        python::def("_shortestPathDijkstra",&pyShortestPathDijkstraTypeFactory,
            python::return_value_policy<python::manage_new_object>() 
        );

    }

    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {   
        // - Dijkstra
        exportShortestPathAlgorithms();
    }

        static ShortestPathDijkstraType * pyShortestPathDijkstraTypeFactory(const Graph & g){
        return new ShortestPathDijkstraType(g);
    }

    static NumpyAnyArray pyShortestPathDistance(
        const ShortestPathDijkstraType & sp,
        FloatNodeArray distanceArray = FloatNodeArray()
    ){
        // reshape output
        distanceArray.reshapeIfEmpty(  IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(sp.graph()));

        // numpy arrays => lemon maps
        FloatNodeArrayMap distanceArrayMap(sp.graph(),distanceArray);

        copyNodeMap(sp.graph(),sp.distances(),distanceArrayMap);

        return distanceArray;
    }

    static float pyShortestPathSingleDist(
        const ShortestPathDijkstraType & sp,
        const PyNode & target
    ){
        return sp.distance(target);
    }


    static NumpyAnyArray pyShortestPathPredecessors(
        const ShortestPathDijkstraType & sp,
        Int32NodeArray predecessorsArray = FloatNodeArray()
    ){
        // reshape output
        predecessorsArray.reshapeIfEmpty(  IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(sp.graph()));

        // numpy arrays => lemon maps
        Int32NodeArrayMap predecessorsArrayMap(sp.graph(),predecessorsArray);

        for(NodeIt n(sp.graph());n!=lemon::INVALID;++n){
            const Node pred = sp.predecessors()[*n];
            predecessorsArrayMap[*n]= (pred!=lemon::INVALID ? sp.graph().id(pred) : -1);
        }
        return predecessorsArray;
    }


    static NumpyAnyArray makeNodeIdPath(
        const ShortestPathDijkstraType & sp,
        PyNode target,
        NumpyArray<1,Singleband<UInt32> > nodeIdPath = NumpyArray<1,Singleband<UInt32> >()
    ){
        typename  ShortestPathDijkstraType::PredecessorsMap predMap = sp.predecessors();
        const Node source = sp.source();
        Node currentNode = target; 
        // comput length of the path
        const size_t length = pathLength(Node(source),Node(target),predMap);
        nodeIdPath.reshapeIfEmpty(typename NumpyArray<1,Singleband<UInt32> >::difference_type(length));
        pathIds(sp.graph(),source,target,predMap,nodeIdPath);
        return nodeIdPath;
        
    }

    static NumpyAnyArray makeNodeCoordinatePath(
        const ShortestPathDijkstraType & sp,
        PyNode target,
        NodeCoorinateArray nodeCoordinates = NodeCoorinateArray()
    ){
        typename  ShortestPathDijkstraType::PredecessorsMap predMap = sp.predecessors();
        const Node source = sp.source();
        // comput length of the path
        const size_t length = pathLength(Node(source),Node(target),predMap);
        nodeCoordinates.reshapeIfEmpty(typename NumpyArray<1,Singleband<UInt32> >::difference_type(length));
        pathCoordinates(sp.graph(),source,target,predMap,nodeCoordinates);
        return nodeCoordinates;
    }

    static void runShortestPath(
        ShortestPathDijkstraType & sp,
        FloatEdgeArray edgeWeightsArray,
        PyNode source,
        PyNode target
    ){
        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(sp.graph(),edgeWeightsArray);

        // run algorithm itself
        sp.run(edgeWeightsArrayMap,source,target);
    }

    static void runShortestPathNoTarget(
        ShortestPathDijkstraType & sp,
        FloatEdgeArray edgeWeightsArray,
        PyNode source
    ){
        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(sp.graph(),edgeWeightsArray);

        // run algorithm itself
        sp.run(edgeWeightsArrayMap,source);
    }

};




} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_SHORTEST_PATH_VISITOR_HXX
