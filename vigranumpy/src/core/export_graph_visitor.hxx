#ifndef EXPORT_GRAPH_VISITOR_HXX
#define EXPORT_GRAPH_VISITOR_HXX
//#define NO_IMPORT_ARRAY

/*boost python before anything else*/
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
/*std*/
#include <sstream>
#include <string>

/*vigra*/
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/graphs.hxx>
#include <vigra/graph_maps.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/python_graph.hxx>


namespace python = boost::python;

namespace vigra{




template<class GRAPH>
class LemonUndirectedGraphCoreVisitor 
:   public boost::python::def_visitor<LemonUndirectedGraphCoreVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;
    typedef LemonUndirectedGraphCoreVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;

    typedef typename Graph::NodeIt              NodeIt;
    typedef typename Graph::EdgeIt              EdgeIt;
    typedef typename Graph::ArcIt               ArcIt;
    //typedef typename Graph::NeighborNodeIt      NeighborNodeIt;
    //typedef EdgeIteratorHolder<Graph> EdgeIteratorHolderType;
    //typedef NodeIteratorHolder<Graph> NodeIteratorHolderType;
    //typedef ArcIteratorHolder<Graph>  ArcIteratorHolderType;

    typedef NodeIteratorHolder<Graph>          PyNodeIteratorHolder;
    typedef EdgeIteratorHolder<Graph>          PyEdgeIteratorHolder;

    typedef NeighbourNodeIteratorHolder<Graph> PyNeighbourNodeIteratorHolder;
    typedef IncEdgeIteratorHolder<Graph>       PyIncEdgeIteratorHolder;

    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;

    typedef std::vector<PyEdge> PyEdgeVector;
    typedef std::vector<PyNode> PyNodeVector;
    typedef std::vector<PyArc > PyArcVector;


        

    LemonUndirectedGraphCoreVisitor(const std::string clsName)
    :clsName_(clsName){

    }
    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {
        
        //import_vigranumpy();

        const std::string edgeHolderClsName = std::string("Edge")+clsName_;
        const std::string nodeHolderClsName = std::string("Node")+clsName_;
        const std::string arcHolderClsName  = std::string("Arc") +clsName_;


        const std::string edgeHolderVectorClsName = std::string("EdgeVector")+clsName_;
        const std::string nodeHolderVectorClsName = std::string("NodeVector")+clsName_;
        const std::string  arcHolderVectorClsName = std::string("ArcVector") +clsName_;

        python::class_<PyEdge>(edgeHolderClsName.c_str(),python::init< >())
        .add_property("id",  &PyEdge::id )
        .add_property("u",  &PyEdge::u )
        .add_property("v",  &PyEdge::v )
        .def("__eq__",&eqToInvalid<PyEdge>)
        .def("__ne__",&neqToInvalid<PyEdge>)
        .def("coord",&PyEdge::intrinsicEdgeCoordinate)
        ;

        python::class_<PyNode>(nodeHolderClsName.c_str(),python::init< >())
        .add_property("id",  &PyNode::id )
        .def("__eq__",&eqToInvalid<PyNode>)
        .def("__ne__",&neqToInvalid<PyNode>)
        .def("coord",&PyNode::intrinsicNodeCoordinate)
        ;

        python::class_< PyArc>( arcHolderClsName.c_str(),python::init< >())
        .add_property("id",  &PyArc::id )
        .def("__eq__",&eqToInvalid<PyArc>)
        .def("__ne__",&neqToInvalid<PyArc>)
        .def("coord",&PyArc::intrinsicArcCoordinate)
        ;


        python::class_<PyEdgeVector>(edgeHolderVectorClsName.c_str(),python::init<>())
        .def(python::vector_indexing_suite< PyEdgeVector >() );
        ;   




        const std::string nodeIteratorHolderClsName = std::string("NodeIteratorHolder")+clsName_;
        python::class_<PyNodeIteratorHolder>(nodeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&PyNodeIteratorHolder::begin,&PyNodeIteratorHolder::end))
        ;

        const std::string edgeIteratorHolderClsName = std::string("EdgeIteratorHolder")+clsName_;
        python::class_<PyEdgeIteratorHolder>(edgeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&PyEdgeIteratorHolder::begin,&PyEdgeIteratorHolder::end))
        ;

        const std::string neighbourNodeIteratorHolderClsName = std::string("NeighbourNodeIteratorHolder")+clsName_;
        python::class_<PyNeighbourNodeIteratorHolder>(neighbourNodeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&PyNeighbourNodeIteratorHolder::begin,&PyNeighbourNodeIteratorHolder::end))
        ;

        const std::string incEdgeIteratorHolderClsName = std::string("IncEdgeIteratorHolder")+clsName_;
        python::class_<PyIncEdgeIteratorHolder>(incEdgeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&PyIncEdgeIteratorHolder::begin,&PyIncEdgeIteratorHolder::end))
        ;
        //const std::string nodeIteratorHolderClsName = std::string("NodeIteratorHolder")+clsName_;
        //python::class_<NodeIteratorHolderType>(nodeIteratorHolderClsName.c_str(),python::no_init)
        //.def("__iter__",python::range(&NodeIteratorHolderType::begin,&NodeIteratorHolderType::end))
        ;

        c

            //special members
            // print graph
            .def("__str__",&asStr)
            .def("__len__",&Graph::edgeNum)


            // basic properties
            .add_property("nodeNum",  &Graph::nodeNum ,"number of nodes within the graph")
            .add_property("edgeNum",  &Graph::edgeNum ,"number of edges within the graph")
            .add_property("arcNum",   &Graph::arcNum  ,"number of arcs (2*edgeNum for undirected graphs)")
            .add_property("maxNodeId",&Graph::maxNodeId,"maximum id of a valid edge in the graph")
            .add_property("maxEdgeId",&Graph::maxEdgeId,"maximum id of a valid node in the graph")
            .add_property("maxArcId" ,&Graph::maxNodeId,"maximum id of a valid arc in the graph")

            // id functions
            .def("id",&nodeId, "get the id of a given node")
            .def("id",&edgeId, "get the id of a given edge")
            .def("id",&arcId , "get the id of a given arc")

            // item from id
            .def("nodeFromId",&nodeFromId,
                (
                    python::arg("id")
                ),
                "get the node descriptor from the given id"
            )
            .def("edgeFromId",&edgeFromId,"get the edge descriptor from the given id")//,python::with_custodian_and_ward_postcall<0,1>())
            .def("arcFromId", &arcFromId ,"get the arc descriptor from the given id")

            // find edges
            .def("findEdge",&findEdge,"find an edge between node u and v")
            .def("findEdge",&findEdgeFromIds,"find the edge between two nodes given their id")
            //  uv source target
            .def("u",&u,"get the u node of an edge")
            .def("v",&v,"geht the v node of an edge")
            .def("uId",&uId)
            .def("vId",&vId)
            .def("uvId",&uvId)
            .def("uvId",&uvIdFromId)

            .def("source",&source)
            .def("target",&target)


            // iterators
            //.def("edgeIter",&edgeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted
            //.def("nodeIter",&nodeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted

            .def("nodeIter",&getNodeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())
            .def("edgeIter",&getEdgeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())
            .def("neighbourNodeIter",&getNeighbourNodeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())
            .def("incEdgeIter",      &getIncEdgeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())

            // intrinsic shape of maps
            .def("intrinsicNodeMapShape",&IntrinsicGraphShape<Graph>::intrinsicNodeMapShape)
            .def("intrinsicEdgeMapShape",&IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape)
            .def("intrinsicArcMapShape" , &IntrinsicGraphShape<Graph>::intrinsicArcMapShape)
            // commented out because TaggedShape is currently not exported to Python
            //.def("taggedNodeMapShape",&TaggedGraphShape<Graph>::taggedNodeMapShape)
            //.def("taggedEdgeMapShape",&TaggedGraphShape<Graph>::taggedEdgeMapShape)
            //.def("taggedArcMapShape" , &TaggedGraphShape<Graph>::taggedArcMapShape)
            .def("axistagsNodeMap",&TaggedGraphShape<Graph>::axistagsNodeMap)
            .def("axistagsEdgeMap",&TaggedGraphShape<Graph>::axistagsEdgeMap)
            .def("axistagsArcMap" , &TaggedGraphShape<Graph>::axistagsArcMap)

            // intrinsic coordinate of node/edge/arc
            .def("intrinsicNodeCoordinate",& GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate)
            .def("intrinsicEdgeCoordinate",&GraphDescriptorToMultiArrayIndex<Graph>::intrinsicEdgeCoordinate)
            .def("intrinsicArcCoordinate",& GraphDescriptorToMultiArrayIndex<Graph>::intrinsicArcCoordinate )


            ///////////////////////////////////
            // vectorized  api
            ///////////////////////////////////
            .def("nodeIds",registerConverters(&itemIds<Node,NodeIt>),( python::arg("out")=python::object() ) )
            .def("edgeIds",registerConverters(&itemIds<Edge,EdgeIt>),( python::arg("out")=python::object() ) )
            .def("arcIds" ,registerConverters(&itemIds<Arc ,ArcIt >),( python::arg("out")=python::object() ) )

            .def("nodeIdMap",registerConverters(&nodeIdMap),( python::arg("out")=python::object() ) )

            .def("findEdges",registerConverters(&findEdges),( python::arg("nodeIdPairs"), python::arg("out")=python::object() ) )

            .def("uIds" ,registerConverters(&uIds), ( python::arg("out")=python::object() ) )
            .def("vIds" ,registerConverters(&vIds), ( python::arg("out")=python::object() ) )
            .def("uvIds",registerConverters(&uvIds),( python::arg("out")=python::object() ) )
            .def("uIds",registerConverters(&uIdsSubset),( python::arg("edgeIds"),python::arg("out")=python::object() ) )
            .def("vIds",registerConverters(&vIdsSubset),( python::arg("edgeIds"),python::arg("out")=python::object() ) )
            .def("uvIds",registerConverters(&uvIdsSubset),( python::arg("edgeIds"),python::arg("out")=python::object() ) )
            


            // these functions are defined on the "FULL SET"
            .def("validEdgeIds",registerConverters(&validIds<Edge,EdgeIt>),( python::arg("out")=python::object() ) )
            .def("validNodeIds",registerConverters(&validIds<Node,NodeIt>),( python::arg("out")=python::object() ) )
            .def("validArcIds" ,registerConverters(&validIds<Arc ,ArcIt >),( python::arg("out")=python::object() ) )




            //.def("dtypetest",registerConverters(&dtypetest<Edge,EdgeIt>),( python::arg("out")=python::object() ) )
            
        ;
    }

    static PyEdgeIteratorHolder getEdgeIteratorHolder(const Graph & self){
        return PyEdgeIteratorHolder(self);
    }

    static PyNodeIteratorHolder getNodeIteratorHolder(const Graph & self){
        return PyNodeIteratorHolder(self);
    }


    static PyIncEdgeIteratorHolder getIncEdgeIteratorHolder(const Graph & self,const PyNode & node){
        return PyIncEdgeIteratorHolder(self,node);
    }


    static PyNeighbourNodeIteratorHolder getNeighbourNodeIteratorHolder(const Graph & self,const PyNode & node){
        return PyNeighbourNodeIteratorHolder(self,node);
    }

    template<class ITEM,class ITEM_IT>
    static NumpyAnyArray dtypetest(   
        const Graph & g, 
        NumpyArray<1,ITEM> out =NumpyArray<1,ITEM>() 
    ){
        typedef GraphItemHelper<Graph,ITEM> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,ITEM>::difference_type(  ItemHelper::itemNum(g)  ));
        size_t  counter=0;
        for(ITEM_IT i(g);i!=lemon::INVALID;++i){
            const ITEM item = *i;
            out(counter)=item;
            ++counter;
        }
        return out;
    }


    static NumpyAnyArray findEdges(
        const Graph & g,
        NumpyArray<2,UInt32> nodeIdPairs, 
        NumpyArray<1,Int32> out =NumpyArray<1,Int32>()  
    ){
        out.reshapeIfEmpty(typename NumpyArray<1,Int32>::difference_type(  nodeIdPairs.shape(0)  ));
        for(MultiArrayIndex i=0; i<nodeIdPairs.shape(0); ++i){
            const Edge e = g.findEdge(
                g.nodeFromId(nodeIdPairs(i,0)),
                g.nodeFromId(nodeIdPairs(i,1))
            );
            out(i) = e==lemon::INVALID ? -1 : g.id(e);
        }
       
        return out;
    }


    static std::string asStr(const Graph &g){
        std::stringstream ss;
        ss<<"Nodes: "<<g.nodeNum()<<" Edges: "<<g.edgeNum()<<" maxNodeId: "<<g.maxNodeId()<<" maxEdgeId: "<<g.maxEdgeId();
        return ss.str();
    }


    static NumpyAnyArray uIds(const Graph & g, NumpyArray<1,UInt32> out =NumpyArray<1,UInt32>() ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  ItemHelper::itemNum(g)  ));
        size_t  counter=0;
        for(EdgeIt i(g);i!=lemon::INVALID;++i){
            out(counter)=g.id(g.u(*i));
            ++counter;
        }
        return out;
    }
    static NumpyAnyArray vIds(const Graph & g, NumpyArray<1,UInt32> out =NumpyArray<1,UInt32>() ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  ItemHelper::itemNum(g)  ));
        size_t  counter=0;
        for(EdgeIt i(g);i!=lemon::INVALID;++i){
            out(counter)=g.id(g.v(*i));
            ++counter;
        }
        return out;
    }

    static NumpyAnyArray uvIds(const Graph & g, NumpyArray<2,UInt32> out =NumpyArray<1,UInt32>() ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<2,UInt32>::difference_type(  ItemHelper::itemNum(g) ,2 ));
        size_t  counter=0;
        for(EdgeIt i(g);i!=lemon::INVALID;++i){
            out(counter,0)=g.id(g.u(*i));
            out(counter,1)=g.id(g.v(*i));
            ++counter;
        }
        return out;
    }

    static NumpyAnyArray uIdsSubset(
        const Graph & g, 
        NumpyArray<1,UInt32> edgeIds,
        NumpyArray<1,UInt32> out =NumpyArray<1,UInt32>() 
    ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  edgeIds.shape(0)));
        for(MultiArrayIndex i=0; i<edgeIds.shape(0); ++i){
            const index_type edgeId=edgeIds(i);
            const Edge edge  = g.edgeFromId(edgeId);
            if(edge!=lemon::INVALID){
                out(i)=g.id(g.u(edge));
            }
        }
        return out;
    }
    static NumpyAnyArray vIdsSubset(
        const Graph & g, 
        NumpyArray<1,UInt32> edgeIds,
        NumpyArray<1,UInt32> out =NumpyArray<1,UInt32>() 
    ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  edgeIds.shape(0)));
        for(MultiArrayIndex i=0; i<edgeIds.shape(0); ++i){
            const index_type edgeId=edgeIds(i);
            const Edge edge  = g.edgeFromId(edgeId);
            if(edge!=lemon::INVALID){
                out(i)=g.id(g.v(edge));
            }
        }
        return out;
    }

    static NumpyAnyArray uvIdsSubset(
        const Graph & g, 
        NumpyArray<1,UInt32> edgeIds,
        NumpyArray<2,UInt32> out =NumpyArray<1,UInt32>() 
    ){
        typedef GraphItemHelper<Graph,Edge> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<2,UInt32>::difference_type(  edgeIds.shape(0) ,2 ));
        for(MultiArrayIndex i=0; i<edgeIds.shape(0); ++i){
            const index_type edgeId=edgeIds(i);
            const Edge edge  = g.edgeFromId(edgeId);
            if(edge!=lemon::INVALID){
                out(i,0)=g.id(g.u(edge));
                out(i,1)=g.id(g.v(edge));
            }
        }
        return out;
    }

    template<class ITEM,class ITEM_IT>
    static NumpyAnyArray validIds(const Graph & g, NumpyArray<1,bool> out =NumpyArray<1,UInt32>() ){
        typedef GraphItemHelper<Graph,ITEM> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  ItemHelper::maxItemId(g)  ));
        std::fill(out.begin(),out.end(),false);
        size_t  counter=0;
        for(ITEM_IT i(g);i!=lemon::INVALID;++i){
            out(g.id(*i))=true;
            ++counter;
        }
        return out;
    }

    template<class ITEM,class ITEM_IT>
    static NumpyAnyArray itemIds(const Graph & g, NumpyArray<1,UInt32> out =NumpyArray<1,UInt32>() ){
        typedef GraphItemHelper<Graph,ITEM> ItemHelper;
        out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(  ItemHelper::itemNum(g)  ));
        size_t  counter=0;
        for(ITEM_IT i(g);i!=lemon::INVALID;++i){
            out(counter)=g.id(*i);
            ++counter;
        }
        return out;
    }


    static  python::tuple uvId(const Graph & self,const PyEdge & e){
        return  python::make_tuple(self.id(self.u(e)),self.id(self.u(e)) );
    }

    static  python::tuple uvIdFromId(const Graph & self,const index_type i){
        const Edge e = self.edgeFromId(i);
        return  python::make_tuple(self.id(self.u(e)),self.id(self.u(e)) );
    }

    static  index_type uId(const Graph & self,const PyEdge & e){
        return self.id(self.u(e));
    }
    static  index_type vId(const Graph & self,const PyEdge & e){
        return self.id(self.v(e));
    }
    static  PyNode u(const Graph & self,const PyEdge & e){
        return PyNode(self,self.u(e));
    }
    static  PyNode v(const Graph & self,const PyEdge & e){
        return PyNode(self,self.v(e));
    }
    static  PyNode source(const Graph & self,const PyArc & a){
        return PyNode(self,self.source(a));
    }
    static  PyNode target(const Graph & self,const PyArc & a){
        return PyNode(self,self.target(a));
    }

    template<class ITEM>
    static bool eqToInvalid(const ITEM &  item,const lemon::Invalid iv){
        return item.graph_==NULL || item==lemon::INVALID;
    }

    template<class ITEM>
    static bool neqToInvalid(const ITEM &  item,const lemon::Invalid iv){
        return item.graph_!=NULL && item!=lemon::INVALID;
    }
    
    static PyNode nodeFromId(const Graph & self,const index_type id){
        return PyNode(self,self.nodeFromId(id));
    }
    static PyEdge edgeFromId(const Graph & self,const index_type id){
        return PyEdge(self,self.edgeFromId(id));
    }
    static PyArc   arcFromId(const Graph & self,const index_type id){
        return PyArc(self,self.arcFromId(id));
    }


    static PyEdge findEdge( const Graph & self ,const PyNode & u , const PyNode & v){
        return PyEdge(self,self.findEdge(u,v));
    }

    static PyEdge findEdgeFromIds( const Graph & self ,const index_type u , const index_type v){
        return PyEdge(self,self.findEdge(self.nodeFromId(u),self.nodeFromId(v)));
    }

    static index_type nodeId( const Graph & self,const PyNode & node ){return  self.id(node);}
    static index_type edgeId( const Graph & self,const PyEdge & edge ){return  self.id(edge);}
    static index_type arcId(  const Graph & self,const PyArc & arc ){return  self.id(arc);}



    static NumpyAnyArray nodeIdMap(
        const Graph & graph,
        typename PyNodeMapTraits<Graph,   UInt32>::Array  idArray
    ){  
        //reshape output
        idArray.reshapeIfEmpty(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(graph));

        // array to lemon map
        typename PyNodeMapTraits<Graph,   UInt32>::Map idArrayMap(graph, idArray);

        for(NodeIt iter(graph);iter!=lemon::INVALID;++iter){
            idArrayMap[*iter]=graph.id(*iter);
        }

        return idArray;
    }

};


template<class GRAPH>
class LemonUndirectedGraphAddItemsVisitor 
:   public boost::python::def_visitor<LemonUndirectedGraphAddItemsVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonUndirectedGraphAddItemsVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;

    typedef typename Graph::NodeIt              NodeIt;
    typedef typename Graph::EdgeIt              EdgeIt;
    typedef typename Graph::ArcIt               ArcIt;
    //typedef typename Graph::NeighborNodeIt      NeighborNodeIt;



    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;

    LemonUndirectedGraphAddItemsVisitor(const std::string clsName)
    :clsName_(clsName){

    }
    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {
        
        //import_vigranumpy();


        c
            // add node
            .def("addNode",&addNode)
            .def("addNode",&addNodeFromId)
            .def("addEdge",&addEdge)
            .def("addEdges",registerConverters(&addEdges),
                (
                    python::arg("edges"),
                    python::arg("out")=python::object()
                )
            )
        ;
    }

    static PyNode addNode(Graph & self){
        return PyNode(self,self.addNode());
    }
    static PyNode addNodeFromId(Graph & self,const index_type id){
        return PyNode(self,self.addNode(id));
    }

    static PyEdge addEdge(Graph & self ,const PyNode & u , const PyNode & v){
        return PyEdge(self,self.addEdge(u,v));
    }

    static NumpyAnyArray addEdges(Graph & self,
        NumpyArray<2,UInt32> edges,
        NumpyArray<1,UInt32> edgeIds  =NumpyArray<1,index_type>()
    ){
        edgeIds.reshapeIfEmpty(typename NumpyArray<1,index_type>::difference_type(edges.shape(0)));
        for(MultiArrayIndex i=0; i<edges.shape(0); ++i){
            const Edge e = self.addEdge(edges(i,0),edges(i,1));
            edgeIds(i)=self.id(e);
        }
        return edgeIds;
    }
};

} // end namespace vigra

#endif // EXPORT_GRAPH_VISITOR_HXX