#ifndef EXPORT_GRAPH_VISITOR_HXX
#define EXPORT_GRAPH_VISITOR_HXX
//#define NO_IMPORT_ARRAY

/*std*/
#include <sstream>
#include <string>

/*vigra*/
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <boost/python.hpp>
#include <vigra/graphs.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/python_graph_generalization.hxx>
namespace python = boost::python;

namespace vigra{



template<class GRAPH>
struct EdgeIteratorHolder{

    typedef typename GraphIteratorAccessor<GRAPH>::EdgeIt const_iterator;
    EdgeIteratorHolder( ):graph_(NULL){}
    EdgeIteratorHolder(const GRAPH & g):graph_(&g){}
    const_iterator begin()const{
        return GraphIteratorAccessor<GRAPH>::edgesBegin(*graph_);
    }
    const_iterator end()const{
        return GraphIteratorAccessor<GRAPH>::edgesEnd(*graph_);
    }

    const GRAPH * graph_;
};


template<class GRAPH>
struct NodeIteratorHolder{

    typedef typename GraphIteratorAccessor<GRAPH>::NodeIt const_iterator;
    NodeIteratorHolder( ):graph_(NULL){}
    NodeIteratorHolder(const GRAPH & g):graph_(&g){}
    const_iterator begin()const{
        return GraphIteratorAccessor<GRAPH>::nodesBegin(*graph_);
    }
    const_iterator end()const{
        return GraphIteratorAccessor<GRAPH>::nodesBegin(*graph_);
    }
    const GRAPH * graph_;
};

template<class GRAPH>
struct ArcIteratorHolder{

    typedef typename GraphIteratorAccessor<GRAPH>::ArcIt const_iterator;
    ArcIteratorHolder( ):graph_(NULL){}
    ArcIteratorHolder(const GRAPH & g):graph_(&g){}
    const_iterator begin()const{
        return GraphIteratorAccessor<GRAPH>::arcsBegin(*graph_);
    }
    const_iterator end()const{
        return GraphIteratorAccessor<GRAPH>::arcsEnd(*graph_);
    }
    const GRAPH * graph_;
};


template<class GRAPH>
class LemonDirectedGraphCoreVisitor 
:   public boost::python::def_visitor<LemonDirectedGraphCoreVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;
    typedef LemonDirectedGraphCoreVisitor<GRAPH> VisitorType;
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


    typedef NeighbourNodeIteratorHolder<Graph> PyNeighbourNodeIteratorHolder;
    typedef IncEdgeIteratorHolder<Graph>       PyIncEdgeIteratorHolder;

    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;

    LemonDirectedGraphCoreVisitor(const std::string clsName)
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
        ;
        python::class_< PyArc>( arcHolderClsName.c_str(),python::init< >())
        .add_property("id",  &PyArc::id )
        .def("__eq__",&eqToInvalid<PyArc>)
        .def("__ne__",&neqToInvalid<PyArc>)
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
            .def("nodeFromId",&nodeFromId,"get the node descriptor from the given id")
            .def("edgeFromId",&edgeFromId,"get the edge descriptor from the given id")//,python::with_custodian_and_ward_postcall<0,1>())
            .def("arcFromId", &arcFromId ,"get the arc descriptor from the given id")

            // find edges
            .def("findEdge",&findEdge)
            .def("findEdge",&findEdgeFromIds)
            //  uv source target
            .def("u",&u)
            .def("v",&v)
            .def("uId",&uId)
            .def("vId",&vId)
            .def("uvId",&uvId)
            .def("uvId",&uvIdFromId)

            .def("source",&source)
            .def("target",&target)


            // iterators
            //.def("edgeIter",&edgeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted
            //.def("nodeIter",&nodeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted

            .def("neighbourNodeIter",&getNeighbourNodeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())
            .def("incEdgeIter",      &getIncEdgeIteratorHolder,python::with_custodian_and_ward_postcall<0,1>())

            // intrinsic shape of maps
            .def("intrinsicNodeMapShape",&IntrinsicGraphShape<Graph>::intrinsicNodeMapShape)
            .def("intrinsicEdgeMapShape",&IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape)
            .def("intrinsicArcMapShape" , &IntrinsicGraphShape<Graph>::intrinsicArcMapShape)
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

            .def("findEdges",registerConverters(&findEdges),( python::arg("nodeIdPairs"), python::arg("out")=python::object() ) )

            .def("uIds" ,registerConverters(&uIds), ( python::arg("out")=python::object() ) )
            .def("vIds" ,registerConverters(&uIds), ( python::arg("out")=python::object() ) )
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
        for(size_t i=0;i<nodeIdPairs.shape(0);++i){
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
        for(size_t i=0;i<edgeIds.shape(0);++i){
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
        for(size_t i=0;i<edgeIds.shape(0);++i){
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
        for(size_t i=0;i<edgeIds.shape(0);++i){
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
        return self.id(self.u(e));
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
    //static EdgeIteratorHolderType edgeHolder(const Graph & self){return EdgeIteratorHolderType(self);}
    //static NodeIteratorHolderType nodeHolder(const Graph & self){return NodeIteratorHolderType(self);}

};


template<class GRAPH>
class LemonDirectedGraphAddItemsVisitor 
:   public boost::python::def_visitor<LemonDirectedGraphAddItemsVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonDirectedGraphAddItemsVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;

    typedef typename Graph::NodeIt              NodeIt;
    typedef typename Graph::EdgeIt              EdgeIt;
    typedef typename Graph::ArcIt               ArcIt;
    //typedef typename Graph::NeighborNodeIt      NeighborNodeIt;
    typedef EdgeIteratorHolder<Graph> EdgeIteratorHolderType;
    typedef NodeIteratorHolder<Graph> NodeIteratorHolderType;
    typedef ArcIteratorHolder<Graph>  ArcIteratorHolderType;


    typedef EdgeHolder<Graph> PyEdge;
    typedef NodeHolder<Graph> PyNode;
    typedef  ArcHolder<Graph> PyArc;

    LemonDirectedGraphAddItemsVisitor(const std::string clsName)
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
        for(size_t i=0;i<edges.shape(0);++i){
            const Edge e = self.addEdge(edges(i,0),edges(i,1));
            edgeIds(i)=self.id(e);
        }
        return edgeIds;
    }
};






template<class GRAPH>
class LemonGraphMapVisitor 
:   public boost::python::def_visitor<LemonGraphMapVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonGraphMapVisitor<GRAPH> VisitorType;
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

    LemonGraphMapVisitor(const std::string clsName)
    :clsName_(clsName){

    }
    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {
        
        


        //c
        //    // add node
        //    .def("addNode",&addNode)
        //    .def("addEdge",&addEdge)
        //;
    }


};




} // end namespace vigra

#endif // EXPORT_GRAPH_VISITOR_HXX