#ifndef EXPORT_GRAPH_VISITOR_HXX
#define EXPORT_GRAPH_VISITOR_HXX
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <boost/python.hpp>
#include <vigra/graphs.hxx>
#include <vigra/graph_generalization.hxx>
namespace python = boost::python;

namespace vigra{


template<class GRAPH>
struct NodeHolder{
    typedef typename GRAPH::Node Node;
    NodeHolder(const lemon::Invalid & iv = lemon::INVALID)
    : item_(lemon::INVALID),
      graph_(NULL)
    {}
    NodeHolder(const GRAPH & g , const Node & item)
    : item_(item),
      graph_(&g)
    {}

    operator Node ()const{
        return item_;
    }
    typename GRAPH::index_type id()const{
        return graph_->id(item_);
    }
    Node item_;
    const GRAPH * graph_;
};



template<class GRAPH>
struct EdgeHolder{

    typedef typename GRAPH::Edge Edge;
    EdgeHolder(const lemon::Invalid & iv = lemon::INVALID)
    : item_(lemon::INVALID),
      graph_(NULL)
    {}
    EdgeHolder(const GRAPH & g , const Edge & item)
    : item_(item),
      graph_(&g)
    {}

    operator Edge ()const{
        return item_;
    }
    typename GRAPH::index_type id()const{
        return graph_->id(item_);
    }

    NodeHolder<GRAPH> u()const{
        return NodeHolder<GRAPH>(*graph_,graph_->u(item_));
    }
    NodeHolder<GRAPH> v()const{
        return NodeHolder<GRAPH>(*graph_,graph_->v(item_));
    }

    Edge item_;
    const GRAPH * graph_; 
};



template<class GRAPH>
struct ArcHolder{
    typedef typename GRAPH::Arc Arc;
    ArcHolder(const lemon::Invalid & iv = lemon::INVALID)
    : item_(lemon::INVALID),
      graph_(NULL)
    {}
    ArcHolder(const GRAPH & g , const Arc & item)
    : item_(item),
      graph_(&g)
    {}

    operator Arc ()const{
        return item_;
    }
    typename GRAPH::index_type id()const{
        return graph_->id(item_);
    }
    Arc item_;
    const GRAPH * graph_;
};

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


template<class GRAPH,class NODE_RET_POLICY,class EDGE_RET_POLICY>
class LemonDirectedGraphCoreVisitor 
:   public boost::python::def_visitor<LemonDirectedGraphCoreVisitor<GRAPH,NODE_RET_POLICY,EDGE_RET_POLICY> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH GraphType;

    typedef LemonDirectedGraphCoreVisitor<GRAPH,NODE_RET_POLICY,EDGE_RET_POLICY> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename GraphType::index_type       index_type;
    typedef typename GraphType::Edge             Edge;
    typedef typename GraphType::Node             Node;
    typedef typename GraphType::Arc              Arc;

    typedef typename GraphType::NodeIt              NodeIt;
    typedef typename GraphType::EdgeIt              EdgeIt;
    typedef typename GraphType::ArcIt               ArcIt;
    //typedef typename GraphType::NeighborNodeIt      NeighborNodeIt;
    typedef EdgeIteratorHolder<GraphType> EdgeIteratorHolderType;
    typedef NodeIteratorHolder<GraphType> NodeIteratorHolderType;
    typedef ArcIteratorHolder<GraphType>  ArcIteratorHolderType;


    typedef EdgeHolder<GraphType> PyEdge;
    typedef NodeHolder<GraphType> PyNode;
    typedef  ArcHolder<GraphType> PyArc;

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




        const std::string edgeIteratorHolderClsName = std::string("EdgeIteratorHolder")+clsName_;
        python::class_<EdgeIteratorHolderType>(edgeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&EdgeIteratorHolderType::begin,&EdgeIteratorHolderType::end))
        ;

        const std::string nodeIteratorHolderClsName = std::string("NodeIteratorHolder")+clsName_;
        python::class_<NodeIteratorHolderType>(nodeIteratorHolderClsName.c_str(),python::no_init)
        .def("__iter__",python::range(&NodeIteratorHolderType::begin,&NodeIteratorHolderType::end))
        ;

        c
            // basic properties
            .add_property("nodeNum",  &GraphType::nodeNum )
            .add_property("edgeNum",  &GraphType::edgeNum )
            .add_property("arcNum",   &GraphType::arcNum  )
            .add_property("maxNodeId",&GraphType::maxNodeId )
            .add_property("maxEdgeId",&GraphType::maxEdgeId )
            .add_property("maxArcId" ,&GraphType::maxNodeId )

            // id functions
            .def("nodeId",&nodeId)
            .def("edgeId",&edgeId)
            .def("arcId" ,&arcId)
            .def("id",&nodeId)
            .def("id",&edgeId)
            .def("id",&arcId)

            // item from id
            .def("nodeFromId",&nodeFromId)
            .def("edgeFromId",&edgeFromId)
            .def("arcFromId", &arcFromId)

            // find edges
            .def("findEdge",&findEdge)

            //  uv source target
            .def("u",&u)
            .def("v",&v)
            .def("source",&source)
            .def("target",&target)
            // iterators
            //.def("edgeIter",&edgeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted
            //.def("nodeIter",&nodeHolder,python::with_custodian_and_ward_postcall<0,1>() )  // graph may not be deleted bevore holder is deleted


            //numpy batch interface
            .def("edgeIds",&edgeIds//,
                //(
                //    python::arg("out")=python::object()
                //)
            )
        ;
    }

    static void edgeIds(
    //static NumpyAnyArray edgeIds(
        const GraphType & g//,
        //NumpyArray<1,index_type> out = NumpyArray<1,index_type>()
    ){
        //out.reshapeIfEmpty(typename NumpyArray<1,index_type>::difference_type(g.edgeNum()));
        size_t  counter=0;
        for(EdgeIt e(g);e!=lemon::INVALID;++e){
            const Edge edge(*e);
            //out(counter)=g.id(edge);
            ++counter;
        }
        //return out;
    }

    static  PyNode u(const GraphType & self,const PyEdge & e){
        return PyNode(self,self.u(e.item_));
    }
    static  PyNode v(const GraphType & self,const PyEdge & e){
        return PyNode(self,self.v(e.item_));
    }
    static  PyNode source(const GraphType & self,const PyArc & a){
        return PyNode(self,self.source(a.item_));
    }
    static  PyNode target(const GraphType & self,const PyArc & a){
        return PyNode(self,self.target(a.item_));
    }

    template<class ITEM>
    static bool eqToInvalid(const ITEM &  item,const lemon::Invalid iv){
        return item.graph_==NULL || item.item_==lemon::INVALID;
    }

    template<class ITEM>
    static bool neqToInvalid(const ITEM &  item,const lemon::Invalid iv){
        return item.graph_!=NULL && item.item_!=lemon::INVALID;
    }
    
    static PyNode nodeFromId(const GraphType & self,const index_type id){
        return PyNode(self,self.nodeFromId(id));
    }
    static PyEdge edgeFromId(const GraphType & self,const index_type id){
        return PyEdge(self,self.edgeFromId(id));
    }
    static PyArc   arcFromId(const GraphType & self,const index_type id){
        return PyArc(self,self.arcFromId(id));
    }


    static PyEdge findEdge( const GraphType & self ,const PyNode & u , const PyNode & v){
        return PyEdge(self,self.findEdge(u.item_,v.item_));
    }
    static index_type nodeId( const GraphType & self,const PyNode & node ){return  self.id(node.item_);}
    static index_type edgeId( const GraphType & self,const PyEdge & edge ){return  self.id(edge.item_);}
    static index_type arcId(  const GraphType & self,const PyNode & node ){return  self.id(node.item_);}
    //static EdgeIteratorHolderType edgeHolder(const GraphType & self){return EdgeIteratorHolderType(self);}
    //static NodeIteratorHolderType nodeHolder(const GraphType & self){return NodeIteratorHolderType(self);}

};

template<class GRAPH>
class MergeGraphAddOnVisitor 
:   public boost::python::def_visitor<MergeGraphAddOnVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH GraphType;
    // Lemon Graph Typedefs
    
    typedef typename GraphType::index_type       index_type;
    typedef typename GraphType::Edge             Edge;
    typedef typename GraphType::Node             Node;
    //typedef typename GraphType::Arc              Arc;



    template <class classT>
    void visit(classT& c) const
    {
        c
            // basic properties
            .add_property("initNodeNum",    &GraphType::initNumberOfNodes,"get the initial number of nodes" )
            .add_property("initEdgeNum",    &GraphType::initNumberOfEdges,"get the initial number of edges" )

            // special iterators
            .def("nodeIdIterator", python::range(
                    &GraphType::nodeIdsBegin  , 
                    &GraphType::nodeIdsEnd
                )
            )
            .def("edgeIdIterator", python::range(
                    &GraphType::nodeIdsBegin  , 
                    &GraphType::nodeIdsEnd
                )
            )
        ;
    }

   //static void foo(X& self);
   // static void bar(X& self);
};


} // end namespace vigra

#endif // EXPORT_GRAPH_VISITOR_HXX