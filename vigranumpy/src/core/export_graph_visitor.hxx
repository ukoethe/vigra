#ifndef EXPORT_GRAPH_VISITOR_HXX
#define EXPORT_GRAPH_VISITOR_HXX

#include <boost/python.hpp>

namespace python = boost::python;

namespace vigra{


template<class GRAPH,class NODE_RET_POLICY,class EDGE_RET_POLICY>
class LemonDirectedGraphCoreVisitor 
:   public boost::python::def_visitor<LemonDirectedGraphCoreVisitor<GRAPH,NODE_RET_POLICY,EDGE_RET_POLICY> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH GraphType;
    // Lemon Graph Typedefs
    
    typedef typename GraphType::index_type       index_type;
    typedef typename GraphType::Edge             Edge;
    typedef typename GraphType::Node             Node;
    //typedef typename GraphType::Arc              Arc;

    typedef typename GraphType::NodeIt              NodeIt;
    typedef typename GraphType::EdgeIt              EdgeIt;
    typedef typename GraphType::NeighborNodeIt      NeighborNodeIt;

    template <class classT>
    void visit(classT& c) const
    {
        c
            // basic properties
            .add_property("edgeNum",    &GraphType::edgeNum )
            .add_property("nodeNum",    &GraphType::nodeNum )
            .add_property("maxEdgeId",  &GraphType::maxEdgeId )
            .add_property("maxEdgeId",  &GraphType::maxEdgeId )
            // id functions
            .def("id",&nodeId)
            .def("id",&edgeId)
            // get edge / node from id
            .def("nodeFromId",&GraphType::nodeFromId,NODE_RET_POLICY() )
            .def("edgeFromId",&GraphType::edgeFromId,EDGE_RET_POLICY() )

            // basic iterators
            .def("nodeIterator", python::range< NODE_RET_POLICY >(
                    &GraphType::nodesBegin  , 
                    &GraphType::nodesEnd
                )
            )
            .def("edgeIterator", python::range< EDGE_RET_POLICY >(
                    &GraphType::edgesBegin  , 
                    &GraphType::edgesEnd
                )
            )
            .def("incEdgeIdIt",&incEdgeIdIt)
        ;
    }

   static index_type nodeId( const GraphType & self,const Node & node ){return  self.id(node);}
   static index_type edgeId( const GraphType & self,const Edge & edge ){return  self.id(edge);}
   static python::object  incEdgeIdIt( const GraphType & self,const Node & node ){
        return python::range< python::objects::default_iterator_call_policies, index_type >(
            boost::bind( boost::mem_fn(&GraphType::neigbourEdgeIdsBegin) , &self , node),
            boost::bind( boost::mem_fn(&GraphType::neigbourEdgeIdsEnd)   , &self , node)
        );
   }

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