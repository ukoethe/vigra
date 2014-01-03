#ifndef EXPORT_GRAPH_VISITOR_HXX
#define EXPORT_GRAPH_VISITOR_HXX

#include <boost/python.hpp>

namespace python = boost::python;

namespace vigra{


template<class GRAPH>
class LemonDirectedGraphCoreVisitor 
:   public boost::python::def_visitor<LemonDirectedGraphCoreVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH GraphType;
    // Lemon Graph Typedefs
    
    typedef typename GraphType::index_type       index_type;
    typedef typename GraphType::Edge             Edge;
    typedef typename GraphType::Node             Node;
    typedef typename GraphType::Arc              Arc;



    template <class classT>
    void visit(classT& c) const
    {
        c
            // basic properties
            .add_property("edgeNum",    &GraphType::edgeNum )
            .add_property("nodeNum",    &GraphType::nodeNum )
            .add_property("maxEdgeId",  &GraphType::maxEdgeId )
            .add_property("maxEdgeId",  &GraphType::maxEdgeId )
            // basic iterators
            .def("nodeIterator", python::range< python::return_internal_reference<> >(
                    &GraphType::nodesBegin  , 
                    &GraphType::nodesEnd
                )
            )
            .def("edgeIterator", python::range< python::return_value_policy<python::return_by_value> >(
                    &GraphType::edgesBegin  , 
                    &GraphType::edgesEnd
                )
            )
        ;
    }

   //static void foo(X& self);
   // static void bar(X& self);
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
    typedef typename GraphType::Arc              Arc;



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