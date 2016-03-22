#ifndef VIGRA_EXPORT_GRID_RAG_2D_VISITOR_HXX
#define VIGRA_EXPORT_GRID_RAG_2D_VISITOR_HXX
//#define NO_IMPORT_ARRAY

/*boost python before anything else*/
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

/*std*/
#include <sstream>
#include <string>

/*vigra*/
#include "vigra/numpy_array.hxx"
#include "vigra/numpy_array_converters.hxx"
#include "vigra/graphs.hxx"
#include "vigra/graph_maps.hxx"
#include "vigra/python_graph.hxx"
#include "vigra/graph_algorithms.hxx"
#include "vigra/metrics.hxx"
#include "vigra/multi_gridgraph.hxx"
#include "vigra/error.hxx"
#include "vigra/graph_rag_project_back.hxx"
#include "vigra/threadpool.hxx"
#include "vigra/grid_rag_2d_features.hxx"
#include "vigra/accumulator.hxx"
#include "vigra/splineimageview.hxx"
namespace python = boost::python;

namespace vigra{



template<class GRAPH>
class LemonGraphGridRag2dVisitor 
:   public boost::python::def_visitor<LemonGraphGridRag2dVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;
    typedef AdjacencyListGraph RagGraph;

    typedef LemonGraphGridRag2dVisitor<GRAPH> VisitorType;
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
    typedef NumpyArray<RagEdgeMapDim +1,Multiband <float > > RagMultiFloatEdgeArray;

    typedef NumpyScalarEdgeMap<RagGraph,RagFloatEdgeArray>         RagFloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagFloatNodeArray>         RagFloatNodeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagUInt32NodeArray>        RagUInt32NodeArrayMap;
    typedef NumpyScalarNodeMap<RagGraph,RagInt32NodeArray>         RagInt32NodeArrayMap;
    typedef NumpyMultibandNodeMap<RagGraph,RagMultiFloatNodeArray> RagMultiFloatNodeArrayMap;



    typedef typename RagGraph:: template EdgeMap< std::vector<Edge> > RagAffiliatedEdges;
    typedef std::vector<Edge> EdgeVec;

    typedef typename GraphDescriptorToMultiArrayIndex<Graph>::IntrinsicNodeMapShape NodeCoordinate;
    typedef NumpyArray<1,NodeCoordinate>  NodeCoorinateArray;

    LemonGraphGridRag2dVisitor(const std::string clsName)
    :clsName_(clsName){

    }


    template <class classT>
    void visit(classT& /*c*/) const
    {   

        python::class_<SortedLines>("GridRag2dSortedLines",python::no_init)
        ;
        
        python::def("_gridRagSortedLines", &pySortedLinesFactory,
                python::with_custodian_and_ward_postcall< 0,1 ,
                    python::with_custodian_and_ward_postcall< 0 ,2,
                        python::return_value_policy<   python::manage_new_object
                > > >()
        );

        // typedef OnTheFlyEdgeMap2<
        //     Graph, typename PyNodeMapTraits<Graph,float>::Map,
        //     MeanFunctor<float>, float
        // > ImplicitEdgeMap;


        python::def("_gridRag2dCurvatureFeatures",
            registerConverters(
                &pyCurvatureFeatures
            ),
            (
                python::arg("rag"),
                python::arg("sortedLines"),
                python::arg("nThreads") = -1,
                python::arg("out")=python::object()
            )
        );

        python::def("_gridRag2dGeometricFeatures",
            registerConverters(
                &pyGeometricFeatures
            ),
            (
                python::arg("rag"),
                python::arg("sortedLines"),
                python::arg("nThreads") = -1,
                python::arg("out")=python::object()
            )
        );
    }


    static SortedLines * pySortedLinesFactory(
        const AdjacencyListGraph & g, 
        NumpyArray<2, UInt32> labels
    ){
        SortedLines * ret = NULL;
        {
            PyAllowThreads _pythread;
            ret = new SortedLines(g, labels);
        }
        return ret;
    }




    static NumpyAnyArray pyCurvatureFeatures(
        const RagGraph      & rag,
        const SortedLines   & sortedLines,
        const int           nThreads,
        NumpyArray<RagEdgeMapDim+1, float> ragEdgeFeaturesArray
    ){
        const size_t NFeatures = 14;
        // resize out
        typename MultiArray<RagEdgeMapDim+1,int>::difference_type outShape;
        for(size_t d=0;d<RagEdgeMapDim;++d){
            outShape[d]=IntrinsicGraphShape<RagGraph>::intrinsicEdgeMapShape(rag)[d];
        }
        outShape[RagEdgeMapDim]= NFeatures;
        ragEdgeFeaturesArray.reshapeIfEmpty(outShape);

        {
            PyAllowThreads _pythread;
            curvature(sortedLines,nThreads, ragEdgeFeaturesArray);
        }


        return ragEdgeFeaturesArray;
    }


    static NumpyAnyArray pyGeometricFeatures(
        const RagGraph      & rag,
        const SortedLines   & sortedLines,
        const int           nThreads,
        NumpyArray<RagEdgeMapDim+1, float> ragEdgeFeaturesArray
    ){
        const size_t NFeatures = 5;
        // resize out
        typename MultiArray<RagEdgeMapDim+1,int>::difference_type outShape;
        for(size_t d=0;d<RagEdgeMapDim;++d){
            outShape[d]=IntrinsicGraphShape<RagGraph>::intrinsicEdgeMapShape(rag)[d];
        }
        outShape[RagEdgeMapDim]= NFeatures;
        ragEdgeFeaturesArray.reshapeIfEmpty(outShape);

        {
            PyAllowThreads _pythread;
            geometricFeatures(sortedLines,nThreads, ragEdgeFeaturesArray);
        }


        return ragEdgeFeaturesArray;
    }


    

 

  
private:
    std::string clsName_;
};






} // end namespace vigra

#endif // VIGRA_EXPORT_GRID_RAG_2D_VISITOR_HXX
