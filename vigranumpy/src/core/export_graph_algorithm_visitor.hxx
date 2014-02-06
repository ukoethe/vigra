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
#include <vigra/python_graph_generalization.hxx>
#include <vigra/graph_algorithms.hxx>
namespace python = boost::python;

namespace vigra{


template<class GRAPH>
class LemonGraphSegmentationAlgotithmVisitor 
:   public boost::python::def_visitor<LemonGraphSegmentationAlgotithmVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;

    typedef LemonGraphSegmentationAlgotithmVisitor<GRAPH> VisitorType;
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

    typedef NumpyArray<IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension,float>  FloatEdgeArray;
    typedef NumpyArray<IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension,UInt32> NodeLabelArray;
    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>    FloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<Graph,NodeLabelArray>    NodeLabelArrayMap;

    LemonGraphSegmentationAlgotithmVisitor(const std::string clsName)
    :clsName_(clsName){

    }
    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {
        
        // free functions
        python::def("watershedsSegmentation",&pyWatershedSegmentation,
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );

    }



    static NumpyAnyArray pyWatershedSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        NodeLabelArray seedsArray,
        NodeLabelArray labelsArray
    ){
        // reize output ? 
        labelsArray.reshapeIfEmpty(seedsArray.shape());

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        NodeLabelArrayMap seedsArrayMap(g,seedsArray);
        NodeLabelArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        watershedsSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

};




} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX