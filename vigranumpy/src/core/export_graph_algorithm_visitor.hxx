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
#include <vigra/metrics.hxx>
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
    const static unsigned int EdgeMapDim = IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension;
    const static unsigned int NodeMapDim = IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension;

    typedef NumpyArray<EdgeMapDim,   Singleband<float > > FloatEdgeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<float > > FloatNodeArray;
    typedef NumpyArray<NodeMapDim,   Singleband<UInt32> > UInt32NodeArray;
    typedef NumpyArray<NodeMapDim +1,Multiband <float > > MultiFloatNodeArray;

    typedef NumpyScalarEdgeMap<Graph,FloatEdgeArray>         FloatEdgeArrayMap;
    typedef NumpyScalarNodeMap<Graph,FloatNodeArray>         FloatNodeArrayMap;
    typedef NumpyScalarNodeMap<Graph,UInt32NodeArray>        UInt32NodeArrayMap;
    typedef NumpyMultibandNodeMap<Graph,MultiFloatNodeArray> MultiFloatNodeArrayMap;

    LemonGraphSegmentationAlgotithmVisitor(const std::string clsName)
    :clsName_(clsName){

    }
    std::string clsName_;
    template <class classT>
    void visit(classT& c) const
    {
        
        // free functions:
        // - node feature distance to edge weights
        // - watersheds-segmentation
        // - felzenwalb-segmentation
        // - labeling


        python::def("nodeFeatureDistToEdgeWeight",registerConverters(&pyNodeFeatureDistToEdgeWeight),
            (
                python::arg("graph"),
                python::arg("nodeFeatures"),
                python::arg("metric"),
                python::arg("out")=python::object()
            ),
            "convert node features to edge weights with the given metric"
        );

        python::def("watershedsSegmentation",registerConverters(&pyWatershedSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("seeds"),
                python::arg("out")=python::object()
            ),
            "Seeded watersheds on a edge weighted graph"
        );

        python::def("felzenszwalbSegmentation",registerConverters(&pyFelzenszwalbSegmentation),
            (
                python::arg("graph"),
                python::arg("edgeWeights"),
                python::arg("nodeSizes")=python::object(),
                python::arg("k")=300.0f,
                python::arg("out")=python::object()
            ),
            "Felzenwalb graph based segmentation"
        );

    }


   


    static NumpyAnyArray pyNodeFeatureDistToEdgeWeight(
        const GRAPH & g,
        const MultiFloatNodeArray & nodeFeaturesArray,
        const std::string & functor,
        FloatEdgeArray edgeWeightsArray
    ){
        edgeWeightsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(g) );

        if(functor=="eucledian" || functor=="norm" || functor=="l2"){
            typedef  metrics::Norm<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        if(functor=="squaredNorm"){
            typedef  metrics::SquaredNorm<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else if (functor=="manhattan" || functor=="l1"){
            typedef  metrics::Manhattan<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else if (functor=="chiSquared"){
            typedef  metrics::ChiSquared<float> DistFunctor;
            DistFunctor f;
            return pyNodeFeatureDistToEdgeWeightT<DistFunctor>(g,nodeFeaturesArray,f,edgeWeightsArray);
        }
        else{
            throw std::runtime_error(
                "distance not supported\n"
                "supported distance types:\n"
                "- eucledian/norm/l2\n"
                "- squaredNorm\n"
                "- manhattan/l1\n"
                "- chiSquared\n"
            );
        }
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




    static NumpyAnyArray pyWatershedSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        UInt32NodeArray seedsArray,
        UInt32NodeArray labelsArray
    ){
        // reize output ? 
        labelsArray.reshapeIfEmpty( IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // numpy arrays => lemon maps
        FloatEdgeArrayMap edgeWeightsArrayMap(g,edgeWeightsArray);
        UInt32NodeArrayMap seedsArrayMap(g,seedsArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        watershedsSegmentation(g,edgeWeightsArrayMap,seedsArrayMap,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

    static NumpyAnyArray pyFelzenszwalbSegmentation(
        const GRAPH & g,
        FloatEdgeArray edgeWeightsArray,
        FloatNodeArray nodeSizesArray,
        const float k,
        UInt32NodeArray labelsArray
    ){
        // reize output ? 
        labelsArray.reshapeIfEmpty(  IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(g) );

        // is size array empty?
        // if empty fill with ones
        if(nodeSizesArray.shape(0)==0){
            nodeSizesArray.reshapeIfEmpty(edgeWeightsArray.shape());
            std::fill(nodeSizesArray.begin(),nodeSizesArray.end(),1.0);
        }

        // numpy arrays => lemon maps
        FloatEdgeArrayMap  edgeWeightsArrayMap(g,edgeWeightsArray);
        FloatNodeArrayMap  nodeSizesArrayMap(g,nodeSizesArray);
        UInt32NodeArrayMap labelsArrayMap(g,labelsArray);

        // call algorithm itself
        felzenszwalbSegmentation(g,edgeWeightsArrayMap,nodeSizesArrayMap,k,labelsArrayMap);

        // retun labels
        return labelsArray;
    }

};




} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_ALGORITHM_VISITOR_HXX