#ifndef VIGRA_EXPORT_GRAPH_HIERARCHICAL_CLUSTERING_VISTITOR_HXX
#define VIGRA_EXPORT_GRAPH_HIERARCHICAL_CLUSTERING_VISTITOR_HXX
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
#include <vigra/merge_graph_adaptor.hxx>
#include <vigra/hierarchical_clustering.hxx>
#include <vigra/timing.hxx>
namespace python = boost::python;

namespace vigra{



template<class GRAPH>
class LemonGraphHierachicalClusteringVisitor 
:   public boost::python::def_visitor<LemonGraphHierachicalClusteringVisitor<GRAPH> >
{
public:

    friend class def_visitor_access;

    typedef GRAPH Graph;
    typedef MergeGraphAdaptor<Graph> MergeGraph;

    typedef LemonGraphHierachicalClusteringVisitor<GRAPH> VisitorType;
    // Lemon Graph Typedefs
    
    typedef typename Graph::index_type       index_type;
    typedef typename Graph::Edge             Edge;
    typedef typename Graph::Node             Node;
    typedef typename Graph::Arc              Arc;
    typedef typename Graph::NodeIt           NodeIt;
    typedef typename Graph::EdgeIt           EdgeIt;
    typedef typename Graph::ArcIt            ArcIt;


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


    typedef cluster_operators::EdgeWeightNodeFeatures<
        MergeGraph,
        FloatEdgeArrayMap,
        FloatEdgeArrayMap,
        MultiFloatNodeArrayMap,
        FloatNodeArrayMap,
        FloatEdgeArrayMap,
        UInt32NodeArrayMap
    > DefaultClusterOperator;


    //typedef cluster_operators::EdgeWeightNodeFeatures2<
    //    MergeGraph,
    //    FloatEdgeArrayMap,
    //    FloatEdgeArrayMap,
    //    MultiFloatNodeArrayMap,
    //    FloatNodeArrayMap,
    //    FloatEdgeArrayMap
    //> NeuroClusterOperator;

    typedef cluster_operators::PythonOperator<MergeGraph> PythonClusterOperator;

   


    LemonGraphHierachicalClusteringVisitor(const std::string clsName)
    :clsName_(clsName){

    }



    void exportMergeGraph()const{
        const std::string mgAdaptorClsName = clsName_ + std::string("MergeGraph");
        python::class_<MergeGraph,boost::noncopyable>(
            mgAdaptorClsName.c_str(),python::init<const Graph &>()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()]
        )
        .def(LemonUndirectedGraphCoreVisitor<MergeGraph>(mgAdaptorClsName))
        .def("inactiveEdgesNode",&pyInactiveEdgesNode)
        .def("graph",&pyMergeGraphsGraph, python::return_internal_reference<>())
        .def("contractEdge",&pyContractEdgeA)
        .def("contractEdge",&pyContractEdgeB)
        .def("hasEdgeId",&pyHasEdgeId)

        .def("graphLabels",registerConverters(&pyCurrentLabeling<MergeGraph>),
            (
                python::arg("out")=python::object()
            )
        )

        ;

        python::def("__mergeGraph",&pyMergeGraphConstructor ,  
            python::with_custodian_and_ward_postcall< 0,1 ,
                    python::return_value_policy<   python::manage_new_object      >  >()  
        )
        ;
    }

    void exportHierarchicalClusteringOperators()const{
        {   
            const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("MinEdgeWeightNodeDistOperator");
            python::class_<DefaultClusterOperator  >(operatorName.c_str(),python::no_init)
            .def("__init__", python::make_constructor(&pyEdgeWeightNodeFeaturesConstructor))
            ;
            python::def("__minEdgeWeightNodeDistOperator",registerConverters(&pyEdgeWeightNodeFeaturesConstructor),
                python::with_custodian_and_ward_postcall< 0,1 ,
                    python::with_custodian_and_ward_postcall< 0 ,2,
                        python::with_custodian_and_ward_postcall< 0 ,3,
                            python::with_custodian_and_ward_postcall< 0 ,4,
                                python::with_custodian_and_ward_postcall< 0 ,5,
                                    python::with_custodian_and_ward_postcall< 0 ,6,
                                        python::with_custodian_and_ward_postcall< 0 ,7,
                                            python::return_value_policy<   python::manage_new_object      
                >  >    >   >   >   >   >   >()  
            );

        }
        //{   
        //    const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("NeuroOperator");
        //    python::class_<NeuroClusterOperator  >(operatorName.c_str(),python::no_init)
        //    .def("__init__", python::make_constructor(&pyNeuroConstructor))
        //    ;
        //    python::def("__neuroOperator",registerConverters(&pyNeuroConstructor),
        //        python::with_custodian_and_ward_postcall< 0,1 ,
        //            python::with_custodian_and_ward_postcall< 0 ,2,
        //                python::with_custodian_and_ward_postcall< 0 ,3,
        //                    python::with_custodian_and_ward_postcall< 0 ,4,
        //                        python::with_custodian_and_ward_postcall< 0 ,5,
        //                            python::with_custodian_and_ward_postcall< 0 ,6,
        //                                python::return_value_policy<   python::manage_new_object      
        //        >  >    >   >   >   >   >()  
        //    );
        //}
        {

            const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("PythonOperator");
            python::class_<PythonClusterOperator  >(operatorName.c_str(),python::no_init)
            .def("__init__", python::make_constructor(&pyPythonOperatorConstructor))//,
                //python::return_value_policy<python::manage_new_object>() )
            ;
            python::def("__pythonClusterOperator",registerConverters(&pyPythonOperatorConstructor),
                python::with_custodian_and_ward_postcall< 0,1 ,
                    python::with_custodian_and_ward_postcall< 0,2 ,
                        python::return_value_policy<   python::manage_new_object      >  >  >()  
            );
        }
    }

    template<class CLUSTER_OPERATOR>
    void exportHierarchicalClustering(const std::string & opClsName)const{
        typedef CLUSTER_OPERATOR ClusterOperator;
        typedef HierarchicalClustering<ClusterOperator> HCluster;

        const std::string clsName = std::string("HierarchicalClustering")+ opClsName;
        python::class_<HCluster,boost::noncopyable>(
            clsName.c_str(),python::init<ClusterOperator &>()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()]
        )
        .def("cluster",&HCluster::cluster)
        .def("reprNodeIds",registerConverters(&pyReprNodeIds<HCluster>))
        .def("resultLabels",registerConverters(&pyResultLabels<HCluster>),
            (
                python::arg("out")=python::object()
            )
        )
        ;

        // free function
        python::def("__hierarchicalClustering",registerConverters(&pyHierarchicalClusteringConstructor<ClusterOperator>),
            python::with_custodian_and_ward_postcall< 0,1 ,
                    python::return_value_policy<   python::manage_new_object      >  >()  
        );
    }

    static const Graph & pyMergeGraphsGraph(const MergeGraph & mg){
        return mg.graph();
    }


    static MergeGraph * pyMergeGraphConstructor(const GRAPH & graph){
        return new MergeGraphAdaptor<GRAPH>(graph);
    }


    template<class MG>
    static NumpyAnyArray pyCurrentLabeling(
        const MG  & mergeGraph,
        UInt32NodeArray   resultArray
    ){
        resultArray.reshapeIfEmpty(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(mergeGraph.graph()));

        UInt32NodeArrayMap resultArrayMap(mergeGraph.graph(),resultArray);
        //std::cout<<"find result labels\n";

        //USETICTOC;

        //TIC;
        for(NodeIt iter(mergeGraph.graph());iter!=lemon::INVALID;++iter ){
            resultArrayMap[*iter]=mergeGraph.reprNodeId(mergeGraph.graph().id(*iter));
        }
        //TOC;
        return resultArray;
    }


    template <class classT>
    void visit(classT& c) const
    {   
        // the merge graph itself and factory functions to get a merge graph
        exportMergeGraph();

        // export the clustering operators
        exportHierarchicalClusteringOperators();

        // export Hierarchical Clustering (for all  cluster operators)
        //{
        //    const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("NeuroOperator");
        //    exportHierarchicalClustering<NeuroClusterOperator>(operatorName);
        //}
        {
            const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("MinEdgeWeightNodeDistOperator");
            exportHierarchicalClustering<DefaultClusterOperator>(operatorName);
        }
        {
            const std::string operatorName = clsName_ + std::string("MergeGraph") + std::string("PythonOperator");
            exportHierarchicalClustering<PythonClusterOperator>(operatorName);
        }

    }



    static NodeHolder<MergeGraph> pyInactiveEdgesNode(
        const MergeGraph & mg,
        const EdgeHolder<MergeGraph> & edge
    ){
        return NodeHolder<MergeGraph>(mg,mg.inactiveEdgesNode(edge));
    }

    static EdgeHolder<MergeGraph> pyReprEdgeID(
        const MergeGraph & mg,
        const EdgeHolder<MergeGraph> & edge
    ){
        return EdgeHolder<MergeGraph>(mg,mg.reprEdge(edge));
    }
    

    static void pyContractEdgeA(
        MergeGraph & mg,
        const EdgeHolder<MergeGraph> & edge
    ){
        mg.contractEdge(edge);
    }

    static void pyContractEdgeB(
        MergeGraph & mg,
        const EdgeHolder<Graph> & graphEdge
    ){
        mg.contractEdge(mg.reprEdge(graphEdge));
    }

    static bool pyHasEdgeId(
        MergeGraph & mg,
        typename MergeGraph::index_type id
    ){
        return mg.hasEdgeId(id);
    }



    template<class CLUSTER_OP>
    static HierarchicalClustering<CLUSTER_OP> * pyHierarchicalClusteringConstructor(
        CLUSTER_OP & clusterOp,
        const size_t nodeNumStopCond,
        const bool buildMergeTreeEncoding


    ){
        typename HierarchicalClustering<CLUSTER_OP>::Parameter param;
        param.nodeNumStopCond_=nodeNumStopCond;
        param.buildMergeTreeEncoding_=buildMergeTreeEncoding;
        param.verbose_=true;
        return new  HierarchicalClustering<CLUSTER_OP>(clusterOp,param);
    }




    static DefaultClusterOperator * 
    pyEdgeWeightNodeFeaturesConstructor(
        MergeGraph &                mergeGraph,
        FloatEdgeArray              edgeIndicatorMapArray,
        FloatEdgeArray              edgeSizeMapArray,
        MultiFloatNodeArray         nodeFeatureMapArray,
        FloatNodeArray              nodeSizeMapArray,
        FloatEdgeArray              edgeMinWeightMapArray,
        UInt32NodeArray             nodeLabelArray,
        const float                 beta,
        const metrics::MetricType   nodeDistType,
        const float                 wardness,
        const float                 gamma
    ){

        FloatEdgeArrayMap       edgeIndicatorMap(mergeGraph.graph(),edgeIndicatorMapArray);
        FloatEdgeArrayMap       edgeSizeMap(mergeGraph.graph(),edgeSizeMapArray);   
        MultiFloatNodeArrayMap  nodeFeatureMap(mergeGraph.graph(),nodeFeatureMapArray);
        FloatNodeArrayMap       nodeSizeMap(mergeGraph.graph(),nodeSizeMapArray);
        FloatEdgeArrayMap       edgeMinWeightMap(mergeGraph.graph(),edgeMinWeightMapArray);
        UInt32NodeArrayMap      nodeLabelMap(mergeGraph.graph(),nodeLabelArray);


        return new DefaultClusterOperator(mergeGraph,
            edgeIndicatorMap,edgeSizeMap,
            nodeFeatureMap, nodeSizeMap,
            edgeMinWeightMap,nodeLabelMap,
            beta,nodeDistType,wardness, gamma
        );
    }






    static PythonClusterOperator * 
    pyPythonOperatorConstructor(
        MergeGraph & mergeGraph,
        python::object object,
        const bool useMergeNodeCallback,
        const bool useMergeEdgesCallback,
        const bool useEraseEdgeCallback
    ){
        return new PythonClusterOperator(mergeGraph,object,useMergeNodeCallback,useMergeEdgesCallback,useEraseEdgeCallback);
    }   



    template<class HCLUSTER>
    static void pyReprNodeIds(
        const HCLUSTER &     hcluster,
        NumpyArray<1,UInt32> labels
    ){
        for(MultiArrayIndex i=0; i<labels.shape(0); ++i)
            labels(i)=hcluster.reprNodeId(labels(i));
    }


    template<class HCLUSTER>
    static NumpyAnyArray pyResultLabels(
        const HCLUSTER  & hcluster,
        UInt32NodeArray   resultArray
    ){
        resultArray.reshapeIfEmpty(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(hcluster.graph()));

        UInt32NodeArrayMap resultArrayMap(hcluster.graph(),resultArray);
        //std::cout<<"find result labels\n";

        //USETICTOC;

        //TIC;
        for(NodeIt iter(hcluster.graph());iter!=lemon::INVALID;++iter ){
            resultArrayMap[*iter]=hcluster.mergeGraph().reprNodeId(hcluster.graph().id(*iter));
        }
        //TOC;
        return resultArray;
    }

    template<class HCLUSTER>
    static python::tuple mergeTreeEncodingAsNumpyArray(const HCLUSTER & hcluster) {
        typedef typename HCLUSTER::MergeTreeEncoding      MergeTreeEncoding;
        typedef typename HCLUSTER::MergeGraphIndexType    MergeGraphIndexType;
        typedef typename HCLUSTER::ValueType              ValueType;
        const MergeTreeEncoding & encoding = hcluster.mergeTreeEndcoding();
        const MergeGraphIndexType numMerges = encoding.size();
        //CPP BUG?!?
        NumpyArray<1,ValueType> w = NumpyArray<1,ValueType>(typename NumpyArray<1,ValueType>::difference_type(numMerges));
        NumpyArray<2,MergeGraphIndexType> indices = NumpyArray<2,MergeGraphIndexType>(typename NumpyArray<2,MergeGraphIndexType>::difference_type(numMerges,3));
        for(MergeGraphIndexType m=0;m<numMerges;++m){
            w(int(m))=encoding[m].w_;
            indices(m,0)=encoding[m].a_;
            indices(m,1)=encoding[m].b_;
            indices(m,2)=encoding[m].r_;
        }
        return python::make_tuple(indices,w);
    } 


    template<class HCLUSTER>
    static python::tuple leafNodeIdsAsNumpyArray(
        const HCLUSTER &            hcluster,
        const typename HCLUSTER::MergeGraphIndexType treeNodeId,
        NumpyArray<1,UInt32>  leafes  = (NumpyArray<1,UInt32>())
    ) {        
        leafes.reshapeIfEmpty( typename NumpyArray<1,UInt32>::difference_type( hcluster.graph().nodeNum()) );
        if(leafes.shape(0)!=hcluster.graph().nodeNum()){
            throw std::runtime_error("out.shape(0) must be equal nodeNum");
        }

        // todo make leafes size check
        const size_t leafNum=hcluster.leafNodeIds(treeNodeId,leafes.begin());
        return python::make_tuple(leafes,leafNum);
    } 

    /*

    template<class T>
    static NumpyAnyArray pyRagProjectNodeFeaturesToBaseGraph(
        const RagGraph &                                         rag,
        const Graph    &                                         graph,
        const typename PyNodeMapTraits<Graph,   UInt32>::Array & labelsWhichGeneratedRagArray,
        const typename PyNodeMapTraits<RagGraph,T     >::Array & ragNodeFeaturesArray,
        const Int32                                              ignoreLabel=-1,
        typename PyNodeMapTraits<Graph,T>::Array                 graphNodeFeaturesArray=typename PyNodeMapTraits<Graph,T>::Array() // out
    ){
        // reshape out  ( last argument (out) will be reshaped if empty, and #channels is taken from second argument)
        reshapeNodeMapIfEmpty(graph,ragNodeFeaturesArray,graphNodeFeaturesArray);
        // numpy arrays => lemon maps 
        typename PyNodeMapTraits<Graph,   UInt32>::Map labelsWhichGeneratedRagArrayMap(graph, labelsWhichGeneratedRagArray);
        typename PyNodeMapTraits<RagGraph,T     >::Map ragNodeFeaturesArrayMap(rag,ragNodeFeaturesArray);
        typename PyNodeMapTraits<Graph,   T     >::Map graphNodeFeaturesArrayMap(graph,graphNodeFeaturesArray);
        // run algorithm
        for(typename Graph::NodeIt iter(graph);iter!=lemon::INVALID;++iter){
            if(ignoreLabel==-1 || static_cast<Int32>(labelsWhichGeneratedRagArrayMap[*iter])!=ignoreLabel)
                graphNodeFeaturesArrayMap[*iter]=ragNodeFeaturesArrayMap[rag.nodeFromId(labelsWhichGeneratedRagArrayMap[*iter])];
            else{
                // ???
                // aks U. Koethe here
            }
        }
        return graphNodeFeaturesArray; // out
    }
    */
private:
    std::string clsName_;
};




} // end namespace vigra

#endif // VIGRA_EXPORT_GRAPH_HIERARCHICAL_CLUSTERING_VISTITOR_HXX
