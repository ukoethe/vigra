/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpygraphs_PyArray_API
#define NO_IMPORT_ARRAY

#include "export_graph_visitor.hxx"

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/merge_graph/merge_graph_adaptor.hxx>
#include <vigra/merge_graph/maps/clustering_operator.hxx>
#include <vigra/python_graph_generalization.hxx>
#include <vigra/hierarchical_clustering.hxx>

namespace python = boost::python;

namespace vigra{

    template<class GRAPH>
    cluster_operators::EdgeWeightNodeFeatures<
        MergeGraphAdaptor<GRAPH>,
        NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > ,
        NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > , 
        NumpyMultibandNodeMap<GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> > >,
        NumpyScalarNodeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float> > ,
        NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > 
    > * 
    pyEdgeWeightNodeFeaturesConstructor(
        MergeGraphAdaptor<GRAPH> & mergeGraph,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeIndicatorMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeSizeMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> >  nodeFeatureMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float >               nodeSizeMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeMinWeightMapArray
    ){
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>               EdgeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float>               NodeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> > NodeMultibandFloatArray;
        typedef NumpyScalarEdgeMap<GRAPH,EdgeFloatArray>                EdgeFloatMap;
        typedef NumpyScalarNodeMap<GRAPH,NodeFloatArray>                NodeFloatMap;
        typedef NumpyMultibandNodeMap<GRAPH,NodeMultibandFloatArray>    NodeMultibandFloatMap;

        typedef cluster_operators::EdgeWeightNodeFeatures<
            MergeGraphAdaptor<GRAPH> ,
            NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > ,
            NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > , 
            NumpyMultibandNodeMap<GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> > >,
            NumpyScalarNodeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float> > ,
            NumpyScalarEdgeMap<   GRAPH, NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float> > 
        > OperatorType;

        // allocate the thin wrappers
        EdgeFloatMap           edgeIndicatorMap(mergeGraph.graph(),edgeIndicatorMapArray);
        EdgeFloatMap           edgeSizeMap(mergeGraph.graph(),edgeSizeMapArray);   
        NodeMultibandFloatMap  nodeFeatureMap(mergeGraph.graph(),nodeFeatureMapArray);
        NodeFloatMap           nodeSizeMap(mergeGraph.graph(),nodeSizeMapArray);
        EdgeFloatMap           edgeMinWeightMap(mergeGraph.graph(),edgeMinWeightMapArray);



        return new OperatorType(mergeGraph,
            edgeIndicatorMap,edgeSizeMap,
            nodeFeatureMap, nodeSizeMap,
            edgeMinWeightMap
        );
    }


    template<class CLUSTER_OPERATOR>
    void defineHierarchicalClustering(const std::string & opClsName){
        typedef CLUSTER_OPERATOR ClusterOperator;
        typedef HierarchicalClustering<ClusterOperator> HCluster;

        const std::string clsName = std::string("HierarchicalClustering")+ opClsName;
        python::class_<HCluster,boost::noncopyable>(
            clsName.c_str(),python::init<ClusterOperator &>()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()]
        );


    }


    template<class GRAPH>
    void defineMergeGraphT(const std::string & clsName){

        typedef GRAPH Graph;
        typedef MergeGraphAdaptor<Graph> MergeGraphAdaptor;
        typedef typename MergeGraphAdaptor::Edge Edge;
        typedef typename MergeGraphAdaptor::Node Node;


        const std::string mgAdaptorClsName = clsName + std::string("MergeGraphAdaptor");
        python::class_<MergeGraphAdaptor,boost::noncopyable>(
            mgAdaptorClsName.c_str(),python::init<const Graph &>()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()]
        )
        .def(LemonDirectedGraphCoreVisitor<MergeGraphAdaptor>(mgAdaptorClsName))
        ;

        // define operator and the cluster class for this operator
        {
            typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                EdgeFloatArray;
            typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float>                NodeFloatArray;
            typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> >  NodeMultibandFloatArray;
            typedef NumpyScalarEdgeMap<GRAPH,EdgeFloatArray>                EdgeFloatMap;
            typedef NumpyScalarNodeMap<GRAPH,NodeFloatArray>                NodeFloatMap;
            typedef NumpyMultibandNodeMap<GRAPH,NodeMultibandFloatArray>    NodeMultibandFloatMap;
            typedef cluster_operators::EdgeWeightNodeFeatures<MergeGraphAdaptor,EdgeFloatMap,EdgeFloatMap,
                NodeMultibandFloatMap,NodeFloatMap,EdgeFloatMap> OperatorType;
            const std::string operatorName = mgAdaptorClsName + std::string("EdgeWeightNodeFeatureOperator");

            python::class_<OperatorType  >(operatorName.c_str(),python::no_init)
            .def("__init__", python::make_constructor(&pyEdgeWeightNodeFeaturesConstructor<Graph>))//,
                //python::return_value_policy<python::manage_new_object>() )
            ;
            defineHierarchicalClustering<OperatorType>(operatorName);
        }


        
    }

    void defineMergeGraphNew(){
        defineMergeGraphT<AdjacencyListGraph>("AdjacencyListGraph");
    }

} 


