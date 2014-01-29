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
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/hierarchical_clustering.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/python_graph_generalization.hxx>

namespace python = boost::python;

namespace vigra{

    template<class GRAPH,class HYPER_EDGE_MAP,class SIZE_ARRAY>
    void pyHyperEdgeSizes(
        const GRAPH &                                                           g,
        const HYPER_EDGE_MAP &                                                  hyperEdges,
        SIZE_ARRAY                                                              sizeArray
    ){
        // convert numpy array to an edge map (API::LEMON)
        NumpyScalarEdgeMap<GRAPH,SIZE_ARRAY> sizeMap(g,sizeArray);

        hyperEdgeSizes(g,hyperEdges,sizeMap);
    }

    template<class RAG,class GRAPH,class LABEL_ARRAY,class SIZE_ARRAY>
    void pyHyperNodeSizes(
        const RAG   &                                                           rag,
        const GRAPH &                                                           graph,
        LABEL_ARRAY                                                             nodeLabelArray, // wrt graph
        SIZE_ARRAY                                                              sizeArray     // wrt rag
    ){
        // convert numpy array to an edge map (API::LEMON)
        NumpyScalarNodeMap<GRAPH,LABEL_ARRAY> labelMap(graph,nodeLabelArray); // wrt graph (for example grid graph)
        NumpyScalarNodeMap<RAG,  SIZE_ARRAY>  sizeMap(rag,sizeArray);  // wrt rag (based on graph)

        hyperNodeSizes(rag,graph,labelMap,sizeMap);
    }

    template<class RAG,class HYPER_EDGE_COORD_MAP,class IMAGE>
    void pyHyperEdgeImageFeatures(
        const RAG & rag,
        const HYPER_EDGE_COORD_MAP & hyperEdgeCoordMap,
        IMAGE  image,
        NumpyArray<
            IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicEdgeMapDimension,
            typename IMAGE::value_type
        > outArray
    ){
        // convert numpy array to an edge map (API::LEMON)
        typedef NumpyArray<
            IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicEdgeMapDimension,
            typename IMAGE::value_type
        > ArrayView;
        NumpyScalarEdgeMap<RAG,ArrayView> outMap(rag,outArray);
        hyperEdgeImageFeatures(rag,hyperEdgeCoordMap,image,outMap);
    }


    template<class RAG,class GRAPH,class LABEL_ARRAY,class IMAGE_ARRAY,class FEATURE_ARRAY>
    void pyHyperNodeImageFeatures(
        const RAG           & rag,
        const GRAPH         & graph,
        const LABEL_ARRAY   & labelArray,
        IMAGE_ARRAY           imageArray,
        FEATURE_ARRAY         featureArray
    ){
        // convert numpy array to an edge map (API::LEMON)
        NumpyScalarNodeMap<GRAPH,   LABEL_ARRAY> labelMap(graph,labelArray);    // wrt graph (for example grid graph)
        NumpyMultibandNodeMap<GRAPH,IMAGE_ARRAY> imageMap(graph,imageArray);    // wrt graph (for example grid graph)
        NumpyMultibandNodeMap<RAG,FEATURE_ARRAY> featureMap(rag,featureArray);  // wrt rag (based on graph)

        hyperNodeImageFeatures(rag,graph,labelMap,imageMap,featureMap);
    }


    template<class HCLUSTER>
    python::tuple mergeTreeEncodingAsNumpyArray(const HCLUSTER & hcluster) {
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
    python::tuple leadNodeIdsAsNumpyArray(
        const HCLUSTER &            hcluster,
        const typename HCLUSTER::MergeGraphIndexType treeNodeId,
        NumpyArray<1,UInt32>  leafes  = NumpyArray<1,UInt32>()
    ) {
        typedef typename HCLUSTER::MergeTreeEncoding      MergeTreeEncoding;
        typedef typename HCLUSTER::MergeGraphIndexType    MergeGraphIndexType;
        typedef typename HCLUSTER::ValueType              ValueType;
        
        leafes.reshapeIfEmpty( typename NumpyArray<1,UInt32>::difference_type( hcluster.graph().nodeNum()) );

        if(leafes.shape(0)!=hcluster.graph().nodeNum()){
            throw std::runtime_error("out.shape(0) must be equal nodeNum");
        }

        // todo make leafes size check
        const size_t leafNum=hcluster.leafNodeIds(treeNodeId,leafes.begin());
        return python::make_tuple(leafes,leafNum);
    } 


    void defineHyperEdgeSizes(){
        typedef GridGraph<2,boost::undirected_tag> GridGraphUndirected2d;
        typedef GridGraph<3,boost::undirected_tag> GridGraphUndirected3d;
        // gridgraph 2d and adjacency list graph
        {
            typedef AdjacencyListGraph::EdgeMap< std::vector<GridGraphUndirected2d::Edge > > HyperEdgeMap;
            typedef NumpyArray<IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicEdgeMapDimension,float>  SizeArray;
            python::def(
                "_hyperEdgeSizes",
                registerConverters(&pyHyperEdgeSizes<AdjacencyListGraph,HyperEdgeMap,SizeArray>)
            );
        }
        // gridgraph 3d and adjacency list graph
        {
            typedef AdjacencyListGraph::EdgeMap< std::vector<GridGraphUndirected3d::Edge > > HyperEdgeMap;
            typedef NumpyArray<IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicEdgeMapDimension,float>  SizeArray;
            python::def(
                "_hyperEdgeSizes",
                registerConverters(&pyHyperEdgeSizes<AdjacencyListGraph,HyperEdgeMap,SizeArray>)
            );
        }
    }
    void defineHyperNodeSizes(){
        typedef GridGraph<2,boost::undirected_tag> GridGraphUndirected2d;
        typedef GridGraph<3,boost::undirected_tag> GridGraphUndirected3d;
        // gridgraph 2d and adjacency list graph
        {
            typedef NumpyArray<IntrinsicGraphShape<GridGraphUndirected2d>::IntrinsicNodeMapDimension,UInt32>  LabelArray;
            typedef NumpyArray<IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicNodeMapDimension,float>      SizeArray;
            python::def(
                "_hyperNodeSizes",
                registerConverters(&pyHyperNodeSizes<AdjacencyListGraph,GridGraphUndirected2d,LabelArray,SizeArray>)
            );
        }
        // gridgraph 3d and adjacency list graph
        {
            typedef NumpyArray<IntrinsicGraphShape<GridGraphUndirected3d>::IntrinsicNodeMapDimension,UInt32>  LabelArray;
            typedef NumpyArray<IntrinsicGraphShape<AdjacencyListGraph>::IntrinsicNodeMapDimension,float>      SizeArray;
            python::def(
                "_hyperNodeSizes",
                registerConverters(&pyHyperNodeSizes<AdjacencyListGraph,GridGraphUndirected3d,LabelArray,SizeArray>)
            );
        }
    }
    void defineHyperEdgeFeatures(){
        typedef GridGraph<2,boost::undirected_tag> GridGraphUndirected2d;
        typedef GridGraph<3,boost::undirected_tag> GridGraphUndirected3d;
        //gridgraph 2d and adjacency list graph
        {
            typedef AdjacencyListGraph::EdgeMap< std::vector<GridGraphUndirected2d::Edge > > HyperEdgeMap;
            typedef NumpyArray<2,float>                     FeatureImage;
            python::def(
                "_hyperEdgeImageFeatures",
                registerConverters(&pyHyperEdgeImageFeatures<AdjacencyListGraph,HyperEdgeMap,FeatureImage>)
            );
        }
        //gridgraph 3d and adjacency list graph
        {
            typedef AdjacencyListGraph::EdgeMap< std::vector<GridGraphUndirected3d::Edge > > HyperEdgeMap;
            typedef NumpyArray<3,float>                     FeatureImage;
            python::def(
                "_hyperEdgeImageFeatures",
                registerConverters(&pyHyperEdgeImageFeatures<AdjacencyListGraph,HyperEdgeMap,FeatureImage>)
            );
        }
    }

    void defineHyperNodeFeatures(){
        typedef GridGraph<2,boost::undirected_tag> GridGraphUndirected2d;
        typedef GridGraph<3,boost::undirected_tag> GridGraphUndirected3d;
        //gridgraph 2d and adjacency list graph
        {
            typedef NumpyArray< 3,Multiband<float> >  ImageArray;
            typedef NumpyArray< 2,UInt32 >            LabelArray;
            typedef NumpyArray< 2,Multiband<float> >  FeatureArray;

            python::def("_hyperNodeImageFeatures",
                registerConverters(
                    &pyHyperNodeImageFeatures<AdjacencyListGraph,GridGraphUndirected2d,LabelArray,ImageArray,FeatureArray>
                )
            );
        }

        //gridgraph 3d and adjacency list graph
        {
            typedef NumpyArray< 4,Multiband<float> >  ImageArray;
            typedef NumpyArray< 3,UInt32 >            LabelArray;
            typedef NumpyArray< 2,Multiband<float> >  FeatureArray;

            python::def("_hyperNodeImageFeatures",
                registerConverters(
                    &pyHyperNodeImageFeatures<AdjacencyListGraph,GridGraphUndirected3d,LabelArray,ImageArray,FeatureArray>
                )
            );
        }
    }

    template<class HCLUSTER>
    void pyReprNodeIds(
        const HCLUSTER &     hcluster,
        NumpyArray<1,UInt32> labels
    ){
        for(size_t i=0;i<labels.shape(0);++i)
            labels(i)=hcluster.reprNodeId(labels(i));
    }

    template<class GRAPH,class HCLUSTER>
    HCLUSTER * pyHierarchicalClusteringConstructor(
        const GRAPH & graph,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeIndicatorMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeSizeMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> >    nodeFeatureMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float >               nodeSizeMapArray,
        NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>                edgeMinWeightMapArray,
        typename HCLUSTER::Parameter & param
    ){
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>             EdgeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float>            NodeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> > NodeMultibandFloatArray;
        typedef NumpyScalarEdgeMap<GRAPH,EdgeFloatArray>          EdgeFloatMap;
        typedef NumpyScalarNodeMap<GRAPH,NodeFloatArray>          NodeFloatMap;
        typedef NumpyMultibandNodeMap<GRAPH,NodeMultibandFloatArray> NodeMultibandFloatMap;
        // allocate the thin wrappers
        EdgeFloatMap           edgeIndicatorMap(graph,edgeIndicatorMapArray);
        EdgeFloatMap           edgeSizeMap(graph,edgeSizeMapArray);   
        NodeMultibandFloatMap  nodeFeatureMap(graph,nodeFeatureMapArray);
        NodeFloatMap           nodeSizeMap(graph,nodeSizeMapArray);
        EdgeFloatMap           edgeMinWeightMap(graph,edgeMinWeightMapArray);
        return new HCLUSTER(
            graph,
            edgeIndicatorMap,
            edgeSizeMap,
            nodeFeatureMap,
            nodeSizeMap,
            edgeMinWeightMap,
            param
        );
    }


    /*
    template<class GRAPH>
    void defineHierarchicalClustering(const std::string clsNamePostfix){
        typedef GRAPH GraphType;
        // python edge map type
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicEdgeMapDimension,float>             EdgeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension,float>             NodeFloatArray;
        typedef NumpyArray< IntrinsicGraphShape<GRAPH>::IntrinsicNodeMapDimension+1,Multiband<float> > NodeMultibandFloatArray;
        typedef NumpyScalarEdgeMap<GRAPH,EdgeFloatArray>          EdgeFloatMap;
        typedef NumpyScalarNodeMap<GRAPH,NodeFloatArray>          NodeFloatMap;
        typedef NumpyMultibandNodeMap<GRAPH,NodeMultibandFloatArray> NodeMultibandFloatMap;
        typedef HierarchicalClustering<
            GraphType,
            EdgeFloatMap,
            EdgeFloatMap,
            NodeMultibandFloatMap,
            NodeFloatMap,
            EdgeFloatMap
        > HierarchicalClusteringType;

        // parameter object
        typedef typename HierarchicalClusteringType::Parameter Parameter;
        const std::string hcParamName = std::string("HierarchicalClusteringParameter")+clsNamePostfix;
        python::class_<Parameter>(hcParamName.c_str(),
            python::init<>()
        )
        .def_readwrite("nodeNumStopCond", &Parameter::nodeNumStopCond_)
        .def_readwrite("buildMergeTree", &Parameter::buildMergeTreeEncoding_)
        .def_readwrite("verbose", &Parameter::verbose_)
        .def_readwrite("beta", &Parameter::beta_)
        .def_readwrite("degree1Fac", &Parameter::degree1Fac_)
        .def_readwrite("nodeDistType", &Parameter::nodeDistType_)
        .def_readwrite("wardness", &Parameter::wardness_)
        ;

        // factory function
        const std::string hcFactoryName = std::string("_hierarchicalClustering")+clsNamePostfix;
        python::def(hcFactoryName.c_str(),
            registerConverters(pyHierarchicalClusteringConstructor<GraphType,HierarchicalClusteringType>),
            python::return_value_policy<python::manage_new_object>() );

        const std::string hcClsName = std::string("HierarchicalClustering")+clsNamePostfix;
        python::class_<HierarchicalClusteringType,boost::noncopyable>(hcClsName.c_str(),python::no_init )
        .def("cluster",&HierarchicalClusteringType::cluster)
        .def("transformInputMaps",&HierarchicalClusteringType::transformInputMaps)
        .def("mergeTreeEncoding",&mergeTreeEncodingAsNumpyArray<HierarchicalClusteringType>)
        .def("mergeGraph",&HierarchicalClusteringType::mergeGraph,python::return_internal_reference<>() )
        .def("leafNodeIds",&leadNodeIdsAsNumpyArray<HierarchicalClusteringType>,
            (
                python::arg("treeNodeId"),
                python::arg("out")=python::object()
            )
        )
        .def("reprNodeIds",registerConverters(&pyReprNodeIds<HierarchicalClusteringType>))
        ;
    }
    */
    void defineGraphAlgorithms(){
        
        //typedef AdjacencyListGraph  AdjacencyListGraph;
        typedef GridGraph<2,boost::undirected_tag> GridGraphUndirected2d;
        typedef GridGraph<3,boost::undirected_tag> GridGraphUndirected3d;


        defineHyperEdgeSizes();
        defineHyperNodeSizes();
        defineHyperEdgeFeatures();
        defineHyperNodeFeatures();
        //defineHierarchicalClustering<AdjacencyListGraph>("AdjacencyListGraph");

    }

} 


