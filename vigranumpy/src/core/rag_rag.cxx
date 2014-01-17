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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyrag_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/rag/rag.hxx>
#include <vigra/hiracical_clustering.hxx>
#include "export_graph_visitor.hxx"

namespace python = boost::python;

namespace vigra{




//// helper functions
template<class RAG,class LABELS>
RAG * makeRag(const LABELS & labeling){
    return new RAG(labeling);
}


// exporter  / definer functions 

template<unsigned int DIM>
void defineRagClass(const std::string & clsName ){

    typedef Rag<DIM,UInt32> RagType;
    typedef typename RagType::InputLabelingView InputLabelingView;
    typedef typename RagType::InLabelType InLabelType;
    typedef NumpyArray<RagType::Dimension ,Singleband < InLabelType > > NumpyLabelArray;

    // rag class itself
    python::class_<RagType>(clsName.c_str(),python::init<>())
    //python::class_<RagType>("Rag",python::init<const NumpyLabelArray & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()] )
    .def(LemonDirectedGraphCoreVisitor<
            RagType,
            python::return_value_policy<python::return_by_value>,
            python::return_value_policy<python::return_by_value>
        >()
    )
    .def("__init__",python::make_constructor(registerConverters(&makeRag<RagType,NumpyLabelArray>) /*,python::with_custodian_and_ward<1,2>()*/ ) )
    .add_property("shape",&RagType::shape)
    .add_property("topolocgialShape",&RagType::topolocgialShape)
    ;
}

template<unsigned int DIM>
void defineRagMaps(const std::string & clsNamePrefix){
    typedef Rag<DIM,UInt32> RagType;
    typedef typename GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
    typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;

    
    const std::string edgeCoordClsName = clsNamePrefix+std::string("EdgeCoordinatesMap");
    python::class_<RagEdgeCoordinatesMap>(edgeCoordClsName.c_str(),python::init<const RagType & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const RagType & */>()] )
    .def("size",&RagEdgeCoordinatesMap::size,"size");
    ;

    const std::string edgeMapFloatClsName = clsNamePrefix+std::string("EdgeFloatMap");
    python::class_<RagEdgeFloatMap>(edgeMapFloatClsName.c_str(),python::init<const RagType & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const RagType & */>()] )
    .def("size",&RagEdgeFloatMap::size,"size");
    ;
}

// .def("__init__",make_constructor(registerConverters(&constructPythonMap<PythonMapType,MergeGraphType>)))

template<unsigned int DIM>
void defineRagFunctions(){

    typedef Rag<DIM,UInt32> RagType;
    typedef typename GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
    typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;
    typedef NumpyArray<RagType::Dimension ,vigra::Singleband < float > > SingleBandFloatImage;
    //typedef NumpyArray<RagType::Dimension ,float > SingleBandFloatImage;
    

    // free functions
    python::def("extractEdgeCoordinates",& extractEdgeCoordinates<DIM,UInt32>,
        "extract the coordinates of the given graph");


    python::def("extractEdgeFeaturesFromImage",
        registerConverters(
            & extractEdgeFeaturesFromImage<
                RagType,
                RagEdgeCoordinatesMap,
                SingleBandFloatImage,
                RagEdgeFloatMap
            >
        )
    );

    python::def("extractEdgeSizeFromCoords",
        & edgeSizeMapFromCoordMap<
            RagType,
            RagEdgeCoordinatesMap,
            RagEdgeFloatMap
        >
    );

    python::def("edgeValueImage",
        registerConverters(
            & edgeValueImage<
                RagType,
                RagEdgeCoordinatesMap,
                RagEdgeFloatMap,
                SingleBandFloatImage
            >
        )
    );
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


template<unsigned int DIM>
void defineHiracicalClustering(const std::string & clsNamePrefix){
    typedef Rag<DIM,UInt32> RagType;
    typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap; 
    typedef HiracicalClustering<RagType,RagEdgeFloatMap,RagEdgeFloatMap> HiracicalClusteringType;


    const std::string hcClsName = clsNamePrefix+std::string("HiracicalClustering");
    python::class_<HiracicalClusteringType,boost::noncopyable>(hcClsName.c_str(),python::init<const RagType & , RagEdgeFloatMap &, RagEdgeFloatMap &>() )
    .def("cluster",&HiracicalClusteringType::cluster)
    .def("transformInputMaps",&HiracicalClusteringType::transformInputMaps)
    .def("mergeTreeEncoding",&mergeTreeEncodingAsNumpyArray<HiracicalClusteringType>)

    .def("leafNodeIds",&leadNodeIdsAsNumpyArray<HiracicalClusteringType>,
        (
            python::arg("treeNodeId"),
            python::arg("out")=python::object()
        )
    )
    ;
}


void defineRag(){

    {
        typedef Rag<2,UInt32> RagType;
        typedef GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
        typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;
        typedef NumpyArray<RagType::Dimension ,vigra::Singleband < float > > SingleBandFloatImage;
        //typedef NumpyArray<RagType::Dimension ,float > SingleBandFloatImage;
        // the rag class itself
        defineRagClass<2>("Rag2d");
        // define all maps
        defineRagMaps<2>("Rag2d");
        // define free functions
        defineRagFunctions<2>();

        // define hiracical clustering
        defineHiracicalClustering<2>("Rag2d");
    }
    
    {
        typedef Rag<3,UInt32> RagType;
        typedef GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
        typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;
        typedef NumpyArray<RagType::Dimension ,vigra::Singleband < float > > SingleBandFloatImage;
        //typedef NumpyArray<RagType::Dimension ,float > SingleBandFloatImage;
        // the rag class itself
        defineRagClass<3>("Rag3d");
        // define all maps
        defineRagMaps<3>("Rag3d");
        // define free functions
        defineRagFunctions<3>();
    }
    
}

} // namespace vigra
