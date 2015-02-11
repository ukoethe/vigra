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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyilastiktools_PyArray_API
//#define NO_IMPORT_ARRAY


/*vigra*/
#include <vigra/ilastiktools/carving.hxx>


/*vigra python */
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "export_graph_visitor.hxx"
#include "export_graph_rag_visitor.hxx"
#include "export_graph_algorithm_visitor.hxx"
#include "export_graph_shortest_path_visitor.hxx"
#include "export_graph_hierarchical_clustering_visitor.hxx"



namespace python = boost::python;

namespace vigra{



} // namespace vigra

using namespace vigra;
using namespace boost::python;


template<unsigned int DIM, class LABEL_TYPE>
void pyAssignLabels(
    GridRag<DIM, LABEL_TYPE> & graph,
    const NumpyArray<DIM, LABEL_TYPE> & labels
){
    graph.assignLabels(labels);
}

template<unsigned int DIM, class LABEL_TYPE>
void pyAssignLabelsFromSerialization(
    GridRag<DIM, LABEL_TYPE> & graph,
    const NumpyArray<DIM, LABEL_TYPE> & labels,
    NumpyArray<1, UInt32> serialization 
){
    graph.assignLabelsFromSerialization(labels, serialization);
}



template<unsigned int DIM, class LABEL_TYPE,
         class FEATURES_IN>
NumpyAnyArray pyAccumulateEdgeFeatures(
    GridRag<DIM, LABEL_TYPE> & graph,
    const NumpyArray<DIM, FEATURES_IN> & featuresIn,
    NumpyArray<1, typename NumericTraits<FEATURES_IN>::RealPromote > out 
){
    typedef TinyVector<MultiArrayIndex, 1>  Shape1;
    Shape1 shape(graph.edgeNum());
    out.reshapeIfEmpty(shape);
    graph.accumulateEdgeFeatures(featuresIn, out);
    return out;
}

template<unsigned int DIM, class LABEL_TYPE>
void pyPreprocessing(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<DIM, LABEL_TYPE> & labels,
    const NumpyArray<DIM, float> & weightArray      
){
    gridSegmentor.preprocessing(labels, weightArray);
}


template<unsigned int DIM, class LABEL_TYPE>
void pyPreprocessingFromSerialization(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<DIM, LABEL_TYPE> & labels,
    const NumpyArray<1, LABEL_TYPE> & serialization,
    const NumpyArray<1, float> & edgeWeights,
    const NumpyArray<1, UInt8> & nodeSeeds,
    const NumpyArray<1, UInt8> & resultSegmentation
){
    gridSegmentor.preprocessingFromSerialization(labels, serialization,
                                                 edgeWeights, nodeSeeds,
                                                 resultSegmentation);
}







template<unsigned int DIM, class LABEL_TYPE>
void pyAddLabels(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<DIM, UInt8> & brushStroke,
    const TinyVector<MultiArrayIndex, DIM> roiBegin,
    const TinyVector<MultiArrayIndex, DIM> roiEnd,
    const UInt8 maxValidLabel
){
    gridSegmentor.addLabels(brushStroke, roiBegin, roiEnd, maxValidLabel);;
}

template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyGetSegmentation(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const TinyVector<MultiArrayIndex, DIM> roiBegin,
    const TinyVector<MultiArrayIndex, DIM> roiEnd,
    NumpyArray<DIM, UInt8>  segmentation
){
    typedef TinyVector<MultiArrayIndex, DIM>  ShapeN;
    ShapeN shape(roiEnd-roiBegin);
    segmentation.reshapeIfEmpty(shape);

    {
        PyAllowThreads _pythread;
        gridSegmentor.getSegmentation(roiBegin, roiEnd, segmentation);
    }
   
    return segmentation;
}


template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyGetSuperVoxelSeg(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, UInt8>  segmentation
){
    typedef TinyVector<MultiArrayIndex, 1>  Shape1;
    Shape1 shape(gridSegmentor.maxNodeId()+1);
    segmentation.reshapeIfEmpty(shape);

    {
        PyAllowThreads _pythread;
        gridSegmentor.getSuperVoxelSeg(segmentation);
    }
   
    return segmentation;
}

template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyGetSuperVoxelSeeds(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, UInt8>  seeds
){
    typedef TinyVector<MultiArrayIndex, 1>  Shape1;
    Shape1 shape(gridSegmentor.maxNodeId()+1);
    seeds.reshapeIfEmpty(shape);

    {
        PyAllowThreads _pythread;
        gridSegmentor.getSuperVoxelSeeds(seeds);
    }
   
    return seeds;
}





template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pySerializeGraph(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, UInt32> serialization 
){
    serialization.reshapeIfEmpty( NumpyArray<1, UInt32>::difference_type(gridSegmentor.graph().serializationSize()));
    gridSegmentor.graph().serialize(serialization.begin());
    return serialization;
}

template<unsigned int DIM, class LABEL_TYPE>
void pyDeserializeGraph(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<1, UInt32> & serialization 
){
    gridSegmentor.graph().clear();
    gridSegmentor.graph().deserialize(serialization.begin(),serialization.end());
}




template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyEdgeWeights(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, float> out
){
    out.reshapeIfEmpty( 
        NumpyArray<1, UInt32>::difference_type(gridSegmentor.edgeNum())
    );
    out = gridSegmentor.edgeWeights();
    return out;
}

template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyNodeSeeds(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, UInt8> out
){
    out.reshapeIfEmpty( 
        NumpyArray<1, UInt32>::difference_type(gridSegmentor.maxNodeId()+1)
    );
    out = gridSegmentor.nodeSeeds();
    return out;
}

template<unsigned int DIM, class LABEL_TYPE>
NumpyAnyArray pyGetResultSegmentation(
    const GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    NumpyArray<1, UInt8> out
){
    out.reshapeIfEmpty( 
        NumpyArray<1, UInt32>::difference_type(gridSegmentor.maxNodeId()+1)
    );
    out = gridSegmentor.resultSegmentation();
    return out;
}



template<unsigned int DIM, class LABEL_TYPE>
void pySetSeeds(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<2, Int64> & fgSeeds,
    const NumpyArray<2, Int64> & bgSeeds
){
    gridSegmentor.setSeeds(fgSeeds, bgSeeds);
}

template<unsigned int DIM, class LABEL_TYPE>
void pySetResulFgObj(
    GridSegmentor<DIM , LABEL_TYPE, float> & gridSegmentor,
    const NumpyArray<1, Int64> & fgNodes
){
    gridSegmentor.setResulFgObj(fgNodes);
}

template<unsigned int DIM, class LABEL_TYPE>
void defineGridRag(const std::string & clsName){


    typedef GridRag<DIM, LABEL_TYPE> Graph;

    python::class_<Graph>(clsName.c_str(),python::init<  >())
        //.def("__init__",python::make_constructor(&pyGridGraphFactory3d<DIM,boost::undirected_tag>))
        .def(LemonUndirectedGraphCoreVisitor<Graph>(clsName))
        .def(LemonGraphAlgorithmVisitor<Graph>(clsName))
        // my functions
        .def("assignLabels",registerConverters(&pyAssignLabels<DIM, LABEL_TYPE>))
        .def("accumulateEdgeFeatures", 
            registerConverters(&pyAccumulateEdgeFeatures<DIM, LABEL_TYPE, float>),
            (
                python::arg("features"),
                python::arg("out") = python::object()
            )
        )

    ;
}




template<unsigned int DIM, class LABEL_TYPE>
void defineGridSegmentor(const std::string & clsName){


    typedef GridSegmentor<DIM, LABEL_TYPE, float> Segmentor;

    python::class_<Segmentor>(clsName.c_str(),python::init<  >())
        .def("preprocessing", 
            registerConverters( & pyPreprocessing<DIM, LABEL_TYPE>),
            (
                python::arg("labels"),
                python::arg("weightArray")
            )
        )
        .def("preprocessingFromSerialization", 
            registerConverters( & pyPreprocessingFromSerialization<DIM, LABEL_TYPE>),
            (
                python::arg("labels"),
                python::arg("serialization"),
                python::arg("edgeWeights"),
                python::arg("nodeSeeds"),
                python::arg("resultSegmentation")
            )
        )
        .def("addSeeds", 
            registerConverters( & pyAddLabels<DIM, LABEL_TYPE>),
            (
                python::arg("brushStroke"),
                python::arg("roiBegin"),
                python::arg("roiEnd"),
                python::arg("maxValidLabel")
            )
        )
        .def("setSeeds", 
            registerConverters( & pySetSeeds<DIM, LABEL_TYPE>),
            (
                python::arg("fgSeeds"),
                python::arg("bgSeeds")
            )
        )
        .def("setResulFgObj", 
            registerConverters( & pySetResulFgObj<DIM, LABEL_TYPE>),
            (
                python::arg("fgNodes")
            )
        )

        .def("getSegmentation", 
            registerConverters( & pyGetSegmentation<DIM, LABEL_TYPE>),
            (
                python::arg("roiBegin"),
                python::arg("roiEnd"),
                python::arg("out") = python::object()
            )
        )
        .def("nodeNum",&Segmentor::nodeNum)
        .def("edgeNum",&Segmentor::edgeNum)
        .def("maxNodeId",&Segmentor::maxNodeId)
        .def("maxEdgeId",&Segmentor::maxEdgeId)
        .def("run",&Segmentor::run)
        .def("clearSeeds",&Segmentor::clearSeeds)
        .def("serializeGraph", registerConverters(&pySerializeGraph<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
        .def("deserializeGraph", registerConverters(&pyDeserializeGraph<DIM, LABEL_TYPE>),
            (
                python::arg("serialization")
            )
        )
        .def("getEdgeWeights",registerConverters(pyEdgeWeights<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
        .def("getNodeSeeds",registerConverters(pyNodeSeeds<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
        .def("getResultSegmentation",registerConverters(pyGetResultSegmentation<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
        .def("getSuperVoxelSeg",registerConverters(pyGetSuperVoxelSeg<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
        .def("getSuperVoxelSeeds",registerConverters(pyGetSuperVoxelSeeds<DIM, LABEL_TYPE>),
            (
                python::arg("out") = python::object()
            )
        )
    ;
}





BOOST_PYTHON_MODULE_INIT(ilastiktools)
{
    import_vigranumpy();

    python::docstring_options doc_options(true, true, false);


    defineGridRag<3, vigra::UInt32>("GridRag_3D_UInt32");
    defineGridSegmentor<3, vigra::UInt32>("GridSegmentor_3D_UInt32");

}


