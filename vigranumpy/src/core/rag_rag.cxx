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
#include <vigra/ucm.hxx>
#include "export_graph_visitor.hxx"

namespace python = boost::python;

namespace vigra{




//// helper functions
template<class RAG,class LABELS>
RAG * makeRag(const LABELS & labeling){
    return new RAG(labeling);
}


// exporter  / definer functions 

void defineRagClass(){

    typedef Rag<2,UInt32> RagType;
    typedef RagType::InputLabelingView InputLabelingView;
    typedef RagType::InLabelType InLabelType;
    typedef NumpyArray<RagType::Dimension ,vigra::Singleband < InLabelType > > NumpyLabelArray;

    // rag class itself
    python::class_<RagType>("Rag",python::init<>())
    //python::class_<RagType>("Rag",python::init<const NumpyLabelArray & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()] )
    .def(LemonDirectedGraphCoreVisitor<
            RagType,
            python::return_value_policy<python::return_by_value>,
            python::return_value_policy<python::return_by_value>
        >()
    )
    .def("__init__",python::make_constructor(registerConverters(&makeRag<RagType,NumpyLabelArray>) /*,python::with_custodian_and_ward<1,2>()*/ ) )
    ;
}

void defineRagMaps(){
    typedef Rag<2,UInt32> RagType;
    typedef GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
    typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;

    
    python::class_<RagEdgeCoordinatesMap>("RagEdgeCoordinatesMap",python::init<const RagType & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const RagType & */>()] )
    .def("size",&RagEdgeCoordinatesMap::size,"size");
    ;

    
    python::class_<RagEdgeFloatMap>("RagEdgeFloatMap",python::init<const RagType & >()[python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const RagType & */>()] )
    .def("size",&RagEdgeFloatMap::size,"size");
    ;
}

// .def("__init__",make_constructor(registerConverters(&constructPythonMap<PythonMapType,MergeGraphType>)))


void defineRag(){

    typedef Rag<2,UInt32> RagType;
    typedef GraphCoordinateTraits<RagType>::EdgeCoordinatesMap RagEdgeCoordinatesMap;
    typedef DenseEdgeReferenceMap<RagType,float> RagEdgeFloatMap;
    typedef NumpyArray<RagType::Dimension ,vigra::Singleband < float > > SingleBandFloatImage;
    //typedef NumpyArray<RagType::Dimension ,float > SingleBandFloatImage;
    // the rag class itself
    defineRagClass();
    // define all maps
    defineRagMaps();

    

    // free functions
    python::def("extractEdgeCoordinates",& extractEdgeCoordinates<2,UInt32>,
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


    python::def("ucmTransform",
        & ucmTransform<
            RagType,
            RagEdgeFloatMap,
            RagEdgeFloatMap
        >
    );
}

} // namespace vigra
