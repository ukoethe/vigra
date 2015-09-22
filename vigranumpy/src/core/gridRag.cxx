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
#include "export_graph_rag_visitor.hxx"
#include "export_graph_algorithm_visitor.hxx"
#include "export_graph_shortest_path_visitor.hxx"
#include "export_graph_hierarchical_clustering_visitor.hxx"

#include "vigra/numpy_array.hxx"
#include "vigra/numpy_array_converters.hxx"
#include "vigra/adjacency_list_graph.hxx"
#include "vigra/graph_features.hxx"



namespace python = boost::python;

namespace vigra{





    template<unsigned int DIM>
    struct ExportFeatureExtractor{


        typedef GridRagFeatureExtractor<DIM, UInt32> GridRagFeatureExtractorType;

        static void exportFeatureExtractor(
            const std::string & clsName
        ){

            python::class_<GridRagFeatureExtractorType>(
                clsName.c_str(),
                python::no_init
            )
            .def("accumulatedFeatures", registerConverters(&accumulatedFeatures<float>),
                (
                    python::arg("data"),
                    python::arg("minVal"),
                    python::arg("maxVal"),
                    python::arg("out") = python::object()
                )
            )
            .def("geometricFeatures", registerConverters(&geometricFeatures),
                (
                    python::arg("out") = python::object()
                )
            )

            ;

            python::def("gridRagFeatureExtractor", registerConverters(&factory),
                python::return_value_policy<python::manage_new_object>()
            );
        }

        // the factory
        static GridRagFeatureExtractorType * factory(
            const AdjacencyListGraph & graph,
            NumpyArray<DIM, UInt32> labels
        ){
            if(graph.edgeNum() != graph.maxEdgeId() + 1){
                throw std::runtime_error("graph.edgeNum() == graph.maxEdgeId() +1 is violated");
            }
            return new GridRagFeatureExtractorType(graph, labels);
        }


        template<class DATA_TYPE>
        static NumpyAnyArray accumulatedFeatures(
            const GridRagFeatureExtractorType & extractor,
            const NumpyArray<DIM, DATA_TYPE>  & data,
            const DATA_TYPE minVal,
            const DATA_TYPE maxVal,
            NumpyArray<2, float> out
        ){


            TinyVector<UInt32,2> outShape(extractor.edgeNum(),100);
            out.reshape(outShape);
            extractor.accumulatedFeatures(data, minVal, maxVal, out);
            return out;
        }

        static NumpyAnyArray geometricFeatures(
            const GridRagFeatureExtractorType & extractor,
            NumpyArray<2, float> out
        ){
            TinyVector<UInt32,2> outShape(extractor.edgeNum(),100);
            out.reshape(outShape);
            extractor.geometricFeatures(out);
            return out;
        }

    };




    template<unsigned int DIM>
    void defineGridRag(){

        std::string clsName;
        if (DIM == 2)
            clsName = "GridRagFeatureExtractor2D";
        if (DIM == 3)
            clsName = "GridRagFeatureExtractor3D";
        if (DIM == 4)
            clsName = "GridRagFeatureExtractor4D";

        ExportFeatureExtractor<DIM>::exportFeatureExtractor(clsName) ;
    }


    template void defineGridRag<2>();
    template void defineGridRag<3>();


} 


