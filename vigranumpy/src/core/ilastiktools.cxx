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


template<unsigned int DIM, class LABEL_TYPE
         class FEATURES_IN>
void pyAccumulateEdgeFeatures(
    GridRag<DIM, LABEL_TYPE> & graph,
    const NumpyArray<DIM, FEATURES_IN> & featuresIn,
    NumpyArray<DIM, typename NumericTraits<FEATURES_IN>::RealPromote > out 
){
    typedef TinyVector<1, DIM>  Shape1;
    Shape1 shape(graph.edgeNum());
    out.reshapeIfEmpty(shape);

    graph.accumulateEdgeFeatures(featuresIn, out);

    return out;
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


BOOST_PYTHON_MODULE_INIT(ilastiktools)
{
    import_vigranumpy();

    python::docstring_options doc_options(true, true, false);


    defineGridRag<3, vigra::UInt32>("GridRag_3D_UInt32");

}


