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

#define WITH_BOOST_GRAPH

/*vigra*/

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/python_graph.hxx>
#include <vigra/graph_rag.hxx>


namespace python = boost::python;

namespace vigra{


    void defineRagOptions();



    template<class GRAPH, class LABEL_TYPE>
    typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType * pyMakeRag(
        const GRAPH & graph,
        typename PyNodeMapTraits<GRAPH,LABEL_TYPE >::Array labelArray,
        AdjacencyListGraph & rag,
        const RagOptions  & options
    ){
        // numpy to lemon map
        typename PyNodeMapTraits<GRAPH ,LABEL_TYPE >::Map labelsMap(graph,labelArray);

        typedef  typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType AffiliatedEdgesType;

        AffiliatedEdgesType *  affiliatedEdges = new AffiliatedEdgesType();
        makeRag(graph, labelsMap, rag, *affiliatedEdges, options);

        return affiliatedEdges;

    }

    template<class GRAPH, class LABEL_TYPE>
    void defineMakeRagTemplate(){

        python::def("_makeRag", registerConverters(&pyMakeRag<GRAPH, LABEL_TYPE>),
            (
                (
                    python::arg("graph"),
                    python::arg("labels"),
                    python::arg("rag"),
                    python::arg("options")= RagOptions()
                )
            ),
            python::return_value_policy<python::manage_new_object>()
        )
        ;

    }

    void defineMakeRag(){

        typedef GridGraph<3,boost::undirected_tag> GridGraph3d;

        defineMakeRagTemplate<GridGraph3d, UInt32 >();

    }

    void defineNewRag(){
        defineRagOptions();
        defineMakeRag();
    }


} 


