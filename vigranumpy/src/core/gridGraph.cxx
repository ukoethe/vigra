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
#include <vigra/graph_algorithms.hxx>
#include <vigra/python_graph_generalization.hxx>


namespace python = boost::python;

namespace vigra{


    template<unsigned int DIM,class DTAG,class LABEL_TYPE>
    void pyGetRag(
        const GridGraph<DIM,DTAG>                                                                & graph,
        NumpyArray<DIM,LABEL_TYPE>                                                                labelMap,
        AdjacencyListGraph                                                                       & rag,
        AdjacencyListGraph:: template EdgeMap< std::vector<typename GridGraph<DIM,DTAG>::Edge> > & hyperEdges,
        const Int64                                                                                ignoreLabel=-1
    ){
        makeRegionAdjacencyGraph(graph,labelMap,rag,hyperEdges,ignoreLabel);
    }



    template<unsigned int DIM>
    void defineGridGraphT(const std::string & clsName){

        typedef GridGraph<DIM,boost::undirected_tag> Graph;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::Node Node;

        typedef typename AdjacencyListGraph:: template EdgeMap< std::vector<Edge> > HyperEdgeMap;
        typedef typename AdjacencyListGraph:: template NodeMap< std::vector<Node> > HyperNodeMap;


        const std::string hyperEdgeMapNamClsName = clsName + std::string("HyperEdgeMap");
        python::class_<HyperEdgeMap>(hyperEdgeMapNamClsName.c_str(),python::init<>())
        ;



        
        python::class_<Graph>(clsName.c_str(),python::init< TinyVector<Int64,DIM> >())
        .def(LemonDirectedGraphCoreVisitor<Graph>(clsName))
        .def("getRegionAdjacencyGraph",registerConverters(&pyGetRag<DIM,boost::undirected_tag,UInt32> ),
           (
               python::arg("labels"),
               python::arg("rag"),
               python::arg("hyperEdges"),
               python::arg("ignoreLabel")=-1
           )
        )
        ;

        
    }

    void defineGridGraph(){
        defineGridGraphT<2>("GridGraphUndirected2d");
        defineGridGraphT<3>("GridGraphUndirected3d");
    }

    void defineInvalid(){
        python::class_<lemon::Invalid>("Invalid",python::init<>())
        ;
    }

} 


