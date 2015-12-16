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
#include "export_graph_visitor.hxx"
#include "export_graph_rag_visitor.hxx"
#include "export_graph_algorithm_visitor.hxx"
#include "export_graph_shortest_path_visitor.hxx"
#include "export_graph_hierarchical_clustering_visitor.hxx"

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/python_graph.hxx>


namespace python = boost::python;

namespace vigra{






    template<unsigned int DIM,class DTAG>
    NodeHolder< GridGraph<DIM,boost::undirected_tag> > pyCoordinateToNode2d(
        const GridGraph<DIM,boost::undirected_tag> & g,
        const typename MultiArray<DIM,int>::difference_type & coordinate
    ){
        typename GridGraph<DIM,boost::undirected_tag>::Node node(coordinate);
        return NodeHolder<  GridGraph<DIM,boost::undirected_tag> >(g,node);
    }

    template<unsigned int DIM,class DTAG>
    GridGraph<DIM,DTAG>  * pyGridGraphFactory2d(
        typename MultiArray<DIM,int>::difference_type shape,
        const bool directNeighborhood
    ){
        return new GridGraph<DIM,DTAG>(shape,directNeighborhood?DirectNeighborhood:IndirectNeighborhood);
    }



    template<unsigned int DIM>
    void defineGridGraphRagSerialization();



    template<unsigned int DIM>
    void defineGridGraphT2d(const std::string & clsName){

        typedef GridGraph<DIM,boost::undirected_tag> Graph;
        typedef typename MultiArray<DIM,int>::difference_type ShapeType;

        
        python::class_<Graph>(clsName.c_str(),python::init< ShapeType >())
        .def("__init__",python::make_constructor(&pyGridGraphFactory2d<DIM,boost::undirected_tag>))
        .def(LemonUndirectedGraphCoreVisitor<Graph>(clsName))
        .def(LemonGraphAlgorithmVisitor<Graph>(clsName))
        .def(LemonGridGraphAlgorithmAddonVisitor<Graph>(clsName))
        .def(LemonGraphShortestPathVisitor<Graph>(clsName))
        .def(LemonGraphRagVisitor<Graph>(clsName))
        .def(LemonGraphHierachicalClusteringVisitor<Graph>(clsName))
        .def("coordinateToNode",pyCoordinateToNode2d<DIM,boost::undirected_tag>)
        ;

       
    }

    void defineGridGraph2d(){
        defineGridGraphT2d<2>("GridGraphUndirected2d");
        defineGridGraphRagSerialization<2>();
        //defineGridGraphT<3>("GridGraphUndirected3d");
    }


} 


