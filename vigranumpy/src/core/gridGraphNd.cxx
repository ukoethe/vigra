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





    template<unsigned int DIM>
    NumpyAnyArray pySerializeAffiliatedEdges(
        const GridGraph<DIM, boost::undirected_tag> & gridGraph,
        const AdjacencyListGraph & rag,
        const typename AdjacencyListGraph:: template EdgeMap< std::vector<  typename GridGraph<DIM, boost::undirected_tag >::Edge   > > & affiliatedEdges,
        NumpyArray<1,UInt32> out
    ){
        // reshape
        const size_t sSize = affiliatedEdgesSerializationSize(gridGraph, rag, affiliatedEdges);
        out.reshapeIfEmpty( typename NumpyArray<1,UInt32>::difference_type(sSize) );

        // do serialization
        serializeAffiliatedEdges(gridGraph, rag, affiliatedEdges, out.begin());

        // return 
        return out;

    }


    template<unsigned int DIM>
    AdjacencyListGraph:: template EdgeMap< std::vector<  typename GridGraph<DIM, boost::undirected_tag >::Edge   > > * 
    pyDeserializeAffiliatedEdges(
        const GridGraph<DIM, boost::undirected_tag> & gridGraph,
        const AdjacencyListGraph & rag,
        NumpyArray<1,UInt32> serialization
    ){

        AdjacencyListGraph:: template EdgeMap< std::vector<  typename GridGraph<DIM, boost::undirected_tag >::Edge   > > * affEdges_
        = new AdjacencyListGraph:: template EdgeMap< std::vector<  typename GridGraph<DIM, boost::undirected_tag >::Edge   > > () ; 

        // do serialization
        deserializeAffiliatedEdges(gridGraph, rag, *affEdges_, serialization.begin(), serialization.end());

        // return 
        return affEdges_;
    }




    template<unsigned int DIM>
    void defineGridGraphRagSerialization(){
        python::def("_serialzieGridGraphAffiliatedEdges",registerConverters(&pySerializeAffiliatedEdges<DIM>),
            (
                python::arg("gridGraph"),
                python::arg("rag"),
                python::arg("affiliatedEdges"),
                python::arg("out") = python::object()
            )
        );

        python::def("_deserialzieGridGraphAffiliatedEdges",registerConverters(&pyDeserializeAffiliatedEdges<DIM>),
            (
                python::arg("gridGraph"),
                python::arg("rag"),
                python::arg("serialization")
            ),
            python::return_value_policy<python::manage_new_object>()
        );
    };

    template void defineGridGraphRagSerialization< 2 >();
    template void defineGridGraphRagSerialization< 3  >();

} 


