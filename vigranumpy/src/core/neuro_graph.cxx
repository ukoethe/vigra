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
#include <vigra/neuro_features.hxx>


namespace python = boost::python;

namespace vigra{




    // python::with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const InputLabelingView & */>()

    template<class SELF>
    void pyAssignEdgeCues(
        SELF  & self,
        const NumpyArray<2, float> & edgeCues
    ){
        self.assignEdgeCues(edgeCues);
    }

    template<class SELF>
    void pyAssignNodeCues(
        SELF  & self,
        const NumpyArray<2, float> & nodeCues
    ){
        self.assignNodeCues(nodeCues);
    }


    template<class SELF>
    void pyAssignEdgeSizes(
        SELF  & self,
        const NumpyArray<1, float> & edgeSizes
    ){
        self.assignEdgeSizes(edgeSizes);
    }

    template<class SELF>
    void pyAssignNodeSizes(
        SELF  & self,
        const NumpyArray<1, float> & nodeSizes
    ){
        self.assignNodeSizes(nodeSizes);
    }


    void defineNeuroGraph(){

        typedef vigra::AdjacencyListGraph Graph;
        typedef NeuroDynamicFeatures PyNeuroDynamicFeatures;
        typedef PyNeuroDynamicFeatures::MergeGraph MergeGraph;

        python::class_<PyNeuroDynamicFeatures>(
            "NeuroDynamicFeatures", 
            python::init<const Graph &, MergeGraph &>()
            [
                python::with_custodian_and_ward<1 , 2,
                    python::with_custodian_and_ward<1 ,3 > 
                >()
            ]
        )
        .def("assignEdgeCues", pyAssignEdgeCues<PyNeuroDynamicFeatures> )
        .def("assignNodeCues", pyAssignNodeCues<PyNeuroDynamicFeatures> )
        .def("assignEdgeSizes", pyAssignEdgeSizes<PyNeuroDynamicFeatures> )
        .def("assignNodeSizes", pyAssignNodeSizes<PyNeuroDynamicFeatures> )

        .def("registerCallbacks", & PyNeuroDynamicFeatures::registerCallbacks)
        ;
        


    }


} 


