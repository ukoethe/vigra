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
#include <vigra/graph_maps.hxx>
#include <vigra/python_graph.hxx>


namespace python = boost::python;

namespace vigra{




    template<class GRAPH, class T_NODE,class FUNCTOR , class OTF_EDGE_MAP>
    OTF_EDGE_MAP * makeImplicitEdgeMap(
        const GRAPH & graph,
        const typename PyNodeMapTraits<GRAPH, T_NODE>::Array & nodeArray
    ){
        // generate lemon compatible map to node array (cheap view)
        typename PyNodeMapTraits<GRAPH,   T_NODE>::Map nodeArrayMap(graph, nodeArray);
        FUNCTOR f;
        OTF_EDGE_MAP * res = new OTF_EDGE_MAP(graph, nodeArrayMap, f);
        return res;

    }



    template<class GRAPH, class T_NODE,class NODE_MAP, class FUNCTOR, class RESULT>
    void defineImplicitEdgeMapT(const std::string & clsName, const std::string & factoryName){

        
        typedef OnTheFlyEdgeMap2<GRAPH, NODE_MAP, FUNCTOR, RESULT> EdgeMap;

        python::class_<EdgeMap>(clsName.c_str(),python::no_init)
        ;



        python::def(factoryName.c_str(),registerConverters(&makeImplicitEdgeMap<GRAPH, T_NODE,FUNCTOR, EdgeMap>),
            python::with_custodian_and_ward_postcall< 0,1 ,
                python::with_custodian_and_ward_postcall< 0 ,2,
                    python::return_value_policy<   python::manage_new_object      
                >   >   >()  
        );
       
    }

    template<int DIM, class T_NODE, class T_RES, class FUNCTOR>
    void defineGridGraphImplicitEdgeMapT(const std::string & clsName, const std::string & factoryName){
            


        typedef GridGraph<DIM, boost::undirected_tag> Graph;
        typedef typename PyNodeMapTraits<Graph, T_NODE>::Map NodeMap;
       
        //typedef OnTheFlyEdgeMap2<Graph, NodeMap, FUNCTOR, T_RES> EdgeMap;
        defineImplicitEdgeMapT<Graph,T_NODE, NodeMap, FUNCTOR,  T_RES>(clsName,factoryName);





    }


    void defineGridGraphImplicitEdgeMap(){

        {
            typedef float NodeValue;
            typedef float EdgeValue;
            typedef MeanFunctor<EdgeValue> Functor;
            defineGridGraphImplicitEdgeMapT<3, NodeValue, EdgeValue,Functor>(
                std::string("ImplicitMEanEdgeMap_3d_float_float"),
                std::string("implicitMeanEdgeMap")
            );
        }
    }

} 


