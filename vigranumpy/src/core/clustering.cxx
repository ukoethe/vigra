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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyclustering_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include "vigra/merge_graph/merge_graph.hxx"
#include "vigra/merge_graph/maps/multi_array_map.hxx"

namespace python = boost::python;

namespace vigra
{


template<class LABEL_TYPE>
void setInitalEdges(
    MergeGraph<LABEL_TYPE> & graph, 
    NumpyArray<2, LABEL_TYPE > edges
){
    const size_t numberOfEdges = edges.shape(0);
    //vigra_precondition(numberOfEdges==graph.initNumberOfEdges(),
    //           "wrong number of edges passed to setInitalEdge");

    for(size_t edgeIndex=0;edgeIndex<numberOfEdges;++edgeIndex){
        graph.setInitalEdge(edgeIndex,edges(edgeIndex,0),edges(edgeIndex,1));
    }
}

template<class LABEL_TYPE>
vigra::NumpyAnyArray stateOfInitalEdges(
    MergeGraph<LABEL_TYPE> & graph, 
    vigra::NumpyArray<1 , LABEL_TYPE> res = vigra::NumpyArray<1,LABEL_TYPE >()
){
    const size_t initNumberOfEdges=graph.initNumberOfEdges();
    res.reshapeIfEmpty(typename vigra::NumpyArray<1 , LABEL_TYPE>::difference_type(initNumberOfEdges));
    graph.stateOfInitalEdges(res.begin(),res.end());
    return res;
}




template<class MGRAPH>
void registerMergeNodeCallBack(
MGRAPH & graph, typename MGRAPH::MergeNodeCallBackType & callback
){
    graph.registerMergeNodeCallBack(callback);
}
template<class MGRAPH>
void registerMergeEdgeCallBack(
MGRAPH & graph, typename MGRAPH::MergeEdgeCallBackType & callback
){
    graph.registerMergeEdgeCallBack(callback);
}
template<class MGRAPH>
void registerEraseEdgeCallBack(
MGRAPH & graph, typename MGRAPH::EraseEdgeCallBackType & callback
){
    graph.registerEraseEdgeCallBack(callback);
}


template<class LABEL_TYPE>
void defineMergeGraph(){
    using namespace python;
    docstring_options doc_options(true, true, false);
    NumpyArrayConverter<NumpyArray<2, double, UnstridedArrayTag> >();
    typedef LABEL_TYPE LabelType;
    typedef MergeGraph<LabelType> MergeGraphType;
    typedef typename MergeGraphType::MergeItemsCallBackType MergeItemsCallBackType;
    typedef typename MergeGraphType::EraseEdgeCallBackType EraseEdgeCallBackType;
    ///////////////////////////////////////////////////////////////////////////////////////
    //  
    ///////////////////////////////////////////////////////////////////////////////////////
    class_<MergeGraphType,boost::noncopyable>("MergeGraph",init<const size_t,const size_t >
        (
            ( 
                arg("numberOfNodes"),
                arg("numberOfEdges")
            )
        ) 
    )
        .def("initNumberOfNodes",&MergeGraphType::initNumberOfNodes,"get the initial number of nodes")
        .def("initNumberOfEdges",&MergeGraphType::initNumberOfEdges,"get the initial number of edges")
        .def("numberOfNodes",&MergeGraphType::numberOfNodes,"get the current number of nodes")
        .def("numberOfEdges",&MergeGraphType::numberOfEdges,"get the current number of edges")
        .def("mergeParallelEdges",&MergeGraphType::mergeParallelEdges,"merge parallel / double edges")
        .def("mergeRegions",&MergeGraphType::mergeRegions,"mergeTwoRegions")
        .def("setInitalEdge",&MergeGraphType::setInitalEdge,
            (
                arg("initEdge"),
                arg("initNode0"),
                arg("initNode1")
            ),
            "set inital edges (bevore andy merging)"
        )
        .def("setInitalEdges",registerConverters(&setInitalEdges<LabelType>),
            (
                arg("edges")
            ),
            "a  initNumberOfEdges x 2 array"
        )

        .def("registerMergeNodeCallBack",&registerMergeNodeCallBack<MergeGraphType>,"register am MergeNode call-back")
        .def("registerMergeEdgeCallBack",&registerMergeEdgeCallBack<MergeGraphType>,"register am MergeEdge call-back")
        .def("registerEraseEdgeCallBack",&registerEraseEdgeCallBack<MergeGraphType>,"register am EraseNode call-back")

        .def("stateOfInitalEdges",vigra::registerConverters(&stateOfInitalEdges<LABEL_TYPE>),
            (
                arg("out")=object()
            ),
            "get edge state remapped on inital edges"
        )
        //.def("activeEdgeLabels",vigra::registerConverters(&activeEdgeLabels),
        //    (
        //        arg("out")=object()
        //    ),
        //    "activeEdgeLabels"
        //)
    ;
    class_<MergeItemsCallBackType >("MergeItemsCallBackType",init<>() ); 
    class_<EraseEdgeCallBackType >("EraseEdgeCallBackType",init<>() ); 
}

// helper function to construct some class from an array
template<class TO_CONSTRUCT,class ARRAY>
TO_CONSTRUCT * constructSumMap(
    ARRAY array
){
    return new TO_CONSTRUCT(array);
}

template<class TO_CONSTRUCT,class ARRAY,class WEIGHT_MAP>
TO_CONSTRUCT * constructWeightedMeanMap(
    ARRAY array,
    const WEIGHT_MAP & weightMap
){
    return new TO_CONSTRUCT(array,weightMap);
}


template<class TO_CONSTRUCT,class MERGE_GRAPH,class ARRAY,class EDGE_MAP,class NODE_MAP>
TO_CONSTRUCT * constructMinWeightEdgeMap(
    const MERGE_GRAPH & graph,
    ARRAY array,
    const EDGE_MAP & edgeMap,
    const NODE_MAP & nodeMap
){
    return new TO_CONSTRUCT(graph,array,edgeMap,nodeMap);
}


template<class LABEL_TYPE,class VALUE_TYPE>
void defineMaps(){
    using namespace python;
    docstring_options doc_options(true, true, false);
    NumpyArrayConverter<NumpyArray<2, double, UnstridedArrayTag> >();
    typedef LABEL_TYPE LabelType;
    typedef VALUE_TYPE ValueType;
    typedef MergeGraph<LabelType> MergeGraphType;

    typedef typename MergeGraphType::MergeItemsCallBackType MergeItemsCallBackType;
    typedef typename MergeGraphType::EraseEdgeCallBackType  EraseEdgeCallBackType;

    //////////////////////////////////////////////////////////////
    // SUM MAP 
    //////////////////////////////////////////////////////////////
    typedef view_maps::SumMap<0,ValueType>  SumMap0;
    typedef typename SumMap0::ArrayViewType SumMap0ArrayView;

    typedef view_maps::SumMap<1,ValueType>  SumMap1;
    typedef typename SumMap1::ArrayViewType SumMap1ArrayView;

    class_<SumMap0,boost::noncopyable>("SumMap0",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructSumMap<SumMap0,SumMap0ArrayView>)))
    .def("mergeCallback",& SumMap0:: template mergeCallback<MergeItemsCallBackType>,"callback to merge function", return_value_policy<return_by_value>())
    ;

    class_<SumMap1,boost::noncopyable>("SumMap1",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructSumMap<SumMap1,SumMap1ArrayView>)))
    .def("mergeCallback",& SumMap1:: template mergeCallback<MergeItemsCallBackType>,"callback to merge function", return_value_policy<return_by_value>())
    ;

    
    //////////////////////////////////////////////////////////////
    // WEIGHTED MEAN MAP 
    //////////////////////////////////////////////////////////////
    typedef view_maps::WeightedMeanMap<0,ValueType,SumMap0>     WeightedMeanMap0;
    typedef typename WeightedMeanMap0::ArrayViewType            WeightedMeanMap0ArrayView;

    typedef view_maps::WeightedMeanMap<1,ValueType,SumMap0>     WeightedMeanMap1;
    typedef typename WeightedMeanMap1::ArrayViewType            WeightedMeanMap1ArrayView;

    class_<WeightedMeanMap0,boost::noncopyable>("WeightedMeanMap0",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructWeightedMeanMap<WeightedMeanMap0,WeightedMeanMap0ArrayView,SumMap0>)))
    .def("mergeCallback",& WeightedMeanMap0:: template mergeCallback<MergeItemsCallBackType>,"callback to merge function")
    ;

    class_<WeightedMeanMap1,boost::noncopyable>("WeightedMeanMap1",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructWeightedMeanMap<WeightedMeanMap1,WeightedMeanMap1ArrayView,SumMap0>)))
    .def("mergeCallback",& WeightedMeanMap1:: template mergeCallback<MergeItemsCallBackType>,"callback to merge function")
    ;
    

    /////////////////////////////////
    // Min weight edge map
    /////////////////////////////////
    typedef view_maps::MinWeightEdgeMap<ValueType,MergeGraphType,WeightedMeanMap0,WeightedMeanMap1> MinWeightEdgeMap;
    typedef typename MinWeightEdgeMap::ArrayViewType             MinWeightEdgeMapArrayView;


    class_<MinWeightEdgeMap,boost::noncopyable>("MinWeightEdgeMap",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructMinWeightEdgeMap<MinWeightEdgeMap,MergeGraphType,MinWeightEdgeMapArrayView,WeightedMeanMap0,WeightedMeanMap1>)))
    .def("mergeEdgeCallback",& MinWeightEdgeMap:: template mergeEdgeCallback<MergeItemsCallBackType>,"callback to merge edges function")
    .def("eraseEdgeCallback",& MinWeightEdgeMap:: template eraseEdgeCallback<EraseEdgeCallBackType>,"callback to erase edges function")
    .def("minWeightEdgeLabel",&MinWeightEdgeMap::minWeightEdgeLabel,"get the label of the edge with min. weight")
    ;
}


void defineClustering()
{
    using namespace python;
    docstring_options doc_options(true, true, false);
    NumpyArrayConverter<NumpyArray<2, double, UnstridedArrayTag> >();


    typedef UInt32 LabelType;
    typedef float  ValueType;
    
    defineMergeGraph<LabelType>();
    defineMaps<LabelType,ValueType>();



}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(clustering)
{
    import_vigranumpy();
    defineClustering();
}
