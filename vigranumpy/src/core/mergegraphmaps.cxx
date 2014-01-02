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
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include "vigra/merge_graph/merge_graph.hxx"
#include "vigra/merge_graph/maps/multi_array_map.hxx"
#include "vigra/merge_graph/maps/python_map.hxx"
namespace python = boost::python;

namespace vigra
{


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


template<class TO_CONSTRUCT,class MERGE_GRAPH>
TO_CONSTRUCT * constructPythonMap(
    const MERGE_GRAPH & graph,
    boost::python::object obj
){
    return new TO_CONSTRUCT(graph,obj);
}

template<class LABEL_TYPE,class VALUE_TYPE>
void defineMergeGraphMapsT(){
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


    /////////////////////////////////
    // python map
    /////////////////////////////////
    typedef PythonGraphMap<ValueType,MergeGraphType>  PythonMapType;

    class_<PythonMapType,boost::noncopyable>("PythonMap",no_init) 
    .def("__init__",make_constructor(registerConverters(&constructPythonMap<PythonMapType,MergeGraphType>)))
    .def("mergeNodeCallback",& PythonMapType:: template mergeNodeCallback<MergeItemsCallBackType>,"callback to merge node function")
    .def("mergeEdgeCallback",& PythonMapType:: template mergeEdgeCallback<MergeItemsCallBackType>,"callback to merge edges function")
    .def("eraseEdgeCallback",& PythonMapType:: template eraseEdgeCallback<EraseEdgeCallBackType>,"callback to erase edges function")
    ;

}


void defineMergeGraphMaps()
{
    using namespace python;
    docstring_options doc_options(true, true, false);
    NumpyArrayConverter<NumpyArray<2, double, UnstridedArrayTag> >();
    typedef UInt32 LabelType;
    typedef float  ValueType;    
    defineMergeGraphMapsT<LabelType,ValueType>();
}


} // namespace vigra