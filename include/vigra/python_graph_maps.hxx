/************************************************************************/
/*                                                                      */
/*    Copyright 2011 by Ullrich Koethe and Thorsten Beier               */
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

#ifndef VIGRA_PYTHON_GRAPH_MAPS_HXX
#define VIGRA_PYTHON_GRAPH_MAPS_HXX

/*vigra*/
//#include <vigra/array_vector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/graphs.hxx>
#include <vigra/numpy_array.hxx>

namespace vigra{

template<class T>
struct IsMultiband : VigraFalseType{
};

template<class T>
struct IsMultiband<Multiband<T> > : VigraTrueType{
};




template<class G,class AV>
class NumpyScalarEdgeMap{

public:
    typedef G  Graph;
    typedef AV ArrayView;
    typedef typename  Graph::Edge                Key;
    typedef typename  ArrayView::value_type      Value;
    typedef typename  ArrayView::reference       Reference;
    typedef typename  ArrayView::const_reference ConstReference;

    typedef typename  Graph::Edge                key_type;
    typedef typename  ArrayView::value_type      value_type;
    typedef typename  ArrayView::reference       reference;
    typedef typename  ArrayView::const_reference const_reference;

    NumpyScalarEdgeMap()
    :   graph_(NULL),
        array_(){
    }

    NumpyScalarEdgeMap(const Graph & graph,ArrayView array)
    :   graph_(&graph),
        array_(array){
    }

    Reference operator[](const Key & key){
        return array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicEdgeCoordinate(*graph_,key)];
    }
    ConstReference operator[](const Key & key)const{
        return   array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicEdgeCoordinate(*graph_,key)];
    }
private:
    const Graph * graph_;
    MultiArrayView<IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension,Value> array_;

};

template<class G,class AV>
class NumpyScalarNodeMap{

public:
    typedef G  Graph;
    typedef AV ArrayView;
    typedef typename  Graph::Node                Key;
    typedef typename  ArrayView::value_type      Value;
    typedef typename  ArrayView::reference       Reference;
    typedef typename  ArrayView::const_reference ConstReference;

    typedef typename  Graph::Node                key_type;
    typedef typename  ArrayView::value_type      value_type;
    typedef typename  ArrayView::reference       reference;
    typedef typename  ArrayView::const_reference const_reference;
    //typedef Value &                                Reference;
    //typedef const Value &                          ConstReference;

    NumpyScalarNodeMap()
    :   graph_(NULL),
        array_(){
    }

    NumpyScalarNodeMap(const Graph & graph,ArrayView array)
    :   graph_(&graph),
        array_(array){
    }

    Reference operator[](const Key & key){
        return array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(*graph_,key)];
    }
    ConstReference operator[](const Key & key)const{
        return   array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(*graph_,key)];
    }
private:
    const Graph * graph_;
    MultiArrayView<IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension,Value> array_;

};


template<class G,class AV>
class NumpyMultibandNodeMap{

public:
    typedef G  Graph;
    typedef AV ArrayView;
    typedef typename  Graph::Node                Key;
    typedef typename  Graph::Node                key_type;

    //typedef typename  ArrayView::value_type      Value;
    //typedef typename  ArrayView::reference       Reference;
    //typedef typename  ArrayView::const_reference ConstReference;

    typedef  MultiArray<1,typename AV::value_type>           Value;
    typedef  MultiArrayView<1,typename AV::value_type>       Reference;
    typedef  MultiArrayView<1,typename AV::value_type> ConstReference;
    typedef  MultiArray<1,typename AV::value_type>           value_type;
    typedef  MultiArrayView<1,typename AV::value_type>       reference;
    typedef  MultiArrayView<1,typename AV::value_type> const_reference;
    //typedef Value &                                Reference;
    //typedef const Value &                          ConstReference;

    NumpyMultibandNodeMap()
    :   graph_(NULL),
        array_(){
    }

    NumpyMultibandNodeMap(const Graph & graph,ArrayView array)
    :   graph_(&graph),
        array_(array){
    }

    Reference operator[](const Key & key){
        return array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(*graph_,key)];
    }
    ConstReference operator[](const Key & key)const{
        return array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicNodeCoordinate(*graph_,key)];
    }
private:
    const Graph * graph_;
    mutable AV array_;

};


template<class G,class AV>
class NumpyMultibandEdgeMap{

public:
    typedef G  Graph;
    typedef AV ArrayView;
    typedef typename  Graph::Edge                Key;
    typedef typename  Graph::Edge                key_type;

    //typedef typename  ArrayView::value_type      Value;
    //typedef typename  ArrayView::reference       Reference;
    //typedef typename  ArrayView::const_reference ConstReference;

    typedef  MultiArray<1,typename AV::value_type>           Value;
    typedef  MultiArrayView<1,typename AV::value_type>       Reference;
    typedef  MultiArrayView<1,typename AV::value_type> ConstReference;
    typedef  MultiArray<1,typename AV::value_type>           value_type;
    typedef  MultiArrayView<1,typename AV::value_type>       reference;
    typedef  MultiArrayView<1,typename AV::value_type> const_reference;
    //typedef Value &                                Reference;
    //typedef const Value &                          ConstReference;

    NumpyMultibandEdgeMap()
    :   graph_(NULL),
        array_(){
    }

    NumpyMultibandEdgeMap(const Graph & graph,ArrayView array)
    :   graph_(&graph),
        array_(array){
    }

    Reference operator[](const Key & key){
        return array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicEdgeCoordinate(*graph_,key)];
    }
    ConstReference operator[](const Key & key)const{
        return   array_[GraphDescriptorToMultiArrayIndex<Graph>::intrinsicEdgeCoordinate(*graph_,key)];
    }
private:
    const Graph * graph_;
    mutable AV array_;

};



// TODO ASK UKOETHE FOR HELP HERE
template<unsigned int G_DIM ,class T,class G,unsigned int OG_DIM>
void reshapeNodeMapIfEmpty(
    const G & graph,
    const NumpyArray<OG_DIM,T> & otherArray,
    NumpyArray<G_DIM ,T> & toReshapeArray
){
    const static unsigned int GraphNodeMapDim = IntrinsicGraphShape<G>::IntrinsicNodeMapDimension;
    const static unsigned int OutShapeLength  = GraphNodeMapDim == G_DIM ? G_DIM  : GraphNodeMapDim +1;
    typedef typename MultiArray<OutShapeLength,int>::difference_type OutShapeType;
    OutShapeType outShape;
    for(size_t d=0;d<GraphNodeMapDim;++d){
        outShape[d]=IntrinsicGraphShape<G>::intrinsicNodeMapShape(graph)[d];
    }
    if( G_DIM == GraphNodeMapDim + 1){

        outShape[GraphNodeMapDim]=otherArray.shape(OG_DIM);
        if(GraphNodeMapDim==1)
            toReshapeArray.reshapeIfEmpty( NumpyArray<G_DIM ,T>::ArrayTraits::taggedShape(outShape,"xc"));
        else if(GraphNodeMapDim==2)
            toReshapeArray.reshapeIfEmpty( NumpyArray<G_DIM ,T>::ArrayTraits::taggedShape(outShape,"xyc"));
        else if(GraphNodeMapDim==3)
            toReshapeArray.reshapeIfEmpty( NumpyArray<G_DIM ,T>::ArrayTraits::taggedShape(outShape,"xyzc"));
        else if(GraphNodeMapDim==4)
            toReshapeArray.reshapeIfEmpty( NumpyArray<G_DIM ,T>::ArrayTraits::taggedShape(outShape,"xyztc"));
        else
            throw std::runtime_error("reshapeNodeMapIfEmpty does onnly support graphs with an intrinsic node map shape <=4");
    }
    else{
        toReshapeArray.reshapeIfEmpty(outShape);
    }
}




template<class G,class T>
struct NumpyNodeMap
: 
    IfBool<
        IsMultiband<T>::value,
        NumpyMultibandNodeMap< G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension+1,T> > , 
        NumpyScalarNodeMap<    G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension  ,T> >
    >::type

{
    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension+1,T> , 
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension  ,T>
    >::type NumpyArrayType;


    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyMultibandNodeMap< G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension+1,T> > , 
        NumpyScalarNodeMap<    G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension  ,T> >
    >::type BaseType;

    NumpyNodeMap(const G & g, NumpyArrayType numpyArray)
    :BaseType(g,numpyArray){
    }

};


template<class G,class T>
struct NumpyEdgeMap
: 
    IfBool<
        IsMultiband<T>::value,
        NumpyMultibandEdgeMap< G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension+1,T> > , 
        NumpyScalarEdgeMap<    G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension  ,T> >
    >::type

{
    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension+1,T> , 
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension  ,T>
    >::type NumpyArrayType;


    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyMultibandEdgeMap< G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension+1,T> > , 
        NumpyScalarEdgeMap<    G ,  NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension  ,T> >
    >::type BaseType;

    NumpyEdgeMap(const G & g, NumpyArrayType numpyArray)
    :BaseType(g,numpyArray){
    }

};



template<class G,class T>
struct PyEdgeMapTraits{
    typedef NumpyEdgeMap<G,T> Map;
    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension+1,T> , 
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicEdgeMapDimension  ,T>
    >::type Array;
};




template<class G,class T>
struct PyNodeMapTraits{
    typedef NumpyNodeMap<G,T> Map;
    typedef typename IfBool<
        IsMultiband<T>::value,
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension+1,T> , 
        NumpyArray<IntrinsicGraphShape<G>::IntrinsicNodeMapDimension  ,T>
    >::type Array;
};


} // end namespace vigra




#endif //VIGRA_PYTHON_GRAPH_MAPS_HXX 