/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 Stefan Schmidt and Ullrich Koethe            */
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

/**
 * This header provides definitions of graph-related types
 * and optionally provides a gateway to popular graph libraries
 * (for now, BGL is supported).
 */

#ifndef VIGRA_PYTHON_GRAPH_HXX
#define VIGRA_PYTHON_GRAPH_HXX

/*boost*/
#include <boost/python.hpp>
#include <boost/iterator/transform_iterator.hpp>

/*vigra*/
#include <vigra/graphs.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graphs.hxx>
#include <vigra/priority_queue.hxx>
#include <vigra/merge_graph_adaptor.hxx>
namespace vigra{






template<class MAP>
struct GraphMapTypeTraits;

template< unsigned int DIM,class T>
struct GraphMapTypeTraits<NumpyArray<DIM,T> >{
    typedef typename NumpyArray<DIM,T>::value_type Value;
    typedef Value &                                Reference;
    typedef const Value  &                         ConstReference;
};


    
template<class GRAPH>
struct NodeHolder :  GRAPH::Node
{
    typedef typename GRAPH::Node Node;
    NodeHolder(const lemon::Invalid & iv = lemon::INVALID)
    : Node(lemon::INVALID),
      graph_(NULL)
    {}
    NodeHolder(const GRAPH & g , const Node & item)
    : Node(item),
      graph_(&g)
    {}

    typename GRAPH::index_type id()const{
        return graph_->id(*this);
    }

    typename GraphDescriptorToMultiArrayIndex<GRAPH>::IntrinsicNodeMapShape
    intrinsicNodeCoordinate()const{
        return GraphDescriptorToMultiArrayIndex<GRAPH>::intrinsicNodeCoordinate(*graph_,*this);
    }

    const GRAPH * graph_;
};



template<class GRAPH>
struct EdgeHolder : GRAPH::Edge
{

    typedef typename GRAPH::Edge Edge;
    EdgeHolder(const lemon::Invalid & iv = lemon::INVALID)
    : Edge(lemon::INVALID),
      graph_(NULL)
    {}
    EdgeHolder(const GRAPH & g , const Edge & item)
    : Edge(item),
      graph_(&g)
    {}

    typename GRAPH::index_type id()const{
        return graph_->id(*this);
    }

    NodeHolder<GRAPH> u()const{
        return NodeHolder<GRAPH>(*graph_,graph_->u(*this));
    }
    NodeHolder<GRAPH> v()const{
        return NodeHolder<GRAPH>(*graph_,graph_->v(*this));
    }

    typename GraphDescriptorToMultiArrayIndex<GRAPH>::IntrinsicEdgeMapShape
    intrinsicEdgeCoordinate()const{
        return GraphDescriptorToMultiArrayIndex<GRAPH>::intrinsicEdgeCoordinate(*graph_,*this);
    }

    const GRAPH * graph_; 
};



template<class GRAPH>
struct ArcHolder: GRAPH::Arc {
    typedef typename GRAPH::Arc Arc;
    ArcHolder(const lemon::Invalid & iv = lemon::INVALID)
    : Arc(lemon::INVALID),
      graph_(NULL)
    {}
    ArcHolder(const GRAPH & g , const Arc & item)
    : Arc(item),
      graph_(&g)
    {}

    typename GRAPH::index_type id()const{
        return graph_->id(*this);
    }

    typename GraphDescriptorToMultiArrayIndex<GRAPH>::IntrinsicArcMapShape
    intrinsicArcCoordinate()const{
        return GraphDescriptorToMultiArrayIndex<GRAPH>::intrinsicArcCoordinate(*graph_,*this);
    }


    const GRAPH * graph_;
};


namespace detail_python_graph{

template<class GRAPH>
struct ArcToTargetNodeHolder{
    typedef typename GRAPH::Node Node;
    typedef typename GRAPH::Arc Arc;
    ArcToTargetNodeHolder(const GRAPH & graph)
    : graph_(&graph){
    }
    NodeHolder<GRAPH> operator()(const Arc & arc)const{
        return NodeHolder<GRAPH>(*graph_,graph_->target(arc));
    }
    const GRAPH * graph_;
};

template<class GRAPH>
struct ArcToEdgeHolder{
    typedef typename GRAPH::Edge Edge;
    typedef typename GRAPH::Arc Arc;
    ArcToEdgeHolder(const GRAPH & graph)
    : graph_(&graph){
    }
    EdgeHolder<GRAPH> operator()(const Arc & arc)const{
        const Edge edge(arc);
        return EdgeHolder<GRAPH>(*graph_,edge);
    }
    const GRAPH * graph_;
};

template<class GRAPH>
struct NodeToNodeHolder{
    typedef typename GRAPH::Node Node;
    NodeToNodeHolder(const GRAPH & graph)
    : graph_(&graph){
    }
    NodeHolder<GRAPH> operator()(const Node & node)const{
        return NodeHolder<GRAPH>(*graph_,node);
    }
    const GRAPH * graph_;
};

template<class GRAPH>
struct EdgeToEdgeHolder{
    typedef typename GRAPH::Edge Edge;
    EdgeToEdgeHolder(const GRAPH & graph)
    : graph_(&graph){
    }
    EdgeHolder<GRAPH> operator()(const Edge & edge)const{
        return EdgeHolder<GRAPH>(*graph_,edge);
    }
    const GRAPH * graph_;
};

} // end namespace detail_python_graph



template<class GRAPH>
struct NodeIteratorHolder{
    typedef typename GRAPH::Node Node;
    typedef typename GRAPH::NodeIt Iter;
    typedef detail_python_graph::NodeToNodeHolder<GRAPH> Transform;
    typedef boost::transform_iterator<Transform ,Iter ,NodeHolder<GRAPH>, NodeHolder<GRAPH> > const_iterator;
    NodeIteratorHolder(const GRAPH & graph,const Node & node = Node(lemon::INVALID) )
    : graph_(&graph),
      node_(node){
    }
    const_iterator begin()const{

        Iter iter = GraphIteratorAccessor<GRAPH>::nodesBegin(*graph_);
        return const_iterator(iter,Transform(*graph_));
    }
    const_iterator end()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::nodesEnd(*graph_);
        return const_iterator(iter,Transform(*graph_));
    }
    const GRAPH * graph_;
    Node node_;
};

template<class GRAPH>
struct EdgeIteratorHolder{
    typedef typename GRAPH::Edge Edge;
    typedef typename GRAPH::EdgeIt Iter;
    typedef detail_python_graph::EdgeToEdgeHolder<GRAPH> Transform;
    typedef boost::transform_iterator<Transform ,Iter ,EdgeHolder<GRAPH>, EdgeHolder<GRAPH> > const_iterator;
    EdgeIteratorHolder(const GRAPH & graph,const Edge & edge = Edge(lemon::INVALID) )
    : graph_(&graph),
      edge_(edge){
    }
    const_iterator begin()const{

        Iter iter = GraphIteratorAccessor<GRAPH>::edgesBegin(*graph_);
        return const_iterator(iter,Transform(*graph_));
    }
    const_iterator end()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::edgesEnd(*graph_);
        return const_iterator(iter,Transform(*graph_));
    }
    const GRAPH * graph_;
    Edge edge_;
};


template<class GRAPH>
struct NeighbourNodeIteratorHolder{
    typedef typename GRAPH::Node Node;
    typedef typename GRAPH::OutArcIt Iter;
    typedef detail_python_graph::ArcToTargetNodeHolder<GRAPH> Transform;
    typedef boost::transform_iterator<Transform ,Iter ,NodeHolder<GRAPH>, NodeHolder<GRAPH> > const_iterator;
    NeighbourNodeIteratorHolder(const GRAPH & graph,const Node & node)
    : graph_(&graph),
      node_(node){
    }
    const_iterator begin()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::outArcBegin(*graph_,node_);
        return const_iterator(iter,Transform(*graph_));
    }
    const_iterator end()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::outArcEnd(*graph_,node_);
        return const_iterator(iter,Transform(*graph_));
    }
    const GRAPH * graph_;
    Node node_;
};


template<class GRAPH>
struct IncEdgeIteratorHolder{
    typedef typename GRAPH::Node Node;
    typedef typename GRAPH::Edge Edge;
    typedef typename GRAPH::OutArcIt Iter;
    typedef detail_python_graph::ArcToEdgeHolder<GRAPH> Transform;
    typedef boost::transform_iterator<Transform ,Iter ,EdgeHolder<GRAPH>, EdgeHolder<GRAPH> > const_iterator;
    IncEdgeIteratorHolder(const GRAPH & graph,const Node & node)
    : graph_(&graph),
      node_(node){
    }
    const_iterator begin()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::outArcBegin(*graph_,node_);
        return const_iterator(iter,Transform(*graph_));
    }
    const_iterator end()const{
        Iter iter = GraphIteratorAccessor<GRAPH>::outArcEnd(*graph_,node_);
        return const_iterator(iter,Transform(*graph_));
    }
    const GRAPH * graph_;
    Node node_;
};

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
    bool any()const{
        return array_.any();
    }
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
    bool any()const{
        return array_.any();
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
    bool any()const{
        return array_.any();
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
    bool any()const{
        return array_.any();
    }
private:
    const Graph * graph_;
    mutable AV array_;

};






// tagged shape for lemon graphs
// edge map / node map / arc map 
template<class G>
class TaggedGraphShape{
public:
    typedef G Graph;
    const static unsigned int ND = IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension;
    const static unsigned int ED = IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension;
    const static unsigned int AD = IntrinsicGraphShape<Graph>::IntrinsicArcMapDimension;
    static TaggedShape  taggedNodeMapShape(const Graph & graph){
        return NumpyArray<ND,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(graph),"n");
    }
    static TaggedShape  taggedEdgeMapShape(const Graph & graph){
        return NumpyArray<ED,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(graph),"e");
    }
    static TaggedShape  taggedArcMapShape(const Graph & graph){
        return NumpyArray<AD,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicArcMapShape(graph),"e");
    }

    static AxisInfo  axistagsNodeMap(const Graph & graph){
<<<<<<< HEAD
        return AxisInfo("n");
    }
    static AxisInfo  axistagsEdgeMap(const Graph & graph){
        return AxisInfo("e");
    }
    static AxisTags  axistagsArcMap(const Graph & graph){
        return AxisInfo("e");
=======
       return AxisInfo("n");
    }
    static AxisInfo  axistagsEdgeMap(const Graph & graph){
       return AxisInfo("e");
    }
    static AxisInfo  axistagsArcMap(const Graph & graph){
       return AxisInfo("e");
>>>>>>> c5dcff5ade5acabe897b222b7f51deec44a8e166
    }
};

// macro to specialize TaggedGraphShape for 
// grid graphs up to 4 dimensions
#define VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO(DIM,tn,te,ta) \
template<class BOOST_DIRECTED_TAG> \
class TaggedGraphShape<GridGraph<DIM,BOOST_DIRECTED_TAG> >{ \
public: \
    typedef GridGraph<DIM,BOOST_DIRECTED_TAG> Graph; \
    const static unsigned int ND = IntrinsicGraphShape<Graph>::IntrinsicNodeMapDimension; \
    const static unsigned int ED = IntrinsicGraphShape<Graph>::IntrinsicEdgeMapDimension; \
    const static unsigned int AD = IntrinsicGraphShape<Graph>::IntrinsicArcMapDimension; \
    static TaggedShape  taggedNodeMapShape(const Graph & graph){ \
       return NumpyArray<ND,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicNodeMapShape(graph),tn); \
    } \
    static TaggedShape  taggedEdgeMapShape(const Graph & graph){  \
       return NumpyArray<ED,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicEdgeMapShape(graph),te);  \
    } \
    static TaggedShape  taggedArcMapShape(const Graph & graph){  \
       return NumpyArray<AD,int>::ArrayTraits::taggedShape(IntrinsicGraphShape<Graph>::intrinsicArcMapShape(graph),ta);  \
    } \
    static AxisInfo  axistagsNodeMap(const Graph & graph){ \
<<<<<<< HEAD
        return AxisInfo(tn); \
    } \
    static AxisInfo  axistagsEdgeMap(const Graph & graph){ \
        return AxisInfo(te); \
    } \
    static AxisTags  axistagsArcMap(const Graph & graph){ \
        return AxisInfo(ta); \
=======
       return AxisInfo(tn); \
    } \
    static AxisInfo  axistagsEdgeMap(const Graph & graph){ \
       return AxisInfo(te); \
    } \
    static AxisInfo  axistagsArcMap(const Graph & graph){ \
       return AxisInfo(ta); \
>>>>>>> c5dcff5ade5acabe897b222b7f51deec44a8e166
    } \
};

VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO(1,"x","xe","xe");
VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO(2,"xy","xye","xye");
VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO(3,"xyz","xyze","xyze");
VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO(4,"xyzt","xyzte","xyzte");

#undef VIGRA_MAKE_TAGGED_GRAPH_SHAPE_MACRO



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


namespace cluster_operators{

template<class MERGE_GRAPH>
class PythonOperator{
    
    typedef PythonOperator<MERGE_GRAPH > SelfType;
public:


    typedef float WeightType;
    typedef MERGE_GRAPH MergeGraph;
    typedef typename MergeGraph::Graph Graph;
    typedef typename Graph::Edge GraphEdge;
    typedef typename Graph::Node GraphNode;
    typedef typename MergeGraph::Edge Edge;
    typedef typename MergeGraph::Node Node;
    typedef typename MergeGraph::EdgeIt EdgeIt;
    typedef typename MergeGraph::NodeIt NodeIt;
    typedef typename MergeGraph::IncEdgeIt IncEdgeIt;
    typedef typename MergeGraph::index_type index_type;
    typedef MergeGraphItemHelper<MergeGraph,Edge> EdgeHelper;
    typedef MergeGraphItemHelper<MergeGraph,Node> NodeHelper;


    typedef NodeHolder<MERGE_GRAPH> NodeHolderType;
    typedef EdgeHolder<MERGE_GRAPH> EdgeHolderType;

    PythonOperator(
        MergeGraph & mergeGraph,
        boost::python::object object,
        const bool useMergeNodeCallback,
        const bool useMergeEdgesCallback,
        const bool useEraseEdgeCallback
    )
    :   mergeGraph_(mergeGraph),
        object_(object)
    {
        if(useMergeNodeCallback)
            mergeGraph_.registerMergeNodeCallBack(*this,& SelfType::mergeNodes);
        if(useMergeEdgesCallback)
            mergeGraph_.registerMergeEdgeCallBack(*this,& SelfType::mergeEdges);
        if(useEraseEdgeCallback)
            mergeGraph_.registerEraseEdgeCallBack(*this,& SelfType::eraseEdge);

    }

    void mergeEdges(const Edge & a,const Edge & b){
        const EdgeHolderType aa(mergeGraph_,a);
        const EdgeHolderType bb(mergeGraph_,b);
        object_.attr("mergeEdges")(aa,bb);
    }
    void mergeNodes(const Node & a,const Node & b){
        const NodeHolderType aa(mergeGraph_,a);
        const NodeHolderType bb(mergeGraph_,b);
        object_.attr("mergeNodes")(aa,bb);
    }
    void eraseEdge(const Edge & e){
        const EdgeHolderType ee(mergeGraph_,e);
        object_.attr("eraseEdge")(ee);
    }
    Edge contractionEdge(){
        EdgeHolderType eh = boost::python::extract<EdgeHolderType>(object_.attr("contractionEdge")());
        return eh;
    }
    WeightType contractionWeight()const{
        return boost::python::extract<WeightType>(object_.attr("contractionWeight")());
    }

    MergeGraph & mergeGraph(){
        return mergeGraph_;
    }
private:
    MergeGraph & mergeGraph_;
    boost::python::object object_;
};

} // end namespace cluster_operators


} // namespace vigra

#endif // VIGRA_PYTHON_GRAPH_HXX
