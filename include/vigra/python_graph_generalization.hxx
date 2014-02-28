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

#ifndef VIGRA_PYTHON_GRAPH_GENERALIZATION_HXX
#define VIGRA_PYTHON_GRAPH_GENERALIZATION_HXX

/*boost*/
#include <boost/iterator/transform_iterator.hpp>

/*vigra*/
#include <vigra/graphs.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/python_graph_maps.hxx>
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


namespace detail_python_graph_generalization{

template<class GRAPH>
struct ArcToTargetNodeHolder{
    typedef typename GRAPH::Node Node;
    typedef typename GRAPH::Arc  Arc;
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
    typedef typename GRAPH::Arc  Arc;
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

}


template<class GRAPH>
struct NodeIteratorHolder{
    typedef typename GRAPH::Node    Node;
    typedef typename GRAPH::NodeIt  Iter;
    typedef detail_python_graph_generalization::NodeToNodeHolder<GRAPH> Transform;
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
    typedef typename GRAPH::Edge    Edge;
    typedef typename GRAPH::EdgeIt  Iter;
    typedef detail_python_graph_generalization::EdgeToEdgeHolder<GRAPH> Transform;
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
    typedef detail_python_graph_generalization::ArcToTargetNodeHolder<GRAPH> Transform;
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
    typedef detail_python_graph_generalization::ArcToEdgeHolder<GRAPH> Transform;
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




} // namespace vigra

#endif // VIGRA_PYTHON_GRAPH_GENERALIZATION_HXX
