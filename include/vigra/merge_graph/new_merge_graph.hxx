    
#ifndef VIGRA_NEW_MERGE_GRAPH_HXX
#define VIGRA_NEW_MERGE_GRAPH_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* boost */
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>

#include <vigra/graphs.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>







/* this project*/
#include "merge_graph_callbacks.hxx"

/* this project (TO BE REFACTORED) */
#include "iterable_partition.hxx"
#include "macros.hxx"

namespace vigra {


template<class MERGE_GRAPH>
class MergeGraphNodeIt
:   public boost::iterator_facade<MergeGraphNodeIt<MERGE_GRAPH>,const typename MERGE_GRAPH::Node,boost::forward_traversal_tag>{

public:
    typedef MERGE_GRAPH Graph;
    typedef typename Graph::Node Node;
    // Invalid constructor & conversion. 
    MergeGraphNodeIt(const lemon::Invalid & invalid = lemon::INVALID)
    :   graph_(NULL),
        nodeIdIt_(),
        node_(){

    }
    MergeGraphNodeIt(const Graph & g)
    :   graph_(&g),
        nodeIdIt_(g.nodeUfd_.begin()),
        node_(){
    }
    MergeGraphNodeIt(const Graph & g,const Node & node)
    :   graph_(&g),
        nodeIdIt_(g.nodeUfd_.iteratorAt(g.id(node))),
        node_(){

    }
    bool isEnd()const{ 
        return graph_==NULL || nodeIdIt_==graph_->nodeUfd_.end();
    }
    bool isBegin()const{
        return graph_!=NULL && nodeIdIt_==graph_->nodeUfd_.begin();
    }
private:
    friend class boost::iterator_core_access;
    
    
    bool equal(const MergeGraphNodeIt<MERGE_GRAPH> & other)const{
        return (isEnd()&&other.isEnd()) || nodeIdIt_==other.nodeIdIt_;
    }
    void increment(){++nodeIdIt_;}
    const Node & dereference()const{
        node_=Node(*nodeIdIt_);
        return node_;
    }
    // members
    const Graph * graph_;
    typename Graph::NodeIdIt nodeIdIt_;
    mutable Node  node_;
};

template<class MERGE_GRAPH>
class MergeGraphEdgeIt
:   public boost::iterator_facade<MergeGraphEdgeIt<MERGE_GRAPH>,const typename MERGE_GRAPH::Edge,boost::forward_traversal_tag>{

public:
    typedef MERGE_GRAPH Graph;
    typedef typename Graph::Edge Edge;
    // Invalid constructor & conversion. 
    MergeGraphEdgeIt(const lemon::Invalid & invalid = lemon::INVALID)
    :   graph_(NULL),
        edgeIdIt_(),
        edge_(){
    }
    MergeGraphEdgeIt(const Graph & g)
    :   graph_(&g),
        edgeIdIt_(g.edgeUfd_.begin()),
        edge_(){

    }
    MergeGraphEdgeIt(const Graph & g,const Edge & node)
    :   graph_(&g),
        edgeIdIt_(g.edgeUfd_.iteratorAt(g.id(node))),
        edge_(){
    }
    bool isEnd()const{ 
        return graph_==NULL || edgeIdIt_==graph_->edgeUfd_.end();
    }
    bool isBegin()const{
        return graph_!=NULL && edgeIdIt_==graph_->edgeUfd_.begin();
    }
private:
    friend class boost::iterator_core_access;
    
    
    bool equal(const MergeGraphEdgeIt<MERGE_GRAPH> & other)const{
        return (isEnd()&&other.isEnd()) || edgeIdIt_==other.edgeIdIt_;
    }
    void increment(){
        ++edgeIdIt_;
    }
    const Edge & dereference()const{
        edge_=Edge(*edgeIdIt_);
        return edge_;
    }
    // members
    const Graph * graph_;
    typename Graph::EdgeIdIt edgeIdIt_;
    mutable Edge  edge_;
};


template<class GRAPH>
class MergeGraphArcIt
: public boost::iterator_facade<
    MergeGraphArcIt<GRAPH>,
    const typename GRAPH::Arc,
    boost::forward_traversal_tag
>

{
public:
    typedef GRAPH Graph;
    typedef typename  Graph::Arc Arc;
    typedef typename  Graph::Edge Edge;
    typedef typename  Graph::EdgeIt EdgeIt;
    MergeGraphArcIt(const lemon::Invalid invalid = lemon::INVALID )
    :   graph_(NULL),
        pos_(),
        inFirstHalf_(false),
        veryEnd_(true),
        arc_(){
    }
    MergeGraphArcIt(const GRAPH & g )
    :   graph_(&g),
        pos_(g),
        inFirstHalf_(true),
        veryEnd_( g.edgeNum()==0 ? true : false),
        arc_(){
    }

    MergeGraphArcIt(const GRAPH & g , const Arc & arc )
    :   graph_(&g),
        pos_(g,arc.edgeId()),
        inFirstHalf_(g.id(arc)<=g.maxEdgeId()),
        veryEnd_(false),
        arc_(){
    }

private:

    bool isEnd()const{
        return veryEnd_ || graph_==NULL;
    }

    bool isBegin()const{
        return graph_!=NULL &&  veryEnd_==false && pos_ == EdgeIt(*graph_);         
    }


    friend class boost::iterator_core_access;

    void increment() {
        if(inFirstHalf_){
            ++pos_;
            if(pos_ == lemon::INVALID  ) {
                pos_ = EdgeIt(*graph_);
                inFirstHalf_=false;
            }
            return;
        }
        else{
            ++pos_;
            if(pos_ == lemon::INVALID){
                veryEnd_=true;
            }
            return;
        }
    
       
    }
    bool equal(MergeGraphArcIt const& other) const{
        return (
            (
                isEnd()==other.isEnd()                  &&
                inFirstHalf_==other.inFirstHalf_ 
            ) &&
            (isEnd() || graph_==NULL || pos_==other.pos_ )
            );
            
    }

    const Arc & dereference() const { 
        //std::cout<<graph_->id(*pos_)<<"\n";
        arc_ = graph_->direct(*pos_,inFirstHalf_);
        return arc_;
    }


    const GRAPH * graph_;
    EdgeIt pos_;
    bool inFirstHalf_;
    bool veryEnd_;
    mutable Arc arc_;
};



template<class GRAPH>
class MergeGraphAdaptor 
:   
    public MergeGraphCallbacks<vigra::Int64>  // CALLBACKS
{

    public:
    typedef vigra::Int64             IdType;
    typedef IdType                   index_type;
    typedef MergeGraphAdaptor<GRAPH> MergeGraphType;


    typedef detail::GenericNode<index_type>  Node;
    typedef detail::GenericEdge<index_type>  Edge;
    typedef detail::GenericArc<index_type>   Arc;

    typedef GRAPH Graph;
    typedef typename Graph::Node GraphNode;
    typedef typename Graph::Edge GraphEdge;
    typedef typename Graph::Node GraphArc;

    typedef detail::GenericNodeImpl<index_type,std::set<index_type> >  NodeStorage;
    typedef detail::GenericEdgeImpl<index_type >                       EdgeStorage;



    private:
        
        typedef std::map<vigra::UInt64 , std::vector<IdType>  > DoubleMap;
        typedef merge_graph_detail::IterablePartition<IdType> UfdType;
        typedef typename UfdType::const_iterator ConstUdfIter;
        typedef ConstUdfIter                                                EdgeIdIt;
        typedef ConstUdfIter                                                NodeIdIt;
        typedef detail::NeighborNodeFilter<MergeGraphType>                  NnFilter;
        typedef detail::IncEdgeFilter<MergeGraphType>                       IncFilter;
        typedef detail::IsInFilter<MergeGraphType>                          InFlter;
        typedef detail::IsOutFilter<MergeGraphType>                         OutFilter;
    public:
        typedef MergeGraphNodeIt<MergeGraphType>                                 NodeIt;
        typedef MergeGraphEdgeIt<MergeGraphType>                                 EdgeIt;
        typedef MergeGraphArcIt<MergeGraphType>                                  ArcIt;
        typedef detail::GenericIncEdgeIt<MergeGraphType,NodeStorage,NnFilter  >  NeighborNodeIt;
        typedef detail::GenericIncEdgeIt<MergeGraphType,NodeStorage,IncFilter >  IncEdgeIt;
        typedef detail::GenericIncEdgeIt<MergeGraphType,NodeStorage,InFlter   >  InArcIt;
        typedef detail::GenericIncEdgeIt<MergeGraphType,NodeStorage,OutFilter >  OutArcIt;
        

           


        
    private:
        MergeGraphAdaptor();                               // non empty-construction
        MergeGraphAdaptor( const MergeGraphAdaptor& other );      // non construction-copyable
        MergeGraphAdaptor& operator=( const MergeGraphAdaptor& ); // non copyable
    public:
        MergeGraphAdaptor(const Graph &  graph);
        void   setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1);

        // query (sizes) 
        size_t edgeNum()const;
        size_t nodeNum()const;
        size_t arcNum()const;

        IdType maxEdgeId()const;
        IdType maxNodeId()const;
        IdType maxArcId()const;


        // query (iterators )
        EdgeIdIt  edgeIdsBegin()const;
        EdgeIdIt  edgeIdsEnd()const;
        NodeIdIt  nodeIdsBegin()const;
        NodeIdIt  nodeIdsEnd()const;




        //  query (get edge / nodes from id)
        Edge  edgeFromId(const IdType index)const;
        Node  nodeFromId(const IdType index)const;
        Arc   arcFromId( const IdType index)const{
            if(index<=maxEdgeId( ))
                return  Arc(index,index);
            else
                return Arc(index, index-maxEdgeId() -1);
        }

        // query ( has edge )
        bool hasEdgeId(const IdType edgeIndex)const;
        bool hasNodeId(const IdType nodeIndex)const;
        bool hasArcId(const IdType  arcId)const{
            return hasEdgeId(arcFromId(arcId).edgeId());
        }


        Edge findEdge(const Node & a,const Node & b)const;
        Arc  findArc(const Node & u,const Node & v)const;


        IdType id(const Edge & edge)const;
        IdType id(const Node & node)const;
        IdType id(const Arc & arc)const;


        size_t degree(const Node & node)const;



        Node  u(const Edge & edge)const;
        Node  v(const Edge & edge)const;

        Node source(const Arc & arc)const{
            if(arc!=lemon::INVALID)
                return direction(arc) ? u(Edge(arc)) : v(Edge(arc));
            else
                return Node(lemon::INVALID);
        }
        Node target(const Arc & arc)const{
            if(arc!=lemon::INVALID)
                return direction(arc) ? v(Edge(arc)) : u(Edge(arc));
            else
                return Node(lemon::INVALID);
        }


        




        // query (w.r.t. inital nodesIds/edgesIds)
        IdType reprEdgeId(const IdType edgeIndex)const;
        IdType reprNodeId(const IdType nodeIndex)const;
        bool stateOfInitalEdge(const IdType initalEdge)const;
        template<class OUT_ITER>
        void stateOfInitalEdges(OUT_ITER begin,OUT_ITER end)const;

        // modification
        void mergeParallelEdges();
        void mergeRegions(const IdType edgeIndex);


        Node oppositeNode(Node const &n, const Edge &e) const {
            const Node uNode = u(e);
            const Node vNode = v(e);
            if(id(uNode)==id(n)){
                return vNode;
            }
            else if(id(vNode)==id(n)){
                return uNode;
            }
            else{
                return Node(lemon::INVALID);
            }
        }


        Arc direct(const Edge & edge,const bool forward)const{
            if(edge!=lemon::INVALID){
                if(forward)
                    return Arc(id(edge),id(edge));
                else
                    return Arc(id(edge)+(maxEdgeId()+1),id(edge));
            }
            else{
                return Arc(lemon::INVALID);
            }
        }
        Arc direct(const Edge & edge,const Node & node)const{
            if(u(edge)==node)
                return direct(edge,true);
            else if(v(edge)==node)
                return direct(edge,false);
            else
                return Arc(lemon::INVALID);
        }

        bool direction(const Arc & arc)const{
            return arc.id()==arc.edgeId();
        }

  
    private:
        typedef std::map<IdType, Node > NodeMap;
        typedef typename NodeMap::const_iterator ConstNodeMapIterator;

        // needs acces to const nodeImpl
        template<class G,class NIMPL,class FILT>
        friend class detail::GenericIncEdgeIt;


        friend class MergeGraphNodeIt<MergeGraphType>;
        friend class MergeGraphArcIt<MergeGraphType>;
        friend class MergeGraphEdgeIt<MergeGraphType>;


        index_type  uId(const index_type edgeId)const;
        index_type  vId(const index_type edgeId)const;
        index_type  graphUId(const index_type edgeId)const;
        index_type  graphVId(const index_type edgeId)const;
        //index_type  uId(const Edge & edge)const{return uId(id(edge));}
        //index_type  vId(const Edge & edge)const{return vId(id(edge));}
        const NodeStorage & nodeImpl(const Node & node)const{
            return nodeVector_[id(node)];
        }
        NodeStorage & nodeImpl(const Node & node){
            return nodeVector_[id(node)];
        }





        void combineDoubleEdges(const std::vector<IdType> & ,const IdType ,const IdType );
        void searchLocalDoubleEdges(const NodeStorage & node , DoubleMap & doubleMap,const IdType relabelFrom,const IdType relabelTo);




        const GRAPH & graph_;
        size_t nInitNodes_;
        size_t nInitEdges_;

        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< NodeStorage >  nodeVector_;
};








template<class GRAPH>
MergeGraphAdaptor<GRAPH>::MergeGraphAdaptor(const GRAPH & graph )
:   MergeGraphCallbacks<typename MergeGraphAdaptor<GRAPH>::IdType >(),
    graph_(graph),
    nInitNodes_(0),
    nInitEdges_(0),
    nodeUfd_(graph.maxNodeId()+1),
    edgeUfd_(graph.maxEdgeId()+1),
    nodeVector_(graph.maxNodeId()+1)
{
    for(index_type possibleNodeId = 0 ; possibleNodeId <= graph_.maxNodeId(); ++possibleNodeId){
        if(graph_.nodeFromId(possibleNodeId)==lemon::INVALID){
            nodeUfd_.eraseElement(possibleNodeId);
        }
        else{
            nodeVector_[possibleNodeId].id_ = -1;
            CGP_ASSERT_OP(nodeUfd_.find(possibleNodeId),==,possibleNodeId);
        }
    }
    
    for(index_type possibleEdgeId = 0 ; possibleEdgeId <= graph_.maxEdgeId(); ++possibleEdgeId){



        const GraphEdge possibleEdge(graph_.edgeFromId(possibleEdgeId));

        
        if(possibleEdge==lemon::INVALID){
            edgeUfd_.eraseElement(possibleEdgeId);
        }
        else{

            CGP_ASSERT_OP(edgeUfd_.find(possibleEdgeId),==,possibleEdgeId);

            const index_type guid = graphUId(possibleEdgeId);
            const index_type gvid = graphVId(possibleEdgeId);

            CGP_ASSERT_OP(guid,<=,graph.maxNodeId());
            CGP_ASSERT_OP(gvid,<=,graph.maxNodeId());

            CGP_ASSERT_OP(guid,>=,0);
            CGP_ASSERT_OP(gvid,>=,0);


            nodeVector_[ guid ].insertEdgeId(possibleEdgeId);
            nodeVector_[ gvid ].insertEdgeId(possibleEdgeId);   
        }
    }
    
}



template<class GRAPH>
inline  typename MergeGraphAdaptor<GRAPH>::Edge
MergeGraphAdaptor<GRAPH>::findEdge  (
    const typename MergeGraphAdaptor<GRAPH>::Node & nodeA,
    const typename MergeGraphAdaptor<GRAPH>::Node & nodeB
)const{


    const std::pair<IdType,bool> result = nodeImpl(nodeA).sharedEdge(nodeImpl(nodeB));
    if (result.second)
        return Edge(result.first);
    else
        return Edge(lemon::INVALID);
}

template<class GRAPH>
inline  typename MergeGraphAdaptor<GRAPH>::Arc
MergeGraphAdaptor<GRAPH>::findArc  (
    const typename MergeGraphAdaptor<GRAPH>::Node & uNode,
    const typename MergeGraphAdaptor<GRAPH>::Node & vNode
)const{
    const Edge edge = findEdge(uNode,vNode);
    if(edge==lemon::INVALID)
        return Arc(lemon::INVALID);
    else
        return  direct(edge,u(edge)==uNode);
}


template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::Node 
MergeGraphAdaptor<GRAPH>::u(const Edge & edge)const{
    return nodeFromId(uId(id(edge)));
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::Node 
MergeGraphAdaptor<GRAPH>::v(const Edge & edge)const{
    return nodeFromId(vId(id(edge)));
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::index_type 
MergeGraphAdaptor<GRAPH>::uId(const index_type edgeId)const{
    return reprNodeId(graphUId(edgeId));
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::index_type 
MergeGraphAdaptor<GRAPH>::vId(const index_type edgeId)const{
    return reprNodeId(graphVId(edgeId));
}



template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::index_type 
MergeGraphAdaptor<GRAPH>::graphUId(const index_type edgeId)const{
    return graph_.id(graph_.u(graph_.edgeFromId(edgeId)));
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::index_type 
MergeGraphAdaptor<GRAPH>::graphVId(const index_type edgeId)const{
    return graph_.id(graph_.v(graph_.edgeFromId(edgeId)));
}


template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::maxEdgeId()const {
    return static_cast<index_type>(edgeUfd_.lastRep());
}
template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::maxNodeId()const {
    return static_cast<index_type>(nodeUfd_.lastRep());
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::maxArcId()const {
    return maxEdgeId()*2 +1 ;
}


template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::id(
    const typename MergeGraphAdaptor<GRAPH>::Edge & edge
)const{
    return edge.id();
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::id(
    const typename MergeGraphAdaptor<GRAPH>::Node & node
)const{
    return node.id();
}
   
template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::id(
    const typename MergeGraphAdaptor<GRAPH>::Arc & arc
)const{
    return arc.id();
}
 


template<class GRAPH>
inline size_t 
MergeGraphAdaptor<GRAPH>::degree(
    const MergeGraphAdaptor<GRAPH>::Node & node
)const{
    return static_cast<size_t>( nodeVector_[id(node)].edgeNum() );
}



template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::EdgeIdIt 
MergeGraphAdaptor<GRAPH>::edgeIdsBegin()const{
    return edgeUfd_.begin();
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::EdgeIdIt 
MergeGraphAdaptor<GRAPH>::edgeIdsEnd()const{
    return edgeUfd_.end();
}


template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::NodeIdIt 
MergeGraphAdaptor<GRAPH>::nodeIdsBegin()const{
    return nodeUfd_.begin();
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::NodeIdIt 
MergeGraphAdaptor<GRAPH>::nodeIdsEnd()const{
    return nodeUfd_.end();
}


template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::Edge 
MergeGraphAdaptor<GRAPH>::edgeFromId(
    const typename MergeGraphAdaptor<GRAPH>::IdType index
)const{
    if (hasEdgeId(index))
        return Edge(index);
    else
        return Edge(lemon::INVALID);
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::Node 
MergeGraphAdaptor<GRAPH>::nodeFromId(
    const typename MergeGraphAdaptor<GRAPH>::IdType index
)const{
    if(hasNodeId(index))
        return Node(index);
    else
        return Node(-1);
}

template<class GRAPH>
inline bool 
MergeGraphAdaptor<GRAPH>::hasEdgeId(
    const typename MergeGraphAdaptor<GRAPH>::IdType edgeIndex
)const{
    if(!edgeUfd_.isErased(edgeIndex)){
        const IdType reprEdgeIndex = reprEdgeId(edgeIndex);
        if(reprEdgeIndex!=edgeIndex){
            return false;
        }
        else{
            const index_type rnid0=  uId(reprEdgeIndex);
            const index_type rnid1=  vId(reprEdgeIndex);
            return rnid0!=rnid1;
        }
    }
    else{
        return false;
    }
}

template<class GRAPH>
inline bool 
MergeGraphAdaptor<GRAPH>::hasNodeId(
    const typename MergeGraphAdaptor<GRAPH>::IdType nodeIndex
)const{
    return !nodeUfd_.isErased(nodeIndex) && nodeUfd_.find(nodeIndex)==nodeIndex;
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::reprEdgeId(
    const typename MergeGraphAdaptor<GRAPH>::IdType edgeIndex
)const{
    return edgeUfd_.find(edgeIndex);
}

template<class GRAPH>
inline typename MergeGraphAdaptor<GRAPH>::IdType 
MergeGraphAdaptor<GRAPH>::reprNodeId(
    const typename MergeGraphAdaptor<GRAPH>::IdType nodeIndex
)const{
    return nodeUfd_.find(nodeIndex);
}

template<class GRAPH>
inline bool MergeGraphAdaptor<GRAPH>::stateOfInitalEdge(
    const typename MergeGraphAdaptor<GRAPH>::IdType initalEdge
)const{
    const index_type rep = reprEdgeId(initalEdge);

    const index_type rnid0=  reprNodeId( graphUId(initalEdge) );
    const index_type rnid1=  reprNodeId( graphVId(initalEdge) );
    return rnid0!=rnid1;
}

template<class GRAPH>
template<class OUT_ITER>
void MergeGraphAdaptor<GRAPH>::stateOfInitalEdges(
    OUT_ITER begin,
    OUT_ITER end
)const{
    
    throw std::runtime_error("REFACTORE ME");
}

template<class GRAPH>
inline void MergeGraphAdaptor<GRAPH>::setInitalEdge(
    const size_t initEdge,
    const size_t initNode0,
    const size_t initNode1
){
    // set up the edges of a given region mapping
    nodeVector_[initNode0].edges_.insert(initEdge);
    nodeVector_[initNode1].edges_.insert(initEdge);
}

template<class GRAPH>
inline size_t MergeGraphAdaptor<GRAPH>::nodeNum()const{
    return nodeUfd_.numberOfSets();
}



template<class GRAPH>
inline size_t MergeGraphAdaptor<GRAPH>::arcNum()const{
    return edgeNum()*2;
}

template<class GRAPH>
inline size_t MergeGraphAdaptor<GRAPH>::edgeNum()const{
    //return dynamicEdges_.size();
    return edgeUfd_.numberOfSets();
}


template<class GRAPH>
void MergeGraphAdaptor<GRAPH>::mergeParallelEdges(){
    typedef typename DoubleMap::const_iterator MapIter;
    DoubleMap pEdgeFinder;





    for(EdgeIdIt eIter = edgeUfd_.begin() ;eIter!=edgeUfd_.end()  ;++eIter){

        const GraphEdge graphEdge(*eIter);
        IdType n0=graphVId(graphEdge);
        IdType n1=graphUId(graphEdge);

        if(n0<n1){
            std::swap(n0,n1);
        }
        const size_t key = n0 + nInitNodes_*n1;
        pEdgeFinder[key].push_back(graph_.id(graphEdge));
    }

    for(MapIter iter=pEdgeFinder.begin();iter!=pEdgeFinder.end();++iter){
        const std::vector<IdType> & dEdges = iter->second;
        CGP_ASSERT_OP(dEdges.size(),!=,0);

        if(dEdges.size()>1){
            //std::cout<<"found double edges "<<dEdges.size()<<"\n";
            const size_t key = iter->first;
            const size_t r1  = key/nInitNodes_;
            const size_t r0  = key - nInitNodes_*r1;
            this->combineDoubleEdges(dEdges,r0,r1);
        }
    }
}

template<class GRAPH>
void MergeGraphAdaptor<GRAPH>::combineDoubleEdges(
    const std::vector<typename MergeGraphAdaptor<GRAPH>::IdType > & toCombine,
    const typename MergeGraphAdaptor<GRAPH>::IdType  r0,
    const typename MergeGraphAdaptor<GRAPH>::IdType  r1
)
{
    //std::set<IdType> toCombineSet(toCombine.begin(),toCombine.end());
    //CGP_ASSERT_OP(toCombine.size(),==,toCombineSet.size());

    const IdType newIndex = edgeUfd_.multiMerge(toCombine.front(),toCombine.begin()+1,toCombine.end());
   

    // update the two region between the double edge 
    const IdType regions[2]={r0,r1};
    for(size_t r=0;r<2;++r){
        const size_t ri=regions[r];
        std::set<IdType> & nodesEdges = nodeVector_[ri].edges_;
        for(size_t i=0;i<toCombine.size();++i){
            if(toCombine[i]!=newIndex){
                const bool found = static_cast<bool>(nodesEdges.find(toCombine[i])!=nodesEdges.end());
                CGP_ASSERT_OP(found,==,true);
                const size_t nErased = nodesEdges.erase(toCombine[i]);
                CGP_ASSERT_OP(nErased,==,1);
            }
        }
    }

    // call the registerd callbacks to merge the edges
    for(size_t i=0;i<toCombine.size();++i){
        if(toCombine[i]!=newIndex){
            this->callMergeEdgeCallbacks(newIndex,toCombine[i]);
        }
    }


    //CGP_ASSERT_OP(dynamicEdges_.size(),==,edgeUfd_.numberOfSets());
}

template<class GRAPH>
void MergeGraphAdaptor<GRAPH>::searchLocalDoubleEdges(
    const MergeGraphAdaptor<GRAPH>::NodeStorage & node , 
    MergeGraphAdaptor<GRAPH>::DoubleMap & doubleMap,
    const typename MergeGraphAdaptor<GRAPH>::IdType relabelFrom,
    const typename MergeGraphAdaptor<GRAPH>::IdType relabelTo
){
    // loop over all edges of the new formed region
    for(
        typename std::set<IdType>::const_iterator  edgeIter = node.edges_.begin();
        edgeIter!=node.edges_.end();
        ++edgeIter
    ){
        const IdType outEdgeIndex = *edgeIter;
        //*edgeIter = reprEdgeId(outEdgeIndex);
        //CGP_ASSERT_OP(outEdgeIndex,!=,edgeIndex);

        const Edge  oldEdge    = this->edgeFromId(reprEdgeId(outEdgeIndex));
        //const IdType oldNodes[2]= {dynamicEdges_[outEdgeIndex].first,dynamicEdges_[outEdgeIndex].second };
        // do the relabling 
        //IdType newNodes[2]={
        //    id(u(oldEdge)) ==relabelFrom ? relabelTo : id(u(oldEdge)) , 
        //    id(v(oldEdge)) ==relabelFrom ? relabelTo : id(v(oldEdge)) , 
        //};

        IdType newNodes[2]={
            uId(id(oldEdge)),vId(id(oldEdge))
        };


        CGP_ASSERT_OP(newNodes[0],!=,newNodes[1]);

        
        if(newNodes[1]<newNodes[0]){
            std::swap(newNodes[1],newNodes[0]);
        }
        const size_t  key = newNodes[0] + newNodes[1]*(this->maxNodeId()+1);
        doubleMap[key].push_back(outEdgeIndex);
    }
}

template<class GRAPH>
void MergeGraphAdaptor<GRAPH>::mergeRegions(
    const typename MergeGraphAdaptor<GRAPH>::IdType toDeleteEdgeIndex
){
    //std::cout<<"merge edge "<<toDeleteEdgeIndex<<"\n";
    const size_t preNumNodes = this->nodeNum();

    // assertions that edge is active and
    // its own repr.
    CGP_ASSERT_OP(reprEdgeId(toDeleteEdgeIndex),==,toDeleteEdgeIndex);
    CGP_ASSERT_OP(hasEdgeId(toDeleteEdgeIndex),==,true);

    const Edge toDeleteEdge = edgeFromId(toDeleteEdgeIndex);
    //const size_t nodes[2]= {dynamicEdges_[toDeleteEdgeIndex].first,dynamicEdges_[toDeleteEdgeIndex].second };
    std::vector<size_t> nodes(2);
    nodes[0]=id(u(toDeleteEdge));
    nodes[1]=id(v(toDeleteEdge));
    CGP_ASSERT_OP(nodes[0],!=,nodes[1]);

    for(size_t n=0;n<2;++n){
        // assertions that node is active and
        // its own repr.
        const size_t  ni=nodes[n];
        CGP_ASSERT_OP(reprNodeId(ni),==,ni);
        CGP_ASSERT_OP(hasNodeId(ni),==,true);
    }


    // merge the two nodes
    nodeUfd_.merge(nodes[0],nodes[1]);
    const IdType newNodeRep    = reprNodeId(nodes[0]);
    const IdType notNewNodeRep =  (newNodeRep == nodes[0] ? nodes[1] : nodes[0] );



    const size_t  edgeSizeRep    = nodeVector_[newNodeRep].edgeNum();
    const size_t  edgeSizeNotRep = nodeVector_[notNewNodeRep].edgeNum();

    // the new region wich is the result of the merge
    NodeStorage & newFormedNode = nodeVector_[newNodeRep];

    // merge the edges of the nodes
    newFormedNode.mergeEdges(nodeVector_[notNewNodeRep]);
    CGP_ASSERT_OP(newFormedNode.edgeNum(),==,edgeSizeRep+edgeSizeNotRep-1);

    // free old regions edge set (not needed but for consistency)
    nodeVector_[notNewNodeRep].clear();

    // delete the edge which has been between those two regions
    // which we merge (since this edge is the one getting deleted)
    newFormedNode.eraseEdge(toDeleteEdgeIndex);
    //dynamicEdges_.erase(toDeleteEdgeIndex);
    CGP_ASSERT_OP(newFormedNode.edgeNum(),==,edgeSizeRep+edgeSizeNotRep-2);


    // bevore processing with merging the edges we call the "merge" of the node maps
    // - we need to do this bevore any "merge" within the nodeMaps such that
    //   we can guarantee that the nodes maps are tidy when the edge-maps mergers
    //   are called
    this->callMergeNodeCallbacks(newNodeRep,notNewNodeRep);


    edgeUfd_.eraseElement(toDeleteEdgeIndex);

    // construct the "DoubleMap"
    // - if an vector in the map has a size >=2 
    //   this means that there are multiple edges
    //   between a pair of regions which needs to be merged
    DoubleMap doubleEdgeMap;
    CGP_ASSERT_OP(notNewNodeRep ,!= , newNodeRep);



    this->searchLocalDoubleEdges(newFormedNode,doubleEdgeMap,notNewNodeRep,newNodeRep);

    // loop over the double map
    // if an vector in the map has a size >=2 
    // this means that there are multiple edges
    // between a pair of regions which needs to be merged
    for( typename DoubleMap::const_iterator dIter = doubleEdgeMap.begin();dIter!=doubleEdgeMap.end();++dIter){

        // if this vector has a size >=2 this means we have multiple
        // edges between 2 regions
        // the 2 regions are encoded in the key (dIter->first)
        // but we do not need them here
        const std::vector<IdType> & edgeVec = dIter->second;
        if(edgeVec.size()>=2){
            CGP_ASSERT_OP(edgeVec.size(),==,2);
            // merge all these edges in the ufd and get the new representative
            //CGP_ASSERT_OP(hasEdgeId(toMergeEdgeIndex),==,true);
            const IdType newEdgeRep = edgeUfd_.multiMerge(edgeVec.front(),edgeVec.begin()+1,edgeVec.end());
            //CGP_ASSERT_OP(hasEdgeId(toMergeEdgeIndex),==,false);
            // delte all edges which are not needed any more
            //  - edgeVec.size() -1 edges will be deleted 
            //  - (all edges except the new representative "newEdgeRep")
            // furthermore  the edge-sets all nodes adjacent to the "newFormedNode"
            // must be visited since they might refere to nodes which are deleted /merged
            
            

            for(size_t td=0;td<edgeVec.size();++td){

                // index of the edge which is considered for deletion
                const IdType toMergeEdgeIndex = edgeVec[td];
                // delte this edge only if it is NOT the new representative edge
                if(toMergeEdgeIndex!=newEdgeRep){

                    // delete the edge from the new formed region
                    newFormedNode.edges_.erase(toMergeEdgeIndex);

                    //  not true any more
                    //CGP_ASSERT_OP(hasEdgeId(toMergeEdgeIndex),==,true);
                    

                    // at least one of the nodes of the edge "toMergeEdgeIndex" must be the "newFormedNode"
                    //  - we want to get the nodes adjacent to the "newFormedNode"
                    //CGP_ASSERT_OP(dynamicEdges_[toMergeEdgeIndex].hasNodeId(newNodeRep),==,true);
                    //const size_t adjacentNodeIndex = edgeFromId(toMergeEdgeIndex).otherNodeId(newNodeRep);

                    const index_type nodeUId = uId(toMergeEdgeIndex);
                    const index_type nodeVId = vId(toMergeEdgeIndex);


                    //const size_t adjacentNodeIndex = id(this->oppositeNode(nodeFromId(newNodeRep),edgeFromId(toMergeEdgeIndex)));  
                    const size_t adjacentNodeIndex = nodeUId == newNodeRep ? nodeVId : nodeUId ;
                    newFormedNode.eraseAndInsert(toMergeEdgeIndex,newEdgeRep);  

                    if(nodeVector_[adjacentNodeIndex].hasEdgeId(toMergeEdgeIndex)){
                        nodeVector_[adjacentNodeIndex].eraseAndInsert(toMergeEdgeIndex,newEdgeRep);  
                    }

                    
                    
                    // finaly delete the unneeded edge
                    //dynamicEdges_.erase(toMergeEdgeIndex);
                    //CGP_ASSERT_OP(hasEdge_OLD(toMergeEdgeIndex),==,false);
                }
            }
            CGP_ASSERT_OP(edgeVec.size(),==,2)
            CGP_ASSERT_OP(edgeVec[0],!=,toDeleteEdgeIndex);
            CGP_ASSERT_OP(edgeVec[1],!=,toDeleteEdgeIndex);
            CGP_ASSERT_OP(edgeVec[0],!=,edgeVec[1]);
            CGP_ASSERT_OP(hasEdgeId(newEdgeRep),==,true);
            // CALL CALLBACKS TO MERGE EDGES
            this->callMergeEdgeCallbacks(newEdgeRep, (newEdgeRep==edgeVec[0] ? edgeVec[1] : edgeVec[0]));
        } 
    }

    // CALL CALLBACKS TO ERASE EDGE
    this->callEraseEdgeCallbacks(toDeleteEdgeIndex);

    CGP_ASSERT_OP(nodeUfd_.numberOfSets(),==,preNumNodes-1);
    CGP_ASSERT_OP(this->nodeNum(),==,preNumNodes-1);
}




} // end namespace vigra



#endif //VIGRA_NEW_MERGE_GRAPH_HXX