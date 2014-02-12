    
#ifndef VIGRA_NEW_MERGE_GRAPH_HXX
#define VIGRA_NEW_MERGE_GRAPH_HXX




/* boost */
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>
//#include <boost/unordered_map.hpp> 

/* std tr1 library */
//#include <tr1/unordered_map>

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

/*parallel std*/
//#include <parallel/algorithm>


/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graphs.hxx>
#include <vigra/graph_helper/dense_map.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>

/* vigra merge graph */
#include "merge_graph_callbacks.hxx"
#include "iterable_partition.hxx"
#include <vigra/random_access_set.hxx>

namespace vigra {

template<class GRAPH,class ITEM>
struct MergeGraphItemHelper;

template<class MG>
struct MergeGraphItemHelper<MG,typename MG::Edge>{
    typedef typename MG::Graph Graph;
    typedef typename MG::index_type index_type ;
    typedef typename MG::Edge Item;
    typedef typename Graph::Edge GraphItem;
    typedef typename MG::EdgeIt ItemIt;


    static index_type maxItemId(const MG & g){
        return g.maxEdgeId();
    }
    static index_type itemNum(const MG & g){
        return g.edgeNum();
    }

    static GraphItem itemToGraphItem(const MG & g,const Item & item){
        const index_type id = g.id(item);
        return g.graph().edgeFromId(id);
    }
};

template<class MG>
struct MergeGraphItemHelper<MG,typename MG::Node>{
    typedef typename MG::Graph Graph;
    typedef typename MG::index_type index_type ;
    typedef typename MG::Node Item;
    typedef typename Graph::Node GraphItem;
    typedef typename MG::NodeIt ItemIt;


    static index_type maxItemId(const MG & g){
        return g.maxNodeId();
    }
    static index_type itemNum(const MG & g){
        return g.nodeNum();
    }
    static GraphItem itemToGraphItem(const MG & g,const Item & item){
        const index_type id = g.id(item);
        return g.graph().nodeFromId(id);
    }
};


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
:   public NewMergeGraphCallbacks<
        detail::GenericNode<vigra::Int64> ,
        detail::GenericEdge<vigra::Int64> 
    > 

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



    

    //typedef  std::set<index_type>   NodeStorageEdgeSet;
    typedef detail::GenericNodeImpl<index_type,false >  NodeStorage;
    typedef detail::GenericEdgeImpl<index_type >        EdgeStorage;



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
        

        template<class T>
        struct EdgeMap : DenseEdgeReferenceMap<MergeGraphType,T> {
            EdgeMap()
            : DenseEdgeReferenceMap<MergeGraphType,T>(){
            }
            EdgeMap(const MergeGraphType & g)
            : DenseEdgeReferenceMap<MergeGraphType,T>(g){
            }
        };


        template<class T>
        struct NodeMap : DenseNodeReferenceMap<MergeGraphType,T> {
            NodeMap()
            : DenseNodeReferenceMap<MergeGraphType,T>(){
            }
            NodeMap(const MergeGraphType & g)
            : DenseNodeReferenceMap<MergeGraphType,T>(g){
            }
        };


        template<class T>
        struct ArcMap : DenseArcReferenceMap<MergeGraphType,T> {
            ArcMap()
            : DenseArcReferenceMap<MergeGraphType,T>(){
            }
            ArcMap(const MergeGraphType & g)
            : DenseArcReferenceMap<MergeGraphType,T>(g){
            }
        };


        
    private:
        MergeGraphAdaptor();                               // non empty-construction
        MergeGraphAdaptor( const MergeGraphAdaptor& other );      // non construction-copyable
        MergeGraphAdaptor& operator=( const MergeGraphAdaptor& ); // non copyable
    public:
        MergeGraphAdaptor(const Graph &  graph);
        //void   setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1);

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
        Arc   arcFromId( const IdType index)const;





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
        void contractEdge(const Edge & edge);


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

    
        // special merge graph members 
        GraphEdge reprGraphEdge(const GraphEdge & edge)const{
            return  graph_.edgeFromId(reprEdgeId(graph_.id(edge)));
        }
        GraphNode reprGraphNode(const GraphNode & node)const{
            return graph_.nodeFromId(reprNodeId(graph_.id(node)));
        }


        Edge reprEdge(const GraphEdge & edge)const{
            return  edgeFromId(reprEdgeId(graph_.id(edge)));
        }
        Node reprNode(const GraphNode & node)const{
            return nodeFromId(reprNodeId(graph_.id(node)));
        }

        const Graph & graph()const{
            return graph_;
        }
        const Graph & graph(){
            return graph_;
        }

        // in which node is a "merged inactive" edge
        Node inactiveEdgesNode(const Edge edge)const{
            return reprNodeId(graphUId(id(edge)));
        }

    private:
        // needs acces to const nodeImpl
        template<class G,class NIMPL,class FILT>
        friend class detail::GenericIncEdgeIt;

        template<class G>
        friend class detail::NeighborNodeFilter;
        template<class G>
        friend class detail::IncEdgeFilter;
        template<class G>
        friend class detail::BackEdgeFilter;
        template<class G>
        friend class detail::IsOutFilter;
        template<class G>
        friend class detail::IsInFilter;
        friend class MergeGraphNodeIt<MergeGraphType>;
        friend class MergeGraphArcIt<MergeGraphType>;
        friend class MergeGraphEdgeIt<MergeGraphType>;

        Edge  edgeFromIdUnsave(const IdType index)const;

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





        //void combineDoubleEdges(const std::vector<IdType> & ,const IdType ,const IdType );
        void searchLocalDoubleEdges(const NodeStorage & node);// , DoubleMap & doubleMap);




        const GRAPH & graph_;
        size_t nInitNodes_;
        size_t nInitEdges_;

        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< NodeStorage >  nodeVector_;

        size_t nDoubleEdges_;
        std::vector<std::pair<index_type,index_type> > doubleEdges_;
};








template<class GRAPH>
MergeGraphAdaptor<GRAPH>::MergeGraphAdaptor(const GRAPH & graph )
:   NewMergeGraphCallbacks<Node,Edge >(),
    graph_(graph),
    nInitNodes_(0),
    nInitEdges_(0),
    nodeUfd_(graph.maxNodeId()+1),
    edgeUfd_(graph.maxEdgeId()+1),
    nodeVector_(graph.maxNodeId()+1),
    nDoubleEdges_(0),
    doubleEdges_(graph_.edgeNum()/2 +1)
{
    for(index_type possibleNodeId = 0 ; possibleNodeId <= graph_.maxNodeId(); ++possibleNodeId){
        if(graph_.nodeFromId(possibleNodeId)==lemon::INVALID){
            nodeUfd_.eraseElement(possibleNodeId);
        }
        else{
            nodeVector_[possibleNodeId].id_ = -1;
            //CGP_ASSERT_OP(nodeUfd_.find(possibleNodeId),==,possibleNodeId);
        }
    }
    
    for(index_type possibleEdgeId = 0 ; possibleEdgeId <= graph_.maxEdgeId(); ++possibleEdgeId){



        const GraphEdge possibleEdge(graph_.edgeFromId(possibleEdgeId));

        
        if(possibleEdge==lemon::INVALID){
            edgeUfd_.eraseElement(possibleEdgeId);
        }
        else{

            //CGP_ASSERT_OP(edgeUfd_.find(possibleEdgeId),==,possibleEdgeId);

            const index_type guid = graphUId(possibleEdgeId);
            const index_type gvid = graphVId(possibleEdgeId);

            //CGP_ASSERT_OP(guid,<=,graph.maxNodeId());
            //CGP_ASSERT_OP(gvid,<=,graph.maxNodeId());

            //CGP_ASSERT_OP(guid,>=,0);
            //CGP_ASSERT_OP(gvid,>=,0);


            nodeVector_[ guid ].insert(gvid,possibleEdgeId);
            nodeVector_[ gvid ].insert(guid,possibleEdgeId);   
        }
    }
    
}



template<class GRAPH>
inline  typename MergeGraphAdaptor<GRAPH>::Edge
MergeGraphAdaptor<GRAPH>::findEdge  (
    const typename MergeGraphAdaptor<GRAPH>::Node & a,
    const typename MergeGraphAdaptor<GRAPH>::Node & b
)const{

    if(a!=b){
        std::pair<index_type,bool> res =  nodeVector_[id(a)].findEdge(id(b));
        if(res.second){
            return Edge(res.first);
        }
    }
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
    typename MergeGraphAdaptor<GRAPH>::Node const & node
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
MergeGraphAdaptor<GRAPH>::edgeFromIdUnsave(
    const typename MergeGraphAdaptor<GRAPH>::IdType index
)const{
    return Edge(index);
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
inline typename MergeGraphAdaptor<GRAPH>::Arc 
MergeGraphAdaptor<GRAPH>::arcFromId(
    const typename MergeGraphAdaptor<GRAPH>::IdType index
)const{
    if(index<=maxEdgeId( ))
        return  Arc(index,index);
    else
        return Arc(index, index-maxEdgeId() -1);
}

template<class GRAPH>
inline bool 
MergeGraphAdaptor<GRAPH>::hasEdgeId(
    const typename MergeGraphAdaptor<GRAPH>::IdType edgeIndex
)const{
    if(edgeIndex<=maxEdgeId() && !edgeUfd_.isErased(edgeIndex)){
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

    return nodeIndex<=maxNodeId() &&  !nodeUfd_.isErased(nodeIndex) && nodeUfd_.find(nodeIndex)==nodeIndex;
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
void MergeGraphAdaptor<GRAPH>::contractEdge(
    const typename MergeGraphAdaptor<GRAPH>::Edge & toDeleteEdge
){
    //std::cout<<"node num "<<nodeNum()<<"\n";
    const index_type toDeleteEdgeIndex = id(toDeleteEdge);
    const index_type nodesIds[2]={id(u(toDeleteEdge)),id(v(toDeleteEdge))};

    // merge the two nodes
    nodeUfd_.merge(nodesIds[0],nodesIds[1]);
    const IdType newNodeRep    = reprNodeId(nodesIds[0]);
    const IdType notNewNodeRep =  (newNodeRep == nodesIds[0] ? nodesIds[1] : nodesIds[0] );

    typename NodeStorage::AdjIt iter=nodeVector_[notNewNodeRep].adjacencyBegin();
    typename NodeStorage::AdjIt end =nodeVector_[notNewNodeRep].adjacencyEnd();
   
    nDoubleEdges_=0;
    for(;iter!=end;++iter){
        const size_t adjToDeadNodeId = iter->nodeId(); 
        if(adjToDeadNodeId!=newNodeRep){

            // REFACTOR ME,  we can make that faster if
            // we do that in set intersect style
            std::pair<index_type,bool> found=nodeVector_[adjToDeadNodeId].findEdge(newNodeRep);


            if(found.second){
                edgeUfd_.merge(iter->edgeId(),found.first);
                
                const index_type edgeA = iter->edgeId();
                const index_type edgeB = found.first;
                const index_type edgeR  = edgeUfd_.find(edgeA);
                const index_type edgeNR = edgeR==edgeA ? edgeB : edgeA; 

                nodeVector_[adjToDeadNodeId].eraseFromAdjacency(notNewNodeRep);

                // refactor me ... this DOES NOT change the key
                nodeVector_[adjToDeadNodeId].eraseFromAdjacency(newNodeRep);
                nodeVector_[adjToDeadNodeId].insert(newNodeRep,edgeR);

                // refactor me .. this DOES NOT change the key
                nodeVector_[newNodeRep].eraseFromAdjacency(adjToDeadNodeId);
                nodeVector_[newNodeRep].insert(adjToDeadNodeId,edgeR);

                doubleEdges_[nDoubleEdges_]=std::pair<index_type,index_type>(edgeR,edgeNR );
                ++nDoubleEdges_;
            }
            else{
                nodeVector_[adjToDeadNodeId].eraseFromAdjacency(notNewNodeRep);
                //nodeVector_[adjToDeadNodeId].eraseFromAdjacency(newNodeRep);
                nodeVector_[adjToDeadNodeId].insert(newNodeRep,iter->edgeId());

                // symetric
                //nodeVector_[newNodeRep].eraseFromAdjacency(adjToDeadNodeId);
                nodeVector_[newNodeRep].insert(adjToDeadNodeId,iter->edgeId());

            }
        }
    }

    //nodeVector_[newNodeRep].merge(nodeVector_[notNewNodeRep]);
    nodeVector_[newNodeRep].eraseFromAdjacency(notNewNodeRep);
    //nodeVector_[newNodeRep].eraseFromAdjacency(newNodeRep); // no self adjacecy
    nodeVector_[notNewNodeRep].clear();
    
    edgeUfd_.eraseElement(toDeleteEdgeIndex);

    //std::cout<<"merge nodes callbacks\n";
    
    this->callMergeNodeCallbacks(Node(newNodeRep),Node(notNewNodeRep));

    //std::cout<<"merge double edge callbacks\n";
    for(size_t de=0;de<nDoubleEdges_;++de){
        this->callMergeEdgeCallbacks(Edge(doubleEdges_[de].first),Edge(doubleEdges_[de].second));
    }
    //std::cout<<"erase edge callbacks\n";
    this->callEraseEdgeCallbacks(Edge(toDeleteEdgeIndex));

    //std::cout<<"and done\n";
}







} // end namespace vigra



#endif //VIGRA_NEW_MERGE_GRAPH_HXX