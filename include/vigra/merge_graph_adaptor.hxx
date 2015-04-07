    
/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
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


#ifndef VIGRA_NEW_MERGE_GRAPH_HXX
#define VIGRA_NEW_MERGE_GRAPH_HXX


/* delegates / callbacks */
#include "delegate/delegate.hxx"

/* std library */
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <map>

/* vigra */
#include "multi_array.hxx"
#include "tinyvector.hxx"
#include "multi_array.hxx"
#include "graphs.hxx"
#include "graph_maps.hxx"
#include "graph_item_impl.hxx"
#include "random_access_set.hxx"
#include "iteratorfacade.hxx"


namespace vigra {

namespace merge_graph_detail {

// ufd  data structure structure for merge graph
// only useful for merge graphs internal usage
template<class T>
class IterablePartition;

// representative element iterator
// for IterablePartition
// only useful for merge graphs internal usage
template<class T>
struct  ConstRepIter
:  public ForwardIteratorFacade<
      ConstRepIter<T>,T,true
   >

{
   typedef IterablePartition<T> IterablePartitionType;
   ConstRepIter(const IterablePartitionType & p,const T cr)
   :  partition_(&p),
      currentRep_(cr){
   }


   ConstRepIter()
   :  partition_(NULL),
      currentRep_()
   {
   }

private:
  friend class vigra::IteratorFacadeCoreAccess;


   bool isBegin()const{
      return partition_!=NULL  && currentRep_==partition_->firstRep();
   }
   bool isEnd()const{
      return  partition_==NULL || currentRep_>partition_->lastRep();
   }

   bool equal(const ConstRepIter & other)const{
      return   (this->isEnd() && other.isEnd() )  || ((this->isEnd()==other.isEnd() ) && this->currentRep_==other.currentRep_);
   }

   void increment(){
      if(partition_->jumpVec_[currentRep_].second==0){
         currentRep_+=1;
      }
      else{
         currentRep_+=partition_->jumpVec_[currentRep_].second;
      }
   }

   void decrement(){
      if(partition_->jumpVec_[currentRep_].first==0){
         //VIGRA_ASSERT_OP(currentRep_,==,partition_->firstRep());
         //currentRep_+=1;
      }
      else{
         currentRep_-=partition_->jumpVec_[currentRep_].first;
      }
   }

   const T & dereference()const{
      return currentRep_;
   }



   const IterablePartitionType * partition_;
   T currentRep_;
   
};



// ufd  data structure structure for merge graph
// only useful for merge graphs internal usage
/// Disjoint set data structure with path compression.
/// \ingroup datastructures
template<class T>
class IterablePartition {
public:
   friend struct ConstRepIter<T>;
   typedef T value_type;
   typedef std::size_t SizeTType;
   IterablePartition();
   IterablePartition(const value_type&);

   // query
   value_type find(const value_type&) const; // without path compression
   value_type find(value_type); // with path compression
   value_type numberOfElements() const;
   value_type numberOfSets() const;
   template<class Iterator> void elementLabeling(Iterator) const;
   template<class Iterator> void representatives(Iterator) const;
   void representativeLabeling(std::map<value_type, value_type>&) const;

   // manipulation
   void reset(const value_type&);
   void merge(value_type, value_type);

   value_type firstRep()const{
      return firstRep_;
   }
   value_type lastRep()const{
      return lastRep_;
   }
   typedef ConstRepIter<T> const_iterator;

   const_iterator begin()const{
      if(numberOfSets_!=0)
         return ConstRepIter<T>(*this,firstRep_);
      else
         return ConstRepIter<T>(*this,lastRep_+1);
   }
   const_iterator end()const{
      return ConstRepIter<T>(*this,lastRep_+1);
   }


   const_iterator iteratorAt(const value_type & rep)const{
      if(numberOfSets_!=0)
         return ConstRepIter<T>(*this,rep);
      else
         return ConstRepIter<T>(*this,lastRep_+1);
   }

   bool isErased(const value_type & value)const{
      return jumpVec_[value].first == -1 && jumpVec_[value].second == -1;
   }

   void eraseElement(const value_type & value,const bool reduceSize=true){
      const T notRep=value;
      const T jumpMinus = jumpVec_[notRep].first;
      const T jumpPlus  = jumpVec_[notRep].second;

      if(jumpMinus==0){
         const T nextRep = notRep+jumpPlus;
         firstRep_=nextRep;
         jumpVec_[nextRep].first=0;
      }
      else if(jumpPlus==0){
         //VIGRA_ASSERT_OP(lastRep_,==,notRep);
         const T prevRep = notRep-jumpMinus;
         lastRep_=prevRep;
         jumpVec_[prevRep].second=0;

      }
      else{
         const T nextRep = notRep+jumpPlus;
         const T prevRep = notRep-jumpMinus;
         jumpVec_[nextRep].first+=jumpVec_[notRep].first;
         jumpVec_[prevRep].second+=jumpVec_[notRep].second;
      }   
      if(reduceSize){
         --numberOfSets_;
      }
      jumpVec_[notRep].first  =-1;
      jumpVec_[notRep].second =-1;
   }

private:
   std::vector<value_type> parents_;
   std::vector<value_type> ranks_;
   std::vector< std::pair< vigra::Int64, vigra::Int64> > jumpVec_;
   value_type firstRep_;
   value_type lastRep_;
   value_type numberOfElements_;
   value_type numberOfSets_;
};


} // end namespa merge graph detail 



// helper classes to generalize
// some functionality for 
// nodes,edges and arcs 
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

// merge graphs LEMON compatible iterator
template<class MERGE_GRAPH>
class MergeGraphNodeIt
:   public ForwardIteratorFacade<MergeGraphNodeIt<MERGE_GRAPH>,typename MERGE_GRAPH::Node,true>{
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
   friend class vigra::IteratorFacadeCoreAccess;
    
    
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

// merge graphs LEMON compatible iterator
template<class MERGE_GRAPH>
class MergeGraphEdgeIt
:   public ForwardIteratorFacade<MergeGraphEdgeIt<MERGE_GRAPH>,typename MERGE_GRAPH::Edge,true>{
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
    friend class vigra::IteratorFacadeCoreAccess;
    
    
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

// merge graphs LEMON compatible iterator
template<class GRAPH>
class MergeGraphArcIt
: public ForwardIteratorFacade<
    MergeGraphArcIt<GRAPH>,typename GRAPH::Arc,true
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
    friend class vigra::IteratorFacadeCoreAccess;
    bool isEnd()const{
        return veryEnd_ || graph_==NULL;
    }

    bool isBegin()const{
        return graph_!=NULL &&  veryEnd_==false && pos_ == EdgeIt(*graph_);         
    }

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


// callbacks of merge graph 
// to update node and edge maps w.r.t. edge contractions
template<class NODE,class EDGE>
class MergeGraphCallbacks{
    public:

        typedef delegate2<void ,const NODE & ,const NODE &>        MergeNodeCallBackType;
        typedef delegate2<void ,const EDGE & ,const EDGE &>        MergeEdgeCallBackType;
        typedef delegate1<void ,const EDGE &>                      EraseEdgeCallBackType;

        MergeGraphCallbacks(){}

        void registerMergeNodeCallBack(MergeNodeCallBackType  f){
            mergeNodeCallbacks_.push_back(f);
        }
        void registerMergeEdgeCallBack(MergeEdgeCallBackType  f){
            mergeEdgeCallbacks_.push_back(f);
        }
        void registerEraseEdgeCallBack(EraseEdgeCallBackType  f){
            eraseEdgeCallbacks_.push_back(f);
        }

    protected:
        void callMergeNodeCallbacks(const NODE & a,const NODE & b){
            for(size_t i=0;i<mergeNodeCallbacks_.size();++i)
                mergeNodeCallbacks_[i](a,b);
        }
        void callMergeEdgeCallbacks(const EDGE & a,const EDGE & b){
            for(size_t i=0;i<mergeEdgeCallbacks_.size();++i)
                mergeEdgeCallbacks_[i](a,b);
        }
        void callEraseEdgeCallbacks(const EDGE & a){
            for(size_t i=0;i<eraseEdgeCallbacks_.size();++i)
                eraseEdgeCallbacks_[i](a);
        }
        void clearCallbacks(){
            mergeNodeCallbacks_.clear();
            mergeEdgeCallbacks_.clear();
            eraseEdgeCallbacks_.clear();
        }
    private:
        std::vector<MergeNodeCallBackType> mergeNodeCallbacks_;
        std::vector<MergeEdgeCallBackType> mergeEdgeCallbacks_;
        std::vector<EraseEdgeCallBackType> eraseEdgeCallbacks_;
};



/** \brief undirected graph adaptor 
      for edge contraction and feature merging
    */
template<class GRAPH>
class MergeGraphAdaptor 
:   public MergeGraphCallbacks<
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
            EdgeMap(): DenseEdgeReferenceMap<MergeGraphType,T>(){
            }
            EdgeMap(const MergeGraphType & g)
            : DenseEdgeReferenceMap<MergeGraphType,T>(g){
            }
            EdgeMap(const MergeGraphType & g,const T & val)
            : DenseEdgeReferenceMap<MergeGraphType,T>(g,val){
            }
        };

        template<class T>
        struct NodeMap : DenseNodeReferenceMap<MergeGraphType,T> {
            NodeMap(): DenseNodeReferenceMap<MergeGraphType,T>(){
            }
            NodeMap(const MergeGraphType & g)
            : DenseNodeReferenceMap<MergeGraphType,T>(g){
            }
            NodeMap(const MergeGraphType & g,const T & val)
            : DenseNodeReferenceMap<MergeGraphType,T>(g,val){
            }
        };

        template<class T>
        struct ArcMap : DenseArcReferenceMap<MergeGraphType,T> {
            ArcMap(): DenseArcReferenceMap<MergeGraphType,T>(){
            }
            ArcMap(const MergeGraphType & g)
            : DenseArcReferenceMap<MergeGraphType,T>(g){
            }
            ArcMap(const MergeGraphType & g,const T & val)
            : DenseArcReferenceMap<MergeGraphType,T>(g,val){
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
        size_t maxDegree()const{
            size_t md=0;
            for(NodeIt it(*this);it!=lemon::INVALID;++it){
                std::max(md, size_t( degree(*it) ) );
            }
            return md;
        }

        void reset();

    private:
        // needs acces to const nodeImpl
        template<class G,class NIMPL,class FILT>
        friend class detail::GenericIncEdgeIt;

        template<class G>
        friend struct detail::NeighborNodeFilter;
        template<class G>
        friend struct detail::IncEdgeFilter;
        template<class G>
        friend struct detail::BackEdgeFilter;
        template<class G>
        friend struct detail::IsOutFilter;
        template<class G>
        friend struct detail::IsInFilter;
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


        const GRAPH & graph_;
        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< NodeStorage >  nodeVector_;

        size_t nDoubleEdges_;
        std::vector<std::pair<index_type,index_type> > doubleEdges_;
};


template<class GRAPH>
MergeGraphAdaptor<GRAPH>::MergeGraphAdaptor(const GRAPH & graph )
:   MergeGraphCallbacks<Node,Edge >(),
    graph_(graph),
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
            nodeVector_[possibleNodeId].id_ = possibleNodeId;
        }
    }
    for(index_type possibleEdgeId = 0 ; possibleEdgeId <= graph_.maxEdgeId(); ++possibleEdgeId){
        const GraphEdge possibleEdge(graph_.edgeFromId(possibleEdgeId));
        if(possibleEdge==lemon::INVALID){
            edgeUfd_.eraseElement(possibleEdgeId);
        }
        else{
            const index_type guid = graphUId(possibleEdgeId);
            const index_type gvid = graphVId(possibleEdgeId);
            nodeVector_[ guid ].insert(gvid,possibleEdgeId);
            nodeVector_[ gvid ].insert(guid,possibleEdgeId);   
        }
    }
    
}


template<class GRAPH>
void MergeGraphAdaptor<GRAPH>::reset  (){

    nodeUfd_.reset(graph_.maxNodeId()+1),
    edgeUfd_.reset(graph_.maxEdgeId()+1),

    this->clearCallbacks();

    // clean nodes_
    for(index_type possibleNodeId = 0 ; possibleNodeId <= graph_.maxNodeId(); ++possibleNodeId){

        nodeVector_[possibleNodeId].clear();
        if(graph_.nodeFromId(possibleNodeId)==lemon::INVALID){
            nodeUfd_.eraseElement(possibleNodeId);
        }
        else{
            nodeVector_[possibleNodeId].id_ = possibleNodeId;
        }
    }

    for(index_type possibleEdgeId = 0 ; possibleEdgeId <= graph_.maxEdgeId(); ++possibleEdgeId){
        const GraphEdge possibleEdge(graph_.edgeFromId(possibleEdgeId));
        if(possibleEdge==lemon::INVALID){
            edgeUfd_.eraseElement(possibleEdgeId);
        }
        else{
            const index_type guid = graphUId(possibleEdgeId);
            const index_type gvid = graphVId(possibleEdgeId);
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


#ifndef DOXYGEN  // doxygen doesn't understand this

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

#endif //DOXYGEN


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
        return Node(lemon::INVALID);
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
inline size_t MergeGraphAdaptor<GRAPH>::nodeNum()const{
    return nodeUfd_.numberOfSets();
}

template<class GRAPH>
inline size_t MergeGraphAdaptor<GRAPH>::arcNum()const{
    return edgeNum()*2;
}

template<class GRAPH>
inline size_t MergeGraphAdaptor<GRAPH>::edgeNum()const{
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



namespace merge_graph_detail {
/// Construct a partition.
template<class T>
IterablePartition<T>::IterablePartition()
: parents_(),
  ranks_(),
  jumpVec_(),
  firstRep_(0),
  lastRep_(0),
  numberOfElements_(0),
  numberOfSets_(0)
{}

/// Construct a partition.
///
/// \param size Number of distinct sets.
///
template<class T>
inline
IterablePartition<T>::IterablePartition
(
   const value_type& size
)
: parents_(static_cast<SizeTType>(size)),
  ranks_(static_cast<SizeTType>(size)),
  jumpVec_(static_cast<SizeTType>(size)),
  firstRep_(0),
  lastRep_(static_cast<SizeTType>(size)-1),
  numberOfElements_(size),
  numberOfSets_(size)
{
   for(T j=0; j<size; ++j) {
      parents_[static_cast<SizeTType>(j)] = j;
   }

   jumpVec_.front().first=0;
   jumpVec_.front().second=1;
   for(T j=1; j<size-1;++j){
      jumpVec_[j].first =1;
      jumpVec_[j].second=1;
   }
   jumpVec_.back().first=1;
   jumpVec_.back().second=0;
}

/// Reset a partition such that each set contains precisely one element
///
/// \param size Number of distinct sets.
///
template<class T>
inline void
IterablePartition<T>::reset
(
   const value_type& size
)
{
   numberOfElements_ = size;
   numberOfSets_ = size;
   ranks_.resize(static_cast<SizeTType>(size));
   parents_.resize(static_cast<SizeTType>(size));
   jumpVec_.resize(static_cast<SizeTType>(size));
   firstRep_=0;
   lastRep_=static_cast<SizeTType>(size)-1;
   for(T j=0; j<size; ++j) {
      ranks_[static_cast<SizeTType>(j)] = 0;
      parents_[static_cast<SizeTType>(j)] = j;
   }

   jumpVec_.front().first=0;
   jumpVec_.front().second=1;
   for(T j=1; j<size-1;++j){
      jumpVec_[j].first =1;
      jumpVec_[j].second=1;
   }
   jumpVec_.back().first=1;
   jumpVec_.back().second=0;
}

/// Find the representative element of the set that contains the given element.
///
/// This constant function does not compress the search path.
///
/// \param element Element.
///
template<class T>
inline typename IterablePartition<T>::value_type
IterablePartition<T>::find
(
   const value_type& element
) const
{
   // find the root
   value_type root = element;
   while(parents_[static_cast<SizeTType>(root)] != root) {
      root = parents_[static_cast<SizeTType>(root)];
   }
   return root;
}

/// Find the representative element of the set that contains the given element.
///
/// This mutable function compresses the search path.
///
/// \param element Element.
///
template<class T>
inline typename IterablePartition<T>::value_type
IterablePartition<T>::find
(
   value_type element // copy to work with
)
{
   // find the root
   value_type root = element;
   while(parents_[static_cast<SizeTType>(root)] != root) {
      root = parents_[static_cast<SizeTType>(root)];
   }
   // path compression
   while(element != root) {
      value_type tmp = parents_[static_cast<SizeTType>(element)];
      parents_[static_cast<SizeTType>(element)] = root;
      element = tmp;
   }
   return root;
}

/// Merge two sets.
///
/// \param element1 Element in the first set.
/// \param element2 Element in the second set.
///
template<class T>
inline void
IterablePartition<T>::merge
(
   value_type element1,
   value_type element2
)
{
   // merge by rank
   element1 = find(element1);
   element2 = find(element2);
   if(element1!=element2){
      T notRep;
      if(ranks_[static_cast<SizeTType>(element1)] < ranks_[static_cast<SizeTType>(element2)]) {
         parents_[static_cast<SizeTType>(element1)] = element2;
         --numberOfSets_;
         //rep=element2;
         notRep=element1;
      }
      else if(ranks_[static_cast<SizeTType>(element1)] > ranks_[static_cast<SizeTType>(element2)]) {
         parents_[static_cast<SizeTType>(element2)] = element1;
         --numberOfSets_;
         //rep=element1;
         notRep=element2;
      }
      else if(element1 != element2) {
         parents_[static_cast<SizeTType>(element2)] = element1;
         ++ranks_[static_cast<SizeTType>(element1)];
         --numberOfSets_;
         //rep=element1;
         notRep=element2;
      }
      this->eraseElement(notRep,false);
   }
}  

template<class T>
inline typename IterablePartition<T>::value_type
IterablePartition<T>::numberOfElements() const
{
   return numberOfElements_;
}

template<class T>
inline typename IterablePartition<T>::value_type
IterablePartition<T>::numberOfSets() const
{
   return numberOfSets_;
}

template<class T>
inline bool operator == (const ConstRepIter<T> & iter,const lemon::Invalid & iv){
    return iter.isEnd();
}
template<class T>
inline bool operator == (const lemon::Invalid & iv , const ConstRepIter<T> & iter){
    return iter.isEnd();
}

template<class T>
inline bool operator != (const ConstRepIter<T> & iter,const lemon::Invalid & iv){
    return !iter.isEnd();
}
template<class T>
inline bool operator != (const lemon::Invalid & iv , const ConstRepIter<T> & iter){
    return !iter.isEnd();
}


} // end namespace merge_graph_detail


} // end namespace vigra



#endif //VIGRA_NEW_MERGE_GRAPH_HXX
