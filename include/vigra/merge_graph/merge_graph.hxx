    
#ifndef VIGRA_MERGE_GRAPH_HXX
#define VIGRA_MERGE_GRAPH_HXX

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
//#include <boost/signals2.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/transform_iterator.hpp>

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>


/* this project*/
#include "merge_graph_node.hxx"
#include "merge_graph_edge.hxx"
#include "merge_graph_iterators.hxx"
#include "merge_graph_callbacks.hxx"
#include "merge_graph_invalid.hxx"

/* this project (TO BE REFACTORED) */
#include "iterable_partition.hxx"
#include "macros.hxx"

namespace vigra {

template<class ID_TYPE>
class MergeGraph 
: 
    public MergeGraphCallbacks<ID_TYPE>  // CALLBACKS
{

    public:
    typedef ID_TYPE                  IdType;
    typedef MergeGraph<ID_TYPE>      MergeGraphType;
    typedef MergeGraphEdge<IdType>   Edge;
    typedef MergeGraphNode<IdType>   Node;
    typedef MergeGraphArc<IdType>    Arc;

    private:
        typedef merge_graph_detail::filter::SmallerThan<IdType>                   BackFilter;
        typedef merge_graph_detail::filter::BackEdgeIdFilter<MergeGraphType>      BackEdgeIdFilter;
        typedef merge_graph_detail::transform::IdToGraphItem<MergeGraphType,Edge> IdToEdgeTransform;
        typedef merge_graph_detail::transform::IdToGraphItem<MergeGraphType,Node> IdToNodeTransform;
        typedef merge_graph_detail::transform::OtherNodeId<MergeGraphType>        OtherNodeIdTransform;
        typedef merge_graph_detail::transform::OtherNode<MergeGraphType>          OtherNodeTransform;
    public:
        
        

        typedef std::map<vigra::UInt64 , std::vector<IdType>  > DoubleMap;
        typedef merge_graph_detail::IterablePartition<IdType> UfdType;
        typedef typename UfdType::const_iterator ConstUdfIter;


        typedef ConstUdfIter                                                                                    EdgeIdIt;
        typedef ConstUdfIter                                                                                    NodeIdIt;
        typedef merge_graph_detail::TransformIter<IdToEdgeTransform,EdgeIdIt,Edge,Edge>                         EdgeIt;
        typedef merge_graph_detail::TransformIter<IdToNodeTransform,NodeIdIt,const Node &,Node >                NodeIt;
  
        typedef typename Node::EdgeIdIt                                                                         NeighborEdgeIdIt;
        typedef merge_graph_detail::TransformIter<IdToEdgeTransform,NeighborEdgeIdIt,Edge,Edge>                 NeighborEdgeIt;
        typedef merge_graph_detail::TransformIter<OtherNodeIdTransform,NeighborEdgeIdIt,IdType,IdType>          NeighborNodeIdIt;
        typedef merge_graph_detail::TransformIter<OtherNodeTransform,NeighborEdgeIdIt,const Node &, Node >      NeighborNodeIt;

        typedef merge_graph_detail::FilterIter< BackFilter,NeighborNodeIdIt>                                    BackNeighborNodeIdIt;
        typedef merge_graph_detail::TransformIter<IdToNodeTransform,BackNeighborNodeIdIt,const Node &,Node >    BackNeighborNodeIt;
        typedef merge_graph_detail::FilterIter< BackFilter,NeighborNodeIdIt>                                    BackNeighborEdgeIdIt;
        typedef merge_graph_detail::TransformIter<IdToEdgeTransform,BackNeighborEdgeIdIt,Edge,Edge>             BackNeighborEdgeIt;

        // redundandet typedefs
        typedef IdType              index_type;
        typedef NeighborNodeIt      adjacency_iterator;
        typedef NeighborNodeIt      neighbor_vertex_iterator;
        typedef BackNeighborNodeIt  back_neighbor_vertex_iterator;
        


    // boost typedefs
    public: 
        typedef boost::directed_tag                 directed_category;
        typedef boost::disallow_parallel_edge_tag   edge_parallel_category;

        typedef Edge edge_descriptor;
        typedef Node vertex_descriptor;

        typedef EdgeIt edge_iterator;
        typedef NodeIt node_iterator;                    

        typedef size_t vertices_size_type;
        typedef size_t edges_size_type;
        typedef size_t degree_size_type;
        
    private:
        MergeGraph();                               // non empty-construction
        MergeGraph( const MergeGraph& other );      // non construction-copyable
        MergeGraph& operator=( const MergeGraph& ); // non copyable
    public:
        MergeGraph(const size_t nNodes,const size_t nEdges);
        void   setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1);

        // query (sizes) 
        size_t numberOfNodes()const;
        size_t numberOfEdges()const;
        size_t numberOfArcs()const;

        size_t arcNum()const;
        size_t edgeNum()const;
        size_t nodeNum()const;

        size_t num_arcs()const;
        size_t num_edges()const;
        size_t num_nodes()const;

        IdType maxArcId()const;
        IdType maxEdgeId()const;
        IdType maxNodeId()const;


        // query (iterators )
        EdgeIdIt  edgeIdsBegin()const;
        EdgeIdIt  edgeIdsEnd()const;
        EdgeIt    edgesBegin()const;
        EdgeIt    edgesEnd()const;
        NodeIdIt  nodeIdsBegin()const;
        NodeIdIt  nodeIdsEnd()const;
        NodeIt    nodesBegin()const;
        NodeIt    nodesEnd()const;


        NeighborEdgeIdIt neigbourEdgeIdsBegin(const Node & node)const{
            return node.edgeIdsBegin();
        }
        NeighborEdgeIdIt neigbourEdgeIdsEnd(const Node & node)const{
            return node.edgeIdsEnd();
        }

        NeighborEdgeIt neigbourEdgesBegin(const Node & node)const{
            return NeighborEdgeIt(node.edgeIdsBegin(),IdToEdgeTransform(*this));
        }
        NeighborEdgeIt neigbourEdgesEnd(const Node & node)const{
            return NeighborEdgeIt(node.edgeIdsEnd(),  IdToEdgeTransform(*this));
        }

        NeighborNodeIdIt neigbourNodeIdsBegin(const Node & node)const{
            return NeighborNodeIdIt(node.edgeIdsBegin(),OtherNodeIdTransform(*this,node.id()));
        }
        NeighborNodeIdIt neigbourNodeIdsEnd(const Node & node)const{
            return NeighborNodeIdIt(node.edgeIdsEnd(),  OtherNodeIdTransform(*this,node.id()));
        }

        NeighborNodeIt neigbourNodesBegin(const Node & node)const{
            return NeighborNodeIt(node.edgeIdsBegin(),OtherNodeTransform(*this,node.id()));
        }
        NeighborNodeIt neigbourNodesEnd(const Node & node)const{
            return NeighborNodeIt(node.edgeIdsEnd(),  OtherNodeTransform(*this,node.id()));
        }

        BackNeighborNodeIdIt backNeigbourNodeIdsBegin(const Node & node)const{
            return BackNeighborNodeIt( BackFilter(node.id()),neigbourNodeIdsBegin(),neigbourNodeIdsEnd() );
        }

        BackNeighborNodeIdIt backNeigbourNodeIdsEnd(const Node & node)const{
            return BackNeighborNodeIt( BackFilter(node.id()),neigbourNodeIdsEnd(),neigbourNodeIdsEnd() );
        }

        BackNeighborNodeIt backNeigbourNodesBegin(const Node & node)const{
            return BackNeighborNodeIt(node.backNeigbourNodeIdsBegin(),OtherNodeTransform(*this,node.id()));
        }
        BackNeighborNodeIt backNeigbourNodesEnd(const Node & node)const{
            return BackNeighborNodeIt(node.backNeigbourNodeIdsEnd(),  OtherNodeTransform(*this,node.id()));
        }

        BackNeighborEdgeIdIt backNeigbourEdgeIdsBegin(const Node & node)const{
            return BackNeighborEdgeIdIt(BackEdgeIdFilter(*this,node.id()),neigbourEdgeIdsBegin(),neigbourEdgeIdsEnd());
        }

        BackNeighborEdgeIdIt backNeigbourEdgeIdsEnd(const Node & node)const{
            return BackNeighborEdgeIdIt(BackEdgeIdFilter(*this,node.id()),neigbourEdgeIdsEnd(),neigbourEdgeIdsEnd());
        }

        BackNeighborEdgeIt backNeigbourEdgesBegin(const Node & node)const{
            return BackNeighborEdgeIt(backNeigbourEdgeIdsBegin(),IdToEdgeTransform(*this));
        }

        BackNeighborEdgeIt backNeigbourEdgesEnd(const Node & node)const{
            return BackNeighborEdgeIt(backNeigbourEdgeIdsEnd(),IdToEdgeTransform(*this));
        }


        neighbor_vertex_iterator get_neighbor_vertex_iterator(const Node & )const;
        neighbor_vertex_iterator get_neighbor_vertex_end_iterator(const Node &)const;
        back_neighbor_vertex_iterator get_back_neighbor_vertex_iterator(const Node & )const;
        back_neighbor_vertex_iterator get_back_neighbor_vertex_end_iterator(const Node &)const;



        // query (inital sizes)
        size_t initNumberOfNodes()const;
        size_t initNumberOfEdges()const;

        //  query (get edge / nodes from id)
        Arc          arcFromId (const IdType index)const;
        Edge         edgeFromId(const IdType index)const;
        const Node & nodeFromId(const IdType index)const;

        // query ( has edge )
        bool hasEdgeId(const IdType edgeIndex)const;
        bool hasNodeId(const IdType nodeIndex)const;

        std::pair<IdType,bool> findEdgeId(const IdType nodeIdA,const IdType nodeIdB)const;
        std::pair<Edge,     bool> findEdge(const IdType nodeIdA,const IdType nodeIdB)const;

        std::pair<IdType,bool> findArcId(const IdType nodeIdA,const IdType nodeIdB)const;
        std::pair<Arc,      bool> findArc(const IdType nodeIdA,const IdType nodeIdB)const;

        IdType id(const Edge & edge)const;
        IdType id(const Node & node)const;


        degree_size_type degree(const Node & node)const;


        const Node & source(const Arc & arc)const;
        const Node & target(const Arc & arc)const;
        const Node & u(const Edge & edge)const;
        const Node & v(const Edge & edge)const;

        // query (w.r.t. inital nodesIds/edgesIds)
        const Edge & getInitalEdge(const IdType index);
        IdType reprEdgeId(const IdType edgeIndex)const;
        IdType reprNodeId(const IdType nodeIndex)const;
        bool stateOfInitalEdge(const IdType initalEdge)const;
        template<class OUT_ITER>
        void stateOfInitalEdges(OUT_ITER begin,OUT_ITER end)const;

        // modification
        void mergeParallelEdges();
        void mergeRegions(const IdType edgeIndex);



  
    private:
        typedef std::map<IdType, Node > NodeMap;
        typedef typename NodeMap::const_iterator ConstNodeMapIterator;


        void combineDoubleEdges(const std::vector<IdType> & ,const IdType ,const IdType );
        void searchLocalDoubleEdges(const Node & node , DoubleMap & doubleMap,const IdType relabelFrom,const IdType relabelTo);

        size_t nInitNodes_;
        size_t nInitEdges_;

        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< Edge >     initEdges_;
        std::vector<Node>  nodeVector_;
};



template<class ID_TYPE>
MergeGraph<ID_TYPE>::MergeGraph(const size_t nNodes,const size_t nEdges)
:   MergeGraphCallbacks<ID_TYPE>(),
    nInitNodes_(nNodes),
    nInitEdges_(nEdges),
    nodeUfd_(nNodes),
    edgeUfd_(nEdges),
    initEdges_(nEdges),
    nodeVector_(nNodes)
{
     // fill nodes
    for(size_t n=0;n<nNodes;++n){
        nodeVector_[n].id_ = n;
        CGP_ASSERT_OP(nodeUfd_.find(n),==,n);
    }
}

template<class ID_TYPE>
inline std::pair<typename MergeGraph<ID_TYPE>::IdType,bool> 
MergeGraph<ID_TYPE>::findEdgeId(
    const IdType nodeIdA,
    const IdType nodeIdB
)const{
   return nodeFromId(nodeIdA).sharedEdge(nodeFromId(nodeIdB));
}

template<class ID_TYPE>
inline std::pair<typename MergeGraph<ID_TYPE>::Edge,     bool> 
MergeGraph<ID_TYPE>::findEdge  (
    const IdType nodeIdA,
    const IdType nodeIdB
)const{
    const std::pair<IdType,bool> result = findEdgeId(nodeIdA,nodeIdB);
    if (result.second)
        return std::pair<Edge,bool>(edgeFromId(result.first),true);
    else
        return std::pair<Edge,bool>(Edge(),false);
}


template<class ID_TYPE>
inline std::pair<typename MergeGraph<ID_TYPE>::IdType,bool> 
MergeGraph<ID_TYPE>::findArcId(
    const IdType nodeIdA,
    const IdType nodeIdB
)const{
    std::pair<IdType,bool> result = findEdgeId(nodeIdA,nodeIdB);

    // check if there is an edge for these nodes
    if(result.second){
        const Edge edge = edgeFromId(result.first);
        // change if other direction
        if(edge[1]==nodeIdA && edge[0]==nodeIdB){
            result.first+=initNumberOfEdges();
        }
    }
    return result;
}

template<class ID_TYPE>
inline std::pair<typename MergeGraph<ID_TYPE>::Arc,     bool> 
MergeGraph<ID_TYPE>::findArc  (
    const IdType nodeIdA,
    const IdType nodeIdB
)const{
    const std::pair<IdType,bool> result = findArcId(nodeIdA,nodeIdB);
    if (result.second)
        return std::pair<Arc,bool>(arcFromId(result.first),true);
    else
        return std::pair<Arc,bool>(Arc(),false);
}


template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Node & 
MergeGraph<ID_TYPE>::source(const Arc & arc)const{
    return nodeFromId(arc[0]);
}

template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Node & 
MergeGraph<ID_TYPE>::target(const Arc & arc)const{
    return nodeFromId(arc[1]);
}

template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Node & 
MergeGraph<ID_TYPE>::u(const Edge & edge)const{
    return nodeFromId(edge[0]);
}

template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Node & 
MergeGraph<ID_TYPE>::v(const Edge & edge)const{
    return nodeFromId(edge[1]);
}
       



template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::maxArcId()const  {
    return edgeUfd_.lastRep()+initNumberOfEdges();
}
template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::maxEdgeId()const {
    return edgeUfd_.lastRep();
}
template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::maxNodeId()const {
    return nodeUfd_.lastRep();
}



template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::id(
    const typename MergeGraph<ID_TYPE>::Edge & edge
)const{
    std::pair<IdType,bool> r = findEdgeId(edge[0],edge[1]);
    if(r.second){
        return r.first;
    }
    else{
        /// TODO
        CGP_ASSERT_OP(false,==,true);
    }
    
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::id(
    const typename MergeGraph<ID_TYPE>::Node & node
)const{
    return node.id();
}
   
template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::degree_size_type 
MergeGraph<ID_TYPE>::degree(
    const MergeGraph<ID_TYPE>::Node & node
)const{
    return static_cast<degree_size_type>(node.numberOfEdges());
}



template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::EdgeIdIt 
MergeGraph<ID_TYPE>::edgeIdsBegin()const{
    return edgeUfd_.begin();
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::EdgeIdIt 
MergeGraph<ID_TYPE>::edgeIdsEnd()const{
    return edgeUfd_.end();
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::EdgeIt 
MergeGraph<ID_TYPE>::edgesBegin()const{
    return EdgeIt(edgeIdsBegin(),IdToEdgeTransform(*this));
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::EdgeIt 
MergeGraph<ID_TYPE>::edgesEnd()const{
    return EdgeIt(edgeIdsEnd(), IdToEdgeTransform(*this));
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::NodeIt 
MergeGraph<ID_TYPE>::nodesBegin()const{
    return NodeIt(nodeIdsBegin(),IdToNodeTransform(*this));
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::NodeIt 
MergeGraph<ID_TYPE>::nodesEnd()const{
    return NodeIt(nodeIdsEnd(),IdToNodeTransform(*this));
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::neighbor_vertex_iterator 
MergeGraph<ID_TYPE>::get_neighbor_vertex_iterator(const Node & node)const{
    return neighbor_vertex_iterator(*this,this->id(node),node,node.edgeIdsBegin());
}


template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::neighbor_vertex_iterator 
MergeGraph<ID_TYPE>::get_neighbor_vertex_end_iterator(const Node & node)const{
    return neighbor_vertex_iterator(*this,this->id(node),node,node.edgeIdsEnd());
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::back_neighbor_vertex_iterator 
MergeGraph<ID_TYPE>::get_back_neighbor_vertex_iterator(const Node & node)const{
    return back_neighbor_vertex_iterator(
        merge_graph_detail::filter::SmallerThan<IdType>(this->id(node)),
        this->get_neighbor_vertex_iterator(),
        this->get_neighbor_vertex_end_iterator()
    );
}


template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::back_neighbor_vertex_iterator 
MergeGraph<ID_TYPE>::get_back_neighbor_vertex_end_iterator(const Node & node)const{
    return back_neighbor_vertex_iterator(
        merge_graph_detail::filter::SmallerThan<IdType>(this->id(node)),
        this->get_neighbor_vertex_end_iterator(),
        this->get_neighbor_vertex_end_iterator()
    );
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::NodeIdIt 
MergeGraph<ID_TYPE>::nodeIdsBegin()const{
    return nodeUfd_.begin();
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::NodeIdIt 
MergeGraph<ID_TYPE>::nodeIdsEnd()const{
    return nodeUfd_.end();
}

template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Edge & 
MergeGraph<ID_TYPE>::getInitalEdge(
    const typename MergeGraph<ID_TYPE>::IdType index
){
    return initEdges_[index];
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::Arc 
MergeGraph<ID_TYPE>::arcFromId(
    const typename MergeGraph<ID_TYPE>::IdType index
)const{

    if(index < initNumberOfEdges()){
        Edge edge = initEdges_[index];
        edge[0]=reprNodeId(edge[0]);
        edge[1]=reprNodeId(edge[1]);
        return Arc(edge[0],edge[1]);
    }
    else{
        Arc arc = initEdges_[index-initNumberOfEdges()];
        Edge edge = initEdges_[index];
        edge[0]=reprNodeId(edge[0]);
        edge[1]=reprNodeId(edge[1]);
        return Arc(edge[1],edge[0]);
    }
}


template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::Edge 
MergeGraph<ID_TYPE>::edgeFromId(
    const typename MergeGraph<ID_TYPE>::IdType index
)const{
    Edge edge = initEdges_[index];
    edge[0]=reprNodeId(edge[0]);
    edge[1]=reprNodeId(edge[1]);
    return edge;
}

template<class ID_TYPE>
inline const typename MergeGraph<ID_TYPE>::Node & 
MergeGraph<ID_TYPE>::nodeFromId(
    const typename MergeGraph<ID_TYPE>::IdType index
)const{
    return nodeVector_[index];
}

template<class ID_TYPE>
inline bool 
MergeGraph<ID_TYPE>::hasEdgeId(
    const typename MergeGraph<ID_TYPE>::IdType edgeIndex
)const{
    const IdType rep = reprEdgeId(edgeIndex);
    if(rep!=edgeIndex){
        return false;
    }
    else{
        const Edge edge=edgeFromId(rep);
        return( edge[0]!=edge[1] );
    }
}

template<class ID_TYPE>
inline bool 
MergeGraph<ID_TYPE>::hasNodeId(
    const typename MergeGraph<ID_TYPE>::IdType nodeIndex
)const{
    return nodeUfd_.find(nodeIndex)==nodeIndex;
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::reprEdgeId(
    const typename MergeGraph<ID_TYPE>::IdType edgeIndex
)const{
    return edgeUfd_.find(edgeIndex);
}

template<class ID_TYPE>
inline typename MergeGraph<ID_TYPE>::IdType 
MergeGraph<ID_TYPE>::reprNodeId(
    const typename MergeGraph<ID_TYPE>::IdType nodeIndex
)const{
    return nodeUfd_.find(nodeIndex);
}

template<class ID_TYPE>
inline bool MergeGraph<ID_TYPE>::stateOfInitalEdge(
    const typename MergeGraph<ID_TYPE>::IdType initalEdge
)const{
    const Edge edge=edgeFromId(initalEdge);
    return (edge[0]!=edge[1]);
}

template<class ID_TYPE>
template<class OUT_ITER>
void MergeGraph<ID_TYPE>::stateOfInitalEdges(
    OUT_ITER begin,
    OUT_ITER end
)const{
    const size_t d = std::distance(begin,end);
    for(size_t ie=0;ie<initNumberOfEdges();++ie){
        const Edge edge=edgeFromId(ie);
        if(edge[0]!=edge[1]){
            begin[ie]=1;
        }
        else{
            begin[ie]=0;
        }
    }
}

template<class ID_TYPE>
inline void MergeGraph<ID_TYPE>::setInitalEdge(
    const size_t initEdge,
    const size_t initNode0,
    const size_t initNode1
){

    // set up inital and dynamic edges
    initEdges_[initEdge].first =initNode0;
    initEdges_[initEdge].second=initNode1;


    // set up the edges of a given region mapping
    nodeVector_[initNode0].edges_.insert(initEdge);
    nodeVector_[initNode1].edges_.insert(initEdge);
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::numberOfNodes()const{
    return nodeUfd_.numberOfSets();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::numberOfEdges()const{
    //return dynamicEdges_.size();
    return edgeUfd_.numberOfSets();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::numberOfArcs()const{
    //return dynamicEdges_.size();
    return numberOfEdges()*2;
}


template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::nodeNum()const{
    return numberOfNodes();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::edgeNum()const{
    return numberOfEdges();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::arcNum()const{
    return numberOfArcs();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::num_nodes()const{
    return numberOfNodes();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::num_edges()const{
    return numberOfEdges();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::num_arcs()const{
    return numberOfArcs();
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::initNumberOfNodes()const{
    return nInitNodes_;
}

template<class ID_TYPE>
inline size_t MergeGraph<ID_TYPE>::initNumberOfEdges()const{
    return nInitEdges_;
}

template<class ID_TYPE>
void MergeGraph<ID_TYPE>::mergeParallelEdges(){
    typedef typename DoubleMap::const_iterator MapIter;
    DoubleMap pEdgeFinder;
    for(size_t e=0;e<nInitEdges_;++e){
        IdType n0=initEdges_[e].first;
        IdType n1=initEdges_[e].second;
        if(n0<n1){
            std::swap(n0,n1);
        }
        const size_t key = n0 + nInitNodes_*n1;
        pEdgeFinder[key].push_back(e);
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

template<class ID_TYPE>
void MergeGraph<ID_TYPE>::combineDoubleEdges(
    const std::vector<ID_TYPE> & toCombine,
    const ID_TYPE r0,
    const ID_TYPE r1
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

template<class ID_TYPE>
void MergeGraph<ID_TYPE>::searchLocalDoubleEdges(
    const MergeGraph<ID_TYPE>::Node & node , 
    MergeGraph<ID_TYPE>::DoubleMap & doubleMap,
    const ID_TYPE relabelFrom,
    const ID_TYPE relabelTo
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
        IdType newNodes[2]={
            oldEdge[0]==relabelFrom ? relabelTo : oldEdge[0] , 
            oldEdge[1]==relabelFrom ? relabelTo : oldEdge[1]
        };
        if(newNodes[1]<newNodes[0]){
            std::swap(newNodes[1],newNodes[0]);
        }
        const size_t  key = newNodes[0] + newNodes[1]*this->initNumberOfNodes();
        doubleMap[key].push_back(outEdgeIndex);
    }
}

template<class ID_TYPE>
void MergeGraph<ID_TYPE>::mergeRegions(
    const typename MergeGraph<ID_TYPE>::IdType toDeleteEdgeIndex
){
    //std::cout<<"merge edge "<<toDeleteEdgeIndex<<"\n";
    const size_t preNumNodes = this->numberOfNodes();

    // assertions that edge is active and
    // its own repr.
    CGP_ASSERT_OP(reprEdgeId(toDeleteEdgeIndex),==,toDeleteEdgeIndex);
    CGP_ASSERT_OP(hasEdgeId(toDeleteEdgeIndex),==,true);

    const Edge toDeleteEdge = edgeFromId(toDeleteEdgeIndex);
    //const size_t nodes[2]= {dynamicEdges_[toDeleteEdgeIndex].first,dynamicEdges_[toDeleteEdgeIndex].second };
    std::vector<size_t> nodes(2);
    nodes[0]=toDeleteEdge[0];
    nodes[1]=toDeleteEdge[1];
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



    const size_t  edgeSizeRep    = nodeVector_[newNodeRep].numberOfEdges();
    const size_t  edgeSizeNotRep = nodeVector_[notNewNodeRep].numberOfEdges();

    // the new region wich is the result of the merge
    Node & newFormedNode = nodeVector_[newNodeRep];

    // merge the edges of the nodes
    newFormedNode.mergeEdges(nodeVector_[notNewNodeRep]);
    CGP_ASSERT_OP(newFormedNode.numberOfEdges(),==,edgeSizeRep+edgeSizeNotRep-1);

    // free old regions edge set (not needed but for consistency)
    nodeVector_[notNewNodeRep].clear();

    // delete the edge which has been between those two regions
    // which we merge (since this edge is the one getting deleted)
    newFormedNode.eraseEdge(toDeleteEdgeIndex);
    //dynamicEdges_.erase(toDeleteEdgeIndex);
    CGP_ASSERT_OP(newFormedNode.numberOfEdges(),==,edgeSizeRep+edgeSizeNotRep-2);


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
                    const size_t adjacentNodeIndex = edgeFromId(toMergeEdgeIndex).otherNodeId(newNodeRep);

                    newFormedNode.eraseAndInsert(toMergeEdgeIndex,newEdgeRep);  
                    nodeVector_[adjacentNodeIndex].eraseAndInsert(toMergeEdgeIndex,newEdgeRep);  
                    
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
    CGP_ASSERT_OP(this->numberOfNodes(),==,preNumNodes-1);
}




} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_HXX