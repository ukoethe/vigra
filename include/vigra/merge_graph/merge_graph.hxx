    
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

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/is_end_transform_iterator.hxx>
#include <vigra/is_end_filter_iterator.hxx>
#include <vigra/graph_iterator_functors.hxx>
#include <vigra/graph_api_helper.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>
#include <vigra/graph_helper/graph_crtp_base.hxx>


/* this project*/
#include "merge_graph_callbacks.hxx"

/* this project (TO BE REFACTORED) */
#include "iterable_partition.hxx"
#include "macros.hxx"

namespace vigra {

template<class ID_TYPE>
class MergeGraphImpl 
:   
    public MergeGraphCallbacks<ID_TYPE>  // CALLBACKS
{

    public:
    typedef ID_TYPE                  IdType;
    typedef IdType                   index_type;
    typedef MergeGraphImpl<ID_TYPE>  MergeGraphType;

    /*
    typedef MergeGraphEdge<IdType>   Edge;
    typedef MergeGraphNode<IdType>   Node;
    */

    typedef detail::GenericNode<index_type>  Node;
    typedef detail::GenericEdge<index_type>  Edge;


    typedef detail::GenericNodeImpl<index_type,std::set<index_type> >  NodeStorageType;
    typedef detail::GenericEdgeImpl<index_type >                       EdgeStorageType;


    private:
        typedef detail::filter::SmallerThan<IdType>                             BackFilter;
        typedef detail::filter::BackEdgeIdFilter<MergeGraphType>                BackEdgeIdFilter;
        typedef detail::transform::IdToItem<IdType,Edge>                        IdToEdgeTransform;
        typedef detail::transform::IdToItem<IdType,Node  >                      IdToNodeTransform;
        typedef detail::transform::OtherNodeId<MergeGraphType>                  OtherNodeIdTransform;
        typedef detail::transform::OtherNode<MergeGraphType>                    OtherNodeTransform;
    public:
        
        

        typedef std::map<vigra::UInt64 , std::vector<IdType>  > DoubleMap;
        typedef merge_graph_detail::IterablePartition<IdType> UfdType;
        typedef typename UfdType::const_iterator ConstUdfIter;


        typedef ConstUdfIter                                                                EdgeIdIt;
        typedef ConstUdfIter                                                                NodeIdIt;
        typedef TransformIter<IdToEdgeTransform,EdgeIdIt,Edge,Edge>                         EdgeIt;
        typedef TransformIter<IdToNodeTransform,NodeIdIt,Node,Node >                        NodeIt;
  
        typedef typename NodeStorageType::EdgeIdIt                                          NeighborEdgeIdIt;
        typedef TransformIter<IdToEdgeTransform,NeighborEdgeIdIt,Edge,Edge>                 NeighborEdgeIt;
        typedef TransformIter<OtherNodeIdTransform,NeighborEdgeIdIt,IdType,IdType>          NeighborNodeIdIt;
        typedef TransformIter<OtherNodeTransform,NeighborEdgeIdIt, Node , Node >            NeighborNodeIt;

        typedef FilterIter< BackFilter,NeighborNodeIdIt>                                    BackNeighborNodeIdIt;
        typedef TransformIter<IdToNodeTransform,BackNeighborNodeIdIt,Node,Node >            BackNeighborNodeIt;
        typedef FilterIter< BackFilter,NeighborNodeIdIt>                                    BackNeighborEdgeIdIt;
        typedef TransformIter<IdToEdgeTransform,BackNeighborEdgeIdIt,Edge,Edge>             BackNeighborEdgeIt;

        // redundandet typedefs
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
        MergeGraphImpl();                               // non empty-construction
        MergeGraphImpl( const MergeGraphImpl& other );      // non construction-copyable
        MergeGraphImpl& operator=( const MergeGraphImpl& ); // non copyable
    public:
        MergeGraphImpl(const size_t nNodes,const size_t nEdges);
        void   setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1);

        // query (sizes) 
        size_t numberOfNodes()const;
        size_t numberOfEdges()const;

        size_t edgeNum()const;
        size_t nodeNum()const;


        size_t num_edges()const;
        size_t num_nodes()const;

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
            return nodeVector_[id(node)].edgeIdsBegin();
        }
        NeighborEdgeIdIt neigbourEdgeIdsEnd(const Node & node)const{
            return nodeVector_[id(node)].edgeIdsEnd();
        }

        // rename me
        NeighborEdgeIt neigbourEdgesBegin(const Node & node)const{
            return NeighborEdgeIt(neigbourEdgeIdsBegin(node),IdToEdgeTransform());
        }
        NeighborEdgeIt neigbourEdgesEnd(const Node & node)const{
            return NeighborEdgeIt(neigbourEdgeIdsEnd(node),  IdToEdgeTransform());
        }



        NeighborEdgeIt incEdgeItBegin(const Node & node)const{
            return NeighborEdgeIt(neigbourEdgeIdsBegin(node),IdToEdgeTransform());
        }
        NeighborEdgeIt incEdgeIt(const Node & node)const{
            return NeighborEdgeIt(neigbourEdgeIdsEnd(node),  IdToEdgeTransform());
        }



        NeighborNodeIdIt neigbourNodeIdsBegin(const Node & node)const{
            return NeighborNodeIdIt(neigbourEdgeIdsBegin(node),OtherNodeIdTransform(*this,node.id()));
        }
        NeighborNodeIdIt neigbourNodeIdsEnd(const Node & node)const{
            return NeighborNodeIdIt(neigbourEdgeIdsEnd(node),  OtherNodeIdTransform(*this,node.id()));
        }

        NeighborNodeIt neigbourNodesBegin(const Node & node)const{
            return NeighborNodeIt(neigbourEdgeIdsBegin(node),OtherNodeTransform(*this,node.id()));
        }
        NeighborNodeIt neigbourNodesEnd(const Node & node)const{
            return NeighborNodeIt(neigbourEdgeIdsEnd(node),  OtherNodeTransform(*this,node.id()));
        }

        BackNeighborNodeIdIt backNeigbourNodeIdsBegin(const Node & node)const{
            return BackNeighborNodeIt( BackFilter(node.id()),neigbourNodeIdsBegin(),neigbourNodeIdsEnd() );
        }

        BackNeighborNodeIdIt backNeigbourNodeIdsEnd(const Node & node)const{
            return BackNeighborNodeIt( BackFilter(node.id()),neigbourNodeIdsEnd(),neigbourNodeIdsEnd() );
        }

        BackNeighborNodeIt backNeigbourNodesBegin(const Node & node)const{
            return BackNeighborNodeIt(backNeigbourNodeIdsBegin(),OtherNodeTransform(*this,node.id()));
        }
        BackNeighborNodeIt backNeigbourNodesEnd(const Node & node)const{
            return BackNeighborNodeIt(backNeigbourNodeIdsEnd(),  OtherNodeTransform(*this,node.id()));
        }

        BackNeighborEdgeIdIt backNeigbourEdgeIdsBegin(const Node & node)const{
            return BackNeighborEdgeIdIt(BackEdgeIdFilter(*this,node.id()),neigbourEdgeIdsBegin(),neigbourEdgeIdsEnd());
        }

        BackNeighborEdgeIdIt backNeigbourEdgeIdsEnd(const Node & node)const{
            return BackNeighborEdgeIdIt(BackEdgeIdFilter(*this,node.id()),neigbourEdgeIdsEnd(),neigbourEdgeIdsEnd());
        }

        BackNeighborEdgeIt backNeigbourEdgesBegin(const Node & node)const{
            return BackNeighborEdgeIt(backNeigbourEdgeIdsBegin(),IdToEdgeTransform());
        }

        BackNeighborEdgeIt backNeigbourEdgesEnd(const Node & node)const{
            return BackNeighborEdgeIt(backNeigbourEdgeIdsEnd(),IdToEdgeTransform());
        }


        neighbor_vertex_iterator get_neighbor_vertex_iterator(const Node & )const;
        neighbor_vertex_iterator get_neighbor_vertex_end_iterator(const Node &)const;
        back_neighbor_vertex_iterator get_back_neighbor_vertex_iterator(const Node & )const;
        back_neighbor_vertex_iterator get_back_neighbor_vertex_end_iterator(const Node &)const;



        // query (inital sizes)
        size_t initNumberOfNodes()const;
        size_t initNumberOfEdges()const;

        //  query (get edge / nodes from id)
        Edge  edgeFromId(const IdType index)const;
        Node  nodeFromId(const IdType index)const;

        // query ( has edge )
        bool hasEdgeId(const IdType edgeIndex)const;
        bool hasNodeId(const IdType nodeIndex)const;

        std::pair<IdType,bool> findEdgeId(const IdType nodeIdA,const IdType nodeIdB)const;
        std::pair<Edge,     bool> findEdge(const IdType nodeIdA,const IdType nodeIdB)const;



        IdType id(const Edge & edge)const;
        IdType id(const Node & node)const;



        degree_size_type degree(const Node & node)const;



        Node  u(const Edge & edge)const;
        Node  v(const Edge & edge)const;

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
                return Node(-1);
            }
        }
  
    private:
        typedef std::map<IdType, Node > NodeMap;
        typedef typename NodeMap::const_iterator ConstNodeMapIterator;


        void combineDoubleEdges(const std::vector<IdType> & ,const IdType ,const IdType );
        void searchLocalDoubleEdges(const NodeStorageType & node , DoubleMap & doubleMap,const IdType relabelFrom,const IdType relabelTo);

        size_t nInitNodes_;
        size_t nInitEdges_;

        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< EdgeStorageType >  initEdges_;
        std::vector< NodeStorageType >  nodeVector_;
};


template<class ID_TYPE>
class MergeGraph 
:   public MergeGraphImpl<ID_TYPE>,
    public detail::ArcHelper
    < 
        MergeGraph<ID_TYPE> ,                               // GRAPH
        typename MergeGraphImpl<ID_TYPE>::index_type,       // INDEX_TYPE
        typename MergeGraphImpl<ID_TYPE>::Edge ,            // EDGE
        typename MergeGraphImpl<ID_TYPE>::Node ,            // NODE
        typename MergeGraphImpl<ID_TYPE>::EdgeIt,           // EDGE_IT 
        typename MergeGraphImpl<ID_TYPE>::NodeIt,           // NODE_IT 
        typename MergeGraphImpl<ID_TYPE>::NeighborEdgeIdIt  // NEIGHBOR_EDGE_ID_IT
    >
{
public:
    typedef MergeGraphImpl<ID_TYPE> MergeGraphImplType;
    typedef detail::ArcHelper
    < 
        MergeGraph<ID_TYPE> ,                               // GRAPH
        typename MergeGraphImpl<ID_TYPE>::index_type,       // INDEX_TYPE
        typename MergeGraphImpl<ID_TYPE>::Edge ,            // EDGE
        typename MergeGraphImpl<ID_TYPE>::Node ,            // NODE
        typename MergeGraphImpl<ID_TYPE>::EdgeIt,           // EDGE_IT 
        typename MergeGraphImpl<ID_TYPE>::NodeIt,           // NODE_IT 
        typename MergeGraphImpl<ID_TYPE>::NeighborEdgeIdIt  // NEIGHBOR_EDGE_ID_IT
    > ArcHelperType;
    MergeGraph(const size_t nNodes,const size_t nEdges)
    : MergeGraphImpl<ID_TYPE>(nNodes,nEdges){
    }

    typename MergeGraphImplType::index_type id(const typename ArcHelperType::Arc & arc)const{
        return this->arcId(arc);
    }

    typename MergeGraphImplType::index_type id(const typename MergeGraphImplType::Node & node)const{
        return static_cast<const MergeGraphImplType *>(this)->id(node);
    }

    typename MergeGraphImplType::index_type id(const typename MergeGraphImplType::Edge & edge)const{
        return static_cast<const MergeGraphImplType *>(this)->id(edge);
    }
};



template<class ID_TYPE>
MergeGraphImpl<ID_TYPE>::MergeGraphImpl(const size_t nNodes,const size_t nEdges)
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
inline std::pair<typename MergeGraphImpl<ID_TYPE>::IdType,bool> 
MergeGraphImpl<ID_TYPE>::findEdgeId(
    const IdType nodeIdA,
    const IdType nodeIdB
)const{
   return nodeFromId(nodeIdA).sharedEdge(nodeFromId(nodeIdB));
}

template<class ID_TYPE>
inline std::pair<typename MergeGraphImpl<ID_TYPE>::Edge,     bool> 
MergeGraphImpl<ID_TYPE>::findEdge  (
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
inline typename MergeGraphImpl<ID_TYPE>::Node 
MergeGraphImpl<ID_TYPE>::u(const Edge & edge)const{
    const index_type rnid = reprNodeId(initEdges_[id(edge)].u());
    return nodeFromId(rnid);
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::Node  
MergeGraphImpl<ID_TYPE>::v(const Edge & edge)const{
    const index_type rnid = reprNodeId(initEdges_[id(edge)].v());
    return nodeFromId(rnid);
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::maxEdgeId()const {
    return edgeUfd_.lastRep();
}
template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::maxNodeId()const {
    return nodeUfd_.lastRep();
}



template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::id(
    const typename MergeGraphImpl<ID_TYPE>::Edge & edge
)const{
    return edge.id();
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::id(
    const typename MergeGraphImpl<ID_TYPE>::Node & node
)const{
    return node.id();
}
   
template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::degree_size_type 
MergeGraphImpl<ID_TYPE>::degree(
    const MergeGraphImpl<ID_TYPE>::Node & node
)const{
    return static_cast<degree_size_type>( nodeVector_[id(node)].edgeNum() );
}



template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::EdgeIdIt 
MergeGraphImpl<ID_TYPE>::edgeIdsBegin()const{
    return edgeUfd_.begin();
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::EdgeIdIt 
MergeGraphImpl<ID_TYPE>::edgeIdsEnd()const{
    return edgeUfd_.end();
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::EdgeIt 
MergeGraphImpl<ID_TYPE>::edgesBegin()const{
    return EdgeIt(edgeIdsBegin(),IdToEdgeTransform());
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::EdgeIt 
MergeGraphImpl<ID_TYPE>::edgesEnd()const{
    return EdgeIt(edgeIdsEnd(), IdToEdgeTransform());
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::NodeIt 
MergeGraphImpl<ID_TYPE>::nodesBegin()const{
    return NodeIt(nodeIdsBegin(),IdToNodeTransform());
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::NodeIt 
MergeGraphImpl<ID_TYPE>::nodesEnd()const{
    return NodeIt(nodeIdsEnd(),IdToNodeTransform());
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::neighbor_vertex_iterator 
MergeGraphImpl<ID_TYPE>::get_neighbor_vertex_iterator(const Node & node)const{
    return neighbor_vertex_iterator(*this,this->id(node),node,node.edgeIdsBegin());
}


template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::neighbor_vertex_iterator 
MergeGraphImpl<ID_TYPE>::get_neighbor_vertex_end_iterator(const Node & node)const{
    return neighbor_vertex_iterator(*this,this->id(node),node,node.edgeIdsEnd());
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::back_neighbor_vertex_iterator 
MergeGraphImpl<ID_TYPE>::get_back_neighbor_vertex_iterator(const Node & node)const{
    return back_neighbor_vertex_iterator(
        BackFilter(this->id(node)),
        this->get_neighbor_vertex_iterator(),
        this->get_neighbor_vertex_end_iterator()
    );
}


template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::back_neighbor_vertex_iterator 
MergeGraphImpl<ID_TYPE>::get_back_neighbor_vertex_end_iterator(const Node & node)const{
    return back_neighbor_vertex_iterator(
        BackFilter(this->id(node)),
        this->get_neighbor_vertex_end_iterator(),
        this->get_neighbor_vertex_end_iterator()
    );
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::NodeIdIt 
MergeGraphImpl<ID_TYPE>::nodeIdsBegin()const{
    return nodeUfd_.begin();
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::NodeIdIt 
MergeGraphImpl<ID_TYPE>::nodeIdsEnd()const{
    return nodeUfd_.end();
}

template<class ID_TYPE>
inline const typename MergeGraphImpl<ID_TYPE>::Edge & 
MergeGraphImpl<ID_TYPE>::getInitalEdge(
    const typename MergeGraphImpl<ID_TYPE>::IdType index
){
    return initEdges_[index];
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::Edge 
MergeGraphImpl<ID_TYPE>::edgeFromId(
    const typename MergeGraphImpl<ID_TYPE>::IdType index
)const{
    /*
    Edge edge = initEdges_[index];
    edge[0]=reprNodeId(edge[0]);
    edge[1]=reprNodeId(edge[1]);
    return edge;
    */
    if (hasEdgeId(index))
        return Edge(index);
    else
        return Edge(-1);
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::Node 
MergeGraphImpl<ID_TYPE>::nodeFromId(
    const typename MergeGraphImpl<ID_TYPE>::IdType index
)const{
    if(hasNodeId(index))
        return Node(index);
    else
        return Node(-1);
}

template<class ID_TYPE>
inline bool 
MergeGraphImpl<ID_TYPE>::hasEdgeId(
    const typename MergeGraphImpl<ID_TYPE>::IdType edgeIndex
)const{
    const IdType reprEdgeIndex = reprEdgeId(edgeIndex);
    if(reprEdgeIndex!=edgeIndex){
        return false;
    }
    else{

        const Node n0 = initEdges_[reprEdgeIndex].u();
        const Node n1 = initEdges_[reprEdgeIndex].v();
        const index_type rnid0=  reprNodeId(id(n0));
        const index_type rnid1=  reprNodeId(id(n1));
        return rnid0!=rnid1;
    }
}

template<class ID_TYPE>
inline bool 
MergeGraphImpl<ID_TYPE>::hasNodeId(
    const typename MergeGraphImpl<ID_TYPE>::IdType nodeIndex
)const{
    return nodeUfd_.find(nodeIndex)==nodeIndex;
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::reprEdgeId(
    const typename MergeGraphImpl<ID_TYPE>::IdType edgeIndex
)const{
    return edgeUfd_.find(edgeIndex);
}

template<class ID_TYPE>
inline typename MergeGraphImpl<ID_TYPE>::IdType 
MergeGraphImpl<ID_TYPE>::reprNodeId(
    const typename MergeGraphImpl<ID_TYPE>::IdType nodeIndex
)const{
    return nodeUfd_.find(nodeIndex);
}

template<class ID_TYPE>
inline bool MergeGraphImpl<ID_TYPE>::stateOfInitalEdge(
    const typename MergeGraphImpl<ID_TYPE>::IdType initalEdge
)const{
    const index_type rep = reprEdgeId(initalEdge);
    const index_type rnid0=  reprNodeId(initEdges_[rep].u());
    const index_type rnid1=  reprNodeId(initEdges_[rep].v());
    return rnid0!=rnid1;
}

template<class ID_TYPE>
template<class OUT_ITER>
void MergeGraphImpl<ID_TYPE>::stateOfInitalEdges(
    OUT_ITER begin,
    OUT_ITER end
)const{
    const size_t d = std::distance(begin,end);
    for(size_t ie=0;ie<initNumberOfEdges();++ie){
        const bool state = stateOfInitalEdge(ie);
        begin[ie]=state;
    }
}

template<class ID_TYPE>
inline void MergeGraphImpl<ID_TYPE>::setInitalEdge(
    const size_t initEdge,
    const size_t initNode0,
    const size_t initNode1
){

    // set up inital and dynamic edges
    initEdges_[initEdge][0] =initNode0;
    initEdges_[initEdge][1] =initNode1;
    initEdges_[initEdge][2] =initEdge;

    // set up the edges of a given region mapping
    nodeVector_[initNode0].edges_.insert(initEdge);
    nodeVector_[initNode1].edges_.insert(initEdge);
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::numberOfNodes()const{
    return nodeUfd_.numberOfSets();
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::numberOfEdges()const{
    //return dynamicEdges_.size();
    return edgeUfd_.numberOfSets();
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::nodeNum()const{
    return numberOfNodes();
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::edgeNum()const{
    return numberOfEdges();
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::num_nodes()const{
    return numberOfNodes();
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::num_edges()const{
    return numberOfEdges();
}



template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::initNumberOfNodes()const{
    return nInitNodes_;
}

template<class ID_TYPE>
inline size_t MergeGraphImpl<ID_TYPE>::initNumberOfEdges()const{
    return nInitEdges_;
}

template<class ID_TYPE>
void MergeGraphImpl<ID_TYPE>::mergeParallelEdges(){
    typedef typename DoubleMap::const_iterator MapIter;
    DoubleMap pEdgeFinder;
    for(size_t e=0;e<nInitEdges_;++e){

        IdType n0=initEdges_[e].u();
        IdType n1=initEdges_[e].v();

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
void MergeGraphImpl<ID_TYPE>::combineDoubleEdges(
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
void MergeGraphImpl<ID_TYPE>::searchLocalDoubleEdges(
    const MergeGraphImpl<ID_TYPE>::NodeStorageType & node , 
    MergeGraphImpl<ID_TYPE>::DoubleMap & doubleMap,
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
        //IdType newNodes[2]={
        //    id(u(oldEdge)) ==relabelFrom ? relabelTo : id(u(oldEdge)) , 
        //    id(v(oldEdge)) ==relabelFrom ? relabelTo : id(v(oldEdge)) , 
        //};

        IdType newNodes[2]={
            id(u(oldEdge)),id(v(oldEdge))
        };


        CGP_ASSERT_OP(newNodes[0],!=,newNodes[1]);

        
        if(newNodes[1]<newNodes[0]){
            std::swap(newNodes[1],newNodes[0]);
        }
        const size_t  key = newNodes[0] + newNodes[1]*this->initNumberOfNodes();
        doubleMap[key].push_back(outEdgeIndex);
    }
}

template<class ID_TYPE>
void MergeGraphImpl<ID_TYPE>::mergeRegions(
    const typename MergeGraphImpl<ID_TYPE>::IdType toDeleteEdgeIndex
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



    const size_t  edgeSizeRep    = nodeVector_[newNodeRep].numberOfEdges();
    const size_t  edgeSizeNotRep = nodeVector_[notNewNodeRep].numberOfEdges();

    // the new region wich is the result of the merge
    NodeStorageType & newFormedNode = nodeVector_[newNodeRep];

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


                    const index_type nodeUId = id(u(toMergeEdgeIndex));
                    const index_type nodeVId = id(v(toMergeEdgeIndex));


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
    CGP_ASSERT_OP(this->numberOfNodes(),==,preNumNodes-1);
}




} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_HXX