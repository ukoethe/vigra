#ifndef VIGRA_RAG_HXX
#define VIGRA_RAG_HXX



// WHY RAG SHOULD TAKE VIEW AND NOT(!) GRAPH WITH LABLED NODE MAP
// - WHO SHOULD EXTRACT EDGE FEATURES 
//    - if rag TAKES graph then we do not know anything about edge coordinates ?
//      we could use a edgeMap or nodeMap of features, but who is going to
//      extract the inital features....anyone who should do does
//      need to know about the RAG and the labeling and also
//      which edge correspondes to which pixels.....and well..this
//      whould need another data structure
//      



/*std*/
#include <vector>
#include  <set>

/*boost*/
#include <boost/foreach.hpp>

/*vigra*/
#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/graphs.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/random_access_set.hxx>
#include <vigra/graph_helper/dense_map.hxx>


#include <vigra/algorithm.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>

#include <vigra/rag/rag_iterators.hxx>

namespace vigra{





    template<unsigned int OFFSET>
    class AdjacencyListGraph
    {
        
    public:
        // public typdedfs
        typedef Int64                                                     index_type;
    private:
        // private typedes which are needed for defining public typedes
        typedef AdjacencyListGraph<OFFSET>                     GraphType;
        typedef RandomAccessSet<index_type>                               NodeStorageEdgeSet;
        typedef detail::GenericNodeImpl<index_type,NodeStorageEdgeSet>    NodeStorage;
        typedef detail::GenericEdgeImpl<index_type >                      EdgeStorage;
        typedef detail::NeighborNodeFilter<GraphType>                       NnFilter;
        typedef detail::IncEdgeFilter<GraphType>                            IncFilter;
        typedef detail::IsInFilter<GraphType>                               InFlter;
        typedef detail::IsOutFilter<GraphType>                              OutFilter;
    public:
        // LEMON API TYPEDEFS (and a few more(NeighborNodeIt))
        typedef detail::GenericNode<index_type>                           Node;
        typedef detail::GenericEdge<index_type>                           Edge;
        typedef detail::GenericArc<index_type>                            Arc;
        typedef detail::RagItemIt<GraphType,Edge>                           EdgeIt;
        typedef detail::RagItemIt<GraphType,Node>                           NodeIt; 
        typedef detail::RagItemIt<GraphType,Arc>                            ArcIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,NnFilter  >  NeighborNodeIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,IncFilter >  IncEdgeIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,InFlter   >  InArcIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,OutFilter >  OutArcIt;



    // public member functions
    public:
        AdjacencyListGraph(const size_t nodes=0,const size_t edges=0);

        index_type edgeNum()const;
        index_type nodeNum()const;
        index_type arcNum()const;

        index_type maxEdgeId()const;
        index_type maxNodeId()const;
        index_type maxArcId()const;

        Arc direct(const Edge & edge,const bool forward)const;
        Arc direct(const Edge & edge,const Node & node)const;
        bool direction(const Arc & arc)const;

        Node u(const Edge & edge)const;
        Node v(const Edge & edge)const;
        Node source(const Arc & arc)const;
        Node target(const Arc & arc)const;
        Node oppositeNode(Node const &n, const Edge &e) const;

        Node baseNode(const IncEdgeIt & iter)const;
        Node baseNode(const OutArcIt & iter)const;

        Node runningNode(const IncEdgeIt & iter)const;
        Node runningNode(const OutArcIt & iter)const;


        // ids 
        index_type id(const Node & node)const;
        index_type id(const Edge & edge)const;
        index_type id(const Arc  & arc )const;

        // get edge / node from id
        Edge edgeFromId(const index_type id)const;
        Node nodeFromId(const index_type id)const;
        Arc  arcFromId(const index_type id)const;


        // find edge
        Edge findEdge(const Node & a,const Node & b)const;
        Arc  findArc(const Node & u,const Node & v)const;


        // iterators
        NodeIt nodesBegin()const;
        NodeIt nodesEnd()const  ;
        EdgeIt edgesBegin()const;
        EdgeIt edgesEnd()const  ;


        Node addNode(){
            const index_type id = nodes_.size()+OFFSET;
            nodes_.push_back(NodeStorage(id));
            return Node(id);
        }

        Edge addEdge(const Node & u , const Node & v){
            const Edge foundEdge  = findEdge(u,v);
            if(foundEdge!=lemon::INVALID){
                return foundEdge;
            }
            else{
                const index_type id = edges_.size()+OFFSET;
                edges_.push_back(EdgeStorage(u.id(),v.id(),id));
                nodeImpl(u).insertEdgeId(id);
                nodeImpl(v).insertEdgeId(id);
                return Edge(id);
            }   
        }

    private:
        // private typedefs
        typedef std::vector<NodeStorage> NodeVector;
        typedef std::vector<EdgeStorage> EdgeVector;


        // needs acces to const nodeImpl
        template<class G,class NIMPL,class FILT>
        friend class detail::GenericIncEdgeIt;


        const NodeStorage & nodeImpl(const Node & node)const{
            return nodes_[id(node)-OFFSET];
        }

        NodeStorage & nodeImpl(const Node & node){
            return nodes_[id(node)-OFFSET];
        }


        // graph
        NodeVector nodes_;
        EdgeVector edges_;

    };



    template<unsigned int OFFSET>
    AdjacencyListGraph<OFFSET>::AdjacencyListGraph(
        size_t reserveNodes,
        size_t reserveEdges
    )
    :   nodes_(),
        edges_(){
        nodes_.reserve(reserveNodes);
        edges_.reserve(reserveEdges);
    }


    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Arc 
    AdjacencyListGraph<OFFSET>::direct(
        const typename AdjacencyListGraph<OFFSET>::Edge & edge,
        const bool forward
    )const{
        if(edge!=lemon::INVALID){
            if(forward)
                return Arc(id(edge),id(edge));
            else
                return Arc(id(edge)+maxEdgeId(),id(edge));
        }
        else
            return Arc(lemon::INVALID);
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Arc 
    AdjacencyListGraph<OFFSET>::direct(
        const typename AdjacencyListGraph<OFFSET>::Edge & edge,
        const typename AdjacencyListGraph<OFFSET>::Node & node
    )const{
        if(u(edge)==node){
            return Arc(id(edge),id(edge));
        }
        else if(v(edge)==node){
            return Arc(id(edge)+maxEdgeId(),id(edge));
        }
        else{
            return Arc(lemon::INVALID);
        }
    }

    template<unsigned int OFFSET>
    inline bool
    AdjacencyListGraph<OFFSET>::direction(
        const typename AdjacencyListGraph<OFFSET>::Arc & arc
    )const{
        return id(arc)<=maxEdgeId();
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::u(
        const typename AdjacencyListGraph<OFFSET>::Edge & edge
    )const{
        return Node(edges_[id(edge)-1].u());
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node
    AdjacencyListGraph<OFFSET>::v(
        const typename AdjacencyListGraph<OFFSET>::Edge & edge
    )const{
        return Node(edges_[id(edge)-1].v());
    }


    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::source(
        const typename AdjacencyListGraph<OFFSET>::Arc & arc
    )const{
        const index_type arcIndex  = id(arc);
        if (arcIndex > maxEdgeId() ){
            const index_type edgeIndex = arc.edgeId();
            const Edge edge = edgeFromId(edgeIndex);
            return v(edge);
        }
        else{
            const index_type edgeIndex = arcIndex;
            const Edge edge = edgeFromId(edgeIndex);
            return u(edge);
        }
    }   


    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::target(
        const typename AdjacencyListGraph<OFFSET>::Arc & arc
    )const{
        const index_type arcIndex  = id(arc);
        if (arcIndex > maxEdgeId() ){
            const index_type edgeIndex = arcIndex-maxEdgeId();
            const Edge edge = edgeFromId(edgeIndex);
            return u(edge);
        }
        else{
            const index_type edgeIndex = arcIndex;
            const Edge edge = edgeFromId(edgeIndex);
            return v(edge);
        }
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node
    AdjacencyListGraph<OFFSET>::oppositeNode(
        const typename AdjacencyListGraph<OFFSET>::Node &n,
        const typename AdjacencyListGraph<OFFSET>::Edge &e
    ) const {
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


    template<unsigned int OFFSET>
    typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::baseNode(
        const typename AdjacencyListGraph<OFFSET>::IncEdgeIt & iter
    )const{
        return u(*iter);
    }

    template<unsigned int OFFSET>
    typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::baseNode(
        const typename AdjacencyListGraph<OFFSET>::OutArcIt & iter 
    )const{
        return source(*iter);
    }


    template<unsigned int OFFSET>
    typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::runningNode(
        const typename AdjacencyListGraph<OFFSET>::IncEdgeIt & iter
    )const{
        return v(*iter);
    }

    template<unsigned int OFFSET>
    typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::runningNode(
        const typename AdjacencyListGraph<OFFSET>::OutArcIt & iter 
    )const{
        return target(*iter);
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::edgeNum()const{
        return edges_.size();
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::nodeNum()const{
        return nodes_.size();
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::arcNum()const{
        return edgeNum()*2;
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::maxEdgeId()const{
        return edgeNum()-1+OFFSET;
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::maxNodeId()const{
        return nodeNum()-1+OFFSET;
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::maxArcId()const{
        return maxEdgeId()*2;
    }

    // ids 
    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::id(
        const typename AdjacencyListGraph<OFFSET>::Node & node
    )const{
        return node.id();
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::id(
        const typename AdjacencyListGraph<OFFSET>::Edge & edge
    )const{
        return edge.id();
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::index_type 
    AdjacencyListGraph<OFFSET>::id(
        const typename AdjacencyListGraph<OFFSET>::Arc & arc
    )const{
        return arc.id();
    }

    // get edge / node from id
    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Edge 
    AdjacencyListGraph<OFFSET>::edgeFromId(
        const typename AdjacencyListGraph<OFFSET>::index_type id
    )const{
        if(id>=OFFSET && id<=maxEdgeId())
            return Edge(id);
        else{
            return Edge(lemon::INVALID);
        }
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Node 
    AdjacencyListGraph<OFFSET>::nodeFromId(
        const typename AdjacencyListGraph<OFFSET>::index_type id
    )const{
        if(id>=OFFSET && id<=maxNodeId())
            return Node(id);
        else{
            return Node(lemon::INVALID);
        }
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::Arc 
    AdjacencyListGraph<OFFSET>::arcFromId(
        const typename AdjacencyListGraph<OFFSET>::index_type id
    )const{
        if(id>=OFFSET && id<=maxArcId()){
            if(id<=maxEdgeId()){
                return Arc(id,id);
            }
            else{
                return Arc(id,id-maxEdgeId());
            }
        }
        else{
            return Arc(lemon::INVALID);
        }
    }

    template<unsigned int OFFSET>
    inline  typename AdjacencyListGraph<OFFSET>::Edge  
    AdjacencyListGraph<OFFSET>::findEdge(
        const typename AdjacencyListGraph<OFFSET>::Node & a,
        const typename AdjacencyListGraph<OFFSET>::Node & b
    )const{
        std::pair<index_type,bool> res =  nodes_[id(a)-1].sharedEdge(nodes_[id(b)-1]);
        if(res.second){
            return Edge(res.first);
        }
        else{
            return Edge(lemon::INVALID);
        }
    }


    template<unsigned int OFFSET>
    inline  typename AdjacencyListGraph<OFFSET>::Arc  
    AdjacencyListGraph<OFFSET>::findArc(
        const typename AdjacencyListGraph<OFFSET>::Node & uNode,
        const typename AdjacencyListGraph<OFFSET>::Node & vNode
    )const{
        const Edge e = findEdge(uNode,vNode);
        if(e==lemon::INVALID){
            return Arc(lemon::INVALID);
        }
        else{
            if(u(e)==uNode)
                return direct(e,true) ;
            else
                return direct(e,false) ;
        }
    }


    // iterators
    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::NodeIt 
    AdjacencyListGraph<OFFSET>::nodesBegin()const{
        return NodeIt(0,nodeNum());
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::NodeIt 
    AdjacencyListGraph<OFFSET>::nodesEnd()const{  
        return NodeIt(nodeNum(),nodeNum());
    }


    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::EdgeIt 
    AdjacencyListGraph<OFFSET>::edgesBegin()const{
        return EdgeIt(0,edgeNum());
    }

    template<unsigned int OFFSET>
    inline typename AdjacencyListGraph<OFFSET>::EdgeIt 
    AdjacencyListGraph<OFFSET>::edgesEnd()const{  
        return EdgeIt(edgeNum(),edgeNum());
    }

}

#endif /*VIGRA_RAG_HXX*/