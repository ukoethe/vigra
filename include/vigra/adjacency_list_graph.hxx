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


    namespace detail_adjacency_list_graph{

        template<class G,class ITEM>
        class ItemIter
        : public boost::iterator_facade<
            ItemIter<G,ITEM>,
            const ITEM,
            boost::forward_traversal_tag
        >
        {

            typedef vigra::GraphItemHelper<G,ITEM> ItemHelper;
            typedef typename G::index_type index_type;

        public:
            ItemIter(const lemon::Invalid & iv = lemon::INVALID)
            :   graph_(NULL),
                id_(-1),
                item_(lemon::INVALID)
            {
            }

            ItemIter(const G & g)
            :   graph_(&g),
                id_(g.zeroStart() ? 0 : 1),
                item_(ItemHelper::itemFromId(*graph_,id_))
            {
                while( !isEnd()  &&  item_==lemon::INVALID ){
                    ++id_;
                    item_ = ItemHelper::itemFromId(*graph_,id_);
                }
            }

            ItemIter(const G & g,const ITEM & item)
            :   graph_(&g),
                id_(g.id(item)),
                item_(item)
            {

            }

        private:

            friend class boost::iterator_core_access;
            bool isEnd( )const{
                return graph_==NULL || id_>ItemHelper::maxItemId(*graph_);
            }
            bool isBegin( )const{
                return graph_!=NULL && (graph_->zeroStart() && id_ == 0 || !graph_->zeroStart() && id_ == 1);
            }

            bool equal(const ItemIter & other) const{
                return   (isEnd() && other.isEnd() ) || (isEnd()==other.isEnd() && (id_ == other.id_) );
            }

            void increment(){
                ++id_;
                item_ = ItemHelper::itemFromId(*graph_,id_);
                while( !isEnd()  &&  item_==lemon::INVALID ){
                    ++id_;
                    item_ = ItemHelper::itemFromId(*graph_,id_);
                }
            }
            const ITEM & dereference()const{
                return item_;
            }
            const G * graph_;
            index_type id_;
            ITEM item_;
        };



        template<class GRAPH>
        class ArcIt
        : public boost::iterator_facade<
            ArcIt<GRAPH>,
            const typename GRAPH::Arc,
            boost::forward_traversal_tag
        >
        {
        public:
            typedef GRAPH Graph;
            typedef typename  Graph::Arc Arc;
            typedef typename  Graph::Edge Edge;
            typedef typename  Graph::EdgeIt EdgeIt;
            ArcIt(const lemon::Invalid invalid = lemon::INVALID )
            :   graph_(NULL),
                pos_(),
                inFirstHalf_(false),
                veryEnd_(true),
                arc_(){
            }
            ArcIt(const GRAPH & g )
            :   graph_(&g),
                pos_(g),
                inFirstHalf_(true),
                veryEnd_( g.edgeNum()==0 ? true : false),
                arc_(){
            }

            ArcIt(const GRAPH & g , const Arc & arc )
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
            bool equal(ArcIt const& other) const{
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

    }


    
    class AdjacencyListGraph
    {
        
    public:
        // public typdedfs
        typedef Int64                                                     index_type;
    private:
        // private typedes which are needed for defining public typedes
        typedef AdjacencyListGraph                                          GraphType;
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


        typedef detail_adjacency_list_graph::ItemIter<GraphType,Edge>    EdgeIt;
        typedef detail_adjacency_list_graph::ItemIter<GraphType,Node>    NodeIt; 
        typedef detail_adjacency_list_graph::ArcIt<GraphType>            ArcIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,NnFilter  >  NeighborNodeIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,IncFilter >  IncEdgeIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,InFlter   >  InArcIt;
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,OutFilter >  OutArcIt;



    // public member functions
    public:
        AdjacencyListGraph(const size_t nodes=0,const size_t edges=0,const bool zeroStart=false);

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
            const index_type id = nodes_.size();
            nodes_.push_back(NodeStorage(id));
            ++nodeNum_;
            return Node(id);
        }

        Edge addEdge(const Node & u , const Node & v){
            const Edge foundEdge  = findEdge(u,v);
            if(foundEdge!=lemon::INVALID){
                return foundEdge;
            }
            else{
                const index_type id = edges_.size();
                edges_.push_back(EdgeStorage(u.id(),v.id(),id));
                nodeImpl(u).insertEdgeId(id);
                nodeImpl(v).insertEdgeId(id);
                ++edgeNum_;
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


        friend class detail_adjacency_list_graph::ItemIter<GraphType,Node>;
        friend class detail_adjacency_list_graph::ItemIter<GraphType,Edge>;


        const NodeStorage & nodeImpl(const Node & node)const{
            return nodes_[node.id()];
        }

        NodeStorage & nodeImpl(const Node & node){
            return nodes_[node.id()];
        }


        bool zeroStart()const{
            return zeroStart_;
        }




        // graph
        NodeVector nodes_;
        EdgeVector edges_;

        size_t nodeNum_;
        size_t edgeNum_;

        bool zeroStart_;
    };



    
    AdjacencyListGraph::AdjacencyListGraph(
        const size_t reserveNodes,
        const size_t reserveEdges,
        const bool   zeroStart
    )
    :   nodes_(),
        edges_(),
        nodeNum_(0),
        edgeNum_(0),
        zeroStart_(zeroStart)
    {
        nodes_.reserve(reserveNodes);
        edges_.reserve(reserveEdges);

        if(!zeroStart_){
            nodes_.push_back(NodeStorage(lemon::INVALID));
            edges_.push_back(EdgeStorage(lemon::INVALID));
        }
    }


    
    inline AdjacencyListGraph::Arc 
    AdjacencyListGraph::direct(
        const AdjacencyListGraph::Edge & edge,
        const bool forward
    )const{
        if(edge!=lemon::INVALID){
            if(forward)
                return Arc(id(edge),id(edge));
            else
                return Arc(id(edge)+maxEdgeId()+1,id(edge));
        }
        else
            return Arc(lemon::INVALID);
    }

    
    inline AdjacencyListGraph::Arc 
    AdjacencyListGraph::direct(
        const AdjacencyListGraph::Edge & edge,
        const AdjacencyListGraph::Node & node
    )const{
        if(u(edge)==node){
            return Arc(id(edge),id(edge));
        }
        else if(v(edge)==node){
            return Arc(id(edge)+maxEdgeId()+1,id(edge));
        }
        else{
            return Arc(lemon::INVALID);
        }
    }

    
    inline bool
    AdjacencyListGraph::direction(
        const AdjacencyListGraph::Arc & arc
    )const{
        return id(arc)<=maxEdgeId();
    }

    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::u(
        const AdjacencyListGraph::Edge & edge
    )const{
        return Node(edges_[id(edge)].u());
    }

    
    inline AdjacencyListGraph::Node
    AdjacencyListGraph::v(
        const AdjacencyListGraph::Edge & edge
    )const{
        return Node(edges_[id(edge)].v());
    }


    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::source(
        const AdjacencyListGraph::Arc & arc
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


    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::target(
        const AdjacencyListGraph::Arc & arc
    )const{
        const index_type arcIndex  = id(arc);
        if (arcIndex > maxEdgeId() ){
            const index_type edgeIndex = arc.edgeId();
            const Edge edge = edgeFromId(edgeIndex);
            return u(edge);
        }
        else{
            const index_type edgeIndex = arcIndex;
            const Edge edge = edgeFromId(edgeIndex);
            return v(edge);
        }
    }

    
    inline AdjacencyListGraph::Node
    AdjacencyListGraph::oppositeNode(
        const AdjacencyListGraph::Node &n,
        const AdjacencyListGraph::Edge &e
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


    
    AdjacencyListGraph::Node 
    AdjacencyListGraph::baseNode(
        const AdjacencyListGraph::IncEdgeIt & iter
    )const{
        return u(*iter);
    }

    
    AdjacencyListGraph::Node 
    AdjacencyListGraph::baseNode(
        const AdjacencyListGraph::OutArcIt & iter 
    )const{
        return source(*iter);
    }


    
    AdjacencyListGraph::Node 
    AdjacencyListGraph::runningNode(
        const AdjacencyListGraph::IncEdgeIt & iter
    )const{
        return v(*iter);
    }

    
    AdjacencyListGraph::Node 
    AdjacencyListGraph::runningNode(
        const AdjacencyListGraph::OutArcIt & iter 
    )const{
        return target(*iter);
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::edgeNum()const{
        return edgeNum_;
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::nodeNum()const{
        return nodeNum_;
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::arcNum()const{
        return edgeNum()*2;
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::maxEdgeId()const{
        return edges_.back().id();
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::maxNodeId()const{
        return nodes_.back().id();
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::maxArcId()const{
        return maxEdgeId()*2+1;
    }

    // ids 
    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::id(
        const AdjacencyListGraph::Node & node
    )const{
        return node.id();
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::id(
        const AdjacencyListGraph::Edge & edge
    )const{
        return edge.id();
    }

    
    inline AdjacencyListGraph::index_type 
    AdjacencyListGraph::id(
        const AdjacencyListGraph::Arc & arc
    )const{
        return arc.id();
    }

    // get edge / node from id
    
    inline AdjacencyListGraph::Edge 
    AdjacencyListGraph::edgeFromId(
        const AdjacencyListGraph::index_type id
    )const{
        if(id<edges_.size())
            return edges_[id].id();
        else
            return Edge(lemon::INVALID);
    }

    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::nodeFromId(
        const AdjacencyListGraph::index_type id
    )const{
        if(id<nodes_.size())
            return nodes_[id].id();
        else
            return Node(lemon::INVALID);
    }

    
    inline AdjacencyListGraph::Arc 
    AdjacencyListGraph::arcFromId(
        const AdjacencyListGraph::index_type id
    )const{
        if(id<=maxEdgeId()){
            const Edge e = edgeFromId(id);
            if(e==lemon::INVALID)
                return Arc(lemon::INVALID);
            else
                return Arc(id,id);
        }
        else{
            const index_type edgeId = id - (maxEdgeId() + 1);
            const Edge e = edgeFromId(edgeId);
            if(e==lemon::INVALID)
                return Arc(lemon::INVALID);
            else
                return Arc(id,edgeId);
        }
    }

    
    inline  AdjacencyListGraph::Edge  
    AdjacencyListGraph::findEdge(
        const AdjacencyListGraph::Node & a,
        const AdjacencyListGraph::Node & b
    )const{
        std::pair<index_type,bool> res =  nodes_[id(a)].sharedEdge(nodes_[id(b)]);
        if(res.second){
            return Edge(res.first);
        }
        else{
            return Edge(lemon::INVALID);
        }
    }


    
    inline  AdjacencyListGraph::Arc  
    AdjacencyListGraph::findArc(
        const AdjacencyListGraph::Node & uNode,
        const AdjacencyListGraph::Node & vNode
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
    
    inline AdjacencyListGraph::NodeIt 
    AdjacencyListGraph::nodesBegin()const{
        return NodeIt(0,nodeNum());
    }

    
    inline AdjacencyListGraph::NodeIt 
    AdjacencyListGraph::nodesEnd()const{  
        return NodeIt(nodeNum(),nodeNum());
    }


    
    inline AdjacencyListGraph::EdgeIt 
    AdjacencyListGraph::edgesBegin()const{
        return EdgeIt(0,edgeNum());
    }

    
    inline AdjacencyListGraph::EdgeIt 
    AdjacencyListGraph::edgesEnd()const{  
        return EdgeIt(edgeNum(),edgeNum());
    }

}

#endif /*VIGRA_RAG_HXX*/