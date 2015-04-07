/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
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


#ifndef VIGRA_ADJACENCY_LIST_GRAPH_HXX
#define VIGRA_ADJACENCY_LIST_GRAPH_HXX

/*std*/
#include <vector>
#include  <set>

/*vigra*/
#include "multi_array.hxx"
#include "multi_gridgraph.hxx"
#include "graphs.hxx"
#include "tinyvector.hxx"
#include "random_access_set.hxx"
#include "graph_maps.hxx"
#include "iteratorfacade.hxx"

#include "algorithm.hxx"
#include "graph_item_impl.hxx"


namespace vigra{

/** \addtogroup GraphDataStructures
*/
//@{

    namespace detail_adjacency_list_graph{

        template<class G,class ITEM>
        class ItemIter
         : public ForwardIteratorFacade<
            ItemIter<G,ITEM>,ITEM,true
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
                id_(0),
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

            friend class vigra::IteratorFacadeCoreAccess;
            bool isEnd( )const{
                return graph_==NULL ||  ItemHelper::itemNum(*graph_)==0 || id_>ItemHelper::maxItemId(*graph_);
            }
            bool isBegin( )const{
                return graph_!=NULL &&  id_ == 0 ;
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
        : public ForwardIteratorFacade<
            ArcIt<GRAPH>,
            typename GRAPH::Arc,true
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

    } // namespace detail_adjacency_list_graph


    /** \brief undirected adjacency list graph in the LEMON API 

    */
    class AdjacencyListGraph
    {
        
    public:
        // public typdedfs
        typedef Int64                                                     index_type;
    private:
        // private typedes which are needed for defining public typedes
        typedef AdjacencyListGraph                                          GraphType;
        typedef detail::GenericNodeImpl<index_type,false>                   NodeStorage;
        typedef detail::GenericEdgeImpl<index_type >                        EdgeStorage;
        typedef detail::NeighborNodeFilter<GraphType>                       NnFilter;
        typedef detail::IncEdgeFilter<GraphType>                            IncFilter;
        typedef detail::IsInFilter<GraphType>                               InFlter;
        typedef detail::IsOutFilter<GraphType>                              OutFilter;
        typedef detail::IsBackOutFilter<GraphType>                          BackOutFilter;
    public:
        // LEMON API TYPEDEFS (and a few more(NeighborNodeIt))

        /// node descriptor
        typedef detail::GenericNode<index_type>                           Node;
        /// edge descriptor
        typedef detail::GenericEdge<index_type>                           Edge;
        /// arc descriptor
        typedef detail::GenericArc<index_type>                            Arc;
        /// edge iterator
        typedef detail_adjacency_list_graph::ItemIter<GraphType,Edge>    EdgeIt;
        /// node iterator
        typedef detail_adjacency_list_graph::ItemIter<GraphType,Node>    NodeIt; 
        /// arc iterator
        typedef detail_adjacency_list_graph::ArcIt<GraphType>            ArcIt;
        
        /// incident edge iterator
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,IncFilter >  IncEdgeIt;
        /// incoming arc iterator
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,InFlter   >  InArcIt;
        /// outgoing arc iterator
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,OutFilter >  OutArcIt;

        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,NnFilter  >  NeighborNodeIt;


        /// outgoing back arc iterator
        typedef detail::GenericIncEdgeIt<GraphType,NodeStorage,BackOutFilter >  OutBackArcIt;


        // BOOST GRAPH API TYPEDEFS
        // - categories (not complete yet)
        typedef boost::directed_tag     directed_category;
        // iterators
        typedef NeighborNodeIt          adjacency_iterator;
        typedef EdgeIt                  edge_iterator;
        typedef NodeIt                  vertex_iterator;
        typedef IncEdgeIt               in_edge_iterator;
        typedef IncEdgeIt               out_edge_iterator;

        // size types
        typedef size_t                  degree_size_type;
        typedef size_t                  edge_size_type;
        typedef size_t                  vertex_size_type;
        // item descriptors
        typedef Edge edge_descriptor;
        typedef Node vertex_descriptor;


        /// default edge map 
        template<class T>
        struct EdgeMap : DenseEdgeReferenceMap<GraphType,T> {
            EdgeMap(): DenseEdgeReferenceMap<GraphType,T>(){
            }
            EdgeMap(const GraphType & g)
            : DenseEdgeReferenceMap<GraphType,T>(g){
            }
            EdgeMap(const GraphType & g,const T & val)
            : DenseEdgeReferenceMap<GraphType,T>(g,val){
            }
        };

        /// default node map 
        template<class T>
        struct NodeMap : DenseNodeReferenceMap<GraphType,T> {
            NodeMap(): DenseNodeReferenceMap<GraphType,T>(){
            }
            NodeMap(const GraphType & g)
            : DenseNodeReferenceMap<GraphType,T>(g){
            }
            NodeMap(const GraphType & g,const T & val)
            : DenseNodeReferenceMap<GraphType,T>(g,val){
            }
        };

        /// default arc map 
        template<class T>
        struct ArcMap : DenseArcReferenceMap<GraphType,T> {
            ArcMap(): DenseArcReferenceMap<GraphType,T>(){
            }
            ArcMap(const GraphType & g)
            : DenseArcReferenceMap<GraphType,T>(g){
            }
            ArcMap(const GraphType & g,const T & val)
            : DenseArcReferenceMap<GraphType,T>(g,val){
            }
        };



    // public member functions
    public:
        /** \brief Constructor.
         
            @param nodes : reserve space for so many nodes
            @param edges : reserve space for so many edges
        */
        AdjacencyListGraph(const size_t nodes=0,const size_t edges=0);

        /** \brief Get the number of edges in this graph (API: LEMON).
        */
        index_type edgeNum()const;

        /** \brief Get the number of nodes in this graph (API: LEMON).
        */
        index_type nodeNum()const;
        /** \brief Get the number of arcs in this graph (API: LEMON).
        */
        index_type arcNum()const;

        /** \brief Get the maximum ID of any edge in this graph (API: LEMON).
        */
        index_type maxEdgeId()const;
        /** \brief Get the maximum ID of any node in this graph (API: LEMON).
        */
        index_type maxNodeId()const;
        /** \brief Get the maximum ID of any edge in arc graph (API: LEMON).
        */
        index_type maxArcId()const;

        /** \brief Create an arc for the given edge \a e, oriented along the 
            edge's natural (<tt>forward = true</tt>) or reversed 
            (<tt>forward = false</tt>) direction (API: LEMON).
        */
        Arc direct(const Edge & edge,const bool forward)const;

        /** \brief Create an arc for the given edge \a e oriented
            so that node \a n is the starting node of the arc (API: LEMON), or
            return <tt>lemon::INVALID</tt> if the edge is not incident to this node.
        */
        Arc direct(const Edge & edge,const Node & node)const;

        /** \brief Return <tt>true</tt> when the arc is looking on the underlying
            edge in its natural (i.e. forward) direction, <tt>false</tt> otherwise (API: LEMON).
        */
        bool direction(const Arc & arc)const;
        /** \brief Get the start node of the given edge \a e (API: LEMON,<br/>
            the boost::graph API provides the free function <tt>boost::source(e, graph)</tt>).
        */
        Node u(const Edge & edge)const;
        /** \brief Get the end node of the given edge \a e (API: LEMON,<br/>
            the boost::graph API provides the free function <tt>boost::target(e, graph)</tt>).
        */
        Node v(const Edge & edge)const;
        /** \brief Get the start node of the given arc \a a (API: LEMON).
        */
        Node source(const Arc & arc)const;
        /** \brief Get the end node of the given arc \a a (API: LEMON).
        */
        Node target(const Arc & arc)const;
        /** \brief Return the opposite node of the given node \a n
            along edge \a e (API: LEMON), or return <tt>lemon::INVALID</tt>
            if the edge is not incident to this node.
        */
        Node oppositeNode(Node const &n, const Edge &e) const;

        /** \brief Return the start node of the edge the given iterator is referring to (API: LEMON).
        */
        Node baseNode(const IncEdgeIt & iter)const;
        /** \brief Return the start node of the edge the given iterator is referring to (API: LEMON).
        */
        Node baseNode(const OutArcIt & iter)const;

        /** \brief Return the end node of the edge the given iterator is referring to (API: LEMON).
        */
        Node runningNode(const IncEdgeIt & iter)const;
        /** \brief Return the end node of the edge the given iterator is referring to (API: LEMON).
        */ 
        Node runningNode(const OutArcIt & iter)const;


        /** \brief Get the ID  for node desciptor \a v (API: LEMON).
        */
        index_type id(const Node & node)const;
        /** \brief Get the ID  for edge desciptor \a v (API: LEMON).
        */
        index_type id(const Edge & edge)const;
        /** \brief Get the ID  for arc desciptor \a v (API: LEMON).
        */
        index_type id(const Arc  & arc )const;

        /** \brief Get edge descriptor for given node ID \a i (API: LEMON).
            Return <tt>Edge(lemon::INVALID)</tt> when the ID does not exist in this graph.
        */
        Edge edgeFromId(const index_type id)const;

        /** \brief Get node descriptor for given node ID \a i (API: LEMON).
            Return <tt>Node(lemon::INVALID)</tt> when the ID does not exist in this graph.
        */
        Node nodeFromId(const index_type id)const;
        /** \brief Get arc descriptor for given node ID \a i (API: LEMON).
            Return <tt>Arc(lemon::INVALID)</tt> when the ID does not exist in this graph.
        */
        Arc  arcFromId(const index_type id)const;


        /** \brief Get a descriptor for the edge connecting vertices \a u and \a v,<br/>or <tt>lemon::INVALID</tt> if no such edge exists (API: LEMON).
        */
        Edge findEdge(const Node & a,const Node & b)const;
        /** \brief Get a descriptor for the arc connecting vertices \a u and \a v,<br/>or <tt>lemon::INVALID</tt> if no such edge exists (API: LEMON).
        */
        Arc  findArc(const Node & u,const Node & v)const;

        /* \brief add a new node to the graph. 
            the next unused id will be assigned to the node
        */
        Node addNode();
        /* \brief add a  node to the graph with a given id.
            If there is  another node with this id, no
            new node will be added.
        */
        Node addNode(const index_type id);


        /*  \brief this will remove any nodes if there are existing nodes (and edges)
            and will add nodes in the range of ids , endId is not included!
        */
        void assignNodeRange(const index_type beginId, const index_type endId);



        /* \brief add an edge to the graph.
            If there is an other edge between u and v no new edge will be added.
        */
        Edge addEdge(const Node & u , const Node & v);
        /* \brief add an edge to the graph.
            If there is an other edge between u and v no new edge will be added.
            If the nodes for the given id's are not in the graph, they will be added.
        */
        Edge addEdge(const index_type u ,const index_type v);

        
        size_t maxDegree()const{
            size_t md=0;
            for(NodeIt it(*this);it!=lemon::INVALID;++it){
                std::max(md, size_t( degree(*it) ) );
            }
            return md;
        }


        ////////////////////////
        // BOOST API
        /////////////////////////
        // - sizes 
        // - iterators
        vertex_iterator  get_vertex_iterator()const;
        vertex_iterator  get_vertex_end_iterator()const  ;
        edge_iterator    get_edge_iterator()const;
        edge_iterator    get_edge_end_iterator()const  ;
        degree_size_type degree(const vertex_descriptor & node)const{
            return nodeImpl(node).numberOfEdges();
        }

        static const bool is_directed = false;

    public:

        void reserveMaxNodeId(const index_type mxid ){
            if(nodeNum()==0 ||  mxid>maxNodeId())
                nodes_.reserve(mxid+1);
        }

        void reserveEdges(const size_t size ){
            if(size>nodeNum())
                edges_.reserve(size);
        }


        void clear(){
            nodeNum_=0;
            edgeNum_=0;
            edges_.clear();
            nodes_.clear();
        }
        size_t serializationSize()const{

            // num edges + num nodes 
            // max edge id  + max node id
            size_t size=4;

            // edge ids
            size+= 2*edgeNum();


            for(NodeIt iter(*this); iter!= lemon::INVALID ; ++iter){
                size+= 2+this->degree(*iter)*2;    
            }

            return size;
        }

        template<class ITER>
        void serialize(ITER outIter) const {

            // sizes of graph
            *outIter = nodeNum(); ++outIter;
            *outIter = edgeNum(); ++outIter;
            *outIter = maxNodeId(); ++outIter;
            *outIter = maxEdgeId(); ++outIter;

            // edges
            for(EdgeIt iter(*this); iter!=lemon::INVALID; ++iter){
                const Edge e(*iter);
                const size_t ui = this->id(this->u(e));
                const size_t vi = this->id(this->v(e));
                *outIter = ui; ++outIter;
                *outIter = vi; ++outIter;
            }



            // node neighbors
            for(NodeIt iter(*this); iter!= lemon::INVALID ; ++iter){
                const Node n(*iter);

                *outIter = this->id(*iter); ++outIter;
                *outIter = this->degree(*iter); ++outIter;

                for(OutArcIt eIter(*this,n); eIter!=lemon::INVALID; ++eIter){
                    const Edge e(*eIter);
                    const Node oNode(this->target(*eIter));

                    const size_t ei = this->id(e);
                    const size_t oni = this->id(oNode);

                    *outIter = ei; ++outIter;
                    *outIter = oni; ++outIter;
                }
            }

        }

        template<class ITER>
        void deserialize(ITER begin, ITER end){


            nodeNum_ = *begin; ++begin;
            edgeNum_ = *begin; ++begin;
            const size_t maxNid = *begin; ++begin;
            const size_t maxEid = *begin; ++begin;

            nodes_.clear();
            edges_.clear();
            nodes_.resize(maxNid+1, NodeStorage());
            edges_.resize(maxEid+1, EdgeStorage());

            // set up edges
            for(size_t eid=0; eid<edgeNum_; ++eid){
                const size_t u = *begin; ++begin;
                const size_t v = *begin; ++begin;
                nodes_[u].setId(u);
                nodes_[v].setId(v);
                edges_[eid]=EdgeStorage(u,v,eid);
            }

            // set up nodes
            for(size_t i=0; i<nodeNum_; ++i){

                const size_t id = *begin; ++begin;
                const size_t nodeDegree=*begin; ++begin;

                NodeStorage & nodeImpl = nodes_[id];
                nodeImpl.setId(id);
                for(size_t d=0; d<nodeDegree; ++d){
                    const size_t ei  = *begin; ++begin;
                    const size_t oni =  *begin; ++begin;
                    nodeImpl.insert(oni, ei);
                }
            }
        }

    private:
        // private typedefs
        typedef std::vector<NodeStorage> NodeVector;
        typedef std::vector<EdgeStorage> EdgeVector;


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
        friend struct detail::IsBackOutFilter;
        template<class G>
        friend struct detail::IsInFilter;


        friend class detail_adjacency_list_graph::ItemIter<GraphType,Node>;
        friend class detail_adjacency_list_graph::ItemIter<GraphType,Edge>;


        const NodeStorage & nodeImpl(const Node & node)const{
            return nodes_[node.id()];
        }

        NodeStorage & nodeImpl(const Node & node){
            return nodes_[node.id()];
        }





        // graph
        NodeVector nodes_;
        EdgeVector edges_;

        size_t nodeNum_;
        size_t edgeNum_;
    };



#ifndef DOXYGEN  // doxygen doesn't like out-of-line definitions

    inline AdjacencyListGraph::AdjacencyListGraph(
        const size_t reserveNodes,
        const size_t reserveEdges
    )
    :   nodes_(),
        edges_(),
        nodeNum_(0),
        edgeNum_(0)
    {
        nodes_.reserve(reserveNodes);
        edges_.reserve(reserveEdges);
    }


    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::addNode(){
        const index_type id = nodes_.size();
        nodes_.push_back(NodeStorage(id));
        ++nodeNum_;
        return Node(id);
    }

    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::addNode(const AdjacencyListGraph::index_type id){
        if(id == nodes_.size()){
            nodes_.push_back(NodeStorage(id));
            ++nodeNum_;
            return Node(id);
        }
        else if((std::size_t)id < nodes_.size()){
            const Node node = nodeFromId(id);
            if(node==lemon::INVALID){
                nodes_[id]=NodeStorage(id);
                ++nodeNum_;
                return Node(id);
            }
            else{
                return node;
            }
        }
        else{
            // refactor me
            while(nodes_.size() < (std::size_t)id){
                nodes_.push_back(NodeStorage(lemon::INVALID));
            }
            nodes_.push_back(NodeStorage(id));
            ++nodeNum_;
            return Node(id);
        }
    }


    inline void 
    AdjacencyListGraph::assignNodeRange(const AdjacencyListGraph::index_type beginId, const AdjacencyListGraph::index_type endId){
        nodes_.clear();
        edges_.clear();
        edgeNum_=0;
        nodeNum_ = endId - beginId;
        nodes_.resize(endId);
        for(index_type i=beginId; i<endId; ++i)
            nodes_[i]=NodeStorage(i);
    }



    inline AdjacencyListGraph::Edge 
    AdjacencyListGraph::addEdge(
        const AdjacencyListGraph::Node & u , 
        const AdjacencyListGraph::Node & v
    ){
        const Edge foundEdge  = findEdge(u,v);
        if(foundEdge!=lemon::INVALID){
            return foundEdge;
        }
        else if(u==lemon::INVALID || v==lemon::INVALID){
            return Edge(lemon::INVALID);
        }
        else{
            const index_type eid  = edges_.size();
            const index_type uid = u.id();
            const index_type vid = v.id();
            edges_.push_back(EdgeStorage(uid,vid,eid));
            nodeImpl(u).insert(vid,eid);
            nodeImpl(v).insert(uid,eid);
            ++edgeNum_;
            return Edge(eid);
        }   
    }

    inline AdjacencyListGraph::Edge 
    AdjacencyListGraph::addEdge(
        const AdjacencyListGraph::index_type u ,
        const AdjacencyListGraph::index_type v
    ){
        const Node uu = addNode(u);
        const Node vv = addNode(v);
        return addEdge(uu,vv);
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


    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::baseNode(
        const AdjacencyListGraph::IncEdgeIt & iter
    )const{
        return u(*iter);
    }

    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::baseNode(
        const AdjacencyListGraph::OutArcIt & iter 
    )const{
        return source(*iter);
    }


    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::runningNode(
        const AdjacencyListGraph::IncEdgeIt & iter
    )const{
        return v(*iter);
    }

    
    inline AdjacencyListGraph::Node 
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
        if((std::size_t)id < edges_.size() && edges_[id].id() != -1)
            return Edge(edges_[id].id());
        else
            return Edge(lemon::INVALID);
    }

    
    inline AdjacencyListGraph::Node 
    AdjacencyListGraph::nodeFromId(
        const AdjacencyListGraph::index_type id
    )const{
        if((std::size_t)id < nodes_.size() && nodes_[id].id() != -1)
            return Node(nodes_[id].id());
        else
            return Node(lemon::INVALID);
    }

    
    inline AdjacencyListGraph::Arc 
    AdjacencyListGraph::arcFromId(
        const AdjacencyListGraph::index_type id
    )const{
        if(id<=maxEdgeId()){
            if(edgeFromId(id)==lemon::INVALID)
                return Arc(lemon::INVALID);
            else
                return Arc(id,id);
        }
        else{
            const index_type edgeId = id - (maxEdgeId() + 1);
            if( edgeFromId(edgeId)==lemon::INVALID)
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
        if(a!=b){
            std::pair<index_type,bool> res =  nodes_[id(a)].findEdge(id(b));
            if(res.second){
                return Edge(res.first);
            }
        }
        return Edge(lemon::INVALID);
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
    
    inline AdjacencyListGraph::vertex_iterator 
    AdjacencyListGraph::get_vertex_iterator()const{
        return NodeIt(0,nodeNum());
    }

    
    inline AdjacencyListGraph::vertex_iterator 
    AdjacencyListGraph::get_vertex_end_iterator()const{  
        return NodeIt(nodeNum(),nodeNum());
    }


    
    inline AdjacencyListGraph::edge_iterator 
    AdjacencyListGraph::get_edge_iterator()const{
        return EdgeIt(0,edgeNum());
    }

    
    inline AdjacencyListGraph::edge_iterator 
    AdjacencyListGraph::get_edge_end_iterator()const{  
        return EdgeIt(edgeNum(),edgeNum());
    }

#endif //DOXYGEN

//@}

} // namespace vigra


// boost free functions specialized for adjacency list graph
namespace boost{

    ////////////////////////////////////
    // functions to get size of the graph
    ////////////////////////////////////
    inline vigra::AdjacencyListGraph::vertex_size_type
    num_vertices(const vigra::AdjacencyListGraph & g){
        return g.nodeNum();
    }
    inline vigra::AdjacencyListGraph::edge_size_type
    num_edges(const vigra::AdjacencyListGraph & g){
        return g.edgeNum();
    }


    ////////////////////////////////////
    // functions to get degrees of nodes
    // (degree / indegree / outdegree)
    ////////////////////////////////////
    inline vigra::AdjacencyListGraph::degree_size_type
    degree(const vigra::AdjacencyListGraph::vertex_descriptor & v , const vigra::AdjacencyListGraph & g){
        return g.degree(v);
    }
    // ??? check if this is the right impl. for undirected graphs
    inline vigra::AdjacencyListGraph::degree_size_type
    in_degree(const vigra::AdjacencyListGraph::vertex_descriptor & v , const vigra::AdjacencyListGraph & g){
        return g.degree(v);
    }
    // ??? check if this is the right impl. for undirected graphs
    inline vigra::AdjacencyListGraph::degree_size_type
    out_degree(const vigra::AdjacencyListGraph::vertex_descriptor & v , const vigra::AdjacencyListGraph & g){
        return g.degree(v);
    }


    ////////////////////////////////////
    // functions to u/v source/target
    ////////////////////////////////////
    inline vigra::AdjacencyListGraph::vertex_descriptor
    source(const vigra::AdjacencyListGraph::edge_descriptor & e , const vigra::AdjacencyListGraph & g){
        return g.u(e);
    }

    inline vigra::AdjacencyListGraph::vertex_descriptor
    target(const vigra::AdjacencyListGraph::edge_descriptor & e , const vigra::AdjacencyListGraph & g){
        return g.v(e);
    }

    ////////////////////////////////////
    // functions to get iterator pairs
    ////////////////////////////////////
    inline  std::pair< vigra::AdjacencyListGraph::vertex_iterator, vigra::AdjacencyListGraph::vertex_iterator >
    vertices(const vigra::AdjacencyListGraph & g ){
        return std::pair< vigra::AdjacencyListGraph::vertex_iterator, vigra::AdjacencyListGraph::vertex_iterator >(
            g.get_vertex_iterator(), g.get_vertex_end_iterator());
    }
    inline  std::pair< vigra::AdjacencyListGraph::edge_iterator, vigra::AdjacencyListGraph::edge_iterator >
    edges(const vigra::AdjacencyListGraph & g ){
        return std::pair< vigra::AdjacencyListGraph::edge_iterator, vigra::AdjacencyListGraph::edge_iterator >(
            g.get_edge_iterator(),g.get_edge_end_iterator());
    }


    inline  std::pair< vigra::AdjacencyListGraph::in_edge_iterator, vigra::AdjacencyListGraph::in_edge_iterator >
    in_edges(const vigra::AdjacencyListGraph::vertex_descriptor & v,  const vigra::AdjacencyListGraph & g ){
        return std::pair< vigra::AdjacencyListGraph::in_edge_iterator, vigra::AdjacencyListGraph::in_edge_iterator >(
            vigra::AdjacencyListGraph::in_edge_iterator(g,v),vigra::AdjacencyListGraph::in_edge_iterator(lemon::INVALID)
        );
    }
    inline  std::pair< vigra::AdjacencyListGraph::out_edge_iterator, vigra::AdjacencyListGraph::out_edge_iterator >
    out_edges(const vigra::AdjacencyListGraph::vertex_descriptor & v,  const vigra::AdjacencyListGraph & g ){
        return std::pair< vigra::AdjacencyListGraph::out_edge_iterator, vigra::AdjacencyListGraph::out_edge_iterator >(
            vigra::AdjacencyListGraph::out_edge_iterator(g,v),vigra::AdjacencyListGraph::out_edge_iterator(lemon::INVALID)
        );
    }

}  // namespace boost

#endif /*VIGRA_ADJACENCY_LIST_GRAPH_HXX*/
