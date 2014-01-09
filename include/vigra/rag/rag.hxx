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

/*vigra*/
#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/graphs.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/random_access_set.hxx>
#include <vigra/is_end_iterator.hxx>
#include <vigra/is_end_transform_iterator.hxx>
#include <vigra/is_end_filter_iterator.hxx>
#include <vigra/is_end_enumeration_iterator.hxx>
#include <vigra/graph_iterator_functors.hxx>
#include <vigra/algorithm.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>
#include <vigra/graph_helper/graph_crtp_base.hxx>

namespace vigra{

    template<class GRAPH,class ITEM>
    struct RagItemIt 
    :   public EnumerationIterator< ITEM > {
    private:
        typedef EnumerationIterator< ITEM > BaseIterType;
    public:
        typedef GRAPH Graph;
        typedef typename Graph::index_type index_type;
        // make this private
        RagItemIt(const index_type pos,const index_type size)
        :   BaseIterType(pos,size){
        }
        // default constructor
        RagItemIt()
        :   BaseIterType(0,0,1){
        }
        RagItemIt(const RagItemIt & other)
        :   BaseIterType(other){
        }

        // Invalid constructor & conversion. 
        RagItemIt(const lemon::Invalid & invalid)
        :   BaseIterType(0,0,1){
        }
        RagItemIt(const Graph & g)
        :   BaseIterType(0, GraphItemHelper<GRAPH,ITEM>::itemNum(g) ,1){
        }
        RagItemIt(const Graph & g,const ITEM & item)
        :   BaseIterType(g.id(item)-1, GraphItemHelper<GRAPH,ITEM>::itemNum(g) , 1 ){
        }
    };

    // specialization for arc
    template<class GRAPH>
    class RagItemIt<GRAPH,typename GRAPH::Arc>
    : public boost::iterator_facade<
        RagItemIt<GRAPH,typename GRAPH::Arc>,
        const typename GRAPH::Arc,
        //boost::forward_traversal_tag 
        boost::forward_traversal_tag
    >
    {
    public:
        typedef GRAPH Graph;
        typedef typename Graph::index_type index_type;
        typedef typename Graph::Arc   Arc;

        RagItemIt(const lemon::Invalid &  invalid = lemon::INVALID)
        :   current_(0),
            maxEdgeId_(0),
            maxArcId_(0){
        }

        RagItemIt(const Graph & g)
        :   current_(1),
            maxEdgeId_(g.maxEdgeId()),
            maxArcId_(g.maxArcId()){
        }
        RagItemIt(const Graph & g,const Arc & arc)
        :   current_(g.id(arc)),
            maxEdgeId_(g.maxEdgeId()),
            maxArcId_(g.maxArcId()){
        }
    private:

        bool isEnd()const{
            return  (maxArcId_==0 && current_==0 )||  current_>maxArcId_;
        }

        bool isBegin()const{
            return current_ == 0 && maxArcId_!=0;
        }



        friend class boost::iterator_core_access;

        void increment() {
            ++current_;
        }
        

        
        bool equal(const RagItemIt & other) const{
            return   (isEnd() && other.isEnd() ) || (current_ == other.current_);
        }

        const Arc & dereference() const { 
            if(current_<=maxEdgeId_){
                arc_ =  Arc(current_,current_);
            }
            else {
                arc_ =  Arc(current_,current_-maxEdgeId_);
            }
            return arc_;
        }




        ptrdiff_t current_;
        ptrdiff_t maxEdgeId_;
        ptrdiff_t maxArcId_;
        mutable Arc arc_;
    };





    template<unsigned int DIM , class IN_LABEL_TYPE>
    class Rag
    {
    private:
        typedef Rag<DIM,IN_LABEL_TYPE> RagType;
        
    public:
        typedef IN_LABEL_TYPE InLabelType;
        typedef MultiArrayView<DIM,InLabelType> InputLabelingView;
        typedef MultiArray<DIM,InLabelType>     InputLabelingArray;

        // graph item typedefs
        typedef Int64 index_type;
        //typedef RagEdgeStorage<RagType> EdgeStorageType;


        typedef detail::GenericNodeImpl<index_type,RandomAccessSet<index_type> > NodeStorageType;
        typedef detail::GenericEdgeImpl<index_type >                             EdgeStorageType;

        typedef detail::GenericNode<index_type>  Node;
        typedef detail::GenericEdge<index_type>  Edge;
        typedef detail::GenericArc<index_type>   Arc;
    private:

        typedef std::vector<NodeStorageType> NodeVector;
        typedef std::vector<EdgeStorageType> EdgeVector;
    public:

        typedef RagItemIt<RagType,Edge>                                                                 EdgeIt;
        typedef RagItemIt<RagType,Node>                                                                 NodeIt;
        typedef RagItemIt<RagType,Arc>                                                                  ArcIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorageType,detail::NeighborNodeFilter<RagType> >  NeighborNodeIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorageType,detail::IncEdgeFilter<RagType>      >  IncEdgeIt;

        typedef detail::GenericIncEdgeIt<RagType,NodeStorageType,detail::IsInFilter<RagType>      >     InArcIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorageType,detail::IsOutFilter<RagType>     >     OutArcIt;

        typedef UInt32 CoordinateElementType;
        typedef TinyVector<CoordinateElementType,4> EdgeCoordinate;
        typedef TinyVector<CoordinateElementType,2> NodeCoordinate;

        typedef std::vector<EdgeCoordinate> EdgeCoordinates;
        typedef std::vector<NodeCoordinate> NodeCoordinates;

        typedef std::vector<EdgeCoordinates>  EdgeCoordinatesVector;
        typedef std::vector<NodeCoordinates>  NodeCoordinatesVector;


    // public member functions
    public:
        Rag();
        Rag(const InputLabelingView & labels);
    



        // sizes
        index_type edgeNum()const;
        index_type nodeNum()const;
        index_type arcNum()const{return edgeNum()*2;}
        
        index_type maxEdgeId()const;
        index_type maxNodeId()const;
        index_type maxArcId()const{
            return maxEdgeId()*2;
        }

        Arc direct(const Edge & edge)const{
            return Arc(id(edge),id(edge));
        }

        Arc direct(const Edge & edge,const Node & node)const{
            if(u(edge)==node){
                std::cout<<"direct  no flip\n";
                return Arc(id(edge),id(edge));
            }
            else if(v(edge)==node){

                std::cout<<"direct   flip\n";
                return Arc(id(edge)+maxEdgeId(),id(edge));
            }
            else{
                std::cout<<"\n\n\n  WAAARNING !!! \n\n\n";
                return Arc(lemon::INVALID);
            }
        }

        // u / v

        Node u(const Edge & edge)const{
            return Node(edges_[id(edge)-1].u());
        }
        Node v(const Edge & edge)const{
            return Node(edges_[id(edge)-1].v());
        }


        Node source(const Arc & arc)const{
            const index_type arcIndex  = id(arc);
            if (arcIndex > maxEdgeId() ){
                std::cout<<"source from  flip\n";
                const index_type edgeIndex = arc.edgeId();
                const Edge edge = edgeFromId(edgeIndex);
                return v(edge);
            }
            else{
                std::cout<<"source from no flip\n";
                const index_type edgeIndex = arcIndex;
                const Edge edge = edgeFromId(edgeIndex);
                return u(edge);
            }
        }

        Node target(const Arc & arc)const{
            const index_type arcIndex  = id(arc);
            if (arcIndex > maxEdgeId() ){
                std::cout<<"target from flip\n";
                const index_type edgeIndex = arcIndex-maxEdgeId();
                const Edge edge = edgeFromId(edgeIndex);
                return u(edge);
            }
            else{
                std::cout<<"target from no flip\n";
                const index_type edgeIndex = arcIndex;
                const Edge edge = edgeFromId(edgeIndex);
                return v(edge);
            }
        }


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

        const NodeStorageType & nodeImpl(const Node & node)const{
        	return nodes_[id(node)-1];
        }


        // ids 
        index_type id(const Node & node)const;
        index_type id(const Edge & edge)const;
        index_type id(const Arc & arc)const{
            return arc.id();
        }

        // get edge / node from id
        Edge  edgeFromId(const index_type id)const;
        Node  nodeFromId(const index_type id)const;
        Arc   arcFromId(const index_type id)const{
            if(id>maxEdgeId()){
                return Arc(id,id);
            }
            else{
                return Arc(id,id-(maxEdgeId()+1));
            }
        }


        // find edge
        Edge findEdge(const Node & a,const Node & b)const;

        // iterators
        NodeIt nodesBegin()const;
        NodeIt nodesEnd()const  ;

        EdgeIt edgesBegin()const;
        EdgeIt edgesEnd()const  ;


        // extractors
        void extractEdgeCoordinates( EdgeCoordinatesVector & coordsVec )const;
        void extractNodeCoordinates( NodeCoordinatesVector & coordsVec )const;

    private:


        // graph
        NodeVector nodes_;
        EdgeVector edges_;

        // invalid edge and invalid node
        Node invalidNode_;
        Edge invalidEdge_;
    
        // superpixel / supervoxel labeling
        InputLabelingView  labeling_;

        UInt64 gridEdges_;
    };



    template<unsigned int DIM , class IN_LABEL_TYPE>
    Rag<DIM,IN_LABEL_TYPE>::Rag()
    :   nodes_(),
        edges_(),
        labeling_(),
        gridEdges_(0){

    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    Rag<DIM,IN_LABEL_TYPE>::Rag(
        const typename Rag<DIM,IN_LABEL_TYPE>::InputLabelingView & labels
    )
    :   nodes_(),
        edges_(),
        invalidNode_(),
        invalidEdge_(),
        labeling_(labels),
        gridEdges_(0)
    {
        
        // get min max label
        const UInt64 minLabel = *std::min_element(labels.begin(),labels.end());
        const UInt64 maxLabel = *std::max_element(labels.begin(),labels.end());
        if(minLabel!=1){
            throw std::runtime_error("minimum label must be 1");
        }

        // allocate space for nodes
        nodes_.resize(maxLabel);

        for(index_type nodeLabel=1;nodeLabel<=maxLabel;++nodeLabel){
            nodes_[nodeLabel-1].id_=nodeLabel;
        }

        // get grid graph from shape
        typedef GridGraph<DIM,boost::undirected_tag> GridGraphType;
        typedef typename GridGraphType::IndexMap GridGraphIndexMapType;
        typedef typename GridGraphType::Node     GridGraphNode;
        typedef typename GridGraphType::NodeIt   GridGraphGraphScanner;
        typedef typename GridGraphType::OutArcIt GridGraphNeighborIterator;

        GridGraphType   g(labeling_.shape());
        GridGraphIndexMapType indexMap = g.indexMap();

        std::set<UInt64> edgeSet;

        // iterate over all nodes (i.e. pixels)
        for (GridGraphGraphScanner node(g); node != lemon::INVALID; ++node){
            const InLabelType label = labeling_[indexMap[*node]];
            // iterate over all neighbors of the current node
            for (GridGraphNeighborIterator arc(g, node); arc != lemon::INVALID; ++arc){
                const GridGraphNode otherNode = g.target(*arc);
                const InLabelType otherLabel = labeling_[indexMap[otherNode]];
                if(otherLabel!=label){
                    const UInt64 l0 = std::min(label,otherLabel);
                    const UInt64 l1 = std::max(label,otherLabel);
                    const UInt64 key = l0 + l1*(maxLabel+1);
                    edgeSet.insert(key);
                    ++gridEdges_;
                }
            }
        }

        // allocate space for edges
        edges_.resize(edgeSet.size());
        typedef typename std::set<UInt64>::const_iterator EdgeSetIt;
        index_type edgeLabel=1;

        // fill graph
        for(EdgeSetIt iter=edgeSet.begin();iter!=edgeSet.end();++iter){
            const UInt64 key = *iter;
            const UInt64 l1 = key / (maxLabel+1);
            const UInt64 l0 = key - l1*(maxLabel+1);
            EdgeStorageType & edge  = edges_[edgeLabel-1];
            edge[0]=l0;
            edge[1]=l1;
            edge[2]=edgeLabel;

            NodeStorageType & n0 = nodes_[l0-1];
            NodeStorageType & n1 = nodes_[l1-1];
            n0.insertEdgeId(edgeLabel);
            n1.insertEdgeId(edgeLabel);
            ++edgeLabel;
        }

    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::edgeNum()const{
        return edges_.size();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::nodeNum()const{
        return nodes_.size();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::maxEdgeId()const{
        return edgeNum();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::maxNodeId()const{
        return nodeNum();
    }


    // ids 
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::id(
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return node.id();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::id(
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge & edge
    )const{
        return edge.id();
    }

    // get edge / node from id
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Edge 
    Rag<DIM,IN_LABEL_TYPE>::edgeFromId(
        const typename Rag<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        return Edge(id);
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::nodeFromId(
        const typename Rag<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        return Node(id);
    }




    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline  typename Rag<DIM,IN_LABEL_TYPE>::Edge  
    Rag<DIM,IN_LABEL_TYPE>::findEdge(
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & a,
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & b
    )const{
        std::pair<index_type,bool> res =  nodes_[id(a)-1].sharedEdge(nodes_[id(b)-1]);
        if(res.second){
            return Edge(res.first);
        }
        else{
            return Edge(-1);
        }
    }


    // iterators
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::NodeIt 
    Rag<DIM,IN_LABEL_TYPE>::nodesBegin()const{
        return NodeIt(0,nodeNum());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::NodeIt 
    Rag<DIM,IN_LABEL_TYPE>::nodesEnd()const{  
        return NodeIt(nodeNum(),nodeNum());
    }


    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::EdgeIt 
    Rag<DIM,IN_LABEL_TYPE>::edgesBegin()const{
        return EdgeIt(0,edgeNum());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::EdgeIt 
    Rag<DIM,IN_LABEL_TYPE>::edgesEnd()const{  
        return EdgeIt(edgeNum(),edgeNum());
    }




    template<unsigned int DIM , class IN_LABEL_TYPE>
    void Rag<DIM,IN_LABEL_TYPE>::extractEdgeCoordinates( 
        typename Rag<DIM,IN_LABEL_TYPE>::EdgeCoordinatesVector & coordsVec 
    )const{
        // resize coords Vec
        coordsVec.resize(edgeNum());
        // get grid graph from shape
        typedef GridGraph<DIM,boost::undirected_tag> GridGraphType;
        typedef typename GridGraphType::IndexMap GridGraphIndexMapType;
        typedef typename GridGraphType::Node     GridGraphNode;
        typedef typename GridGraphType::NodeIt   GridGraphGraphScanner;
        typedef typename GridGraphType::OutArcIt GridGraphNeighborIterator;
        typedef typename GridGraphIndexMapType::value_type GridGraphCoordType;
        GridGraphType   g(labeling_.shape());
        GridGraphIndexMapType indexMap = g.indexMap();



        // iterate over all nodes (i.e. pixels)
        for (GridGraphGraphScanner node(g); node != lemon::INVALID; ++node){

            const GridGraphCoordType & nodeCoord = indexMap[*node];
            const InLabelType label = labeling_[nodeCoord];

            // iterate over all neighbors of the current node
            bool is_local_maximum = true;
            for (GridGraphNeighborIterator arc(g, node); arc != lemon::INVALID; ++arc){
                const GridGraphNode otherNode = g.target(*arc);
                const GridGraphCoordType & otherNodeCoord = indexMap[*otherNode];
                const InLabelType otherLabel = labeling_[otherNodeCoord];

                if(otherLabel!=label){
                    // get the edge index
                    Edge  e = findEdge(nodeFromId(label),nodeFromId(otherLabel));
                    if(e==lemon::INVALID){
                        throw std::runtime_error("internal error");
                    }
                    const index_type edgeIndex = this->id(e);
                    const EdgeCoordinate edgeCoord( nodeCoord[0],nodeCoord[1],otherNodeCoord[0],otherNodeCoord[1]);
                    coordsVec[edgeIndex].push_back(edgeCoord);
                }
            }
        }
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    void Rag<DIM,IN_LABEL_TYPE>::extractNodeCoordinates( 
        typename Rag<DIM,IN_LABEL_TYPE>::NodeCoordinatesVector & coordsVec 
    )const{

    }

}

#endif /*VIGRA_RAG_HXX*/