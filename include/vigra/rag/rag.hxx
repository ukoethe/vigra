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
#include <vigra/graph_helper/graph_iterator_helper.hxx>
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
        :   BaseIterType(0,0){
        }
        RagItemIt(const RagItemIt & other)
        :   BaseIterType(other){
        }

        // Invalid constructor & conversion. 
        RagItemIt(const lemon::Invalid & invalid)
        :   BaseIterType(0,0){
        }
        RagItemIt(const Graph & g)
        :   BaseIterType(0, GraphItemHelper<GRAPH,ITEM>::itemNum(g) ){
        }
        RagItemIt(const Graph & g,const ITEM & item)
        :   BaseIterType(g.id(item), GraphItemHelper<GRAPH,ITEM>::itemNum(g) ){
        }
    };





    template<unsigned int DIM , class IN_LABEL_TYPE>
    class RagImpl
    {
    private:
        typedef RagImpl<DIM,IN_LABEL_TYPE> RagType;
        typedef RagImpl<DIM,IN_LABEL_TYPE> RagImplType;
        
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
    private:

        typedef detail::filter::SmallerThan<index_type>                 BackFilter;
        typedef detail::filter::BackEdgeIdFilter<RagType>               BackEdgeIdFilter;
        typedef detail::transform::IdToItem<index_type,Edge  >          IdToEdgeTransform;
        typedef detail::transform::IdToItem<index_type,Node  >          IdToNodeTransform;
        typedef detail::transform::OtherNodeId<RagType>                 OtherNodeIdTransform;
        typedef detail::transform::OtherNode<RagType>                   OtherNodeTransform;
        typedef std::vector<NodeStorageType> NodeVector;
        typedef std::vector<EdgeStorageType> EdgeVector;
    public:

        typedef RagItemIt<RagImplType,Edge>                                                      EdgeIt;
        typedef RagItemIt<RagImplType,Node>                                                      NodeIt;

        typedef typename NodeStorageType::EdgeIdIt                                          NeighborEdgeIdIt;
        friend class detail::GenericIncEdgeIt<RagImplType,typename NodeStorageType::EdgeIdIt > ;



        typedef TransformIter<IdToEdgeTransform,NeighborEdgeIdIt,Edge,Edge>                 NeighborEdgeIt;
        typedef TransformIter<OtherNodeIdTransform,NeighborEdgeIdIt,index_type,index_type>  NeighborNodeIdIt;
        typedef TransformIter<OtherNodeTransform,NeighborEdgeIdIt,Node, Node >              NeighborNodeIt;

        typedef UInt32 CoordinateElementType;
        typedef TinyVector<CoordinateElementType,4> EdgeCoordinate;
        typedef TinyVector<CoordinateElementType,2> NodeCoordinate;

        typedef std::vector<EdgeCoordinate> EdgeCoordinates;
        typedef std::vector<NodeCoordinate> NodeCoordinates;

        typedef std::vector<EdgeCoordinates>  EdgeCoordinatesVector;
        typedef std::vector<NodeCoordinates>  NodeCoordinatesVector;


    // public member functions
    public:
        RagImpl();
        RagImpl(const InputLabelingView & labels);
    



        // sizes
        index_type edgeNum()const;
        index_type nodeNum()const;
        
        index_type maxEdgeId()const;
        index_type maxNodeId()const;


        // u / v

        Node u(const Edge & edge)const{
            return Node(edges_[id(edge)].u());
        }
        Node v(const Edge & edge)const{
            return Node(edges_[id(edge)].v());
        }


        // ids 
        index_type id(const Node & node)const;
        index_type id(const Edge & edge)const;


        // get edge / node from id
        Edge  edgeFromId(const index_type id)const;
        Node  nodeFromId(const index_type id)const;

        // find edge
        Edge findEdge(const Node & a,const Node & b)const;

        // iterators
        NodeIt nodesBegin()const;
        NodeIt nodesEnd()const  ;

        EdgeIt edgesBegin()const;
        EdgeIt edgesEnd()const  ;


        NeighborEdgeIdIt neigbourEdgeIdsPos(const Node & node,const Edge & edge)const{
            return nodes_[id(node)].find(id(edge));
        }
        NeighborEdgeIdIt neigbourEdgeIdsBegin(const Node & node)const;
        NeighborEdgeIdIt neigbourEdgeIdsEnd(const Node & node)const;

        NeighborEdgeIt neigbourEdgesBegin(const Node & node)const;
        NeighborEdgeIt neigbourEdgesEnd(const Node & node)const;


        NeighborNodeIdIt neigbourNodeIdsBegin(const Node & node)const;
        NeighborNodeIdIt neigbourNodeIdsEnd(const Node & node)const;

        NeighborNodeIt neigbourNodesBegin(const Node & node)const;
        NeighborNodeIt neigbourNodesEnd(const Node & node)const;


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
    class Rag
    :   public RagImpl<DIM,IN_LABEL_TYPE>,
        public detail::ArcHelper
        < 
            Rag<DIM,IN_LABEL_TYPE> ,                                // GRAPH
            typename RagImpl<DIM,IN_LABEL_TYPE>::index_type,        // INDEX_TYPE
            typename RagImpl<DIM,IN_LABEL_TYPE>::Edge ,             // EDGE
            typename RagImpl<DIM,IN_LABEL_TYPE>::Node ,             // NODE
            typename RagImpl<DIM,IN_LABEL_TYPE>::EdgeIt,            // EDGE_IT 
            typename RagImpl<DIM,IN_LABEL_TYPE>::NodeIt,            // NODE_IT 
            typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIdIt   // NEIGHBOR_EDGE_ID_IT

        >
    {
        typedef RagImpl<DIM,IN_LABEL_TYPE> RagImplType;
        typedef detail::ArcHelper
        < 
            Rag<DIM,IN_LABEL_TYPE> ,                                // GRAPH
            typename RagImpl<DIM,IN_LABEL_TYPE>::index_type,        // INDEX_TYPE
            typename RagImpl<DIM,IN_LABEL_TYPE>::Edge ,             // EDGE
            typename RagImpl<DIM,IN_LABEL_TYPE>::Node ,             // NODE
            typename RagImpl<DIM,IN_LABEL_TYPE>::EdgeIt,            // EDGE_IT 
            typename RagImpl<DIM,IN_LABEL_TYPE>::NodeIt,            // NODE_IT 
            typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIdIt   // NEIGHBOR_EDGE_ID_IT
        > ArcHelperType;

        public:
            Rag()
            :   RagImpl<DIM,IN_LABEL_TYPE>(){
            }
            Rag(const typename RagImpl<DIM,IN_LABEL_TYPE>::InputLabelingView & labels)
            :   RagImpl<DIM,IN_LABEL_TYPE>(labels){
            }

            typename RagImplType::index_type id(const typename ArcHelperType::Arc & arc)const{
                return this->arcId(arc);
            }
            typename RagImplType::index_type id(const typename RagImplType::Node & node)const{
                return static_cast<const RagImplType *>(this)->id(node);
            }

            typename RagImplType::index_type id(const typename RagImplType::Edge & edge)const{
                return static_cast<const RagImplType *>(this)->id(edge);
            }

    
    };

    template<unsigned int DIM , class IN_LABEL_TYPE>
    RagImpl<DIM,IN_LABEL_TYPE>::RagImpl()
    :   nodes_(),
        edges_(),
        labeling_(),
        gridEdges_(0){

    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    RagImpl<DIM,IN_LABEL_TYPE>::RagImpl(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::InputLabelingView & labels
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

        // allocate space for nodes
        nodes_.resize(maxLabel+1);

        for(index_type nodeIndex=0;nodeIndex<=maxLabel;++nodeIndex){
            nodes_[nodeIndex].id_=nodeIndex;
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
        index_type edgeIndex=0;

        // fill graph
        for(EdgeSetIt iter=edgeSet.begin();iter!=edgeSet.end();++iter){
            const UInt64 key = *iter;
            const UInt64 l1 = key / (maxLabel+1);
            const UInt64 l0 = key - l1*(maxLabel+1);
            EdgeStorageType & edge  = edges_[edgeIndex];
            edge[0]=l0;
            edge[1]=l1;
            edge[2]=edgeIndex;

            NodeStorageType & n0 = nodes_[l0];
            NodeStorageType & n1 = nodes_[l1];
            n0.insertEdgeId(edgeIndex);
            n1.insertEdgeId(edgeIndex);
            ++edgeIndex;
        }

    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::edgeNum()const{
        return edges_.size();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::nodeNum()const{
        return nodes_.size();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::maxEdgeId()const{
        return edgeNum()-1;
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::maxNodeId()const{
        return nodeNum()-1;
    }


    // ids 
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::id(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return node.id();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::index_type 
    RagImpl<DIM,IN_LABEL_TYPE>::id(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Edge & edge
    )const{
        return edge.id();
    }

    // get edge / node from id
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::Edge 
    RagImpl<DIM,IN_LABEL_TYPE>::edgeFromId(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        return Edge(id);
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::Node 
    RagImpl<DIM,IN_LABEL_TYPE>::nodeFromId(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        return Node(id);
    }




    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline  typename RagImpl<DIM,IN_LABEL_TYPE>::Edge  
    RagImpl<DIM,IN_LABEL_TYPE>::findEdge(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & a,
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & b
    )const{
        std::pair<index_type,bool> res =  nodes_[id(a)].sharedEdge(nodes_[id(b)]);
        if(res.second){
            return Edge(res.first);
        }
        else{
            return Edge(-1);
        }
    }


    // iterators
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NodeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::nodesBegin()const{
        return NodeIt(0,nodeNum());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NodeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::nodesEnd()const{  
        return NodeIt(nodeNum(),nodeNum());
    }


    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::EdgeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::edgesBegin()const{
        return EdgeIt(0,edgeNum());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::EdgeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::edgesEnd()const{  
        return EdgeIt(edgeNum(),edgeNum());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIdIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourEdgeIdsBegin(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return nodes_[id(node)].edgeIdsBegin();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIdIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourEdgeIdsEnd(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return nodes_[id(node)].edgeIdsEnd();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourEdgesBegin(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborEdgeIt(neigbourEdgeIdsBegin(),IdToEdgeTransform());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborEdgeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourEdgesEnd(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborEdgeIt(neigbourEdgeIdsEnd(),  IdToEdgeTransform());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborNodeIdIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourNodeIdsBegin(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborNodeIdIt(neigbourEdgeIdsBegin(),OtherNodeIdTransform(*this,node.id()));
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborNodeIdIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourNodeIdsEnd(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborNodeIdIt(neigbourEdgeIdsEnd(),  OtherNodeIdTransform(*this,node.id()));
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborNodeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourNodesBegin(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborNodeIt(neigbourEdgeIdsBegin(),OtherNodeTransform(*this,node.id()));
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename RagImpl<DIM,IN_LABEL_TYPE>::NeighborNodeIt 
    RagImpl<DIM,IN_LABEL_TYPE>::neigbourNodesEnd(
        const typename RagImpl<DIM,IN_LABEL_TYPE>::Node & node
    )const{
        return NeighborNodeIt(neigbourEdgeIdsBegin(),  OtherNodeTransform(*this,node.id()));
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    void RagImpl<DIM,IN_LABEL_TYPE>::extractEdgeCoordinates( 
        typename RagImpl<DIM,IN_LABEL_TYPE>::EdgeCoordinatesVector & coordsVec 
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
    void RagImpl<DIM,IN_LABEL_TYPE>::extractNodeCoordinates( 
        typename RagImpl<DIM,IN_LABEL_TYPE>::NodeCoordinatesVector & coordsVec 
    )const{

    }

}

#endif /*VIGRA_RAG_HXX*/