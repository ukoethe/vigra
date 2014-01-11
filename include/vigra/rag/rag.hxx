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

//#include <vigra/is_end_iterator.hxx>
//#include <vigra/is_end_transform_iterator.hxx>
//#include <vigra/is_end_filter_iterator.hxx>
//#include <vigra/is_end_enumeration_iterator.hxx>
//#include <vigra/graph_iterator_functors.hxx>
//#include <vigra/graph_helper/graph_crtp_base.hxx>


#include <vigra/algorithm.hxx>
#include <vigra/graph_helper/graph_item_impl.hxx>

#include <vigra/rag/rag_iterators.hxx>

namespace vigra{





    template<unsigned int DIM , class IN_LABEL_TYPE>
    class Rag
    {
        
    public:
        // public typdedfs
        typedef IN_LABEL_TYPE                                             InLabelType;
        typedef MultiArrayView<DIM,InLabelType>                           InputLabelingView;
        typedef MultiArray<DIM,InLabelType>                               InputLabelingArray;
        typedef Int64                                                     index_type;
        const static unsigned int Dimension = DIM;
    private:
        // private typedes which are needed for defining public typedes
        typedef Rag<DIM,IN_LABEL_TYPE>                                    RagType;
        typedef RandomAccessSet<index_type>                               NodeStorageEdgeSet;
        typedef detail::GenericNodeImpl<index_type,NodeStorageEdgeSet>    NodeStorage;
        typedef detail::GenericEdgeImpl<index_type >                      EdgeStorageType;
        typedef detail::NeighborNodeFilter<RagType>                       NnFilter;
        typedef detail::IncEdgeFilter<RagType>                            IncFilter;
        typedef detail::IsInFilter<RagType>                               InFlter;
        typedef detail::IsOutFilter<RagType>                              OutFilter;
    public:
        // LEMON API TYPEDEFS (and a few more(NeighborNodeIt))
        typedef detail::GenericNode<index_type>                           Node;
        typedef detail::GenericEdge<index_type>                           Edge;
        typedef detail::GenericArc<index_type>                            Arc;
        typedef detail::RagItemIt<RagType,Edge>                           EdgeIt;
        typedef detail::RagItemIt<RagType,Node>                           NodeIt;
        typedef detail::RagItemIt<RagType,Arc>                            ArcIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorage,NnFilter  >  NeighborNodeIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorage,IncFilter >  IncEdgeIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorage,InFlter   >  InArcIt;
        typedef detail::GenericIncEdgeIt<RagType,NodeStorage,OutFilter >  OutArcIt;

        // RAG's addition typedefs 
        typedef UInt32                                                    CoordinateElementType;
        typedef TinyVector<CoordinateElementType,DIM>                     EdgeCoordinate;
        typedef TinyVector<CoordinateElementType,DIM  >                   NodeCoordinate;
        typedef std::vector<EdgeCoordinate>                               EdgeCoordinates;
        typedef std::vector<NodeCoordinate>                               NodeCoordinates;


        // predefined map types 
        typedef DenseEdgeReferenceMap<RagType,EdgeCoordinates>            EdgeCoordinatesMap;
        typedef DenseNodeReferenceMap<RagType,NodeCoordinates>            NodeCoordinatesMap;


    // public member functions
    public:
        Rag();
        Rag(const InputLabelingView & labels);

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


        // extractors
        void extractEdgeCoordinates( EdgeCoordinatesMap & coordsVec )const;
        void extractNodeCoordinates( NodeCoordinatesMap & coordsVec )const;

    private:
        // private typedefs
        typedef std::vector<NodeStorage> NodeVector;
        typedef std::vector<EdgeStorageType> EdgeVector;


        // needs acces to const nodeImpl
        template<class G,class NIMPL,class FILT>
        friend class detail::GenericIncEdgeIt;


        const NodeStorage & nodeImpl(const Node & node)const{
            return nodes_[id(node)-1];
        }


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

            NodeStorage & n0 = nodes_[l0-1];
            NodeStorage & n1 = nodes_[l1-1];
            n0.insertEdgeId(edgeLabel);
            n1.insertEdgeId(edgeLabel);
            ++edgeLabel;
        }

    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Arc 
    Rag<DIM,IN_LABEL_TYPE>::direct(
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge & edge,
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

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Arc 
    Rag<DIM,IN_LABEL_TYPE>::direct(
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge & edge,
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & node
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

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline bool
    Rag<DIM,IN_LABEL_TYPE>::direction(
        const typename Rag<DIM,IN_LABEL_TYPE>::Arc & arc
    )const{
        return id(arc)<=maxEdgeId();
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::u(
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge & edge
    )const{
        return Node(edges_[id(edge)-1].u());
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node
    Rag<DIM,IN_LABEL_TYPE>::v(
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge & edge
    )const{
        return Node(edges_[id(edge)-1].v());
    }


    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::source(
        const typename Rag<DIM,IN_LABEL_TYPE>::Arc & arc
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


    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::target(
        const typename Rag<DIM,IN_LABEL_TYPE>::Arc & arc
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

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node
    Rag<DIM,IN_LABEL_TYPE>::oppositeNode(
        const typename Rag<DIM,IN_LABEL_TYPE>::Node &n,
        const typename Rag<DIM,IN_LABEL_TYPE>::Edge &e
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


    template<unsigned int DIM , class IN_LABEL_TYPE>
    typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::baseNode(
        const typename Rag<DIM,IN_LABEL_TYPE>::IncEdgeIt & iter
    )const{
        return u(*iter);
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::baseNode(
        const typename Rag<DIM,IN_LABEL_TYPE>::OutArcIt & iter 
    )const{
        return source(*iter);
    }


    template<unsigned int DIM , class IN_LABEL_TYPE>
    typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::runningNode(
        const typename Rag<DIM,IN_LABEL_TYPE>::IncEdgeIt & iter
    )const{
        return v(*iter);
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::runningNode(
        const typename Rag<DIM,IN_LABEL_TYPE>::OutArcIt & iter 
    )const{
        return target(*iter);
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
    Rag<DIM,IN_LABEL_TYPE>::arcNum()const{
        return edgeNum()*2;
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

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::maxArcId()const{
        return maxEdgeId()*2;
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

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::index_type 
    Rag<DIM,IN_LABEL_TYPE>::id(
        const typename Rag<DIM,IN_LABEL_TYPE>::Arc & arc
    )const{
        return arc.id();
    }

    // get edge / node from id
    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Edge 
    Rag<DIM,IN_LABEL_TYPE>::edgeFromId(
        const typename Rag<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        if(id>=1 && id<=maxEdgeId())
            return Edge(id);
        else{
            return Edge(lemon::INVALID);
        }
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Node 
    Rag<DIM,IN_LABEL_TYPE>::nodeFromId(
        const typename Rag<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        if(id>=1 && id<=maxNodeId())
            return Node(id);
        else{
            return Node(lemon::INVALID);
        }
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline typename Rag<DIM,IN_LABEL_TYPE>::Arc 
    Rag<DIM,IN_LABEL_TYPE>::arcFromId(
        const typename Rag<DIM,IN_LABEL_TYPE>::index_type id
    )const{
        if(id>=1 && id<=maxArcId()){
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
            return Edge(lemon::INVALID);
        }
    }


    template<unsigned int DIM , class IN_LABEL_TYPE>
    inline  typename Rag<DIM,IN_LABEL_TYPE>::Arc  
    Rag<DIM,IN_LABEL_TYPE>::findArc(
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & uNode,
        const typename Rag<DIM,IN_LABEL_TYPE>::Node & vNode
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
        typename Rag<DIM,IN_LABEL_TYPE>::EdgeCoordinatesMap & coordMap 
    )const{

        // get min max label (todo , remove this / refactore minLabel-maxLabel)
        const UInt64 minLabel = *std::min_element(labeling_.begin(),labeling_.end());
        const UInt64 maxLabel = *std::max_element(labeling_.begin(),labeling_.end());
        if(minLabel!=1){
            throw std::runtime_error("minimum label must be 1");
        }


        // get grid graph types
        typedef GridGraph<DIM,boost::undirected_tag> GridGraphType;
        typedef typename GridGraphType::IndexMap     GridGraphIndexMapType;
        typedef typename GridGraphType::Node         GridGraphNode;
        typedef typename GridGraphType::NodeIt       GridGraphGraphScanner;
        typedef typename GridGraphType::OutArcIt     GridGraphNeighborIterator;

        // construct grid graph
        GridGraphType   g(labeling_.shape());
        GridGraphIndexMapType indexMap = g.indexMap();



        // iterate over all nodes (i.e. pixels)
        for (GridGraphGraphScanner node(g); node != lemon::INVALID; ++node){


            const InLabelType label = labeling_[indexMap[*node]];
            // iterate over all neighbors of the current node
            for (GridGraphNeighborIterator arc(g, node); arc != lemon::INVALID; ++arc){
                const GridGraphNode otherNode = g.target(*arc);
                const InLabelType otherLabel = labeling_[indexMap[otherNode]];

                if(otherLabel!=label){
                    // find the the edge
                    const Edge edge = findEdge(nodeFromId(label),nodeFromId(otherLabel));
                    // remove when unit tests are finished
                    if(edge==lemon::INVALID){
                        throw std::runtime_error("internal error");
                    }
                    // add coordinates
                    EdgeCoordinate e1=indexMap[*node];
                    EdgeCoordinate e2=indexMap[otherNode];
                    EdgeCoordinates & coords = coordMap[edge];
                    coords.push_back(e1);
                    coords.push_back(e2);
                }
            }
        }
    }

    template<unsigned int DIM , class IN_LABEL_TYPE>
    void Rag<DIM,IN_LABEL_TYPE>::extractNodeCoordinates( 
        typename Rag<DIM,IN_LABEL_TYPE>::NodeCoordinatesMap & coordsVec 
    )const{

    }



    /// FREE FUNCTIONS / TRAIT CLASSES

    template<class G>
    struct GraphCoordinateTraits{
        typedef typename G::EdgeCoordinatesMap EdgeCoordinatesMap;
        typedef typename G::NodeCoordinatesMap NodeCoordinatesMap;
    };



    template<class G>
    void extractEdgeCoordinates(const G & g, typename GraphCoordinateTraits<G>::EdgeCoordinatesMap & coordMap){
        throw std::runtime_error("need to be specialized for GRAPH");
    }

    template<unsigned int DIM,class IN_LABEL_TYPE>
    void extractEdgeCoordinates(
        const Rag<DIM,IN_LABEL_TYPE> & g, typename GraphCoordinateTraits< Rag<DIM,IN_LABEL_TYPE> >::EdgeCoordinatesMap & coordMap
    ){
        g.extractEdgeCoordinates(coordMap);
    }




    template<class G,class COORD_MAP, class IMAGE_TYPE,class OUT_MAP>
    void extractEdgeFeaturesFromImage(
        const G &             g,
        const COORD_MAP &     coordMap,
        const IMAGE_TYPE &    image,
        OUT_MAP &             outMap
    ){

        // typedefs of graph
        typedef typename G::EdgeIt EdgeIt;
        typedef typename G::Edge   Edge;

        // typedefs of coord map
        typedef typename COORD_MAP::ConstReference      ConstCordVecRef;
        typedef typename COORD_MAP::Value::value_type   Coordinate;

        // typedefs of the image
        typedef typename IMAGE_TYPE::value_type ImageValueType;

        // iterate over all edges
        for(EdgeIt e(g);e!=lemon::INVALID;++e){
            const Edge edge = *e;
            // get a const reference to the coordinates for the edge
            ConstCordVecRef coordVec = coordMap[edge];


            const size_t coordNum = std::distance(coordVec.begin(),coordVec.end());
            // iterate over all coordinates
            size_t c=0;
            ImageValueType accVal=image[*coordVec.begin()];
            BOOST_FOREACH( const Coordinate & coord, coordVec )
            {   
                if(c==0)
                    accVal=image[coord];
                else
                    accVal+=image[coord];
                ++c;
            }
            accVal/=coordNum;
            outMap[edge]=accVal;

        }
    }

}

#endif /*VIGRA_RAG_HXX*/