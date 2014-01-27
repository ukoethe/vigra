#ifndef VIGRA_MERGE_GRAPH_VIEW_MAP_HXX
#define VIGRA_MERGE_GRAPH_VIEW_MAP_HXX

/* vigra - merge graph */
#include <vigra/priority_queue.hxx>
#include <vigra/merge_graph/distance.hxx>
namespace vigra{

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


template<
    class MERGE_GRAPH,
    class EDGE_INDICATOR_MAP,
    class EDGE_SIZE_MAP,
    class NODE_FEATURE_MAP,
    class NODE_SIZE_MAP,
    class MIN_WEIGHT_MAP
>
class HClusterMap{
    typedef typename EDGE_INDICATOR_MAP::Value ValueType;
    typedef MERGE_GRAPH MergeGraph;
    typedef typename MergeGraph::Graph Graph;




public:

    static const size_t CHI_SQUARED_DISTANCE=0;
    static const size_t NORM_SQUARED_DISTANCE=1;
    static const size_t NORM_DISTANCE=2;
  
    typedef typename Graph::Edge GraphEdge;
    typedef typename Graph::Node GraphNode;
    typedef typename MergeGraph::Edge Edge;
    typedef typename MergeGraph::Node Node;
    typedef typename MergeGraph::EdgeIt EdgeIt;
    typedef typename MergeGraph::NodeIt NodeIt;
    typedef typename MergeGraph::IncEdgeIt IncEdgeIt;
    typedef typename MergeGraph::index_type index_type;
    typedef MergeGraphItemHelper<MergeGraph,Edge> EdgeHelper;
    typedef MergeGraphItemHelper<MergeGraph,Node> NodeHelper;


    typedef typename EDGE_INDICATOR_MAP::Reference EdgeIndicatorReference;
    typedef typename NODE_FEATURE_MAP::Reference NodeFeatureReference;
    HClusterMap(
        const MergeGraph & mergeGraph,
        EDGE_INDICATOR_MAP  edgeIndicatorMap,
        EDGE_SIZE_MAP       edgeSizeMap,
        NODE_FEATURE_MAP    nodeFeatureMap,
        NODE_SIZE_MAP       nodeSizeMap,
        MIN_WEIGHT_MAP      minWeightEdgeMap,
        const size_t        nodeDistType,
        const ValueType     beta,
        const ValueType     degree1Fac,
        const ValueType     wardness
    )
    :   mergeGraph_(mergeGraph),
        edgeIndicatorMap_(edgeIndicatorMap),
        edgeSizeMap_(edgeSizeMap),
        nodeFeatureMap_(nodeFeatureMap),
        nodeSizeMap_(nodeSizeMap),
        minWeightEdgeMap_(minWeightEdgeMap),
        pq_(mergeGraph.maxEdgeId()+1),
        nodeDistType_(nodeDistType),
        beta_(beta),
        degree1Fac_(degree1Fac),
        wardness_(wardness)
    {
        for(EdgeIt e(mergeGraph);e!=lemon::INVALID;++e){
            const Edge edge = *e;
            const GraphEdge graphEdge=EdgeHelper::itemToGraphItem(mergeGraph_,edge);
            const index_type edgeId = mergeGraph_.id(edge);
            const ValueType currentWeight = this->getEdgeWeight(edge);
            pq_.push(edgeId,currentWeight);
            minWeightEdgeMap_[graphEdge]=currentWeight;
        }
    }

    void mergeEdges(const Edge & a,const Edge & b){
        // update features / weigts etc
        const GraphEdge aa=EdgeHelper::itemToGraphItem(mergeGraph_,a);
        const GraphEdge bb=EdgeHelper::itemToGraphItem(mergeGraph_,b);
        EdgeIndicatorReference va=edgeIndicatorMap_[aa];
        EdgeIndicatorReference vb=edgeIndicatorMap_[bb];
        va*=edgeSizeMap_[aa];
        vb*=edgeSizeMap_[bb];
        va+=vb;
        edgeSizeMap_[aa]+=edgeSizeMap_[bb];
        va/=(edgeSizeMap_[aa]);
        vb/=edgeSizeMap_[bb];
        // delete b from pq
        pq_.deleteItem(b.id());
    }
    void mergeNodes(const Node & a,const Node & b){
        const GraphNode aa=NodeHelper::itemToGraphItem(mergeGraph_,a);
        const GraphNode bb=NodeHelper::itemToGraphItem(mergeGraph_,b);
        NodeFeatureReference va=nodeFeatureMap_[aa];
        NodeFeatureReference vb=nodeFeatureMap_[bb];
        va*=nodeSizeMap_[aa];
        vb*=nodeSizeMap_[bb];
        va+=vb;
        nodeSizeMap_[aa]+=nodeSizeMap_[bb];
        va/=(nodeSizeMap_[aa]);
        vb/=nodeSizeMap_[bb];
    }
    void eraseEdge(const Edge & edge){
        // delete edge from pq
        pq_.deleteItem(edge.id());
        // get the new region the edge is in
        // (since the edge is no any more an active edge)
        const Node newNode = mergeGraph_.inactiveEdgesNode(edge);

        // iterate over all edges of this node
        for (IncEdgeIt e(mergeGraph_,newNode);e!=lemon::INVALID;++e){
            const Edge      incEdge(*e);
            const GraphEdge incGraphEdge = EdgeHelper::itemToGraphItem(mergeGraph_,incEdge);
            // compute the new weight for this edge 
            // (this should involve region differences)
            const ValueType newWeight = getEdgeWeight(incEdge);
            // change the weight in pq by repushing
            pq_.push(incEdge.id(),newWeight);
            // remember edge weight
            minWeightEdgeMap_[incGraphEdge]=newWeight;
        } 
    }
    Edge minWeightEdge(){
        index_type minLabel = pq_.top();
        while(mergeGraph_.hasEdgeId(minLabel)==false){
            pq_.deleteItem(minLabel);
            index_type minLabel = pq_.top();
        }
        return Edge(minLabel);
    }
    ValueType minWeight()const{
        return pq_.topPriority();
    }
private:
    ValueType getEdgeWeight(const Edge & e){

        

        const Node u  = mergeGraph_.u(e);
        const Node v  = mergeGraph_.v(e);
        const GraphEdge ee=EdgeHelper::itemToGraphItem(mergeGraph_,e);
        const GraphNode uu=NodeHelper::itemToGraphItem(mergeGraph_,u);
        const GraphNode vv=NodeHelper::itemToGraphItem(mergeGraph_,v);

        const size_t degU = mergeGraph_.degree(u);
        const size_t degV = mergeGraph_.degree(v);

        const ValueType degFac     = degU==1 || degV==1 ? degree1Fac_ : 1.0;
        const ValueType wardFacRaw = 1.0 / ( 1.0/nodeSizeMap_[uu] + 1.0/nodeSizeMap_[vv]  );
        const ValueType wardFac    = (wardFacRaw*wardness_) + (1.0-wardness_);

        const ValueType fromEdgeIndicator = edgeIndicatorMap_[ee];
        ValueType fromNodeDist;
        if(nodeDistType_==NORM_DISTANCE){
            distances::Norm<ValueType> nodeDistFunctor;      
            fromNodeDist= nodeDistFunctor(nodeFeatureMap_[uu],nodeFeatureMap_[vv]);
        }
        else if(nodeDistType_==NORM_SQUARED_DISTANCE){
            distances::SquaredNorm<ValueType> nodeDistFunctor;      
            fromNodeDist= nodeDistFunctor(nodeFeatureMap_[uu],nodeFeatureMap_[vv]);
        }
        else if(nodeDistType_==CHI_SQUARED_DISTANCE){
            distances::ChiSquared<ValueType> nodeDistFunctor;      
            fromNodeDist= nodeDistFunctor(nodeFeatureMap_[uu],nodeFeatureMap_[vv]);
        }
        else{
            throw std::runtime_error("wrong distance type");
        }
        const ValueType totalWeight  = ((1.0-beta_)*fromEdgeIndicator + beta_*fromNodeDist)*degFac*wardFac;
        return totalWeight;
    }


    const MergeGraph &  mergeGraph_;
    EDGE_INDICATOR_MAP  edgeIndicatorMap_;
    EDGE_SIZE_MAP       edgeSizeMap_;
    NODE_FEATURE_MAP    nodeFeatureMap_;
    NODE_SIZE_MAP       nodeSizeMap_;
    MIN_WEIGHT_MAP      minWeightEdgeMap_;
    vigra::ChangeablePriorityQueue< ValueType > pq_;

    size_t        nodeDistType_;
    ValueType     beta_;
    ValueType     degree1Fac_;
    ValueType     wardness_;
};



} // end namespace vigra
#endif // VIGRA_MERGE_GRAPH_VIEW_MAP_HXX