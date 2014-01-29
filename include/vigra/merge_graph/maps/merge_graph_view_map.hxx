

#ifndef VIGRA_MERGE_GRAPH_VIEW_MAP_HXX
#define VIGRA_MERGE_GRAPH_VIEW_MAP_HXX

/* vigra - merge graph */
#include <vigra/merge_graph/min_indexed_pq.hxx>

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







template<class MERGE_GRAPH,class FEATURE_MAP,class WEIGHT_MAP,class MG_ITEM_TYPE>
class WeightedMeanMap{


    typedef MERGE_GRAPH MergeGraph;
    typedef MG_ITEM_TYPE MergeGraphItem;
    typedef MergeGraphItemHelper<MergeGraph,MergeGraphItem> ItemHelper;
    typedef FEATURE_MAP FeatureMap;
    typedef WEIGHT_MAP  WeightMap;

public:
    typedef typename MergeGraph::Edge           Edge;
    typedef typename ItemHelper::Item           Key;
    typedef typename FeatureMap::Value          Value;
    typedef typename FeatureMap::Reference      Reference;
    typedef typename FeatureMap::ConstReference ConstReference;

    typedef typename WeightMap::Value          WeightValue;
    typedef typename WeightMap::Reference      WeightReference;
    typedef typename WeightMap::ConstReference WeightConstReference;

    WeightedMeanMap(const MergeGraph & mg,FeatureMap & featureMap,WeightMap & weightMap)
    :   mg_(mg),
        featureMap_(featureMap),
        weightMap_(weightMap){
    }

    Reference operator[](const Key & key){
        return featureMap_[ItemHelper::itemToGraphItem(mg_,key)];
    }
    const ConstReference operator[](const Key & key)const{
        return featureMap_[ItemHelper::itemToGraphItem(mg_,key)];
    }

    void merge(const Key &  a , const Key & b){
        Reference va=this->operator[](a);
        Reference vb=this->operator[](b);
        va*=weight(a);
        vb*=weight(b);
        va+=vb;
        weight(a)+=weight(b);
        va/=(weight(a));
        vb/=weight(b);
    }
private:

    WeightReference weight(const Key & key){
        return weightMap_[ItemHelper::itemToGraphItem(mg_,key)];
    }
    const WeightConstReference weight(const Key & key)const{
        return weightMap_[ItemHelper::itemToGraphItem(mg_,key)];
    }


    const MergeGraph & mg_;
    FeatureMap & featureMap_;
    WeightMap  & weightMap_;

};



template<class MERGE_GRAPH,class EDGE_MAP>
class MinWeightEdgeViewMap 
{
private:
    MinWeightEdgeViewMap();                                           // non empty-construction
    MinWeightEdgeViewMap( const MinWeightEdgeViewMap& other );        // non construction-copyable
    MinWeightEdgeViewMap& operator=( const MinWeightEdgeViewMap& );   // non copyable
public:
    typedef MERGE_GRAPH MergeGraphType;
    typedef typename MergeGraphType::Edge   Edge;
    typedef typename MergeGraphType::EdgeIt EdgeIt;
    typedef typename MergeGraphType::index_type index_type;
    typedef typename EDGE_MAP::Value Value;

    void operator[](const Edge & edge)const{
        return edgeMap_[edge];
    }


    void mergeEdges(const Edge & a,const Edge & b){
        edgeMap_.merge(a,b);
        pq_.deleteValue(b.id());
        changePqWeight(a.id(),edgeMap_[a]);
    } 
    void eraseEdge(const Edge & edge){
        pq_.deleteValue(edge.id());
    }

    template<class CB>
    CB eraseEdgeCallback(){
        return  boost::bind(boost::mem_fn(&MinWeightEdgeViewMap<MERGE_GRAPH,EDGE_MAP>::eraseEdge), this , _1);
    }

    template<class CB>
    CB mergeEdgeCallback(){
        return  boost::bind(boost::mem_fn(&MinWeightEdgeViewMap<MERGE_GRAPH,EDGE_MAP>::mergeEdges), this , _1,_2);
    }

    MinWeightEdgeViewMap(const MergeGraphType & mergeGraph,EDGE_MAP & edgeMap) 
    :   mergeGraph_(mergeGraph),
        edgeMap_(edgeMap),
        pq_(mergeGraph.maxEdgeId()+1){
        for(EdgeIt e(mergeGraph);e!=lemon::INVALID;++e){
            const Edge edge = *e;
            const index_type edgeId = mergeGraph_.id(edge);
            pq_.insert(edgeId,edgeMap_[edge]);
        }
    }

    Edge minWeightEdge(){
        index_type minLabel = pq_.minIndex();
        while(mergeGraph_.hasEdgeId(minLabel)==false){
            pq_.deleteValue(minLabel);
            index_type minLabel = pq_.minIndex();
        }
        return Edge(minLabel);
    }

    Value minWeight()const{
        return pq_.minValue();
    }


private:
    void changePqWeight(const index_type l,const typename EDGE_MAP::Value newWeigt){
        pq_.changeValue(l,newWeigt);
    }
    const MergeGraphType & mergeGraph_;
    EDGE_MAP & edgeMap_;
    vigra::MinIndexedPQ<typename EDGE_MAP::Value>     pq_;
};













}

#endif // VIGRA_MERGE_GRAPH_VIEW_MAP_HXX