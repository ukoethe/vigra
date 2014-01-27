#ifndef VIGRA_HIERARCHICAL_CLUSTERING_HXX
#define VIGRA_HIERARCHICAL_CLUSTERING_HXX

//#include <valgrind/callgrind.h>

/*std*/
#include <queue>          
#include <iomanip>

/*vigra*/

#include <vigra/merge_graph/merge_graph_adaptor.hxx>
#include <vigra/merge_graph/maps/merge_graph_view_map.hxx>
namespace vigra{      


    template<
        class G,
        class EDGE_INDICATOR_MAP,  // edge indicator (gradmag)
        class EDGE_SIZE_MAP,       // size / length of edges/faces/boundaries
        class NODE_FEATURE_MAP,    // node features (as an rgb triple for each node/region)
        class NODE_SIZE_MAP,       // size of nodes as superpixel size
        class MIN_WEIGHT_OUT_MAP   // map to store the current out weight
    >
    class HierarchicalClustering{

    public:


        typedef EDGE_INDICATOR_MAP                  EdgeIndicatorMap;
        typedef EDGE_SIZE_MAP                       EdgeSizeMap;
        typedef NODE_FEATURE_MAP                    NodeFeatureMap;
        typedef NODE_SIZE_MAP                       NodeSizeMap;
        typedef MIN_WEIGHT_OUT_MAP                  MinWeightOutMap;

        typedef G                                       Graph;
        typedef typename Graph::Edge                    Edge;
        typedef typename Graph::Node                    Node;
        typedef typename Graph::EdgeIt                  EdgeIt;
        typedef typename Graph::NodeIt                  NodeIt;
        typedef MergeGraphAdaptor<Graph>                MergeGraph;
        typedef typename MergeGraph::Edge               MergeGraphEdge;
        typedef typename EdgeIndicatorMap::Value    ValueType;
        typedef typename  EdgeSizeMap::Value        SizeType;



        typedef HClusterMap<
            MergeGraph,
            EdgeIndicatorMap,
            EdgeSizeMap,
            NodeFeatureMap,
            NodeSizeMap,
            MinWeightOutMap
        >  MgMinWeightOperator;
        typedef typename MergeGraph::index_type MergeGraphIndexType;
        struct Parameter{
            Parameter(
                const size_t      nodeNumStopCond = 1,
                const bool        buildMergeTree  = true,
                const bool        verbose         = false,
                const ValueType   beta            = 0.5,
                const ValueType   degree1Fac      = 1.0,
                const ValueType   wardness        = 0.0,
                const size_t      nodeDist = MgMinWeightOperator::CHI_SQUARED_DISTANCE
            )
            :   nodeNumStopCond_ (nodeNumStopCond),
                buildMergeTreeEncoding_(buildMergeTree),
                verbose_(verbose),
                nodeDistType_(nodeDist),
                beta_(beta),
                degree1Fac_(degree1Fac),
                wardness_(wardness)
            {                
            }

            size_t nodeNumStopCond_;
            bool   buildMergeTreeEncoding_;
            bool   verbose_;
            size_t nodeDistType_;
            ValueType beta_;
            ValueType degree1Fac_;
            ValueType wardness_;
        };


        struct MergeItem{
            MergeItem(
                const MergeGraphIndexType  a,
                const MergeGraphIndexType  b,
                const MergeGraphIndexType  r,
                const ValueType            w
            ):
            a_(a),b_(b),r_(r),w_(w){
            }
            MergeGraphIndexType a_;
            MergeGraphIndexType b_;
            MergeGraphIndexType r_;
            ValueType           w_;
        };

        typedef std::vector<MergeItem> MergeTreeEncoding;


        HierarchicalClustering(
            const Graph      & graph,
            EdgeIndicatorMap    edgeIndicatorMap,
            EdgeSizeMap         edgeSizeMap,
            NodeFeatureMap      nodeFeatureMap,
            NodeSizeMap         nodeSizeMap,
            MinWeightOutMap     minWeightOutMap,
            Parameter          parameter =  Parameter()
        )
        :   graph_(graph),
            edgeIndicatorMap_(edgeIndicatorMap),
            edgeSizeMap_(edgeSizeMap),
            nodeFeatureMap_(nodeFeatureMap),
            nodeSizeMap_(nodeSizeMap),
            minWeightOutMap_(minWeightOutMap),
            mergeGraph_(graph),
            minWeightOperator_(mergeGraph_,edgeIndicatorMap,edgeSizeMap,nodeFeatureMap,nodeSizeMap,minWeightOutMap,
                parameter.nodeDistType_,parameter.beta_,parameter.degree1Fac_,parameter.wardness_),
            param_(parameter),
            timestamp_(graph.maxNodeId()+1),
            toTimeStamp_(graph.maxNodeId()+1),
            timeStampIndexToMergeIndex_(graph.maxNodeId()+1),
            mergeTreeEndcoding_()
        {
            // this can be be made smater since user can pass
            // stoping condition based on nodeNum
            mergeTreeEndcoding_.reserve(graph.nodeNum()*2);

            for(MergeGraphIndexType nodeId=0;nodeId<=mergeGraph_.maxNodeId();++nodeId){
                toTimeStamp_[nodeId]=nodeId;
            }
        }

        void cluster(){
            mergeGraph_.registerMergeNodeCallBack(minWeightOperator_,& MgMinWeightOperator::mergeNodes);
            mergeGraph_.registerMergeEdgeCallBack(minWeightOperator_,& MgMinWeightOperator::mergeEdges);
            mergeGraph_.registerEraseEdgeCallBack(minWeightOperator_,& MgMinWeightOperator::eraseEdge);
            if(param_.verbose_)
                std::cout<<"\n"; 
            while(mergeGraph_.nodeNum()>param_.nodeNumStopCond_ && mergeGraph_.edgeNum()>0){
                

                const Edge edgeToRemove = minWeightOperator_.minWeightEdge();
                if(param_.buildMergeTreeEncoding_){
                    const MergeGraphIndexType uid = mergeGraph_.id(mergeGraph_.u(edgeToRemove)); 
                    const MergeGraphIndexType vid = mergeGraph_.id(mergeGraph_.v(edgeToRemove));  
                    const ValueType w             = minWeightOperator_.minWeight();
                    // do the merge 
                    mergeGraph_.contractEdge( edgeToRemove);
                    const MergeGraphIndexType aliveNodeId = mergeGraph_.hasNodeId(uid) ? uid : vid;
                    const MergeGraphIndexType deadNodeId  = aliveNodeId==vid ? uid : vid;
                    timeStampIndexToMergeIndex_[timeStampToIndex(timestamp_)]=mergeTreeEndcoding_.size();
                    mergeTreeEndcoding_.push_back(MergeItem( toTimeStamp_[aliveNodeId],toTimeStamp_[deadNodeId],timestamp_,w));
                    toTimeStamp_[aliveNodeId]=timestamp_;
                    timestamp_+=1;
                }
                else{
                    // do the merge 
                    mergeGraph_.contractEdge( edgeToRemove );
                }
                if(param_.verbose_ && mergeGraph_.nodeNum()%10==0)
                    std::cout<<"\rNodes: "<<std::setw(10)<<mergeGraph_.nodeNum()<<std::flush;
                
            }
            if(param_.verbose_)
                std::cout<<"\n"; 
        }

        void transformInputMaps(){
            for(EdgeIt e(graph_);e!=lemon::INVALID;++e){
                const Edge reprEdge =  mergeGraph_.reprGraphEdge(*e);
                if(reprEdge!=*e){
                    edgeIndicatorMap_[*e]=edgeIndicatorMap_[reprEdge];
                    edgeSizeMap_[*e]     =edgeSizeMap_[reprEdge];
                    minWeightOutMap_[*e] =minWeightOutMap_[reprEdge];
                }
            }
            for(NodeIt n(graph_);n!=lemon::INVALID;++n){
                const Node reprNode = mergeGraph_.reprGraphNode(*n);
                if(reprNode!=*n){
                    nodeFeatureMap_[*n]=nodeFeatureMap_[reprNode];
                    nodeSizeMap_[*n]   =nodeSizeMap_[reprNode];
                }
            }
        }

        const MergeTreeEncoding & mergeTreeEndcoding()const{
            return mergeTreeEndcoding_;
        }

        template<class OUT_ITER>
        size_t leafNodeIds(const MergeGraphIndexType treeNodeId, OUT_ITER begin)const{
            if(treeNodeId<=graph_.maxNodeId()){
                *begin=treeNodeId;
                ++begin;
                return 1;
            }
            else{
                size_t leafNum=0;
                std::queue<MergeGraphIndexType>     queue;
                queue.push(treeNodeId);
                
                while(!queue.empty()){

                    const MergeGraphIndexType id = queue.front();
                    queue.pop();
                    const MergeGraphIndexType mergeIndex = timeStampToMergeIndex(id);
                    const MergeGraphIndexType ab[]= { mergeTreeEndcoding_[mergeIndex].a_, mergeTreeEndcoding_[mergeIndex].b_};

                    for(size_t i=0;i<2;++i){
                        if(ab[i]<=graph_.maxNodeId()){
                            *begin=ab[i];
                            ++begin;
                            ++leafNum;
                        }
                        else{
                            queue.push(ab[i]);
                        }
                    }
                }
                return leafNum;
            }
        }


        const Graph & graph()const{
            return graph_;
        }
        const MergeGraph & mergeGraph()const{
            return mergeGraph_;
        }

        const MergeGraphIndexType reprNodeId(const MergeGraphIndexType id)const{
            return mergeGraph_.reprNodeId(id);
        }
    private:
        // privat methods
        MergeGraphIndexType timeStampToIndex(const MergeGraphIndexType timestamp)const{
            return timestamp- graph_.maxNodeId();
        }


        MergeGraphIndexType timeStampToMergeIndex(const MergeGraphIndexType timestamp)const{
            return timeStampIndexToMergeIndex_[timeStampToIndex(timestamp)];
        }

        // input graph and input maps (maps will be modified inplace)
        const Graph &  graph_;
        EdgeIndicatorMap edgeIndicatorMap_;
        EdgeSizeMap      edgeSizeMap_;
        NodeFeatureMap   nodeFeatureMap_;
        NodeSizeMap      nodeSizeMap_; 
        MinWeightOutMap  minWeightOutMap_;

        // merge graph / merge graph maps
        MergeGraph mergeGraph_;
        MgMinWeightOperator minWeightOperator_;

        // parameter object
        Parameter          param_;

        // timestamp
        MergeGraphIndexType timestamp_;
        std::vector<MergeGraphIndexType> toTimeStamp_;
        std::vector<MergeGraphIndexType> timeStampIndexToMergeIndex_;
        // data which can reconstruct the merge tree
        MergeTreeEncoding mergeTreeEndcoding_;


};

}



#endif // VIGRA_HIERARCHICAL_CLUSTERING_HXX