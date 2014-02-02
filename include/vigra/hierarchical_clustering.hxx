#ifndef VIGRA_HIERARCHICAL_CLUSTERING_HXX
#define VIGRA_HIERARCHICAL_CLUSTERING_HXX

/*std*/
#include <queue>          
#include <iomanip>


namespace vigra{      





    template< class CLUSTER_OPERATOR>
    class HierarchicalClustering{

    public:
        typedef CLUSTER_OPERATOR                        ClusterOperator;
        typedef typename ClusterOperator::MergeGraph    MergeGraph;
        typedef typename MergeGraph::Graph              Graph;
        typedef typename Graph::Edge                    Edge;
        typedef typename Graph::Node                    Node;
        typedef typename Graph::EdgeIt                  EdgeIt;
        typedef typename Graph::NodeIt                  NodeIt;
        typedef typename MergeGraph::Edge               MergeGraphEdge;
        typedef typename CLUSTER_OPERATOR::WeightType   ValueType;
        typedef typename MergeGraph::index_type         MergeGraphIndexType;

        struct Parameter{
            Parameter(
                const size_t      nodeNumStopCond = 1,
                const bool        buildMergeTree  = true,
                const bool        verbose         = false
            )
            :   nodeNumStopCond_ (nodeNumStopCond),
                buildMergeTreeEncoding_(buildMergeTree),
                verbose_(verbose){                
            }
            size_t nodeNumStopCond_;
            bool   buildMergeTreeEncoding_;
            bool   verbose_;
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
            ClusterOperator & clusterOperator,
            const Parameter & parameter = Parameter()

        )
        :  
            clusterOperator_(clusterOperator),
            param_(parameter),
            mergeGraph_(clusterOperator_.mergeGraph()),
            graph_(mergeGraph_.graph()),
            timestamp_(graph_.maxNodeId()+1),
            toTimeStamp_(graph_.maxNodeId()+1),
            timeStampIndexToMergeIndex_(graph_.maxNodeId()+1),
            mergeTreeEndcoding_()
        {
            // this can be be made smater since user can pass
            // stoping condition based on nodeNum
            mergeTreeEndcoding_.reserve(graph_.nodeNum()*2);

            for(MergeGraphIndexType nodeId=0;nodeId<=mergeGraph_.maxNodeId();++nodeId){
                toTimeStamp_[nodeId]=nodeId;
            }
        }

        void cluster(){
            // ClusterOperator does registration by itself
            //mergeGraph_.registerMergeNodeCallBack(clusterOperator_,& MgMinWeightOperator::mergeNodes);
            //mergeGraph_.registerMergeEdgeCallBack(clusterOperator_,& MgMinWeightOperator::mergeEdges);
            //mergeGraph_.registerEraseEdgeCallBack(clusterOperator_,& MgMinWeightOperator::eraseEdge);
            if(param_.verbose_)
                std::cout<<"\n"; 
            while(mergeGraph_.nodeNum()>param_.nodeNumStopCond_ && mergeGraph_.edgeNum()>0){
                

                const Edge edgeToRemove = clusterOperator_.contractionEdge();
                if(param_.buildMergeTreeEncoding_){
                    const MergeGraphIndexType uid = mergeGraph_.id(mergeGraph_.u(edgeToRemove)); 
                    const MergeGraphIndexType vid = mergeGraph_.id(mergeGraph_.v(edgeToRemove));  
                    const ValueType w             = clusterOperator_.contractionWeight();
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

        //void transformInputMaps(){
        //    for(EdgeIt e(graph_);e!=lemon::INVALID;++e){
        //        const Edge reprEdge =  mergeGraph_.reprGraphEdge(*e);
        //        if(reprEdge!=*e){
        //            edgeIndicatorMap_[*e]=edgeIndicatorMap_[reprEdge];
        //            edgeSizeMap_[*e]     =edgeSizeMap_[reprEdge];
        //            minWeightOutMap_[*e] =minWeightOutMap_[reprEdge];
        //        }
        //    }
        //    for(NodeIt n(graph_);n!=lemon::INVALID;++n){
        //        const Node reprNode = mergeGraph_.reprGraphNode(*n);
        //        if(reprNode!=*n){
        //            nodeFeatureMap_[*n]=nodeFeatureMap_[reprNode];
        //            nodeSizeMap_[*n]   =nodeSizeMap_[reprNode];
        //        }
        //    }
        //}

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

        ClusterOperator & clusterOperator_;
        Parameter          param_;
        MergeGraph & mergeGraph_;
        const Graph  & graph_;
        // parameter object
        

        // timestamp
        MergeGraphIndexType timestamp_;
        std::vector<MergeGraphIndexType> toTimeStamp_;
        std::vector<MergeGraphIndexType> timeStampIndexToMergeIndex_;
        // data which can reconstruct the merge tree
        MergeTreeEncoding mergeTreeEndcoding_;


    };
}



#endif // VIGRA_HIERARCHICAL_CLUSTERING_HXX