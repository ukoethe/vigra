#ifndef VIGRA_HIERARCHICAL_CLUSTERING_HXX
#define VIGRA_HIERARCHICAL_CLUSTERING_HXX

//#include <valgrind/callgrind.h>

/*std*/
#include <queue>          

/*vigra*/

#include <vigra/merge_graph/merge_graph_adaptor.hxx>
#include <vigra/merge_graph/maps/merge_graph_view_map.hxx>
namespace vigra{      


    template<class G,class EDGE_INDICATOR_MAP,class EDGE_SIZE_MAP>
    class HierarchicalClustering{

    public:


        typedef EDGE_INDICATOR_MAP               EdgeIndicatorMap;
        typedef EDGE_SIZE_MAP                    EdgeSizeMap;
        typedef G                                Graph;
        typedef typename Graph::Edge             Edge;
        typedef typename Graph::Node             Node;
        typedef typename Graph::EdgeIt           EdgeIt;
        typedef MergeGraphAdaptor<Graph>         MergeGraph;
        typedef typename MergeGraph::Edge        MergeGraphEdge;
        typedef typename EdgeIndicatorMap::Value ValueType;
        typedef typename  EdgeSizeMap::Value     SizeType;
        typedef WeightedMeanMap<MergeGraph,EdgeIndicatorMap,EdgeSizeMap,MergeGraphEdge>     MgWeightedMeanView;
        typedef MinWeightEdgeViewMap<MergeGraph,MgWeightedMeanView>                         MgMinWeightOperator;
        typedef typename MergeGraph::index_type MergeGraphIndexType;


        struct Parameter{
            Parameter(
                const size_t nodeNumStopCond = 1,
                const bool   buildMergeTree  = true
            )
            :   nodeNumStopCond_ (nodeNumStopCond){                
            }

            size_t nodeNumStopCond_;
            bool   buildMergeTreeEncoding_;
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
            EdgeIndicatorMap & edgeIndicatorMap,
            EdgeSizeMap      & edgeSizeMap,
            Parameter          parameter =  Parameter()
        )
        :   graph_(graph),
            edgeIndicatorMap_(edgeIndicatorMap),
            edgeSizeMap_(edgeSizeMap),
            mergeGraph_(graph),
            mgWeightedMeanView_(mergeGraph_,edgeIndicatorMap,edgeSizeMap),
            minWeightOperator_(mergeGraph_,mgWeightedMeanView_),
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
            //mg.registerMergeEdgeCallBack(mgWeightedMeanView,& MgWeightedMeanView::merge);
            mergeGraph_.registerMergeEdgeCallBack(minWeightOperator_,& MgMinWeightOperator::mergeEdges);
            mergeGraph_.registerEraseEdgeCallBack(minWeightOperator_,& MgMinWeightOperator::eraseEdge);
            while(mergeGraph_.nodeNum()>1){
                if(mergeGraph_.nodeNum()%100==0)
                    std::cout<<"Nodes: "<<mergeGraph_.nodeNum()<<"\n";

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
                
            }
        }

        void transformInputMaps(){
            for(EdgeIt e(graph_);e!=lemon::INVALID;++e){
                const Edge reprEdge =  mergeGraph_.reprGraphEdge(*e);
                if(reprEdge!=*e){
                    edgeIndicatorMap_[*e]=edgeIndicatorMap_[reprEdge];
                    edgeSizeMap_[*e]     =edgeSizeMap_[reprEdge];
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
        EdgeIndicatorMap & edgeIndicatorMap_;
        EdgeSizeMap      & edgeSizeMap_;

        // merge graph / merge graph maps
        MergeGraph mergeGraph_;
        MgWeightedMeanView  mgWeightedMeanView_;
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