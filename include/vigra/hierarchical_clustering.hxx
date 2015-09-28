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



#ifndef VIGRA_HIERARCHICAL_CLUSTERING_HXX
#define VIGRA_HIERARCHICAL_CLUSTERING_HXX



/*std*/
#include <queue>          
#include <iomanip>

/*vigra*/
#include "priority_queue.hxx"
#include "metrics.hxx" 

namespace vigra{      

namespace cluster_operators{







template<
        class MERGE_GRAPH,
        class EDGE_INDICATOR_MAP,
        class EDGE_SIZE_MAP,
        class NODE_SIZE_MAP,
        class MIN_WEIGHT_MAP
    >
    class EdgeWeightedUcm{
        
        typedef EdgeWeightedUcm<
            MERGE_GRAPH,
            EDGE_INDICATOR_MAP,
            EDGE_SIZE_MAP,
            NODE_SIZE_MAP,
            MIN_WEIGHT_MAP
        > SelfType;
    public:


        typedef typename EDGE_INDICATOR_MAP::Value ValueType;
        typedef ValueType WeightType;
        typedef MERGE_GRAPH MergeGraph;
        typedef typename MergeGraph::Graph Graph;
        typedef typename Graph::Edge BaseGraphEdge;
        typedef typename Graph::Node BaseGraphNode;
        typedef typename MergeGraph::Edge Edge;
        typedef typename MergeGraph::Node Node;
        typedef typename MergeGraph::EdgeIt EdgeIt;
        typedef typename MergeGraph::NodeIt NodeIt;
        typedef typename MergeGraph::IncEdgeIt IncEdgeIt;
        typedef typename MergeGraph::index_type index_type;
        typedef MergeGraphItemHelper<MergeGraph,Edge> EdgeHelper;
        typedef MergeGraphItemHelper<MergeGraph,Node> NodeHelper;

        typedef typename MergeGraph::MergeNodeCallBackType MergeNodeCallBackType;
        typedef typename MergeGraph::MergeEdgeCallBackType MergeEdgeCallBackType;
        typedef typename MergeGraph::EraseEdgeCallBackType EraseEdgeCallBackType;

        typedef typename EDGE_INDICATOR_MAP::Reference EdgeIndicatorReference;
        /// \brief construct cluster operator
        EdgeWeightedUcm(
            MergeGraph & mergeGraph,
            EDGE_INDICATOR_MAP edgeIndicatorMap,
            EDGE_SIZE_MAP edgeSizeMap,
            NODE_SIZE_MAP nodeSizeMap,
            MIN_WEIGHT_MAP minWeightEdgeMap,
            const ValueType wardness=1.0
        )
        :   mergeGraph_(mergeGraph),
            edgeIndicatorMap_(edgeIndicatorMap),
            edgeSizeMap_(edgeSizeMap),
            nodeSizeMap_(nodeSizeMap),
            minWeightEdgeMap_(minWeightEdgeMap),
            pq_(mergeGraph.maxEdgeId()+1),
            wardness_(wardness)
        {



            MergeNodeCallBackType cbMn(MergeNodeCallBackType:: template from_method<SelfType,&SelfType::mergeNodes>(this));
            MergeEdgeCallBackType cbMe(MergeEdgeCallBackType:: template from_method<SelfType,&SelfType::mergeEdges>(this));
            EraseEdgeCallBackType cbEe(EraseEdgeCallBackType:: template from_method<SelfType,&SelfType::eraseEdge>(this));

            mergeGraph_.registerMergeNodeCallBack(cbMn);
            mergeGraph_.registerMergeEdgeCallBack(cbMe);
            mergeGraph_.registerEraseEdgeCallBack(cbEe);

          

            for(EdgeIt e(mergeGraph_);e!=lemon::INVALID;++e){
                const Edge edge = *e;
                const BaseGraphEdge graphEdge=EdgeHelper::itemToGraphItem(mergeGraph_,edge);
                const index_type edgeId = mergeGraph_.id(edge);
                const ValueType currentWeight = this->getEdgeWeight(edge);
                pq_.push(edgeId,currentWeight);
                minWeightEdgeMap_[graphEdge]=currentWeight;
            }

        }

        void setWardness(const float w){
            wardness_ = w;
        }
        void resetMgAndPq(){
            mergeGraph_.reset();


            MergeNodeCallBackType cbMn(MergeNodeCallBackType:: template from_method<SelfType,&SelfType::mergeNodes>(this));
            MergeEdgeCallBackType cbMe(MergeEdgeCallBackType:: template from_method<SelfType,&SelfType::mergeEdges>(this));
            EraseEdgeCallBackType cbEe(EraseEdgeCallBackType:: template from_method<SelfType,&SelfType::eraseEdge>(this));

            mergeGraph_.registerMergeNodeCallBack(cbMn);
            mergeGraph_.registerMergeEdgeCallBack(cbMe);
            mergeGraph_.registerEraseEdgeCallBack(cbEe);


            pq_.reset();
            for(EdgeIt e(mergeGraph_);e!=lemon::INVALID;++e){
                const Edge edge = *e;
                const BaseGraphEdge graphEdge=EdgeHelper::itemToGraphItem(mergeGraph_,edge);
                const index_type edgeId = mergeGraph_.id(edge);
                const ValueType currentWeight = this->getEdgeWeight(edge);
                pq_.push(edgeId,currentWeight);
                minWeightEdgeMap_[graphEdge]=currentWeight;
            }
        }

        /// \brief will be called via callbacks from mergegraph
        void mergeEdges(const Edge & a,const Edge & b){
            // update features / weigts etc
            const BaseGraphEdge aa=EdgeHelper::itemToGraphItem(mergeGraph_,a);
            const BaseGraphEdge bb=EdgeHelper::itemToGraphItem(mergeGraph_,b);
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

        /// \brief will be called via callbacks from mergegraph
        void mergeNodes(const Node & a,const Node & b){
            const BaseGraphNode aa=NodeHelper::itemToGraphItem(mergeGraph_,a);
            const BaseGraphNode bb=NodeHelper::itemToGraphItem(mergeGraph_,b);
            nodeSizeMap_[aa]+=nodeSizeMap_[bb];
        }

        /// \brief will be called via callbacks from mergegraph
        void eraseEdge(const Edge & edge){

            //std::cout<<"start to erase edge "<<mergeGraph_.id(edge)<<"\n";
            // delete edge from pq
            pq_.deleteItem(edge.id());
            // get the new region the edge is in
            // (since the edge is no any more an active edge)
            //std::cout<<"get the new node  \n";
            const Node newNode = mergeGraph_.inactiveEdgesNode(edge);
            //std::cout<<"new node "<<mergeGraph_.id(newNode)<<"\n";


            // iterate over all edges of this node
            for (IncEdgeIt e(mergeGraph_,newNode);e!=lemon::INVALID;++e){

                //std::cout<<"get inc edge\n";
                const Edge incEdge(*e);

                //std::cout<<"get inc graph edge\n";
                const BaseGraphEdge incGraphEdge = EdgeHelper::itemToGraphItem(mergeGraph_,incEdge);

                //std::cout<<"get inc edge weight"<<counter<<"\n";
                // compute the new weight for this edge
                // (this should involve region differences)
                const ValueType newWeight =    getEdgeWeight(incEdge);


                // change the weight in pq by repushing

                //std::cout<<"push\n";
                pq_.push(incEdge.id(),newWeight);
                minWeightEdgeMap_[incGraphEdge]=newWeight;

            }
            //std::cout<<"done\n";
        }

        /// \brief get the edge which should be contracted next
        Edge contractionEdge(){
            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                eraseEdge(Edge(minLabel));
                minLabel = pq_.top();
            }
            return Edge(minLabel);
        }

        /// \brief get the edge weight of the edge which should be contracted next
        WeightType contractionWeight(){
            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                eraseEdge(Edge(minLabel));
                minLabel = pq_.top();
            }
            return pq_.topPriority();

        }


        /// \brief get a reference to the merge
        MergeGraph & mergeGraph(){
            return mergeGraph_;
        }

        bool done(){

            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                eraseEdge(Edge(minLabel));
                minLabel = pq_.top();
            }
            return mergeGraph_.edgeNum()==0 || mergeGraph_.nodeNum()==1;
        }

    private:
        ValueType getEdgeWeight(const Edge & e){
            
            const Node u = mergeGraph_.u(e);
            const Node v = mergeGraph_.v(e);

            const size_t dU = mergeGraph_.degree(u);
            const size_t dV = mergeGraph_.degree(u);
            const BaseGraphEdge ee=EdgeHelper::itemToGraphItem(mergeGraph_,e);
            const BaseGraphNode uu=NodeHelper::itemToGraphItem(mergeGraph_,u);
            const BaseGraphNode vv=NodeHelper::itemToGraphItem(mergeGraph_,v);

            const float sizeU = nodeSizeMap_[uu] ;
            const float sizeV = nodeSizeMap_[vv] ;

            const ValueType wardFac = 2.0 / ( 1.0/std::pow(sizeU,wardness_) + 1/std::pow(sizeV,wardness_) );
            //const ValueType wardFac = (wardFacRaw*wardness_) + (1.0-wardness_);

            const ValueType fromEdgeIndicator = edgeIndicatorMap_[ee];
            const ValueType totalWeight = fromEdgeIndicator*wardFac;
            return totalWeight;
        }


        MergeGraph & mergeGraph_;
        EDGE_INDICATOR_MAP edgeIndicatorMap_;
        EDGE_SIZE_MAP edgeSizeMap_;
        NODE_SIZE_MAP nodeSizeMap_;
        MIN_WEIGHT_MAP minWeightEdgeMap_;
        vigra::ChangeablePriorityQueue< ValueType > pq_;
        ValueType wardness_;;
    };





















    /// \brief  This Cluster Operator is a MONSTER.
    /// It can to really a lot.
    ///
    /// Each edge has a single scalar weight w_e.
    /// Each node has a feature vector f_n.
    /// (all f_n's have the same length).
    /// Edges and nodes have a length / size
    /// 
    /// The total edge weight is computed via a complicated formula
    /// 
    /// The main idea is the following.
    /// Use a  mixture between the edge weights w_e,
    /// and node based edge weights which are computed
    /// via a metric which measures the 'difference' between
    /// the u/v feature vectors f_n.
    ///
    /// Furthermore a 'Ward'-like regularization can be applied.
    /// This is useful if one  have clusters with sizes
    /// in the same magnitude (or 'similar' sizes).
    /// The amount of 'ward'-regularization is controlled
    /// with the 'wardness' parameter.
    ///
    /// Also labels (in the sense of seeds) can be attached to get a 'watershed-ish'
    /// behavior (nodes with different labels will never be merged)
    /// The '0'-Label is used to indicate that there is no label at all.
    /// If certain connected regions share the same seed/label
    /// it is not guaranteed that they will merge. But a certain prior / multiplier
    /// must be specified. The total weight of an edge where the u/v node have
    /// the same label is multiplied with this very multiplier.
    template<
        class MERGE_GRAPH,
        class EDGE_INDICATOR_MAP,
        class EDGE_SIZE_MAP,
        class NODE_FEATURE_MAP,
        class NODE_SIZE_MAP,
        class MIN_WEIGHT_MAP,
        class NODE_LABEL_MAP
    >
    class EdgeWeightNodeFeatures{
        
        typedef EdgeWeightNodeFeatures<
            MERGE_GRAPH,
            EDGE_INDICATOR_MAP,
            EDGE_SIZE_MAP,
            NODE_FEATURE_MAP,
            NODE_SIZE_MAP,
            MIN_WEIGHT_MAP,
            NODE_LABEL_MAP
        > SelfType;
    public:


        typedef typename EDGE_INDICATOR_MAP::Value ValueType;
        typedef ValueType WeightType;
        typedef MERGE_GRAPH MergeGraph;
        typedef typename MergeGraph::Graph Graph;
        typedef typename Graph::Edge BaseGraphEdge;
        typedef typename Graph::Node BaseGraphNode;
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
        /// \brief construct cluster operator
        EdgeWeightNodeFeatures(
            MergeGraph & mergeGraph,
            EDGE_INDICATOR_MAP edgeIndicatorMap,
            EDGE_SIZE_MAP edgeSizeMap,
            NODE_FEATURE_MAP nodeFeatureMap,
            NODE_SIZE_MAP nodeSizeMap,
            MIN_WEIGHT_MAP minWeightEdgeMap,
            NODE_LABEL_MAP nodeLabelMap,
            const ValueType beta,
            const metrics::MetricType metricType,
            const ValueType wardness=1.0,
            const ValueType gamma = 10000000.0,
            const ValueType sameLabelMultiplier = 0.8
        )
        :   mergeGraph_(mergeGraph),
            edgeIndicatorMap_(edgeIndicatorMap),
            edgeSizeMap_(edgeSizeMap),
            nodeFeatureMap_(nodeFeatureMap),
            nodeSizeMap_(nodeSizeMap),
            minWeightEdgeMap_(minWeightEdgeMap),
            nodeLabelMap_(nodeLabelMap),
            pq_(mergeGraph.maxEdgeId()+1),
            beta_(beta),
            wardness_(wardness),
            gamma_(gamma),
            sameLabelMultiplier_(sameLabelMultiplier),
            metric_(metricType)
        {
            typedef typename MergeGraph::MergeNodeCallBackType MergeNodeCallBackType;
            typedef typename MergeGraph::MergeEdgeCallBackType MergeEdgeCallBackType;
            typedef typename MergeGraph::EraseEdgeCallBackType EraseEdgeCallBackType;


            MergeNodeCallBackType cbMn(MergeNodeCallBackType:: template from_method<SelfType,&SelfType::mergeNodes>(this));
            MergeEdgeCallBackType cbMe(MergeEdgeCallBackType:: template from_method<SelfType,&SelfType::mergeEdges>(this));
            EraseEdgeCallBackType cbEe(EraseEdgeCallBackType:: template from_method<SelfType,&SelfType::eraseEdge>(this));

            mergeGraph_.registerMergeNodeCallBack(cbMn);
            mergeGraph_.registerMergeEdgeCallBack(cbMe);
            mergeGraph_.registerEraseEdgeCallBack(cbEe);

          

            for(EdgeIt e(mergeGraph);e!=lemon::INVALID;++e){
                const Edge edge = *e;
                const BaseGraphEdge graphEdge=EdgeHelper::itemToGraphItem(mergeGraph_,edge);
                const index_type edgeId = mergeGraph_.id(edge);
                const ValueType currentWeight = this->getEdgeWeight(edge);
                pq_.push(edgeId,currentWeight);
                minWeightEdgeMap_[graphEdge]=currentWeight;
            }

        }

        /// \brief will be called via callbacks from mergegraph
        void mergeEdges(const Edge & a,const Edge & b){
            // update features / weigts etc
            const BaseGraphEdge aa=EdgeHelper::itemToGraphItem(mergeGraph_,a);
            const BaseGraphEdge bb=EdgeHelper::itemToGraphItem(mergeGraph_,b);
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

        /// \brief will be called via callbacks from mergegraph
        void mergeNodes(const Node & a,const Node & b){
            const BaseGraphNode aa=NodeHelper::itemToGraphItem(mergeGraph_,a);
            const BaseGraphNode bb=NodeHelper::itemToGraphItem(mergeGraph_,b);
            NodeFeatureReference va=nodeFeatureMap_[aa];
            NodeFeatureReference vb=nodeFeatureMap_[bb];
            va*=nodeSizeMap_[aa];
            vb*=nodeSizeMap_[bb];
            va+=vb;
            nodeSizeMap_[aa]+=nodeSizeMap_[bb];
            va/=(nodeSizeMap_[aa]);
            vb/=nodeSizeMap_[bb];


            // update labels
            const UInt32 labelA = nodeLabelMap_[aa];
            const UInt32 labelB = nodeLabelMap_[bb];

            if(labelA!=0 && labelB!=0 && labelA!=labelB){
                throw std::runtime_error("both nodes have labels");
            }
            else{
                const UInt32 newLabel  = std::max(labelA, labelB);
                nodeLabelMap_[aa] = newLabel;
            }
        }

        /// \brief will be called via callbacks from mergegraph
        void eraseEdge(const Edge & edge){

            //std::cout<<"start to erase edge "<<mergeGraph_.id(edge)<<"\n";
            // delete edge from pq
            pq_.deleteItem(edge.id());
            // get the new region the edge is in
            // (since the edge is no any more an active edge)
            //std::cout<<"get the new node  \n";
            const Node newNode = mergeGraph_.inactiveEdgesNode(edge);
            //std::cout<<"new node "<<mergeGraph_.id(newNode)<<"\n";


            // iterate over all edges of this node
            for (IncEdgeIt e(mergeGraph_,newNode);e!=lemon::INVALID;++e){

                //std::cout<<"get inc edge\n";
                const Edge incEdge(*e);

                //std::cout<<"get inc graph edge\n";
                const BaseGraphEdge incGraphEdge = EdgeHelper::itemToGraphItem(mergeGraph_,incEdge);

                //std::cout<<"get inc edge weight"<<counter<<"\n";
                // compute the new weight for this edge
                // (this should involve region differences)
                const ValueType newWeight =    getEdgeWeight(incEdge);


                // change the weight in pq by repushing

                //std::cout<<"push\n";
                pq_.push(incEdge.id(),newWeight);
                minWeightEdgeMap_[incGraphEdge]=newWeight;

            }
            //std::cout<<"done\n";
        }

        /// \brief get the edge which should be contracted next
        Edge contractionEdge(){
            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteItem(minLabel);
                minLabel = pq_.top();
            }
            return Edge(minLabel);
        }

        /// \brief get the edge weight of the edge which should be contracted next
        WeightType contractionWeight(){
            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteItem(minLabel);
                minLabel = pq_.top();
            }
            return pq_.topPriority();

        }


        /// \brief get a reference to the merge
        MergeGraph & mergeGraph(){
            return mergeGraph_;
        }

        bool done(){

            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteItem(minLabel);
                minLabel = pq_.top();
            }
            const ValueType p =  pq_.topPriority();

            return p>= gamma_;
        }

    private:
        ValueType getEdgeWeight(const Edge & e){
            
            const Node u = mergeGraph_.u(e);
            const Node v = mergeGraph_.v(e);

            const size_t dU = mergeGraph_.degree(u);
            const size_t dV = mergeGraph_.degree(u);
            const BaseGraphEdge ee=EdgeHelper::itemToGraphItem(mergeGraph_,e);
            const BaseGraphNode uu=NodeHelper::itemToGraphItem(mergeGraph_,u);
            const BaseGraphNode vv=NodeHelper::itemToGraphItem(mergeGraph_,v);

            const float sizeU = std::min(nodeSizeMap_[uu] , float(std::pow(50.f,3)));
            const float sizeV = std::min(nodeSizeMap_[vv] , float(std::pow(50.f,3)));

            const ValueType wardFacRaw = 1.0 / ( 1.0/std::sqrt(sizeU) + 1.0/std::sqrt(sizeV) );
            const ValueType wardFac = (wardFacRaw*wardness_) + (1.0-wardness_);

            const ValueType fromEdgeIndicator = edgeIndicatorMap_[ee];
            ValueType fromNodeDist = metric_(nodeFeatureMap_[uu],nodeFeatureMap_[vv]);
            ValueType totalWeight = ((1.0-beta_)*fromEdgeIndicator + beta_*fromNodeDist)*wardFac;


            const UInt32 labelA = nodeLabelMap_[uu];
            const UInt32 labelB = nodeLabelMap_[vv];

            if(labelA!=0 && labelB!=0){
                if(labelA == labelB){
                    totalWeight*=sameLabelMultiplier_;
                }
                else{
                    totalWeight += gamma_;
                }
            }
            return totalWeight;
        }


        MergeGraph & mergeGraph_;
        EDGE_INDICATOR_MAP edgeIndicatorMap_;
        EDGE_SIZE_MAP edgeSizeMap_;
        NODE_FEATURE_MAP nodeFeatureMap_;
        NODE_SIZE_MAP nodeSizeMap_;
        MIN_WEIGHT_MAP minWeightEdgeMap_;
        NODE_LABEL_MAP nodeLabelMap_;
        vigra::ChangeablePriorityQueue< ValueType > pq_;
        ValueType beta_;
        ValueType wardness_;
        ValueType gamma_;
        ValueType sameLabelMultiplier_;
        metrics::Metric<float> metric_;
    };



} // end namespace cluster_operators



    /// \brief  do hierarchical clustering with a given cluster operator
    template< class CLUSTER_OPERATOR>
    class HierarchicalClustering{

    public:
        typedef CLUSTER_OPERATOR                        ClusterOperator;
        typedef typename ClusterOperator::MergeGraph    MergeGraph;
        typedef typename MergeGraph::Graph              Graph;
        typedef typename Graph::Edge                    BaseGraphEdge;
        typedef typename Graph::Node                    BaseGraphNode;
        typedef typename MergeGraph::Edge               Edge;
        typedef typename MergeGraph::Node               Node;
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

        /// \brief construct HierarchicalClustering from clusterOperator and an optional parameter object
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
            toTimeStamp_(),
            timeStampIndexToMergeIndex_(),
            mergeTreeEndcoding_()
        {
            if(param_.buildMergeTreeEncoding_){
                // this can be be made smater since user can pass
                // stoping condition based on nodeNum
                mergeTreeEndcoding_.reserve(graph_.nodeNum()*2);
                toTimeStamp_.resize(graph_.maxNodeId()+1);
                timeStampIndexToMergeIndex_.resize(graph_.maxNodeId()+1);
                for(MergeGraphIndexType nodeId=0;nodeId<=mergeGraph_.maxNodeId();++nodeId){
                    toTimeStamp_[nodeId]=nodeId;
                }
            }



        }

        /// \brief start the clustering
        void cluster(){
            if(param_.verbose_)
                std::cout<<"\n"; 
            while(mergeGraph_.nodeNum()>param_.nodeNumStopCond_ && mergeGraph_.edgeNum()>0 && !clusterOperator_.done()){
                
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
                    //std::cout<<"constract\n";
                    // do the merge 
                    mergeGraph_.contractEdge( edgeToRemove );
                }
                if(param_.verbose_ && mergeGraph_.nodeNum()%1==0){
                    std::cout<<"\rNodes: "<<std::setw(10)<<mergeGraph_.nodeNum()<<std::flush;
                }
                
            }
            if(param_.verbose_)
                std::cout<<"\n"; 
        }

        /// \brief get the encoding of the merge tree
        const MergeTreeEncoding & mergeTreeEndcoding()const{
            return mergeTreeEndcoding_;
        }

        template<class EDGE_MAP>
        void ucmTransform(EDGE_MAP & edgeMap)const{
            typedef typename Graph::EdgeIt  BaseGraphEdgeIt;

            for(BaseGraphEdgeIt iter(graph()); iter!=lemon::INVALID; ++iter ){
                const BaseGraphEdge edge=*iter;
                edgeMap[edge] = edgeMap[mergeGraph().reprGraphEdge(edge)];
            }
        }

        /// \brief get the node id's which are the leafes of a treeNodeId
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

        /// \brief get the graph the merge graph is based on
        const Graph & graph()const{
            return graph_;
        }

        /// \brief get the merge graph
        const MergeGraph & mergeGraph()const{
            return mergeGraph_;
        }   

        /// \brief get the representative node id
        const MergeGraphIndexType reprNodeId(const MergeGraphIndexType id)const{
            return mergeGraph_.reprNodeId(id);
        }
    private:

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
