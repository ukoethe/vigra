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
class EdgeWeightedUcm
{
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
        for (IncEdgeIt e(mergeGraph_,newNode);e!=lemon::INVALID;++e)
        {
            //std::cout<<"get inc edge\n";
            const Edge incEdge(*e);

            //std::cout<<"get inc graph edge\n";
            const BaseGraphEdge incGraphEdge = EdgeHelper::itemToGraphItem(mergeGraph_,incEdge);

            //std::cout<<"get inc edge weight"<<counter<<"\n";
            // compute the new weight for this edge
            // (this should involve region differences)
            const ValueType newWeight =    getEdgeWeight(incEdge);

            // change the weight in pq by repushing
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
    /// It can really do a lot.
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
            bool done = false;
            const BaseGraphEdge aa=EdgeHelper::itemToGraphItem(mergeGraph_,a);
            const BaseGraphEdge bb=EdgeHelper::itemToGraphItem(mergeGraph_,b);
            if(!isLifted_.empty()){
                auto isLiftedA =  isLifted_[mergeGraph_.graph().id(aa)];
                auto isLiftedB =  isLifted_[mergeGraph_.graph().id(bb)];
                if(isLiftedA && isLiftedB){
                    pq_.deleteItem(b.id());
                    done = true;
                }
                isLifted_[mergeGraph_.graph().id(aa)] = isLiftedA && isLiftedB;
            }
            if(!done){
                
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
            //std::cout<<"mg e"<<mergeGraph_.edgeNum()<<" mg n"<<mergeGraph_.nodeNum()<<" cw"<< this->contractionWeight()<<"\n";
            if(!isLifted_.empty()){
                if(isLifted_[minLabel])
                    throw std::runtime_error("use lifted edges only if you are DerThorsten or know what you are doing\n");
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


        template<class ITER>
        void setLiftedEdges(ITER idsBegin, ITER idsEnd){
            if(isLifted_.size()<mergeGraph_.graph().maxEdgeId()+1){
                isLifted_.resize(mergeGraph_.graph().maxEdgeId()+1,false);
                std::fill(isLifted_.begin(), isLifted_.end(), false);
            }
            while(idsBegin!=idsEnd){
                isLifted_[*idsBegin] = true;

                const ValueType currentWeight = this->getEdgeWeight(Edge(*idsBegin));
                pq_.push(*idsBegin,currentWeight);
                minWeightEdgeMap_[mergeGraph_.graph().edgeFromId(*idsBegin)]=currentWeight;


                ++idsBegin;
            }
        }
    private:
        ValueType getEdgeWeight(const Edge & e){
            const BaseGraphEdge ee=EdgeHelper::itemToGraphItem(mergeGraph_,e);
            if(!isLifted_.empty() && isLifted_[mergeGraph_.graph().id(ee)]){
                //std::cout<<"found lifted edge\n";
                return 10000000.0;// std::numeric_limits<ValueType>::infinity();
            }

            const Node u = mergeGraph_.u(e);
            const Node v = mergeGraph_.v(e);
           
            const BaseGraphNode uu=NodeHelper::itemToGraphItem(mergeGraph_,u);
            const BaseGraphNode vv=NodeHelper::itemToGraphItem(mergeGraph_,v);

            const float sizeU = nodeSizeMap_[uu];
            const float sizeV = nodeSizeMap_[vv];


            const ValueType wardFac = 2.0 / ( 1.0/std::pow(sizeU,wardness_) + 1.0/std::pow(sizeV,wardness_) );

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

        std::vector<bool> isLifted_;
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
            graph_(mergeGraph_.graph())
        {
        }

        /// \brief start the clustering
        void cluster(){
            if(param_.verbose_)
                std::cout<<"\n";
            while(mergeGraph_.nodeNum()>param_.nodeNumStopCond_ && mergeGraph_.edgeNum()>0 && !clusterOperator_.done()){

                const Edge edgeToRemove = clusterOperator_.contractionEdge();
                           
                // const MergeGraphIndexType uid = mergeGraph_.id(mergeGraph_.u(edgeToRemove));
                // const MergeGraphIndexType vid = mergeGraph_.id(mergeGraph_.v(edgeToRemove));
                
                mergeGraph_.contractEdge( edgeToRemove );


                if(param_.verbose_ && mergeGraph_.nodeNum()%1==0){
                    std::cout<<"\rNodes: "<<std::setw(10)<<mergeGraph_.nodeNum()<<std::flush;
                }

            }
            if(param_.verbose_)
                std::cout<<"\n";
        }


        template<class EDGE_MAP>
        void ucmTransform(EDGE_MAP & edgeMap)const{
            typedef typename Graph::EdgeIt  BaseGraphEdgeIt;

            for(BaseGraphEdgeIt iter(graph()); iter!=lemon::INVALID; ++iter ){
                const BaseGraphEdge edge=*iter;
                edgeMap[edge] = edgeMap[mergeGraph().reprGraphEdge(edge)];
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



        ClusterOperator & clusterOperator_;
        Parameter          param_;
        MergeGraph & mergeGraph_;
        const Graph  & graph_;
        // parameter object



    };


}

#endif // VIGRA_HIERARCHICAL_CLUSTERING_HXX
