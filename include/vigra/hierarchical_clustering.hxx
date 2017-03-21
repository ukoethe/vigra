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
#include "merge_graph_adaptor.hxx"

namespace vigra{

/** \addtogroup GraphDataStructures
*/
//@{

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
            const ValueType wardness=static_cast<ValueType>(1.0),
            const ValueType gamma = static_cast<ValueType>(10000000.0),
            const ValueType sameLabelMultiplier = static_cast<ValueType>(0.8)
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
            metric_(metricType),
            useStopWeight_(false),
            stopWeight_()
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
                const bool isLiftedA =  isLifted_[mergeGraph_.graph().id(aa)];
                const bool isLiftedB =  isLifted_[mergeGraph_.graph().id(bb)];
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
            if(useStopWeight_){
                if(p >= stopWeight_){
                    return true;
                }
            }
            return p>= gamma_;
        }

        template<class ITER>
        void setLiftedEdges(ITER idsBegin, ITER idsEnd){
            if(isLifted_.size()<std::size_t(mergeGraph_.graph().maxEdgeId()+1)){
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

        void enableStopWeight(const ValueType stopWeight){
            useStopWeight_ = true;
            stopWeight_ = stopWeight;
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


            const ValueType wardFac = 2.0 / ( 1.0/std::pow(sizeU,wardness_) + 1/std::pow(sizeV,wardness_) );

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
        bool useStopWeight_;
        ValueType stopWeight_;
    };



} // end namespace cluster_operators

/** \brief  Options object for hierarchical clustering.

  <b>\#include</b> \<vigra/hierarchical_clustering.hxx\><br/>
  Namespace: vigra

  This class allows to set various parameters of \ref hierarchicalClustering().
  See there for usage examples.
*/
class ClusteringOptions
{
  public:

    ClusteringOptions(
        const size_t      nodeNumStopCond = 1,
        const bool        buildMergeTree  = false,
        const bool        verbose         = false)
    : nodeNumStopCond_ (nodeNumStopCond)
    , maxMergeWeight_(NumericTraits<double>::max())
    , nodeFeatureImportance_(0.5)
    , sizeImportance_(1.0)
    , nodeFeatureMetric_(metrics::ManhattanMetric)
    , buildMergeTreeEncoding_(buildMergeTree)
    , verbose_(verbose)
    {}

        /** Stop merging when the number of clusters reaches this threshold.

            Default: 1 (don't stop early)
        */
    ClusteringOptions & minRegionCount(size_t count)
    {
        nodeNumStopCond_ = count;
        return *this;
    }

        /** Stop merging when the weight of the cheapest edge exceeds this threshold.

            Default: infinity (don't stop early)
        */
    ClusteringOptions & maxMergeDistance(double val)
    {
        maxMergeWeight_ = val;
        return *this;
    }

        /** Importance of node features relative to edge weights.

            Must be between 0 and 1, with 0 meaning that node features
            are ignored, and 1 meaning that edge weights are ignored.

            Default: 0.5 (equal importance)
        */
    ClusteringOptions & nodeFeatureImportance(double val)
    {
        vigra_precondition(0.0 <= val && val <= 1.0,
            "ClusteringOptions::nodePropertyImportance(val): 0 <= val <= 1 required.");
        nodeFeatureImportance_ = val;
        return *this;
    }

        /** Importance of node size.

            Must be between 0 and 1, with 0 meaning that size is ignored,
            and 1 meaning that the algorithm prefers to keep cluster sizes
            balanced.

            Default: 1.0 (prefer like-sized clusters)
        */
    ClusteringOptions & sizeImportance(double val)
    {
        vigra_precondition(0.0 <= val && val <= 1.0,
            "ClusteringOptions::sizeImportance(val): 0 <= val <= 1 required.");
        sizeImportance_ = val;
        return *this;
    }


        /** Metric to be used when transforming node features into cluster distances.

            The cluster (= node) distance is the respective norm of the difference
            vector between the corresponding node feature vectors.

            Default: metrics::ManhattanMetric (L1-norm of the feature difference)
        */
    ClusteringOptions & nodeFeatureMetric(metrics::MetricType metric)
    {
        nodeFeatureMetric_ = metric;
        return *this;
    }

    ClusteringOptions & buildMergeTreeEncoding(bool val=true)
    {
        buildMergeTreeEncoding_ = val;
        return *this;
    }

        /** Display progress information.

            Default: false
        */
    ClusteringOptions & verbose(bool val=true)
    {
        verbose_ = val;
        return *this;
    }

    size_t nodeNumStopCond_;
    double maxMergeWeight_;
    double nodeFeatureImportance_;
    double sizeImportance_;
    metrics::MetricType nodeFeatureMetric_;
    bool   buildMergeTreeEncoding_;
    bool   verbose_;
};

// \brief  do hierarchical clustering with a given cluster operator
template< class CLUSTER_OPERATOR>
class HierarchicalClusteringImpl
{
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

    typedef ClusteringOptions Parameter;

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

    /// \brief construct HierarchicalClusteringImpl from clusterOperator and an optional parameter object
    HierarchicalClusteringImpl(
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


/********************************************************/
/*                                                      */
/*                hierarchicalClustering                */
/*                                                      */
/********************************************************/

/** \brief Reduce the number of nodes in a graph by iteratively contracting
    the cheapest edge.

    <b> Declarations:</b>

    \code
    namespace vigra {
        template <class GRAPH,
                  class EDGE_WEIGHT_MAP,  class EDGE_LENGTH_MAP,
                  class NODE_FEATURE_MAP, class NOSE_SIZE_MAP,
                  class NODE_LABEL_MAP>
        void
        hierarchicalClustering(GRAPH const & graph,
                               EDGE_WEIGHT_MAP const & edgeWeights, EDGE_LENGTH_MAP const & edgeLengths,
                               NODE_FEATURE_MAP const & nodeFeatures, NOSE_SIZE_MAP const & nodeSizes,
                               NODE_LABEL_MAP & labelMap,
                               ClusteringOptions options = ClusteringOptions());
    }
    \endcode

    Hierarchical clustering is a simple and versatile image segmentation
    algorithm that typically operates either directly on the pixels (e.g. on
    a \ref vigra::GridGraph) or on a region adjacency graph over suitable
    superpixels (e.g. on an \ref vigra::AdjacencyListGraph). The graph is
    passed to the function in its first argument. After clustering is completed,
    the parameter \a labelMap contains a mapping from original node IDs to
    the ID of the cluster each node belongs to. Cluster IDs correspond to
    the ID of an arbitrarily chosen representative node within each cluster,
    i.e. they form a sparse subset of the original IDs.

    Properties of the graph's edges and nodes are provided in the property maps
    \a edgeWeights, \a edgeLengths, \a nodeFeatures, and \a nodeSizes. These maps
    are indexed by edge or node ID and return the corresponding feature. Features
    must by arithmetic scalars or, in case of node features, scalars or vectors
    of scalars (precisely: objects that provide <tt>begin()</tt> and <tt>end()</tt>
    to create an STL range). Edge weights are typically derived from an edge
    indicator such as the gradient magnitude, and node features are either the
    responses of a filter family (when clustering on the pixel grid), or region
    statistics as computed by \ref FeatureAccumulators (when clustering on
    superpixels).

    In each step, the algorithm merges the two nodes \f$u\f$ and \f$v\f$ whose
    cluster distance is smallest, where the cluster distance is defined as

    \f[
        d_{uv} = \left( (1-\beta) w_{uv} + \beta || f_u - f_v ||_M \right)
                 \cdot \frac{2}{s_u^{-\omega} + s_v^{-\omega}}
    \f]

    with \f$ w_{uv} \f$ denoting the weight of edge \f$uv\f$, \f$f_u\f$ and \f$f_v\f$
    being the node features (possibly vectors to be compared with metric \f$M\f$),
    and \f$s_u\f$ and \f$s_v\f$ the corresponding node sizes. The metric is defined
    in the option object by calling \ref vigra::ClusteringOptions::nodeFeatureMetric()
    and must be selected from the tags defined in \ref vigra::metrics::MetricType.

    The parameters \f$0 \le \beta \le 1\f$ and \f$0 \le \omega \le 1\f$ control the
    relative influence of the inputs: With \f$\beta = 0\f$, the node features are
    ignored, whereas with \f$\beta = 1\f$ the edge weights are ignored. Similarly,
    with \f$\omega = 0\f$, the node size is ignored, whereas with \f$\omega = 1\f$,
    cluster distances are scaled by the harmonic mean of the cluster sizes, making
    the merging of small clusters more favorable. The parameters are defined in the
    option object by calling \ref vigra::ClusteringOptions::nodeFeatureImportance() and
    \ref vigra::ClusteringOptions::sizeImportance() respectively.

    After each merging step, the features of the resulting cluster \f$z\f$ and the weights
    of its outgoing edges are updated by mean of the corresponding properties of the original
    clusters \f$u\f$ and \f$v\f$, weighted by the respective node sizes \f$s_z\f$ and
    edge lengths \f$l_{zy}\f$:

    \f{eqnarray*}{
        s_z & = & s_u + s_v \\
        f_z & = & \frac{s_u f_u + s_v f_v}{s_z} \\
        l_{zy} & = & l_{uy} + l_{vy} \textrm{ for all nodes }y\textrm{ connected to }u\textrm{ or }v \\
        w_{zy} & = & \frac{l_{uy} w_{uy} + l_{vy} w_{vy}}{l_{zy}}
    \f}

    Clustering normally stops when only one cluster remains. This default can be overridden
    by the option object parameters \ref vigra::ClusteringOptions::minRegionCount()
    and \ref vigra::ClusteringOptions::maxMergeDistance() to stop at a particular number of
    clusters or a particular cluster distance respectively.

    <b> Usage:</b>

    <b>\#include</b> \<vigra/hierarchical_clustering.hxx\><br>
    Namespace: vigra

    A fully worked example can be found in <a href="graph_agglomerative_clustering_8cxx-example.html">graph_agglomerative_clustering.cxx</a>
*/
doxygen_overloaded_function(template <...> void hierarchicalClustering)

template <class GRAPH,
          class EDGE_WEIGHT_MAP,  class EDGE_LENGTH_MAP,
          class NODE_FEATURE_MAP, class NOSE_SIZE_MAP,
          class NODE_LABEL_MAP>
void
hierarchicalClustering(GRAPH const & graph,
                       EDGE_WEIGHT_MAP const & edgeWeights, EDGE_LENGTH_MAP const & edgeLengths,
                       NODE_FEATURE_MAP const & nodeFeatures, NOSE_SIZE_MAP const & nodeSizes,
                       NODE_LABEL_MAP & labelMap,
                       ClusteringOptions options = ClusteringOptions())
{
    typedef typename NODE_LABEL_MAP::Value LabelType;
    typedef MergeGraphAdaptor<GRAPH> MergeGraph;
    typedef typename GRAPH::template EdgeMap<float>     EdgeUltrametric;
    typedef typename GRAPH::template NodeMap<LabelType> NodeSeeds;

    MergeGraph mergeGraph(graph);

    // create property maps to store the computed ultrametric and
    // to provide optional cannot-link constraints;
    // we don't use these options here and therefore leave the maps empty
    EdgeUltrametric edgeUltrametric(graph);
    NodeSeeds nodeSeeds(graph);

    // create an operator that stores all property maps needed for
    // hierarchical clustering and updates them after every merge step
    typedef cluster_operators::EdgeWeightNodeFeatures<
        MergeGraph,
        EDGE_WEIGHT_MAP,
        EDGE_LENGTH_MAP,
        NODE_FEATURE_MAP,
        NOSE_SIZE_MAP,
        EdgeUltrametric,
        NodeSeeds>
    MergeOperator;

    MergeOperator mergeOperator(mergeGraph,
                                edgeWeights, edgeLengths,
                                nodeFeatures, nodeSizes,
                                edgeUltrametric, nodeSeeds,
                                options.nodeFeatureImportance_,
                                options.nodeFeatureMetric_,
                                options.sizeImportance_,
                                options.maxMergeWeight_);

    typedef HierarchicalClusteringImpl<MergeOperator> Clustering;

    Clustering clustering(mergeOperator, options);
    clustering.cluster();

    for(typename GRAPH::NodeIt node(graph); node != lemon::INVALID; ++node)
    {
        labelMap[*node] = mergeGraph.reprNodeId(graph.id(*node));
    }
}

//@}

} // namespace vigra

#endif // VIGRA_HIERARCHICAL_CLUSTERING_HXX
