/************************************************************************/
/*                                                                      */
/*     Copyright 2009-2010 by Ullrich Koethe and Hans Meine             */
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
#ifndef VIGRA_MERGE_GRAPH_CLUSTERING_OPERATOR_HXX
#define VIGRA_MERGE_GRAPH_CLUSTERING_OPERATOR_HXX

/*vigra*/
#include <vigra/priority_queue.hxx>
#include <vigra/merge_graph/distance.hxx>

namespace vigra{



    namespace cluster_operators{

    template<
    class MERGE_GRAPH,
    class EDGE_INDICATOR_MAP,
    class EDGE_SIZE_MAP,
    class NODE_FEATURE_MAP,
    class NODE_SIZE_MAP,
    class MIN_WEIGHT_MAP
    >
    class EdgeWeightNodeFeatures{
        
        typedef EdgeWeightNodeFeatures<
            MERGE_GRAPH,
            EDGE_INDICATOR_MAP,
            EDGE_SIZE_MAP,
            NODE_FEATURE_MAP,
            NODE_SIZE_MAP,
            MIN_WEIGHT_MAP
        > SelfType;
    public:

        static const size_t CHI_SQUARED_DISTANCE=0;
        static const size_t NORM_SQUARED_DISTANCE=1;
        static const size_t NORM_DISTANCE=2;


        typedef typename EDGE_INDICATOR_MAP::Value ValueType;
        typedef ValueType WeightType;
        typedef MERGE_GRAPH MergeGraph;
        typedef typename MergeGraph::Graph Graph;
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
        EdgeWeightNodeFeatures(
            MergeGraph & mergeGraph,
            EDGE_INDICATOR_MAP edgeIndicatorMap,
            EDGE_SIZE_MAP edgeSizeMap,
            NODE_FEATURE_MAP nodeFeatureMap,
            NODE_SIZE_MAP nodeSizeMap,
            MIN_WEIGHT_MAP minWeightEdgeMap,
            const ValueType beta,
            const size_t nodeDistType,
            const ValueType wardness=1.0
        )
        :   mergeGraph_(mergeGraph),
            edgeIndicatorMap_(edgeIndicatorMap),
            edgeSizeMap_(edgeSizeMap),
            nodeFeatureMap_(nodeFeatureMap),
            nodeSizeMap_(nodeSizeMap),
            minWeightEdgeMap_(minWeightEdgeMap),
            pq_(mergeGraph.maxEdgeId()+1),
            beta_(beta),
            nodeDistType_(nodeDistType),
            wardness_(wardness)
        {

            mergeGraph_.registerMergeNodeCallBack(*this,& SelfType::mergeNodes);
            mergeGraph_.registerMergeEdgeCallBack(*this,& SelfType::mergeEdges);
            mergeGraph_.registerEraseEdgeCallBack(*this,& SelfType::eraseEdge);


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

            //std::cout<<"start to erase edge "<<mergeGraph_.id(edge)<<"\n";
            // delete edge from pq
            pq_.deleteItem(edge.id());
            // get the new region the edge is in
            // (since the edge is no any more an active edge)
            //std::cout<<"get the new node  \n";
            const Node newNode = mergeGraph_.inactiveEdgesNode(edge);
            //std::cout<<"new node "<<mergeGraph_.id(newNode)<<"\n";

            size_t counter=0;
            // iterate over all edges of this node
            for (IncEdgeIt e(mergeGraph_,newNode);e!=lemon::INVALID;++e){

                //std::cout<<"get inc edge\n";
                const Edge incEdge(*e);

                //std::cout<<"get inc graph edge\n";
                const GraphEdge incGraphEdge = EdgeHelper::itemToGraphItem(mergeGraph_,incEdge);

                //std::cout<<"get inc edge weight"<<counter<<"\n";
                // compute the new weight for this edge
                // (this should involve region differences)
                const ValueType newWeight = getEdgeWeight(incEdge);
                // change the weight in pq by repushing

                //std::cout<<"push\n";
                pq_.push(incEdge.id(),newWeight);
                // remember edge weight

                //std::cout<<"set new\n";
                minWeightEdgeMap_[incGraphEdge]=newWeight;
                ++counter;
            }
            //std::cout<<"done\n";
        }
        Edge contractionEdge(){
            index_type minLabel = pq_.top();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteItem(minLabel);
                index_type minLabel = pq_.top();
            }
            return Edge(minLabel);
        }
        WeightType contractionWeight()const{
            return pq_.topPriority();
        }
        MergeGraph & mergeGraph(){
            return mergeGraph_;
        }
    private:
        ValueType getEdgeWeight(const Edge & e){

            //std::cout<<"raw id ?"<<e.id()<<"\n";
            //std::cout<<"is E Invalid ?"<<bool(e==lemon::INVALID)<<"\n";
            //std::cout<<"GET THE IEDGE WEIGHT for edge<<"<< mergeGraph_.id(e) <<"<<\n";
            //std::cout<<"1\n";
            const Node u = mergeGraph_.u(e);
            const Node v = mergeGraph_.v(e);

            //std::cout<<"2\n";
            const GraphEdge ee=EdgeHelper::itemToGraphItem(mergeGraph_,e);
            const GraphNode uu=NodeHelper::itemToGraphItem(mergeGraph_,u);
            const GraphNode vv=NodeHelper::itemToGraphItem(mergeGraph_,v);

            //std::cout<<"3\n";
            const ValueType wardFacRaw = 1.0 / ( 1.0/std::log(nodeSizeMap_[uu]) + 1.0/std::log(nodeSizeMap_[vv]) );
            const ValueType wardFac = (wardFacRaw*wardness_) + (1.0-wardness_);

            //std::cout<<"4\n";
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
            //std::cout<<"5\n";
            const ValueType totalWeight = ((1.0-beta_)*fromEdgeIndicator + beta_*fromNodeDist)*wardFac;
            return totalWeight;
        }


        MergeGraph & mergeGraph_;
        EDGE_INDICATOR_MAP edgeIndicatorMap_;
        EDGE_SIZE_MAP edgeSizeMap_;
        NODE_FEATURE_MAP nodeFeatureMap_;
        NODE_SIZE_MAP nodeSizeMap_;
        MIN_WEIGHT_MAP minWeightEdgeMap_;
        vigra::ChangeablePriorityQueue< ValueType > pq_;

        ValueType beta_;
        size_t nodeDistType_;
        ValueType wardness_;
    };






    } // end namespace cluster_operators


} // end namespace vigra

#endif