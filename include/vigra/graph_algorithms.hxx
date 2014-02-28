/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier Ullrich Koeth                   */
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

/**
 * This header provides definitions of graph-related algorithms
 */

#ifndef VIGRA_GRAPH_MAP_ALGORITHMS_HXX
#define VIGRA_GRAPH_MAP_ALGORITHMS_HXX

/*std*/
#include <algorithm>
#include <vector>
#include <functional>

/*boost*/
 #include <boost/iterator/transform_iterator.hpp>

/*vigra*/
#include <vigra/graphs.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/priority_queue.hxx>
#include <vigra/union_find.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_helper/on_the_fly_edge_map.hxx>

#define LEMON_UNDIRECTED_GRAPH_TYPEDEFS(GRAPH_CLS,PREFIX,POSTFIX) \
    typedef typename GRAPH_CLS::Edge        PREFIX##Edge      ## POSTFIX; \
    typedef typename GRAPH_CLS::Node        PREFIX##Node      ## POSTFIX; \
    typedef typename GRAPH_CLS::Arc         PREFIX##Arc       ## POSTFIX; \
    typedef typename GRAPH_CLS::EdgeIt      PREFIX##EdgeIt##POSTFIX; \
    typedef typename GRAPH_CLS::NodeIt      PREFIX##NodeIt    ## POSTFIX; \
    typedef typename GRAPH_CLS::ArcIt       PREFIX##ArcIt     ## POSTFIX; \
    typedef typename GRAPH_CLS::OutArcIt    PREFIX##OutArcIt  ## POSTFIX; \
    typedef typename GRAPH_CLS::InArcIt     PREFIX##InArcIt   ## POSTFIX


namespace vigra{




    namespace detail_graph_algorithms{
        template <class GRAPH_MAP,class COMPERATOR>
        struct GraphItemCompare
        {
            
            GraphItemCompare(const GRAPH_MAP & map,const COMPERATOR & comperator)
            : map_(map),
              comperator_(comperator){

            }

            template<class KEY>
            bool operator()(const KEY & a, const KEY & b) const{
                return comperator_(map_[a],map_[b]);
            }

            const GRAPH_MAP & map_;
            const COMPERATOR & comperator_;
        };
    }

    template<class GRAPH,class WEIGHTS,class COMPERATOR>
    void edgeSort(
        const GRAPH   & g,
        const WEIGHTS & weights,
        const COMPERATOR  & comperator,
        std::vector<typename GRAPH::Edge> & sortedEdges
    ){
        sortedEdges.resize(g.edgeNum());
        size_t c=0;
        for(typename GRAPH::EdgeIt e(g);e!=lemon::INVALID;++e){
            sortedEdges[c]=*e;
            ++c;
        }
        detail_graph_algorithms::GraphItemCompare<WEIGHTS,COMPERATOR> edgeComperator(weights,comperator);
        std::sort(sortedEdges.begin(),sortedEdges.end(),edgeComperator);
    }


    template<class GRAPH,class MAP>
    class EdgeMapIteratorHelper{
        public:
            typedef typename GraphMapTypeTraits<MAP>::Reference      Reference;
            typedef typename GraphMapTypeTraits<MAP>::ConstReference ConstReference;
            typedef typename GraphMapTypeTraits<MAP>::Value          Value;
        private:
            struct Transform{
                
                Transform(MAP & map)
                : map_(&map){
                }
                template<class ITEM>
                Reference operator()(const ITEM & item)const{
                    return  map_->operator[](item);
                }
                mutable MAP * map_;
            };  
            struct ConstTransform{
                
                ConstTransform(const MAP & map)
                : map_(&map){
                }
                template<class ITEM>
                ConstReference operator()(const ITEM & item)const{
                    return  map_->operator[](item);
                }
                const MAP * map_;
            }; 
        public:
            typedef  boost::transform_iterator< Transform,      typename GRAPH::EdgeIt,Reference      ,Value> iterator;
            typedef  boost::transform_iterator< ConstTransform, typename GRAPH::EdgeIt,ConstReference ,Value> const_iterator;
            static iterator
            begin(const GRAPH & g, MAP & map){
                Transform f(map);
                typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesBegin(g);
                return iterator(iter,f);
            }
            static iterator
            end(const GRAPH & g, MAP & map){
                Transform f(map);
                typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesEnd(g);
                return iterator(iter,f);
            }
            static const_iterator
            begin(const GRAPH & g, const MAP & map){
                ConstTransform f(map);
                typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesBegin(g);
                return const_iterator(iter,f);
            }
            static const_iterator
            end(const GRAPH & g,const MAP & map){
                ConstTransform f(map);
                typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesEnd(g);
                return const_iterator(iter,f);
            }
        private:
    };

    template<class GRAPH,class MAP>
    class NodeMapIteratorHelper{
        public:
            typedef typename GraphMapTypeTraits<MAP>::Reference      Reference;
            typedef typename GraphMapTypeTraits<MAP>::ConstReference ConstReference;
            typedef typename GraphMapTypeTraits<MAP>::Value          Value;
        private:
            struct Transform{
                
                Transform(MAP & map)
                : map_(&map){
                }
                template<class ITEM>
                Reference operator()(const ITEM & item)const{
                    return  map_->operator[](item);
                }
                mutable MAP * map_;
            };  
            struct ConstTransform{
                
                ConstTransform(const MAP & map)
                : map_(&map){
                }
                template<class ITEM>
                ConstReference operator()(const ITEM & item)const{
                    return  map_->operator[](item);
                }
                const MAP * map_;
            }; 
        public:
            typedef  boost::transform_iterator< Transform,      typename GRAPH::NodeIt,Reference      ,Value> iterator;
            typedef  boost::transform_iterator< ConstTransform, typename GRAPH::NodeIt,ConstReference ,Value> const_iterator;
            static iterator
            begin(const GRAPH & g, MAP & map){
                Transform f(map);
                typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesBegin(g);
                return iterator(iter,f);
            }
            static iterator
            end(const GRAPH & g, MAP & map){
                Transform f(map);
                typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesEnd(g);
                return iterator(iter,f);
            }
            static const_iterator
            begin(const GRAPH & g, const MAP & map){
                ConstTransform f(map);
                typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesBegin(g);
                return const_iterator(iter,f);
            }
            static const_iterator
            end(const GRAPH & g,const MAP & map){
                ConstTransform f(map);
                typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesEnd(g);
                return const_iterator(iter,f);
            }
        private:
    };


    template<class G,class A,class B>
    void copyNodeMap(const G & g,const A & a ,B & b){
        std::copy(NodeMapIteratorHelper<G,A>::begin(g,a),NodeMapIteratorHelper<G,A>::end(g,a), NodeMapIteratorHelper<G,B>::begin(g,b));
    }
    template<class G,class A,class B>
    void copyEdgeMap(const G & g,const A & a ,B & b){
        std::copy(EdgeMapIteratorHelper<G,A>::begin(g,a),EdgeMapIteratorHelper<G,A>::end(g,a), EdgeMapIteratorHelper<G,B>::begin(g,b));
    }
    template<class G,class A,class T>
    void fillNodeMap(const G & g, A & a ,const T & value){
        std::fill(NodeMapIteratorHelper<G,A>::begin(g,a),NodeMapIteratorHelper<G,A>::end(g,a), value);
    }
    template<class G,class A,class T>
    void fillEdgeMap(const G & g,A & a ,const T & value){
        std::fill(EdgeMapIteratorHelper<G,A>::begin(g,a),EdgeMapIteratorHelper<G,A>::end(g,a), value);
    }


    template<
        class GRAPH_IN,
        class GRAPH_IN_NODE_LABEL_MAP
    >
    void makeRegionAdjacencyGraph(
        GRAPH_IN                   graphIn,
        GRAPH_IN_NODE_LABEL_MAP    labels,
        AdjacencyListGraph & rag,
        typename AdjacencyListGraph:: template EdgeMap< std::vector<typename GRAPH_IN::Edge> > & affiliatedEdges,
        const Int64   ignoreLabel=-1
    ){
        rag=AdjacencyListGraph();
        typedef GRAPH_IN_NODE_LABEL_MAP LabelMap;
        typedef typename GraphMapTypeTraits<GRAPH_IN_NODE_LABEL_MAP>::Value LabelType;
        typedef GRAPH_IN GraphIn;
        typedef AdjacencyListGraph GraphOut;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(GraphIn, , GraphIn);
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(GraphOut, ,GraphOut);
        // iterate over all labels in the node map to find min max
        typedef NodeMapIteratorHelper<GraphIn,LabelMap> NodeIterHelper;
        
        for(NodeItGraphIn iter(graphIn);iter!=lemon::INVALID;++iter){
            const LabelType l=labels[*iter];
            if(ignoreLabel==-1 || static_cast<Int64>(l)!=ignoreLabel)
                rag.addNode(l);
        }

        // add al edges
        for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
            const EdgeGraphIn edge(*e);
            const LabelType lu = labels[graphIn.u(edge)];
            const LabelType lv = labels[graphIn.v(edge)];
            if(  lu!=lv && ( ignoreLabel==-1 || (static_cast<Int64>(lu)!=ignoreLabel  && static_cast<Int64>(lv)!=ignoreLabel) )  ){
                // if there is an edge between lu and lv no new edge will be added
                rag.addEdge( rag.nodeFromId(lu),rag.nodeFromId(lv));
            }
        }
        // SET UP HYPEREDGES
        affiliatedEdges.assign(rag);
        // add edges
        for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
            const EdgeGraphIn edge(*e);
            const LabelType lu = labels[graphIn.u(edge)];
            const LabelType lv = labels[graphIn.v(edge)];
            //std::cout<<"edge between ?? "<<lu<<" "<<lv<<"\n";
            if(  lu!=lv && ( ignoreLabel==-1 || (static_cast<Int64>(lu)!=ignoreLabel  && static_cast<Int64>(lv)!=ignoreLabel) )  ){
                //std::cout<<"find edge between "<<lu<<" "<<lv<<"\n";
                EdgeGraphOut ragEdge= rag.findEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
                //std::cout<<"invalid?"<<bool(ragEdge==lemon::INVALID)<<" id "<<rag.id(ragEdge)<<"\n";
                affiliatedEdges[ragEdge].push_back(edge);
                //std::cout<<"write done\n";
            }
        }
    }


    template<class GRAPH,class WEIGHT_TYPE>
    class ShortestPathDijkstra{
    public:
        typedef GRAPH Graph;

        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );

        typedef WEIGHT_TYPE WeightType;
        typedef ChangeablePriorityQueue<WeightType>           PqType;
        typedef typename Graph:: template NodeMap<Node>       PredecessorsMap;
        typedef typename Graph:: template NodeMap<WeightType> DistanceMap;
            

        ShortestPathDijkstra(const Graph & g)
        :   graph_(g),
            pq_(g.maxNodeId()+1),
            predMap_(g),
            distMap_(g)
        {
        }

        template<class WEIGHTS>
        void run(const WEIGHTS & weights,const Node & source,const Node & target = lemon::INVALID){
            source_=source;
            target_=target;

            this->initializeMaps();

            bool finished=false;
            while(!pq_.empty() ){ //&& !finished){
                const Node topNode(graph_.nodeFromId(pq_.top()));
                if(topNode==target_){
                    finished=true;
                    break;
                }
                if(predMap_[topNode]!=lemon::INVALID || topNode==source_ ){
                    pq_.pop();
                    // loop over all neigbours
                    for(OutArcIt outArcIt(graph_,topNode);outArcIt!=lemon::INVALID;++outArcIt){
                        const Node otherNode = graph_.target(*outArcIt);
                        const size_t otherNodeId = graph_.id(otherNode);
    
                        if(pq_.contains(otherNodeId)){
                            const Edge edge(*outArcIt);
                            const WeightType currentDist     = distMap_[otherNode];
                            const WeightType alternativeDist = distMap_[topNode]+weights[edge];
                            if(alternativeDist<currentDist){
                                pq_.push(otherNodeId,alternativeDist);
                                distMap_[otherNode]=alternativeDist;
                                predMap_[otherNode]=topNode;
                            }
                        }

                    }
                }
                else{
                    finished=true;
                    break;
                }
                if(finished){
                    break;
                }
            }
        }

        const Graph & graph()const{
            return graph_;
        }

        const Node & source()const{
            return source_;
        }
        const Node & target()const{
            return target_;
        }

        bool hasTarget()const{
            return target_!=lemon::INVALID;
        }

        const PredecessorsMap & predecessors()const{
            return predMap_;
        }
        const DistanceMap & distances()const{
            return distMap_;
        }

        WeightType distance(const Node & target)const{
            return distMap_[target];
        }


    private:

        void initializeMaps(){
            for(NodeIt n(graph_);n!=lemon::INVALID;++n){
                const Node node(*n);
                pq_.push(graph_.id(node),std::numeric_limits<WeightType>::infinity() );
                distMap_[node]=std::numeric_limits<WeightType>::infinity();
                predMap_[node]=lemon::INVALID;
            }
            distMap_[source_]=static_cast<WeightType>(0.0);
            pq_.push(graph_.id(source_),0.0);
        }

        const Graph  & graph_;
        PqType  pq_;
        PredecessorsMap predMap_;
        DistanceMap     distMap_;

        Node source_;
        Node target_;
    };


    template<class NODE,class PREDECESSORS>
    size_t pathLength(
        const NODE source,
        const NODE target,
        const PREDECESSORS & predecessors
    ){
        if(predecessors[target]==lemon::INVALID)
            return 0;
        else{
            NODE currentNode = target;
            size_t length=1;
            while(currentNode!=source){
                currentNode=predecessors[currentNode];
                length+=1;
            }
            return length;
        }
    }



    template<class GRAPH,class T>
    struct ZeroHeuristc{
        typedef typename GRAPH::Node Node;
        T operator()(const Node & a ,const Node & b )const{
            return static_cast<T>(0.0);
        }
    };


    template<class GRAPH,class WEIGHTS,class PREDECESSORS,class DISTANCE,class HEURSTIC>
    void shortestPathAStar(
        const GRAPH         &           graph,
        const typename GRAPH::Node &    source,
        const typename GRAPH::Node &    target,
        const WEIGHTS       &           weights,
        PREDECESSORS        &           predecessors,
        DISTANCE            &           distance,
        const HEURSTIC      &           heuristic
    ){

        typedef GRAPH                       Graph;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );

        typedef typename WEIGHTS::value_type     WeightType;
        typedef typename DISTANCE::value_type    DistanceType;

        typename  GRAPH:: template NodeMap<bool> closedSet(graph);
        vigra::ChangeablePriorityQueue<typename WEIGHTS::value_type> estimatedDistanceOpenSet(graph.maxNodeId()+1);
        // initialize
        for(NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            closedSet[node]=false;
            distance[node]=std::numeric_limits<DistanceType>::infinity();
            predecessors[node]=lemon::INVALID;
        }
        // distance and estimated distance for start node
        distance[source]=static_cast<DistanceType>(0.0);
        estimatedDistanceOpenSet.push(graph.id(source),heuristic(source,target));

        // while any nodes left in openSet
        while(!estimatedDistanceOpenSet.empty()){

            // get the node with the lpwest estimated distance in the open set
            const Node current = graph.nodeFromId(estimatedDistanceOpenSet.top());

            // reached target?
            if(current==target)
                break;

            // remove current from openSet
            // add current to closedSet
            estimatedDistanceOpenSet.pop();
            closedSet[current]=true;

            // iterate over neigbours of current
            for(OutArcIt outArcIt(graph,current);outArcIt!=lemon::INVALID;++outArcIt){

                // get neigbour node and id
                const Node neighbour = graph.target(*outArcIt);
                const size_t neighbourId = graph.id(neighbour);

                // if neighbour is not yet in closedSet
                if(!closedSet[neighbour]){

                    // get edge between current and neigbour
                    const Edge edge(*outArcIt);

                    // get tentative score
                    const DistanceType tenativeScore = distance[current] + weights[edge];

                    // neighbour NOT in openSet OR tentative score better than the current distance
                    if(!estimatedDistanceOpenSet.contains(neighbourId) || tenativeScore < distance[neighbour] ){
                        // set predecessors and distance
                        predecessors[neighbour]=current;
                        distance[neighbour]=tenativeScore;

                        // update the estimated cost from neighbour to target
                        // ( and neigbour will be (re)-added to openSet)
                        estimatedDistanceOpenSet.push(neighbourId,distance[neighbour]+heuristic(neighbour,target));
                    }
                }
            }
        }
    }
    


    template<class GRAPH,class WEIGHTS,class NODE_LABEL_MAP>
    void minimumSpanningTreeSegmentation(
        const GRAPH     & graph,
        WEIGHTS         & weights,
        const bool useWeightThreshold,
        const bool useNodeThreshold,
        typename WEIGHTS::Value weightThreshold,
        const  size_t nodeNumThreshold,
        NODE_LABEL_MAP     &  nodeLabeling
    ){  

        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );
        typedef typename WEIGHTS::Value WeightType;
        typedef EdgeMapIteratorHelper<GRAPH,WEIGHTS>     WeightIterHelper;
        typedef detail::Partition<size_t> UfdType;
    
        // sort the edges by their weights
        std::vector<Edge> sortedEdges;
        std::less<WeightType> comperator;
        edgeSort(graph,weights,comperator,sortedEdges);

        UfdType ufd(graph.maxNodeId()+1);
        const size_t unusedSets = ufd.numberOfSets()-graph.nodeNum();

        for(size_t i=0;i<sortedEdges.size();++i){
            const Edge e=sortedEdges[i];
            const size_t uId=graph.id(graph.u(e));
            const size_t vId=graph.id(graph.v(e));

            if(useWeightThreshold && weightThreshold <= weights[e])
                break;
            if(useNodeThreshold && (ufd.numberOfSets()-unusedSets)<=nodeNumThreshold )
                break; 
            if(ufd.find(uId)!=ufd.find(vId)){
                ufd.merge(uId,vId);
            }
        }
        for(typename  GRAPH::NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            nodeLabeling[node]=ufd.find(graph.id(node));
        }
    }  

    namespace detail_watersheds_segmentation{

    struct IdentityFunctor{
        template<class L, class T>
        T operator()(const L label,const T  priority)const{
            return priority;
        }
    };

    template<class PRIORITY_TYPE,class LABEL_TYPE>
    struct CarvingFunctor{
        CarvingFunctor(const LABEL_TYPE backgroundLabel,const PRIORITY_TYPE & factor)
        :   backgroundLabel_(backgroundLabel),
            factor_(factor){
        }
        PRIORITY_TYPE operator()(const LABEL_TYPE label,const PRIORITY_TYPE  priority)const{
            return (label==backgroundLabel_ ? priority*factor_ : priority);
        }
        LABEL_TYPE     backgroundLabel_;
        PRIORITY_TYPE  factor_;
    };


    template<class GRAPH,class EDGE_WEIGHTS,class SEEDS,class PRIORITY_MANIP_FUNCTOR,class LABELS>
    void watershedsSegmentationImpl(
        const GRAPH & g,
        const EDGE_WEIGHTS      & edgeWeights,
        const SEEDS             & seeds,
        PRIORITY_MANIP_FUNCTOR  & priorManipFunctor,
        LABELS                  & labels
    ){  
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );
        typedef typename EDGE_WEIGHTS::Value WeightType;
        typedef typename SEEDS::Value   SeedType;
        typedef typename LABELS::Value  LabelType;
        typedef typename Graph:: template NodeMap<bool>    NodeBoolMap;
        typedef NodeMapIteratorHelper<Graph,NodeBoolMap  > InPQHelper;
        typedef NodeMapIteratorHelper<Graph,SEEDS  >       SeedsHelper;
        typedef NodeMapIteratorHelper<Graph,LABELS >       LabelsHelper;
        typedef PriorityQueue<Node,WeightType,true> PQ;

        PQ pq;
        NodeBoolMap inPQ(g);
        copyNodeMap(g,seeds,labels);
        fillNodeMap(g,inPQ,false);

        bool anySeed=false;
        for(NodeIt n(g);n!=lemon::INVALID;++n){
            const Node node(*n);
            if(labels[node]!=static_cast<LabelType>(0)){
                anySeed=true;
                for(OutArcIt a(g,node);a!=lemon::INVALID;++a){
                    const Edge edge(*a);
                    const Node neigbour=g.target(*a);
                    //std::cout<<"n- node "<<g.id(neigbour)<<"\n";
                    if(labels[neigbour]==static_cast<LabelType>(0) && !inPQ[neigbour]){
                        const WeightType priority = priorManipFunctor(labels[node],edgeWeights[edge]);
                        pq.push(neigbour,priority);
                        inPQ[neigbour]=true;
                    }
                }
            }
        }


        if(anySeed){

            while(!pq.empty()){
                const Node node       = pq.top();
                const LabelType label = labels[node]; 
                //std::cout<<"node "<<g.id(node)<<" with label "<<label<<"\n";
                if(label!=0){
                    throw std::runtime_error("this should not happen 0");
                }

                pq.pop();
                bool moreThanOneLabel = false;
                LabelType labelFound  = 0 ;
                for(OutArcIt a(g,node);a!=lemon::INVALID;++a){
                    const Edge edge(*a);
                    const Node neigbour=g.target(*a);
                    if(labels[neigbour]!=static_cast<LabelType>(0)){
                        if(labelFound==0){
                            labelFound=labels[neigbour];
                        }
                        else{
                            moreThanOneLabel=true;
                            break;
                        }
                    }
                }
                if(labelFound!=0 && !moreThanOneLabel ){
                    labels[node]=labelFound;
                    for(OutArcIt a(g,node);a!=lemon::INVALID;++a){
                        const Edge edge(*a);
                        const Node neigbour=g.target(*a);
                        if(labels[neigbour]==static_cast<LabelType>(0)){
                            if(!inPQ[neigbour]){
                                const WeightType priority = priorManipFunctor(labelFound,edgeWeights[edge]);
                                pq.push(neigbour,edgeWeights[edge]);
                                inPQ[neigbour]=true;
                            }
                        }
                    }
                }
            }


            for(NodeIt n(g);n!=lemon::INVALID;++n){
                const Node node(*n);
                if(labels[node]==static_cast<LabelType>(0)){

                    WeightType minWeight       = std::numeric_limits<WeightType>::infinity();
                    LabelType  minWeightLabel  = static_cast<LabelType>(0);
                    for(OutArcIt a(g,node);a!=lemon::INVALID;++a){
                        const Edge edge(*a);
                        const Node neigbour=g.target(*a);
                        const WeightType priority = priorManipFunctor(labels[neigbour],edgeWeights[edge]);
                        if(labels[neigbour]!=0 && priority<minWeight){
                            minWeight=priority;
                            minWeightLabel=labels[neigbour];
                        }
                    }
                    if(minWeightLabel==0){
                        throw std::runtime_error("this should not happen 2");
                    }
                    labels[node]=minWeightLabel;
                }
            }
        }
    }

    } // end namespace detail_watersheds_segmentation

    template<class GRAPH,class EDGE_WEIGHTS,class SEEDS,class LABELS>
    void watershedsSegmentation(
        const GRAPH & g,
        const EDGE_WEIGHTS & edgeWeights,
        const SEEDS        & seeds,
        LABELS             & labels
    ){  
        detail_watersheds_segmentation::IdentityFunctor f;
        detail_watersheds_segmentation::watershedsSegmentationImpl(g,edgeWeights,seeds,f,labels);
    }   
    

    template<class GRAPH,class EDGE_WEIGHTS,class SEEDS,class LABELS>
    void carvingSegmentation(
        const GRAPH                         & g,
        const EDGE_WEIGHTS                  & edgeWeights,
        const SEEDS                         & seeds,
        const typename LABELS::Value        backgroundLabel,
        const typename EDGE_WEIGHTS::Value  backgroundBias,
        LABELS                      & labels
    ){
        typedef typename EDGE_WEIGHTS::Value WeightType;
        typedef typename LABELS::Value       LabelType;
        detail_watersheds_segmentation::CarvingFunctor<WeightType,LabelType> f(backgroundLabel,backgroundBias);
        detail_watersheds_segmentation::watershedsSegmentationImpl(g,edgeWeights,seeds,f,labels);
    }



    template< class GRAPH , class EDGE_WEIGHTS, class NODE_SIZE,class NODE_LABEL_MAP>
    void felzenszwalbSegmentation(
        const GRAPH &         graph,
        const EDGE_WEIGHTS &  edgeWeights,
        const NODE_SIZE    &  nodeSizes,
        float           k,
        NODE_LABEL_MAP     &  nodeLabeling,
        const int             nodeNumStopCond = -1
    ){
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );
        typedef typename EDGE_WEIGHTS::Value WeightType;
        typedef typename EDGE_WEIGHTS::Value NodeSizeType;
        typedef typename Graph:: template NodeMap<WeightType>   NodeIntDiffMap;
        typedef typename Graph:: template NodeMap<NodeSizeType> NodeSizeAccMap;

        typedef NodeMapIteratorHelper<GRAPH,NODE_SIZE  >      NodeSizeMapHelper;
        typedef NodeMapIteratorHelper<GRAPH,NodeSizeAccMap  > NodeAccSizeMapHelper;
        typedef NodeMapIteratorHelper<GRAPH,NodeIntDiffMap  > NodeIntDiffMapHelper;
        typedef detail::Partition<size_t> UfdType;

        // initalize node size map  and internal diff map
        NodeIntDiffMap internalDiff(graph);
        NodeSizeAccMap nodeSizeAcc(graph);  
        copyNodeMap(graph,nodeSizes,nodeSizeAcc);
        fillNodeMap(graph,internalDiff,static_cast<WeightType>(0.0));



        // initlaize internal node diff map

        // sort the edges by their weights
        std::vector<Edge> sortedEdges;
        std::less<WeightType> comperator;
        edgeSort(graph,edgeWeights,comperator,sortedEdges);

        // make the ufd
        UfdType ufd(graph.maxNodeId()+1);

        size_t nodeNum = graph.nodeNum();


        while(true){
            // iterate over edges is the sorted order
            for(size_t i=0;i<sortedEdges.size();++i){
                const Edge e  = sortedEdges[i];
                const size_t rui = ufd.find(graph.id(graph.u(e)));
                const size_t rvi = ufd.find(graph.id(graph.v(e)));
                const Node   ru  = graph.nodeFromId(rui);
                const Node   rv  = graph.nodeFromId(rvi);
                if(rui!=rvi){

                    //check if to merge or not ?
                    const WeightType   w         = edgeWeights[e];
                    const NodeSizeType sizeRu    = nodeSizeAcc[ru];
                    const NodeSizeType sizeRv    = nodeSizeAcc[rv];
                    const WeightType tauRu       = static_cast<WeightType>(k)/static_cast<WeightType>(sizeRu);
                    const WeightType tauRv       = static_cast<WeightType>(k)/static_cast<WeightType>(sizeRv);
                    const WeightType minIntDiff  = std::min(internalDiff[ru]+tauRu,internalDiff[rv]+tauRv);
                    if(w<=minIntDiff){
                        // do merge
                        ufd.merge(rui,rvi);
                        --nodeNum;
                        // update size and internal difference
                        const size_t newRepId = ufd.find(rui);
                        const Node newRepNode = graph.nodeFromId(newRepId);
                        internalDiff[newRepNode]=w;
                        nodeSizeAcc[newRepNode] = sizeRu+sizeRv;
                    }
                }
                if(nodeNum==nodeNumStopCond){
                    break;
                }
            }
            if(nodeNumStopCond==-1){
                break;
            }
            else{
                if(nodeNum>nodeNumStopCond){
                    k*=1.2;
                }
                else{
                    break;
                }
            }
        }
        for(typename  GRAPH::NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            nodeLabeling[node]=ufd.find(graph.id(node));
        }
    } 




    namespace detail_graph_smoothing{

    template<
        class GRAPH, 
        class NODE_FEATURES_IN,
        class EDGE_WEIGHTS,
        class WEIGHTS_TO_SMOOTH_FACTOR,
        class NODE_FEATURES_OUT
    >
    void graphSmoothingImpl(
        const GRAPH & g,
        const NODE_FEATURES_IN   & nodeFeaturesIn,
        const EDGE_WEIGHTS       & edgeWeights,
        WEIGHTS_TO_SMOOTH_FACTOR & weightsToSmoothFactor,
        NODE_FEATURES_OUT        & nodeFeaturesOut

    ){
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRAPH_TYPEDEFS(Graph, , );
        typedef typename NODE_FEATURES_IN::ConstReference NodeFeatureInConstRef;
        typedef typename NODE_FEATURES_IN::Value          NodeFeatureInValue;
        typedef typename NODE_FEATURES_OUT::Reference     NodeFeatureOutRef;
        typedef typename EDGE_WEIGHTS::ConstReference SmoothFactorType;


        fillNodeMap(g,nodeFeaturesOut,0.0);

        for(NodeIt n(g);n!=lemon::INVALID;++n){

            const Node node(*n);

            NodeFeatureInValue    featIn  = nodeFeaturesIn[node];
            NodeFeatureOutRef     featOut = nodeFeaturesOut[node];

            featOut=0;
            float weightSum = 0.0;
            size_t degree    = 0;
            for(OutArcIt a(g,node);a!=lemon::INVALID;++a){
                const Edge edge(*a);
                const Node neigbour(g.target(*a));
                SmoothFactorType smoothFactor= weightsToSmoothFactor(edgeWeights[edge]);

                NodeFeatureInValue neighbourFeat = nodeFeaturesIn[neigbour];
                neighbourFeat*=smoothFactor;
                featOut += neighbourFeat;
                weightSum+=smoothFactor;
                ++degree;
            }
            // fixme..set me to right type 
            featIn*=static_cast<float>(degree);
            weightSum+=static_cast<float>(degree);
            featOut+=featIn;
            featOut/=weightSum;
        }
    }

    template<class T>
    struct ExpSmoothFactor{
        ExpSmoothFactor(const T lambda,const T edgeThreshold,const T scale)
        :   lambda_(lambda),
            edgeThreshold_(edgeThreshold),
            scale_(scale){
        }
        T operator()(const T weight){
            return weight> edgeThreshold_ ? 0 :  std::exp(-1.0*lambda_*weight)*scale_;
        }
        T lambda_;
        T edgeThreshold_;
        T scale_;
    };



    }

    template<class GRAPH, class NODE_FEATURES_IN,class EDGE_INDICATOR,class NODE_FEATURES_OUT>
    void graphSmoothing(
        const GRAPH & g,
        const NODE_FEATURES_IN  & nodeFeaturesIn,
        const EDGE_INDICATOR    & edgeIndicator,
        const float lambda,
        const float edgeThreshold,
        const float scale,
        NODE_FEATURES_OUT       & nodeFeaturesOut
    ){
        detail_graph_smoothing::ExpSmoothFactor<float> functor(lambda,edgeThreshold,scale);
        detail_graph_smoothing::graphSmoothingImpl(g,nodeFeaturesIn,edgeIndicator,functor,nodeFeaturesOut);
    }


    template<class GRAPH, class NODE_FEATURES_IN,class EDGE_INDICATOR,class NODE_FEATURES_OUT>
    void recursiveGraphSmoothing(
        const GRAPH & g,
        const NODE_FEATURES_IN   & nodeFeaturesIn,
        const EDGE_INDICATOR     & edgeIndicator,
        const float lambda,
        const float edgeThreshold,
        const float scale,
        size_t                    iterations,
        NODE_FEATURES_OUT       & nodeFeaturesBuffer,
        NODE_FEATURES_OUT       & nodeFeaturesOut
    ){

        iterations = std::max(size_t(1),iterations);
        // initial run
        graphSmoothing(g,nodeFeaturesIn,edgeIndicator,lambda,edgeThreshold,scale,nodeFeaturesOut);
        iterations -=1;

        bool outAsIn=true;
        for(size_t i=0;i<iterations;++i){
            if(outAsIn){
                graphSmoothing(g,nodeFeaturesOut,edgeIndicator,lambda,edgeThreshold,scale,nodeFeaturesBuffer);
                outAsIn=false;
            }
            else{
                graphSmoothing(g,nodeFeaturesBuffer,edgeIndicator,lambda,edgeThreshold,scale,nodeFeaturesOut);
                outAsIn=true;
            }
        }
        if(!outAsIn){
            copyNodeMap(g,nodeFeaturesBuffer,nodeFeaturesOut);
        }
    }

} // namespace vigra

#endif // VIGRA_GRAPH_MAP_ALGORITHMS_HXX
