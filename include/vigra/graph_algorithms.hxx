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
#include <vigra/numpy_array.hxx>
#include <vigra/priority_queue.hxx>
#include <vigra/union_find.hxx>
#include <vigra/adjacency_list_graph.hxx>


#define LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(GRAPH_CLS,PREFIX,POSTFIX) \
    typedef typename GRAPH_CLS::Edge        PREFIX##Edge      ## POSTFIX; \
    typedef typename GRAPH_CLS::Node        PREFIX##Node      ## POSTFIX; \
    typedef typename GRAPH_CLS::Arc         PREFIX##Arc       ## POSTFIX; \
    typedef typename GRAPH_CLS::EdgeIt      PREFIX##EdgeIt##POSTFIX; \
    typedef typename GRAPH_CLS::NodeIt      PREFIX##NodeIt    ## POSTFIX; \
    typedef typename GRAPH_CLS::ArcIt       PREFIX##ArcIt     ## POSTFIX; \
    typedef typename GRAPH_CLS::OutArcIt    PREFIX##OutArcIt  ## POSTFIX; \
    typedef typename GRAPH_CLS::OutArcIt    PREFIX##OutArcIt  ## POSTFIX


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



    template<class GRAPH,class HYPER_EDGE_MAP,class HYPER_EDGE_SIZE_MAP>
    void hyperEdgeSizes(
        const GRAPH & g,
        const HYPER_EDGE_MAP & hyperEdges,
        HYPER_EDGE_SIZE_MAP  & out
    ){
        for(typename  GRAPH::EdgeIt e(g);e!=lemon::INVALID;++e){
            const size_t size = hyperEdges[*e].size();
            out[*e]=static_cast<typename HYPER_EDGE_SIZE_MAP::Value>(size);
        }
    }

    template<class RAG,class GRAPH,class NODE_LABEL_MAP,class HYPER_NODE_SIZE_MAP>
    void hyperNodeSizes(
        const RAG   & rag,
        const GRAPH & g,
        const NODE_LABEL_MAP & nodeLabelMap,
        HYPER_NODE_SIZE_MAP  & out
    ){
        typedef NodeMapIteratorHelper<GRAPH,NODE_LABEL_MAP> LabelIterHelper;
        typedef NodeMapIteratorHelper<RAG,HYPER_NODE_SIZE_MAP> OutIterHelper;
        typedef typename NODE_LABEL_MAP::Value LabelType;
        typedef typename HYPER_NODE_SIZE_MAP::Value SizeType;
        //std::fill(OutIterHelper::begin(rag,out),OutIterHelper::end(rag,out),static_cast<SizeType>(0));

        typename LabelIterHelper::const_iterator labelIter    = LabelIterHelper::begin(g,nodeLabelMap);
        typename LabelIterHelper::const_iterator labelIterEnd = LabelIterHelper::end(g,nodeLabelMap);

        std::vector<SizeType> labelCount(rag.maxNodeId()+1);
        while(labelIter!=labelIterEnd){
            const LabelType l = *labelIter;
            if(l<=rag.maxNodeId()){
                labelCount[l]+=static_cast<SizeType>(1);
            }
            ++labelIter;
        }
        for(typename RAG::NodeIt n(rag);n!=lemon::INVALID;++n){
            out[*n]=labelCount[rag.id(*n)];
        }
    }

    template<class RAG,class GRAPH,class NODE_LABEL_MAP,class IMAGE_MAP,class HYPER_NODE_FEATURE_MAP>
    void hyperNodeImageFeatures(
        const RAG   & rag,
        const GRAPH & g,
        const NODE_LABEL_MAP & nodeLabelMap,
        // if image map is const type there is error...
        IMAGE_MAP & imageMap,
        HYPER_NODE_FEATURE_MAP & out
    ){
        typedef NodeMapIteratorHelper<GRAPH,NODE_LABEL_MAP>       LabelIterHelper;
        typedef NodeMapIteratorHelper<RAG,HYPER_NODE_FEATURE_MAP> OutIterHelper;
        typedef typename NODE_LABEL_MAP::Value LabelType;
        std::fill(OutIterHelper::begin(rag,out),OutIterHelper::end(rag,out),0.0);
        std::vector<size_t> labelCount(rag.maxNodeId()+1);

        for(typename  GRAPH::NodeIt n(g);n!=lemon::INVALID;++n){
            const LabelType label = nodeLabelMap[*n];
            const typename RAG::Node ragNode = rag.nodeFromId(label);
            if(ragNode!=lemon::INVALID){
                labelCount[label]+=1;
                out[ragNode]+=imageMap[*n];
            }
        }
        for(typename  RAG::NodeIt n(rag);n!=lemon::INVALID;++n){
            out[*n]/=labelCount[rag.id(*n)];
        }    
    }

    template<class RAG,class GRAPH,class HYPER_EDGE_COORDINATE_MAP,class IMAGE,class HYPER_EDGE_FEATURE_MAP>
    void hyperEdgeImageFeatures(
        const RAG & rag,
        const GRAPH & graph,
        const HYPER_EDGE_COORDINATE_MAP & hyperEdges,
        const IMAGE & image,
        HYPER_EDGE_FEATURE_MAP  & out
    ){
        typedef typename HYPER_EDGE_COORDINATE_MAP::ConstReference HyperEdgeCoordVec;
        typedef typename HYPER_EDGE_COORDINATE_MAP::Value::value_type TinyVectorType;
        typedef HYPER_EDGE_FEATURE_MAP OutMap;
        typedef typename OutMap::Value Value;
        typedef typename IMAGE::difference_type ImageCoord;

        for(typename  RAG::EdgeIt e(rag);e!=lemon::INVALID;++e){
            HyperEdgeCoordVec hyperEdgeCoords = hyperEdges[*e];
            out[*e]=Value(0);
            const size_t nEdges = hyperEdgeCoords.size();
            for(size_t i=0;i<nEdges;++i){
                out[*e]+=image[graph.u( hyperEdgeCoords[i]) ];
                out[*e]+=image[graph.v( hyperEdgeCoords[i]) ];
            }
            // todo: replace me with functors normalization / outpu
            out[*e]/=static_cast<Value>(2*hyperEdgeCoords.size());
        }
    }



    template<
        class GRAPH_IN,
        class GRAPH_IN_NODE_LABEL_MAP
    >
    void makeRegionAdjacencyGraph(
        GRAPH_IN                   graphIn,
        GRAPH_IN_NODE_LABEL_MAP    labels,
        AdjacencyListGraph & rag,
        typename AdjacencyListGraph:: template EdgeMap< std::vector<typename GRAPH_IN::Edge> > & hyperEdges,
        const Int64   ignoreLabel=-1
    ){
        typedef GRAPH_IN_NODE_LABEL_MAP LabelMap;
        typedef typename GraphMapTypeTraits<GRAPH_IN_NODE_LABEL_MAP>::Value LabelType;
        typedef GRAPH_IN GraphIn;
        typedef AdjacencyListGraph GraphOut;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(GraphIn, , GraphIn);
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(GraphOut, ,GraphOut);
        typedef typename GraphOut:: template EdgeMap< std::vector<EdgeGraphIn> > HyperEdgeMap;
        typedef typename GraphOut:: template NodeMap< std::vector<NodeGraphIn> > HyperNodeMap;
        // iterate over all labels in the node map to find min max
        typedef NodeMapIteratorHelper<GraphIn,LabelMap> NodeIterHelper;
        const LabelType  minLabel = *std::min_element(NodeIterHelper::begin(graphIn,labels),NodeIterHelper::end(graphIn,labels));
        const LabelType  maxLabel = *std::max_element(NodeIterHelper::begin(graphIn,labels),NodeIterHelper::end(graphIn,labels));

        if(ignoreLabel==0 && rag.zeroStart())
            throw std::runtime_error("if ignoreLabel == 0, rag.zeroStart() must be false");
        if( minLabel==1 && rag.zeroStart() )
            throw std::runtime_error("if minLabel==1, rag.zeroStart() must be false");
        if( ignoreLabel==1  && !rag.zeroStart())
            throw std::runtime_error("if ignoreLabel == 1, rag.zerStart() must be true");
        if( ignoreLabel==1  && minLabel!=0)
            throw std::runtime_error("if ignoreLabel == 1, minLabel must be 0");


        if(ignoreLabel!=0 || minLabel!= 1 || rag.zeroStart()){
            throw std::runtime_error("WORK IN PROGRESS LIMITATION VIOLATED  => this should evaluate to false (ignoreLabel!=0 || minLabel!= 1 || rag.zeroStart() ) ");
        }



        size_t numberOfNodes = maxLabel; // TODO
        size_t reserveEdges  = 0;

        // SET UP RAG GAPPH
        // add nodes
        for(size_t n=0;n<numberOfNodes;++n){
            rag.addNode();
        }
        // add al edges
        for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
            const EdgeGraphIn edge(*e);
            const LabelType lu = labels[graphIn.u(edge)];
            const LabelType lv = labels[graphIn.v(edge)];
            if(lu!=lv){
                // if there is an edge between lu and lv no new edge will be added
                rag.addEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
            }
        }

        // SET UP HYPEREDGES


        hyperEdges.assign(rag);
        // add edges
        for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
            const EdgeGraphIn edge(*e);
            const LabelType lu = labels[graphIn.u(edge)];
            const LabelType lv = labels[graphIn.v(edge)];
            if(lu!=lv){
                EdgeGraphOut ragEdge= rag.findEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
                hyperEdges[ragEdge].push_back(edge);
            }
        }
    }
        
    template<class GRAPH,class WEIGHTS,class PREDECESSORS,class DISTANCE>
    void shortestPathDijkstra(
        const GRAPH         &           graph,
        const typename GRAPH::Node &    source,
        const WEIGHTS       &           weights,
        PREDECESSORS        &           predecessors,
        DISTANCE            &           distance,
        const typename GRAPH::Node &    target  = lemon::INVALID
    ){

        typedef GRAPH                       Graph;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
        typedef typename WEIGHTS::value_type     WeightType;
        typedef typename DISTANCE::value_type    DistanceType;
        const size_t maxNodeId = graph.maxNodeId();
        vigra::ChangeablePriorityQueue<typename WEIGHTS::value_type> pq(maxNodeId+1);

        for(NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            pq.push(graph.id(node),std::numeric_limits<WeightType>::infinity() );
            distance[node]=std::numeric_limits<DistanceType>::infinity();
            predecessors[node]=lemon::INVALID;
        }

        distance[source]=static_cast<DistanceType>(0.0);
        pq.push(graph.id(source),0.0);

        bool finished=false;
        while(!pq.empty() && !finished){
            const WeightType minDist = pq.topPriority();
            if(minDist < std::numeric_limits<DistanceType>::infinity()){
                const Node topNode(graph.nodeFromId(pq.top()));
                pq.pop();
                // loop over all neigbours
                for(OutArcIt outArcIt(graph,topNode);outArcIt!=lemon::INVALID;++outArcIt){
                    const Node otherNode = graph.target(*outArcIt);
                    const size_t otherNodeId = graph.id(otherNode);

                    if(pq.contains(otherNodeId)){
                        const Edge edge(*outArcIt);
                        const DistanceType currentDist      = distance[otherNode];
                        const DistanceType alternativeDist  = distance[topNode]+weights[edge];
                        if(alternativeDist<currentDist){
                            pq.push(otherNodeId,alternativeDist);
                            distance[otherNode]=alternativeDist;
                            predecessors[otherNode]=topNode;
                        }
                    }
                    if(target==otherNode){
                        finished=true;
                        break;
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




    template<class GRAPH,class WEIGHTS,class PREDECESSORS,class DISTANCE>
    void shortestPathDijkstraFast(
        const GRAPH         &                                        graph,
        const typename GRAPH::Node &                                 source,
        const WEIGHTS       &                                        weights,
        PREDECESSORS        &                                        predecessors,
        DISTANCE            &                                        distance,
        vigra::ChangeablePriorityQueue<typename WEIGHTS::value_type> pq,
        const typename GRAPH::Node &    target  = lemon::INVALID,
        const bool fill=true
    ){

        typedef GRAPH                       Graph;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
        typedef typename WEIGHTS::value_type     WeightType;
        typedef typename DISTANCE::value_type    DistanceType;
        const size_t maxNodeId = graph.maxNodeId();
        //vigra::ChangeablePriorityQueue<typename WEIGHTS::value_type> pq(maxNodeId+1);

        for(NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            pq.push(graph.id(node),std::numeric_limits<WeightType>::infinity() );
            distance[node]=std::numeric_limits<DistanceType>::infinity();
            predecessors[node]=lemon::INVALID;
        }

        distance[source]=static_cast<DistanceType>(0.0);
        pq.push(graph.id(source),0.0);

        bool finished=false;
        while(!pq.empty() && !finished){
            const WeightType minDist = pq.topPriority();
            if(minDist < std::numeric_limits<DistanceType>::infinity()){
                const Node topNode(graph.nodeFromId(pq.top()));
                pq.pop();
                // loop over all neigbours
                for(OutArcIt outArcIt(graph,topNode);outArcIt!=lemon::INVALID;++outArcIt){
                    const Node otherNode = graph.target(*outArcIt);
                    const size_t otherNodeId = graph.id(otherNode);

                    if(pq.contains(otherNodeId)){
                        const Edge edge(*outArcIt);
                        const DistanceType currentDist      = distance[otherNode];
                        const DistanceType alternativeDist  = distance[topNode]+weights[edge];
                        if(alternativeDist<currentDist){
                            pq.push(otherNodeId,alternativeDist);
                            distance[otherNode]=alternativeDist;
                            predecessors[otherNode]=topNode;
                        }
                    }
                    if(target==otherNode){
                        finished=true;
                        break;
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
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );

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
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
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

    template<class GRAPH,class EDGE_WEIGHTS,class SEEDS,class LABELS>
    void watershedsSegmentation(
        const GRAPH & g,
        const EDGE_WEIGHTS & edgeWeights,
        const SEEDS        & seeds,
        LABELS             & labels
    ){  
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
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
                        pq.push(neigbour,edgeWeights[edge]);
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
                        if(labels[neigbour]!=0 && edgeWeights[edge]<minWeight){
                            minWeight=edgeWeights[edge];;
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
    

    template<class GRAPH,class EDGE_WEIGHTS,class SEEDS,class LABELS>
    void carvingSegmentation(
        const GRAPH & g,
        const EDGE_WEIGHTS & edgeWeights,
        const SEEDS        & seeds,
        LABELS             & labels
    ){
        
    }



    template< class GRAPH , class EDGE_WEIGHTS, class NODE_SIZE,class NODE_LABEL_MAP>
    void felzenszwalbSegmentation(
        const GRAPH &         graph,
        const EDGE_WEIGHTS &  edgeWeights,
        const NODE_SIZE    &  nodeSizes,
        const float           k,
        NODE_LABEL_MAP     &  nodeLabeling
    ){
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
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
                    // update size and internal difference
                    const size_t newRepId = ufd.find(rui);
                    const Node newRepNode = graph.nodeFromId(newRepId);
                    internalDiff[newRepNode]=w;
                    nodeSizeAcc[newRepNode] = sizeRu+sizeRv;
                }
            }
        }
        for(typename  GRAPH::NodeIt n(graph);n!=lemon::INVALID;++n){
            const Node node(*n);
            nodeLabeling[node]=ufd.find(graph.id(node));
        }
    } 






    template<class GRAPH, class NODE_FEATURES,class NODE_SIZES,class EDGE_WEIGHTS,class NODE_FEATURES_OUT>
    void graphSmoothing(
        const GRAPH & g,
        const NODE_FEATURES & nodeFeatures,
        const NODE_FEATURES & nodeSizes,
        const EDGE_WEIGHTS  & edgeWeights,
        NODE_FEATURES_OUT   & nodeFeaturesOut

    ){
        typedef GRAPH Graph;
        LEMON_UNDIRECTED_GRPAPH_TYPEDEFS(Graph, , );
        typedef typename NODE_FEATURES::Value NodeFeaturesValueType;
        for(NodeIt n(g);n!=lemon::INVALID;++n){
            const Node node(*n);
            NodeFeaturesValueType nodeFeature = nodeFeatures[n];

        }
    }

} // namespace vigra

#endif // VIGRA_GRAPH_MAP_ALGORITHMS_HXX
