#ifndef VIGRA_NEW_RAG_HXX
#define VIGRA_NEW_RAG_HXX

#include <vigra/graph_generalization.hxx>
#include <vigra/graph_map_algorithms.hxx>
#include <vigra/numpy_array.hxx>

namespace vigra{





template<
    class GRAPH_IN,
    class GRAPH_IN_NODE_LABEL_MAP
>
void makeRegionAdjacencyGraph(
    GRAPH_IN                   graphIn,
    GRAPH_IN_NODE_LABEL_MAP    labels,
    typename AdjacencyListGraph:: template EdgeMap< std::vector<typename GRAPH_IN::Edge> > & hyperEdges,
    typename AdjacencyListGraph:: template NodeMap< std::vector<typename GRAPH_IN::Node> > & hyperNodes,
    const Int64   ignoreLabel=-1
){
    typedef GRAPH_IN_NODE_LABEL_MAP LabelMap;
    typedef typename GraphMapTypeTraits<GRAPH_IN_NODE_LABEL_MAP>::Value LabelType;
    typedef GRAPH_IN GraphIn;
    typedef typename GraphIn::Edge EdgeGraphIn;
    typedef typename GraphIn::Node NodeGraphIn;
    typedef typename GraphIn::EdgeIt EdgeItInGraph;
    typedef typename GraphIn::NodeIt NodeItInGraph;
    
    typedef AdjacencyListGraph      GraphOut;
    typedef typename GraphOut::Edge EdgeGraphOut;
    typedef typename GraphOut::Node NodeGraphOut;
    typedef typename GraphOut::EdgeIt EdgeItOutGraph;
    typedef typename GraphOut::NodeIt NodeItOutGraph;
    typedef typename GraphOut:: template EdgeMap< std::vector<EdgeGraphIn> > HyperEdgeMap;
    typedef typename GraphOut:: template NodeMap< std::vector<NodeGraphIn> > HyperNodeMap;

    



    // iterate over all labels in the node map to find min max
    typedef NodeMapIteratorHelper<GraphIn,LabelMap> NodeIterHelper;
    const LabelType  minLabel = *std::min_element(NodeIterHelper::begin(graphIn,labels),NodeIterHelper::end(graphIn,labels));
    //const LabelType  maxLabel = *std::max_element(NodeIterHelper::begin(graphIn,labels),NodeIterHelper::end(graphIn,labels));


    size_t numberOfNodes = 0; // TODO
    size_t reserveEdges  = 0;

    // SET UP RAG GAPPH
    AdjacencyListGraph ragGraph(numberOfNodes,reserveEdges);
    // add nodes
    for(size_t n=0;n<numberOfNodes;++n){
        ragGraph.addNode();
    }
    // add al edges
    for(EdgeItInGraph e(graphIn);e!=lemon::INVALID;++e){
        const EdgeGraphIn edge(*e);
        const LabelType lu = labels[graphIn.u(edge)];
        const LabelType lv = labels[graphIn.v(edge)];
        if(lu!=lv){
            // if there is an edge between lu and lv no new edge will be added
            ragGraph.addEdge(ragGraph.nodeFromId(lu),ragGraph.nodeFromId(lv));
        }
    }

    // SET UP HYPEREDGES


    hyperEdges.assign(ragGraph);
    hyperNodes.assign(ragGraph);

    // add edges
    for(EdgeItInGraph e(graphIn);e!=lemon::INVALID;++e){
        const EdgeGraphIn edge(*e);
        const LabelType lu = labels[graphIn.u(edge)];
        const LabelType lv = labels[graphIn.v(edge)];
        if(lu!=lv){
            EdgeGraphOut ragEdge= ragGraph.findEdge(ragGraph.nodeFromId(lu),ragGraph.nodeFromId(lv));
            hyperEdges[ragEdge].push_back(edge);
        }
    }

    // add  nodes
    for(NodeItInGraph n(graphIn);n!=lemon::INVALID;++n){
        const NodeGraphIn node(*n);
        const LabelType ln = labels[node];
        const NodeGraphOut ragNode =ragGraph.nodeFromId(ln);
        hyperNodes[ragNode].push_back(node);
    }

}

} // end namespace vigra

#endif // VIGRA_NEW_RAG_HXX