#ifndef VIGRA_NEW_RAG_HXX
#define VIGRA_NEW_RAG_HXX

namespace vigra{





template<
    class GRAPH_IN,
    class GRAPH_IN_NODE_LABEL_MAP,
>
void makeRegionAdjacencyGraph(
    GRAPH_IN                   graphIn,
    GRAPH_IN_NODE_LABEL_MAP    labels
){

    typedef typename GRAPH_IN_NODE_LABEL_MAP::Value LabelType;
    typedef GRAPH_IN GraphIn;
    typedef typename GraphIn::Edge EdgeGraphIn;
    typedef typename GraphIn::Node NodeGraphIn;
    typedef typename GraphIn::EdgeIt EdgeItInGraph;
    typedef typename GraphIn::NodeIt NodeItInGraph;
    
    typedef AdjacencyListGraph<1> GraphOut;
    typedef typename GraphOut::Edge EdgeGraphOut;
    typedef typename GraphOut::Node NodeGraphOut;
    typedef typename GraphOut::EdgeIt EdgeItOutGraph;
    typedef typename GraphOut::NodeIt NodeItOutGraph;
    typedef typename GraphOut:: template EdgeMap< std::vector<EdgeGraphIn> > HyperEdgeMap;
    typedef typename GraphOut:: template NodeMap< std::vector<NodeGraphIn> > HyperNodeMap;

    
    // iterate over all labels in the node map to find min max
    std::pair<LabelType,LabelType> minMaxLabel  = nodeMapMinMap(graphIn,labels);

    size_t numberOfNodes = 0; // TODO
    size_t reserveEdges  = 0;

    // SET UP RAG GAPPH
    AdjacencyListGraph<1> ragGraph(numberOfNodes,reserveEdges);
    // add nodes
    for(size_t n=0;n<numberOfNodes;++n){
        ragGraph.addNdoe();
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


    HyperEdgeMap hyperEdges(ragGraph);
    HyperNodeMap hyperNodes(ragGraph);

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

#endif // VIGRA_NEW_RAG_HXX