    
#ifndef VIGRA_MERGE_GRAPH_ITEM_HXX
#define VIGRA_MERGE_GRAPH_ITEM_HXX


/* this project*/
#include "merge_graph.hxx"


namespace vigra {


namespace merge_graph_detail{

    template<class MERGE_GRAPH,class NODE_OR_EDGE>
    struct GraphItem;


    template<class MERGE_GRAPH>
    struct GraphItem<MERGE_GRAPH,typename MERGE_GRAPH::Edge> {
        typedef MERGE_GRAPH MergeGraphType;
        typedef  typename MERGE_GRAPH::Edge result;
        typedef  typename MergeGraphType::EdgeIdIt  GraphItemIdIt;

    };


    template<class MERGE_GRAPH>
    struct GraphItem<MERGE_GRAPH,typename MERGE_GRAPH::Node> {
        typedef MERGE_GRAPH MergeGraphType;
        typedef  const typename MERGE_GRAPH::Node & result;
        typedef  typename MergeGraphType::NodeIdIt  GraphItemIdIt;

    };


}





} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_ITEM_HXX