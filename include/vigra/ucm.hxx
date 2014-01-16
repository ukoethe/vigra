#ifndef VIGRA_UCM_HXX
#define VIGRA_UCM_HXX

#include <valgrind/callgrind.h>

#include <vigra/merge_graph/new_merge_graph.hxx>
#include <vigra/merge_graph/maps/multi_array_map.hxx>
#include <vigra/merge_graph/maps/merge_graph_view_map.hxx>
namespace vigra{      


    template<class G,class EDGE_INDICATOR_MAP,class EDGE_SIZE_MAP>
    void ucmTransform(
        const G &             g,
        EDGE_INDICATOR_MAP &  edgeIndicator,
        EDGE_SIZE_MAP      &  edgeSizeMap
    ){

        CALLGRIND_START_INSTRUMENTATION;

        // typedefs graph and merge graph
        typedef typename G::Edge     Edge;
        typedef typename G::EdgeIt   EdgeIt;
        typedef MergeGraphAdaptor<G> MergeGraph;
        typedef typename MergeGraph::Edge MergeGraphEdge;

        // typedefs edge indactor map
        typedef typename EDGE_INDICATOR_MAP::Value ValueType;
        //typedef edge size map 
        typedef typename  EDGE_SIZE_MAP::Value SizeType;
        typedef view_maps::SumMap<0,SizeType>  MgEdgeSizeMap;

        typedef WeightedMeanMap<MergeGraph,EDGE_INDICATOR_MAP,EDGE_SIZE_MAP,Edge> MgWeightedMeanView;
        typedef MinWeightEdgeViewMap<MergeGraph,MgWeightedMeanView>               MgMinWeightOperator;

        MergeGraph mg(g);

        MgWeightedMeanView  mgWeightedMeanView(mg,edgeIndicator,edgeSizeMap);
        MgMinWeightOperator minWeightOperator(mg,mgWeightedMeanView);


       



        std::cout<<"Register callbacks\n";

        //mg.registerMergeEdgeCallBack(mgWeightedMeanView,& MgWeightedMeanView::merge);
        mg.registerMergeEdgeCallBack(minWeightOperator,& MgMinWeightOperator::mergeEdges);
        mg.registerEraseEdgeCallBack(minWeightOperator,& MgMinWeightOperator::eraseEdge);


        while(mg.nodeNum()>1){
            if(mg.nodeNum()%1000==0)
                std::cout<<mg.nodeNum()<<"\n";
                mg.mergeRegions( mg.id(minWeightOperator.minWeightEdge()));
        }

        for(EdgeIt e(g);e!=lemon::INVALID;++e){
            const Edge reprEdge =  mg.reprGraphEdge(*e);
            if(reprEdge!=*e){
                edgeIndicator[*e]=edgeIndicator[reprEdge];
                edgeSizeMap[*e]=edgeSizeMap[reprEdge];
            }
        }

        CALLGRIND_STOP_INSTRUMENTATION;
    }



}



#endif // VIGRA_UCM_HXX