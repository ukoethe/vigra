#ifndef VIGRA_UCM_HXX
#define VIGRA_UCM_HXX

#include <vigra/merge_graph/new_merge_graph.hxx>
#include <vigra/merge_graph/maps/multi_array_map.hxx>


namespace vigra{

	template<class G,class EDGE_INDICATOR_MAP,class EDGE_OUT_UCM_MAP>
    void ucmTransform(
        const G &                   g,
        const EDGE_INDICATOR_MAP &  edgeIndicator,
        EDGE_OUT_UCM_MAP   &  edgeOutUcmMap
    ){
        
        typedef MergeGraphAdaptor<G> MergeGraph;
        typedef typename EDGE_INDICATOR_MAP::Value ValueType;
        typedef view_maps::SumMap<0,ValueType> UcmWeightMap;
        typedef view_maps::WeightedMeanMap<0,ValueType,UcmWeightMap > UcmEdgeIndicatorMap;
        typedef view_maps::MinWeightEdgeMapSimple<ValueType,MergeGraph,UcmEdgeIndicatorMap> MinMap;

        MergeGraph mg(g);
        std::cout<<"START UCM TRANSFORM \n";
        std::cout<<"nodes "<<mg.nodeNum();
        std::cout<<"edges "<<mg.edgeNum();

        MultiArray<1,ValueType> sizeDummy(typename MultiArray<1,ValueType>::difference_type(g.maxEdgeId()+1)   );
       		

        UcmWeightMap 	     mgEdgeSizeMap(sizeDummy);
       	UcmEdgeIndicatorMap  mgEdgeIndicatorMap(mgEdgeSizeMap,mgEdgeSizeMap);
       	MinMap 				 minEdgeMap(mg,edgeOutUcmMap,mgEdgeIndicatorMap);

       	
       	std::cout<<"Register callbacks\n";

       	mg.registerMergeEdgeCallBack(mgEdgeSizeMap,     & UcmWeightMap::merge);
       	mg.registerMergeEdgeCallBack(mgEdgeIndicatorMap,& UcmEdgeIndicatorMap::merge);
       	mg.registerMergeEdgeCallBack(minEdgeMap,        & MinMap::mergeEdges);
       	mg.registerEraseEdgeCallBack(minEdgeMap,        & MinMap::eraseEdge);


       	while(mg.nodeNum()>1){
       		if(mg.nodeNum()%100==0)
       			std::cout<<mg.nodeNum()<<"\n";
       		mg.mergeRegions( minEdgeMap.minWeightEdgeLabel());
       	}
    }


}



#endif // VIGRA_UCM_HXX