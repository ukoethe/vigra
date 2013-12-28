#ifndef VIGRA_ACCUMULATOR_CHAIN_MAP
#define VIGRA_ACCUMULATOR_CHAIN_MAP

/*std*/
#include <vector>



/*vigra*/
#include "graph_map_base.hxx"
#include <vigra/multi_array.hxx>
#include <vigra/accumulator.hxx>



namespace vigra{


namespace acc_map_detail{



}



/** \brief an Graph map which stores an vigra::AccumulatorChainArray 

<b>\#include</b> \<vigra/merge_graph/maps/accumulator_chain_map.hxx\> <br/>
Namespace: vigra
*/
template<
    class MERGE_GRAPH,
    unsigned int DIM,
    class T,
    class SELECTED,
    bool DYNAMIC=false
>
class AccumulatorChainMap
:   public GraphMapBase<MERGE_GRAPH>

{
public:
    // public typedefs
    typedef MERGE_GRAPH MergeGraphType;
    typedef typename MergeGraphType::LabelType LabelType;
    typedef T ValueType;
    typedef SELECTED SelectedStatisticsType;
    typedef vigra::CoupledArrays<DIM, ValueType ,LabelType> CoupledArraysType;
    typedef vigra::acc::AccumulatorChainArray<CoupledArraysType,SelectedStatisticsType> AccumulatorChainArrayType;



    typedef vigra::MultiArrayView<DIM,LabelType> LabelArrayType;
    typedef vigra::MultiArrayView<DIM,ValueType> ValueTypeArrayType;



    //typedef MergeGraphType::LabelType LabelType;
    //AccumulatorChainMap(){}

    
    AccumulatorChainMap(MergeGraphType & graph):GraphMapBase<MergeGraphType>(graph){}
    AccumulatorChainMap(MergeGraphType & , ValueTypeArrayType & ,LabelArrayType & );

    /**
    * merge two accumulator chains inplace
    * @param a will store the result, a=merge(a,b)
    * @param b will be merged in a
    * @return void/nothing
    */
    void merge(const LabelType a,const LabelType b);


    template<class TAG>
    struct TagResult
    {
        typedef typename acc::LookupTag<TAG, AccumulatorChainArrayType>::result_type result_type;
    };


    template<class TAG>
    typename acc::LookupTag<TAG, AccumulatorChainArrayType>::result_type
    get(const LabelType label)const{
        return acc::getAccumulator<TAG>(accChainArray_, label).get();
    }

    const AccumulatorChainArrayType & accChainArray()const{
        return accChainArray_;
    }

private:
    AccumulatorChainArrayType accChainArray_;
};




template<class ACC_CHAIN_MAP,class TAG>
class AccumulatorChainMapTagView{
    public:
        typedef ACC_CHAIN_MAP AccumulatorChainMapType;
        typedef typename AccumulatorChainMapType::LabelType LabelType;
        typedef typename AccumulatorChainMapType::AccumulatorChainArrayType AccumulatorChainArrayType;
        typedef typename acc::LookupTag<TAG, AccumulatorChainArrayType>::result_type result_type;

        AccumulatorChainMapTagView(const AccumulatorChainMapType & accChainMap)
        : accChainMap_(accChainMap){

        }

        result_type operator[](const LabelType label)const{
            return accChainMap_. template  get<TAG>(label);
        }
    private:
        const AccumulatorChainMapType & accChainMap_;
};




template<class MERGE_GRAPH,unsigned int DIM,class T, class SELECTED,bool DYNAMIC>
void AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::merge(
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelType a,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelType b
){
    accChainArray_.merge(a,b);
}



template<class MERGE_GRAPH,unsigned int DIM,class T, class SELECTED,bool DYNAMIC>
AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::AccumulatorChainMap(
    MERGE_GRAPH & graph,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::ValueTypeArrayType   & data,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelArrayType       & labels
)
:   GraphMapBase<MERGE_GRAPH>(graph),
    accChainArray_()
{
    accChainArray_.setMaxRegionLabel(graph.initNumberOfNodes()-1);
    vigra::acc::extractFeatures(data,labels,accChainArray_);
}


}

#endif /*VIGRA_ACCUMULATOR_CHAIN_MAP*/