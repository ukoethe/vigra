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


/*
using namespace vigra::acc;
vigra::MultiArray<3, double> data(...);
vigra::MultiArray<3, int> labels(...);
AccumulatorChainArray<CoupledArrays<3, double, int>,
Select<DataArg<1>, LabelArg<2>, // in which array to look (coordinates are always arg 0)
Mean, Variance, //per-region statistics over values
Coord<Mean>, Coord<Variance>, //per-region statistics over coordinates
Global<Mean>, Global<Variance> > > //global statistics
a;
a.ignoreLabel(0); //statistics will not be computed for region 0 (e.g. background)
extractFeatures(data, labels, a);
int regionlabel = ...;
std::cout << get<Mean>(a, regionlabel) << std::endl; //get Mean of region with label 'regionlabel'
*/


/** \brief Base class for any GraphMap (EdgeMap/NodeMap) which uses vigra::AccumulatorChain as ValueType

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
private:
    AccumulatorChainArrayType accChainArrayType_;
};





template<class MERGE_GRAPH,unsigned int DIM,class T, class SELECTED,bool DYNAMIC>
void AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::merge(
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelType a,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelType b
){
    accChainArrayType_.merge(a,b);
}



template<class MERGE_GRAPH,unsigned int DIM,class T, class SELECTED,bool DYNAMIC>
AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::AccumulatorChainMap(
    MERGE_GRAPH & graph,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::ValueTypeArrayType   & data,
    typename  AccumulatorChainMap<MERGE_GRAPH,DIM,T,SELECTED,DYNAMIC>::LabelArrayType       & labels
)
:   GraphMapBase<MERGE_GRAPH>(graph),
    accChainArrayType_()
{
    std::cout<<"extract Features\n";
    vigra::acc::extractFeatures(data,labels,accChainArrayType_);
}


}

#endif /*VIGRA_ACCUMULATOR_CHAIN_MAP*/