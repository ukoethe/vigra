#ifndef VIGRA_GRAPH_MAP_BASE
#define VIGRA_GRAPH_MAP_BASE


namespace vigra{

template<class DYNAMIC_GRAPH>
class GraphMapBase{
public:
    typedef DYNAMIC_GRAPH DynamicGraphType;
    GraphMapBase(DynamicGraphType &);


    const DynamicGraphType & graph() const ;
protected:
    DynamicGraphType & graph();

    DynamicGraphType & graph_;


};  





// GRAPH BASE

template<class DYNAMIC_GRAPH>
inline GraphMapBase<DYNAMIC_GRAPH>::GraphMapBase(
    typename GraphMapBase<DYNAMIC_GRAPH>::DynamicGraphType & graph
)
:   graph_(graph){

}   


template<class DYNAMIC_GRAPH>
inline const typename GraphMapBase<DYNAMIC_GRAPH>::DynamicGraphType & 
GraphMapBase<DYNAMIC_GRAPH>::graph()const{
    return graph_;
}

template<class DYNAMIC_GRAPH>
inline typename GraphMapBase<DYNAMIC_GRAPH>::DynamicGraphType & 
GraphMapBase<DYNAMIC_GRAPH>::graph(){
    return graph_;
}



}


#endif /*VIGRA_GRAPH_MAP_BASE*/