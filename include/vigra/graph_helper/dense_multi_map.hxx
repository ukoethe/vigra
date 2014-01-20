#ifndef VIGRA_DENSE_REFERENCE_MULTI_MAP_HXX
#define VIGRA_DENSE_REFERENCE_MULTI_MAP_HXX

/*vigra*/
//#include <vigra/array_vector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graphs.hxx>

namespace vigra{


namespace detail_dense_multi_map{

template<unsigned int DIM,class T>
typename MultiArray<DIM+1,T>::difference_type shapeHelper(
    const size_t  maxItemId,
    const typename MultiArray<DIM,T>::difference_type & valueShape
){
    typename MultiArray<DIM+1,T>::difference_type resultShape;
    resultShape[0]=maxItemId+1;
    for(size_t i=0;i<valueShape.size();++i){
        resultShape[i+1]=valueShape[i];
    }
    return resultShape;
}

} // end detail_dense_multi_map



template<unsigned int DIM,class T,class KEY>
class DenseReferenceMultiMap
: public MultiArray<DIM+1,T>
{
public:
    typedef KEY                             Key;
    typedef MultiArray<DIM,T>               Value;
    typedef MultiArrayView<DIM,T>           Reference;
    typedef MultiArrayView<DIM,const T>     ConstReference;

    typedef typename MultiArray<DIM  ,T>::difference_type ValueShapeType;
    typedef typename MultiArray<DIM+1,T>::difference_type TotalShapeType;

    DenseReferenceMultiMap(const Key & maxKey, const ValueShapeType valueShape)
    : MultiArray<DIM+1,T>(detail_dense_multi_map::shapeHelper(maxKey,valueShape)){

    }

    ConstReference operator[](const KEY & key)const{
        //return this->operator[](key.id());
        return this->bindInner(key);
    }
    Reference operator[](const KEY & key){
        return this->bindInner(key);
    }

    size_t size()const{
        return this->shape(0);
    }
private:
    // NONE
};


template<class GRAPH,class ITEM,unsigned int DIM,class T>
class DenseGraphItemReferenceMultiMap
: public DenseReferenceMultiMap<DIM,T,ITEM>
{
    typedef GRAPH Graph;
    typedef ITEM  Item;
    typedef DenseReferenceMultiMap<DIM,T,ITEM> DenseReferenceMultiMapType;
    typedef GraphItemHelper<Graph,ITEM> ItemHelper;
    typedef typename ItemHelper::ItemIt ItemIt;

public:
    DenseGraphItemReferenceMultiMap(const Graph & g,typename DenseReferenceMultiMapType::ValueShapeType & valueShape)
    :   DenseReferenceMultiMapType(ItemHelper::maxItemId(g),valueShape){

    }
};


template<class GRAPH,unsigned int DIM,class T>
class DenseNodeReferenceMultiMap
: public DenseGraphItemReferenceMultiMap<GRAPH,typename GRAPH::Node,DIM,T>
{
    typedef typename GRAPH::Node Node;
    typedef DenseGraphItemReferenceMultiMap<GRAPH,Node,DIM,T> DenseGraphItemReferenceMultiMapType;
    public:
    DenseNodeReferenceMultiMap(const GRAPH & g,typename DenseGraphItemReferenceMultiMapType::ValueShapeType & valueShape)
    : DenseGraphItemReferenceMultiMapType(g,valueShape){
    }
};


template<class GRAPH,unsigned int DIM,class T>
class DenseEdgeReferenceMultiMap
: public DenseGraphItemReferenceMultiMap<GRAPH,typename GRAPH::Edge,DIM,T>
{
    typedef typename GRAPH::Edge Edge;
    typedef DenseGraphItemReferenceMultiMap<GRAPH,Edge,DIM,T> DenseGraphItemReferenceMultiMapType;
    public:
    DenseEdgeReferenceMultiMap(const GRAPH & g,typename DenseGraphItemReferenceMultiMapType::ValueShapeType & valueShape)
    : DenseGraphItemReferenceMultiMapType(g,valueShape){
    }
};


template<class GRAPH,unsigned int DIM,class T>
class DenseArcReferenceMultiMap
: public DenseGraphItemReferenceMultiMap<GRAPH,typename GRAPH::Arc,DIM,T>
{
    typedef typename GRAPH::Arc Arc;
    typedef DenseGraphItemReferenceMultiMap<GRAPH,Arc,DIM,T> DenseGraphItemReferenceMultiMapType;
    public:
    DenseArcReferenceMultiMap(const GRAPH & g,typename DenseGraphItemReferenceMultiMapType::ValueShapeType & valueShape)
    : DenseGraphItemReferenceMultiMapType(g,valueShape){
    }
};


}

#endif //VIGRA_DENSE_REFERENCE_MULTI_MAP_HXX 