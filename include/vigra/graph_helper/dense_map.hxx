#ifndef VIGRA_DENSE_REFERENCE_MAP_HXX
#define VIGRA_DENSE_REFERENCE_MAP_HXX

/*vigra*/
//#include <vigra/array_vector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graphs.hxx>

namespace vigra{

template<class T,class KEY,class REF,class CREF>
class DenseReferenceMap
: public MultiArray<1,T>
{
public:
	typedef KEY  Key;
	typedef T 	 Value;
	typedef REF  Reference;
	typedef CREF ConstReference;

	typedef typename MultiArray<1,T>::difference_type Shape1Type;

	DenseReferenceMap()
	: MultiArray<1,T>(){
	}
	DenseReferenceMap(const size_t maxKey)
	: MultiArray<1,T>(Shape1Type(maxKey+1)){
	}
	DenseReferenceMap(const size_t maxKey,ConstReference  value)
	: MultiArray<1,T>(Shape1Type(maxKey+1),value){
	}




	ConstReference operator[](const KEY & key)const{
		//return this->operator[](key.id());
		return MultiArray<1,T>::operator()(key.id());
	}
	Reference operator[](const KEY & key){
		return MultiArray<1,T>::operator()(key.id());
	}

	size_t size()const{
		return this->shape(0);
	}
protected:
	void assign(const size_t maxKey){
		this->reshape(Shape1Type(maxKey+1));
	}
private:
	// NONE
};


template<class GRAPH,class ITEM,class T,class REF,class CREF>
class DenseGraphItemReferenceMap
: public DenseReferenceMap<T,ITEM,REF,CREF>
{
	typedef GRAPH Graph;
	typedef ITEM  Item;
	typedef DenseReferenceMap<T,ITEM,REF,CREF> DenseReferenceMapType;
	typedef GraphItemHelper<Graph,ITEM> ItemHelper;
	typedef typename ItemHelper::ItemIt ItemIt;

public:
	DenseGraphItemReferenceMap()
	:	DenseReferenceMapType(){

	}
	DenseGraphItemReferenceMap(const Graph & g)
	:	DenseReferenceMapType(ItemHelper::maxItemId(g)){

	}
	DenseGraphItemReferenceMap(const Graph & g,typename DenseReferenceMapType::ConstReference value)
	:	DenseReferenceMapType(ItemHelper::maxItemId(g),value){

	}
	void assign(const Graph & g){
		DenseReferenceMapType::assign(ItemHelper::maxItemId(g));
	}
};


template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseNodeReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Node,T,REF,CREF>
{
	typedef typename GRAPH::Node Node;
	typedef DenseGraphItemReferenceMap<GRAPH,Node,T,REF,CREF> DenseGraphItemReferenceMapType;
	public:
		DenseNodeReferenceMap()
		: DenseGraphItemReferenceMapType(){
		}
		DenseNodeReferenceMap(const GRAPH & g)
		: DenseGraphItemReferenceMapType(g){
		}
		DenseNodeReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
		: DenseGraphItemReferenceMapType(g,value){
		}
};

template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseEdgeReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Edge,T,REF,CREF>
{
	typedef typename GRAPH::Edge Edge;
	typedef DenseGraphItemReferenceMap<GRAPH,Edge,T,REF,CREF> DenseGraphItemReferenceMapType;
	public:
		DenseEdgeReferenceMap()
		: DenseGraphItemReferenceMapType(){
		}
		DenseEdgeReferenceMap(const GRAPH & g)
		: DenseGraphItemReferenceMapType(g){
		}
		DenseEdgeReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
		: DenseGraphItemReferenceMapType(g,value){
		}
};

template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseArcReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Arc,T,REF,CREF>
{
	typedef typename GRAPH::Arc Arc;
	typedef DenseGraphItemReferenceMap<GRAPH,Arc,T,REF,CREF> DenseGraphItemReferenceMapType;
	public:
		DenseArcReferenceMap()
		: DenseGraphItemReferenceMapType(){
		}
		DenseArcReferenceMap(const GRAPH & g)
		: DenseGraphItemReferenceMapType(g){
		}
		DenseArcReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
		: DenseGraphItemReferenceMapType(g,value){
		}
};

}

#endif //VIGRA_DENSE_REFERENCE_MAP_HXX 