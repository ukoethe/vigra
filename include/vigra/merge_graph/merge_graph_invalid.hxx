#ifndef VIGRA_MERGE_GRAPH_INVALID_HXX
#define VIGRA_MERGE_GRAPH_INVALID_HXX

#include <vigra/graphs.hxx>

//#include "merge_graph.hxx"
#include "merge_graph_node.hxx"
#include "merge_graph_edge.hxx"
#include "merge_graph_iterators.hxx"
#include "iterable_partition.hxx"



namespace vigra {



namespace merge_graph_detail{
	// Iterable partiton
	template<class T>
	inline bool operator == (const ConstRepIter<T> & iter,const lemon::Invalid & iv){
	    return iter.isEnd();
	}
	template<class T>
	inline bool operator == (const lemon::Invalid & iv , const ConstRepIter<T> & iter){
	    return iter.isEnd();
	}

	template<class T>
	inline bool operator != (const ConstRepIter<T> & iter,const lemon::Invalid & iv){
	    return !iter.isEnd();
	}
	template<class T>
	inline bool operator != (const lemon::Invalid & iv , const ConstRepIter<T> & iter){
	    return !iter.isEnd();
	}


	template<class ITER>
	inline bool operator == (const IsEndIter<ITER> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class ITER>
	inline bool operator == (const lemon::Invalid & iv,const IsEndIter<ITER> & iter){
    	return iter.isEnd();
	}
	template<class ITER>
	inline bool operator != (const IsEndIter<ITER> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class ITER>
	inline bool operator != (const lemon::Invalid & iv,const IsEndIter<ITER> & iter){
    	return !iter.isEnd();
	}




	template<class FILTER,class ITER>
	inline bool operator == (const FilterIter<FILTER,ITER> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator == (const lemon::Invalid & iv,const FilterIter<FILTER,ITER> & iter){
    	return iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator != (const FilterIter<FILTER,ITER> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator != (const lemon::Invalid & iv,const FilterIter<FILTER,ITER> & iter){
    	return !iter.isEnd();
	}


	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator == (const TransformIter<TRANSFORM,ITER,REF,VAL> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator == (const lemon::Invalid & iv,const TransformIter<TRANSFORM,ITER,REF,VAL> & iter){
    	return iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator != (const TransformIter<TRANSFORM,ITER,REF,VAL> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator != (const lemon::Invalid & iv,const TransformIter<TRANSFORM,ITER,REF,VAL> & iter){
    	return !iter.isEnd();
	}

}







} // namespace lemon


#endif // VIGRA_MERGE_GRAPH_INVALID_HXX