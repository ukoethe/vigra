#ifndef VIGRA_MERGE_GRAPH_INVALID_HXX
#define VIGRA_MERGE_GRAPH_INVALID_HXX

#include <vigra/graphs.hxx>

//#include "merge_graph.hxx"
#include "merge_graph_node.hxx"
#include "merge_graph_edge.hxx"
#include "merge_graph_iterators.hxx"
#include "iterable_partition.hxx"



namespace vigra {



namespace detail_merge_graph{
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
}


// merge graph item iterator (enge and node iterator)
template<class MERGE_GRAPH,class GRAPH_ITEM_TYPE>
inline bool operator == (const MergeGraphItemIterator<MERGE_GRAPH,GRAPH_ITEM_TYPE> & iter,const lemon::Invalid & iv){
    return iter.isEnd();
}
template<class MERGE_GRAPH,class GRAPH_ITEM_TYPE>
inline bool operator == (const lemon::Invalid & iv , const MergeGraphItemIterator<MERGE_GRAPH,GRAPH_ITEM_TYPE> & iter){
    return iter.isEnd();
}

template<class MERGE_GRAPH,class GRAPH_ITEM_TYPE>
inline bool operator != (const MergeGraphItemIterator<MERGE_GRAPH,GRAPH_ITEM_TYPE> & iter,const lemon::Invalid & iv){
    return !iter.isEnd();
}
template<class MERGE_GRAPH,class GRAPH_ITEM_TYPE>
inline bool operator != (const lemon::Invalid & iv , const MergeGraphItemIterator<MERGE_GRAPH,GRAPH_ITEM_TYPE> & iter){
    return !iter.isEnd();
}





} // namespace lemon


#endif // VIGRA_MERGE_GRAPH_INVALID_HXX