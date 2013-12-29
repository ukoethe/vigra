#ifndef VIGRA_MULTI_ARRAY_MAP_HXX
#define VIGRA_MULTI_ARRAY_MAP_HXX

/* boost */
#include <boost/function.hpp>


/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/accumulator.hxx>

/* vigra - merge graph */
#include <vigra/merge_graph/min_indexed_pq.hxx>


namespace vigra{
namespace view_maps {


	template<unsigned int DIM,class T>
	class SumMap;


	template<class T>
	class SumMap<0,T> : public MultiArrayView<1,T> {
	private:
        SumMap();                               // non empty-construction
        SumMap( const SumMap& other );      	// non construction-copyable
        SumMap& operator=( const SumMap& ); 	// non copyable
	public:
		typedef MultiArrayView<1,T> ArrayViewType;
		typedef UInt32 LabelType;
		typedef T value_type;

		//const value_type & operator[](const LabelType label)const{
		//	return this->operator(label);
		//}

		void merge(const LabelType a,const LabelType b){
			this->operator()(a)+=this->operator()(b);
		} 

		template<class CB>
		CB mergeCallback(){
			CB cb;
			cb = boost::bind(boost::mem_fn(&SumMap<0,T>::merge), this , _1,_2);
			return cb;
		}

		SumMap(ArrayViewType & array) 
		:	ArrayViewType(array){
		}
	};

	template<unsigned int DIM, class T>
	class SumMap : public MultiArrayView<DIM+1,T> {
	private:
        SumMap();                               // non empty-construction
        SumMap( const SumMap& other );      	// non construction-copyable
        SumMap& operator=( const SumMap& ); 	// non copyable
	public:
		typedef MultiArrayView<DIM+1,T> ArrayViewType;
		typedef UInt32 LabelType;
		typedef T value_type;

		//const value_type & operator[](const LabelType label)const{
		//	return this->operator(label);
		//}

		void merge(const LabelType a,const LabelType b){
			this->bindInner(a)+=this->bindInner(b);
		} 

		template<class CB>
		CB mergeCallback(){
			CB cb;
			cb = boost::bind(boost::mem_fn(&SumMap<DIM,T>::merge), this , _1,_2);
			return cb;
		}

		SumMap(ArrayViewType & array) 
		:	ArrayViewType(array){
		}
	};



	template<unsigned int DIM,class T,class WEIGHT_MAP>
	class WeightedMeanMap;


	template<class T,class WEIGHT_MAP>
	class WeightedMeanMap<0,T,WEIGHT_MAP> : public MultiArrayView<1,T> {
	private:
        WeightedMeanMap();                               		// non empty-construction
        WeightedMeanMap( const WeightedMeanMap& other );      	// non construction-copyable
        WeightedMeanMap& operator=( const WeightedMeanMap& ); 	// non copyable
	public:
		typedef WEIGHT_MAP WeightMapType;
		typedef MultiArrayView<1,T> ArrayViewType;
		typedef UInt32 LabelType;
		typedef T value_type;

		void merge(const LabelType a,const LabelType b){
			const T va=this->operator()(a);
			const T vb=this->operator()(b);
			const T wa=weightMap_(a);
			const T wb=weightMap_(b);
			this->operator()(a)=(va*wa + vb*wb)/(wa+wb);
		} 


		template<class CB>
		CB mergeCallback(){
			return  boost::bind(boost::mem_fn(&WeightedMeanMap<0,T,WEIGHT_MAP>::merge), this , _1,_2);
		}

		WeightedMeanMap(ArrayViewType & array,const WeightMapType & weightMap) 
		:	ArrayViewType(array),
			weightMap_(weightMap){
		}
	private:
		const WeightMapType & weightMap_;
	};


	template<unsigned int DIM,class T,class WEIGHT_MAP>
	class WeightedMeanMap : public MultiArrayView<DIM+1,T> {
	private:
        WeightedMeanMap();                               		// non empty-construction
        WeightedMeanMap( const WeightedMeanMap& other );      	// non construction-copyable
        WeightedMeanMap& operator=( const WeightedMeanMap& ); 	// non copyable
	public:
		typedef WEIGHT_MAP WeightMapType;
		typedef MultiArrayView<DIM+1,T> ArrayViewType;
		typedef UInt32 LabelType;
		typedef T value_type;

		void merge(const LabelType a,const LabelType b){
			const T wa=weightMap_(a);
			const T wb=weightMap_(b);
			vigra::MultiArrayView<DIM,T> va = this->bindInner(a);
			vigra::MultiArrayView<DIM,T> 	 vb = this->bindInner(b);
			//a=(a*wa+b*wb)
			va*=wa;
			vb*=wb;
			va+=vb;
			va/=(wa+wb);
			vb/=wb;
		} 


		template<class CB>
		CB mergeCallback(){
			return  boost::bind(boost::mem_fn(&WeightedMeanMap<DIM,T,WEIGHT_MAP>::merge), this , _1,_2);
		}

		WeightedMeanMap(ArrayViewType & array,const WeightMapType & weightMap) 
		:	ArrayViewType(array),
			weightMap_(weightMap){
		}
	private:
		const WeightMapType & weightMap_;
	};





	template<class T,class MERGE_GRAPH,class EDGE_MAP,class NODE_MAP>
	class MinWeightEdgeMap : public MultiArrayView<1,T> {
	private:
        MinWeightEdgeMap();                               			// non empty-construction
        MinWeightEdgeMap( const MinWeightEdgeMap& other );      	// non construction-copyable
        MinWeightEdgeMap& operator=( const MinWeightEdgeMap& ); 	// non copyable
	public:
		typedef MERGE_GRAPH MergeGraphType;
		typedef EDGE_MAP    EdgeMapType;
		typedef NODE_MAP	NodeMapType;
		typedef MultiArrayView<1,T> ArrayViewType;
		typedef UInt32 LabelType;
		typedef T value_type;

		void merge(const LabelType a,const LabelType b){

		} 
		void eraseEdge(const LabelType label){

		}


		template<class CB>
		CB mergeCallback(){
			
		}

		MinWeightEdgeMap(const MergeGraphType & mergeGraph,ArrayViewType & array,const EdgeMapType & edgeMap,const NodeMapType & nodeMap) 
		:	ArrayViewType(array),
			mergeGraph_(mergeGraph),
			edgeMap_(edgeMap),
			nodeMap_(nodeMap){

		}
	private:
		//T 


		const MergeGraphType & mergeGraph_;
		const EdgeMapType 	 & edgeMap_;
		const NodeMapType    & nodeMap_;
	};


} // end namespace view_maps
} // end namespace vigra

#endif // VIGRA_MULTI_ARRAY_MAP_HXX