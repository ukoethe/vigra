#ifndef VIGRA_PYTHON_MAP
#define VIGRA_PYTHON_MAP

#include <boost/python.hpp>

namespace vigra{

	template<class T,class MERGE_GRAPH>
	class PythonGraphMap  {
	private:
        PythonGraphMap();                               			// non empty-construction
        PythonGraphMap( const PythonGraphMap& other );      	// non construction-copyable
        PythonGraphMap& operator=( const PythonGraphMap& ); 	// non copyable
	public:
		typedef MERGE_GRAPH MergeGraphType;
		typedef UInt32 LabelType;
		typedef T value_type;

		void mergeEdges(const LabelType a,const LabelType b){
			object_.attr("mergeEdges")(a,b);
		} 
		void mergeNodes(const LabelType a,const LabelType b){
			object_.attr("mergeNodes")(a,b);
		} 
		void eraseEdge(const LabelType label){
			object_.attr("eraseEdge")(label);
		}

		template<class CB>
		CB eraseEdgeCallback(){
			return  boost::bind(boost::mem_fn(&PythonGraphMap<T,MERGE_GRAPH>::eraseEdge), this , _1);
		}

		template<class CB>
		CB mergeEdgeCallback(){
			return  boost::bind(boost::mem_fn(&PythonGraphMap<T,MERGE_GRAPH>::mergeEdges), this , _1,_2);
		}

		template<class CB>
		CB mergeNodeCallback(){
			return  boost::bind(boost::mem_fn(&PythonGraphMap<T,MERGE_GRAPH>::mergeNodes), this , _1,_2);
		}

		PythonGraphMap(const MergeGraphType & mergeGraph,boost::python::object obj) 
		:	object_(obj),
			mergeGraph_(mergeGraph){

		}

	private:
		boost::python::object object_;
		const MergeGraphType & mergeGraph_;
	};

}

#endif