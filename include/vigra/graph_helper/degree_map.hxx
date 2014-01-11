#ifndef VIGRA_DEGREE_MAP_HXX
#define VIGRA_DEGREE_MAP_HXX


namespace vigra{
namespace detail{

	template<class GRAPH>
	class DegreeMap{

	public:
		typedef GRAPH Graph;
		typedef typename Graph::Node Node;
		typedef typename Graph::index_type index_type;

		DegreeMap(const Graph graph)
		: graph_(graph){
		}
		index_type operator[](const Node & node)const{
			return graph_.degree(node);
		}
		index_type operator[](const Node & node){
			return graph_.degree(node);
		}


	private:

		const Graph & graph_;
	};	
}
}


#endif //VIGRA_DEGREE_MAP_HXX 