#ifndef VIGRA_DEGREE_MAP_HXX
#define VIGRA_DEGREE_MAP_HXX

/*vigra*/ 
#include <vigra/graphs.hxx>

namespace vigra{
namespace detail{

	template<class GRAPH>
	class DegreeMap{

	public:
		typedef GRAPH Graph;
		typedef typename Graph::Node Node;
		typedef typename Graph::index_type index_type;

		typedef Node Key;
		typedef index_type Value;
		typedef index_type ConstReference;
		

		DegreeMap(const Graph graph)
		: graph_(graph){
		}
		index_type operator[](const Node & node)const{
			// use the boost free function
			// => will work with all our graphs
			return boost::degree(node,*graph_);
		}
		index_type operator[](const Node & node){
			// use the boost free function
			// => will work with all our graphs
			return boost::degree(node,*graph_);
		}


	private:

		const Graph * graph_;
	};	
}


#endif //VIGRA_DEGREE_MAP_HXX 