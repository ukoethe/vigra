#ifndef VIIST_NEIGHBORS_HXX_
#define VISIT_NEIGHBORS_HXX_

#include "label_multi_range.hxx"

namespace vigra
{

template <class Range, class Visitor>
void visitNeighbors(Range range, NeighborhoodType neighborhood, Visitor visitor)
{
    enum{ N = Range::dimension };
    label_multi_range_detail::RangePropertyMap<Range> rpm(range); //FIXME: refactor RangePropertyMap into its own module

    typedef typename multi_range_traits::multi_size_type<Range>::type Shape;
    Shape shape = multi_length(range);

    typedef GridGraph<N, undirected_tag> Graph;
    typedef typename Graph::edge_iterator EdgeIterator;
    typedef typename Graph::Node Node;
    
    Graph graph(shape, neighborhood);
    for(EdgeIterator it = graph.get_edge_iterator(); it != graph.get_edge_end_iterator(); ++it)
    {
        Node u = graph.u(*it);
        Node v = graph.v(*it);
        visitor(rpm[u], rpm[v], v - u);
    }
}

}

#endif

