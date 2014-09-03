#ifndef VIGRA_VISIT_BORDER_HXX_
#define VIGRA_VISIT_BORDER_HXX_

#include <vigra/multi_array.hxx>

namespace vigra
{

namespace visit_border_detail
{

template <unsigned int K>
struct visit_border_impl
{
    template <unsigned int N, class Data, class S1,
                              class Label, class S2,
              class Shape, class Visitor>
    static void exec(const MultiArrayView<N, Data, S1>& u_data, MultiArrayView<N, Label, S2> u_labels,
                     const MultiArrayView<N, Data, S1>& v_data, MultiArrayView<N, Label, S2> v_labels,
                     const Shape& difference, NeighborhoodType neighborhood, Visitor visitor)
    {
        static const unsigned int D = K - 1;
        typedef visit_border_impl<D> next;
        
        if(difference[D] == -1)
        {
            MultiArrayIndex data_last = v_data.shape(D) - 1;
            MultiArrayIndex labels_last = v_labels.shape(D) - 1;
            next::exec(u_data.bindAt(D, 0), u_labels.bindAt(D, 0),
                       v_data.bindAt(D, data_last), v_labels.bindAt(D, labels_last),
                       difference, neighborhood, visitor);
        }
        else if(difference[D] == 1)
        {
            MultiArrayIndex data_last = u_data.shape(D) - 1;
            MultiArrayIndex labels_last = u_labels.shape(D) - 1;
            next::exec(u_data.bindAt(D, data_last), u_labels.bindAt(D, labels_last),
                       v_data.bindAt(D, 0), v_labels.bindAt(D, 0),
                       difference, neighborhood, visitor);
        }
        else if(difference[D] == 0)
        {
            vigra_precondition(u_data.shape(D) == u_labels.shape(D) &&
                               v_data.shape(D) == v_labels.shape(D) &&
                               u_data.shape(D) == v_data.shape(D),
                               "block shape mismatch");
            next::exec(u_data, u_labels, v_data, v_labels, difference, neighborhood, visitor);
        }
        else
        {
            vigra_precondition(false, "invalid block difference");
        }
    }
};

template <>
struct visit_border_impl<0>
{
    template <class Data, class S1,
              class Label, class S2,
              class Shape, class Visitor>
    static void exec(const MultiArrayView<0, Data, S1>& u_data, MultiArrayView<0, Label, S2> u_labels,
                     const MultiArrayView<0, Data, S1>& v_data, MultiArrayView<0, Label, S2> v_labels,
                     const Shape& difference, NeighborhoodType neighborhood, Visitor visitor)
    {
        visitor(u_data(0), u_labels(0), v_data(0), v_labels(0));
    }
    template <unsigned int N, class Data, class S1,
                              class Label, class S2,
              class Shape, class Visitor>
    static void exec(const MultiArrayView<N, Data, S1>& u_data, MultiArrayView<N, Label, S2> u_labels,
                     const MultiArrayView<N, Data, S1>& v_data, MultiArrayView<N, Label, S2> v_labels,
                     const Shape& difference, NeighborhoodType neighborhood, Visitor visitor)
    {
        if(neighborhood == DirectNeighborhood)
        {
            typedef typename MultiArrayView<N, Data, S1>::const_iterator DataIterator;
            typedef typename MultiArrayView<N, Label, S2>::iterator LabelsIterator;

            DataIterator u_data_it = u_data.begin();
            LabelsIterator u_labels_it = u_labels.begin();

            DataIterator v_data_it = v_data.begin();
            LabelsIterator v_labels_it = v_labels.begin();

            for( ; u_data_it != u_data.end(); ++u_data_it, ++u_labels_it, ++v_data_it, ++v_labels_it)
            {
                visitor(*u_data_it, *u_labels_it, *v_data_it, *v_labels_it);
            }
        } 
        else if(neighborhood == IndirectNeighborhood)
        {
            typedef GridGraph<N, undirected_tag> Graph;
            typedef typename Graph::NodeIt GraphScanner;
            typedef typename Graph::OutArcIt NeighborIterator;
            
            Graph graph(u_data.shape(), neighborhood);
            for(GraphScanner node(graph); node != lemon::INVALID; ++node)
            {
                visitor(u_data[*node], u_labels[*node], v_data[*node], v_labels[*node]);
                for(NeighborIterator arc(graph, *node); arc != lemon::INVALID; ++arc)
                {
                    visitor(u_data[*node], u_labels[*node], v_data[graph.target(*arc)], v_labels[graph.target(*arc)]);
                }
            }
        }
    }
};

}
template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Shape, class Visitor>
void visitBorder(const MultiArrayView<N, Data, S1>& u_data, MultiArrayView<N, Label, S2>& u_labels,
                 const MultiArrayView<N, Data, S1>& v_data, MultiArrayView<N, Label, S2>& v_labels,
                 const Shape& difference, NeighborhoodType neighborhood, Visitor visitor)
{
    visit_border_detail::visit_border_impl<N>::exec(u_data, u_labels,
                                                    v_data, v_labels,
                                                    difference, neighborhood, visitor);
}

}

#endif

