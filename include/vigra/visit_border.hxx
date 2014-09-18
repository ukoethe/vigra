/************************************************************************/
/*                                                                      */
/*     Copyright 2013-2014 by Martin Bidlingmaier and Ullrich Koethe    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_VISIT_BORDER_HXX
#define VIGRA_VISIT_BORDER_HXX

#include "multi_array.hxx"

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
                     const Shape& block_difference, NeighborhoodType neighborhood, Visitor visitor)
    {
        static const unsigned int D = K - 1;
        typedef visit_border_impl<D> next;
        
        if(block_difference[D] == -1)
        {
            MultiArrayIndex last = v_data.shape(D) - 1;
            next::exec(u_data.bindAt(D, 0), u_labels.bindAt(D, 0),
                       v_data.bindAt(D, last), v_labels.bindAt(D, last),
                       block_difference, neighborhood, visitor);
        }
        else if(block_difference[D] == 1)
        {
            MultiArrayIndex last = u_data.shape(D) - 1;
            next::exec(u_data.bindAt(D, last), u_labels.bindAt(D, last),
                       v_data.bindAt(D, 0), v_labels.bindAt(D, 0),
                       block_difference, neighborhood, visitor);
        }
        else if(block_difference[D] == 0)
        {
            next::exec(u_data, u_labels, v_data, v_labels, block_difference, neighborhood, visitor);
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
                     const Shape& block_difference, NeighborhoodType neighborhood, Visitor visitor)
    {
        visitor(u_data(0), u_labels(0), v_data(0), v_labels(0), block_difference);
    }
    template <unsigned int N, class Data, class S1,
                              class Label, class S2,
              class Shape, class Visitor>
    static void exec(const MultiArrayView<N, Data, S1>& u_data, MultiArrayView<N, Label, S2> u_labels,
                     const MultiArrayView<N, Data, S1>& v_data, MultiArrayView<N, Label, S2> v_labels,
                     const Shape& block_difference, NeighborhoodType neighborhood, Visitor visitor)
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
                visitor(*u_data_it, *u_labels_it, *v_data_it, *v_labels_it, block_difference);
            }
        } 
        else if(neighborhood == IndirectNeighborhood)
        {
            typedef GridGraph<N, undirected_tag> Graph;
            typedef typename Graph::NodeIt GraphScanner;
            typedef typename Graph::OutArcIt NeighborIterator;
            
            static const int global_dim_number = Shape::static_size;
            TinyVector<unsigned int, N> dim_mapping; // mapping of every local dimension to their actual global dimension indices
            int local_dims_pos = 0;
            int global_dims_pos = 0;
            for( ; global_dims_pos != global_dim_number; ++global_dims_pos)
            {
                if(block_difference[global_dims_pos] == 0)
                {
                    vigra_assert(local_dims_pos != N, "");
                    dim_mapping[local_dims_pos] = global_dims_pos;
                    ++local_dims_pos;
                }
            }
            vigra_assert(local_dims_pos == N, "");

            Graph graph(u_data.shape(), neighborhood);
            Shape pixel_difference = block_difference;
            for(GraphScanner node(graph); node != lemon::INVALID; ++node)
            {
                // compare neighbors that have have equal coordinates in all unbound dimensions
                // their pixel-level difference is exactly block_difference
                visitor(u_data[*node], u_labels[*node], v_data[*node], v_labels[*node], block_difference);
                // now let unbound dimensions vary
                for(NeighborIterator arc(graph, *node); arc != lemon::INVALID; ++arc)
                {
                    for(int i = 0; i != N; ++i)
                        pixel_difference[dim_mapping[i]] = graph.target(*arc)[i] - (*node)[i];
                    visitor(u_data[*node], u_labels[*node], v_data[graph.target(*arc)], v_labels[graph.target(*arc)], pixel_difference);
                }
            }
        }
    }
};

} // namespace visit_border_detail

template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Shape, class Visitor>
inline void 
visitBorder(const MultiArrayView<N, Data, S1>& u_data, MultiArrayView<N, Label, S2>& u_labels,
            const MultiArrayView<N, Data, S1>& v_data, MultiArrayView<N, Label, S2>& v_labels,
            const Shape& difference, NeighborhoodType neighborhood, Visitor visitor)
{
    vigra_precondition(u_data.shape() == u_labels.shape() && v_data.shape() == v_labels.shape(),
                       "differing block shapes");
    visit_border_detail::visit_border_impl<N>::exec(u_data, u_labels,
                                                    v_data, v_labels,
                                                    difference, neighborhood, visitor);
}

} // namespace vigra

#endif // VIGRA_VISIT_BORDER_HXX
