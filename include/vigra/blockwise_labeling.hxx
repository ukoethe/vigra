#ifndef VIGRA_BLOCKWISE_LABELING_HXX_
#define VIGRA_BLOCKWISE_LABELING_HXX_

#include <vigra/multi_gridgraph.hxx>
#include <vigra/multi_labeling.hxx>
#include <vigra/union_find.hxx>

#include <vigra/visit_border.hxx>
#include <vigra/blockify.hxx>

namespace vigra
{

namespace blockwise_labeling_detail
{

template <class Equal, class Label>
struct BorderVisitor
{
    Label u_label_offset;
    Label v_label_offset;
    UnionFindArray<Label>* global_unions;
    Equal* equal;
    
    template <class Data>
    void operator()(const Data& u_data, Label& u_label, const Data& v_data, Label& v_label)
    {
        if((*equal)(u_data, v_data))
        {
            global_unions->makeUnion(u_label + u_label_offset, v_label + v_label_offset);       
        }
    }
};

template <class DataBlocksIterator, class LabelBlocksIterator, class Equal>
typename LabelBlocksIterator::value_type::value_type
blockwiseLabeling(DataBlocksIterator data_blocks_begin, DataBlocksIterator data_blocks_end,
                  LabelBlocksIterator label_blocks_begin, DataBlocksIterator label_blocks_end,
                  NeighborhoodType neighborhood, Equal equal,
                  bool with_background, typename DataBlocksIterator::value_type::value_type background_value)
{
    typedef typename LabelBlocksIterator::value_type::value_type Label;
    typedef typename DataBlocksIterator::shape_type Shape;

    Shape blocks_shape = data_blocks_begin.shape();
    vigra_precondition(blocks_shape == label_blocks_begin.shape(), "shapes of blocks of blocks do not match");

    static const unsigned int Dimensions = DataBlocksIterator::dimension + 1;
    MultiArray<Dimensions, Label> label_offsets(label_blocks_begin.shape());
    
    // mapping stage: label each block and save number of labels assigned in blocks before the current block in label_offsets
    Label unmerged_label_number;
    {
        DataBlocksIterator data_blocks_it = data_blocks_begin;
        LabelBlocksIterator label_blocks_it = label_blocks_begin;
        typename MultiArray<Dimensions, Label>::iterator offsets_it = label_offsets.begin();
        Label current_offset = 0;
        for( ; data_blocks_it != data_blocks_end; ++data_blocks_it, ++label_blocks_it, ++offsets_it)
        {
            vigra_assert(label_blocks_it != label_blocks_end && offsets_it != label_offsets.end(), "");
            *offsets_it = current_offset;
            if(with_background)
            {
                current_offset += 1 + labelMultiArrayWithBackground(*data_blocks_it, *label_blocks_it,
                                                                    neighborhood, background_value, equal);
            }
            else
            {
                current_offset += labelMultiArray(*data_blocks_it, *label_blocks_it,
                                                  neighborhood, equal);
            }
        }
        unmerged_label_number = current_offset;
        if(!with_background)
            ++unmerged_label_number;
    }
    // reduce stage: merge adjacent labels if the region overlaps
    UnionFindArray<Label> global_unions(unmerged_label_number);
    if(with_background)
    {
        // merge all labels that refer to background
        for(typename MultiArray<Dimensions, Label>::iterator offsets_it = label_offsets.begin();
                offsets_it != label_offsets.end();
                ++offsets_it)
        {
            global_unions.makeUnion(0, *offsets_it);
        }
    }
    
    typedef GridGraph<Dimensions, undirected_tag> Graph;
    typedef typename Graph::edge_iterator EdgeIterator;
    Graph blocks_graph(blocks_shape);
    for(EdgeIterator it = blocks_graph.get_edge_iterator(); it != blocks_graph.get_edge_end_iterator(); ++it)
    {
        Shape u = blocks_graph.u(*it);
        Shape v = blocks_graph.v(*it);
        Shape difference = v - u;
        
        Label u_label_offset = label_offsets[u];
        Label v_label_offset = label_offsets[v];

        BorderVisitor<Equal, Label> border_visitor;
        border_visitor.u_label_offset = u_label_offset;
        border_visitor.v_label_offset = v_label_offset;
        border_visitor.global_unions = &global_unions;
        border_visitor.equal = &equal;
        visitBorder(data_blocks_begin[u], label_blocks_begin[u],
                    data_blocks_begin[v], label_blocks_begin[v],
                    difference, neighborhood, border_visitor);
    }
    // write global labels in each block
    Label last_label = global_unions.makeContiguous();
    {
        LabelBlocksIterator label_blocks_it = label_blocks_begin;
        typename MultiArray<Dimensions, Label>::iterator offsets_it = label_offsets.begin();
        for( ; label_blocks_it != label_blocks_end; ++label_blocks_it, ++offsets_it)
        {
            vigra_assert(offsets_it != label_offsets.end(), "");
            typedef typename LabelBlocksIterator::reference LabelBlockRef;
            LabelBlockRef label_block = *label_blocks_it;
            typedef typename LabelBlocksIterator::value_type::iterator LabelIterator;
            
            for(LabelIterator label_it = label_block.begin(); label_it != label_block.end(); ++label_it)
            {
                *label_it = global_unions.findLabel(*label_it + *offsets_it);
            }
        }
    }
    return last_label; 
}


const MultiArrayIndex default_block_side_length = 128;

}

template <unsigned int N, class T, class S1,
                          class Label, class S2>
Label labelMultiArrayBlockwise(const MultiArrayView<N, T, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               NeighborhoodType neighborhood = DirectNeighborhood,
                               const typename MultiArrayView<N, T, S1>::difference_type& block_shape
                                       = typename MultiArrayView<N, T, S1>::difference_type(blockwise_labeling_detail::default_block_side_length))
{
    MultiArray<N, MultiArrayView<N, T, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);

    return blockwise_labeling_detail::blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                                        label_blocks.begin(), label_blocks.end(),
                                                        neighborhood, std::equal_to<T>(),
                                                        false, T());
}

template <unsigned int N, class T, class S1,
                          class Label, class S2, class Equal>
Label labelMultiArrayBlockwise(const MultiArrayView<N, T, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               NeighborhoodType neighborhood,
                               const typename MultiArrayView<N, T, S1>::difference_type& block_shape,
                               Equal equal)

{
    MultiArray<N, MultiArrayView<N, T, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);

    return blockwise_labeling_detail::blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                                        label_blocks.begin(), label_blocks.end(),
                                                        neighborhood, equal,
                                                        false, T());
}

template <unsigned int N, class T, class S1,
                          class Label, class S2>
Label labelMultiArrayWithBackgroundBlockwise(const MultiArrayView<N, T, S1>& data,
                                             MultiArrayView<N, Label, S2> labels,
                                             NeighborhoodType neighborhood = DirectNeighborhood,
                                             T background_value = T(),
                                             const typename MultiArrayView<N, T, S1>::difference_type& block_shape
                                                     = typename MultiArrayView<N, T, S1>::difference_type(blockwise_labeling_detail::default_block_side_length))
{
    MultiArray<N, MultiArrayView<N, T, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);

    return blockwise_labeling_detail::blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                                        label_blocks.begin(), label_blocks.end(),
                                                        neighborhood, std::equal_to<T>(),
                                                        true, background_value); 
}

template <unsigned int N, class T, class S1,
                          class Label, class S2,
          class Equal>
Label labelMultiArrayWithBackgroundBlockwise(const MultiArrayView<N, T, S1>& data,
                                             MultiArrayView<N, Label, S2> labels,
                                             NeighborhoodType neighborhood,
                                             T background_value,
                                             const typename MultiArrayView<N, T, S1>::difference_type& block_shape,
                                             Equal equal)
{
    MultiArray<N, MultiArrayView<N, T, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);

    return blockwise_labeling_detail::blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                                        label_blocks.begin(), label_blocks.end(),
                                                        neighborhood, equal,
                                                        true, background_value); 
}

}

#endif

