#ifndef LABEL_MULTI_ARRAY_BLOCKWISE_HXX_
#define LABEL_MULTI_ARRAY_BLOCKWISE_HXX_

#include "label_multi_range.hxx"
#include "visit_neighbors.hxx"
#include "visit_border.hxx"

#include <vigra/range/multi_range_traits.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/range/multi_length.hxx>


namespace vigra
{

namespace label_multi_array_blockwise_detail
{

template <class LabelBlocksRange>
struct label_type
{
    typedef typename multi_range_traits::decayed_value_type<
        typename multi_range_traits::decayed_value_type<LabelBlocksRange>::type>::type type;
};

}

template <class DataBlocksRange, class LabelBlocksRange, class Equal>
typename label_multi_array_blockwise_detail::label_type<LabelBlocksRange>::type
label_multi_range_blockwise(DataBlocksRange data_blocks, LabelBlocksRange label_blocks,
                            NeighborhoodType neighborhood, Equal equal)
{
    using namespace label_multi_array_blockwise_detail;
    typedef typename multi_range_traits::multi_size_type<DataBlocksRange>::type Shape;
    
    Shape blocks_shape = multi_length(data_blocks);
    vigra_assert(blocks_shape == multi_length(label_blocks), "");
    
    typedef typename label_type<LabelBlocksRange>::type Label;
    enum{ N = DataBlocksRange::dimension };

    MultiArray<N, Label> label_offsets(blocks_shape);
    
    auto labels_offset_blocks = multi_zip(label_blocks, multi_range(label_offsets));
    auto blocks = multi_zip(data_blocks,  labels_offset_blocks);
    typedef typename multi_range_traits::value_type<
        decltype(blocks)>::type Block; 
    // BlockType is something like 
    // std::pair<DataBlock::actual_value_type, std::pair<LabelBlock::actual_value_type, LabelType>>

    for_each_flattened(blocks, [&](Block block)
    {
        auto& data_block = block.first;
        auto& label_block = block.second.first;
        auto& label_offset = block.second.second;
        label_offset = labelMultiRange(data_block, label_block, neighborhood, equal);
    });

    Label current_offset = 0;
    for_each_flattened(multi_range(label_offsets), [&current_offset](Label& offset)
    {
        Label current_offset_temp = current_offset;
        current_offset += offset;
        offset = current_offset_temp;
    });

    Label last_unmerged_label = current_offset;

    UnionFindArray<Label> global_unions(last_unmerged_label + 1);
    
    visitNeighbors(blocks, neighborhood, [&global_unions, &neighborhood, &equal](Block u, Block v, Shape difference)
    {
        auto& u_data_block = u.first;
        auto& u_local_label_block = u.second.first;
        auto& u_label_offset = u.second.second;
        
        auto& v_data_block = v.first;
        auto& v_local_label_block = v.second.first;
        auto& v_label_offset = v.second.second;


        auto u_label_block = multi_map(u_local_label_block, [&u_label_offset](Label& local_label)
        {
            return local_label + u_label_offset;
        });
        auto v_label_block = multi_map(v_local_label_block, [&v_label_offset](Label& local_label)
        {
            return local_label + v_label_offset;
        });
        
        auto u_block = multi_zip(u_data_block, u_label_block);
        auto v_block = multi_zip(v_data_block , v_label_block);
        
        typedef typename multi_range_traits::value_type<decltype(u_block)>::type DataLabelPair;
        
        visitBorder(u_block, v_block, difference, neighborhood, [&global_unions, &equal](DataLabelPair a, DataLabelPair b)
        {
            if(equal(a.first, b.first))
            {
                global_unions.makeUnion(a.second, b.second);
            }
        });
    });

    Label last_label = global_unions.makeContiguous();

    for_each_flattened(blocks, [&global_unions](Block block)
    {
        auto& local_label_block = block.second.first;
        auto& label_offset = block.second.second;
        
        typedef typename multi_range_traits::value_type<
            typename std::decay<decltype(local_label_block)>::type>::type Label;

        for_each_flattened(local_label_block, [&global_unions, &label_offset](Label label)
        {
            label = global_unions.findLabel(label + label_offset);
        });
    });
    
    return last_label;
}
                                 
}

#endif

