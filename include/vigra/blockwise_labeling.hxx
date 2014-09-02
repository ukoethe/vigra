#ifndef VIGRA_BLOCKWISE_LABELING_HXX_
#define VIGRA_BLOCKWISE_LABELING_HXX_

#include <vigra/multi_gridgraph.hxx>
#include <vigra/multi_labeling.hxx>

namespace vigra
{

namespace blockwise_labeling_detail
{

template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Visitor>
void visitBorder(const MultiArrayView<N, Data, S1>& data, const MultiArrayView<N, Label, S2>& labels, Visitor visitor)
{

}



template <class DataBlocksIterator, class LabelBlocksIterator, class Equal, class BlocksShape>
typename LabelBlocksIterator::value_type::value_type
blockwiseLabeling(DataBlocksIterator data_blocks_begin, DataBlocksIterator data_blocks_end,
                  LabelBlocksIterator label_blocks_begin, DataBlocksIterator label_blocks_end,
                  NeighborhoodType neighborhood, Equal equal,
                  BlocksShape blocks_shape) // FIXME
{
    typedef typename LabelBlocksIterator::value_type::value_type Label;

    #ifdef VIGRA_CHECK_BOUNDS
        vigra_precondition(label_blocks_it != label_blocks_end, "less label blocks than data blocks");
    #endif

    std::vector<Label> label_offsets(data_blocks_end - data_blocks_begin);
    
    // mapping stage: label each block and save number of labels assigned in blocks before the current block in label_offsets
    DataBlocksIterator data_blocks_it = data_blocks_begin;
    LabelBlocksIterator label_blocks_it = label_blocks_begin;
    typename std::vector<Label>::iterator offsets_it = label_offsets.begin();
    Label current_offset = 0;
    for( ; data_blocks_it != data_blocks_end; ++data_blocks_it, ++label_blocks_it, ++offsets_it)
    {
        *offsets_it = current_offset;
        current_offset += labelMultiArray(*data_blocks_it, *label_blocks_it, neighborhood, equal);
    }

    // reduce stage: merge adjacent labels if the region overlaps
    static const unsigned int Dimensions = DataBlocksIterator::dimension + 1;
    typedef GridGraph<Dimensions, undirected_tag> Graph;
    typedef typename Graph::edge_iterator EdgeIterator;
    
    Graph blocks_graph(blocks_shape);
    for(EdgeIterator it = blocks_graph.get_edge_iterator(); it != blocks_graph.get_edge_end_iterator(); ++it)
    {
        BlocksShape u = blocks_graph.u(*it);
        BlocksShape v = blocks_graph.v(*it);
        BlocksShape difference = v - u;

        typedef typename DataBlocksIterator::value_type::reference DataRef;
        typedef typename LabelBlocksIterator::value_type::reference LabelRef;
        
        struct BorderVisitor
        {
            void operator()(DataRef u_data, LabelRef u_label, DataRef v_data, DataRef v_label)
            {

            }
        } border_visitor;
        /*
        visitBorder(data_blocks_begin[u], data_blocks_begin[u],
                    data_blocks_begin[v], data_blocks_begin[v],
                    difference, neighborhood, border_visitor);                
        */
    }

    

    return 0; // FIXME
}

}

}


#endif

