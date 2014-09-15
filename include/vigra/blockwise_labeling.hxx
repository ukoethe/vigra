#ifndef VIGRA_BLOCKWISE_LABELING_HXX_
#define VIGRA_BLOCKWISE_LABELING_HXX_

#include <vigra/multi_gridgraph.hxx>
#include <vigra/multi_labeling.hxx>
#include <vigra/union_find.hxx>
#include <vigra/multi_array_chunked.hxx>

#include <vigra/visit_border.hxx>
#include <vigra/blockify.hxx>

namespace vigra
{

class LabelOptions;

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

template <class DataBlocksIterator, class LabelBlocksIterator, class Equal, class Mapping>
typename LabelBlocksIterator::value_type::value_type
blockwiseLabeling(DataBlocksIterator data_blocks_begin, DataBlocksIterator data_blocks_end,
                  LabelBlocksIterator label_blocks_begin, LabelBlocksIterator label_blocks_end,
                  NeighborhoodType neighborhood, Equal equal,
                  const typename DataBlocksIterator::value_type::value_type* background_value,
                  Mapping& mapping)
{
    typedef typename LabelBlocksIterator::value_type::value_type Label;
    typedef typename DataBlocksIterator::shape_type Shape;

    Shape blocks_shape = data_blocks_begin.shape();
    vigra_precondition(blocks_shape == label_blocks_begin.shape() &&
                       blocks_shape == mapping.shape(),
                       "shapes of blocks of blocks do not match");

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
            if(background_value)
            {
                current_offset += 1 + labelMultiArrayWithBackground(*data_blocks_it, *label_blocks_it,
                                                                    neighborhood, *background_value, equal);
            }
            else
            {
                current_offset += labelMultiArray(*data_blocks_it, *label_blocks_it,
                                                  neighborhood, equal);
            }
        }
        unmerged_label_number = current_offset;
        if(!background_value)
            ++unmerged_label_number;
    }
    
    // reduce stage: merge adjacent labels if the region overlaps
    UnionFindArray<Label> global_unions(unmerged_label_number);
    if(background_value)
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
    Graph blocks_graph(blocks_shape, neighborhood);
    for(EdgeIterator it = blocks_graph.get_edge_iterator(); it != blocks_graph.get_edge_end_iterator(); ++it)
    {
        Shape u = blocks_graph.u(*it);
        Shape v = blocks_graph.v(*it);
        Shape difference = v - u;
        
        BorderVisitor<Equal, Label> border_visitor;
        border_visitor.u_label_offset = label_offsets[u];
        border_visitor.v_label_offset = label_offsets[v];
        border_visitor.global_unions = &global_unions;
        border_visitor.equal = &equal;
        visitBorder(data_blocks_begin[u], label_blocks_begin[u],
                    data_blocks_begin[v], label_blocks_begin[v],
                    difference, neighborhood, border_visitor);
    }

    // fill mapping (local labels) -> (global labels)
    Label last_label = global_unions.makeContiguous();
    {
        typename MultiArray<Dimensions, Label>::iterator offsets_it = label_offsets.begin();
        Label offset = *offsets_it;
        ++offsets_it;
        typename Mapping::iterator mapping_it = mapping.begin();
        for( ; offsets_it != label_offsets.end(); ++offsets_it, ++mapping_it)
        {
            mapping_it->clear();
            Label next_offset = *offsets_it;
            if(background_value)
            {
                for(Label current_label = offset; current_label != next_offset; ++current_label)
                {
                    mapping_it->push_back(global_unions.findLabel(current_label));
                }
            }
            else
            {
                mapping_it->push_back(0); // local labels start at 1
                for(Label current_label = offset + 1; current_label != next_offset + 1; ++current_label)
                {
                    mapping_it->push_back(global_unions.findLabel(current_label));
                }
            }

            offset = next_offset;
        }
        // last block:
        // instead of next_offset, use last_label+1 
        mapping_it->clear();
        if(background_value)
        {
            for(Label current_label = offset; current_label != unmerged_label_number; ++current_label)
            {
                mapping_it->push_back(global_unions.findLabel(current_label));    
            }
        }
        else
        {
            mapping_it->push_back(0);
            for(Label current_label = offset + 1; current_label != unmerged_label_number; ++current_label)
            {
                mapping_it->push_back(global_unions.findLabel(current_label));
            }
        }
    }
    return last_label; 
}


template <class LabelBlocksIterator, class MappingIterator>
void toGlobalLabels(LabelBlocksIterator label_blocks_begin, LabelBlocksIterator label_blocks_end,
                    MappingIterator mapping_begin, MappingIterator mapping_end)
{
    typedef typename LabelBlocksIterator::value_type LabelBlock;
    for( ; label_blocks_begin != label_blocks_end; ++label_blocks_begin, ++mapping_begin)
    {
        vigra_assert(mapping_begin != mapping_end, "");
        for(typename LabelBlock::iterator labels_it = label_blocks_begin->begin();
            labels_it != label_blocks_begin->end();
            ++labels_it)
        {
            vigra_assert(*labels_it < mapping_begin->size(), "");
            *labels_it = (*mapping_begin)[*labels_it];
        }
    }
}


const MultiArrayIndex default_block_side_length = 128;

template <class T>
const T* getBackground(const LabelOptions& options);
template <class T>
const T* getBlockShape(const LabelOptions& options);
NeighborhoodType getNeighborhood(const LabelOptions& options);

}

class LabelOptions
{
private:
    struct type_erasure_base
    {
        virtual ~type_erasure_base()
        {}
    };
    template <class T>
    struct type_erasure
      : type_erasure_base
    {
        type_erasure(const T& contained_obj)
          : obj(contained_obj)
        {}
        T obj;
    };
    
    std::auto_ptr<type_erasure_base> background_value_;
    std::auto_ptr<type_erasure_base> block_shape_;
    NeighborhoodType neighborhood_;
public:
    LabelOptions()
      : neighborhood_(DirectNeighborhood)
    {}
private:
    LabelOptions(const LabelOptions&); // deleted
    LabelOptions& operator=(const LabelOptions&); // deleted
public:
    template <class T>
    LabelOptions& background(const T& background_value)
    {
        vigra_precondition(background_value_.get() == 0, "background set twice");
        background_value_ = std::auto_ptr<type_erasure_base>(new type_erasure<T>(background_value));
        return *this;
    }
    template <int N>
    LabelOptions& blockShape(const TinyVector<MultiArrayIndex, N>& block_shape)
    {
        vigra_precondition(block_shape_.get() == 0, "block shape set twice");
        block_shape_ = std::auto_ptr<type_erasure_base>(new type_erasure<TinyVector<MultiArrayIndex, N> >(block_shape));
        return *this;
    }

    LabelOptions& neighborhood(NeighborhoodType n)
    {
        neighborhood_ = n;
        return *this;
    }

    
    template <class T>
    friend const T* blockwise_labeling_detail::getBackground(const LabelOptions& options);
    template <class T>
    friend const T* blockwise_labeling_detail::getBlockShape(const LabelOptions& options);
    friend NeighborhoodType blockwise_labeling_detail::getNeighborhood(const LabelOptions& options);
};

namespace blockwise_labeling_detail
{

template <class T>
const T* getBackground(const LabelOptions& options)
{
    if(options.background_value_.get() == 0)
        return 0;

    LabelOptions::type_erasure<T>* background = dynamic_cast<LabelOptions::type_erasure<T>*>(options.background_value_.get());
    vigra_precondition(background != 0, "background value type and data type do not match");
    return &background->obj;
}
template <class T>
const T* getBlockShape(const LabelOptions& options)
{
    if(options.block_shape_.get() == 0)
        return 0;

    LabelOptions::type_erasure<T>* block_shape = dynamic_cast<LabelOptions::type_erasure<T>*>(options.block_shape_.get());
    vigra_precondition(block_shape != 0, "shape type has invalid dimension");
    return &block_shape->obj;
}
NeighborhoodType getNeighborhood(const LabelOptions& options)
{
    return options.neighborhood_;
}

}

template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Equal, class S3>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               const LabelOptions& options,
                               Equal equal, MultiArrayView<N, std::vector<Label>, S3>& mapping)
{
    using namespace blockwise_labeling_detail;
    TinyVector<MultiArrayIndex, N> block_shape(default_block_side_length);
    const TinyVector<MultiArrayIndex, N>* options_block_shape = getBlockShape<TinyVector<MultiArrayIndex, N> >(options);
    if(options_block_shape)
        block_shape = *options_block_shape;
    const Data* background_value = getBackground<Data>(options);
    NeighborhoodType neighborhood = getNeighborhood(options);
    
    MultiArray<N, MultiArrayView<N, Data, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);
    return blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                             label_blocks.begin(), label_blocks.end(),
                             neighborhood, equal, background_value, mapping);
}
template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Equal>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels, const LabelOptions& options, Equal equal) {
    using namespace blockwise_labeling_detail;
    TinyVector<MultiArrayIndex, N> block_shape(default_block_side_length);
    const TinyVector<MultiArrayIndex, N>* options_block_shape = getBlockShape<TinyVector<MultiArrayIndex, N> >(options);
    if(options_block_shape)
        block_shape = *options_block_shape;
    const Data* background_value = getBackground<Data>(options);
    NeighborhoodType neighborhood = getNeighborhood(options);
    
    MultiArray<N, MultiArrayView<N, Data, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);
    MultiArray<N, std::vector<Label> > mapping(data_blocks.shape());
    Label last_label = blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                         label_blocks.begin(), label_blocks.end(),
                                         neighborhood, equal, background_value, mapping);

    // replace local labels by global labels
    toGlobalLabels(label_blocks.begin(), label_blocks.end(), mapping.begin(), mapping.end());
    return last_label;
}
template <unsigned int N, class Data, class S1,
                          class Label, class S2>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               const LabelOptions& options = LabelOptions())
{
    return labelMultiArrayBlockwise(data, labels, options, std::equal_to<Data>());
}



template <unsigned int N, class Data, class Label, class Equal, class S3>
Label labelMultiArrayBlockwise(const ChunkedArray<N, Data>& data,
                               ChunkedArray<N, Label>& labels,
                               const LabelOptions& options,
                               Equal equal, MultiArrayView<N, std::vector<Label>, S3>& mapping)
{    
    using namespace blockwise_labeling_detail;
    const TinyVector<MultiArrayIndex, N>* options_block_shape = getBlockShape<TinyVector<MultiArrayIndex, N> >(options);
    vigra_precondition(options_block_shape == 0, "block shape not supported for chunked arrays, uses chunk size per default");
    
    typedef typename ChunkedArray<N, Data>::shape_type Shape;

    const Data* background_value = getBackground<Data>(options);
    NeighborhoodType neighborhood = getNeighborhood(options);
    
    typedef typename ChunkedArray<N, Data>::chunk_const_iterator DataChunkIterator;
    typedef typename ChunkedArray<N, Label>::chunk_iterator LabelChunkIterator;

    DataChunkIterator data_chunks_begin = data.chunk_begin(Shape(0), data.shape());
    LabelChunkIterator label_chunks_begin = labels.chunk_begin(Shape(0), labels.shape());
    
    return blockwiseLabeling(data_chunks_begin, data_chunks_begin.getEndIterator(),
                             label_chunks_begin, label_chunks_begin.getEndIterator(),
                             neighborhood, equal, background_value, mapping);
}
template <unsigned int N, class Data, class Label, class Equal>
Label labelMultiArrayBlockwise(const ChunkedArray<N, Data>& data,
                               ChunkedArray<N, Label>& labels,
                               const LabelOptions& options, Equal equal)
{   
    using namespace blockwise_labeling_detail;
    MultiArray<N, std::vector<Label> > mapping(data.chunkArrayShape());
    Label result = labelMultiArrayBlockwise(data, labels, options, equal, mapping);
    typedef typename ChunkedArray<N, Data>::shape_type Shape;
    toGlobalLabels(labels.chunk_begin(Shape(0), data.shape()), labels.chunk_end(Shape(0), data.shape()), mapping.begin(), mapping.end());
    return result;
}
template <unsigned int N, class Data, class Label>
Label labelMultiArrayBlockwise(const ChunkedArray<N, Data>& data,
                               ChunkedArray<N, Label>& labels,
                               const LabelOptions& options = LabelOptions())
{
    return labelMultiArrayBlockwise(data, labels, options, std::equal_to<Data>());
}
 
    
}

#endif

