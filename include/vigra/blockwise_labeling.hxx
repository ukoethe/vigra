/************************************************************************/
/*                                                                      */
/*    Copyright 2013-2014 by Martin Bidlingmaier, Ullrich Koethe        */
/*                                                and Kevin Kiefer      */
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

#ifndef VIGRA_BLOCKWISE_LABELING_HXX
#define VIGRA_BLOCKWISE_LABELING_HXX

#include <algorithm>
#include <type_traits>

#include "threadpool.hxx"
#include "counting_iterator.hxx"
#include "multi_gridgraph.hxx"
#include "multi_labeling.hxx"
#include "multi_blockwise.hxx"
#include "union_find.hxx"
#include "multi_array_chunked.hxx"
#include "metaprogramming.hxx"

#include "visit_border.hxx"
#include "blockify.hxx"

namespace vigra
{
/** \addtogroup Labeling
*/
//@{

    /** Options object for labelMultiArrayBlockwise().

        It is simply a subclass of both \ref vigra::LabelOptions
        and \ref vigra::BlockwiseOptions. See there for
        detailed documentation.
    */
template <unsigned int N>
class BlockwiseLabelOptions
: public LabelOptions
, public BlockwiseOptions<N>
{
public:
    typedef typename BlockwiseOptions<N>::Shape Shape;

    // reimplement setter functions to allow chaining

    template <class T>
    BlockwiseLabelOptions<N>& ignoreBackgroundValue(const T& value)
    {
        LabelOptions::ignoreBackgroundValue(value);
        return *this;
    }

    BlockwiseLabelOptions<N> & neighborhood(NeighborhoodType n)
    {
        LabelOptions::neighborhood(n);
        return *this;
    }

    BlockwiseLabelOptions<N> & blockShape(const Shape & shape)
    {
        BlockwiseOptions<N>::blockShape(shape);
        return *this;
    }

    template <class T>
    BlockwiseLabelOptions<N> & blockShape(const TinyVector<T, N> & shape)
    {
        BlockwiseOptions<N>::blockShape(shape);
        return *this;
    }

    BlockwiseLabelOptions<N> & numThreads(const int n)
    {
        BlockwiseOptions<N>::numThreads(n);
        return *this;
    }
};

namespace blockwise_labeling_detail
{

template <class Equal, class Label>
struct BorderVisitor
{
    Label u_label_offset;
    Label v_label_offset;
    UnionFindArray<Label>* global_unions;
    Equal* equal;

    template <class Data, class Shape>
    void operator()(const Data& u_data, Label& u_label, const Data& v_data, Label& v_label, const Shape& diff)
    {
        if(labeling_equality::callEqual(*equal, u_data, v_data, diff))
        {
            global_unions->makeUnion(u_label + u_label_offset, v_label + v_label_offset);
        }
    }
};

// needed by MSVC
template <class LabelBlocksIterator>
struct BlockwiseLabelingResult
{
    typedef typename LabelBlocksIterator::value_type::value_type type;
};

template <class DataBlocksIterator, class LabelBlocksIterator,
          class Equal, class Mapping, unsigned int N>
typename BlockwiseLabelingResult<LabelBlocksIterator>::type
blockwiseLabeling(DataBlocksIterator data_blocks_begin, DataBlocksIterator data_blocks_end,
                  LabelBlocksIterator label_blocks_begin, LabelBlocksIterator label_blocks_end,
                  BlockwiseLabelOptions<N> const & options,
                  Equal equal,
                  Mapping& mapping)
{
    typedef typename LabelBlocksIterator::value_type::value_type Label;
    typedef typename DataBlocksIterator::shape_type Shape;

    Shape blocks_shape = data_blocks_begin.shape();
    vigra_precondition(blocks_shape == label_blocks_begin.shape() &&
                       blocks_shape == mapping.shape(),
                       "shapes of blocks of blocks do not match");
    vigra_precondition(std::distance(data_blocks_begin,data_blocks_end) == std::distance(label_blocks_begin,label_blocks_end),
                       "the sizes of input ranges are different");

    static const unsigned int Dimensions = DataBlocksIterator::dimension + 1;
    MultiArray<Dimensions, Label> label_offsets(label_blocks_begin.shape());

    bool has_background = options.hasBackgroundValue();

    // mapping stage: label each block and save number of labels assigned in blocks before the current block in label_offsets
    Label unmerged_label_number;
    {
        DataBlocksIterator data_blocks_it = data_blocks_begin;
        LabelBlocksIterator label_blocks_it = label_blocks_begin;
        typename MultiArray<Dimensions, Label>::iterator offsets_it = label_offsets.begin();
        Label current_offset = 0;
        // a la OPENMP_PRAGMA FOR

        auto d = std::distance(data_blocks_begin, data_blocks_end);


        std::vector<Label> nSeg(d);
        //std::vector<int> ids(d);
        //std::iota(ids.begin(), ids.end(), 0 );

        parallel_foreach(options.getNumThreads(), d,
            [&](const int /*threadId*/, const uint64_t i){
                Label resVal = labelMultiArray(data_blocks_it[i], label_blocks_it[i],
                                               options, equal);
                if(has_background) // FIXME: reversed condition?
                    ++resVal;
                nSeg[i] = resVal;
            }
        );

        for(int i=0; i<d;++i){
            offsets_it[i] = current_offset;
            current_offset+=nSeg[i];
        }


        unmerged_label_number = current_offset;
        if(!has_background)
            ++unmerged_label_number;
    }

    // reduce stage: merge adjacent labels if the region overlaps
    UnionFindArray<Label> global_unions(unmerged_label_number);
    if(has_background)
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
    Graph blocks_graph(blocks_shape, options.getNeighborhood());
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
                    difference, options.getNeighborhood(), border_visitor);
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
            if(has_background)
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
        if(has_background)
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


// Version for ChunkedArrays. Uses the commit/checkout API internally
template <unsigned int M, unsigned int N1, unsigned int N2,
          class Data, class Label, class TAG1, class TAG2, class Equal, class Mapping>
Label chunkedBlockwiseLabeling(const ChunkedArray<N1,Data,TAG1> & data,
                               ChunkedArray<N2,Label,TAG2> & labels,
                               const MultiBlocking<M, MultiArrayIndex> & blocking,
                               const BlockwiseLabelOptions<M> & options,
                               Equal equal,
                               Mapping& mapping)
{
    typedef TinyVector<MultiArrayIndex,M> Shape;
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;
    typedef typename Blocking::BlockIter BlockIter;

    Shape blocks_shape = data.template minimalShape<M>() / options.getBlockShape();
    Shape reminder = data.template minimalShape<M>() % options.getBlockShape();
    for (unsigned i = 0; i < M; ++i)
        if (reminder[i] != 0)
            blocks_shape[i] += 1;

    MultiArray<M, Label> label_offsets(blocks_shape);
    bool has_background = options.hasBackgroundValue();

    // mapping stage: label each block and save number of labels assigned in blocks before the current block in label_offsets
    Label unmerged_label_number;
    {
        BlockIter data_blocks_it = blocking.blockBegin();
        BlockIter label_blocks_it = blocking.blockBegin();
        typename MultiArray<M, Label>::iterator offsets_it = label_offsets.begin();

        Label current_offset = 0;
        auto d = std::distance(blocking.blockBegin(), blocking.blockEnd());
        std::vector<Label> nSeg(d);

        parallel_foreach(options.getNumThreads(), d,
            [&](const int /*threadId*/, const uint64_t i){
                auto block = blocking.blockBegin()[i];

                MultiArray<M,Data> data_cache(block.size());
                MultiArray<M,Label> label_cache(block.size());

                data.checkoutSubarray(block.begin(), data_cache);

                Label resVal = labelMultiArray(data_cache, label_cache, options, equal);

                labels.commitSubarray(block.begin(), label_cache);

                if(has_background) // FIXME: reversed condition?
                    ++resVal;
                nSeg[i] = resVal;
            }
        );

        for(int i=0; i<d;++i){
            offsets_it[i] = current_offset;
            current_offset+=nSeg[i];
        }


        unmerged_label_number = current_offset;
        if(!has_background)
            ++unmerged_label_number;
    }

    // reduce stage: merge adjacent labels if the region overlaps
    UnionFindArray<Label> global_unions(unmerged_label_number);
    if(has_background)
    {
        // merge all labels that refer to background
        for(typename MultiArray<M, Label>::iterator offsets_it = label_offsets.begin();
                offsets_it != label_offsets.end();
                ++offsets_it)
        {
            global_unions.makeUnion(0, *offsets_it);
        }
    }


    typedef GridGraph<M, undirected_tag> Graph;
    typedef typename Graph::edge_iterator EdgeIterator;
    Graph blocks_graph(blocks_shape, options.getNeighborhood());
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

        auto block_u = blocking.blockDescToBlock(u);
        auto block_v = blocking.blockDescToBlock(v);

        MultiArray<M,Data> data_cache_u(block_u.size());
        MultiArray<M,Data> data_cache_v(block_v.size());

        MultiArray<M,Label> label_cache_u(block_u.size());
        MultiArray<M,Label> label_cache_v(block_v.size());

        data.checkoutSubarray(block_u.begin(), data_cache_u);
        data.checkoutSubarray(block_v.begin(), data_cache_v);

        labels.checkoutSubarray(block_u.begin(), label_cache_u);
        labels.checkoutSubarray(block_v.begin(), label_cache_v);

        visitBorder(data_cache_u, label_cache_u,
                    data_cache_v, label_cache_v,
                    difference, options.getNeighborhood(), border_visitor);

        labels.commitSubarray(block_u.begin(), label_cache_u);
        labels.commitSubarray(block_v.begin(), label_cache_v);
    }

    // fill mapping (local labels) -> (global labels)
    Label last_label = global_unions.makeContiguous();
    {
        typename MultiArray<M, Label>::iterator offsets_it = label_offsets.begin();
        Label offset = *offsets_it;
        ++offsets_it;
        typename Mapping::iterator mapping_it = mapping.begin();
        for( ; offsets_it != label_offsets.end(); ++offsets_it, ++mapping_it)
        {
            mapping_it->clear();
            Label next_offset = *offsets_it;
            if(has_background)
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
        if(has_background)
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


template <unsigned int M, unsigned int N, class Label, class TAG, class MappingIterator>
void chunkedToGlobalLabels(ChunkedArray<N,Label,TAG> & labels,
                           const MultiBlocking<M, MultiArrayIndex> & blocking,
                           MappingIterator mapping_begin,
                           MappingIterator mapping_end)
{
    for(auto blockIt = blocking.blockBegin(); blockIt != blocking.blockEnd(); ++blockIt, ++mapping_begin)
    {
        vigra_assert(mapping_begin != mapping_end, "");

        auto block = *blockIt;
        MultiArray<M,Label> label_cache(block.size());
        labels.checkoutSubarray(block.begin(), label_cache);

        for(auto labels_it = label_cache.begin(); labels_it != label_cache.end(); ++labels_it)
        {
            vigra_assert(*labels_it < mapping_begin->size(), "");
            *labels_it = (*mapping_begin)[*labels_it];
        }

        labels.commitSubarray(block.begin(), label_cache);
    }
}


} // namespace blockwise_labeling_detail

/*************************************************************/
/*                                                           */
/*                      labelMultiArrayBlockwise             */
/*                                                           */
/*************************************************************/

/** \weakgroup ParallelProcessing
    \sa labelMultiArrayBlockwise <B>(...)</B>
*/

/** \brief Connected components labeling for MultiArrays and ChunkedArrays.

    <b> Declarations:</b>

    \code
    namespace vigra {
        // assign local labels and generate mapping (local labels) -> (global labels) for each chunk
        template <unsigned int N, class Data, class S1,
                                  class Label, class S2,
                  class Equal, class S3>
        Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                                       MultiArrayView<N, Label, S2> labels,
                                       const BlockwiseLabelOptions& options,
                                       Equal equal,
                                       MultiArrayView<N, std::vector<Label>, S3>& mapping);

        // assign local labels and generate mapping (local labels) -> (global labels) for each chunk
        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class EqualityFunctor>
        Label labelMultiArrayBlockwise(const ChunkedArray<N, Data>& data,
                                       ChunkedArray<N, Label>& labels,
                                       const BlockwiseLabelOptions& options,
                                       Equal equal,
                                       MultiArrayView<N, std::vector<Label>, S3>& mapping);

        // assign global labels
        template <unsigned int N, class Data, class S1,
                                  class Label, class S2,
                  class Equal>
        Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                                       MultiArrayView<N, Label, S2> labels,
                                       const BlockwiseLabelOptions& options,
                                       Equal equal);

        // assign global labels
        template <unsigned int N, class T, class S1,
                                  class Label, class S2,
                  class EqualityFunctor = std::equal_to<T> >
        Label labelMultiArrayBlockwise(const ChunkedArray<N, Data>& data,
                                       ChunkedArray<N, Label>& labels,
                                       const BlockwiseLabelOptions& options = BlockwiseLabelOptions(),
                                       Equal equal = std::equal_to<T>());
    }
    \endcode

    The resulting labeling is equivalent to a labeling by \ref labelMultiArray, that is, the connected components are the same but may have different ids.
    \ref NeighborhoodType and background value (if any) can be specified with the LabelOptions object.
    If the \a mapping parameter is provided, each chunk is labeled seperately and contiguously (starting at one, zero for background),
    with \a mapping containing a mapping of local labels to global labels for each chunk.
    Thus, the shape of 'mapping' has to be large enough to hold each chunk coordinate.

    Return: the number of regions found (=largest global region label)

    <b> Usage: </b>

    <b>\#include </b> \<vigra/blockwise_labeling.hxx\><br>
    Namespace: vigra

    \code
    Shape3 shape = Shape3(10);
    Shape3 chunk_shape = Shape3(4);
    ChunkedArrayLazy<3, int> data(shape, chunk_shape);
    // fill data ...

    ChunkedArrayLazy<3, size_t> labels(shape, chunk_shape);

    MultiArray<3, std::vector<size_t> > mapping(Shape3(3)); // there are 3 chunks in each dimension

    labelMultiArrayBlockwise(data, labels, LabelOptions().neighborhood(DirectNeighborhood).background(0),
                             std::equal_to<int>(), mapping);

    // check out chunk in the middle
    MultiArray<3, size_t> middle_chunk(Shape3(4));
    labels.checkoutSubarray(Shape3(4), middle_chunk);

    // print number of non-background labels assigned in middle_chunk
    cout << mapping[Shape3(1)].size() << endl;

    // get local label for voxel
    // this may be the same value assigned to different component in another chunk
    size_t local_label = middle_chunk[Shape3(2)];
    // get global label for voxel
    // if another voxel has the same label, they are in the same connected component albeit they may be in different chunks
    size_t global_label = mapping[Shape3(1)][local_label
    \endcode
    */
doxygen_overloaded_function(template <...> unsigned int labelMultiArrayBlockwise)

template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Equal, class S3>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               const BlockwiseLabelOptions<N>& options,
                               Equal equal,
                               MultiArrayView<N, std::vector<Label>, S3>& mapping)
{
    using namespace blockwise_labeling_detail;

    typedef typename MultiArrayShape<N>::type Shape;
    Shape block_shape(options.getBlockShape());

    MultiArray<N, MultiArrayView<N, Data, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);
    return blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                             label_blocks.begin(), label_blocks.end(),
                             options, equal, mapping);
}

template <unsigned int N, class Data, class S1,
                          class Label, class S2,
          class Equal>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               const BlockwiseLabelOptions<N>& options,
                               Equal equal)
{
    using namespace blockwise_labeling_detail;

    typedef typename MultiArrayShape<N>::type Shape;
    Shape block_shape(options.getBlockShape());

    MultiArray<N, MultiArrayView<N, Data, S1> > data_blocks = blockify(data, block_shape);
    MultiArray<N, MultiArrayView<N, Label, S2> > label_blocks = blockify(labels, block_shape);
    MultiArray<N, std::vector<Label> > mapping(data_blocks.shape());
    Label last_label = blockwiseLabeling(data_blocks.begin(), data_blocks.end(),
                                         label_blocks.begin(), label_blocks.end(),
                                         options, equal, mapping);

    // replace local labels by global labels
    toGlobalLabels(label_blocks.begin(), label_blocks.end(), mapping.begin(), mapping.end());
    return last_label;
}

template <unsigned int N, class Data, class S1,
                          class Label, class S2>
Label labelMultiArrayBlockwise(const MultiArrayView<N, Data, S1>& data,
                               MultiArrayView<N, Label, S2> labels,
                               const BlockwiseLabelOptions<N>& options = BlockwiseLabelOptions<N>())
{
    return labelMultiArrayBlockwise(data, labels, options, std::equal_to<Data>());
}



template <unsigned int M, unsigned int N1, unsigned int N2, class Data, class Label, class Equal, class S3, class TAG1, class TAG2>
typename std::enable_if<(N1 > M && N2 > M), Label>::type labelMultiArrayBlockwise(
        const ChunkedArray<N1, Data, TAG1>& data,
        ChunkedArray<N2, Label, TAG2>& labels,
        const BlockwiseLabelOptions<M>& options,
        Equal equal,
        MultiArrayView<M, std::vector<Label>, S3> mapping)
{
    using namespace blockwise_labeling_detail;

    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(data.template minimalShape<M>() == labels.template minimalShape<M>(),
            "Zhe data shape must be equal to the destination shape.");
    vigra_precondition(allGreaterEqual(data.template minimalShape<M>(), options.getBlockShape()),
            "The data shape must be equal or greater then the block shape.");

    const Blocking blocking(data.template minimalShape<M>(), options.getBlockShape());

    return chunkedBlockwiseLabeling(data, labels, blocking, options, equal, mapping);
}


template <unsigned int M, unsigned int N1, unsigned int N2, class Data, class Label, class Equal, class TAG1, class TAG2>
typename std::enable_if<(N1 >= M && N2 >= M), Label>::type labelMultiArrayBlockwise(
        const ChunkedArray<N1, Data, TAG1>& data,
        ChunkedArray<N2, Label, TAG2>& labels,
        const BlockwiseLabelOptions<M>& options,
        Equal equal)
{
    using namespace blockwise_labeling_detail;

    typedef TinyVector<MultiArrayIndex, M> Shape;
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(data.template minimalShape<M>() == labels.template minimalShape<M>(),
            "The data shape must be equal to the destination shape.");
    vigra_precondition(allGreaterEqual(data.template minimalShape<M>(), options.getBlockShape()),
            "The data shape must be equal or greater then the block shape.");

    const Blocking blocking(data.template minimalShape<M>(), options.getBlockShape());

    Shape blockArrayShape = data.template minimalShape<M>() / options.getBlockShape();
    Shape reminder = data.template minimalShape<M>() % options.getBlockShape();
    for (unsigned i = 0; i < M; ++i)
        if (reminder[i] != 0)
            blockArrayShape[i] += 1;

    MultiArray<M, std::vector<Label> > mapping(blockArrayShape);
    Label result = chunkedBlockwiseLabeling(data, labels, blocking, options, equal, mapping);

    chunkedToGlobalLabels(labels, blocking, mapping.begin(), mapping.end());
    return result;
}


template <unsigned int M, unsigned int N1, unsigned int N2, class Data, class Label, class TAG1, class TAG2>
typename std::enable_if<(N1 >= M && N2 >= M), Label>::type inline labelMultiArrayBlockwise(
        const ChunkedArray<N1, Data, TAG1>& data,
        ChunkedArray<N2, Label, TAG2>& labels,
        const BlockwiseLabelOptions<M>& options)
{
    return labelMultiArrayBlockwise(data, labels, options, std::equal_to<Data>());
}

//@}

} // namespace vigra

#endif // VIGRA_BLOCKWISE_LABELING_HXX
