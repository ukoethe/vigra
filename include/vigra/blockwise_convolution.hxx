#ifndef VIGRA_BLOCKWISE_CONVOLUTION_HXX_
#define VIGRA_BLOCKWISE_CONVOLUTION_HXX_

#include <vigra/overlapped_blocks.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/blockify.hxx>
#include <vigra/multi_array.hxx>

namespace vigra
{

namespace blockwise_convolution_detail
{

template <class DataArray, class OutputBlocksIterator, class KernelIterator>
void convolveImpl(const Overlaps<DataArray>& overlaps, OutputBlocksIterator output_blocks_begin, KernelIterator kit)
{
    static const unsigned int N = DataArray::actual_dimension;
    typedef typename MultiArrayShape<N>::type Shape;
    typedef typename OutputBlocksIterator::value_type OutputBlock;

    Shape shape = overlaps.shape();
    vigra_assert(shape == output_blocks_begin.shape(), "");

    MultiCoordinateIterator<N> it(shape);
    MultiCoordinateIterator<N> end = it.getEndIterator();
    for( ; it != end; ++it)
    {
        OutputBlock output_block = output_blocks_begin[*it];
        OverlappingBlock<DataArray> data_block = overlaps[*it];
        separableConvolveMultiArray(data_block.block, output_block, kit, data_block.inner_bounds.first, data_block.inner_bounds.second);
    }
}

template <class Shape, class KernelIterator>
std::pair<Shape, Shape> kernelOverlap(KernelIterator kit)
{
    Shape before;
    Shape after;

    for(unsigned int i = 0; i != Shape::static_size; ++i, ++kit)
    {
        // FIXME: is this correct?
        before[i] = kit->right();
        after[i] = -kit->left();
    }
    return std::make_pair(before, after);
}

}


template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class KernelIterator>
void separableConvolveBlockwise(MultiArrayView<N, T1, S1> source, MultiArrayView<N, T2, S2> dest, KernelIterator kit,
                                const typename MultiArrayView<N, T1, S1>::difference_type& block_shape =
                                     typename MultiArrayView<N, T1, S1>::difference_type(128))
{
    using namespace blockwise_convolution_detail;

    typedef typename MultiArrayView<N, T1, S1>::difference_type Shape;
    
    Shape shape = source.shape();
    vigra_precondition(shape == dest.shape(), "shape mismatch of source and destination");
    
    std::pair<Shape, Shape> overlap = kernelOverlap<Shape, KernelIterator>(kit);
    Overlaps<MultiArrayView<N, T2, S2> > overlaps(source, block_shape, overlap.first, overlap.second);

    MultiArray<N, MultiArrayView<N, T2, S2> > destination_blocks = blockify(dest, block_shape);
    
    convolveImpl(overlaps, destination_blocks.begin(), kit);
}
template <unsigned int N, class T1, class S1,
                          class T2, class S2,
          class T3>
void separableConvolveBlockwise(MultiArrayView<N, T1, S1> source, MultiArrayView<N, T2, S2> dest, const Kernel1D<T3>& kernel,
                                const typename MultiArrayView<N, T1, S1>::difference_type& block_shape =
                                     typename MultiArrayView<N, T1, S1>::difference_type(128))
{
    std::vector<Kernel1D<T3> > kernels(N, kernel);
    separableConvolveBlockwise(source, dest, kernels.begin(), block_shape);
}
template <unsigned int N, class T1, class T2, class KernelIterator>
void separableConvolveBlockwise(const ChunkedArray<N, T1>& source, ChunkedArray<N, T2>& destination, KernelIterator kit)
{
    using namespace blockwise_convolution_detail;

    typedef typename ChunkedArray<N, T1>::shape_type Shape;
    
    Shape shape = source.shape();
    vigra_precondition(shape == destination.shape(), "shape mismatch of source and destination");

    std::pair<Shape, Shape> overlap = kernelOverlap<Shape, KernelIterator>(kit);
    Shape block_shape = source.chunkShape();
    vigra_precondition(block_shape == destination.chunkShape(), "chunk shapes do not match");
    Overlaps<ChunkedArray<N, T1> > overlaps(source, block_shape, overlap.first, overlap.second);
    
    convolveImpl(overlaps, destination.chunk_begin(Shape(0), shape), kit);
}
template <unsigned int N, class T1, class T2, class T>
void separableConvolveBlockwise(const ChunkedArray<N, T1>& source, ChunkedArray<N, T2>& destination, const Kernel1D<T>& kernel)
{
    std::vector<Kernel1D<T> > kernels(N, kernel);
    separableConvolveBlockse(source, destination, kernels.begin());
}


}

#endif

