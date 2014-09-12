#ifndef VIGRA_BLOCKWISE_WATERSHEDS_HXX_
#define VIGRA_BLOCKWISE_WATERSHEDS_HXX_

#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/blockify.hxx>

#include <limits>

namespace vigra
{

namespace blockwise_watersheds_detail
{

template <unsigned int N, class T, class S>
struct OverlappingBlock
{
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;
    
    MultiArrayView<N, T, S> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <class Shape>
bool within(const Shape& coordinates, const std::pair<Shape, Shape>& bounds)
{
    for(unsigned int i = 0; i != Shape::static_size; ++i)
    {
        if(coordinates[i] < bounds.first[i] || coordinates[i] >= bounds.second[i])
            return false;
    }
    return true;
}

template <unsigned int N, class T, class S>
class OverlapsGenerator
{
private:
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;

    MultiArrayView<N, T, S> arr;
    Shape block_shape;
    
    std::pair<Shape, Shape> blockBoundsAt(const Shape& coordinates) const
    {
        Shape block_begin;
        for(int i = 0; i != N; ++i)
        {
            block_begin[i] = coordinates[i] * block_shape[i];
            vigra_assert(block_begin[i] < arr.shape(i), "block coordinates out of bounds");
        }
        Shape block_end;
        for(int i = 0; i != N; ++i)
        {
            block_end[i] = std::min(block_begin[i] + block_shape[i], arr.shape(i));
        }
        return std::make_pair(block_begin, block_end);
    }
    std::pair<Shape, Shape> overlapBoundsAt(const std::pair<Shape, Shape>& block_bounds) const
    {
        Shape overlapped_block_begin = block_bounds.first;
        Shape overlapped_block_end = block_bounds.second;
        for(int i = 0; i != N; ++i)
        {
            if(overlapped_block_begin[i] != 0)
                --overlapped_block_begin[i];
            if(overlapped_block_end[i] != arr.shape(i))
                ++overlapped_block_end[i];
        }
        return std::make_pair(overlapped_block_begin, overlapped_block_end);
    }
public:
    OverlapsGenerator(MultiArrayView<N, T, S> arr, const Shape& block_shape)
      : arr(arr),
        block_shape(block_shape)
    {}
    OverlappingBlock<N, T, S> operator[](const Shape& coordinates) const
    {
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds);
        
        OverlappingBlock<N, T, S> result;
        result.block = arr.subarray(overlap_bounds.first, overlap_bounds.second);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);
        return result;
    }
    Shape shape() const
    {
        Shape result;
        for(int i = 0; i != N; ++i)
        {
            result[i] = arr.shape(i) / block_shape[i];
            if(block_shape[i] * result[i] != arr.shape(i))
                ++result[i];
        }
        return result;
    }
};

template <class T>
struct FlaggedUnsigned
{
    typedef T DataType;
    static const T max_label = std::numeric_limits<T>::max() >> 1;
    static const T flag_bit = ~max_label;
    static const int flag_shift_bytes = std::numeric_limits<T>::digits - 1;
    
    T data;
    
    static T max()
    {
        return max_label;
    }
    static FlaggedUnsigned construct(T data, bool flag)
    {
        vigra_assert(flag <= 1, "flag out of bounds");
        vigra_assert(data <= max_label, "data out of bounds");
        return FlaggedUnsigned{data | (static_cast<T>(flag) << flag_shift_bytes)};
    }
    
    bool getFlag() const
    {
        return (data & flag_bit) != 0;
    }
    T getData() const
    {
        return data & max_label;
    }
};

template <unsigned int N, class T, class S,
          class DirectionsBlocksIterator>
void prepareBlockwiseWatersheds(OverlapsGenerator<N, T, S> overlaps,
                                DirectionsBlocksIterator directions_blocks_begin)
{
    typedef typename MultiArrayShape<N>::type Shape;
    typedef typename DirectionsBlocksIterator::value_type DirectionsBlock;
    Shape shape = overlaps.shape();
    vigra_assert(shape == directions_blocks_begin.shape(), "");
    
    MultiCoordinateIterator<N> it(shape);
    MultiCoordinateIterator<N> end = it.getEndIterator();
    for( ; it != end; ++it)
    {
        DirectionsBlock directions_block = directions_blocks_begin[*it];
        OverlappingBlock<N, T, S> data_block = overlaps[*it];
        
        typedef GridGraph<N, undirected_tag> Graph;
        typedef typename Graph::NodeIt GraphScanner;
        typedef typename Graph::OutArcIt NeighborIterator;
        
        Graph graph(data_block.block.shape());
        for(GraphScanner node(graph); node != lemon::INVALID; ++node)
        {
            if(within(*node, data_block.inner_bounds))
            {
                typedef typename MultiArrayView<N, T, S>::value_type Data;
                Data lowest_neighbor = data_block.block[*node];
                
                typedef typename DirectionsBlock::value_type Direction;
                Direction lowest_neighbor_direction = graph.maxDegree();
                
                for(NeighborIterator arc(graph, *node); arc != lemon::INVALID; ++arc)
                {
                    const Shape& neighbor_coordinates = graph.target(*arc);
                    Data neighbor_data = data_block.block[neighbor_coordinates];
                    if(neighbor_data < lowest_neighbor)
                    {
                        lowest_neighbor = neighbor_data;
                        lowest_neighbor_direction = arc.neighborIndex();
                    }
                }
                directions_block[*node - data_block.inner_bounds.first] = lowest_neighbor_direction;
            }
        }
    }
}

}

}

#endif

