#ifndef VIGRA_BLOCKWISE_WATERSHEDS_HXX_
#define VIGRA_BLOCKWISE_WATERSHEDS_HXX_

#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/blockify.hxx>
#include <vigra/blockwise_labeling.hxx>

#include <limits>

namespace vigra
{

namespace blockwise_watersheds_detail
{

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

template <class Shape>
std::pair<Shape, Shape> blockBoundsAt(const Shape& coordinates, const Shape& global_shape, const Shape& block_shape)
{
    Shape block_begin;
    for(int i = 0; i != Shape::static_size; ++i)
    {
        block_begin[i] = coordinates[i] * block_shape[i];
        vigra_assert(block_begin[i] < global_shape[i], "block coordinates out of bounds");
    }
    Shape block_end;
    for(int i = 0; i != Shape::static_size; ++i)
    {
        block_end[i] = std::min(block_begin[i] + block_shape[i], global_shape[i]);
    }
    return std::make_pair(block_begin, block_end);
}
template <class Shape>
std::pair<Shape, Shape> overlapBoundsAt(const std::pair<Shape, Shape>& block_bounds, const Shape& global_shape)
{
    Shape overlapped_block_begin = block_bounds.first;
    Shape overlapped_block_end = block_bounds.second;
    for(int i = 0; i != Shape::static_size; ++i)
    {
        if(overlapped_block_begin[i] != 0)
            --overlapped_block_begin[i];
        if(overlapped_block_end[i] != global_shape[i])
            ++overlapped_block_end[i];
    }
    return std::make_pair(overlapped_block_begin, overlapped_block_end);
}
template <class Shape>
Shape blocksShape(const Shape& global_shape, const Shape& block_shape)
{
    Shape result;
    for(int i = 0; i != Shape::static_size; ++i)
    {
        result[i] = global_shape[i] / block_shape[i];
        if(block_shape[i] * result[i] != global_shape[i])
            ++result[i];
    }
    return result;

}


template <unsigned int N, class T, class S>
struct OverlappingBlock
{
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;
    
    MultiArrayView<N, T, S> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int N, class T, class S>
class OverlapsGenerator
{
private:
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;

    MultiArrayView<N, T, S> arr;
    Shape block_shape;
public:
    OverlapsGenerator(MultiArrayView<N, T, S> arr, const Shape& block_shape)
      : arr(arr),
        block_shape(block_shape)
    {}
    OverlappingBlock<N, T, S> operator[](const Shape& coordinates) const
    {
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, arr.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, arr.shape());
        
        OverlappingBlock<N, T, S> result;
        result.block = arr.subarray(overlap_bounds.first, overlap_bounds.second);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);
        return result;
    }
    Shape shape() const
    {
        return blocksShape(arr.shape(), block_shape);
    }
};

template <unsigned int N, class T>
struct OverlappingBlock<N, T, ChunkedArrayTag>
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    MultiArray<N, T> block; // not a view, overlaps are checked out
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int N, class T>
class OverlapsGenerator<N, T, ChunkedArrayTag>
{
private:
    typedef typename MultiArrayView<N, T, ChunkedArrayTag>::difference_type Shape;
    
    const ChunkedArray<N, T>& arr;
    Shape block_shape;
public:
    OverlapsGenerator(const ChunkedArray<N, T>& array, const Shape& block_shape)
      : arr(array),
        block_shape(block_shape)
    {}
    OverlappingBlock<N, T, ChunkedArrayTag> operator[](const Shape& coordinates) const
    {
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, arr.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, arr.shape());
        
        OverlappingBlock<N, T, ChunkedArrayTag> result;
        result.block.reshape(overlap_bounds.second - overlap_bounds.first);
        arr.checkoutSubarray(overlap_bounds.first, result.block);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);
        return result;
    }
    Shape shape() const
    {
        return blocksShape(arr.shape(), block_shape);
    }
};

template <unsigned int N, class T, class S,
          class DirectionsBlocksIterator>
void prepareBlockwiseWatersheds(OverlapsGenerator<N, T, S> overlaps,
                                DirectionsBlocksIterator directions_blocks_begin,
                                NeighborhoodType neighborhood)
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
        
        Graph graph(data_block.block.shape(), neighborhood);
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

template <unsigned int N>
struct UnionFindWatershedEquality
{
    GridGraph<N, undirected_tag>* graph;
    
    template <class Shape>
    bool operator()(unsigned short u, const unsigned short v, const Shape& diff) const
    {
        return (u == graph->maxDegree() && v == graph->maxDegree()) ||
               (u != graph->maxDegree() && graph->neighborOffset(u) == diff) ||
               (v != graph->maxDegree() && graph->neighborOffset(graph->oppositeIndex(v)) == diff);
    }

    struct WithDiffTag
    {};
};

}

template <unsigned int N, class Data, class S1,
                          class Label, class S2>
Label unionFindWatershedBlockwise(MultiArrayView<N, Data, S1> data,
                                  MultiArrayView<N, Label, S2> labels,
                                  NeighborhoodType neighborhood = DirectNeighborhood,
                                  const typename MultiArrayView<N, Data, S1>::difference_type& block_shape = 
                                          typename MultiArrayView<N, Data, S1>::difference_type(128))
{
    using namespace blockwise_watersheds_detail;

    typedef typename MultiArrayView<N, Data, S1>::difference_type Shape;
    Shape shape = data.shape();
    vigra_precondition(shape == labels.shape(), "shapes of data and labels do not match");
    
    MultiArray<N, unsigned short> directions(shape);
    
    MultiArray<N, MultiArrayView<N, unsigned short> > directions_blocks = blockify(directions, block_shape);

    OverlapsGenerator<N, Data, S1> overlaps(data, block_shape);
    prepareBlockwiseWatersheds(overlaps, directions_blocks.begin(), neighborhood);
    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    UnionFindWatershedEquality<N> equal = {&graph};
    return labelMultiArrayBlockwise(directions, labels, LabelOptions().neighborhood(neighborhood).blockShape(block_shape), equal);
}

template <unsigned int N, class Data,
                          class Label>
Label unionFindWatershedBlockwise(const ChunkedArray<N, Data>& data,
                                  ChunkedArray<N, Label>& labels,
                                  NeighborhoodType neighborhood,
                                  ChunkedArray<N, unsigned short>& directions)
{
    using namespace blockwise_watersheds_detail;
    
    typedef typename ChunkedArray<N, Data>::shape_type Shape;
    Shape shape = data.shape();
    vigra_precondition(shape == labels.shape() && shape == directions.shape(), "shapes of data and labels do not match");
    Shape chunk_shape = data.chunkShape();
    vigra_precondition(chunk_shape == labels.chunkShape() && chunk_shape == directions.chunkShape(), "chunk shapes do not match");
    
    OverlapsGenerator<N, Data, ChunkedArrayTag> overlaps(data, chunk_shape);
    
    prepareBlockwiseWatersheds(overlaps, directions.chunk_begin(Shape(0), shape), neighborhood);
    
    GridGraph<N, undirected_tag> graph(shape, neighborhood);
    UnionFindWatershedEquality<N> equal = {&graph};
    return labelMultiArrayBlockwise(directions, labels, LabelOptions().neighborhood(neighborhood), equal);
}
template <unsigned int N, class Data,
                          class Label>
Label unionFindWatershedBlockwise(const ChunkedArray<N, Data>& data,
                                  ChunkedArray<N, Label>& labels,
                                  NeighborhoodType neighborhood)
{
    ChunkedArrayLazy<N, unsigned short> directions(data.shape(), data.chunk_shape());
    return unionFindWatershedBlockwise(data, labels, neighborhood, directions);
}

}

#endif

