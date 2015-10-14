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

#ifndef VIGRA_BLOCKWISE_WATERSHEDS_HXX
#define VIGRA_BLOCKWISE_WATERSHEDS_HXX

#include "multi_array.hxx"
#include "multi_gridgraph.hxx"
#include "blockify.hxx"
#include "blockwise_labeling.hxx"
#include "overlapped_blocks.hxx"

#include <limits>

namespace vigra
{

/** \addtogroup SeededRegionGrowing
*/
//@{

namespace blockwise_watersheds_detail
{

template <class DataArray, class DirectionsBlocksIterator>
void prepareBlockwiseWatersheds(const Overlaps<DataArray>& overlaps,
                                DirectionsBlocksIterator directions_blocks_begin,
                                NeighborhoodType neighborhood)
{
    static const unsigned int N = DataArray::actual_dimension;
    typedef typename MultiArrayShape<N>::type Shape;
    typedef typename DirectionsBlocksIterator::value_type DirectionsBlock;
    Shape shape = overlaps.shape();
    vigra_assert(shape == directions_blocks_begin.shape(), "");
    
    MultiCoordinateIterator<N> it3(shape);
    MultiCoordinateIterator<N> end = it3.getEndIterator();

    // trivial prarallel. )(halo 1)
    
    #pragma omp parallel for
    for(MultiCoordinateIterator<N> it=it3; it < end; ++it)
    {
        DirectionsBlock directions_block = directions_blocks_begin[*it];
        OverlappingBlock<DataArray> data_block = overlaps[*it];
        
        typedef GridGraph<N, undirected_tag> Graph;
        typedef typename Graph::NodeIt GraphScanner;
        typedef typename Graph::OutArcIt NeighborIterator;
        
        Graph graph(data_block.block.shape(), neighborhood);
        for(GraphScanner node(graph); node != lemon::INVALID; ++node)
        {
            if(within(*node, data_block.inner_bounds))
            {
                typedef typename DataArray::value_type Data;
                Data lowest_neighbor = data_block.block[*node];
                
                typedef typename DirectionsBlock::value_type Direction;
                Direction lowest_neighbor_direction = std::numeric_limits<unsigned short>::max();
                
                for(NeighborIterator arc(graph, *node); arc != lemon::INVALID; ++arc)
                {
                    Shape neighbor_coordinates = graph.target(*arc);
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
struct UnionFindWatershedsEquality
{
    // FIXME: this graph object shouldn't be necessary, most functions (and state) of graph are not used
    // this probably needs some refactoring in GridGraph
    GridGraph<N, undirected_tag>* graph;

    template <class Shape>
    bool operator()(unsigned short u, const unsigned short v, const Shape& diff) const
    {
        static const unsigned short plateau_id = std::numeric_limits<unsigned short>::max();
        return (u == plateau_id && v == plateau_id) ||
               (u != plateau_id && graph->neighborOffset(u) == diff) ||
               (v != plateau_id && graph->neighborOffset(graph->oppositeIndex(v)) == diff);
    }

    struct WithDiffTag
    {};
};

} // namespace blockwise_watersheds_detail

template <unsigned int N, class Data, class S1,
                          class Label, class S2>
Label unionFindWatershedsBlockwise(MultiArrayView<N, Data, S1> data,
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

    Overlaps<MultiArrayView<N, Data, S1> > overlaps(data, block_shape, Shape(1), Shape(1));
    prepareBlockwiseWatersheds(overlaps, directions_blocks.begin(), neighborhood);
    GridGraph<N, undirected_tag> graph(data.shape(), neighborhood);
    UnionFindWatershedsEquality<N> equal = {&graph};
    return labelMultiArrayBlockwise(directions, labels, LabelOptions().neighborhood(neighborhood).blockShape(block_shape), equal);
}

/*************************************************************/
/*                                                           */
/*                      unionFindWatershedsBlockwise         */
/*                                                           */
/*************************************************************/

/** \brief Blockise union-find watersheds transform for ChunkedArrays.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template <unsigned int N, class Data, class Label>
        Label unionFindWatershedsBlockwise(const ChunkedArray<N, Data>& data,
                                          ChunkedArray<N, Label>& labels,
                                          NeighborhoodType neighborhood = DirectNeighborhood);

        // provide temporary directions storage
        template <unsigned int N, class Data, class Label>
        Label unionFindWatershedsBlockwise(const ChunkedArray<N, Data>& data,
                                          ChunkedArray<N, Label>& labels,
                                          NeighborhoodType neighborhood,
                                          ChunkedArray<N, unsigned short>& temporary_storage);
    }
    \endcode
    
    The resulting labeling is equivalent to a labeling by \ref watershedsUnionFind, that is,
    the components are the same but may have different ids.
    If \a temporary_storage is provided, this array is used for intermediate result storage.
    Otherwise, a newly created \ref vigra::ChunkedArrayLazy is used.

    Return: the number of labels assigned (=largest label, because labels start at one)
    
    <b> Usage: </b>

    <b>\#include </b> \<vigra/blockwise_watersheds.hxx\><br>
    Namespace: vigra

    \code
    Shape3 shape = Shape3(10);
    Shape3 chunk_shape = Shape3(4);
    ChunkedArrayLazy<3, int> data(shape, chunk_shape);
    // fill data ...
    
    ChunkedArrayLazy<3, size_t> labels(shape, chunk_shape);
    
    unionFindWatershedsBlockwise(data, labels, IndirectNeighborhood);
    \endcode
    */
doxygen_overloaded_function(template <...> unsigned int unionFindWatershedsBlockwise)

template <unsigned int N, class Data, class Label>
Label unionFindWatershedsBlockwise(const ChunkedArray<N, Data>& data,
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
    
    Overlaps<ChunkedArray<N, Data> > overlaps(data, data.chunkShape(), Shape(1), Shape(1));
    
    prepareBlockwiseWatersheds(overlaps, directions.chunk_begin(Shape(0), shape), neighborhood);
    
    GridGraph<N, undirected_tag> graph(shape, neighborhood);
    UnionFindWatershedsEquality<N> equal = {&graph};
    return labelMultiArrayBlockwise(directions, labels, LabelOptions().neighborhood(neighborhood), equal);
}

template <unsigned int N, class Data,
                          class Label>
inline Label 
unionFindWatershedsBlockwise(const ChunkedArray<N, Data>& data,
                                   ChunkedArray<N, Label>& labels,
                                   NeighborhoodType neighborhood = DirectNeighborhood)
{
    ChunkedArrayLazy<N, unsigned short> directions(data.shape(), data.chunkShape());
    return unionFindWatershedsBlockwise(data, labels, neighborhood, directions);
}

//@}

} // namespace vigra

#endif // VIGRA_BLOCKWISE_WATERSHEDS_HXX
