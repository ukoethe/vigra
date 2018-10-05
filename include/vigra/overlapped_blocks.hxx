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

#ifndef VIGRA_OVERLAPPED_BLOCKS_HXX
#define VIGRA_OVERLAPPED_BLOCKS_HXX

#include <utility>
#include <algorithm>

#include <vigra/multi_array.hxx>
#include <vigra/multi_array_chunked.hxx>

namespace vigra
{

namespace overlapped_blocks_detail
{

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
std::pair<Shape, Shape> overlapBoundsAt(const std::pair<Shape, Shape>& block_bounds, const Shape& global_shape,
                                        const Shape& overlap_before, const Shape& overlap_after)
{
    Shape overlapped_block_begin = block_bounds.first;
    Shape overlapped_block_end = block_bounds.second;
    for(int i = 0; i != Shape::static_size; ++i)
    {
        if(overlapped_block_begin[i] >= overlap_before[i])
            overlapped_block_begin[i] -= overlap_before[i];
        else
            overlapped_block_begin[i] = 0;

        if(overlapped_block_end[i] <= global_shape[i] - overlap_after[i])
            overlapped_block_end[i] += overlap_after[i];
        else
            overlapped_block_end[i] = global_shape[i];
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

} // namespace overlapped_blocks_detail

template <class Shape>
inline bool
within(const Shape& coordinates, const std::pair<Shape, Shape>& bounds)
{
    return allLessEqual(bounds.first, coordinates) && allLess(coordinates, bounds.second);
}

template <class ArrayType>
struct OverlappingBlock;

template <class ArrayType>
class Overlaps;

template <unsigned int N, class T, class S>
struct OverlappingBlock<MultiArrayView<N, T, S> >
{
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;

    MultiArrayView<N, T, S> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int N, class T, class S>
class Overlaps<MultiArrayView<N, T, S> >
{
private:
    typedef MultiArrayView<N, T, S> View;
    typedef typename View::difference_type Shape;

    View view;
    Shape block_shape;
    Shape overlap_before;
    Shape overlap_after;
public:
    Overlaps(View view, const Shape& block_shape, const Shape& overlap_before, const Shape& overlap_after)
    : view(view),
      block_shape(block_shape),
      overlap_before(overlap_before),
      overlap_after(overlap_after)
    {}
    OverlappingBlock<View> operator[](const Shape& coordinates) const
    {
        using namespace overlapped_blocks_detail;
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, view.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, view.shape(), overlap_before, overlap_after);

        OverlappingBlock<View> result;
        result.block = view.subarray(overlap_bounds.first, overlap_bounds.second);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);
        return result;
    }
    Shape shape() const
    {
        using namespace overlapped_blocks_detail;
        return blocksShape(view.shape(), block_shape);
    }
};

template <unsigned int N, class T, class TAG>
struct OverlappingBlock<ChunkedArray<N, T, TAG> >
{
    typedef typename MultiArrayShape<N>::type Shape;

    MultiArray<N, T> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int N, class T, class TAG>
class Overlaps<ChunkedArray<N, T, TAG> >
{
private:
    typedef ChunkedArray<N, T, TAG> Array;
    typedef typename MultiArrayShape<N>::type Shape;

    const Array& array;
    Shape block_shape;
    Shape overlap_before;
    Shape overlap_after;
public:
    Overlaps(const Array& array, const Shape& block_shape, const Shape& overlap_before, const Shape& overlap_after)
    : array(array),
      block_shape(block_shape),
      overlap_before(overlap_before),
      overlap_after(overlap_after)
    {}

    OverlappingBlock<Array> operator[](const Shape& coordinates) const
    {
        using namespace overlapped_blocks_detail;
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, array.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, array.shape(), overlap_before, overlap_after);

        OverlappingBlock<Array> result;
        result.block.reshape(overlap_bounds.second - overlap_bounds.first);
        array.checkoutSubarray(overlap_bounds.first, result.block);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);

        return result;
    }
    Shape shape() const
    {
        using namespace overlapped_blocks_detail;
        return blocksShape(array.shape(), block_shape);
    }
};



template <unsigned int M, class ArrayType>
struct OverlappingBlockN;

template <unsigned int M, class ArrayType>
class OverlapsN;

template <unsigned int N, class T, class S>
struct OverlappingBlockN<N, MultiArrayView<N, T, S> >
{
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;

    MultiArrayView<N, T, S> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int N, class T, class TAG>
struct OverlappingBlockN<N, ChunkedArray<N, T, TAG> >
{
    typedef typename MultiArrayShape<N>::type Shape;

    MultiArray<N, T> block;
    std::pair<Shape, Shape> inner_bounds;
};


template <unsigned int N, class T, class S>
class OverlapsN<N, MultiArrayView<N, T, S> >
{
private:
    typedef MultiArrayView<N, T, S> View;
    typedef typename View::difference_type Shape;

    View view;
    Shape block_shape;
    Shape overlap_before;
    Shape overlap_after;
public:
    OverlapsN(View view, const Shape& block_shape, const Shape& overlap_before, const Shape& overlap_after)
    : view(view),
      block_shape(block_shape),
      overlap_before(overlap_before),
      overlap_after(overlap_after)
    {}
    OverlappingBlockN<N, View> operator[](const Shape& coordinates) const
    {
        using namespace overlapped_blocks_detail;
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, view.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, view.shape(), overlap_before, overlap_after);

        OverlappingBlockN<N, View> result;
        result.block = view.subarray(overlap_bounds.first, overlap_bounds.second);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);
        return result;
    }
    Shape shape() const
    {
        using namespace overlapped_blocks_detail;
        return blocksShape(view.shape(), block_shape);
    }
};

template <unsigned int N, class T, class TAG>
class OverlapsN<N, ChunkedArray<N, T, TAG> >
{
private:
    typedef ChunkedArray<N, T, TAG> Array;
    typedef typename MultiArrayShape<N>::type Shape;

    const Array& array;
    Shape block_shape;
    Shape overlap_before;
    Shape overlap_after;
public:
    OverlapsN(const Array& array, const Shape& block_shape, const Shape& overlap_before, const Shape& overlap_after)
    : array(array),
      block_shape(block_shape),
      overlap_before(overlap_before),
      overlap_after(overlap_after)
    {}

    OverlappingBlockN<N, Array> operator[](const Shape& coordinates) const
    {
        using namespace overlapped_blocks_detail;
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, array.shape(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, array.shape(), overlap_before, overlap_after);

        OverlappingBlockN<N, Array> result;
        result.block.reshape(overlap_bounds.second - overlap_bounds.first);
        array.checkoutSubarray(overlap_bounds.first, result.block);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);

        return result;
    }
    Shape shape() const
    {
        using namespace overlapped_blocks_detail;
        return blocksShape(array.shape(), block_shape);
    }
};


template <unsigned int M, unsigned int N, class T, class TAG>
struct OverlappingBlockN<M, ChunkedArray<N, T, TAG> >
{
    typedef typename MultiArrayShape<M>::type Shape;

    MultiArray<M, T> block;
    std::pair<Shape, Shape> inner_bounds;
};

template <unsigned int M, unsigned int N, class T, class TAG>
class OverlapsN<M, ChunkedArray<N, T, TAG> >
{
private:
    typedef ChunkedArray<N, T, TAG> Array;
    typedef typename MultiArrayShape<M>::type Shape;

    const Array& array;
    Shape block_shape;
    Shape overlap_before;
    Shape overlap_after;
public:
    OverlapsN(const Array& array, const Shape& block_shape, const Shape& overlap_before, const Shape& overlap_after)
    : array(array),
      block_shape(block_shape),
      overlap_before(overlap_before),
      overlap_after(overlap_after)
    {}

    OverlappingBlockN<M, Array> operator[](const Shape& coordinates) const
    {
        using namespace overlapped_blocks_detail;
        std::pair<Shape, Shape> block_bounds = blockBoundsAt(coordinates, array.template minimalShape<M>(), block_shape);
        std::pair<Shape, Shape> overlap_bounds = overlapBoundsAt(block_bounds, array.template minimalShape<M>(), overlap_before, overlap_after);

        OverlappingBlockN<M, Array> result;
        result.block.reshape(overlap_bounds.second - overlap_bounds.first);
        array.checkoutSubarray(overlap_bounds.first, result.block);
        result.inner_bounds = std::make_pair(block_bounds.first - overlap_bounds.first, block_bounds.second - overlap_bounds.first);

        return result;
    }
    Shape shape() const
    {
        using namespace overlapped_blocks_detail;
        return blocksShape(array.template minimalShape<M>(), block_shape);
    }
};


} // namespace vigra

#endif // VIGRA_OVERLAPPED_BLOCKS_HXX
