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

#ifndef VIGRA_BLOCKIFY_HXX
#define VIGRA_BLOCKIFY_HXX

#include "multi_array.hxx"

namespace vigra
{

namespace blockify_detail
{

template <unsigned int CurrentDimensions>
struct blockify_impl
{
    // for CurrentDimension >= 1
    template <unsigned int N, class T, class S, class Shape>
    static void make(MultiArrayView<N, T, S>& source, 
                     MultiArrayView<N, MultiArrayView<N, T, S> >& blocks,
                     Shape current_block_begin,
                     Shape current_block_end,
                     Shape current_block_pos,
                     Shape block_shape)
    {
        typedef typename Shape::value_type size_type;
        enum{ n = CurrentDimensions - 1};

        size_type blocks_extend = blocks.shape(n);
        
        vigra_assert(blocks_extend != 0, "");
        for(current_block_pos[n] = 0, current_block_begin[n] = 0, current_block_end[n] = block_shape[n];
            current_block_pos[n] != blocks_extend - 1;
            ++current_block_pos[n],
                current_block_begin[n] += block_shape[n], 
                current_block_end[n] += block_shape[n])
        {
            blockify_impl<n>::make(source, blocks, current_block_begin, current_block_end, current_block_pos, block_shape);
        }
        current_block_end[n] = source.shape(n);
        blockify_impl<n>::make(source, blocks, current_block_begin, current_block_end, current_block_pos, block_shape);
    }
};

template <>
struct blockify_impl<0>
{
    template <unsigned int N, class T, class S, class Shape>
    static void make(MultiArrayView<N, T, S>& source, 
                     MultiArrayView<N, MultiArrayView<N, T, S> >& blocks,
                     Shape current_block_begin,
                     Shape current_block_end,
                     Shape current_block_pos,
                     Shape block_shape)
    {
        blocks[current_block_pos] = source.subarray(current_block_begin, current_block_end);
    }
};

} // namespace blockify_detail

template <unsigned int N, class T, class S>
MultiArray<N, MultiArrayView<N, T, S> >
blockify(MultiArrayView<N, T, S> source, typename MultiArrayView<N, T, S>::difference_type block_shape)
{
    using namespace blockify_detail;
    typedef typename MultiArrayView<N, T, S>::difference_type Shape;

    Shape blocks_shape;
    for(unsigned int n = 0; n != N; ++n)
    {
        blocks_shape[n] = source.shape(n) / block_shape[n];
        if(blocks_shape[n] * block_shape[n] != source.shape(n))
            ++blocks_shape[n];
    }
    MultiArray<N, MultiArrayView<N, T, S> > blocks(blocks_shape);
    if(source.size() == 0)
        return blocks;
    Shape a;
    Shape b;
    Shape c;
    blockify_impl<N>::make(source, blocks, a, b, c, block_shape);
    return blocks;
}

} // namespace vigra

#endif // VIGRA_BLOCKIFY_HXX
