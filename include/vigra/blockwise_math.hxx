/************************************************************************/
/*                                                                      */
/*     Copyright 2012-2014 by Ullrich Koethe and Kevin Kiefer           */
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


#ifndef VIGRA_MULTI_ARRAY_BLOCKWISE_MATH_HXX
#define VIGRA_MULTI_ARRAY_BLOCKWISE_MATH_HXX

#include "multi_array.hxx"
#include "multi_array_chunked.hxx"
#include "multi_blockwise.hxx"
#include "threading.hxx"
#include "matrix.hxx"
#include "linear_algebra.hxx"
#include "error.hxx"

namespace vigra {


template <unsigned int N, class T1, class T2, class T3, class S1, class S2, class S3>
void add(MultiArrayView<N,T1,S1> const & in1,
         MultiArrayView<N,T2,S2> const & in2,
         MultiArrayView<N,T3,S3> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == in2.shape(),
                        "In add(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.shape() == out.shape(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> in2_subarray = in2.subarray(block.begin(), block.end());
            MultiArrayView<N,T3,S3> out_subarray = out.subarray(block.begin(), block.end());

            linalg::add(in1_subarray, in2_subarray, out_subarray);
        },
        blocking.numBlocks()
    );
}


template <unsigned int N, class T1, class T2, class S1, class S2>
void add(MultiArrayView<N,T1,S1> const & in1,
         T2 value,
         MultiArrayView<N,T2,S2> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == out.shape(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> out_subarray = out.subarray(block.begin(), block.end());

            out_subarray = linalg::operator+(in1_subarray, value);
        },
        blocking.numBlocks()
    );
}


template <unsigned int N, class T1, class T2, class T3, class S1, class S2, class S3>
void sub(MultiArrayView<N,T1,S1> const & in1,
         MultiArrayView<N,T2,S2> const & in2,
         MultiArrayView<N,T3,S3> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == in2.shape(),
                        "In sub(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.shape() == out.shape(),
                        "In sub(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In sub(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> in2_subarray = in2.subarray(block.begin(), block.end());
            MultiArrayView<N,T3,S3> out_subarray = out.subarray(block.begin(), block.end());

            linalg::sub(in1_subarray, in2_subarray, out_subarray);
        },
        blocking.numBlocks()
    );
}


template <unsigned int N, class T1, class T2, class S1, class S2>
void sub(MultiArrayView<N,T1,S1> const & in1,
         T2 value,
         MultiArrayView<N,T2,S2> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == out.shape(),
                        "In sub(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In sub(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> out_subarray = out.subarray(block.begin(), block.end());

            out_subarray = linalg::operator-(in1_subarray, value);
        },
        blocking.numBlocks()
    );
}


template <unsigned int N, class T1, class T2, class T3, class S1, class S2, class S3>
void pmul(MultiArrayView<N,T1,S1> const & in1,
         MultiArrayView<N,T2,S2> const & in2,
         MultiArrayView<N,T3,S3> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == in2.shape(),
                        "In pmul(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.shape() == out.shape(),
                        "In pmul(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In pmul(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> in2_subarray = in2.subarray(block.begin(), block.end());
            MultiArrayView<N,T3,S3> out_subarray = out.subarray(block.begin(), block.end());

            linalg::pmul(in1_subarray, in2_subarray, out_subarray);
        },
        blocking.numBlocks()
    );
}


template <unsigned int N, class T1, class T2, class S1, class S2>
void smul(MultiArrayView<N,T1,S1> const & in1,
         T2 value,
         MultiArrayView<N,T2,S2> & out,
         BlockwiseOptions<N> const & options)
{
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;

    vigra_precondition(in1.shape() == out.shape(),
                        "In smul(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.shape(), options.getBlockShape()),
                        "In smul(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.shape(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArrayView<N,T1,S1> in1_subarray = in1.subarray(block.begin(), block.end());
            MultiArrayView<N,T2,S2> out_subarray = out.subarray(block.begin(), block.end());

            linalg::smul(in1_subarray, value, out_subarray);
        },
        blocking.numBlocks()
    );
}






template <unsigned int M, unsigned int N1, unsigned int N2, unsigned int N3,
          class T1, class T2, class T3, class TAG1, class TAG2, class TAG3>
void add(ChunkedArray<N1,T1,TAG1> const & in1,
         ChunkedArray<N2,T2,TAG2> const & in2,
         ChunkedArray<N3,T3,TAG3> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == in2.template minimalShape<M>(),
                        "In add(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> in2_cache(block.size());
            MultiArray<M,T3> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);
            in2.checkoutSubarray(block.begin(), in2_cache);

            linalg::add(in1_cache, in2_cache, out_cache);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


template <unsigned int M, unsigned int N1, unsigned int N2, class T1, class T2, class TAG1, class TAG2>
void add(ChunkedArray<N1,T1,TAG1> const & in1,
         T1 value,
         ChunkedArray<N2,T2,TAG2> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);

            out_cache = linalg::operator+(in1_cache, value);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


template <unsigned int M, unsigned int N1, unsigned int N2, unsigned int N3,
          class T1, class T2, class T3, class TAG1, class TAG2, class TAG3>
void sub(ChunkedArray<N1,T1,TAG1> const & in1,
         ChunkedArray<N2,T2,TAG2> const & in2,
         ChunkedArray<N3,T3,TAG3> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == in2.template minimalShape<M>(),
                        "In add(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> in2_cache(block.size());
            MultiArray<M,T3> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);
            in2.checkoutSubarray(block.begin(), in2_cache);

            linalg::sub(in1_cache, in2_cache, out_cache);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


template <unsigned int M, unsigned int N1, unsigned int N2, class T1, class T2, class TAG1, class TAG2>
void sub(ChunkedArray<N1,T1,TAG1> const & in1,
         T1 value,
         ChunkedArray<N2,T2,TAG2> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In add(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In add(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);

            out_cache = linalg::operator-(in1_cache, value);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


template <unsigned int M, unsigned int N1, unsigned int N2, unsigned int N3,
          class T1, class T2, class T3, class TAG1, class TAG2, class TAG3>
void pmul(ChunkedArray<N1,T1,TAG1> const & in1,
         ChunkedArray<N2,T2,TAG2> const & in2,
         ChunkedArray<N3,T3,TAG3> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == in2.template minimalShape<M>(),
                        "In pmul(...): Shapes of input arrays don't match.");
    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In pmul(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In pmul(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> in2_cache(block.size());
            MultiArray<M,T3> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);
            in2.checkoutSubarray(block.begin(), in2_cache);

            linalg::pmul(in1_cache, in2_cache, out_cache);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


template <unsigned int M, unsigned int N1, unsigned int N2, class T1, class T2, class TAG1, class TAG2>
void smul(ChunkedArray<N1,T1,TAG1> const & in1,
         T1 value,
         ChunkedArray<N2,T2,TAG2> & out,
         BlockwiseOptions<M> const & options)
{
    typedef MultiBlocking<M, MultiArrayIndex> Blocking;

    vigra_precondition(in1.template minimalShape<M>() == out.template minimalShape<M>(),
                        "In smul(...): Shapes of input and output array don't match.");
    vigra_precondition(allGreaterEqual(in1.template minimalShape<M>(), options.getBlockShape()),
                        "In smul(...): The array shapes must be greater or equal then the block shape.");

    const Blocking blocking(in1.template minimalShape<M>(), options.getBlockShape());

    parallel_foreach(options.getNumThreads(), blocking.blockBegin(), blocking.blockEnd(),
        [&](const int /*threadId*/, const typename Blocking::Block block){
            MultiArray<M,T1> in1_cache(block.size());
            MultiArray<M,T2> out_cache(block.size());

            in1.checkoutSubarray(block.begin(), in1_cache);

            linalg::smul(in1_cache, value, out_cache);

            out.commitSubarray(block.begin(), out_cache);
        },
        blocking.numBlocks()
    );
}


} // END NAMESPACE VIGRA

#endif // VIGRA_MULTI_ARRAY_BLOCKWISE_MATH_HXX
