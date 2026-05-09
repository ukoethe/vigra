/************************************************************************/
/*                                                                      */
/*     Copyright 2016 by Ullrich Koethe and Kevin Kiefer                */
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

#ifndef VIGRA_CHUNKED_ARRAY_BLOCKWISE_CONVOLUTION_HXX
#define VIGRA_CHUNKED_ARRAY_BLOCKWISE_CONVOLUTION_HXX

#include "multi_array_chunked.hxx"
#include "multi_blocking.hxx"
#include "multi_blockwise.hxx"
#include "multi_convolution.hxx"
#include "threadpool.hxx"


namespace vigra {

namespace chunked_blockwise {

/**
    Helper function to create blockwise parallel filters for chunked arrays.
*/
template <unsigned int N, class T1, class T2, class FUNC, class C>
void chunkedBlockwiseCaller(
        const ChunkedArray<N,T1> & source,
        ChunkedArray<N,T2> & dest,
        FUNC & func,
        const MultiBlocking<N,C> & blocking,
        const typename MultiBlocking<N,C>::Shape & borderWidth,
        const BlockwiseConvolutionOptions<N> & opt
)
{
    typedef typename MultiBlocking<N,C>::BlockWithBorder BlockWithBorder;

    auto begin = blocking.blockWithBorderBegin(borderWidth);
    auto end = blocking.blockWithBorderEnd(borderWidth);

    parallel_foreach(opt.getNumThreads(), begin, end,
        [&](const int /*threadId*/, const BlockWithBorder bwb)
        {
            // copy input of the block into a new allocated array
            MultiArray<N,T1> tmp(bwb.border().end() - bwb.border().begin());
            source.checkoutSubarray(bwb.border().begin(), tmp);

            // get the output as new allocated array
            MultiArray<N,T2> tmpOut(bwb.core().end() - bwb.core().begin());

            func(tmp, tmpOut, bwb.localCore().begin(), bwb.localCore().end());

            // copy output into destination
            dest.commitSubarray(bwb.core().begin(), tmpOut);
        },
        blocking.numBlocks()
    );
}

} // END NAMESPACE chunked_blockwise

/**
    Overload the functions listed below for chunked arrays.
    NOTE: Even if the MultiArrayView version may work in place
    the ChunkedArray overload does not.
*/
#define VIGRA_CHUNKED_BLOCKWISE(FUNCTOR, FUNCTION, ORDER, USES_OUTER_SCALE)                 \
template <unsigned int N, class T1, class T2>                                               \
void FUNCTION(                                                                              \
        const ChunkedArray<N,T1> & source,                                                  \
        ChunkedArray<N,T2> & dest,                                                          \
        BlockwiseConvolutionOptions<N> const & opt                                          \
)                                                                                           \
{                                                                                           \
    typedef MultiBlocking<N, MultiArrayIndex> Blocking;                                     \
    typedef typename Blocking::Shape Shape;                                                 \
\
    CompareChunkedArrays<N,T1,T2>(source, dest);                                            \
\
    const Shape border = blockwise::getBorder(opt, ORDER, USES_OUTER_SCALE);                \
    const Blocking blocking(source.shape(), opt.template getBlockShapeN<N>());              \
\
    BlockwiseConvolutionOptions<N> subOpt(opt);                                             \
    subOpt.subarray(Shape(0), Shape(0));                                                    \
\
    blockwise::FUNCTOR<N> func(subOpt);                                                     \
    chunked_blockwise::chunkedBlockwiseCaller(source, dest, func, blocking, border, opt);   \
}

// Reuse the blockwise functors from \<vigra/multi_blockwise.hxx\>
VIGRA_CHUNKED_BLOCKWISE(GaussianSmoothFunctor,                   gaussianSmoothMultiArray,                   0, false);
VIGRA_CHUNKED_BLOCKWISE(GaussianGradientFunctor,                 gaussianGradientMultiArray,                 1, false);
VIGRA_CHUNKED_BLOCKWISE(SymmetricGradientFunctor,                symmetricGradientMultiArray,                1, false);
VIGRA_CHUNKED_BLOCKWISE(GaussianDivergenceFunctor,               gaussianDivergenceMultiArray,               1, false);
VIGRA_CHUNKED_BLOCKWISE(HessianOfGaussianFunctor,                hessianOfGaussianMultiArray,                2, false);
VIGRA_CHUNKED_BLOCKWISE(HessianOfGaussianEigenvaluesFunctor,     hessianOfGaussianEigenvaluesMultiArray,     2, false);
VIGRA_CHUNKED_BLOCKWISE(HessianOfGaussianFirstEigenvalueFunctor, hessianOfGaussianFirstEigenvalueMultiArray, 2, false);
VIGRA_CHUNKED_BLOCKWISE(HessianOfGaussianLastEigenvalueFunctor,  hessianOfGaussianLastEigenvalueMultiArray,  2, false);
VIGRA_CHUNKED_BLOCKWISE(LaplacianOfGaussianFunctor,              laplacianOfGaussianMultiArray,              2, false);
VIGRA_CHUNKED_BLOCKWISE(GaussianGradientMagnitudeFunctor,        gaussianGradientMagnitudeMultiArray,        1, false);
VIGRA_CHUNKED_BLOCKWISE(StructureTensorFunctor,                  structureTensorMultiArray,                  1, true);

#undef VIGRA_CHUNKED_BLOCKWISE


// Alternative name for backward compatibility.
template <unsigned int N, class T1, class T2>
inline void gaussianGradientMagnitude(
        ChunkedArray<N,T1> const & source,
        ChunkedArray<N,T2> & dest,
        BlockwiseConvolutionOptions<N> const & opt
)
{
    gaussianGradientMagnitudeMultiArray(source, dest, opt);
}

} // END NAMESPACE vigra

#endif // VIGRA_CHUNKED_ARRAY_BLOCKWISE_CONVOLUTION_HXX
