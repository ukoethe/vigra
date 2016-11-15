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

#include <vigra/multi_array_chunked.hxx>
#include <vigra/multi_blockwise.hxx>
#include <vigra/multi_chunked.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/multi_shape.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/unittest.hxx>

#include <iostream>
#include "utils.hxx"

using namespace std;
using namespace vigra;

#define MAKE_TEST_VIEW(TEST_NAME, FUNCTION_NAME, USES_OUTER_SCALE)                              \
template <unsigned int N, class T1, class T2>                                                   \
void TEST_NAME()                                                                                \
{                                                                                               \
    int max = 256;                                                                              \
    TinyVector<MultiArrayIndex,N> shape(300, 50, 300);                                                   \
    TinyVector<MultiArrayIndex,N> chunk_shape(32);                                              \
    TinyVector<MultiArrayIndex,N-1> block_shape(60);                                            \
    TinyVector<MultiArrayIndex,N> view_start(50);                                               \
    TinyVector<MultiArrayIndex,N> view_shape(100,1,100);                                        \
    TinyVector<MultiArrayIndex,N> view_steps(2);                                                \
    TinyVector<MultiArrayIndex,N> stride = view_steps * detail::defaultStride(shape);           \
\
    MultiArray<N,T1> source_array(shape);                                                       \
    MultiArray<N,T2> dest_array(shape);                                                         \
    auto source_array_view = MultiArrayView<N,T1>(view_shape, stride,                           \
            &source_array(view_start[0], view_start[1], view_start[2])).bind<1>(0);             \
    auto dest_array_view = MultiArrayView<N,T1>(view_shape, stride,                             \
            &dest_array(view_start[0], view_start[1], view_start[2])).bind<1>(0);               \
\
    ChunkedArrayLazy<N,T1> source_lazy(shape, chunk_shape);                                     \
    ChunkedArrayLazy<N,T2> dest_lazy(shape, chunk_shape);                                       \
    ChunkedArrayView<N,T1> source_lazy_view(source_lazy, view_start, view_shape, view_steps);   \
    ChunkedArrayView<N,T2> dest_lazy_view(dest_lazy, view_start, view_shape, view_steps);       \
\
    ChunkedArrayCompressed<N,T1> source_comp(shape, chunk_shape);                               \
    ChunkedArrayCompressed<N,T2> dest_comp(shape, chunk_shape);                                 \
    ChunkedArrayView<N,T1> source_comp_view(source_comp, view_start, view_shape, view_steps);   \
    ChunkedArrayView<N,T2> dest_comp_view(dest_comp, view_start, view_shape, view_steps);       \
\
    fillRandom(source_array.begin(), source_array.end(), max);                                  \
    source_lazy_view.commitSubarray(TinyVector<int,N-1>(0), source_array_view);                 \
    source_comp_view.commitSubarray(TinyVector<int,N-1>(0), source_array_view);                 \
\
    BlockwiseConvolutionOptions<N-1> opt;                                                       \
    opt.stdDev(2.0);                                                                            \
    opt.blockShape(block_shape);                                                                \
    if (opt.getNumThreads() > 4)                                                                \
        opt.numThreads(4);                                                                      \
    if (USES_OUTER_SCALE)                                                                       \
        opt.outerScale(1.0);                                                                    \
\
    FUNCTION_NAME(source_array_view, dest_array_view, opt);                                     \
    FUNCTION_NAME(source_lazy_view, dest_lazy_view, opt);                                       \
    FUNCTION_NAME(source_comp_view, dest_comp_view, opt);                                       \
\
    compare(dest_lazy.begin(), dest_lazy.end(), dest_array.begin());                            \
    compare(dest_comp.begin(), dest_comp.end(), dest_array.begin());                            \
}

#define MAKE_TEST_ROI(TEST_NAME, FUNCTION_NAME, USES_OUTER_SCALE)                               \
template <unsigned int N, class T1, class T2>                                                   \
void TEST_NAME()                                                                                \
{                                                                                               \
    int max = 256;                                                                              \
    TinyVector<MultiArrayIndex,N> shape(300);                                                   \
    TinyVector<MultiArrayIndex,N> chunk_shape(32);                                              \
    TinyVector<MultiArrayIndex,N-1> block_shape(60);                                            \
    TinyVector<MultiArrayIndex,N> roi_start(50);                                                \
    TinyVector<MultiArrayIndex,N> roi_shape(100,1,100);                                         \
\
    MultiArray<N,T1> source_array(shape);                                                       \
    MultiArray<N,T2> dest_array(shape);                                                         \
    auto source_array_roi = source_array.subarray(roi_start, roi_start + roi_shape).bind<1>(0); \
    auto dest_array_roi = dest_array.subarray(roi_start, roi_start + roi_shape).bind<1>(0);     \
\
    ChunkedArrayLazy<N,T1> source_lazy(shape, chunk_shape);                                     \
    ChunkedArrayLazy<N,T2> dest_lazy(shape, chunk_shape);                                       \
    ChunkedArrayROI<N,T1> source_lazy_roi(source_lazy, roi_start, roi_shape);                   \
    ChunkedArrayROI<N,T2> dest_lazy_roi(dest_lazy, roi_start, roi_shape);                       \
\
    ChunkedArrayCompressed<N,T1> source_comp(shape, chunk_shape);                               \
    ChunkedArrayCompressed<N,T2> dest_comp(shape, chunk_shape);                                 \
    ChunkedArrayROI<N,T1> source_comp_roi(source_comp, roi_start, roi_shape);                   \
    ChunkedArrayROI<N,T2> dest_comp_roi(dest_comp, roi_start, roi_shape);                       \
\
    fillRandom(source_array.begin(), source_array.end(), max);                                  \
    source_lazy_roi.commitSubarray(TinyVector<int,N-1>(0), source_array_roi);                   \
    source_comp_roi.commitSubarray(TinyVector<int,N-1>(0), source_array_roi);                   \
\
    BlockwiseConvolutionOptions<N-1> opt;                                                       \
    opt.stdDev(2.0);                                                                            \
    opt.blockShape(block_shape);                                                                \
    if (opt.getNumThreads() > 4)                                                                \
        opt.numThreads(4);                                                                      \
    if (USES_OUTER_SCALE)                                                                       \
        opt.outerScale(1.0);                                                                    \
\
    FUNCTION_NAME(source_array_roi, dest_array_roi, opt);                                       \
    FUNCTION_NAME(source_lazy_roi, dest_lazy_roi, opt);                                         \
    FUNCTION_NAME(source_comp_roi, dest_comp_roi, opt);                                         \
\
    compare(dest_lazy.begin(), dest_lazy.end(), dest_array.begin());                            \
    compare(dest_comp.begin(), dest_comp.end(), dest_array.begin());                            \
}

MAKE_TEST_VIEW(gaussianSmoothTestImpl,                   gaussianSmoothMultiArray,                   false)
MAKE_TEST_VIEW(gaussianGradientMagnitudeTestImpl,        gaussianGradientMagnitudeMultiArray,        false)
MAKE_TEST_VIEW(laplacianOfGaussianTestImpl,              laplacianOfGaussianMultiArray,              false)
MAKE_TEST_VIEW(hessianOfGaussianFirstEigenvalueTestImpl, hessianOfGaussianFirstEigenvalueMultiArray, false)
MAKE_TEST_VIEW(hessianOfGaussianLastEigenvalueTestImpl,  hessianOfGaussianLastEigenvalueMultiArray,  false)
MAKE_TEST_ROI(gaussianGradientTestImpl,                 gaussianGradientMultiArray,                 false)
MAKE_TEST_ROI(symmetricGradientTestImpl,                symmetricGradientMultiArray,                false)
MAKE_TEST_ROI(hessianOfGaussianEigenvaluesTestImpl,     hessianOfGaussianEigenvaluesMultiArray,     false)
MAKE_TEST_ROI(gaussianDivergenceTestImpl,               gaussianDivergenceMultiArray,               false)
MAKE_TEST_ROI(hessianOfGaussianTestImpl,                hessianOfGaussianMultiArray,                false)
MAKE_TEST_ROI(structureTensorTestImpl,                  structureTensorMultiArray,                  true)

#undef MAKE_TEST_VIEW
#undef MAKE_TEST_ROI

struct ChunkedConvolutionTestSuite : public test_suite
{
    ChunkedConvolutionTestSuite()
    : test_suite("chunked blockwise convolution test")
    {
        typedef float T1;
        typedef float T2;
        constexpr int N = 3;
        constexpr int M = 2;

        auto gaussianSmoothTest                     = gaussianSmoothTestImpl<N,T1,T2>;
        auto gaussianGradientMagnitudeTest          = gaussianGradientMagnitudeTestImpl<N,T1,T2>;
        auto laplacianOfGaussianTest                = laplacianOfGaussianTestImpl<N,T1,T2>;
        auto hessianOfGaussianFirstEigenvalueTest   = hessianOfGaussianFirstEigenvalueTestImpl<N,T1,T2>;
        auto hessianOfGaussianLastEigenvalueTest    = hessianOfGaussianLastEigenvalueTestImpl<N,T1,T2>;
        auto gaussianGradientTest                   = gaussianGradientTestImpl<N,T1,TinyVector<T2,M> >;
        auto symmetricGradientTest                  = symmetricGradientTestImpl<N,T1,TinyVector<T2,M> >;
        auto hessianOfGaussianEigenvaluesTest       = hessianOfGaussianEigenvaluesTestImpl<N,T1,TinyVector<T2,M> >;
        auto gaussianDivergenceTest                 = gaussianDivergenceTestImpl<N,TinyVector<T1,M>,T2>;
        auto hessianOfGaussianTest                  = hessianOfGaussianTestImpl<N,T1,TinyVector<T2,int(M*(M+1)/2)> >;
        auto structureTensorTest                    = structureTensorTestImpl<N,T1,TinyVector<T2,int(M*(M+1)/2)> >;

        add(testCase(gaussianSmoothTest));
        add(testCase(gaussianGradientMagnitudeTest));
        add(testCase(laplacianOfGaussianTest));
        add(testCase(hessianOfGaussianFirstEigenvalueTest));
        add(testCase(hessianOfGaussianLastEigenvalueTest));
        add(testCase(gaussianGradientTest));
        add(testCase(symmetricGradientTest));
        add(testCase(hessianOfGaussianEigenvaluesTest));
        add(testCase(gaussianDivergenceTest));
        add(testCase(hessianOfGaussianTest));
        add(testCase(structureTensorTest));
    }
};

int main(int argc, char** argv)
{
    ChunkedConvolutionTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    cout << test.report() << endl;

    return failed != 0;
}
