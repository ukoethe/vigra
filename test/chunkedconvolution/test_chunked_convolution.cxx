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
#include <vigra/tinyvector.hxx>
#include <vigra/unittest.hxx>

#include <iostream>
#include "utils.hxx"

using namespace std;
using namespace vigra;

#define MAKE_TEST(TEST_NAME, FUNCTION_NAME, USES_OUTER_SCALE)           \
template <unsigned int N, class T1, class T2>                           \
void TEST_NAME()                                                        \
{                                                                       \
    int max = 256;                                                      \
    TinyVector<int,N> shape(512,512);                                   \
    TinyVector<int,N> chunk_shape(32,32);                               \
    TinyVector<int,N> block_shape(30,30);                               \
\
    MultiArray<N,T1> source_array(shape);                               \
    MultiArray<N,T2> dest_array(shape);                                 \
\
    ChunkedArrayLazy<N,T1> source_lazy(shape, chunk_shape);             \
    ChunkedArrayLazy<N,T2> dest_lazy(shape, chunk_shape);               \
\
    ChunkedArrayCompressed<N,T1> source_comp(shape, chunk_shape);       \
    ChunkedArrayCompressed<N,T2> dest_comp(shape, chunk_shape);         \
\
    fillRandom(source_array.begin(), source_array.end(), max);          \
    source_lazy.commitSubarray(TinyVector<int,N>(0,0), source_array);   \
    source_comp.commitSubarray(TinyVector<int,N>(0,0), source_array);   \
\
    BlockwiseConvolutionOptions<N> opt;                                 \
    opt.stdDev(2.0);                                                    \
    opt.blockShape(block_shape);                                        \
    if (opt.getNumThreads() > 4)                                        \
        opt.numThreads(4);                                              \
    if (USES_OUTER_SCALE)                                               \
        opt.outerScale(1.0, 1.0);                                       \
\
    FUNCTION_NAME(source_array, dest_array, opt);                       \
    FUNCTION_NAME(source_lazy, dest_lazy, opt);                         \
    FUNCTION_NAME(source_comp, dest_comp, opt);                         \
\
    compare(dest_lazy.begin(), dest_lazy.end(), dest_array.begin());    \
    compare(dest_comp.begin(), dest_comp.end(), dest_array.begin());    \
}

MAKE_TEST(gaussianSmoothTestImpl,                   gaussianSmoothMultiArray,                   false)
MAKE_TEST(gaussianGradientMagnitudeTestImpl,        gaussianGradientMagnitudeMultiArray,        false)
MAKE_TEST(laplacianOfGaussianTestImpl,              laplacianOfGaussianMultiArray,              false)
MAKE_TEST(hessianOfGaussianFirstEigenvalueTestImpl, hessianOfGaussianFirstEigenvalueMultiArray, false)
MAKE_TEST(hessianOfGaussianLastEigenvalueTestImpl,  hessianOfGaussianLastEigenvalueMultiArray,  false)
MAKE_TEST(gaussianGradientTestImpl,                 gaussianGradientMultiArray,                 false)
MAKE_TEST(symmetricGradientTestImpl,                symmetricGradientMultiArray,                false)
MAKE_TEST(hessianOfGaussianEigenvaluesTestImpl,     hessianOfGaussianEigenvaluesMultiArray,     false)
MAKE_TEST(gaussianDivergenceTestImpl,               gaussianDivergenceMultiArray,               false)
MAKE_TEST(hessianOfGaussianTestImpl,                hessianOfGaussianMultiArray,                false)
MAKE_TEST(structureTensorTestImpl,                  structureTensorMultiArray,                  true)

struct ChunkedConvolutionTestSuite : public test_suite
{
    ChunkedConvolutionTestSuite()
    : test_suite("chunked blockwise convolution test")
    {
        typedef float T1;
        typedef float T2;
        constexpr int N = 2;

        auto gaussianSmoothTest                     = gaussianSmoothTestImpl<N,T1,T2>;
        auto gaussianGradientMagnitudeTest          = gaussianGradientMagnitudeTestImpl<N,T1,T2>;
        auto laplacianOfGaussianTest                = laplacianOfGaussianTestImpl<N,T1,T2>;
        auto hessianOfGaussianFirstEigenvalueTest   = hessianOfGaussianFirstEigenvalueTestImpl<N,T1,T2>;
        auto hessianOfGaussianLastEigenvalueTest    = hessianOfGaussianLastEigenvalueTestImpl<N,T1,T2>;
        auto gaussianGradientTest                   = gaussianGradientTestImpl<N,T1,TinyVector<T2,N> >;
        auto symmetricGradientTest                  = symmetricGradientTestImpl<N,T1,TinyVector<T2,N> >;
        auto hessianOfGaussianEigenvaluesTest       = hessianOfGaussianEigenvaluesTestImpl<N,T1,TinyVector<T2,N> >;
        auto gaussianDivergenceTest                 = gaussianDivergenceTestImpl<N,TinyVector<T1,N>,T2>;
        auto hessianOfGaussianTest                  = hessianOfGaussianTestImpl<N,T1,TinyVector<T2,int(N*(N+1)/2)> >;
        auto structureTensorTest                    = structureTensorTestImpl<N,T1,TinyVector<T2,int(N*(N+1)/2)> >;

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
