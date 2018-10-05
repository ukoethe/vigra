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
#include <vigra/multi_array.hxx>
#include <vigra/multi_blockwise.hxx>
#include <vigra/multi_shape.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/unittest.hxx>
#include <vigra/error.hxx>

#include <cstddef>
#include <iostream>
#include "utils.hxx"

using namespace vigra;

void chunkedViewTest_basic()
{
    const int N = 3;
    typedef float T;
    Shape3 backend_shape(100, 60, 40);
    Shape3 chunk_shape(32, 16, 16);
    Shape3 view_shape(30, 20, 15);
    Shape3 view_start(10, 5, 3);
    Shape3 view_steps(2);
    Shape3 null_shape(0);

    ChunkedArrayCompressed<N, T> chunked_array1(backend_shape, chunk_shape);
    ChunkedArrayCompressed<N, T> chunked_array2(backend_shape, chunk_shape);

    MultiArray<N, T> array(backend_shape);
    fillRandom(array.begin(), array.end(), 256);
    chunked_array1.commitSubarray(null_shape, array);
    chunked_array2.commitSubarray(null_shape, array);

    ChunkedArrayView<N, T> chunked_view(chunked_array1, view_start, view_shape, view_steps);
    MultiArray<N, T> view(view_shape);
    fillRandom(view.begin(), view.end(), 256);

    chunked_view.commitSubarray(null_shape, view);
    for (auto k = 0; k < view_shape[2]; ++k)
        for (auto j = 0; j < view_shape[1]; ++j)
            for (auto i = 0; i < view_shape[0]; ++i)
                chunked_array2.setItem(view_start + Shape3(i, j, k)*view_steps,
                                       chunked_view.getItem(Shape3(i, j, k)));

    compare(chunked_array1.begin(), chunked_array1.end(), chunked_array2.begin());
}

void chunkedViewTest_minimal()
{
    const int N = 5;
    typedef float T;
    Shape5 backend_shape(50, 40, 30, 20, 10);
    Shape5 chunk_shape(16, 8, 4, 2, 2);
    Shape5 view_start(5, 4, 3, 2, 1);
    Shape5 view_shape(20, 1, 10, 8, 1);
    Shape5 view_steps(2);

    const int minimalDimension = 3;
    Shape3 minimalShape(20, 10, 8);
    Shape3 minimalChunkShape(16, 4, 2);

    ChunkedArrayCompressed<N, T> chunked_array(backend_shape, chunk_shape);
    ChunkedArrayView<N, T> chunked_view(chunked_array, view_start, view_shape, view_steps);

    vigra_assert(minimalDimension == chunked_view.minimalDimension(),
            "ERROR in chunkedviewTest_minimalShape(): minimalDimension is wrong");
    vigra_assert(minimalShape == chunked_view.template minimalShape<minimalDimension>(),
            "ERROR in chunkedviewTest_minimalShape(): minimalShape is wrong");
    vigra_assert(minimalChunkShape == chunked_view.template minimalChunkShape<minimalDimension>(),
            "ERROR in chunkedviewTest_minimalShape(): minimalChunkShape is wrong");

    MultiArray<minimalDimension, T> minimal_view(minimalShape);
    fillRandom(minimal_view.begin(), minimal_view.end(), 256);

    chunked_view.commitSubarray(Shape3(0), minimal_view);

    MultiArray<N, T> view(view_shape);
    chunked_view.checkoutSubarray(Shape5(0), view);

    compare(view.begin(), view.end(), minimal_view.begin());
}

void chunkedViewTest_bindChannel()
{
    const unsigned int N = 3;
    typedef float T;
    typedef TinyVector<T, N> T_VEC;
    Shape3 array_shape(100, 60, 40);
    Shape3 chunk_shape(32, 16, 16);
    Shape3 view_start(10, 5, 2);
    Shape3 view_shape(20, 10, 10);
    Shape3 view_steps(3);

    ChunkedArrayCompressed<N, T_VEC> chunked_array(array_shape, chunk_shape);
    ChunkedArrayView<N, T_VEC> chunked_view(chunked_array, view_start,
                                            view_shape, view_steps);
    MultiArray<N, T_VEC> array(array_shape);
    auto view = array.subarray(view_start, view_start + view_shape);

    auto chunked_view_b = chunked_view.bindElementChannel(1);
    auto view_b = view.bindElementChannel(1);
    fillRandom(view_b.begin(), view_b.end(), 256);

    chunked_view_b.checkoutSubarray(Shape3(0), view_b);
    compare(array.begin(), array.end(), chunked_array.begin());
}

void chunkedViewTest_assignment()
{
    const unsigned int N = 3;
    const unsigned int M = 2;
    typedef float T;
    Shape3 array_shape_3D(100, 60, 40);
    Shape3 chunk_shape_3D(32, 16, 16);
    Shape3 view_start(10, 6, 4);
    Shape3 view_shape(30, 1, 10);
    Shape3 view_steps(2);
    Shape2 array_shape_2D(30, 10);
    Shape2 chunk_shape_2D(8, 4);
    Shape2 block_shape(16, 8);

    BlockwiseOptions<M> opt;
    opt.blockShape(block_shape);

    ChunkedArrayCompressed<N, T> chunked_array_3D(array_shape_3D, chunk_shape_3D);
    ChunkedArrayCompressed<M, T> chunked_array_2D(array_shape_2D, chunk_shape_2D);
    ChunkedArrayView<N, T> view(chunked_array_3D, view_start, view_shape, view_steps);

    fillRandom(chunked_array_2D.begin(), chunked_array_2D.end(), 256);
    view.assign(chunked_array_2D, opt);

    for (auto j = 0; j < array_shape_2D[1]; ++j)
        for (auto i = 0; i < array_shape_2D[0]; ++i)
        {
            shouldEqual(chunked_array_2D.getItem(Shape2(i, j)),
                        chunked_array_3D.getItem(view_start + Shape3(i, 0, j)*view_steps));
        }
}

void chunkedViewTest_overlaps()
{
    const unsigned int N = 2;
    typedef float T;
    Shape2 backend_shape(400);
    Shape2 roi_shape_1(100);
    Shape2 roi_start_1(50);
    Shape2 roi_steps_1(3);
    Shape2 roi_shape_2(50);
    Shape2 roi_start_2(100);
    Shape2 roi_steps_2(2);
    Shape2 roi_shape_3(30);
    Shape2 roi_start_3(300);
    Shape2 roi_steps_3(2);

    ChunkedArrayCompressed<N, T> array_1(backend_shape);
    ChunkedArrayCompressed<N, T> array_2(backend_shape);
    ChunkedArrayROI<N, T> roi_11(array_1, roi_start_1, roi_shape_1);
    ChunkedArrayROI<N, T> roi_12(array_1, roi_start_2, roi_shape_2);
    ChunkedArrayROI<N, T> roi_13(array_1, roi_start_3, roi_shape_3);
    ChunkedArrayROI<N, T> roi_21(array_2, roi_start_1, roi_shape_1);

    vigra_assert(roi_11.overlaps(roi_12), "Failure in overlaps test 12o11");
    vigra_assert(!roi_11.overlaps(roi_21), "Failure in overlaps test !11o21");
    vigra_assert(!roi_13.overlaps(roi_12), "Failure in overlaps test !13o12");
}

struct ChunkedRoiTestSuite
: public test_suite
{
    ChunkedRoiTestSuite()
    : test_suite("ChunkedArrayView test")
    {
        add(testCase(chunkedViewTest_basic));
        add(testCase(chunkedViewTest_minimal));
        add(testCase(chunkedViewTest_bindChannel));
        add(testCase(chunkedViewTest_assignment));
        add(testCase(chunkedViewTest_overlaps));
    }
};

int main(int argc, char** argv)
{
    ChunkedRoiTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return failed != 0;
}
