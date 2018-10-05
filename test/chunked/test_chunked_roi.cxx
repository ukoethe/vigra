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

void chunkedRoiTest_basic()
{
    const int N = 3;
    typedef float T;
    Shape3 backend_shape(100, 60, 40);
    Shape3 chunk_shape(32, 16, 16);
    Shape3 roi_shape(60, 30, 25);
    Shape3 roi_start(10, 5, 3);
    Shape3 null_shape(0);

    ChunkedArrayCompressed<N, T> chunked_array1(backend_shape, chunk_shape);
    ChunkedArrayCompressed<N, T> chunked_array2(backend_shape, chunk_shape);

    MultiArray<N, T> array(backend_shape);
    fillRandom(array.begin(), array.end(), 256);
    chunked_array1.commitSubarray(null_shape, array);
    chunked_array2.commitSubarray(null_shape, array);

    ChunkedArrayROI<N, T> chunked_roi(chunked_array1, roi_start, roi_shape);
    MultiArray<N, T> roi(roi_shape);
    fillRandom(roi.begin(), roi.end(), 256);

    chunked_roi.commitSubarray(null_shape, roi);
    chunked_array2.commitSubarray(roi_start, roi);

    compare(chunked_array1.begin(), chunked_array1.end(), chunked_array2.begin());
}

void chunkedRoiTest_minimal()
{
    const int N = 5;
    typedef float T;
    Shape5 backend_shape(50, 40, 30, 20, 10);
    Shape5 chunk_shape(16, 8, 4, 2, 2);
    Shape5 roi_start(5, 4, 3, 2, 1);
    Shape5 roi_shape(25, 1, 15, 10, 1);

    const int minimalDimension = 3;
    Shape3 minimalShape(25, 15, 10);
    Shape3 minimalChunkShape(16, 4, 2);

    ChunkedArrayCompressed<N, T> chunked_array(backend_shape, chunk_shape);
    ChunkedArrayROI<N, T> chunked_roi(chunked_array, roi_start, roi_shape);

    vigra_assert(minimalDimension == chunked_roi.minimalDimension(),
            "ERROR in chunkedRoiTest_minimalShape(): minimalDimension is wrong");
    vigra_assert(minimalShape == chunked_roi.template minimalShape<minimalDimension>(),
            "ERROR in chunkedRoiTest_minimalShape(): minimalShape is wrong");
    vigra_assert(minimalChunkShape == chunked_roi.template minimalChunkShape<minimalDimension>(),
            "ERROR in chunkedRoiTest_minimalShape(): minimalChunkShape is wrong");

    MultiArray<minimalDimension, T> minimal_roi(minimalShape);
    fillRandom(minimal_roi.begin(), minimal_roi.end(), 256);

    chunked_roi.commitSubarray(Shape3(0), minimal_roi);

    MultiArray<N, T> roi(roi_shape);
    chunked_array.checkoutSubarray(roi_start, roi);

    compare(roi.begin(), roi.end(), minimal_roi.begin());
}

void chunkedRoiTest_bindChannel()
{
    const unsigned int N = 3;
    typedef float T;
    typedef TinyVector<T, N> T_VEC;
    Shape3 array_shape(100, 60, 40);
    Shape3 chunk_shape(32, 16, 16);

    ChunkedArrayCompressed<N, T_VEC> chunked_array(array_shape, chunk_shape);
    MultiArray<N, T_VEC> array(array_shape);

    auto chunked_roi = chunked_array.bindElementChannel(1);
    auto roi = array.bindElementChannel(1);
    fillRandom(roi.begin(), roi.end(), 256);

    chunked_roi.checkoutSubarray(Shape3(0), roi);
    compare(array.begin(), array.end(), chunked_array.begin());
}

void chunkedRoiTest_assignment()
{
    const unsigned int N = 3;
    const unsigned int M = 2;
    typedef float T;
    Shape3 array_shape_3D(100, 60, 40);
    Shape3 chunk_shape_3D(32, 16, 16);
    Shape3 roi_start(10, 6, 4);
    Shape3 roi_shape(50, 1, 20);
    Shape2 array_shape_2D(50, 20);
    Shape2 chunk_shape_2D(8, 4);
    Shape2 block_shape(16, 8);
    BlockwiseOptions<M> opt;
    opt.blockShape(block_shape);

    ChunkedArrayCompressed<N, T> chunked_array_3D(array_shape_3D, chunk_shape_3D);
    ChunkedArrayCompressed<M, T> chunked_array_2D(array_shape_2D, chunk_shape_2D);
    ChunkedArrayROI<N, T> roi(chunked_array_3D, roi_start, roi_shape);

    fillRandom(chunked_array_2D.begin(), chunked_array_2D.end(), 256);
    roi.assign(chunked_array_2D, opt);

    for (auto j = 0; j < array_shape_2D[1]; ++j)
        for (auto i = 0; i < array_shape_2D[0]; ++i)
        {
            shouldEqual(chunked_array_2D.getItem(Shape2(i, j)),
                        chunked_array_3D.getItem(roi_start + Shape3(i, 0, j)));
        }
}

void chunkedRoiTest_overlaps()
{
    const unsigned int N = 2;
    typedef float T;
    Shape2 backend_shape(400);
    Shape2 roi_shape_1(200);
    Shape2 roi_start_1(50);
    Shape2 roi_shape_2(150);
    Shape2 roi_start_2(100);
    Shape2 roi_shape_3(80);
    Shape2 roi_start_3(270);

    ChunkedArrayCompressed<N, T> array_1(backend_shape);
    ChunkedArrayCompressed<N, T> array_2(backend_shape);
    ChunkedArrayROI<N, T> roi_11(array_1, roi_start_1, roi_shape_1);
    ChunkedArrayROI<N, T> roi_12(array_1, roi_start_2, roi_shape_2);
    ChunkedArrayROI<N, T> roi_13(array_1, roi_start_3, roi_shape_3);
    ChunkedArrayROI<N, T> roi_21(array_2, roi_start_1, roi_shape_1);

    vigra_assert(roi_12.overlaps(roi_11), "Failure in overlaps test 12o11");
    vigra_assert(!roi_11.overlaps(roi_13), "Failure in overlaps test !11o13");
    vigra_assert(!roi_21.overlaps(roi_12), "Failure in overlaps test !21o12");
}

struct ChunkedRoiTestSuite
: public test_suite
{
    ChunkedRoiTestSuite()
    : test_suite("ChunkedArrayROI test")
    {
        add(testCase(chunkedRoiTest_basic));
        add(testCase(chunkedRoiTest_minimal));
        add(testCase(chunkedRoiTest_bindChannel));
        add(testCase(chunkedRoiTest_assignment));
        add(testCase(chunkedRoiTest_overlaps));
    }
};

int main(int argc, char** argv)
{
    ChunkedRoiTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return failed != 0;
}
