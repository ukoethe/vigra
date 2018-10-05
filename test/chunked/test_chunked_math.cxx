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
#include <vigra/blockwise_math.hxx>
#include <vigra/linear_algebra.hxx>
#include <vigra/multi_blockwise.hxx>
#include <vigra/unittest.hxx>
#include <vigra/error.hxx>

#include <cstddef>
#include <iostream>
#include "utils.hxx"

using namespace vigra;

void test_add()
{
    const unsigned int N = 2;
    typedef float T;
    Shape2 shape(200, 300);
    Shape2 chunk_shape(32);

    T value = 666;

    MultiArray<N, T> v_array_in1(shape);
    MultiArray<N, T> v_array_in2(shape);
    MultiArray<N, T> v_array_out(shape);

    fillRandom(v_array_in1.begin(), v_array_in1.end(), 256);
    fillRandom(v_array_in2.begin(), v_array_in2.end(), 256);

    ChunkedArrayCompressed<N, T> c_array_in1(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_in2(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_out(shape, chunk_shape);

    c_array_in1.commitSubarray(Shape2(0), v_array_in1);
    c_array_in2.commitSubarray(Shape2(0), v_array_in2);

    MultiArray<N, T> m_array_in1 = v_array_in1;
    MultiArray<N, T> m_array_in2 = v_array_in2;
    MultiArray<N, T> m_array_out(shape);

    BlockwiseOptions<N> opt;
    opt.blockShape(40);

    linalg::add(v_array_in1, v_array_in2, v_array_out);
    add(c_array_in1, c_array_in2, c_array_out, opt);
    add(m_array_in1, m_array_in2, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());

    v_array_out = linalg::operator+(v_array_in1, value);
    add(c_array_in1, value, c_array_out, opt);
    add(m_array_in1, value, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());
}

void test_sub()
{
    const unsigned int N = 2;
    typedef float T;
    Shape2 shape(200, 300);
    Shape2 chunk_shape(32);

    T value = 666;

    MultiArray<N, T> v_array_in1(shape);
    MultiArray<N, T> v_array_in2(shape);
    MultiArray<N, T> v_array_out(shape);

    fillRandom(v_array_in1.begin(), v_array_in1.end(), 256);
    fillRandom(v_array_in2.begin(), v_array_in2.end(), 256);

    ChunkedArrayCompressed<N, T> c_array_in1(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_in2(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_out(shape, chunk_shape);

    c_array_in1.commitSubarray(Shape2(0), v_array_in1);
    c_array_in2.commitSubarray(Shape2(0), v_array_in2);

    MultiArray<N, T> m_array_in1 = v_array_in1;
    MultiArray<N, T> m_array_in2 = v_array_in2;
    MultiArray<N, T> m_array_out(shape);

    BlockwiseOptions<N> opt;
    opt.blockShape(40);

    linalg::sub(v_array_in1, v_array_in2, v_array_out);
    sub(c_array_in1, c_array_in2, c_array_out, opt);
    sub(m_array_in1, m_array_in2, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());

    v_array_out = linalg::operator-(v_array_in1, value);
    sub(c_array_in1, value, c_array_out, opt);
    sub(m_array_in1, value, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());
}

void test_mul()
{
    const unsigned int N = 2;
    typedef float T;
    Shape2 shape(200, 300);
    Shape2 chunk_shape(32);

    T value = 666;

    MultiArray<N, T> v_array_in1(shape);
    MultiArray<N, T> v_array_in2(shape);
    MultiArray<N, T> v_array_out(shape);

    fillRandom(v_array_in1.begin(), v_array_in1.end(), 256);
    fillRandom(v_array_in2.begin(), v_array_in2.end(), 256);

    ChunkedArrayCompressed<N, T> c_array_in1(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_in2(shape, chunk_shape);
    ChunkedArrayCompressed<N, T> c_array_out(shape, chunk_shape);

    c_array_in1.commitSubarray(Shape2(0), v_array_in1);
    c_array_in2.commitSubarray(Shape2(0), v_array_in2);

    MultiArray<N, T> m_array_in1 = v_array_in1;
    MultiArray<N, T> m_array_in2 = v_array_in2;
    MultiArray<N, T> m_array_out(shape);

    BlockwiseOptions<N> opt;
    opt.blockShape(40);

    linalg::pmul(v_array_in1, v_array_in2, v_array_out);
    pmul(c_array_in1, c_array_in2, c_array_out, opt);
    pmul(m_array_in1, m_array_in2, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());

    linalg::smul(v_array_in1, value, v_array_out);
    smul(c_array_in1, value, c_array_out, opt);
    smul(m_array_in1, value, m_array_out, opt);

    compare(c_array_out.begin(), c_array_out.end(), v_array_out.begin());
    compare(m_array_out.begin(), m_array_out.end(), v_array_out.begin());
}

using namespace vigra;
struct ChunkedRoiTestSuite
: public test_suite
{
    ChunkedRoiTestSuite()
    : test_suite("chunked math test")
    {
        add(testCase(test_add));
        add(testCase(test_sub));
        add(testCase(test_mul));
    }
};

int main(int argc, char** argv)
{
    ChunkedRoiTestSuite test;
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return failed != 0;
}
