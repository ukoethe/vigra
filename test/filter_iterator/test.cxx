/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
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
#include <vigra/filter_iterator.hxx>
#include <vigra/unittest.hxx>
#include <vector>
#include <numeric>
#include <iterator>

using namespace vigra;

struct FilterIteratorTests
{
    template <typename ITER0, typename ITER1>
    void test_filter_read_mod2(ITER0 in_begin, ITER0 in_end, ITER1 expected_begin)
    {
        typedef typename std::iterator_traits<ITER0>::value_type value_type;
        std::vector<value_type> out;
        auto filter = [](int x) { return x % 2 == 0; };
        auto begin = make_filter_iterator(filter, in_begin, in_end);
        auto end = make_filter_iterator(filter, in_end, in_end);
        for (auto it = begin; it != end; ++it)
        {
            out.push_back(*it);
        }
        shouldEqualSequence(out.begin(), out.end(), expected_begin);
    }

    template <typename ITER0, typename ITER1>
    void test_filter_write_mod2(ITER0 in_begin, ITER0 in_end, ITER1 expected_begin)
    {
        typedef typename std::iterator_traits<ITER0>::value_type value_type;
        std::vector<value_type> out(in_begin, in_end);
        auto filter = [](int x) { return x % 2 == 0; };
        auto begin = make_filter_iterator(filter, out.begin(), out.end());
        auto end = make_filter_iterator(filter, out.end(), out.end());
        for (auto it = begin; it != end; ++it)
        {
            *it += 100;
        }
        shouldEqualSequence(out.begin(), out.end(), expected_begin);
    }

    void test_filter_iterator_read()
    {
        {
            int in[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
            int out_expected[] = {0, 2, 4, 6, 8};
            test_filter_read_mod2(in, in+10, out_expected);
        }
        {
            int in[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
            int out_expected[] = {0, 2, 4, 6, 8};
            std::vector<int> v_in(in, in+9);
            test_filter_read_mod2(v_in.begin(), v_in.end(), out_expected);
        }
        {
            int in[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
            int out_expected[] = {2, 4, 6, 8};
            std::vector<int> v_in(in, in+9);
            test_filter_read_mod2(v_in.cbegin(), v_in.cend(), out_expected);
        }
        {
            int in[] = {1, 2, 3, 4, 5, 6, 7, 8};
            int out_expected[] = {2, 4, 6, 8};
            test_filter_read_mod2(in, in+8, out_expected);
        }
    }

    void test_filter_iterator_write()
    {
        {
            int in[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
            int out_expected[] = {100, 1, 102, 3, 104, 5, 106, 7, 108, 9};
            test_filter_write_mod2(in, in+10, out_expected);
        }
        {
            int in[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
            int out_expected[] = {100, 1, 102, 3, 104, 5, 106, 7, 108};
            test_filter_write_mod2(in, in+9, out_expected);
        }
        {
            int in[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
            int out_expected[] = {1, 102, 3, 104, 5, 106, 7, 108, 9};
            test_filter_write_mod2(in, in+9, out_expected);
        }
        {
            int in[] = {1, 2, 3, 4, 5, 6, 7, 8};
            int out_expected[] = {1, 102, 3, 104, 5, 106, 7, 108};
            test_filter_write_mod2(in, in+8, out_expected);
        }
    }
};

struct FilterIteratorTestSuite : public test_suite
{
    FilterIteratorTestSuite()
        :
        test_suite("FilterIterator test")
    {
        add(testCase(&FilterIteratorTests::test_filter_iterator_read));
        add(testCase(&FilterIteratorTests::test_filter_iterator_write));
    }
};

int main(int argc, char** argv)
{
    FilterIteratorTestSuite filter_iterator_test;
    int failed = filter_iterator_test.run(testsToBeExecuted(argc, argv));
    std::cout << filter_iterator_test.report() << std::endl;
    return (failed != 0);
}
