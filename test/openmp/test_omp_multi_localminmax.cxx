/************************************************************************/
/*                                                                      */
/*              Copyright 2012-2013 by Ullrich Koethe                   */
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

#define VIGRA_CHECK_BOUNDS
#include "vigra/unittest.hxx"
#include <vigra/multi_shape.hxx>
#include <vigra/multi_iterator.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/multi_localminmax.hxx>
#include <vigra/algorithm.hxx>
#include <vigra/timing.hxx>

#ifdef WITH_BOOST_GRAPH
#  include <boost/graph/graph_concepts.hpp>
#endif


using namespace vigra;

template <unsigned int N>
struct GridGraphAlgorithmTests
{
    int max_threads;

    typedef typename MultiArrayShape<N>::type Shape;

    GridGraphAlgorithmTests()
#ifdef OPENMP
    : max_threads(std::max(omp_get_num_procs() / 2, 2))
#else
    : max_threads(1)
#endif
    {
#ifdef OPENMP
        omp_set_num_threads(max_threads);
#endif
    }

    template <class DirectedTag, NeighborhoodType NType>
    void testLocalMinMax()
    {
        typedef GridGraph<N, DirectedTag> Graph;

        Graph g(Shape(3), NType);
        typename Graph::template NodeMap<int> src(g), dest(g);

        src[Shape(1)] = 1;

        should(1 == boost_graph::localMinMaxGraph(g, src, dest, 1, -9999, std::greater<int>()));

        shouldEqualSequence(src.begin(), src.end(), dest.begin());

        dest.init(0);

        should(1 == lemon_graph::localMinMaxGraph(g, src, dest, 1, -9999, std::greater<int>()));

        shouldEqualSequence(src.begin(), src.end(), dest.begin());
    }
};

template <unsigned int N>
struct GridgraphTestSuiteN
: public vigra::test_suite
{
    GridgraphTestSuiteN()
    : vigra::test_suite((std::string("Gridgraph Test Dimension ") + vigra::asString(N)).c_str())
    {
        add(testCase((&GridGraphAlgorithmTests<N>::template testLocalMinMax<undirected_tag, DirectNeighborhood>)));
    }
};

struct GridgraphTestSuite
: public vigra::test_suite
{
    GridgraphTestSuite()
#ifdef WITH_BOOST_GRAPH
    : vigra::test_suite("Gridgraph BGL Test")
#else
    : vigra::test_suite("Gridgraph Test")
#endif
    {
        USETICTOC

        TIC
        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<2>));
        TOC

        TIC
        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<3>));
        TOC

        TIC
        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<4>));
        TOC
    }
};

int main(int argc, char **argv)
{

    GridgraphTestSuite gridgraphTest;

    int failed = gridgraphTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << gridgraphTest.report() << std::endl;

    return (failed != 0);
}
