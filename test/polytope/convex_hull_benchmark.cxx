/************************************************************************/
/*                                                                      */
/*             Copyright 2011-2012 by Ullrich Koethe                    */
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

#include <iostream>
#include <chrono>

#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/polytope.hxx>
#include <vigra/accumulator.hxx>

namespace chrono = std::chrono;

namespace vigra
{

using namespace acc;

struct ConvexHullBenchmark
{

    typedef chrono::steady_clock   clock_type;

    template <unsigned int N>
    void testTypical()
    {
        std::cout << "# Benchmark for typical case with dim = " << N << ". ";
        std::cout << "All time measures in ms." << std::endl;
        std::cout << "# size, min, mean, max, std" << std::endl;
        for (int size = 1; size < 10000; size *= 2)
        {
            ArrayVector<double> times = constructPolytopeTypical<N>(size, 32);
            AccumulatorChain<
                    MultiArrayIndex,
                    Select<Minimum, Maximum, Mean, StdDev> > acc_chain;
            extractFeatures(times.begin(), times.end(), acc_chain);
            std::cout
                    << size << ", "
                    << get<Minimum>(acc_chain) << ", "
                    << get<Mean   >(acc_chain) << ", "
                    << get<Maximum>(acc_chain) << ", "
                    << get<StdDev >(acc_chain) << std::endl;
        }
    }

    template <unsigned int N>
    void testWorst()
    {
        std::cout << "# Benchmark for worst case with dim = " << N << ". ";
        std::cout << "All time measures in ms." << std::endl;
        std::cout << "# size, min, mean, max, std" << std::endl;
        for (int size = 1; size < 10000; size *= 2)
        {
            ArrayVector<double> times = constructPolytopeWorst<N>(size, 32);
            AccumulatorChain<
                    MultiArrayIndex,
                    Select<Minimum, Maximum, Mean, StdDev> > acc_chain;
            extractFeatures(times.begin(), times.end(), acc_chain);
            std::cout
                    << size << ", "
                    << get<Minimum>(acc_chain) << ", "
                    << get<Mean   >(acc_chain) << ", "
                    << get<Maximum>(acc_chain) << ", "
                    << get<StdDev >(acc_chain) << std::endl;
        }
    }

    template <unsigned int N>
    ArrayVector<double> constructPolytopeTypical(int size, int iterations) const
    {
        ArrayVector<double> ret;
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            ret.push_back(constructPolytopeTypical<N>(size));
        }
        return ret;
    }

    template <unsigned int N>
    double constructPolytopeTypical(int size) const
    {
        clock_type::time_point start = clock_type::now();
        // Construct the base polytope
        ConvexPolytope<N, double> poly;
        TinyVector<double, N> vec;
        poly.addVertex(vec);
        for (int n = 0; n < N; n++)
        {
            vec[n] = 1.;
            if (n > 0)
            {
                vec[n-1] = 0.;
            }
            poly.addVertex(vec);
        }
        poly.close();

        // Add the vertices
        for (int n = 0; n < size; n++)
        {
            do
            {
                for (int dim = 0; dim < N; dim++)
                {
                    vec[dim] = (2*rand() - 1)/static_cast<double>(RAND_MAX);
                }
            }
            while (vec.magnitude() > 1.);
            poly.addExtremeVertex(vec);
        }
        clock_type::time_point stop = clock_type::now();
        return chrono::duration_cast<chrono::microseconds>(stop - start).count();
    }

    template <unsigned int N>
    ArrayVector<double> constructPolytopeWorst(int size, int iterations) const
    {
        ArrayVector<double> ret;
        for (int iteration = 0; iteration < iterations; iteration++)
        {
            ret.push_back(constructPolytopeWorst<N>(size));
        }
        return ret;
    }

    template <unsigned int N>
    double constructPolytopeWorst(int size) const
    {
        clock_type::time_point start = clock_type::now();
        // Construct the base polytope
        ConvexPolytope<N, double> poly;
        TinyVector<double, N> vec;
        poly.addVertex(vec);
        for (int n = 0; n < N; n++)
        {
            vec[n] = 1.;
            if (n > 0)
            {
                vec[n-1] = 0.;
            }
            poly.addVertex(vec);
        }
        poly.close();

        // Add the vertices
        for (int n = 0; n < size; n++)
        {
            for (int dim = 0; dim < N; dim++)
            {
                vec[dim] = (2*rand() - 1)/static_cast<double>(RAND_MAX);
            }
            vec /= norm(vec);
            poly.addExtremeVertex(vec);
        }
        clock_type::time_point stop = clock_type::now();
        return chrono::duration_cast<chrono::microseconds>(stop - start).count();
    }
};

struct ConvexHullBenchmarkSuite : public test_suite
{
    ConvexHullBenchmarkSuite()
    : test_suite("ConvexHullBenchmarkSuite")
    {
        // add(testCase(&ConvexHullBenchmark::testTypical<2>));
        // add(testCase(&ConvexHullBenchmark::testTypical<3>));
        // add(testCase(&ConvexHullBenchmark::testTypical<4>));
        add(testCase(&ConvexHullBenchmark::testWorst<2>));
        add(testCase(&ConvexHullBenchmark::testWorst<3>));
        add(testCase(&ConvexHullBenchmark::testWorst<4>));
    }
};

} // namespace vigra

int main(int argc, char** argv)
{
    vigra::ConvexHullBenchmarkSuite benchmark;
    const int failed = benchmark.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << benchmark.report() << std::endl;

    return failed != 0;
}
