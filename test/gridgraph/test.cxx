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
#include "unittest.hxx"
#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph_neighborhoods.hxx>
#include <vigra/algorithm.hxx>

using namespace vigra;

struct NeighborhoodTests
{
    NeighborhoodTests()
    {}
    
    template <unsigned int N>
    void testDirectNeighborhood()
    {
        typedef typename MultiArrayShape<N>::type Shape;
        
        ArrayVector<ArrayVector<Shape> > neighborhood;
        ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
        ArrayVector<ArrayVector<int> > neighborIndexLookup;
        detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, DirectNeighborhood);
        
        int neighborCount = 2*N;
        shouldEqual(neighborhood[0].size(), neighborCount);
        
        Shape pos, neg, strides = cumprod(Shape(3)) / 3;
        for(int k=0; k<neighborCount; ++k)
        {
            shouldEqual(sum(abs(neighborhood[0][k])), 1); // test that it is a direct neighbor 
            
            if(k < neighborCount/2)
            {
                should(dot(strides, neighborhood[0][k]) < 0); // check that causal neighbors are first
                neg += neighborhood[0][k];                    // check that all causal neighbors are found
            }
            else
            {
                should(dot(strides, neighborhood[0][k]) > 0); // check that anti-causal neighbors are last
                pos += neighborhood[0][k];                    // check that all anti-causal neighbors are found
            }
            
            shouldEqual(neighborhood[0][k], -neighborhood[0][neighborCount-1-k]); // check index of opposite neighbor
        }
        
        shouldEqual(pos, Shape(1));   // check that all causal neighbors were found
        shouldEqual(neg, Shape(-1));  // check that all anti-causal neighbors were found
        
        // check border flags
        MultiArray<N, int> a(Shape(3));
        typename MultiArrayView<N, int>::const_iterator ai = static_cast<MultiArrayView<N, int> const &>(a).begin();
        
        shouldEqual(neighborExists.size(), (int)pow(2.0, (int)N*2));
        
        for(int k=0; k<a.size(); ++k, ++ai)
        {
            int borderType = ai.borderType();
            shouldEqual(neighborExists[borderType].size(), neighborCount);
            
            for(int j=0; j<neighborCount; ++j)
            {
                shouldEqual(a.isInside(ai.point()+neighborhood[0][j]), neighborExists[borderType][j]);
            }
        }
    }
    
    template <unsigned int N>
    void testIndirectNeighborhood()
    {
        typedef typename MultiArrayShape<N>::type Shape;
        
        ArrayVector<ArrayVector<Shape> > neighborhood;
        ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
        ArrayVector<ArrayVector<int> > neighborIndexLookup;
        detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, IndirectNeighborhood);
        
        MultiArray<N, int> a(Shape(3));
        Shape center(1), strides = cumprod(Shape(3)) / 3;
        a[center] = 1;              
        
        int neighborCount = (int)pow(3.0, (int)N) - 1;
        shouldEqual(neighborhood[0].size(), neighborCount);
        
        for(int k=0; k<neighborCount; ++k)
        {
            shouldEqual(abs(neighborhood[0][k]).maximum(), 1); // check that offset is at most 1 in any direction
                 
            if(k < neighborCount/2)
                should(dot(strides, neighborhood[0][k]) < 0); // check that causal neighbors are first
            else
                should(dot(strides, neighborhood[0][k]) > 0); // check that anti-causal neighbors are last

            shouldEqual(neighborhood[0][k], -neighborhood[0][neighborCount-1-k]); // check index of opposite neighbor
            
            a[center+neighborhood[0][k]] += 1;  // check that all neighbors are found
        }
        
          // check that all neighbors are found
         int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
         a.minmax(&min, &max);
        
        shouldEqual(min, 1);
        shouldEqual(max, 1);
        
        // check border flags
        typename MultiArrayView<N, int>::const_iterator ai = static_cast<MultiArrayView<N, int> const &>(a).begin();
        
        shouldEqual(neighborExists.size(), (int)pow(2.0, (int)N*2));
        
        for(int k=0; k<a.size(); ++k, ++ai)
        {
            int borderType = ai.borderType();
            shouldEqual(neighborExists[borderType].size(), neighborCount);
            
            for(int j=0; j<neighborCount; ++j)
            {
                shouldEqual(a.isInside(ai.point()+neighborhood[0][j]), neighborExists[borderType][j]);
            }
        }
    }
};

struct GridgraphTestSuite
: public vigra::test_suite
{
    GridgraphTestSuite()
    : vigra::test_suite("Gridgraph Test")
    {
        add(testCase(&NeighborhoodTests::testDirectNeighborhood<2>));
        add(testCase(&NeighborhoodTests::testDirectNeighborhood<3>));
        add(testCase(&NeighborhoodTests::testIndirectNeighborhood<2>));
        add(testCase(&NeighborhoodTests::testIndirectNeighborhood<3>));
    }
};

int main(int argc, char **argv)
{

    GridgraphTestSuite gridgraphTest;

    int failed = gridgraphTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << gridgraphTest.report() << std::endl;

    return (failed != 0);
}
