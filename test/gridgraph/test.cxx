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
#include <vigra/multi_shape.hxx>
#include <vigra/multi_iterator.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph.hxx>
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
        
        ArrayVector<Shape> neighborOffsets;
        ArrayVector<ArrayVector<bool> > neighborExists;
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, DirectNeighborhood);
        
        int neighborCount = 2*N;
        shouldEqual(neighborOffsets.size(), neighborCount);
        
        Shape pos, neg, strides = cumprod(Shape(3)) / 3;
        for(int k=0; k<neighborCount; ++k)
        {
            shouldEqual(sum(abs(neighborOffsets[k])), 1); // test that it is a direct neighbor 
            
            if(k < neighborCount/2)
            {
                should(dot(strides, neighborOffsets[k]) < 0); // check that causal neighbors are first
                neg += neighborOffsets[k];                    // check that all causal neighbors are found
            }
            else
            {
                should(dot(strides, neighborOffsets[k]) > 0); // check that anti-causal neighbors are last
                pos += neighborOffsets[k];                    // check that all anti-causal neighbors are found
            }
            
            shouldEqual(neighborOffsets[k], -neighborOffsets[neighborCount-1-k]); // check index of opposite neighbor
        }
        
        shouldEqual(pos, Shape(1));   // check that all causal neighbors were found
        shouldEqual(neg, Shape(-1));  // check that all anti-causal neighbors were found
        
        shouldEqual(neighborExists.size(), (int)pow(2.0, (int)N*2));
        MultiArray<1, unsigned char> checkNeighborCodes(Shape1(neighborExists.size()), (unsigned char)0);

        // check neighborhoods at ROI border
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        MultiArray<N, int> a(Shape(3));
        typedef typename MultiArray<N, int>::view_type View;
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            View va = a.subarray(Shape(), *i+Shape(1)); 
            
            // check neighborhood of all pixels
            typename View::iterator vi = va.begin(), viend = vi.getEndIterator();
            for(; vi != viend; ++vi)
            {
                int borderType = vi.borderType();
                
                shouldEqual(neighborExists[borderType].size(), neighborCount);
                checkNeighborCodes[borderType] = 1;
                
                for(int k=0; k<neighborCount; ++k)
                {
                    // check that neighbors are correctly marked as inside or outside in neighborExists
                    shouldEqual(va.isInside(vi.point()+neighborOffsets[k]), neighborExists[borderType][k]);
                }
            }
        }
        
        should(checkNeighborCodes.all()); // check that all possible neighborhoods have been tested
    }
    
    template <unsigned int N>
    void testIndirectNeighborhood()
    {
        typedef typename MultiArrayShape<N>::type Shape;
        
        ArrayVector<Shape> neighborOffsets;
        ArrayVector<ArrayVector<bool> > neighborExists;
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, IndirectNeighborhood);
        
        MultiArray<N, int> a(Shape(3));
        Shape center(1), strides = cumprod(Shape(3)) / 3;
        a[center] = 1;              
        
        int neighborCount = (int)pow(3.0, (int)N) - 1;
        shouldEqual(neighborOffsets.size(), neighborCount);
        
        for(int k=0; k<neighborCount; ++k)
        {
            shouldEqual(abs(neighborOffsets[k]).maximum(), 1); // check that offset is at most 1 in any direction
                 
            if(k < neighborCount/2)
                should(dot(strides, neighborOffsets[k]) < 0); // check that causal neighbors are first
            else
                should(dot(strides, neighborOffsets[k]) > 0); // check that anti-causal neighbors are last

            shouldEqual(neighborOffsets[k], -neighborOffsets[neighborCount-1-k]); // check index of opposite neighbor
            
            a[center+neighborOffsets[k]] += 1;  // check that all neighbors are found
        }
        
          // check that all neighbors are found
         int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
         a.minmax(&min, &max);
        
        shouldEqual(min, 1);
        shouldEqual(max, 1);
        
        shouldEqual(neighborExists.size(), (int)pow(2.0, (int)N*2));
        MultiArray<1, unsigned char> checkNeighborCodes(Shape1(neighborExists.size()), (unsigned char)0);

        // check neighborhoods at ROI border
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        typedef typename MultiArray<N, int>::view_type View;
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            View va = a.subarray(Shape(), *i +Shape(1));
            
            // check neighborhood of all pixels
            typename View::iterator vi = va.begin(), viend = vi.getEndIterator();
            for(; vi != viend; ++vi)
            {
                int borderType = vi.borderType();
                
                shouldEqual(neighborExists[borderType].size(), neighborCount);
                checkNeighborCodes[borderType] = 1;
                
                for(int k=0; k<neighborCount; ++k)
                {
                    // check that neighbors are correctly marked as inside or outside in neighborExists
                    shouldEqual(va.isInside(vi.point()+neighborOffsets[k]), neighborExists[borderType][k]);
                }
            }
        }
        
        should(checkNeighborCodes.all()); // check that all possible neighborhoods have been tested
    }
    
    template <unsigned int N, NeighborhoodType NType>
    void testNeighborhoodIterator()
    {
        typedef typename MultiArrayShape<N>::type Shape;
        
        ArrayVector<Shape> neighborOffsets;
        ArrayVector<ArrayVector<bool> > neighborExists;
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        
        ArrayVector<ArrayVector<Shape> > relativeOffsets, backOffsets, forwardOffsets;
        ArrayVector<ArrayVector<MultiArrayIndex> > neighborIndices, backIndices, forwardIndices;
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, neighborIndices, true, true);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, backOffsets, backIndices, true, false);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, forwardOffsets, forwardIndices, false, true);
        
        // check neighborhoods at ROI border
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        MultiArray<N, int> a(Shape(3));
        typedef typename MultiArray<N, int>::view_type View;
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            View va = a.subarray(Shape(), *i+Shape(1)); 
            
            // check neighborhood of all pixels
            typename View::iterator vi = va.begin(), viend = vi.getEndIterator();
            for(; vi != viend; ++vi)
            {
                int borderType = vi.borderType();
                
                {
                    GridGraphNeighborIterator<N> ni(relativeOffsets[borderType], neighborIndices[borderType], vi.point()),
                                                 nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                }
                
                {
                    GridGraphNeighborIterator<N> ni(backOffsets[borderType], backIndices[borderType], vi.point()),
                                                 nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size()/2; ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                }
                
                {
                    GridGraphNeighborIterator<N> ni(forwardOffsets[borderType], forwardIndices[borderType], vi.point()),
                                                 nend = ni.getEndIterator();
                    
                    for(int k=neighborExists[borderType].size()/2; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                }
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
        add(testCase((&NeighborhoodTests::testNeighborhoodIterator<2, DirectNeighborhood>)));
        add(testCase((&NeighborhoodTests::testNeighborhoodIterator<3, DirectNeighborhood>)));
        add(testCase((&NeighborhoodTests::testNeighborhoodIterator<2, IndirectNeighborhood>)));
        add(testCase((&NeighborhoodTests::testNeighborhoodIterator<3, IndirectNeighborhood>)));
    }
};

int main(int argc, char **argv)
{

    GridgraphTestSuite gridgraphTest;

    int failed = gridgraphTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << gridgraphTest.report() << std::endl;

    return (failed != 0);
}
