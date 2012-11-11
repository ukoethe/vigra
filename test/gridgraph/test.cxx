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

template <unsigned int N>
struct NeighborhoodTests
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    ArrayVector<Shape> neighborOffsets;
    ArrayVector<ArrayVector<bool> > neighborExists;
    ArrayVector<ArrayVector<Shape> > relativeOffsets, backOffsets, forwardOffsets;
    ArrayVector<ArrayVector<GridGraphEdgeDescriptor<N> > > edgeDescrOffsets, backEdgeDescrOffsets, forwardEdgeDescrOffsets;
    ArrayVector<ArrayVector<MultiArrayIndex> > neighborIndices, backIndices, forwardIndices;

    NeighborhoodTests()
    {}
    
    void testVertexIterator()
    {
        typedef typename MultiArrayShape<N>::type Shape;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i+Shape(1);
            MultiCoordinateIterator<N> vi(s), viend = vi.getEndIterator();
            
            MultiArray<N, int> vertex_map(s);

            for(; vi != viend; ++vi)
            {
                should(vi.isValid() && !vi.atEnd());
                vertex_map[*vi] += 1;
            }
            
            should(!vi.isValid() && vi.atEnd());
            
            // check that all vertice are found
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            vertex_map.minmax(&min, &max);
            
            shouldEqual(min, 1);
            shouldEqual(max, 1);
        }
    }   
    
    void testDirectNeighborhood()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, DirectNeighborhood);
        
        static const unsigned int neighborCount = 2*N;
        shouldEqual(neighborOffsets.size(), neighborCount);
        shouldEqual((GridGraphMaxDegree<N, DirectNeighborhood>::value), neighborCount);
        shouldEqual(gridGraphMaxDegree(N, DirectNeighborhood), neighborCount);
        
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
        
        shouldEqual(neighborExists.size(), (MetaPow<2, 2*N>::value));
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
    
    void testIndirectNeighborhood()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, IndirectNeighborhood);
        
        MultiArray<N, int> a(Shape(3));
        Shape center(1), strides = cumprod(Shape(3)) / 3;
        a[center] = 1;              
        
        static const unsigned int neighborCount = MetaPow<3, N>::value - 1;
        shouldEqual(neighborOffsets.size(), neighborCount);
        shouldEqual((GridGraphMaxDegree<N, IndirectNeighborhood>::value), neighborCount);
        shouldEqual(gridGraphMaxDegree(N, IndirectNeighborhood), neighborCount);
        
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
        
        shouldEqual(neighborExists.size(), (MetaPow<2, 2*N>::value));
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
    
    template <NeighborhoodType NType>
    void testNeighborhoodIterator()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, true, true, true);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, backOffsets, backEdgeDescrOffsets, backIndices, true, true, false);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, forwardOffsets, forwardEdgeDescrOffsets, forwardIndices, true, false, true);
        
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
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                
                {
                    GridGraphNeighborIterator<N> ni(backOffsets[borderType], backIndices[borderType], vi.point()),
                                                 nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size()/2; ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                
                {
                    GridGraphNeighborIterator<N> ni(forwardOffsets[borderType], forwardIndices[borderType], vi.point()),
                                                 nend = ni.getEndIterator();
                    
                    for(int k=neighborExists[borderType].size()/2; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point()+neighborOffsets[k], *ni);
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testOutEdgeIteratorDirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, true, true, true);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, backOffsets, backEdgeDescrOffsets, backIndices, true, true, false);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, forwardOffsets, forwardEdgeDescrOffsets, forwardIndices, true, false, true);

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
                    GridGraphOutEdgeIterator<N> ni(edgeDescrOffsets[borderType], neighborIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point(), ni->vertexDescriptor());
                            shouldEqual(k, ni->edgeIndex());
                            shouldEqual(k, ni.neighborIndex());
                            should(!ni->isReversed());
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                        
                {
                    GridGraphOutEdgeIterator<N> ni(backEdgeDescrOffsets[borderType], backIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size()/2; ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point(), ni->vertexDescriptor());
                            shouldEqual(k, ni->edgeIndex());
                            shouldEqual(k, ni.neighborIndex());
                            should(!ni->isReversed());
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                
                {
                    GridGraphOutEdgeIterator<N> ni(forwardEdgeDescrOffsets[borderType], forwardIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=neighborExists[borderType].size()/2; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(vi.point(), ni->vertexDescriptor());
                            shouldEqual(k, ni->edgeIndex());
                            shouldEqual(k, ni.neighborIndex());
                            should(!ni->isReversed());
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testOutEdgeIteratorUndirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, false, true, true);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, backOffsets, backEdgeDescrOffsets, backIndices, false, true, false);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, forwardOffsets, forwardEdgeDescrOffsets, forwardIndices, false, false, true);

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
                    GridGraphOutEdgeIterator<N> ni(edgeDescrOffsets[borderType], neighborIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(k, ni.neighborIndex());
                            if(k < neighborExists[borderType].size() / 2)
                            {
                                shouldEqual(vi.point(), ni->vertexDescriptor());
                                shouldEqual(k, ni->edgeIndex());
                                should(!ni->isReversed());
                            }
                            else
                            {
                                shouldEqual(vi.point()+neighborOffsets[k], ni->vertexDescriptor());
                                shouldEqual(k, (int)neighborOffsets.size() - ni->edgeIndex() - 1);
                                should(ni->isReversed());
                            }
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                        
                {
                    GridGraphOutEdgeIterator<N> ni(backEdgeDescrOffsets[borderType], backIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=0; k<neighborExists[borderType].size()/2; ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(k, ni.neighborIndex());
                            shouldEqual(vi.point(), ni->vertexDescriptor());
                            shouldEqual(k, ni->edgeIndex());
                            should(!ni->isReversed());
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
                
                {
                    GridGraphOutEdgeIterator<N> ni(forwardEdgeDescrOffsets[borderType], forwardIndices[borderType], vi.point()),
                                                nend = ni.getEndIterator();
                    
                    for(int k=neighborExists[borderType].size()/2; k<neighborExists[borderType].size(); ++k)
                    {
                        if(neighborExists[borderType][k])
                        {
                            should(ni.isValid() && !ni.atEnd());
                            shouldEqual(k, ni.neighborIndex());
                            shouldEqual(vi.point()+neighborOffsets[k], ni->vertexDescriptor());
                            shouldEqual(k, (int)neighborOffsets.size() - ni->edgeIndex() - 1);
                            should(ni->isReversed());
                            ++ni;
                        }
                    }
                    should(ni == nend);
                    should(ni.atEnd() && !ni.isValid());
                }
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testEdgeIteratorDirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, true, true, true);

        // check neighborhoods at ROI border
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        typedef typename MultiArrayShape<N+1>::type EdgeMapShape;
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            
            EdgeMapShape es(gridGraphMaxDegree(N, NType));
            es.template subarray<0,N>() = s;
            MultiArray<N+1, int> edge_map(es);
            
            GridGraphEdgeIterator<N> ni(edgeDescrOffsets, neighborIndices, s),
                                     nend = ni.getEndIterator();
            
            for(; ni != nend; ++ni)
            {
                should(ni.isValid() && !ni.atEnd());
                edge_map[*ni] += 1;
            }
            
            should(!ni.isValid() && ni.atEnd());
                
            // check neighborhood of all pixels
            MultiCoordinateIterator<N> vi(s), viend = vi.getEndIterator();
            for(; vi != viend; ++vi)
            {
                int borderType = vi.borderType();
                es.template subarray<0,N>() = *vi;
                
                for(es[N]=0; es[N]<(MultiArrayIndex)neighborExists[borderType].size(); ++es[N])
                {
                    if(neighborExists[borderType][es[N]])
                        shouldEqual(edge_map[es], 1);
                    else
                        shouldEqual(edge_map[es], 0);
                }
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testEdgeIteratorUndirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, backOffsets, backEdgeDescrOffsets, backIndices, false, true, false);

        // check neighborhoods at ROI border
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();
        typedef typename MultiArrayShape<N+1>::type EdgeMapShape;
        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            
            EdgeMapShape es(gridGraphMaxDegree(N, NType) / 2);
            es.template subarray<0,N>() = s;
            MultiArray<N+1, int> edge_map(es);
            
            GridGraphEdgeIterator<N> ni(backEdgeDescrOffsets, backIndices, s),
                                     nend = ni.getEndIterator();
            
            for(; ni != nend; ++ni)
            {
                should(ni.isValid() && !ni.atEnd());
                edge_map[*ni] += 1;
            }
            
            should(!ni.isValid() && ni.atEnd());
                
            // check neighborhood of all pixels
            MultiCoordinateIterator<N> vi(s), viend = vi.getEndIterator();
            for(; vi != viend; ++vi)
            {
                int borderType = vi.borderType();
                es.template subarray<0,N>() = *vi;
                
                for(es[N]=0; es[N]<(MultiArrayIndex)neighborExists[borderType].size()/2; ++es[N])
                {
                    if(neighborExists[borderType][es[N]])
                        shouldEqual(edge_map[es], 1);
                    else
                        shouldEqual(edge_map[es], 0);
                }
            }
        }
    }
};

template <unsigned int N>
struct GridGraphTests
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    template <class DirectedTag, NeighborhoodType NType>
    void testCounts()
    {
        static const bool directed = IsSameType<DirectedTag, vigragraph::directed_tag>::value;
        
        GridGraph<N, DirectedTag> g1(Shape(1), NType);
        
        shouldEqual(g1.num_vertices(), 1);
        shouldEqual(g1.num_edges(), 0);
        shouldEqual(g1.maxDegree(), gridGraphMaxDegree(N, NType));
        
        GridGraph<N, DirectedTag> g3(Shape(3), NType);
        
        // int expectedEdgeCount = directed
                                  // ? 2*gridGraphMaxDegree(N, NType)
                                  // : gridGraphMaxDegree(N, NType);
        shouldEqual(g3.num_vertices(), (MetaPow<3, N>::value));
        // shouldEqual(g3.num_edges(), 0);
        shouldEqual(g3.maxDegree(), gridGraphMaxDegree(N, NType));
    }
};

struct GridgraphTestSuite
: public vigra::test_suite
{
    GridgraphTestSuite()
    : vigra::test_suite("Gridgraph Test")
    {
        add(testCase(&NeighborhoodTests<2>::testVertexIterator));
        add(testCase(&NeighborhoodTests<3>::testVertexIterator));
        
        add(testCase(&NeighborhoodTests<2>::testDirectNeighborhood));
        add(testCase(&NeighborhoodTests<3>::testDirectNeighborhood));
        add(testCase(&NeighborhoodTests<2>::testIndirectNeighborhood));
        add(testCase(&NeighborhoodTests<3>::testIndirectNeighborhood));
        
        add(testCase(&NeighborhoodTests<2>::testNeighborhoodIterator<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testNeighborhoodIterator<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<2>::testNeighborhoodIterator<IndirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testNeighborhoodIterator<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<2>::testOutEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testOutEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<2>::testOutEdgeIteratorDirected<IndirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testOutEdgeIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<2>::testOutEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testOutEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<2>::testOutEdgeIteratorUndirected<IndirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testOutEdgeIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<2>::testEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<2>::testEdgeIteratorDirected<IndirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testEdgeIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<2>::testEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<2>::testEdgeIteratorUndirected<IndirectNeighborhood>));
        add(testCase(&NeighborhoodTests<3>::testEdgeIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase((&GridGraphTests<2>::testCounts<vigragraph::directed_tag, IndirectNeighborhood>)));
    }
};

int main(int argc, char **argv)
{

    GridgraphTestSuite gridgraphTest;

    int failed = gridgraphTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << gridgraphTest.report() << std::endl;

    return (failed != 0);
}
