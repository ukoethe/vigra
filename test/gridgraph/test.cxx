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
            
            shouldEqual(edge_map.template sum<int>(), gridGraphEdgeCount(s, NType, true));
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
            
            shouldEqual(edge_map.template sum<int>(), gridGraphEdgeCount(s, NType, false));
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
        using namespace vigragraph;
        static const bool directed = IsSameType<DirectedTag, vigragraph::directed_tag>::value;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            GridGraph<N, DirectedTag> g(s, NType);
        
            shouldEqual(g.shape(), s);
            shouldEqual(g.isDirected(), directed);
            shouldEqual(g.num_vertices(), prod(s));
            shouldEqual(num_vertices(g), prod(s));
            shouldEqual(g.num_edges(), gridGraphEdgeCount(s, NType, directed));
            shouldEqual(num_edges(g), gridGraphEdgeCount(s, NType, directed));
            shouldEqual(g.maxDegree(), gridGraphMaxDegree(N, NType));
        }
    }
    
    template <class DirectedTag, NeighborhoodType NType>
    void testVertexIterator()
    {
        using namespace vigragraph;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            MultiArray<N, int> vertexMap(s);
            int count = 0;
        
            typename Graph::vertex_iterator j = g.get_vertex_iterator(), 
                                            end = g.get_vertex_end_iterator();
                                            
            should(j == vertices(g).first);
            should(end == vertices(g).second);
            for(; j != end; ++j, ++count)
            {
                should(j.isValid() && !j.atEnd());
                shouldEqual(g.out_degree(j), g.out_degree(*j));
                shouldEqual(g.out_degree(j), out_degree(*j, g));
                shouldEqual(g.forward_degree(j), g.forward_degree(*j));
                shouldEqual(g.back_degree(j), g.back_degree(*j));
                shouldEqual(g.forward_degree(j) + g.back_degree(j), g.out_degree(j));
                shouldEqual(g.in_degree(j), g.out_degree(j));
                put(vertexMap, *j, get(vertexMap, *j) + 1); // same as: vertexMap[*j] += 1;
            }
            should(!j.isValid() && j.atEnd());
            
            // check that all vertices are found exactly once
            shouldEqual(count, g.num_vertices());
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            vertexMap.minmax(&min, &max);
            
            shouldEqual(min, 1);
            shouldEqual(max, 1);
        }
    }
    
    template <class DirectedTag, NeighborhoodType NType>
    void testNeighborIterator()
    {
        using namespace vigragraph;
        
        static const bool directed = IsSameType<DirectedTag, vigragraph::directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            MultiArray<N, int> vertexMap(s);
            MultiArray<N+1, int> edgeMap(g.edge_propmap_shape());
            
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), directed ? g.maxDegree() : g.halfMaxDegree());
            
            int totalCount = 0;
        
            typename Graph::vertex_iterator j = g.get_vertex_iterator(), 
                                            end = g.get_vertex_end_iterator();
            for(; j != end; ++j)
            {
                typename Graph::neighbor_vertex_iterator n = g.get_neighbor_vertex_iterator(j), 
                                                         nend = g.get_neighbor_vertex_end_iterator(j);
                typename Graph::out_edge_iterator        e = g.get_out_edge_iterator(j), 
                                                         eend = g.get_out_edge_end_iterator(j);

                should(n == g.get_neighbor_vertex_iterator(*j));
                should(nend == g.get_neighbor_vertex_end_iterator(*j));
                should(n == adjacent_vertices(*j, g).first);
                should(nend == adjacent_vertices(*j, g).second);
                should(n == adjacent_vertices_at_iterator(j, g).first);
                should(nend == adjacent_vertices_at_iterator(j, g).second);
                should(e == g.get_out_edge_iterator(*j));
                should(eend == g.get_out_edge_end_iterator(*j));
                should(e == out_edges(*j, g).first);
                should(eend == out_edges(*j, g).second);

                int count = 0;
                for(; n != nend; ++n, ++e, ++count)
                {
                    should(*n != *j);
                    should(n.isValid() && !n.atEnd());
                    should(e.isValid() && !e.atEnd());
                    should(e != eend);
                    shouldEqual(source(*e, g), *j);
                    shouldEqual(target(*e, g), *n);
                    vertexMap[*n] += 1;
                    edgeMap[*e] += 1;
                }
                should(!n.isValid() && n.atEnd());
                should(!e.isValid() && e.atEnd());
                should(e == eend);
                shouldEqual(count, g.out_degree(j));
                
                totalCount += count;
            }
            
            // check that all neighbors are found
            if(!directed)
                totalCount /= 2;
            shouldEqual(totalCount, g.num_edges());
            
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            edgeMap.minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, g.num_edges() ? !directed ? 2 : 1 : 0);

            j = g.get_vertex_iterator();
            for(; j != end; ++j)
            {
                shouldEqual(vertexMap[*j], g.in_degree(j));
                if(directed)
                    shouldEqual(edgeMap.bindInner(*j).template sum<int>(), g.out_degree(j));
                else
                    shouldEqual(edgeMap.bindInner(*j).template sum<int>(), 2*g.back_degree(j));
            }
        }
    }
    
    template <class DirectedTag, NeighborhoodType NType>
    void testEdgeIterator()
    {
        using namespace vigragraph;
        
        static const bool directed = IsSameType<DirectedTag, vigragraph::directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            MultiArray<N, int> sourceVertexMap(s),targetVertexMap(s);
            MultiArray<N+1, int> edgeMap(g.edge_propmap_shape());
            
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), directed ? g.maxDegree() : g.halfMaxDegree());
            
            typename Graph::edge_iterator e = g.get_edge_iterator(),
                                          eend = g.get_edge_end_iterator();

            should(e == edges(g).first);
            should(eend == edges(g).second);
            
            int count = 0;
            for(; e != eend; ++e, ++count)
            {
                should(e.isValid() && !e.atEnd());
                put(edgeMap, *e, get(edgeMap, *e) + 1); // same as: edgeMap[*e] += 1;
                sourceVertexMap[source(*e, g)] += 1;
                targetVertexMap[target(*e, g)] += 1;
            }
            should(!e.isValid() && e.atEnd());
            
            // check that all neighbors are found
            shouldEqual(count, g.num_edges());
            
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            edgeMap.minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, g.num_edges() ? 1 : 0);
            
            MultiCoordinateIterator<N> j(s), end = j.getEndIterator();
            for(; j != end; ++j)
            {
                if(directed)
                {
                    shouldEqual(edgeMap.bindInner(*j).template sum<int>(), g.out_degree(j));
                    shouldEqual(sourceVertexMap[*j], g.out_degree(j));
                    shouldEqual(targetVertexMap[*j], g.out_degree(j));
                }
                else
                {
                    shouldEqual(edgeMap.bindInner(*j).template sum<int>(), g.back_degree(j));
                    shouldEqual(sourceVertexMap[*j], g.back_degree(j));
                    shouldEqual(targetVertexMap[*j], g.forward_degree(j));
                }
            }
        }
    }
};

template <unsigned int N>
struct GridgraphTestSuiteN
: public vigra::test_suite
{
    GridgraphTestSuiteN()
    : vigra::test_suite((std::string("Gridgraph Test Dimension ") + vigra::asString(N)).c_str())
    {
        add(testCase(&NeighborhoodTests<N>::testVertexIterator));
        
        add(testCase(&NeighborhoodTests<N>::testDirectNeighborhood));
        add(testCase(&NeighborhoodTests<N>::testIndirectNeighborhood));
        
        add(testCase(&NeighborhoodTests<N>::template testNeighborhoodIterator<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testNeighborhoodIterator<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testOutEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testOutEdgeIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testOutEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testOutEdgeIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testEdgeIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testEdgeIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testEdgeIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testEdgeIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase((&GridGraphTests<N>::template testCounts<vigragraph::directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testCounts<vigragraph::undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testCounts<vigragraph::directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testCounts<vigragraph::undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testVertexIterator<vigragraph::directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<vigragraph::undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<vigragraph::directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<vigragraph::undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<vigragraph::directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<vigragraph::undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<vigragraph::directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<vigragraph::undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<vigragraph::directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<vigragraph::undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<vigragraph::directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<vigragraph::undirected_tag, DirectNeighborhood>)));
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
        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<2>));
        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<3>));
//        add(VIGRA_TEST_SUITE(GridgraphTestSuiteN<4>));
    }
};

int main(int argc, char **argv)
{

    GridgraphTestSuite gridgraphTest;

    int failed = gridgraphTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << gridgraphTest.report() << std::endl;

    return (failed != 0);
}
