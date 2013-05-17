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
#include <vigra/multi_localminmax.hxx>
#include <vigra/algorithm.hxx>

#ifdef WITH_BOOST_GRAPH
#  include <boost/graph/graph_concepts.hpp>
#endif


using namespace vigra;

template <unsigned int N>
struct NeighborhoodTests
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    ArrayVector<Shape> neighborOffsets;
    ArrayVector<ArrayVector<bool> > neighborExists;
    ArrayVector<ArrayVector<Shape> > relativeOffsets, backOffsets, forwardOffsets;
    ArrayVector<ArrayVector<GridGraphArcDescriptor<N> > > edgeDescrOffsets, backEdgeDescrOffsets, forwardEdgeDescrOffsets;
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
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, backIndices, true);
        
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
                    GridGraphNeighborIterator<N> ni(relativeOffsets[borderType], 
                                                    neighborIndices[borderType], 
                                                    backIndices[borderType], 
                                                    vi.point()),
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
                    GridGraphNeighborIterator<N, true> ni(relativeOffsets[borderType], 
                                                          neighborIndices[borderType], 
                                                          backIndices[borderType], 
                                                          vi.point()),
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
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testOutArcIteratorDirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, backIndices, true);

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
                    GridGraphOutArcIterator<N, false> ni(edgeDescrOffsets[borderType], 
                                                         neighborIndices[borderType], 
                                                         backIndices[borderType], 
                                                         vi.point()),
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
                    GridGraphOutArcIterator<N, true> ni(edgeDescrOffsets[borderType], 
                                                        neighborIndices[borderType], 
                                                        backIndices[borderType], 
                                                        vi.point()),
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
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testOutArcIteratorUndirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, backIndices, false);

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
                    GridGraphOutArcIterator<N> ni(edgeDescrOffsets[borderType], 
                                                  neighborIndices[borderType], 
                                                  backIndices[borderType], 
                                                  vi.point()),
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
                    GridGraphOutArcIterator<N, true> ni(edgeDescrOffsets[borderType], 
                                                        neighborIndices[borderType], 
                                                        backIndices[borderType], 
                                                        vi.point()),
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
            }
        }
    }
    
    template <NeighborhoodType NType>
    void testArcIteratorDirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, backIndices, true);

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
            
            GridGraphArcIterator<N, false> ni(edgeDescrOffsets, neighborIndices, backIndices, s),
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
    void testArcIteratorUndirected()
    {
        detail::makeArrayNeighborhood(neighborOffsets, neighborExists, NType);
        detail::computeNeighborOffsets(neighborOffsets, neighborExists, relativeOffsets, edgeDescrOffsets, neighborIndices, backIndices, false);

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
            
            GridGraphArcIterator<N, true> ni(edgeDescrOffsets, neighborIndices, backIndices, s),
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
    void testBasics()
    {
        using namespace boost;
        typedef GridGraph<N, DirectedTag> G;
        
        static const bool directed = IsSameType<DirectedTag, directed_tag>::value;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            G g(s, NType);
            
#ifdef WITH_LEMON
            shouldEqual(directed, !lemon::UndirectedTagIndicator<G>::value);
#endif
            
            should(typename G::vertex_descriptor(lemon::INVALID) == lemon::INVALID);
            should(typename G::edge_descriptor(lemon::INVALID) == lemon::INVALID);
        
            shouldEqual(g.shape(), s);
            shouldEqual(g.isDirected(), directed);
            shouldEqual(g.num_vertices(), prod(s));
            shouldEqual(num_vertices(g), prod(s));
            shouldEqual(g.num_edges(), gridGraphEdgeCount(s, NType, directed));
            shouldEqual(num_edges(g), gridGraphEdgeCount(s, NType, directed));
            shouldEqual(g.maxDegree(), gridGraphMaxDegree(N, NType));
            
            shouldEqual(g.num_vertices(), g.nodeNum());
            shouldEqual(g.num_edges(), g.edgeNum());
            shouldEqual(g.arcNum(), directed ? g.edgeNum() : 2*g.edgeNum());
        }
#ifdef WITH_BOOST_GRAPH
        BOOST_CONCEPT_ASSERT(( boost::GraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::BidirectionalGraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::AdjacencyGraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::EdgeListGraphConcept<G> ));
        BOOST_CONCEPT_ASSERT(( boost::AdjacencyMatrixConcept<G> ));
#endif
    }
    
    template <class DirectedTag, NeighborhoodType NType>
    void testVertexIterator()
    {
        using namespace boost;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            typename Graph::template NodeMap<int> vertexMap(g);
            int count = 0;
        
            typename Graph::vertex_iterator j = g.get_vertex_iterator(), 
                                            end = g.get_vertex_end_iterator();
            typename Graph::NodeIt lj(g);
                                            
            should(j == vertices(g).first);
            should(end == vertices(g).second);
            should(j == typename Graph::vertex_iterator(g));
            for(; j != end; ++j, ++lj, ++count)
            {
                should(j.isValid() && !j.atEnd());
                should(j == lj);
                should(j != lemon::INVALID);
                should(!(j == lemon::INVALID));
                
                shouldEqual(j.scanOrderIndex(), g.id(j));
                shouldEqual(j.scanOrderIndex(), g.id(*j));
                shouldEqual(*j, g.nodeFromId(g.id(*j)));
                
                shouldEqual(g.out_degree(j), g.out_degree(*j));
                shouldEqual(g.out_degree(j), out_degree(*j, g));
                shouldEqual(g.forward_degree(j), g.forward_degree(*j));
                shouldEqual(g.back_degree(j), g.back_degree(*j));
                shouldEqual(g.forward_degree(j) + g.back_degree(j), g.out_degree(j));
                shouldEqual(g.in_degree(j), g.out_degree(j));
                shouldEqual(g.degree(j), g.isDirected() ? 2*g.out_degree(j) : g.out_degree(j));
                put(vertexMap, *j, get(vertexMap, *j) + 1); // same as: vertexMap[*j] += 1;
            }
            should(!j.isValid() && j.atEnd());
            should(j == lj);
            should(j == lemon::INVALID);
            should(!(j != lemon::INVALID));
            
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
        using namespace boost;
        
        static const bool directed = IsSameType<DirectedTag, directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        typedef typename Graph::Node Node;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            MultiArray<N, int> vertexMap(s);
            typename Graph::template EdgeMap<int> edgeMap(g);
            
            typename Graph::template ArcMap<int> arcIdMap(g);            
            typename Graph::template EdgeMap<int> edgeIdMap(g);            

            linearSequence(arcIdMap.begin(), arcIdMap.end());
            linearSequence(edgeIdMap.begin(), edgeIdMap.end());
            
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), g.maxUniqueDegree());
            shouldEqual((arcIdMap.shape().template subarray<0, N>()), s);
            shouldEqual(arcIdMap.shape(N), g.maxDegree());
            shouldEqual(edgeMap.shape(), edgeIdMap.shape());
            
            int totalCount = 0;
        
            typename Graph::vertex_iterator j = g.get_vertex_iterator(), 
                                            end = g.get_vertex_end_iterator();
            for(; j != end; ++j)
            {
                typename Graph::neighbor_vertex_iterator n = g.get_neighbor_vertex_iterator(j), 
                                                         nend = g.get_neighbor_vertex_end_iterator(j);
                typename Graph::out_edge_iterator        e = g.get_out_edge_iterator(j), 
                                                         eend = g.get_out_edge_end_iterator(j);
                typename Graph::in_edge_iterator         i = g.get_in_edge_iterator(j), 
                                                         iend = g.get_in_edge_end_iterator(j);
                typename Graph::IncEdgeIt                le(g, j);
                typename Graph::OutArcIt                 la(g, j);
                typename Graph::InArcIt                  li(g, j);

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
                should(i == in_edges(*j, g).first);
                should(iend == in_edges(*j, g).second);

                int count = 0;
                for(; n != nend; ++n, ++e, ++i, ++la, ++le, ++li, ++count)
                {
                    should(*n != *j);
                    should(n.isValid() && !n.atEnd());
                    should(n != lemon::INVALID);
                    should(!(n == lemon::INVALID));
                    should(e.isValid() && !e.atEnd());
                    should(e == la);
                    should(e == le);
                    should(e != eend);
                    should(e != lemon::INVALID);
                    should(!(e == lemon::INVALID));
                    should(i != lemon::INVALID);
                    should(!(i == lemon::INVALID));
                    should(i == li);
                    
                    shouldEqual(source(*e, g), *j);
                    shouldEqual(target(*e, g), *n);
                    vertexMap[*n] += 1;
                    edgeMap[*e] += 1;
                    
                    shouldEqual(source(*i, g), *n);
                    shouldEqual(target(*i, g), *j);
                    
                    typename Graph::Arc a  = g.findArc(*j, *n),
                                        oa = g.findArc(*n, *j);
                    should(a == *la);
                    shouldEqual(source(oa, g), *n);
                    shouldEqual(target(oa, g), *j);
                    shouldEqual(g.baseNode(la), *j);
                    shouldEqual(g.runningNode(la), *n);
                    shouldEqual(g.oppositeArc(a), oa);
                    
                    typename Graph::Edge ge = g.findEdge(*j, *n),
                                         oe = g.findEdge(*n, *j);
                    should(ge == *le);
                    shouldEqual(g.baseNode(le), *j);
                    shouldEqual(g.runningNode(le), *n);
                    
                    typename Graph::Node u = g.u(oe), 
                                         v = g.v(oe);
                    if(directed)
                    {
                        should(g.direction(*la));
                        shouldEqual(g.direct(*le, true), *la);
                        shouldEqual(g.direct(*le, false), g.oppositeArc(*la));
                    }
                    else
                    {
                        should(oe == *le);
                        if(g.id(n) < g.id(j))
                        {
                            std::swap(u, v);
                            should(g.direction(*la));
                            shouldEqual(g.direct(*le, true), *la);
                            shouldEqual(g.direct(*le, false), g.oppositeArc(*la));
                        }
                        else
                        {
                            should(!g.direction(*la));
                            shouldEqual(g.direct(*le, false), *la);
                            shouldEqual(g.direct(*le, true), g.oppositeArc(*la));
                        }
                    }

                    shouldEqual(u, *n);
                    shouldEqual(v, *j);
                    shouldEqual(g.direct(*le, v), *la);
                    shouldEqual(g.direct(*le, u), g.oppositeArc(*la));
                    shouldEqual(g.oppositeNode(u, *le), v);
                    shouldEqual(g.oppositeNode(v, *le), u);
                    shouldEqual(g.oppositeNode(lemon::INVALID, *le), Node(lemon::INVALID));
                    
                    shouldEqual(arcIdMap[*la], g.id(la));
                    shouldEqual(edgeIdMap[*la], g.id((typename Graph::Edge &)*la));
                    shouldEqual(edgeIdMap[*le], g.id(*le));
                    shouldEqual(*la, g.arcFromId(g.id(*la)));
                    shouldEqual(*le, g.edgeFromId(g.id(*le)));
                }
                should(!n.isValid() && n.atEnd());
                should(n == lemon::INVALID);
                should(!(n != lemon::INVALID));
                should(!e.isValid() && e.atEnd());
                should(e == eend);
                should(e == la);
                should(e == le);
                should(e == lemon::INVALID);
                should(!(e != lemon::INVALID));
                should(i == li);
                should(i == iend);
                should(i == lemon::INVALID);
                should(!(i != lemon::INVALID));
                
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
    void testBackNeighborIterator()
    {
        using namespace boost;
        
        static const bool directed = IsSameType<DirectedTag, directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            MultiArray<N, int> vertexMap(s);
            typename Graph::template EdgeMap<int> edgeMap(g);
            typename Graph::template ArcMap<int> arcMap(g);
            
            typename Graph::template ArcMap<int> arcIdMap(g);            
            typename Graph::template EdgeMap<int> edgeIdMap(g);            

            linearSequence(arcIdMap.begin(), arcIdMap.end());
            linearSequence(edgeIdMap.begin(), edgeIdMap.end());
            
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), g.maxUniqueDegree());
            shouldEqual((arcIdMap.shape().template subarray<0, N>()), s);
            shouldEqual(arcIdMap.shape(N), g.maxDegree());
            shouldEqual(edgeMap.shape(), edgeIdMap.shape());
            
            int totalCount = 0;
        
            typename Graph::vertex_iterator j = g.get_vertex_iterator(), 
                                            end = g.get_vertex_end_iterator();
            for(; j != end; ++j)
            {
                typename Graph::back_neighbor_vertex_iterator n = g.get_back_neighbor_vertex_iterator(j), 
                                                              nend = g.get_back_neighbor_vertex_end_iterator(j);
                typename Graph::out_back_edge_iterator        e = g.get_out_back_edge_iterator(j), 
                                                              eend = g.get_out_back_edge_end_iterator(j);
                typename Graph::IncBackEdgeIt                 le(g, j);
                typename Graph::OutBackArcIt                  la(g, j);

                should(n == g.get_back_neighbor_vertex_iterator(*j));
                should(nend == g.get_back_neighbor_vertex_end_iterator(*j));
                should(n == back_adjacent_vertices(*j, g).first);
                should(nend == back_adjacent_vertices(*j, g).second);
                should(e == g.get_out_back_edge_iterator(*j));
                should(eend == g.get_out_back_edge_end_iterator(*j));
                should(e == out_back_edges(*j, g).first);
                should(eend == out_back_edges(*j, g).second);

                int count = 0;
                for(; n != nend; ++n, ++e, ++la, ++le, ++count)
                {
                    should(*n != *j);
                    should(n.isValid() && !n.atEnd());
                    should(n != lemon::INVALID);
                    should(!(n == lemon::INVALID));
                    should(e.isValid() && !e.atEnd());
                    should(e == la);
                    should(e == le);
                    should(e != eend);
                    should(e != lemon::INVALID);
                    should(!(e == lemon::INVALID));
                    
                    shouldEqual(source(*e, g), *j);
                    shouldEqual(target(*e, g), *n);
                    vertexMap[*n] += 1;
                    edgeMap[*e] += 1;
                    arcMap[*la] += 1;
                    
                    typename Graph::Arc a  = g.findArc(*j, *n),
                                        oa = g.findArc(*n, *j);
                    should(a == *la);
                    shouldEqual(source(oa, g), *n);
                    shouldEqual(target(oa, g), *j);
                    shouldEqual(g.baseNode(la), *j);
                    shouldEqual(g.runningNode(la), *n);
                    shouldEqual(g.oppositeArc(a), oa);
                    
                    typename Graph::Edge ge = g.findEdge(*j, *n),
                                         oe = g.findEdge(*n, *j);
                    should(ge == *le);
                    shouldEqual(g.baseNode(le), *j);
                    shouldEqual(g.runningNode(le), *n);
                    
                    typename Graph::Node u = g.u(oe), 
                                         v = g.v(oe);
                    if(!directed)
                    {
                        should(oe == *le);
                        if(g.id(n) < g.id(j))
                            std::swap(u, v);
                    }
                    shouldEqual(u, *n);
                    shouldEqual(v, *j);
                    
                    shouldEqual(arcIdMap[*la], g.id(la));
                    shouldEqual(edgeIdMap[*la], g.id((typename Graph::Edge &)*la));
                    shouldEqual(edgeIdMap[*le], g.id(*le));
                    shouldEqual(*la, g.arcFromId(g.id(*la)));
                    shouldEqual(*le, g.edgeFromId(g.id(*le)));
                }
                should(!n.isValid() && n.atEnd());
                should(n == lemon::INVALID);
                should(!(n != lemon::INVALID));
                should(!e.isValid() && e.atEnd());
                should(e == eend);
                should(e == la);
                should(e == le);
                should(e == lemon::INVALID);
                should(!(e != lemon::INVALID));
                
                shouldEqual(count, g.back_degree(j));
                
                totalCount += count;
            }
            
            // check that all neighbors are found
            if(directed)
                totalCount *= 2;
            shouldEqual(totalCount, g.num_edges());
            
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            edgeMap.minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, g.num_edges() ? 1 : 0);

            typename Graph::edge_propmap_shape_type estart, estop;
            estop.template subarray<0,N>() = s;
            estop[N] = g.maxDegree() / 2;
            shouldEqualSequence(edgeMap.begin(), edgeMap.end(), arcMap.subarray(estart, estop).begin());
            
            estart[N] = g.maxDegree() / 2;
            estop[N] = g.maxDegree();
            min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            arcMap.subarray(estart, estop).minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, 0);

            j = g.get_vertex_iterator();
            for(; j != end; ++j)
            {
                shouldEqual(vertexMap[*j], g.forward_degree(j));
                shouldEqual(edgeMap.bindInner(*j).template sum<int>(), g.back_degree(j));
            }
        }
    }
    
    template <class DirectedTag, NeighborhoodType NType>
    void testEdgeIterator()
    {
        using namespace boost;
        
        static const bool directed = IsSameType<DirectedTag, directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            MultiArray<N, int> sourceVertexMap(s),targetVertexMap(s);
            MultiArray<N+1, int> edgeMap(g.edge_propmap_shape()),
                                 edgeIdMap(g.edge_propmap_shape());
            
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), g.maxUniqueDegree());
            
            linearSequence(edgeIdMap.begin(), edgeIdMap.end());
            
            typename Graph::edge_iterator e = g.get_edge_iterator(),
                                          eend = g.get_edge_end_iterator();
            typename Graph::EdgeIt el(g);

            should(e == edges(g).first);
            should(eend == edges(g).second);
            
            int count = 0;
            int maxEdgeId = -1;
            for(; e != eend; ++e, ++el, ++count)
            {
                should(e.isValid() && !e.atEnd());
                should(e == el);
                should(e != lemon::INVALID);
                should(!(e == lemon::INVALID));

                put(edgeMap, *e, get(edgeMap, *e) + 1); // same as: edgeMap[*e] += 1;
                sourceVertexMap[source(*e, g)] += 1;
                targetVertexMap[target(*e, g)] += 1;
                
                shouldEqual(edgeIdMap[*e], g.id(el));
                shouldEqual(*el, g.edgeFromId(g.id(*el)));

                if(maxEdgeId < g.id(el))
                    maxEdgeId = g.id(el);
            }
            should(!e.isValid() && e.atEnd());
            should(e == el);
            should(e == lemon::INVALID);
            should(!(e != lemon::INVALID));
            
            shouldEqual(maxEdgeId, g.maxEdgeId());
            
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
    
    template <class DirectedTag, NeighborhoodType NType>
    void testArcIterator()
    {
        using namespace boost;
        
        static const bool directed = IsSameType<DirectedTag, directed_tag>::value;
        
        typedef GridGraph<N, DirectedTag> Graph;
        
        MultiCoordinateIterator<N> i(Shape(3)), iend = i.getEndIterator();        
        for(; i != iend; ++i)
        {
            // create all possible array shapes from 1**N to 3**N
            Shape s = *i + Shape(1);
            Graph g(s, NType);
            
            typename Graph::template NodeMap<int> sourceVertexMap(g), targetVertexMap(g);
            typename Graph::template EdgeMap<int> edgeMap(g);
            typename Graph::template ArcMap<int> arcMap(g);
            
            shouldEqual(sourceVertexMap.shape(), s);
            shouldEqual(targetVertexMap.shape(), s);
            shouldEqual((edgeMap.shape().template subarray<0, N>()), s);
            shouldEqual(edgeMap.shape(N), g.maxUniqueDegree());
            shouldEqual((arcMap.shape().template subarray<0, N>()), s);
            shouldEqual(arcMap.shape(N), g.maxDegree());
            
            typename Graph::ArcIt e(g);

            int count = 0;
            int maxArcId = -1;
            for(; e != lemon::INVALID; ++e, ++count)
            {
                should(e.isValid() && !e.atEnd());
                should(e != lemon::INVALID);
                should(!(e == lemon::INVALID));

                sourceVertexMap[source(*e, g)] += 1;
                targetVertexMap[target(*e, g)] += 1;
                edgeMap[*e] += 1;
                arcMap[*e] += 1;
                
                shouldEqual(*e, g.arcFromId(g.id(*e)));
                
                if(maxArcId < g.id(e))
                    maxArcId = g.id(e);
            }
            should(!e.isValid() && e.atEnd());
            should(e == lemon::INVALID);
            should(!(e != lemon::INVALID));
            
            shouldEqual(maxArcId, g.maxArcId());
            
            // check that all neighbors are found
            shouldEqual(count, g.arcNum());
            shouldEqual(count, directed ? g.edgeNum() : 2*g.edgeNum());
            
            int min = NumericTraits<int>::max(), max = NumericTraits<int>::min();
            edgeMap.minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, g.arcNum() ? directed ? 1 : 2 : 0);
            
            min = NumericTraits<int>::max();
            max = NumericTraits<int>::min();
            arcMap.minmax(&min, &max);
            shouldEqual(min, 0);
            shouldEqual(max, g.arcNum() ? 1 : 0);
            
            MultiCoordinateIterator<N> j(s), end = j.getEndIterator();
            for(; j != end; ++j)
            {
                shouldEqual(edgeMap.bindInner(*j).template sum<int>(), directed ? g.out_degree(j) : 2*g.back_degree(j));
                shouldEqual(arcMap.bindInner(*j).template sum<int>(), g.out_degree(j));
                shouldEqual(sourceVertexMap[*j], g.out_degree(j));
                shouldEqual(targetVertexMap[*j], g.out_degree(j));
            }
        }
    }
};

template <unsigned int N>
struct GridGraphAlgorithmTests
{
    typedef typename MultiArrayShape<N>::type Shape;
    
    template <class DirectedTag, NeighborhoodType NType>
    void testLocalMinMax()
    {
        typedef GridGraph<N, DirectedTag> Graph;
        
        Graph g(Shape(3), NType);
        typename Graph::template NodeMap<int> src(g), dest(g);
        
        src[Shape(1)] = 1;
        
        boost_graph::localMinMaxGraph(g, src, dest, 1, -9999, std::greater<int>());
        
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
        add(testCase(&NeighborhoodTests<N>::testVertexIterator));
        
        add(testCase(&NeighborhoodTests<N>::testDirectNeighborhood));
        add(testCase(&NeighborhoodTests<N>::testIndirectNeighborhood));
        
        add(testCase(&NeighborhoodTests<N>::template testNeighborhoodIterator<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testNeighborhoodIterator<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testOutArcIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testOutArcIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testOutArcIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testOutArcIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testArcIteratorDirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testArcIteratorDirected<IndirectNeighborhood>));
        
        add(testCase(&NeighborhoodTests<N>::template testArcIteratorUndirected<DirectNeighborhood>));
        add(testCase(&NeighborhoodTests<N>::template testArcIteratorUndirected<IndirectNeighborhood>));
        
        add(testCase((&GridGraphTests<N>::template testBasics<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBasics<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBasics<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBasics<undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testVertexIterator<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testVertexIterator<undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testNeighborIterator<undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testBackNeighborIterator<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBackNeighborIterator<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBackNeighborIterator<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testBackNeighborIterator<undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testEdgeIterator<undirected_tag, DirectNeighborhood>)));
        
        add(testCase((&GridGraphTests<N>::template testArcIterator<directed_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testArcIterator<undirected_tag, IndirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testArcIterator<directed_tag, DirectNeighborhood>)));
        add(testCase((&GridGraphTests<N>::template testArcIterator<undirected_tag, DirectNeighborhood>)));
        
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
