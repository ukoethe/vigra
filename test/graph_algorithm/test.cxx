/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
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
#include "vigra/unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/adjacency_list_graph.hxx"
#include "vigra/graph_algorithms.hxx"
#include "vigra/multi_resize.hxx"

using namespace vigra;

struct GraphAlgorithmTest{


    typedef vigra::AdjacencyListGraph            GraphType;
    typedef GraphType::Node                      Node;
    typedef GraphType::Edge                      Edge;
    typedef GraphType::Arc                       Arc;
    typedef GraphType::EdgeIt                    EdgeIt;
    typedef GraphType::NodeIt                    NodeIt;
    typedef GraphType::ArcIt                     ArcIt;
    typedef GraphType::IncEdgeIt                 IncEdgeIt;
    typedef GraphType::InArcIt                   InArcIt;
    typedef GraphType::OutArcIt                  OutArcIt;
    typedef GraphType::NeighborNodeIt            NeighborNodeIt;

    GraphAlgorithmTest(){       

    }

    template <class Graph>
    void testShortestPathImpl(Graph const & g)
    {
        typedef ShortestPathDijkstra<Graph,float> Sp;
        typedef typename Sp::PredecessorsMap PredMap;
        typedef typename Sp::DistanceMap     DistMap;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;

        //   1 | 2
        //   _   _ 
        //   3 | 4 
        
        typename Graph::NodeIt node(g);
        const Node n1=*node++;
        const Node n2=*node++;
        const Node n3=*node++;
        const Node n4=*node;
        const Edge e12= g.findEdge(n1,n2);
        const Edge e13= g.findEdge(n1,n3);
        const Edge e24= g.findEdge(n2,n4);
        const Edge e34= g.findEdge(n3,n4);

        typename Graph::template EdgeMap<float> ew(g);
        ew[e12]=10.0;
        ew[e13]=2.0;
        ew[e24]=3.0;
        ew[e34]=4.0;
        {
            Sp pf(g);
            pf.run(ew,n1,n2);

            should(pf.source() == n1);
            should(pf.target() == n2);

            const PredMap & pmap = pf.predecessors();
            const DistMap & dmap = pf.distances();

            should(pmap[n2]==n4);
            should(pmap[n4]==n3);
            should(pmap[n3]==n1);
            should(pmap[n1]==n1);

            shouldEqual(pf.discoveryOrder().size(), 4);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n3);
            shouldEqual(pf.discoveryOrder()[2], n4);
            shouldEqual(pf.discoveryOrder()[3], n2);

            shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
            shouldEqualTolerance(dmap[n3],2.0f , 0.00001);
            shouldEqualTolerance(dmap[n4],6.0f , 0.00001);
            shouldEqualTolerance(dmap[n2],9.0f , 0.00001);

            pf.run(ew,n1,n2, 8.0f);

            should(pf.source() == n1);
            should(pf.target() == lemon::INVALID);

            shouldEqual(pf.discoveryOrder().size(), 3);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n3);
            shouldEqual(pf.discoveryOrder()[2], n4);
        }
        {
            Sp pf(g);
            pf.run(ew,n1);

            const PredMap & pmap = pf.predecessors();
            const DistMap & dmap = pf.distances();

            should(pf.source() == n1);
            should(pf.target() == n2);

            should(pmap[n2]==n4);
            should(pmap[n4]==n3);
            should(pmap[n3]==n1);
            should(pmap[n1]==n1);

            shouldEqual(pf.discoveryOrder().size(), 4);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n3);
            shouldEqual(pf.discoveryOrder()[2], n4);
            shouldEqual(pf.discoveryOrder()[3], n2);

            shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
            shouldEqualTolerance(dmap[n3],2.0f , 0.00001);
            shouldEqualTolerance(dmap[n4],6.0f , 0.00001);
            shouldEqualTolerance(dmap[n2],9.0f , 0.00001);

            pf.run(ew,n1, lemon::INVALID, 8.0f);

            should(pf.source() == n1);
            should(pf.target() == n4); // n2 is now unreachable within maxDistance = 8.0

            should(pmap[n2]==lemon::INVALID);
            should(pmap[n4]==n3);
            should(pmap[n3]==n1);
            should(pmap[n1]==n1);

            shouldEqual(pf.discoveryOrder().size(), 3);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n3);
            shouldEqual(pf.discoveryOrder()[2], n4);
        }
    }

    template <class Graph>
    void testShortestPathWithROIImpl(Graph const & g)
    {
        typedef ShortestPathDijkstra<Graph,float> Sp;
        typedef typename Sp::PredecessorsMap PredMap;
        typedef typename Sp::DistanceMap     DistMap;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;

        //   1 | 2
        //   _   _ 
        //   3 | 4 
        
        typename Graph::NodeIt node(g);
        const Node n1=*node++;
        const Node n2=*node++;
        const Node n3=*node++;
        const Node n4=*node;
        const Edge e12= g.findEdge(n1,n2);
        const Edge e13= g.findEdge(n1,n3);
        const Edge e24= g.findEdge(n2,n4);
        const Edge e34= g.findEdge(n3,n4);

        typename Graph::template EdgeMap<float> ew(g);
        ew[e12]=10.0;
        ew[e13]=2.0;
        ew[e24]=3.0;
        ew[e34]=4.0;
        {
            Sp pf(g);

            // ROI = entire graph
            pf.run(Node(0), g.shape(), ew, n1, n2);

            should(pf.source() == n1);
            should(pf.target() == n2);

            const PredMap & pmap = pf.predecessors();
            const DistMap & dmap = pf.distances();

            should(pmap[n2]==n4);
            should(pmap[n4]==n3);
            should(pmap[n3]==n1);
            should(pmap[n1]==n1);

            shouldEqual(pf.discoveryOrder().size(), 4);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n3);
            shouldEqual(pf.discoveryOrder()[2], n4);
            shouldEqual(pf.discoveryOrder()[3], n2);

            shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
            shouldEqualTolerance(dmap[n3],2.0f , 0.00001);
            shouldEqualTolerance(dmap[n4],6.0f , 0.00001);
            shouldEqualTolerance(dmap[n2],9.0f , 0.00001);

            // ROI = top half
            pf.run(n1, n2+Node(1), ew, n1, n2);

            should(pf.source() == n1);
            should(pf.target() == n2);

            should(pmap[n2]==n1);
            should(pmap[n1]==n1);

            shouldEqual(pf.discoveryOrder().size(), 2);
            shouldEqual(pf.discoveryOrder()[0], n1);
            shouldEqual(pf.discoveryOrder()[1], n2);

            shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
            shouldEqualTolerance(dmap[n2],10.0f , 0.00001);

            // ROI = top half, maxWeight less then weight of e12
            pf.run(n1, n2+Node(1), ew, n1, n2, 8.0f);

            should(pf.source() == n1);
            should(pf.target() == lemon::INVALID);

            shouldEqual(pf.discoveryOrder().size(), 1);
            shouldEqual(pf.discoveryOrder()[0], n1);
        }
        {
            Sp pf(g);
            // ROI = bottom half
            pf.run(n3, g.shape(), ew, n4, n3);

            should(pf.source() == n4);
            should(pf.target() == n3);

            const PredMap & pmap = pf.predecessors();
            const DistMap & dmap = pf.distances();

            should(pmap[n3]==n4);
            should(pmap[n4]==n4);

            shouldEqual(pf.discoveryOrder().size(), 2);
            shouldEqual(pf.discoveryOrder()[0], n4);
            shouldEqual(pf.discoveryOrder()[1], n3);

            shouldEqualTolerance(dmap[n4],0.0f , 0.00001);
            shouldEqualTolerance(dmap[n3],4.0f , 0.00001);

            // ROI = bottom half, maxWeight less then weight of e34
            pf.run(n3, g.shape(), ew, n4, n3, 2.0);

            should(pf.source() == n4);
            should(pf.target() == lemon::INVALID);

            shouldEqual(pf.discoveryOrder().size(), 1);
            shouldEqual(pf.discoveryOrder()[0], n4);
        }
    }

    void testShortestPathAdjacencyListGraph()
    {
        GraphType g(0,0);
        const Node n1=g.addNode(1);
        const Node n2=g.addNode(2);
        const Node n3=g.addNode(3);
        const Node n4=g.addNode(4);
        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);

        testShortestPathImpl(g);
    }

    void testShortestPathGridGraph()
    {
        GridGraph<2> g(Shape2(2,2), DirectNeighborhood);

        testShortestPathImpl(g);
        testShortestPathWithROIImpl(g);
    }

    void testRegionAdjacencyGraph(){
        {
            GraphType g(0,0);
            const Node n1=g.addNode(1);
            const Node n2=g.addNode(2);
            const Node n3=g.addNode(3);
            const Node n4=g.addNode(4);
            const Node n5=g.addNode(5);
            const Edge e12= g.addEdge(n1,n2);
            g.addEdge(n1,n3);
            g.addEdge(n2,n4);
            const Edge e34= g.addEdge(n3,n4);
            const Edge e45= g.addEdge(n4,n5);

            //   1 | 2
            //   _   _ 
            //   3 | 4 | 5 

            //labeling
            //   1 | 7
            //   _   _ 
            //   1 | 7 | 3 

            GraphType::NodeMap<int> labels(g);
            labels[n1]=1;
            labels[n2]=7;
            labels[n3]=1;
            labels[n4]=7;
            labels[n5]=3;

            GraphType rag;
            GraphType::EdgeMap< std::vector<Edge> > affEdges;

            makeRegionAdjacencyGraph(g,labels,rag,affEdges);

            shouldEqual(rag.nodeNum(),3);
            shouldEqual(rag.edgeNum(),2);

            const Node rn1 = rag.nodeFromId(1);
            const Node rn7 = rag.nodeFromId(7);
            const Node rn3 = rag.nodeFromId(3);

            should(rag.nodeFromId(0)==lemon::INVALID);
            should(rag.nodeFromId(2)==lemon::INVALID);
            should(rag.nodeFromId(4)==lemon::INVALID);
            should(rag.nodeFromId(5)==lemon::INVALID);
            should(rag.nodeFromId(6)==lemon::INVALID);
            should(rag.nodeFromId(8)==lemon::INVALID);

            should(rn1!=lemon::INVALID);
            should(rn7!=lemon::INVALID);
            should(rn3!=lemon::INVALID);

            const Edge re17 = rag.findEdge(rn1,rn7);
            const Edge re73 = rag.findEdge(rn7,rn3);

            should(re17!=lemon::INVALID);
            should(re73!=lemon::INVALID);

            should(rag.findEdge(rn1,rn3)==lemon::INVALID);


            shouldEqual(affEdges[re17].size(),2);
            shouldEqual(affEdges[re73].size(),1);

            should(affEdges[re73][0]==e45);
            should(affEdges[re17][0]==e12 || affEdges[re17][1]==e12 );
            should(affEdges[re17][0]==e34 || affEdges[re17][1]==e34 );
        }   
    }


    void testEdgeSort(){
        {
            GraphType g(0,0);
            const Node n1=g.addNode(1);
            const Node n2=g.addNode(2);
            const Node n3=g.addNode(3);
            const Node n4=g.addNode(4);
            const Node n5=g.addNode(5);
            const Edge e12= g.addEdge(n1,n2);
            const Edge e13= g.addEdge(n1,n3);
            const Edge e24= g.addEdge(n2,n4);
            const Edge e34= g.addEdge(n3,n4);
            const Edge e45= g.addEdge(n4,n5);

            //   1 | 2
            //   _   _ 
            //   3 | 4 | 5 

            GraphType::EdgeMap<float> ew(g);

            ew[e12]=2.0;
            ew[e13]=1.0;
            ew[e24]=5.0;
            ew[e34]=4.0;
            ew[e45]=3.0;


            std::vector<Edge> edgeVec;

            std::less<float> l;
            edgeSort(g,ew,l,edgeVec);

            shouldEqual(edgeVec.size(),g.edgeNum());
            should(edgeVec[0]==e13);
            should(edgeVec[1]==e12);
            should(edgeVec[2]==e45);
            should(edgeVec[3]==e34);
            should(edgeVec[4]==e24);


            std::greater<float> gr;
            edgeSort(g,ew,gr,edgeVec);

            shouldEqual(edgeVec.size(),g.edgeNum());
            should(edgeVec[4]==e13);
            should(edgeVec[3]==e12);
            should(edgeVec[2]==e45);
            should(edgeVec[1]==e34);
            should(edgeVec[0]==e24);

        }
    }

    void testEdgeWeightComputation()
    {
        MultiArray<2, double> nodeMap(Shape2(3,2), LinearSequence);
        MultiArray<2, double> interpolated(Shape2(5,3));
        resizeImageLinearInterpolation(nodeMap, interpolated);

        GridGraph<2> g(nodeMap.shape(), IndirectNeighborhood);
        GridGraph<2>::EdgeMap<double> edgeMap1(g), edgeMap2(g);

        edgeWeightsFromNodeWeights(g, nodeMap, edgeMap1);
        edgeWeightsFromInterpolatedImage(g, interpolated, edgeMap2);

        shouldEqual(edgeMap1.shape(), Shape3(3,2,4));

        double ref[] = {
            0.0, 0.0, 0.0, 0.0, 2.0, 3.0,
            0.0, 0.0, 0.0, 1.5, 2.5, 3.5,
            0.0, 0.0, 0.0, 2.0, 3.0, 0.0,
            0.0, 0.5, 1.5, 0.0, 3.5, 4.5
        };

        shouldEqualSequence(edgeMap1.begin(), edgeMap1.end(), ref);
        shouldEqualSequence(edgeMap2.begin(), edgeMap2.end(), ref);

        edgeWeightsFromNodeWeights(g, nodeMap, edgeMap1, true);
        edgeWeightsFromInterpolatedImage(g, interpolated, edgeMap2, true);

        double ref2[] = {
            0.0, 0.0, 0.0, 0.0, 2.0*M_SQRT2, 3.0*M_SQRT2,
            0.0, 0.0, 0.0, 1.5, 2.5, 3.5,
            0.0, 0.0, 0.0, 2.0*M_SQRT2, 3.0*M_SQRT2, 0.0,
            0.0, 0.5, 1.5, 0.0, 3.5, 4.5
        };

        shouldEqualSequence(edgeMap1.begin(), edgeMap1.end(), ref2);
        shouldEqualSequence(edgeMap2.begin(), edgeMap2.end(), ref2);
    }
};


 
struct GraphAlgorithmTestSuite
: public vigra::test_suite
{
    GraphAlgorithmTestSuite()
    : vigra::test_suite("GraphAlgorithmTestSuite")
    {   
        add( testCase( &GraphAlgorithmTest::testShortestPathAdjacencyListGraph));
        add( testCase( &GraphAlgorithmTest::testShortestPathGridGraph));
        add( testCase( &GraphAlgorithmTest::testRegionAdjacencyGraph));
        add( testCase( &GraphAlgorithmTest::testEdgeSort));
        add( testCase( &GraphAlgorithmTest::testEdgeWeightComputation));
    }
};

int main(int argc, char ** argv)
{
    GraphAlgorithmTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

