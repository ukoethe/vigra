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

    void testShortestPath(){

        {
            typedef ShortestPathDijkstra<GraphType,float> Sp;
            typedef Sp::PredecessorsMap PredMap;
            typedef Sp::DistanceMap     DistMap;
            GraphType g(0,0);
            const Node n1=g.addNode(1);
            const Node n2=g.addNode(2);
            const Node n3=g.addNode(3);
            const Node n4=g.addNode(4);
            const Edge e12= g.addEdge(n1,n2);
            const Edge e13= g.addEdge(n1,n3);
            const Edge e24= g.addEdge(n2,n4);
            const Edge e34= g.addEdge(n3,n4);

            //   1 | 2
            //   _   _ 
            //   3 | 4 
        
            GraphType::EdgeMap<float> ew(g);
            ew[e12]=10.0;
            ew[e13]=2.0;
            ew[e24]=3.0;
            ew[e34]=4.0;
            {
                Sp pf(g);
                pf.run(ew,n1,n2);
                const PredMap & pmap = pf.predecessors();
                const DistMap & dmap = pf.distances();

                should(pmap[n2]==n4);
                should(pmap[n4]==n3);
                should(pmap[n3]==n1);

                shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
                shouldEqualTolerance(dmap[n3],2.0f , 0.00001);
                shouldEqualTolerance(dmap[n4],6.0f , 0.00001);
                shouldEqualTolerance(dmap[n2],9.0f , 0.00001);
            }
            {
                Sp pf(g);
                pf.run(ew,n1);
                const PredMap & pmap = pf.predecessors();
                const DistMap & dmap = pf.distances();

                should(pmap[n2]==n4);
                should(pmap[n4]==n3);
                should(pmap[n3]==n1);

                shouldEqualTolerance(dmap[n1],0.0f , 0.00001);
                shouldEqualTolerance(dmap[n3],2.0f , 0.00001);
                shouldEqualTolerance(dmap[n4],6.0f , 0.00001);
                shouldEqualTolerance(dmap[n2],9.0f , 0.00001);
            }
        }

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




};


 
struct GraphAlgorithmTestSuite
: public vigra::test_suite
{
    GraphAlgorithmTestSuite()
    : vigra::test_suite("GraphAlgorithmTestSuite")
    {   
        


        add( testCase( &GraphAlgorithmTest::testShortestPath));
        add( testCase( &GraphAlgorithmTest::testRegionAdjacencyGraph));
        add( testCase( &GraphAlgorithmTest::testEdgeSort));

    }
};

int main(int argc, char ** argv)
{
    GraphAlgorithmTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

