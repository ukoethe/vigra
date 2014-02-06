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

using namespace vigra;






struct AdjacencyListGraphOffTest{


    typedef vigra::AdjacencyListGraph                     GraphType;
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
    AdjacencyListGraphOffTest(){       

    }

    void adjGraphSimpleTestStart0(){
        GraphType g(0,0,true);

        // add nodes
        shouldEqual(g.nodeNum(),0);
        shouldEqual(g.edgeNum(),0);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)==lemon::INVALID);
        should(g.nodeFromId(2)==lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n1 = g.addNode();
        shouldEqual(g.nodeNum(),1);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),0);
        should(g.nodeFromId(0)!=lemon::INVALID);
        should(g.nodeFromId(1)==lemon::INVALID);
        should(g.nodeFromId(2)==lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n2 = g.addNode();
        shouldEqual(g.nodeNum(),2);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),1);
        should(g.nodeFromId(0)!=lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)==lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n3 = g.addNode();
        shouldEqual(g.nodeNum(),3);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),2);
        should(g.nodeFromId(0)!=lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);


        Node n4 = g.addNode();
        shouldEqual(g.nodeNum(),4);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),3);
        should(g.nodeFromId(0)!=lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)!=lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        shouldEqual(g.edgeNum(),0);

        should(  g.findEdge(n1,n2) == lemon::INVALID  );
        should(  g.findEdge(n1,n3) == lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        // add edges
        // SET UP THIS GRAPH 
        // 1 |3
        // __ __
        // 2 |4

        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)==lemon::INVALID);
        should(g.edgeFromId(2)==lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);


        Edge e12  = g.addEdge(n1,n2);
        shouldEqual(g.edgeNum(),1);
        shouldEqual(g.maxEdgeId(),0);
        should(g.u(e12)==n1);
        should(g.v(e12)==n2);
        should(g.edgeFromId(0)!=lemon::INVALID);
        should(g.edgeFromId(1)==lemon::INVALID);
        should(g.edgeFromId(2)==lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) == lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e13  = g.addEdge(n1,n3);
        shouldEqual(g.edgeNum(),2);
        shouldEqual(g.maxEdgeId(),1);
        should(g.u(e13)==n1);
        should(g.v(e13)==n3);
        should(g.edgeFromId(0)!=lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)==lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e24  = g.addEdge(n2,n4);
        shouldEqual(g.edgeNum(),3);  
        shouldEqual(g.maxEdgeId(),2);
        should(g.u(e24)==n2);
        should(g.v(e24)==n4);
        should(g.edgeFromId(0)!=lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) != lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e34  = g.addEdge(n3,n4);
        shouldEqual(g.edgeNum(),4);
        shouldEqual(g.maxEdgeId(),3);
        should(g.u(e34)==n3);
        should(g.v(e34)==n4);
        should(g.edgeFromId(0)!=lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)!=lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) != lemon::INVALID  );
        should(  g.findEdge(n3,n4) != lemon::INVALID  );

        // WE HAVE THIS THIS GRAPH 
        // 1 |3
        // __ __
        // 2 |4
    }
    void adjGraphSimpleTestStart1(){
        GraphType g(0,0,false);

        // add nodes
        shouldEqual(g.nodeNum(),0);
        shouldEqual(g.edgeNum(),0);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)==lemon::INVALID);
        should(g.nodeFromId(2)==lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n1 = g.addNode();
        shouldEqual(g.nodeNum(),1);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),1);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)==lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n2 = g.addNode();
        shouldEqual(g.nodeNum(),2);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),2);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)==lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        Node n3 = g.addNode();
        shouldEqual(g.nodeNum(),3);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),3);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)!=lemon::INVALID);
        should(g.nodeFromId(4)==lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);


        Node n4 = g.addNode();
        shouldEqual(g.nodeNum(),4);
        shouldEqual(g.edgeNum(),0);
        shouldEqual(g.maxNodeId(),4);
        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)!=lemon::INVALID);
        should(g.nodeFromId(4)!=lemon::INVALID);
        should(g.nodeFromId(5)==lemon::INVALID);

        shouldEqual(g.edgeNum(),0);

        should(  g.findEdge(n1,n2) == lemon::INVALID  );
        should(  g.findEdge(n1,n3) == lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        // add edges
        // SET UP THIS GRAPH 
        // 1 |3
        // __ __
        // 2 |4

        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)==lemon::INVALID);
        should(g.edgeFromId(2)==lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);


        Edge e12  = g.addEdge(n1,n2);
        shouldEqual(g.edgeNum(),1);
        shouldEqual(g.maxEdgeId(),1);
        should(g.u(e12)==n1);
        should(g.v(e12)==n2);
        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)==lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) == lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e13  = g.addEdge(n1,n3);
        shouldEqual(g.edgeNum(),2);
        shouldEqual(g.maxEdgeId(),2);
        should(g.u(e13)==n1);
        should(g.v(e13)==n3);
        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)==lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) == lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e24  = g.addEdge(n2,n4);
        shouldEqual(g.edgeNum(),3);  
        shouldEqual(g.maxEdgeId(),3);
        should(g.u(e24)==n2);
        should(g.v(e24)==n4);
        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)!=lemon::INVALID);
        should(g.edgeFromId(4)==lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) != lemon::INVALID  );
        should(  g.findEdge(n3,n4) == lemon::INVALID  );

        Edge e34  = g.addEdge(n3,n4);
        shouldEqual(g.edgeNum(),4);
        shouldEqual(g.maxEdgeId(),4);
        should(g.u(e34)==n3);
        should(g.v(e34)==n4);
        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)!=lemon::INVALID);
        should(g.edgeFromId(4)!=lemon::INVALID);
        should(g.edgeFromId(5)==lemon::INVALID);
        should(  g.findEdge(n1,n2) != lemon::INVALID  );
        should(  g.findEdge(n1,n3) != lemon::INVALID  );
        should(  g.findEdge(n1,n4) == lemon::INVALID  );
        should(  g.findEdge(n2,n3) == lemon::INVALID  );
        should(  g.findEdge(n2,n4) != lemon::INVALID  );
        should(  g.findEdge(n3,n4) != lemon::INVALID  );

        // WE HAVE THIS THIS GRAPH 
        // 1 |3
        // __ __
        // 2 |4
    }

    void adjGraphEdgeItTestStart0()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g(0,0,true);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);


        

        {
            EdgeIt begin(g);
            EdgeIt invalid(lemon::INVALID);

            should(begin!=lemon::INVALID);

            
            std::vector<Edge> edgeVec(begin,invalid);
            shouldEqual(4,edgeVec.size());
            shouldEqual(0,g.id(edgeVec[0]));
            shouldEqual(1,g.id(edgeVec[1]));
            shouldEqual(2,g.id(edgeVec[2]));
            shouldEqual(3,g.id(edgeVec[3]));
        }
        {
            EdgeIt begin(g);
            should(begin!=lemon::INVALID);

            EdgeIt empty;
            std::vector<Edge> edgeVec(begin,empty);
            shouldEqual(4,edgeVec.size());
            shouldEqual(0,g.id(edgeVec[0]));
            shouldEqual(1,g.id(edgeVec[1]));
            shouldEqual(2,g.id(edgeVec[2]));
            shouldEqual(3,g.id(edgeVec[3]));
        }
        {
            EdgeIt begin(g,g.edgeFromId(1));
            should(begin!=lemon::INVALID);

            EdgeIt empty;
            std::vector<Edge> edgeVec(begin,empty);
            shouldEqual(3,edgeVec.size());
            shouldEqual(1,g.id(edgeVec[0]));
            shouldEqual(2,g.id(edgeVec[1]));
            shouldEqual(3,g.id(edgeVec[2]));
        }
        {
            EdgeIt begin(g,g.edgeFromId(1));
            EdgeIt end(g,g.edgeFromId(2));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);    
            should(begin!=end);

            shouldEqual(std::distance(begin,end),1);
            std::vector<Edge> edgeVec(begin,end);
            shouldEqual(1,edgeVec.size());
            shouldEqual(1,g.id(edgeVec[0]));
        }

        {
            EdgeIt begin(g,g.edgeFromId(1));
            EdgeIt end(g,g.edgeFromId(3));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);

            std::vector<Edge> edgeVec(begin,end);
            shouldEqual(2,edgeVec.size());
            shouldEqual(1,g.id(edgeVec[0]));
            shouldEqual(2,g.id(edgeVec[1]));
        }
    }
    void adjGraphEdgeItTestStart1()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);


        

        {
            EdgeIt begin(g);
            EdgeIt invalid(lemon::INVALID);

            should(begin!=lemon::INVALID);

            
            std::vector<Edge> edgeVec(begin,invalid);
            shouldEqual(4,edgeVec.size());
            shouldEqual(1,g.id(edgeVec[0]));
            shouldEqual(2,g.id(edgeVec[1]));
            shouldEqual(3,g.id(edgeVec[2]));
            shouldEqual(4,g.id(edgeVec[3]));
        }
        {
            EdgeIt begin(g);
            should(begin!=lemon::INVALID);

            EdgeIt empty;
            std::vector<Edge> edgeVec(begin,empty);
            shouldEqual(4,edgeVec.size());
            shouldEqual(1,g.id(edgeVec[0]));
            shouldEqual(2,g.id(edgeVec[1]));
            shouldEqual(3,g.id(edgeVec[2]));
            shouldEqual(4,g.id(edgeVec[3]));
        }
        {
            EdgeIt begin(g,g.edgeFromId(2));
            should(begin!=lemon::INVALID);

            EdgeIt empty;
            std::vector<Edge> edgeVec(begin,empty);
            shouldEqual(3,edgeVec.size());
            shouldEqual(2,g.id(edgeVec[0]));
            shouldEqual(3,g.id(edgeVec[1]));
            shouldEqual(4,g.id(edgeVec[2]));
        }
        {
            EdgeIt begin(g,g.edgeFromId(2));
            EdgeIt end(g,g.edgeFromId(3));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);    
            should(begin!=end);

            shouldEqual(std::distance(begin,end),1);
            std::vector<Edge> edgeVec(begin,end);
            shouldEqual(1,edgeVec.size());
            shouldEqual(2,g.id(edgeVec[0]));
        }

        {
            EdgeIt begin(g,g.edgeFromId(2));
            EdgeIt end(g,g.edgeFromId(4));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);

            std::vector<Edge> edgeVec(begin,end);
            shouldEqual(2,edgeVec.size());
            shouldEqual(2,g.id(edgeVec[0]));
            shouldEqual(3,g.id(edgeVec[1]));
        }
    }


    void adjGraphNodeItTestStart0()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g(0,0,true);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);

        {
            NodeIt begin(g);
            NodeIt invalid(lemon::INVALID);

            should(begin!=lemon::INVALID);            
            std::vector<Node> nodeVec(begin,invalid);
            shouldEqual(4,nodeVec.size());
            shouldEqual(0,g.id(nodeVec[0]));
            shouldEqual(1,g.id(nodeVec[1]));
            shouldEqual(2,g.id(nodeVec[2]));
            shouldEqual(3,g.id(nodeVec[3]));
        }
        {
            NodeIt begin(g);
            should(begin!=lemon::INVALID);

            NodeIt empty;
            std::vector<Node> nodeVec(begin,empty);
            shouldEqual(4,nodeVec.size());
            shouldEqual(0,g.id(nodeVec[0]));
            shouldEqual(1,g.id(nodeVec[1]));
            shouldEqual(2,g.id(nodeVec[2]));
            shouldEqual(3,g.id(nodeVec[3]));
        }
        {
            NodeIt begin(g,g.nodeFromId(1));
            should(begin!=lemon::INVALID);

            NodeIt empty;
            std::vector<Node> nodeVec(begin,empty);
            shouldEqual(3,nodeVec.size());
            shouldEqual(1,g.id(nodeVec[0]));
            shouldEqual(2,g.id(nodeVec[1]));
            shouldEqual(3,g.id(nodeVec[2]));
        }
        {
            NodeIt begin(g,g.nodeFromId(1));
            NodeIt end(g,g.nodeFromId(2));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);    
            should(begin!=end);

            shouldEqual(std::distance(begin,end),1);
            std::vector<Node> nodeVec(begin,end);
            shouldEqual(1,nodeVec.size());
            shouldEqual(1,g.id(nodeVec[0]));
        }

        {
            NodeIt begin(g,g.nodeFromId(1));
            NodeIt end(g,g.nodeFromId(3));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);

            std::vector<Node> nodeVec(begin,end);
            shouldEqual(2,nodeVec.size());
            shouldEqual(1,g.id(nodeVec[0]));
            shouldEqual(2,g.id(nodeVec[1]));
        }
    }
    void adjGraphNodeItTestStart1()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);

        {
            NodeIt begin(g);
            NodeIt invalid(lemon::INVALID);

            should(begin!=lemon::INVALID);            
            std::vector<Node> nodeVec(begin,invalid);
            shouldEqual(4,nodeVec.size());
            shouldEqual(1,g.id(nodeVec[0]));
            shouldEqual(2,g.id(nodeVec[1]));
            shouldEqual(3,g.id(nodeVec[2]));
            shouldEqual(4,g.id(nodeVec[3]));
        }
        {
            NodeIt begin(g);
            should(begin!=lemon::INVALID);

            NodeIt empty;
            std::vector<Node> nodeVec(begin,empty);
            shouldEqual(4,nodeVec.size());
            shouldEqual(1,g.id(nodeVec[0]));
            shouldEqual(2,g.id(nodeVec[1]));
            shouldEqual(3,g.id(nodeVec[2]));
            shouldEqual(4,g.id(nodeVec[3]));
        }
        {
            NodeIt begin(g,g.nodeFromId(2));
            should(begin!=lemon::INVALID);

            NodeIt empty;
            std::vector<Node> nodeVec(begin,empty);
            shouldEqual(3,nodeVec.size());
            shouldEqual(2,g.id(nodeVec[0]));
            shouldEqual(3,g.id(nodeVec[1]));
            shouldEqual(4,g.id(nodeVec[2]));
        }
        {
            NodeIt begin(g,g.nodeFromId(2));
            NodeIt end(g,g.nodeFromId(3));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);    
            should(begin!=end);

            shouldEqual(std::distance(begin,end),1);
            std::vector<Node> nodeVec(begin,end);
            shouldEqual(1,nodeVec.size());
            shouldEqual(2,g.id(nodeVec[0]));
        }

        {
            NodeIt begin(g,g.nodeFromId(2));
            NodeIt end(g,g.nodeFromId(4));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);

            std::vector<Node> nodeVec(begin,end);
            shouldEqual(2,nodeVec.size());
            shouldEqual(2,g.id(nodeVec[0]));
            shouldEqual(3,g.id(nodeVec[1]));
        }
    }
    void adjGraphTestStart0()
    {

        GraphType g(0,0,true);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);




        
        // assert basic sizes
        should(g.edgeNum()==4);
        should(g.nodeNum()==4);
        should(g.arcNum()==8);
        should(g.maxEdgeId()==3);
        should(g.maxNodeId()==3);
        should(g.maxArcId()==7);



        NodeIt nbegin(g);
        NodeIt nend(lemon::INVALID);

        EdgeIt ebegin(g);
        EdgeIt eend(lemon::INVALID);


        ArcIt abegin(g);
        ArcIt aend(lemon::INVALID);

        std::vector<Node> allNodes(nbegin,nend);
        should(allNodes.size()==4);

        std::vector<Edge> allEdges(ebegin,eend);
        should(allEdges.size()==4);


        should(0==g.id( allEdges[0] ) );
        should(1==g.id( allEdges[1] ) );
        should(2==g.id( allEdges[2] ) );
        should(3==g.id( allEdges[3] ) );
        

        std::vector<Arc>  allArcs( abegin,aend);

        should(allArcs.size() ==8);


        shouldEqual(g.id(allArcs[0]),0);
        shouldEqual(g.id(allArcs[1]),1);
        shouldEqual(g.id(allArcs[2]),2);
        shouldEqual(g.id(allArcs[3]),3);
        shouldEqual(g.id(allArcs[4]),4);
        shouldEqual(g.id(allArcs[5]),5);
        shouldEqual(g.id(allArcs[6]),6);
        shouldEqual(g.id(allArcs[7]),7);

        shouldEqual(allArcs[0].edgeId(),0);
        shouldEqual(allArcs[1].edgeId(),1);
        shouldEqual(allArcs[2].edgeId(),2);
        shouldEqual(allArcs[3].edgeId(),3);
        shouldEqual(allArcs[4].edgeId(),0);
        shouldEqual(allArcs[5].edgeId(),1);
        shouldEqual(allArcs[6].edgeId(),2);
        shouldEqual(allArcs[7].edgeId(),3);
        
        shouldEqual(allArcs[0].id(),0);
        shouldEqual(allArcs[1].id(),1);
        shouldEqual(allArcs[2].id(),2);
        shouldEqual(allArcs[3].id(),3);
        shouldEqual(allArcs[4].id(),4);
        shouldEqual(allArcs[5].id(),5);
        shouldEqual(allArcs[6].id(),6);
        shouldEqual(allArcs[7].id(),7);


        shouldEqual( g.id(g.source(allArcs[0])),g.id(g.u(allEdges[0])));
        shouldEqual( g.id(g.source(allArcs[1])),g.id(g.u(allEdges[1])));
        shouldEqual( g.id(g.source(allArcs[2])),g.id(g.u(allEdges[2])));
        shouldEqual( g.id(g.source(allArcs[3])),g.id(g.u(allEdges[3])));
        shouldEqual( g.id(g.target(allArcs[0])),g.id(g.v(allEdges[0])));
        shouldEqual( g.id(g.target(allArcs[1])),g.id(g.v(allEdges[1])));
        shouldEqual( g.id(g.target(allArcs[2])),g.id(g.v(allEdges[2])));
        shouldEqual( g.id(g.target(allArcs[3])),g.id(g.v(allEdges[3])));

        shouldEqual( g.id(g.source(allArcs[0])),g.id(g.target(allArcs[0+4])));
        shouldEqual( g.id(g.source(allArcs[1])),g.id(g.target(allArcs[1+4])));
        shouldEqual( g.id(g.source(allArcs[2])),g.id(g.target(allArcs[2+4])));
        shouldEqual( g.id(g.source(allArcs[3])),g.id(g.target(allArcs[3+4])));
        shouldEqual( g.id(g.target(allArcs[0])),g.id(g.source(allArcs[0+4])));
        shouldEqual( g.id(g.target(allArcs[1])),g.id(g.source(allArcs[1+4])));
        shouldEqual( g.id(g.target(allArcs[2])),g.id(g.source(allArcs[2+4])));
        shouldEqual( g.id(g.target(allArcs[3])),g.id(g.source(allArcs[3+4])));

        shouldEqual( g.id(g.source(allArcs[4])),g.id(g.v(allEdges[0])));
        shouldEqual( g.id(g.source(allArcs[5])),g.id(g.v(allEdges[1])));
        shouldEqual( g.id(g.source(allArcs[6])),g.id(g.v(allEdges[2])));
        shouldEqual( g.id(g.source(allArcs[7])),g.id(g.v(allEdges[3])));
        shouldEqual( g.id(g.target(allArcs[4])),g.id(g.u(allEdges[0])));
        shouldEqual( g.id(g.target(allArcs[5])),g.id(g.u(allEdges[1])));
        shouldEqual( g.id(g.target(allArcs[6])),g.id(g.u(allEdges[2])));
        shouldEqual( g.id(g.target(allArcs[7])),g.id(g.u(allEdges[3])));
    }
    void adjGraphTestStart1()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

       g.addEdge(n1,n2);
       g.addEdge(n1,n3);
       g.addEdge(n2,n4);
       g.addEdge(n3,n4);




        
        // assert basic sizes
        should(g.edgeNum()==4);
        should(g.nodeNum()==4);
        should(g.arcNum()==8);
        should(g.maxEdgeId()==4);
        should(g.maxNodeId()==4);
        should(g.maxArcId()==9);

        should( g.findEdge(g.nodeFromId(1),g.nodeFromId(3) )!=lemon::INVALID);
        should( g.findEdge(g.nodeFromId(1),g.nodeFromId(2) )!=lemon::INVALID);
        should( g.findEdge(g.nodeFromId(2),g.nodeFromId(4) )!=lemon::INVALID);
        should( g.findEdge(g.nodeFromId(3),g.nodeFromId(4) )!=lemon::INVALID);
        should( g.findEdge(g.nodeFromId(1),g.nodeFromId(4) )==lemon::INVALID);
        should( g.findEdge(g.nodeFromId(2),g.nodeFromId(3) )==lemon::INVALID);

        NodeIt nbegin(g);
        NodeIt nend(lemon::INVALID);

        EdgeIt ebegin(g);
        EdgeIt eend(lemon::INVALID);


        ArcIt abegin(g);
        ArcIt aend(lemon::INVALID);

        std::vector<Node> allNodes(nbegin,nend);
        should(allNodes.size()==4);

        std::vector<Edge> allEdges(ebegin,eend);
        should(allEdges.size()==4);


        should(1==g.id( allEdges[0] ) );
        should(2==g.id( allEdges[1] ) );
        should(3==g.id( allEdges[2] ) );
        should(4==g.id( allEdges[3] ) );
        

        std::vector<Arc>  allArcs( abegin,aend);
        should(allArcs.size() ==8);

        shouldEqual(g.id(allArcs[0]),1);
        shouldEqual(g.id(allArcs[1]),2);
        shouldEqual(g.id(allArcs[2]),3);
        shouldEqual(g.id(allArcs[3]),4);
        shouldEqual(g.id(allArcs[4]),6);
        shouldEqual(g.id(allArcs[5]),7);
        shouldEqual(g.id(allArcs[6]),8);
        shouldEqual(g.id(allArcs[7]),9);

        shouldEqual(allArcs[0].edgeId(),1);
        shouldEqual(allArcs[1].edgeId(),2);
        shouldEqual(allArcs[2].edgeId(),3);
        shouldEqual(allArcs[3].edgeId(),4);
        shouldEqual(allArcs[4].edgeId(),1);
        shouldEqual(allArcs[5].edgeId(),2);
        shouldEqual(allArcs[6].edgeId(),3);
        shouldEqual(allArcs[7].edgeId(),4);
        
        shouldEqual(allArcs[0].id(),1);
        shouldEqual(allArcs[1].id(),2);
        shouldEqual(allArcs[2].id(),3);
        shouldEqual(allArcs[3].id(),4);
        shouldEqual(allArcs[4].id(),6);
        shouldEqual(allArcs[5].id(),7);
        shouldEqual(allArcs[6].id(),8);
        shouldEqual(allArcs[7].id(),9);
        

        should( g.id(g.source(allArcs[0]))==g.id(g.u(allEdges[0])));
        should( g.id(g.source(allArcs[1]))==g.id(g.u(allEdges[1])));
        should( g.id(g.source(allArcs[2]))==g.id(g.u(allEdges[2])));
        should( g.id(g.source(allArcs[3]))==g.id(g.u(allEdges[3])));

        should( g.id(g.target(allArcs[0]))==g.id(g.v(allEdges[0])));
        should( g.id(g.target(allArcs[1]))==g.id(g.v(allEdges[1])));
        should( g.id(g.target(allArcs[2]))==g.id(g.v(allEdges[2])));
        should( g.id(g.target(allArcs[3]))==g.id(g.v(allEdges[3])));


        should( g.id(g.source(allArcs[0]))==g.id(g.target(allArcs[0+4])));
        should( g.id(g.source(allArcs[1]))==g.id(g.target(allArcs[1+4])));
        should( g.id(g.source(allArcs[2]))==g.id(g.target(allArcs[2+4])));
        should( g.id(g.source(allArcs[3]))==g.id(g.target(allArcs[3+4])));

        should( g.id(g.target(allArcs[0]))==g.id(g.source(allArcs[0+4])));
        should( g.id(g.target(allArcs[1]))==g.id(g.source(allArcs[1+4])));
        should( g.id(g.target(allArcs[2]))==g.id(g.source(allArcs[2+4])));
        should( g.id(g.target(allArcs[3]))==g.id(g.source(allArcs[3+4])));


        should( g.id(g.source(allArcs[4]))==g.id(g.v(allEdges[0])));
        should( g.id(g.source(allArcs[5]))==g.id(g.v(allEdges[1])));
        should( g.id(g.source(allArcs[6]))==g.id(g.v(allEdges[2])));
        should( g.id(g.source(allArcs[7]))==g.id(g.v(allEdges[3])));

        should( g.id(g.target(allArcs[4]))==g.id(g.u(allEdges[0])));
        should( g.id(g.target(allArcs[5]))==g.id(g.u(allEdges[1])));
        should( g.id(g.target(allArcs[6]))==g.id(g.u(allEdges[2])));
        should( g.id(g.target(allArcs[7]))==g.id(g.u(allEdges[3])));
    }

    void adjGraphArcTest()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);
        
        // check sources and targets of arcs which are just the "natural edges"
        should(  g.source(g.arcFromId(1)) == g.u(g.edgeFromId(1)) );
        should(  g.source(g.arcFromId(2)) == g.u(g.edgeFromId(2)) );
        should(  g.source(g.arcFromId(3)) == g.u(g.edgeFromId(3)) );
        should(  g.source(g.arcFromId(4)) == g.u(g.edgeFromId(4)) );

        should(  g.target(g.arcFromId(1)) == g.v(g.edgeFromId(1)) );
        should(  g.target(g.arcFromId(2)) == g.v(g.edgeFromId(2)) );
        should(  g.target(g.arcFromId(3)) == g.v(g.edgeFromId(3)) );
        should(  g.target(g.arcFromId(4)) == g.v(g.edgeFromId(4)) );




        // check sources and targets of arcs which are flipped "natural edges"
        should(  g.source(g.arcFromId(6)) == g.v(g.edgeFromId(1)) );
        should(  g.source(g.arcFromId(7)) == g.v(g.edgeFromId(2)) );
        should(  g.source(g.arcFromId(8)) == g.v(g.edgeFromId(3)) );
        should(  g.source(g.arcFromId(9)) == g.v(g.edgeFromId(4)) );
        should(  g.target(g.arcFromId(6)) == g.u(g.edgeFromId(1)) );
        should(  g.target(g.arcFromId(7)) == g.u(g.edgeFromId(2)) );
        should(  g.target(g.arcFromId(8)) == g.u(g.edgeFromId(3)) );
        should(  g.target(g.arcFromId(9)) == g.u(g.edgeFromId(4)) );

        // check that arcs are convertible to edges
        should(Edge(g.arcFromId(1))==g.edgeFromId(1));
        should(Edge(g.arcFromId(2))==g.edgeFromId(2));
        should(Edge(g.arcFromId(3))==g.edgeFromId(3));
        should(Edge(g.arcFromId(4))==g.edgeFromId(4));

        should(Edge(g.arcFromId(6))==g.edgeFromId(1));
        should(Edge(g.arcFromId(7))==g.edgeFromId(2));
        should(Edge(g.arcFromId(8))==g.edgeFromId(3));
        should(Edge(g.arcFromId(9))==g.edgeFromId(4));

    }

    void adjGraphArcItTest()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        g.addEdge(n1,n2);
        g.addEdge(n1,n3);
        g.addEdge(n2,n4);
        g.addEdge(n3,n4);
        {
            ArcIt begin(g);
            ArcIt invalid(lemon::INVALID);
            should(begin!=lemon::INVALID);
            shouldEqual(std::distance(begin,invalid),8);            
            std::vector<Arc> arcVec(begin,invalid);
            shouldEqual(8,arcVec.size());
            shouldEqual(1,g.id(arcVec[0]));
            shouldEqual(2,g.id(arcVec[1]));
            shouldEqual(3,g.id(arcVec[2]));
            shouldEqual(4,g.id(arcVec[3]));
            shouldEqual(6,g.id(arcVec[4]));
            shouldEqual(7,g.id(arcVec[5]));
            shouldEqual(8,g.id(arcVec[6]));
            shouldEqual(9,g.id(arcVec[7]));
        }
        {
            ArcIt begin(g);
            should(begin!=lemon::INVALID);

            ArcIt empty;
            std::vector<Arc> arcVec(begin,empty);
            shouldEqual(8,arcVec.size());
            shouldEqual(1,g.id(arcVec[0]));
            shouldEqual(2,g.id(arcVec[1]));
            shouldEqual(3,g.id(arcVec[2]));
            shouldEqual(4,g.id(arcVec[3]));
            shouldEqual(6,g.id(arcVec[4]));
            shouldEqual(7,g.id(arcVec[5]));
            shouldEqual(8,g.id(arcVec[6]));
            shouldEqual(9,g.id(arcVec[7]));
        }
        {
            ArcIt begin(g,g.arcFromId(2));
            should(begin!=lemon::INVALID);

            ArcIt empty;
            std::vector<Arc> arcVec(begin,empty);
            shouldEqual(7,arcVec.size());
            shouldEqual(2,g.id(arcVec[0]));
            shouldEqual(3,g.id(arcVec[1]));
            shouldEqual(4,g.id(arcVec[2]));
            shouldEqual(6,g.id(arcVec[3]));
            shouldEqual(7,g.id(arcVec[4]));
            shouldEqual(8,g.id(arcVec[5]));
            shouldEqual(9,g.id(arcVec[6]));
        }
        {
            ArcIt begin(g,g.arcFromId(2));
            ArcIt end(g,g.arcFromId(3));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);    
            should(begin!=end);

            shouldEqual(std::distance(begin,end),1);
            std::vector<Arc> arcVec(begin,end);
            shouldEqual(1,arcVec.size());
            shouldEqual(2,g.id(arcVec[0]));
        }

        {
            ArcIt begin(g,g.arcFromId(2));
            ArcIt end(g,g.arcFromId(4));

            should(begin!=lemon::INVALID);
            should(end!=lemon::INVALID);

            std::vector<Arc> arcVec(begin,end);
            shouldEqual(2,arcVec.size());
            shouldEqual(2,g.id(arcVec[0]));
            shouldEqual(3,g.id(arcVec[1]));
        }
    }

    void adjGraphFindEdgeTest(){


        {

        GraphType g(0,0,false);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);


        e12 = g.findEdge(n1,n2);
        e13 = g.findEdge(n1,n3);
        e24 = g.findEdge(n2,n4);
        e34 = g.findEdge(n3,n4);



        should(e12!=lemon::INVALID);
        should(e13!=lemon::INVALID);
        should(e24!=lemon::INVALID);
        should(e34!=lemon::INVALID);


        should(g.findEdge(n2,n3)==lemon::INVALID);
        should(g.findEdge(n3,n2)==lemon::INVALID);
        should(g.findEdge(n1,n4)==lemon::INVALID);
        should(g.findEdge(n4,n1)==lemon::INVALID);

        }

        {
        GraphType g(0,0,true);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);


        e12 = g.findEdge(n1,n2);
        e13 = g.findEdge(n1,n3);
        e24 = g.findEdge(n2,n4);
        e34 = g.findEdge(n3,n4);



        should(e12!=lemon::INVALID);
        should(e13!=lemon::INVALID);
        should(e24!=lemon::INVALID);
        should(e34!=lemon::INVALID);


        should(g.findEdge(n2,n3)==lemon::INVALID);
        should(g.findEdge(n3,n2)==lemon::INVALID);
        should(g.findEdge(n1,n4)==lemon::INVALID);
        should(g.findEdge(n4,n1)==lemon::INVALID);
        }
    }


    void adjGraphIncEdgeItTestStart0()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g(0,0,true);

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);

        // get edges


        {
            IncEdgeIt a(g,n1);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e13)!=eSet.end());
            should(eSet.find(e12)!=eSet.end());
            should(eSet.find(e34)==eSet.end());
            should(eSet.find(e24)==eSet.end());
        }
        {
            IncEdgeIt a(g,n2);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e12)!=eSet.end());
            should(eSet.find(e24)!=eSet.end());
            should(eSet.find(e34)==eSet.end());
            should(eSet.find(e13)==eSet.end());
        }
        {
            IncEdgeIt a(g,n3);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e13)!=eSet.end());
            should(eSet.find(e34)!=eSet.end());
            should(eSet.find(e12)==eSet.end());
            should(eSet.find(e24)==eSet.end());
        }
        {
            IncEdgeIt a(g,n4);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e24)!=eSet.end());
            should(eSet.find(e34)!=eSet.end());
            should(eSet.find(e12)==eSet.end());
            should(eSet.find(e13)==eSet.end());
        }
    }
    void adjGraphIncEdgeItTestStart1()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);

        // get edges


        {
            IncEdgeIt a(g,n1);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e13)!=eSet.end());
            should(eSet.find(e12)!=eSet.end());
            should(eSet.find(e34)==eSet.end());
            should(eSet.find(e24)==eSet.end());
        }
        {
            IncEdgeIt a(g,n2);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e12)!=eSet.end());
            should(eSet.find(e24)!=eSet.end());
            should(eSet.find(e34)==eSet.end());
            should(eSet.find(e13)==eSet.end());
        }
        {
            IncEdgeIt a(g,n3);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e13)!=eSet.end());
            should(eSet.find(e34)!=eSet.end());
            should(eSet.find(e12)==eSet.end());
            should(eSet.find(e24)==eSet.end());
        }
        {
            IncEdgeIt a(g,n4);
            IncEdgeIt b(lemon::INVALID);
            should(a!=b);
            should(b==lemon::INVALID);
            should(a!=lemon::INVALID);
            shouldEqual(std::distance(a,b),2);

            std::set<Edge> eSet(a,b);
            shouldEqual(eSet.size(),2);
            should(eSet.find(e24)!=eSet.end());
            should(eSet.find(e34)!=eSet.end());
            should(eSet.find(e12)==eSet.end());
            should(eSet.find(e13)==eSet.end());
        }
    }

    void adjGraphInArcItTest()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);

        n1=g.nodeFromId(1);
        n2=g.nodeFromId(2);
        n3=g.nodeFromId(3);
        n4=g.nodeFromId(4);

        e12 = g.findEdge(n1,n2);
        e13 = g.findEdge(n1,n3);
        e24 = g.findEdge(n2,n4);
        e34 = g.findEdge(n3,n4);

        // get incoming arcs
        {
            InArcIt a(g,n1);
            InArcIt b(lemon::INVALID);
            // this node has NO in arcs !
            should(a==b);
            shouldEqual(std::distance(a,b),0);
        }
        {
            InArcIt a(g,n2);
            InArcIt b(lemon::INVALID);
            // this node has 1 incoming arc
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
            Arc arc(*a);
            Node source = g.source(arc);
            Node target = g.target(arc);
            shouldEqual(g.id(target),g.id(n2));
            shouldEqual(g.id(source),g.id(n1));
        }
        {
            InArcIt a(g,n3);
            InArcIt b(lemon::INVALID);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
            Arc arc(*a);
            Node source = g.source(arc);
            Node target = g.target(arc);
            shouldEqual(g.id(target),g.id(n3));
            shouldEqual(g.id(source),g.id(n1));
        }
        {
            InArcIt a(g,n4);
            InArcIt b(lemon::INVALID);
            // this node has 1 incoming arc
            should(a!=b);
            shouldEqual(std::distance(a,b),2);
            Arc arc1(*a); ++a;
            Arc arc2(*a); ++a;
            should(a==b);
            should(a==lemon::INVALID);
            Node source1 = g.source(arc1);
            Node target1 = g.target(arc1);
            Node source2 = g.source(arc2);
            Node target2 = g.target(arc2);
            shouldEqual(g.id(target1),g.id(n4));
            shouldEqual(g.id(target2),g.id(n4));
            should(source1 !=source2 );
            should(source1==n2 || source2==n2);
            should(source1==n3 || source2==n3);

        }
    
    }


    void adjGraphOutArcItTest()
    {
        // 1 |3
        // __ __
        // 2 |4
        // create g
        GraphType g;

        Node n1=g.addNode();
        Node n2=g.addNode();
        Node n3=g.addNode();
        Node n4=g.addNode();

        Edge e12 = g.addEdge(n1,n2);
        Edge e13 = g.addEdge(n1,n3);
        Edge e24 = g.addEdge(n2,n4);
        Edge e34 = g.addEdge(n3,n4);

        n1=g.nodeFromId(1);
        n2=g.nodeFromId(2);
        n3=g.nodeFromId(3);
        n4=g.nodeFromId(4);

        e12 = g.findEdge(n1,n2);
        e13 = g.findEdge(n1,n3);
        e24 = g.findEdge(n2,n4);
        e34 = g.findEdge(n3,n4);


        {
            OutArcIt a(g,n1);
            OutArcIt b(lemon::INVALID);
            should(a!=b);
            shouldEqual(std::distance(a,b),2);
            Arc arc1(*a); ++a;
            Arc arc2(*a); ++a;
            should(a==b);
            should(a==lemon::INVALID);
            Node source1 = g.source(arc1);
            Node target1 = g.target(arc1);
            Node source2 = g.source(arc2);
            Node target2 = g.target(arc2);
            shouldEqual(g.id(source1),g.id(n1));
            shouldEqual(g.id(source2),g.id(n1));
            should(target1 !=target2 );
            should(target1==n2 || target2==n2);
            should(target1==n3 || target2==n3);
        }
        {
            OutArcIt a(g,n2);
            OutArcIt b(lemon::INVALID);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
            Arc arc(*a);
            Node source = g.source(arc);
            Node target = g.target(arc);
            shouldEqual(g.id(source),g.id(n2));
            shouldEqual(g.id(target),g.id(n4));
            
        }
        {
            OutArcIt a(g,n3);
            OutArcIt b(lemon::INVALID);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
            Arc arc(*a);
            Node source = g.source(arc);
            Node target = g.target(arc);
            shouldEqual(g.id(target),g.id(n4));
            shouldEqual(g.id(source),g.id(n3));
        }
        {
            OutArcIt a(g,n4);
            OutArcIt b(lemon::INVALID);
            should(a==lemon::INVALID);
            should(a==b);
            shouldEqual(std::distance(a,b),0);
        }
    
    }

};


 
struct AdjacencyListGraphTestSuite
: public vigra::test_suite
{
    AdjacencyListGraphTestSuite()
    : vigra::test_suite("AdjacencyListGraphTestSuite")
    {   
        


        add( testCase( &AdjacencyListGraphOffTest::adjGraphSimpleTestStart0));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphSimpleTestStart1));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphNodeItTestStart0));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphNodeItTestStart1));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphEdgeItTestStart0));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphEdgeItTestStart1));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphFindEdgeTest));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphIncEdgeItTestStart0));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphIncEdgeItTestStart1));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphTestStart0));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphTestStart1));

        add( testCase( &AdjacencyListGraphOffTest::adjGraphArcTest));
        add( testCase( &AdjacencyListGraphOffTest::adjGraphArcItTest));
        //add( testCase( &AdjacencyListGraphOffTest::adjGraphInArcItTest));
        //add( testCase( &AdjacencyListGraphOffTest::adjGraphOutArcItTest));


    }
};

int main(int argc, char ** argv)
{
    AdjacencyListGraphTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

