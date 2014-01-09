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
#include "vigra/rag/rag.hxx"
#include "vigra/rag/rag.hxx"

using namespace vigra;



struct EnumerationIteratorTest{



    void emptyConstructorTest(){
        typedef EnumerationIterator<Int32> IterType;

        IterType invalidIter = lemon::INVALID;
        should(invalidIter.isEnd());
        should(invalidIter==lemon::INVALID);


        IterType emptyIter;
        should(emptyIter.isEnd());
        should(emptyIter==lemon::INVALID);

        IterType zeroIter(0,0);
        should(zeroIter.isEnd());
        should(zeroIter==lemon::INVALID);

        IterType endIter(3,3);
        should(endIter.isEnd());
        should(endIter==lemon::INVALID);
    }


    void fillTest(){
        typedef EnumerationIterator<Int32> IterType;

        {
            IterType iter(0,3);
            std::vector<Int32> vec(iter,IterType(3,3));
            should(vec.size()==3);


            should(vec[0]==0);
            should(vec[1]==1);
            should(vec[2]==2);
        }
        {
            IterType iter(0,0);
            std::vector<Int32> vec(iter,iter);
            should(vec.size()==0);
        }

    }

    void basicTest(){
        typedef EnumerationIterator<Int32> IterType;

        IterType iter(0,3);
        IterType endIter(3,3);

        shouldEqual(std::distance(iter,endIter),3);


        should(endIter.isEnd());
        should(endIter==lemon::INVALID);


        should(iter.isBegin());
        should(iter!=endIter);
        should(!iter.isEnd());
        should(iter!=lemon::INVALID);
        should(*iter==0);

        ++iter;
        should(iter!=endIter);
        should(!iter.isEnd());
        should(iter!=lemon::INVALID);
        should(*iter==1);

        ++iter;
        should(iter!=endIter);
        should(!iter.isEnd());
        should(iter!=lemon::INVALID);
        should(*iter==2);

        ++iter;
        should(iter==endIter);
        should(iter.isEnd());
        should(iter==lemon::INVALID);

    }
    void distanceTest(){
        typedef EnumerationIterator<Int32> IterType;

        {
            IterType a(0,3);
            IterType b(3,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),3);
        }
        {
            IterType a(0,3);
            IterType b(2,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),2);
        }
        {
            IterType a(0,3);
            IterType b(1,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
        }
        {
            IterType a(0,3);
            IterType b(0,3);
            should(a==b);
            shouldEqual(std::distance(a,b),0);
        }


        {
            IterType a(1,3);
            IterType b(3,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),2);
        }
        {
            IterType a(2,3);
            IterType b(3,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
        }
        {
            IterType a(3,3);
            IterType b(3,3);
            should(a==b);
            shouldEqual(std::distance(a,b),0);
        }

        {
            IterType a(2,3);
            IterType b(3,3);
            should(a!=b);
            shouldEqual(std::distance(a,b),1);
        }

    }
};



template<class IN_LABEL_TYPE>
struct Rag2Test{

    typedef IN_LABEL_TYPE                        InLabelType;
    typedef vigra::Rag<2,InLabelType>            RagType;
    typedef typename RagType::Node               Node;
    typedef typename RagType::Edge               Edge;
    typedef typename RagType::Arc                Arc;
    typedef typename RagType::EdgeIt             EdgeIt;
    typedef typename RagType::NodeIt             NodeIt;
    typedef typename RagType::ArcIt              ArcIt;
    typedef typename RagType::IncEdgeIt          IncEdgeIt;
    typedef typename RagType::InArcIt            InArcIt;
    typedef typename RagType::OutArcIt           OutArcIt;

    typedef typename RagType::NeighborNodeIt     NeighborNodeIt;
    typedef typename RagType::InputLabelingView  InputLabelingView;
    typedef typename RagType::InputLabelingArray InputLabelingArray;
    Rag2Test(){

    }



    void ragTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4

        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;


        // create rag
        RagType rag(labels);


        
        // assert basic sizes
        should(rag.edgeNum()==4);
        should(rag.nodeNum()==4);
        should(rag.arcNum()==8);
        should(rag.maxEdgeId()==4);
        should(rag.maxNodeId()==4);
        should(rag.maxArcId()==8);

        should( rag.findEdge(rag.nodeFromId(1),rag.nodeFromId(3) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(1),rag.nodeFromId(2) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(2),rag.nodeFromId(4) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(3),rag.nodeFromId(4) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(1),rag.nodeFromId(4) )==lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(2),rag.nodeFromId(3) )==lemon::INVALID);

        NodeIt nbegin(rag);
        NodeIt nend(lemon::INVALID);

        EdgeIt ebegin(rag);
        EdgeIt eend(lemon::INVALID);


        ArcIt abegin(rag);
        ArcIt aend(lemon::INVALID);

        std::vector<Node> allNodes(nbegin,nend);
        should(allNodes.size()==4);

        std::vector<Edge> allEdges(ebegin,eend);
        should(allEdges.size()==4);


        should(1==rag.id( allEdges[0] ) );
        should(2==rag.id( allEdges[1] ) );
        should(3==rag.id( allEdges[2] ) );
        should(4==rag.id( allEdges[3] ) );
        

        std::vector<Arc>  allArcs( abegin,aend);
        should(allArcs.size() ==8);
        for(size_t a=0;a<allArcs.size();++a){
            should(a+1==rag.id( allArcs[a] ) );
        }
        

        should( rag.id(rag.source(allArcs[0]))==rag.id(rag.u(allEdges[0])));
        should( rag.id(rag.source(allArcs[1]))==rag.id(rag.u(allEdges[1])));
        should( rag.id(rag.source(allArcs[2]))==rag.id(rag.u(allEdges[2])));
        should( rag.id(rag.source(allArcs[3]))==rag.id(rag.u(allEdges[3])));

        should( rag.id(rag.target(allArcs[0]))==rag.id(rag.v(allEdges[0])));
        should( rag.id(rag.target(allArcs[1]))==rag.id(rag.v(allEdges[1])));
        should( rag.id(rag.target(allArcs[2]))==rag.id(rag.v(allEdges[2])));
        should( rag.id(rag.target(allArcs[3]))==rag.id(rag.v(allEdges[3])));


        should( rag.id(rag.source(allArcs[0]))==rag.id(rag.target(allArcs[0+4])));
        should( rag.id(rag.source(allArcs[1]))==rag.id(rag.target(allArcs[1+4])));
        should( rag.id(rag.source(allArcs[2]))==rag.id(rag.target(allArcs[2+4])));
        should( rag.id(rag.source(allArcs[3]))==rag.id(rag.target(allArcs[3+4])));

        should( rag.id(rag.target(allArcs[0]))==rag.id(rag.source(allArcs[0+4])));
        should( rag.id(rag.target(allArcs[1]))==rag.id(rag.source(allArcs[1+4])));
        should( rag.id(rag.target(allArcs[2]))==rag.id(rag.source(allArcs[2+4])));
        should( rag.id(rag.target(allArcs[3]))==rag.id(rag.source(allArcs[3+4])));


        should( rag.id(rag.source(allArcs[4]))==rag.id(rag.v(allEdges[0])));
        should( rag.id(rag.source(allArcs[5]))==rag.id(rag.v(allEdges[1])));
        should( rag.id(rag.source(allArcs[6]))==rag.id(rag.v(allEdges[2])));
        should( rag.id(rag.source(allArcs[7]))==rag.id(rag.v(allEdges[3])));

        should( rag.id(rag.target(allArcs[4]))==rag.id(rag.u(allEdges[0])));
        should( rag.id(rag.target(allArcs[5]))==rag.id(rag.u(allEdges[1])));
        should( rag.id(rag.target(allArcs[6]))==rag.id(rag.u(allEdges[2])));
        should( rag.id(rag.target(allArcs[7]))==rag.id(rag.u(allEdges[3])));
    }
    void ragArcTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);
        
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
        should(  g.source(g.arcFromId(5)) == g.v(g.edgeFromId(1)) );
        should(  g.source(g.arcFromId(6)) == g.v(g.edgeFromId(2)) );
        should(  g.source(g.arcFromId(7)) == g.v(g.edgeFromId(3)) );
        should(  g.source(g.arcFromId(8)) == g.v(g.edgeFromId(4)) );

        should(  g.target(g.arcFromId(5)) == g.u(g.edgeFromId(1)) );
        should(  g.target(g.arcFromId(6)) == g.u(g.edgeFromId(2)) );
        should(  g.target(g.arcFromId(7)) == g.u(g.edgeFromId(3)) );
        should(  g.target(g.arcFromId(8)) == g.u(g.edgeFromId(4)) );

        // check that arcs are convertible to edges
        should(Edge(g.arcFromId(1))==g.edgeFromId(1));
        should(Edge(g.arcFromId(2))==g.edgeFromId(2));
        should(Edge(g.arcFromId(3))==g.edgeFromId(3));
        should(Edge(g.arcFromId(4))==g.edgeFromId(4));
        should(Edge(g.arcFromId(5))==g.edgeFromId(1));
        should(Edge(g.arcFromId(6))==g.edgeFromId(2));
        should(Edge(g.arcFromId(7))==g.edgeFromId(3));
        should(Edge(g.arcFromId(8))==g.edgeFromId(4));

    }

    void ragEdgeItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 1 |3
        // __ __
        // 2 |4

        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        

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
    void ragNodeItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);
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
    void ragArcItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);
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
            shouldEqual(5,g.id(arcVec[4]));
            shouldEqual(6,g.id(arcVec[5]));
            shouldEqual(7,g.id(arcVec[6]));
            shouldEqual(8,g.id(arcVec[7]));
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
            shouldEqual(5,g.id(arcVec[4]));
            shouldEqual(6,g.id(arcVec[5]));
            shouldEqual(7,g.id(arcVec[6]));
            shouldEqual(8,g.id(arcVec[7]));
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
            shouldEqual(5,g.id(arcVec[3]));
            shouldEqual(6,g.id(arcVec[4]));
            shouldEqual(7,g.id(arcVec[5]));
            shouldEqual(8,g.id(arcVec[6]));
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

    void ragFindEdgeTest(){
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        const Node n1=g.nodeFromId(1);
        const Node n2=g.nodeFromId(2);
        const Node n3=g.nodeFromId(3);
        const Node n4=g.nodeFromId(4);

        const Edge e12 = g.findEdge(n1,n2);
        const Edge e13 = g.findEdge(n1,n3);
        const Edge e24 = g.findEdge(n2,n4);
        const Edge e34 = g.findEdge(n3,n4);



        should(e12!=lemon::INVALID);
        should(e13!=lemon::INVALID);
        should(e24!=lemon::INVALID);
        should(e34!=lemon::INVALID);


        should(g.findEdge(n2,n3)==lemon::INVALID);
        should(g.findEdge(n3,n2)==lemon::INVALID);
        should(g.findEdge(n1,n4)==lemon::INVALID);
        should(g.findEdge(n4,n1)==lemon::INVALID);

    }

    void ragUVOrderTest(){
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        const Node n1=g.nodeFromId(1);
        const Node n2=g.nodeFromId(2);
        const Node n3=g.nodeFromId(3);
        const Node n4=g.nodeFromId(4);

        const Edge e12 = g.findEdge(n1,n2);
        const Edge e13 = g.findEdge(n1,n3);
        const Edge e24 = g.findEdge(n2,n4);
        const Edge e34 = g.findEdge(n3,n4);


        // for rag id(u(edge)) < id(v(edge));
        should(g.id(g.u(e12))<g.id(g.v(e12)));
        should(g.id(g.u(e12))<g.id(g.v(e12)));
        should(g.id(g.u(e12))<g.id(g.v(e12)));
        should(g.id(g.u(e12))<g.id(g.v(e12)));
    }

    void ragIncEdgeItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        Node n1=g.nodeFromId(1);
        Node n2=g.nodeFromId(2);
        Node n3=g.nodeFromId(3);
        Node n4=g.nodeFromId(4);

        Edge e12 = g.findEdge(n1,n2);
        Edge e13 = g.findEdge(n1,n3);
        Edge e24 = g.findEdge(n2,n4);
        Edge e34 = g.findEdge(n3,n4);

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


    void ragInArcItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        Node n1=g.nodeFromId(1);
        Node n2=g.nodeFromId(2);
        Node n3=g.nodeFromId(3);
        Node n4=g.nodeFromId(4);

        Edge e12 = g.findEdge(n1,n2);
        Edge e13 = g.findEdge(n1,n3);
        Edge e24 = g.findEdge(n2,n4);
        Edge e34 = g.findEdge(n3,n4);

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


    void ragOutArcItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 1 |3
        // __ __
        // 2 |4
        labels(0,0)=1;
        labels(0,1)=2;
        labels(1,0)=3;
        labels(1,1)=4;
        // create rag
        RagType g(labels);

        Node n1=g.nodeFromId(1);
        Node n2=g.nodeFromId(2);
        Node n3=g.nodeFromId(3);
        Node n4=g.nodeFromId(4);

        Edge e12 = g.findEdge(n1,n2);
        Edge e13 = g.findEdge(n1,n3);
        Edge e24 = g.findEdge(n2,n4);
        Edge e34 = g.findEdge(n3,n4);


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


 
struct RagTestSuite
: public vigra::test_suite
{
    RagTestSuite()
    : vigra::test_suite("RagTestSuite")
    {   
        

        add( testCase( &EnumerationIteratorTest::emptyConstructorTest));
        add( testCase( &EnumerationIteratorTest::distanceTest) );
        add( testCase( &EnumerationIteratorTest::basicTest) );
        add( testCase( &EnumerationIteratorTest::fillTest));

        add( testCase( &Rag2Test<vigra::UInt32>::ragTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragArcTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragFindEdgeTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragUVOrderTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragEdgeItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragNodeItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragArcItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragIncEdgeItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragInArcItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragOutArcItTest));


    }
};

int main(int argc, char ** argv)
{
    RagTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

