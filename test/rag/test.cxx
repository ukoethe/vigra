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
    typedef typename RagType::NeighborNodeIt     NeighborNodeIt;
    typedef typename RagType::InputLabelingView  InputLabelingView;
    typedef typename RagType::InputLabelingArray InputLabelingArray;
    Rag2Test(){

    }



    void ragTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));

        // 0 |2
        // __ __
        // 1 |3

        labels(0,0)=0;
        labels(0,1)=1;
        labels(1,0)=2;
        labels(1,1)=3;


        // create rag
        RagType rag(labels);


        
        // assert basic sizes
        should(rag.edgeNum()==4);
        should(rag.nodeNum()==4);
        should(rag.arcNum()==8);
        should(rag.maxEdgeId()==3);
        should(rag.maxNodeId()==3);
        should(rag.maxArcId()==7);

        should( rag.findEdge(rag.nodeFromId(0),rag.nodeFromId(2) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(0),rag.nodeFromId(1) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(1),rag.nodeFromId(3) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(2),rag.nodeFromId(3) )!=lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(0),rag.nodeFromId(3) )==lemon::INVALID);
        should( rag.findEdge(rag.nodeFromId(1),rag.nodeFromId(2) )==lemon::INVALID);


        std::vector<Node> allNodes(rag.nodesBegin(),rag.nodesEnd());
        should(allNodes.size()==4);

        std::vector<Edge> allEdges(rag.edgesBegin(),rag.edgesEnd());
        should(allEdges.size()==4);


        should(0==rag.id( allEdges[0] ) );
        should(1==rag.id( allEdges[1] ) );
        should(2==rag.id( allEdges[2] ) );
        should(3==rag.id( allEdges[3] ) );
        

        std::vector<Arc>  allArcs( rag.arcsBegin() ,rag.arcsEnd());
        should(allArcs.size() ==8);
        for(size_t a=0;a<allArcs.size();++a){
            should(a==rag.id( allArcs[a] ) );
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


    void ragEdgeItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 0 |2
        // __ __
        // 1 |3
        labels(0,0)=0;
        labels(0,1)=1;
        labels(1,0)=2;
        labels(1,1)=3;
        // create rag
        RagType g(labels);

        

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
    void ragNodeItTest()
    {
        // create labels
        InputLabelingArray labels(typename InputLabelingArray::difference_type(2,2));
        // 0 |2
        // __ __
        // 1 |3
        labels(0,0)=0;
        labels(0,1)=1;
        labels(1,0)=2;
        labels(1,1)=3;
        // create rag
        RagType g(labels);
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
        add( testCase( &Rag2Test<vigra::UInt32>::ragEdgeItTest));
        add( testCase( &Rag2Test<vigra::UInt32>::ragNodeItTest));

    }
};

int main(int argc, char ** argv)
{
    RagTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

