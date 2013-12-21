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
#include "vigra/merge_graph/merge_graph.hxx"
#include "vigra/multi_array.hxx"

using namespace vigra;

template<class LABEL_TYPE>
struct MergeGraphTest
{
    typedef LABEL_TYPE LabelType;
    typedef vigra::MergeGraph<LabelType> MergeGraphType;

    MergeGraphTest()
    {

    }
    
    void constructorTest()
    {
        MergeGraphType graph(4,5);
        // 2x2 grid graph
        //  0 | 1
        //  _   _ 
        //  2 | 3

        // set inital edges (with one double edge)
        graph.setInitalEdge(0,0,1);
        graph.setInitalEdge(1,2,3);
        graph.setInitalEdge(2,0,2);
        graph.setInitalEdge(3,1,3);
        graph.setInitalEdge(4,1,3);
        should(graph.numberOfNodes() == 4);
        should(graph.numberOfEdges() == 5);

        // check edge exists
        should(graph.hasEdge(0));
        should(graph.hasEdge(1));
        should(graph.hasEdge(2));
        should(graph.hasEdge(3));
        should(graph.hasEdge(4));
        // check edges nodes
        should(graph.getEdge(0)[0]==0);
        should(graph.getEdge(0)[1]==1);
        should(graph.getEdge(1)[0]==2);
        should(graph.getEdge(1)[1]==3);
        should(graph.getEdge(2)[0]==0);
        should(graph.getEdge(2)[1]==2);
        should(graph.getEdge(3)[0]==1);
        should(graph.getEdge(3)[1]==3);
        should(graph.getEdge(4)[0]==1);
        should(graph.getEdge(4)[1]==3);

        // check nodes exist
        should(graph.hasNode(0));
        should(graph.hasNode(1));
        should(graph.hasNode(2));
        should(graph.hasNode(3));
        
        // check the number edges for each node
        should(graph.getNode(0).numberOfEdges( )==2);
        should(graph.getNode(1).numberOfEdges( )==3);
        should(graph.getNode(2).numberOfEdges( )==2);
        should(graph.getNode(3).numberOfEdges( )==3);

        // check edges
        should(graph.getNode(0).hasEdge(0));
        should(graph.getNode(0).hasEdge(2));

        should(graph.getNode(1).hasEdge(0));
        should(graph.getNode(1).hasEdge(3));
        should(graph.getNode(1).hasEdge(4));

        should(graph.getNode(2).hasEdge(1));
        should(graph.getNode(2).hasEdge(2));

        should(graph.getNode(3).hasEdge(1));
        should(graph.getNode(3).hasEdge(3));
        should(graph.getNode(3).hasEdge(4));


        // merge merge Parallel Edges 
        graph.mergeParallelEdges();
        should(graph.numberOfNodes() == 4);
        should(graph.numberOfEdges() == 4);
        should(graph.reprEdge(3) == graph.reprEdge(4));

        // check the number edges for each node
        // (has changed since we merge edges)
        should(graph.getNode(0).numberOfEdges( )==2);
        should(graph.getNode(1).numberOfEdges( )==2);
        should(graph.getNode(2).numberOfEdges( )==2);
        should(graph.getNode(3).numberOfEdges( )==2);

        // check edges
        should(graph.getNode(0).hasEdge(0));
        should(graph.getNode(0).hasEdge(2));

        should(graph.getNode(1).hasEdge(0));
        should(graph.getNode(1).hasEdge(graph.reprEdge(3)));
        should(graph.getNode(1).hasEdge(graph.reprEdge(4)));

        should(graph.getNode(2).hasEdge(1));
        should(graph.getNode(2).hasEdge(2));

        should(graph.getNode(3).hasEdge(1));
        should(graph.getNode(3).hasEdge(graph.reprEdge(3)));
        should(graph.getNode(3).hasEdge(graph.reprEdge(4)));

        // check which edge is the deleted
        LabelType deletedEdge = graph.reprEdge(3)==3 ? 4 : 3;
        should(!graph.hasEdge(deletedEdge));

    }
    void mergeTest()
    {
        MergeGraphType graph(3,3);
        // triangle graph

        // set inital edges (with one double edge)
        graph.setInitalEdge(0,0,1);
        graph.setInitalEdge(1,0,2);
        graph.setInitalEdge(2,1,2);
        should(graph.numberOfNodes() == 3);
        should(graph.numberOfEdges() == 3);

        // merge merge Parallel Edges  (there are non)
        graph.mergeParallelEdges();
        should(graph.numberOfNodes() == 3);
        should(graph.numberOfEdges() == 3);



        should(graph.reprEdge(0)==0);
        should(graph.reprEdge(1)==1);
        should(graph.reprEdge(2)==2);

        for(size_t e=0;e<3;++e){
            std::cout<<"be="<<e<<" ar="<<graph.reprEdge(e)<<"\n";
        }

        // merge edge 0 
        graph.mergeRegions(0);
        should(graph.numberOfNodes() == 2);
        should(graph.numberOfEdges() == 1);

        for(size_t e=0;e<3;++e){
            std::cout<<"e="<<e<<" br="<<graph.reprEdge(e)<<"\n";
        }

        for(size_t e=0;e<3;++e){
            std::cout<<"n="<<e<<" br="<<graph.reprNode(e)<<"\n";
        }

        should(graph.reprNode(0)==graph.reprNode(1));
        should(graph.reprNode(0)!=graph.reprNode(2));
        should(graph.reprNode(1)!=graph.reprNode(2));

        should(graph.reprEdge(1)==graph.reprEdge(2));
        should(!graph.reprEdge(0)==graph.reprEdge(1));


        for(size_t e=0;e<3;++e){
            std::cout<<"e="<<e<<" cr="<<graph.reprEdge(e)<<"\n";
        }

        should(graph.reprEdge(0)==0);
        should(graph.reprEdge(1)!=0);
        should(graph.reprEdge(2)!=0);

        should(graph.reprEdge(0)!=graph.reprEdge(2));
    }
};

        
struct MergeGraphTestSuite
: public vigra::test_suite
{
    MergeGraphTestSuite()
    : vigra::test_suite("MergeGraphTestSuite")
    {   
        add( testCase( &MergeGraphTest<vigra::UInt32>::constructorTest));
        add( testCase( &MergeGraphTest<vigra::UInt64>::constructorTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::constructorTest));
        add( testCase( &MergeGraphTest<vigra::Int64>::constructorTest));

        add( testCase( &MergeGraphTest<vigra::UInt32>::mergeTest));
        add( testCase( &MergeGraphTest<vigra::UInt64>::mergeTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::mergeTest));
        add( testCase( &MergeGraphTest<vigra::Int64>::mergeTest));
    }
};

int main(int argc, char ** argv)
{
    MergeGraphTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

