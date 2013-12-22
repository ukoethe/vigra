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
    
    void mergeSimpleDoubleEdgeTest()
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

        {
        std::vector<LabelType> activeEdgeVec(graph.edgesBegin(),graph.edgesEnd());
        should(activeEdgeVec.size() == 3);
        should(activeEdgeVec[0] == 0);
        should(activeEdgeVec[1] == 1);
        should(activeEdgeVec[2] == 2);
        }

        std::vector<LabelType> activeNodes(graph.nodesBegin(),graph.nodesEnd());
        should(activeNodes.size() == 3);
        should(activeNodes[0] == 0);
        should(activeNodes[1] == 1);
        should(activeNodes[2] == 2);

        should(graph.numberOfNodes() == 3);
        should(graph.numberOfEdges() == 3);

        // merge merge Parallel Edges  (there are non)
        graph.mergeParallelEdges();
        should(graph.numberOfNodes() == 3);
        should(graph.numberOfEdges() == 3);



        should(graph.reprEdge(0)==0);
        should(graph.reprEdge(1)==1);
        should(graph.reprEdge(2)==2);



        // merge edge 0 
        graph.mergeRegions(0);
        should(graph.numberOfNodes() == 2);
        should(graph.numberOfEdges() == 1);
        {
        std::vector<LabelType> activeEdgeVec(graph.edgesBegin(),graph.edgesEnd());
        should(activeEdgeVec.size() == 1);
        }

        should(graph.reprNode(0)==graph.reprNode(1));
        should(graph.reprNode(0)!=graph.reprNode(2));
        should(graph.reprNode(1)!=graph.reprNode(2));

        should(graph.reprEdge(1)==graph.reprEdge(2));
        should(!graph.reprEdge(0)==graph.reprEdge(1));

        should(graph.reprEdge(0)==0);
        should(graph.reprEdge(1)!=0);
        should(graph.reprEdge(2)!=0);

        should(graph.reprEdge(0)!=graph.reprEdge(2));
    }
    void chainTest()
    {

        const size_t nChainNodes  = 10;
        const size_t nChainEdges  = nChainNodes-1;
        // 0 - 1 - 2 - .. nChainNodes-1
        MergeGraphType graph(nChainNodes,nChainEdges);
        // triangle graph

        // set inital edges (without any double edge)
        for(size_t e=0;e<nChainNodes-1;++e){
            graph.setInitalEdge(e,e,e+1);
        }
        should(graph.numberOfNodes() == nChainNodes);
        should(graph.numberOfEdges() == nChainEdges);

        // remove edges from 0 to nChainNodes -1
        for(size_t e=0;e<nChainNodes-1;++e){
            should(graph.numberOfEdges()==nChainEdges-e);
            should(graph.numberOfNodes()==nChainNodes-e);
            // check that edge is there 
            should(graph.hasEdge(e));
            // fist node is rep of e, second node still untouched e+1
            should(graph.getEdge(e)[0]==graph.reprNode(e));
            should(graph.getEdge(e)[1]==e+1);

            // remove the edge
            graph.mergeRegions(e);

            should(!graph.hasEdge(e));
            should(graph.numberOfEdges()==nChainEdges-e-1);
            should(graph.numberOfNodes()==nChainNodes-e-1);
        }
    }

    void testGrid()
    {
        // 3x3
        //
        //   0 | 1 | 2 
        //   _  _   _
        //   3 | 4 | 5
        //   _   _   _
        //   6 | 7 | 8
        const size_t nChainNodes  = 9;
        const size_t nChainEdges  = 12;
        // 0 - 1 - 2 - .. nChainNodes-1
        MergeGraphType graph(nChainNodes,nChainEdges);
        graph.setInitalEdge(0, 0,1);
        graph.setInitalEdge(1, 1,2);

        graph.setInitalEdge(2, 3,4);
        graph.setInitalEdge(3, 4,5);

        graph.setInitalEdge(4, 6,7);
        graph.setInitalEdge(5, 7,8);

        graph.setInitalEdge(6, 0,3);
        graph.setInitalEdge(7, 1,4);
        graph.setInitalEdge(8, 2,5);

        graph.setInitalEdge(9 , 3,6);
        graph.setInitalEdge(10, 4,7);
        graph.setInitalEdge(11, 5,8);

        should(graph.numberOfNodes()==9);
        should(graph.numberOfEdges()==12);
        

        // check inital values bevore any merges
        bool edgeStateTrue[12]  ={1,1,1,1,1,1,1,1,1,1,1,1};
        bool edgeStateCheck[12] ={0,0,0,0,0,0,0,0,0,0,0,0};


        graph.stateOfInitalEdges(edgeStateCheck,edgeStateCheck+12);
        shouldEqualSequence(edgeStateTrue,edgeStateTrue+12,edgeStateCheck);


        // merge edge between 3-4
        // this will reduce the number of active edges by 1:
        // edge 2 will dissaper
        graph.mergeRegions(2);
        should(graph.numberOfNodes()==8);
        should(graph.numberOfEdges()==11);
        graph.stateOfInitalEdges(edgeStateCheck,edgeStateCheck+12);
        edgeStateTrue[2]=0;
        shouldEqualSequence(edgeStateTrue,edgeStateTrue+12,edgeStateCheck);





    }

};


template<class LABEL_TYPE>
struct PartitonTest
{
    typedef LABEL_TYPE LabelType;
    typedef vigra::detail_merge_graph::Partition<LabelType> PartitionType;
    typedef std::set<LabelType> SetType;
    typedef std::vector<LabelType> VecType;
    PartitonTest()
    {

    }

    void trueReps(const PartitionType ufd,SetType &set){
        set.clear();
        for(LabelType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
    }

    void trueReps(const PartitionType ufd,SetType &set,SetType & c){
        set.clear();
        for(LabelType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
        for(typename  SetType::const_iterator iter=c.begin();iter!=c.end();++iter){
            const LabelType toRemove=*iter;
            should(set.find(toRemove)!=set.end());
            set.erase(toRemove);
        }

    }


    void testReps(const PartitionType ufd,VecType & vec){
        vec.clear();
        vec.assign(ufd.begin(),ufd.end());
    }




    void iteratorTest1(){
        PartitionType ufd(6);
        SetType trueRep;
        VecType testRep;

        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

            

        ufd.merge(0,1);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
        
        ufd.merge(0,2);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
        
        ufd.merge(0,3);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(3,3);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(4,5);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(3,5);
        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
    }

    void iteratorTest2(){
        PartitionType ufd(6);
        SetType trueRep;
        VecType testRep;

        trueReps(ufd,trueRep);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        SetType erased;
        erased.insert(0);
        ufd.eraseElement(0);

        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        ufd.merge(1,2);
        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());

        LabelType rep12 = ufd.find(1);
        erased.insert(rep12);
        ufd.eraseElement(rep12);
        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
    }
    
};

        
struct MergeGraphTestSuite
: public vigra::test_suite
{
    MergeGraphTestSuite()
    : vigra::test_suite("MergeGraphTestSuite")
    {   

        add( testCase( &PartitonTest<vigra::UInt32>::iteratorTest1));
        add( testCase( &PartitonTest<vigra::UInt32>::iteratorTest2));

        add( testCase( &MergeGraphTest<vigra::UInt32>::mergeSimpleDoubleEdgeTest));
        add( testCase( &MergeGraphTest<vigra::UInt32>::mergeTest));
        add( testCase( &MergeGraphTest<vigra::UInt32>::chainTest));
        add( testCase( &MergeGraphTest<vigra::UInt32>::testGrid));


    }
};

int main(int argc, char ** argv)
{
    MergeGraphTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

