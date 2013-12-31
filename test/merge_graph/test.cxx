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

#include "vigra/merge_graph/merge_graph.hxx"
#include "vigra/merge_graph/min_indexed_pq.hxx"
#include "vigra/merge_graph/maps/accumulator_map.hxx"


using namespace vigra;

template<class ID_TYPE>
struct MergeGraphTest
{
    typedef ID_TYPE IdType;
    typedef vigra::MergeGraph<IdType> MergeGraphType;
    typedef typename MergeGraphType::Node Node;
    typedef typename MergeGraphType::Edge Edge;
    typedef std::vector<IdType> Lvec;
    typedef std::set<IdType>    Lset;
    typedef std::vector<Edge>       EVec;
    typedef typename MergeGraphType::EdgeIdIt EdgeIdIt;
    typedef typename MergeGraphType::EdgeIt EdgeIt;
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
        should(graph.hasEdgeId(0));
        should(graph.hasEdgeId(1));
        should(graph.hasEdgeId(2));
        should(graph.hasEdgeId(3));
        should(graph.hasEdgeId(4));
        // check edges nodes
        should(graph.edgeFromId(0)[0]==0);
        should(graph.edgeFromId(0)[1]==1);
        should(graph.edgeFromId(1)[0]==2);
        should(graph.edgeFromId(1)[1]==3);
        should(graph.edgeFromId(2)[0]==0);
        should(graph.edgeFromId(2)[1]==2);
        should(graph.edgeFromId(3)[0]==1);
        should(graph.edgeFromId(3)[1]==3);
        should(graph.edgeFromId(4)[0]==1);
        should(graph.edgeFromId(4)[1]==3);

        // check nodes exist
        should(graph.hasNodeId(0));
        should(graph.hasNodeId(1));
        should(graph.hasNodeId(2));
        should(graph.hasNodeId(3));
        
        // check the number edges for each node
        should(graph.nodeFromId(0).numberOfEdges( )==2);
        should(graph.nodeFromId(1).numberOfEdges( )==3);
        should(graph.nodeFromId(2).numberOfEdges( )==2);
        should(graph.nodeFromId(3).numberOfEdges( )==3);

        // check edges
        should(graph.nodeFromId(0).hasEdgeId(0));
        should(graph.nodeFromId(0).hasEdgeId(2));

        should(graph.nodeFromId(1).hasEdgeId(0));
        should(graph.nodeFromId(1).hasEdgeId(3));
        should(graph.nodeFromId(1).hasEdgeId(4));

        should(graph.nodeFromId(2).hasEdgeId(1));
        should(graph.nodeFromId(2).hasEdgeId(2));

        should(graph.nodeFromId(3).hasEdgeId(1));
        should(graph.nodeFromId(3).hasEdgeId(3));
        should(graph.nodeFromId(3).hasEdgeId(4));


        // merge merge Parallel Edges 
        graph.mergeParallelEdges();
        should(graph.numberOfNodes() == 4);
        should(graph.numberOfEdges() == 4);
        should(graph.reprEdgeId(3) == graph.reprEdgeId(4));

        // check the number edges for each node
        // (has changed since we merge edges)
        should(graph.nodeFromId(0).numberOfEdges( )==2);
        should(graph.nodeFromId(1).numberOfEdges( )==2);
        should(graph.nodeFromId(2).numberOfEdges( )==2);
        should(graph.nodeFromId(3).numberOfEdges( )==2);

        // check edges
        should(graph.nodeFromId(0).hasEdgeId(0));
        should(graph.nodeFromId(0).hasEdgeId(2));

        should(graph.nodeFromId(1).hasEdgeId(0));
        should(graph.nodeFromId(1).hasEdgeId(graph.reprEdgeId(3)));
        should(graph.nodeFromId(1).hasEdgeId(graph.reprEdgeId(4)));

        should(graph.nodeFromId(2).hasEdgeId(1));
        should(graph.nodeFromId(2).hasEdgeId(2));

        should(graph.nodeFromId(3).hasEdgeId(1));
        should(graph.nodeFromId(3).hasEdgeId(graph.reprEdgeId(3)));
        should(graph.nodeFromId(3).hasEdgeId(graph.reprEdgeId(4)));

        // check which edge is the deleted
        IdType deletedEdge = graph.reprEdgeId(3)==3 ? 4 : 3;
        should(!graph.hasEdgeId(deletedEdge));

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
        std::vector<IdType> activeEdgeVec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        should(activeEdgeVec.size() == 3);
        should(activeEdgeVec[0] == 0);
        should(activeEdgeVec[1] == 1);
        should(activeEdgeVec[2] == 2);
        }

        std::vector<IdType> activeNodes(graph.nodeIdsBegin(),graph.nodeIdsEnd());
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



        should(graph.reprEdgeId(0)==0);
        should(graph.reprEdgeId(1)==1);
        should(graph.reprEdgeId(2)==2);



        // merge edge 0 
        graph.mergeRegions(0);
        should(graph.numberOfNodes() == 2);
        should(graph.numberOfEdges() == 1);
        {
        std::vector<IdType> activeEdgeVec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        should(activeEdgeVec.size() == 1);
        }

        should(graph.reprNodeId(0)==graph.reprNodeId(1));
        should(graph.reprNodeId(0)!=graph.reprNodeId(2));
        should(graph.reprNodeId(1)!=graph.reprNodeId(2));

        should(graph.reprEdgeId(1)==graph.reprEdgeId(2));
        should(!graph.reprEdgeId(0)==graph.reprEdgeId(1));

        should(graph.reprEdgeId(0)==0);
        should(graph.reprEdgeId(1)!=0);
        should(graph.reprEdgeId(2)!=0);

        should(graph.reprEdgeId(0)!=graph.reprEdgeId(2));
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
            should(graph.hasEdgeId(e));
            // fist node is rep of e, second node still untouched e+1
            should(graph.edgeFromId(e)[0]==graph.reprNodeId(e));
            should(graph.edgeFromId(e)[1]==e+1);

            // remove the edge
            graph.mergeRegions(e);

            should(!graph.hasEdgeId(e));
            should(graph.numberOfEdges()==nChainEdges-e-1);
            should(graph.numberOfNodes()==nChainNodes-e-1);
        }
    }

    void gridTest()
    {
        // 3x3
        //
        //   0 | 1 | 2 
        //   _   _   _
        //   3 | 4 | 5
        //   _   _   _
        //   6 | 7 | 8
        Lset nodeSet;
        Lvec nodeVec;
        Lset edgeSet;
        Lvec edgeVec;
        EVec explicitEdges;

        const size_t nChainNodes  = 9;
        const size_t nChainEdges  = 12;
        // 0 - 1 - 2 - .. nChainNodes-1
        MergeGraphType graph(nChainNodes,nChainEdges);
        graph.setInitalEdge(0, 0,1);  const size_t e01=0;
        graph.setInitalEdge(1, 1,2);  const size_t e12=1;

        graph.setInitalEdge(2, 3,4);  const size_t e34=2;
        graph.setInitalEdge(3, 4,5);  const size_t e45=3;

        graph.setInitalEdge(4, 6,7);  const size_t e67=4;
        graph.setInitalEdge(5, 7,8);  const size_t e78=5;

        graph.setInitalEdge(6, 0,3);  const size_t e03=6;
        graph.setInitalEdge(7, 1,4);  const size_t e14=7;
        graph.setInitalEdge(8, 2,5);  const size_t e25=8;

        graph.setInitalEdge(9 , 3,6); const size_t e36=9;
        graph.setInitalEdge(10, 4,7); const size_t e47=10;
        graph.setInitalEdge(11, 5,8); const size_t e58=11;

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
        //
        //   0 | 1 | 2 
        //   _   _   _
        //   3   4 | 5
        //   _   _   _
        //   6 | 7 | 8
        graph.mergeRegions(e34);
        should(graph.numberOfNodes()==8);
        should(graph.numberOfEdges()==11);
        graph.stateOfInitalEdges(edgeStateCheck,edgeStateCheck+12);
        edgeStateTrue[e34]=0;
        shouldEqualSequence(edgeStateTrue,edgeStateTrue+12,edgeStateCheck);
        should(graph.reprNodeId(3)==graph.reprNodeId(4));
        const IdType rep34 = graph.reprNodeId(3);
        const IdType del34 = (rep34 == 3 ? 4: 3);
        const Node & n34=graph.nodeFromId(graph.reprNodeId(3));
        should(graph.hasNodeId(rep34));
        should(!graph.hasNodeId(del34));
        should(!graph.hasEdgeId(e34));
        should(n34.numberOfEdges()==5);
        should(n34.hasEdgeId(graph.reprEdgeId(e03)));  // 0-3
        should(n34.hasEdgeId(graph.reprEdgeId(e14)));  // 1-4
        should(n34.hasEdgeId(graph.reprEdgeId(e36)));  // 3-6
        should(n34.hasEdgeId(graph.reprEdgeId(e47)));  // 4-7
        should(n34.hasEdgeId(graph.reprEdgeId(e45)));  // 4-5
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==8);
        should(nodeVec.size()==8);
        should(nodeSet.find(rep34)!=nodeSet.end());
        should(nodeSet.find(del34)==nodeSet.end());
        should(nodeSet.find(0)!=nodeSet.end());
        should(nodeSet.find(1)!=nodeSet.end());
        should(nodeSet.find(2)!=nodeSet.end());
        should(nodeSet.find(5)!=nodeSet.end());
        should(nodeSet.find(6)!=nodeSet.end());
        should(nodeSet.find(7)!=nodeSet.end());
        // check representative edges
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==11);
        should(edgeVec.size()==11);
        should(edgeSet.find(e34)==edgeSet.end());
        should(edgeSet.find(e01)!=edgeSet.end());
        should(edgeSet.find(e12)!=edgeSet.end());
        should(edgeSet.find(e03)!=edgeSet.end());
        should(edgeSet.find(e14)!=edgeSet.end());
        should(edgeSet.find(e25)!=edgeSet.end());
        should(edgeSet.find(e45)!=edgeSet.end());
        should(edgeSet.find(e36)!=edgeSet.end());
        should(edgeSet.find(e47)!=edgeSet.end());
        should(edgeSet.find(e58)!=edgeSet.end());
        should(edgeSet.find(e67)!=edgeSet.end());
        should(edgeSet.find(e78)!=edgeSet.end());

        // merge edge between 6-7
        // this will reduce the number of active edges by 2:
        // edge e67 will dissaper and 3-6 4-7 will merge
        //
        //   0 | 1 | 2 
        //   _   _   _
        //   3   4 | 5
        //   _   _   _
        //   6   7 | 8
        graph.mergeRegions(e67);
        should(graph.numberOfNodes()==7);
        should(graph.numberOfEdges()==9);
        graph.stateOfInitalEdges(edgeStateCheck,edgeStateCheck+12);
        edgeStateTrue[e67]=0;
        shouldEqualSequence(edgeStateTrue,edgeStateTrue+12,edgeStateCheck);
        should(graph.reprNodeId(6)==graph.reprNodeId(7));
        const IdType rep67 = graph.reprNodeId(6);
        const IdType del67 = (rep67 == 6 ? 7: 6);
        const Node & n67=graph.nodeFromId(rep67);
        should(graph.reprEdgeId(e36)==graph.reprEdgeId(e47));
        should(n67.numberOfEdges()==2);
        should(n67.hasEdgeId(graph.reprEdgeId(e36)));  // 0-3
        should(n67.hasEdgeId(graph.reprEdgeId(e47)));  // 1-4
        should(n67.hasEdgeId(graph.reprEdgeId(e78)));  // 3-6
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==7);
        should(nodeVec.size()==7);
        should(nodeSet.find(rep67)!=nodeSet.end());
        should(nodeSet.find(del67)==nodeSet.end());
        should(graph.hasNodeId(rep67));
        should(!graph.hasNodeId(del67));
        // check representative edges
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());




        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==9);
        should(edgeVec.size()==9);
        const size_t rep36_47 = graph.reprEdgeId(e36);
        const size_t del36_47 = (rep36_47 == e36 ? e47 : e36);
        should(rep36_47 == e36 || rep36_47 == e47 );
        should(rep36_47 != del36_47);
        should(graph.hasEdgeId(rep36_47));
        should(!graph.hasEdgeId(del36_47));
        should(edgeSet.find(rep36_47)!=edgeSet.end());
        should(edgeSet.find(del36_47)==edgeSet.end());


        // merge edge between 3-6  (and 4-7)
        // this will reduce the number of active edges by 1:
        // edge  3-6 4-7 will dissapear
        //
        //   0 | 1 | 2 
        //   _   _   _
        //   3   4 | 5
        //           _
        //   6   7 | 8
        graph.mergeRegions(graph.reprEdgeId(e36));
        should(graph.numberOfNodes()==6);
        should(graph.numberOfEdges()==8);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==6);
        should(nodeVec.size()==6);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        explicitEdges = EVec(graph.edgesBegin(),graph.edgesEnd());



        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==8);
        should(explicitEdges.size()==8);
        should(edgeVec.size()==8);

        {
            EdgeIdIt tmp     = graph.edgeIdsBegin();
            ++ tmp;
            should(tmp!=graph.edgeIdsBegin());
            should(tmp!=lemon::INVALID);
            -- tmp;
            should(tmp==graph.edgeIdsBegin());


            should(graph.edgeIdsBegin()!=lemon::INVALID);
            should(graph.edgeIdsEnd()  ==lemon::INVALID);

            should(graph.edgesBegin()!=lemon::INVALID);
            should(graph.edgesEnd()==lemon::INVALID);

            should(graph.nodeIdsBegin()!=lemon::INVALID);
            should(graph.nodeIdsEnd()  ==lemon::INVALID);

            should(graph.nodesBegin()!=lemon::INVALID);
            should(graph.nodesEnd()==lemon::INVALID);
        }


        // merge edge between 5-8  
        // this will reduce the number of active edges by 2:
        // edge  5-8 will dissapear and 4-5 7-8 will be merged
        //
        //   0 | 1 | 2 
        //   _   _   _
        //   3   4 | 5
        //            
        //   6   7 | 8
        graph.mergeRegions(graph.reprEdgeId(e58));
        should(graph.numberOfNodes()==5);
        should(graph.numberOfEdges()==6);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==5);
        should(nodeVec.size()==5);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==6);
        should(edgeVec.size()==6);


        // merge edge between 0-1 
        // this will reduce the number of active edges by 2:
        // edge  0-1 will dissapear and 0-3 1-4 will be merged
        //
        //   0   1 | 2 
        //   _   _   _
        //   3   4 | 5
        //            
        //   6   7 | 8
        graph.mergeRegions(graph.reprEdgeId(e01));
        should(graph.numberOfNodes()==4);
        should(graph.numberOfEdges()==4);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==4);
        should(nodeVec.size()==4);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==4);
        should(edgeVec.size()==4);


        // merge edge between 1-2
        // this will reduce the number of active edges by 1:
        // edge  1-2 will dissaper
        //   0   1   2 
        //   _   _   _
        //   3   4 | 5
        //            
        //   6   7 | 8
        graph.mergeRegions(graph.reprEdgeId(e12));
        should(graph.numberOfNodes()==3);
        should(graph.numberOfEdges()==3);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==3);
        should(nodeVec.size()==3);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==3);
        should(edgeVec.size()==3);


        // merge edge between 4-5
        // this will reduce the number of active edges by 2:
        // edge 4-5 will dissaper and the rest of edges
        // will merge to a single edeg
        //   0   1   2 
        //   _   _   _
        //   3   4   5
        //            
        //   6   7   8
        graph.mergeRegions(graph.reprEdgeId(e45));
        should(graph.numberOfNodes()==2);
        should(graph.numberOfEdges()==1);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==2);
        should(nodeVec.size()==2);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());
        should(edgeSet.size()==1);
        should(edgeVec.size()==1);


        // merge edge between 0-3
        // this will reduce the number of active edges by 1:
        //
        //   0   1   2 
        //           
        //   3   4   5
        //            
        //   6   7   8
        graph.mergeRegions(graph.reprEdgeId(e03));
        should(graph.numberOfNodes()==1);
        should(graph.numberOfEdges()==0);
        // check representatives nodes
        nodeSet=Lset(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        nodeVec=Lvec(graph.nodeIdsBegin(),graph.nodeIdsEnd());
        shouldEqualSequence(nodeSet.begin(),nodeSet.end(),nodeVec.begin());
        should(nodeSet.size()==1);
        should(nodeVec.size()==1);
        edgeSet=Lset(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        edgeVec=Lvec(graph.edgeIdsBegin(),graph.edgeIdsEnd());
        shouldEqualSequence(edgeSet.begin(),edgeSet.end(),edgeVec.begin());

        should(edgeSet.size()==0);
        should(edgeVec.size()==0);

    }

    void accumulatorChaiMapTest()
    {
        MergeGraphType graph(4,4);
        // 2x2 grid graph
        //  0 | 1
        //  _   _ 
        //  2 | 3

        // set inital edges (with one double edge)
        graph.setInitalEdge(0,0,1);
        graph.setInitalEdge(1,2,3);
        graph.setInitalEdge(2,0,2);
        graph.setInitalEdge(3,1,3);
        should(graph.numberOfNodes() == 4);
        should(graph.numberOfEdges() == 4);
        typedef  float MapValueType;
        typedef vigra::acc::Select< vigra::acc::DataArg<1>, vigra::acc::LabelArg<2>,vigra::acc::Sum, vigra::acc::Mean ,vigra::acc::Count> AccChain;
        typedef vigra::AccumulatorChainMap<MergeGraphType,2,MapValueType,AccChain > AccChainMapType;

        typedef vigra::MultiArray<2,IdType> LabelArrayType;
        typedef vigra::MultiArray<2,MapValueType>     ValueTypeArrayType;

        LabelArrayType     labels(typename LabelArrayType::difference_type(2,2));
        ValueTypeArrayType data(typename ValueTypeArrayType::difference_type(2,2));


        // view map 

        typedef AccumulatorChainMapTagView<AccChainMapType,vigra::acc::Sum> SumView;



        labels(0,0)=0;
        labels(1,0)=1;
        labels(0,1)=2;
        labels(1,1)=3;


        data(0,0)=1;
        data(1,0)=2;
        data(0,1)=3;
        data(1,1)=4;

        AccChainMapType nodeMap(graph,data,labels);
        SumView sumView(nodeMap);

        // register callbacks
        graph.registerMergeNodeCallBack(nodeMap,& AccChainMapType::merge);


        //graph.mergeRegions(0);
        should(vigra::acc::get<vigra::acc::Sum>(nodeMap.accChainArray(),graph.reprNodeId(0))==1);
        should(vigra::acc::get<vigra::acc::Sum>(nodeMap.accChainArray(),graph.reprNodeId(1))==2);



        should(sumView[graph.reprNodeId(0)]==1);
        should(sumView[graph.reprNodeId(1)]==2);
        graph.mergeRegions(0);
        should(vigra::acc::get<vigra::acc::Sum>(nodeMap.accChainArray(),graph.reprNodeId(0))==3);
        should(nodeMap.get<vigra::acc::Sum>(graph.reprNodeId(0))==3);
        should(sumView[graph.reprNodeId(0)]==3);
    }

};


template<class ID_TYPE>
struct PartitonTest
{
    typedef ID_TYPE IdType;
    typedef vigra::detail_merge_graph::IterablePartition<IdType> PartitionType;
    typedef std::set<IdType> SetType;
    typedef std::vector<IdType> VecType;
    PartitonTest()
    {

    }

    void trueReps(const PartitionType ufd,SetType &set){
        set.clear();
        for(IdType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
    }

    void trueReps(const PartitionType ufd,SetType &set,SetType & c){
        set.clear();
        for(IdType i=0;i<ufd.numberOfElements();++i){
            if(ufd.find(i)==i){
                set.insert(i);
            }
        }
        for(typename  SetType::const_iterator iter=c.begin();iter!=c.end();++iter){
            const IdType toRemove=*iter;
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

        IdType rep12 = ufd.find(1);
        erased.insert(rep12);
        ufd.eraseElement(rep12);
        trueReps(ufd,trueRep,erased);
        testReps(ufd,testRep);
        shouldEqualSequence(trueRep.begin(),trueRep.end(),testRep.begin());
    }
    
};


struct MinIndexedPQTest
{
    typedef vigra::MinIndexedPQ<float> PqType;

    MinIndexedPQTest()
    {

    }


    void constructorTest(){
        PqType pq(10);
        should(pq.size()==0);
        for(size_t i=0;i<10;++i){
            should(!pq.contains(i));
        }
    }

    
};
        
struct MergeGraphTestSuite
: public vigra::test_suite
{
    MergeGraphTestSuite()
    : vigra::test_suite("MergeGraphTestSuite")
    {   

        add( testCase( &PartitonTest<vigra::Int32>::iteratorTest1));
        add( testCase( &PartitonTest<vigra::Int32>::iteratorTest2));

        add( testCase( &MergeGraphTest<vigra::Int32>::mergeSimpleDoubleEdgeTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::mergeTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::chainTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::gridTest));
        add( testCase( &MergeGraphTest<vigra::Int32>::accumulatorChaiMapTest));


        add( testCase( &MinIndexedPQTest::constructorTest));
    }
};

int main(int argc, char ** argv)
{
    MergeGraphTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
