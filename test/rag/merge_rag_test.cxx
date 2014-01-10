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
#include "vigra/merge_graph/new_merge_graph.hxx"
using namespace vigra;




template<class IN_LABEL_TYPE>
struct Rag2MergeGraphTest{

    typedef IN_LABEL_TYPE                        InLabelType;
    typedef vigra::Rag<2,InLabelType>            RagType;
    typedef vigra::MergeGraphAdaptor<RagType>    MergeGraphType;
    typedef typename MergeGraphType::index_type  index_type;
    typedef typename MergeGraphType::IdType      IdType;
    typedef typename RagType::Node               RagNode;
    typedef typename RagType::Edge               RagEdge;
    typedef typename RagType::Arc                RagArc;
    typedef typename RagType::EdgeIt             RagEdgeIt;
    typedef typename RagType::NodeIt             RagNodeIt;
    typedef typename RagType::ArcIt              RagArcIt;
    typedef typename RagType::IncEdgeIt          RagIncEdgeIt;
    typedef typename RagType::InArcIt            RagInArcIt;
    typedef typename RagType::OutArcIt           RagOutArcIt;


    typedef typename MergeGraphType::Node        Node;
    typedef typename MergeGraphType::Edge        Edge;
    typedef typename MergeGraphType::Arc         Arc;
    typedef typename MergeGraphType::EdgeIt      EdgeIt;
    typedef typename MergeGraphType::NodeIt      NodeIt;
    typedef typename MergeGraphType::ArcIt       ArcIt;
    typedef typename MergeGraphType::IncEdgeIt   IncEdgeIt;
    typedef typename MergeGraphType::InArcIt     InArcIt;
    typedef typename MergeGraphType::OutArcIt    OutArcIt;

    typedef typename RagType::NeighborNodeIt     NeighborNodeIt;
    typedef typename RagType::InputLabelingView  InputLabelingView;
    typedef typename RagType::InputLabelingArray InputLabelingArray;


    RagType rag2x2_;

    bool nodeHasEdge(const MergeGraphType & g,const index_type n,const index_type & e)const{
        const Edge edge  = g.edgeFromId(e);
        const IdType u   = g.id(g.u(edge));
        const IdType v   = g.id(g.v(edge));
        return n==u || n==v;
    }

    bool nodeHasEdge(const MergeGraphType & g,const Node & node,const index_type & e)const{
        const Edge edge  = g.edgeFromId(e);
        const IdType u   = g.id(g.u(edge));
        const IdType v   = g.id(g.v(edge));
        const IdType n   = g.id(node);
        return n==u || n==v;
    }
    bool edgeHasNode(const MergeGraphType & g, Edge edge, Node node){
        if(  g.id(g.u(edge)) == g.id(node) ){
            return true;
        }
        else if(g.id(g.v(edge)) == g.id(node) ){
            return true;
        }
        else{
            return false;
        }
    }

    bool edgeHasNode(const MergeGraphType & g, Edge edge, IdType node){
        if(  g.id(g.u(edge)) == node ){
            return true;
        }
        else if(g.id(g.v(edge)) == node){
            return true;
        }
        else{
            return false;
        }
    }




    size_t nodeDegree(const MergeGraphType & g ,index_type n)const{
        return g.degree(g.nodeFromId(n));
    }
    size_t nodeDegree(const MergeGraphType & g ,const Node & n)const{
        return g.degree(n);
    }


    Rag2MergeGraphTest(){
        //std::cout<<"Rag2MergeGraphTest \n";
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
        rag2x2_=g;
    }



    void ragTest()
    {
        //std::cout<<"ragTest \n";
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

        // create merge graph adpator
        MergeGraphType g(rag);



        //std::cout<<"basic sizes \n";
        // assert basic sizes

        should(g.edgeNum()==4);
        should(g.nodeNum()==4);
        should(g.arcNum()==8);

        //std::cout<<"max  ids \n";

        should(g.maxEdgeId()==4);
        should(g.maxNodeId()==4);
        should(g.maxArcId()==g.maxEdgeId()*2+1);

        //std::cout<<"find edges \n";

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

        std::cout<<"fill node vec \n";
        std::vector<Node> allNodes(nbegin,nend);
        should(allNodes.size()==4);


        std::cout<<"fill edge vec \n";
        std::vector<Edge> allEdges(ebegin,eend);
        should(allEdges.size()==4);


        should(1==g.id( allEdges[0] ) );
        should(2==g.id( allEdges[1] ) );
        should(3==g.id( allEdges[2] ) );
        should(4==g.id( allEdges[3] ) );
        
     
        std::vector<Arc>  allArcs( abegin,aend);
        should(allArcs.size() ==8);
          
        should(1==g.id( allArcs[0] ) );
        should(2==g.id( allArcs[1] ) );
        should(3==g.id( allArcs[2] ) );
        should(4==g.id( allArcs[3] ) );

        should(1+g.maxEdgeId()+1==g.id( allArcs[4] ) );
        should(2+g.maxEdgeId()+1==g.id( allArcs[5] ) );
        should(3+g.maxEdgeId()+1==g.id( allArcs[6] ) );
        should(4+g.maxEdgeId()+1==g.id( allArcs[7] ) );

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
    void ragArcTest()
    {
        //std::cout<<"ragArcTest \n";
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
        MergeGraphType g(rag);
        
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

    void ragEdgeItTest()
    {
        //std::cout<<"ragEdgeItTest \n";
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
        MergeGraphType g(rag);

        

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
        //std::cout<<"ragNodeItTest \n";
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
        MergeGraphType g(rag);
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
        //std::cout<<"ragArcItTest \n";
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
        MergeGraphType g(rag);
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

            shouldEqual(1+g.maxEdgeId()+1,g.id(arcVec[4]));
            shouldEqual(2+g.maxEdgeId()+1,g.id(arcVec[5]));
            shouldEqual(3+g.maxEdgeId()+1,g.id(arcVec[6]));
            shouldEqual(4+g.maxEdgeId()+1,g.id(arcVec[7]));
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

            shouldEqual(1+g.maxEdgeId()+1,g.id(arcVec[4]));
            shouldEqual(2+g.maxEdgeId()+1,g.id(arcVec[5]));
            shouldEqual(3+g.maxEdgeId()+1,g.id(arcVec[6]));
            shouldEqual(4+g.maxEdgeId()+1,g.id(arcVec[7]));
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

            shouldEqual(1+g.maxEdgeId()+1,g.id(arcVec[3]));
            shouldEqual(2+g.maxEdgeId()+1,g.id(arcVec[4]));
            shouldEqual(3+g.maxEdgeId()+1,g.id(arcVec[5]));
            shouldEqual(4+g.maxEdgeId()+1,g.id(arcVec[6]));
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

        //std::cout<<"ragFindEdgeTest \n";

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
        MergeGraphType g(rag);

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

        //std::cout<<"ragUVOrderTest \n";
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
        MergeGraphType g(rag);

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
        //std::cout<<"ragIncEdgeItTest \n";
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
        MergeGraphType g(rag);

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
        //std::cout<<"ragInArcItTest \n";
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
        MergeGraphType g(rag);

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
        //std::cout<<"ragOutArcItTest \n";
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
        MergeGraphType g(rag);

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


    void ragMergeChainTest()
    {

        const size_t nChainNodes  = 6;
        const size_t nChainEdges  = nChainNodes-1;

        InputLabelingArray labels(typename InputLabelingArray::difference_type(nChainNodes,1));

        for(index_type i=0;i<nChainNodes;++i){
            labels(i)=i+1;
        }
        RagType rag(labels);

        
        should(rag.edgeFromId(0)==lemon::INVALID);
        should(rag.edgeFromId(6)==lemon::INVALID);
        should(rag.edgeFromId(1)!=lemon::INVALID);
        should(rag.edgeFromId(2)!=lemon::INVALID);
        should(rag.edgeFromId(3)!=lemon::INVALID);
        should(rag.edgeFromId(4)!=lemon::INVALID);
        should(rag.edgeFromId(5)!=lemon::INVALID);





        MergeGraphType g(rag);
        should(g.edgeFromId(0)==lemon::INVALID);
        should(g.edgeFromId(6)==lemon::INVALID);
        should(g.edgeFromId(1)!=lemon::INVALID);
        should(g.edgeFromId(2)!=lemon::INVALID);
        should(g.edgeFromId(3)!=lemon::INVALID);
        should(g.edgeFromId(4)!=lemon::INVALID);
        should(g.edgeFromId(5)!=lemon::INVALID);

        should(g.nodeFromId(0)==lemon::INVALID);
        should(g.nodeFromId(7)==lemon::INVALID);

        should(g.nodeFromId(1)!=lemon::INVALID);
        should(g.nodeFromId(2)!=lemon::INVALID);
        should(g.nodeFromId(3)!=lemon::INVALID);
        should(g.nodeFromId(4)!=lemon::INVALID);
        should(g.nodeFromId(5)!=lemon::INVALID);
        should(g.nodeFromId(6)!=lemon::INVALID);



        should(g.nodeNum() == nChainNodes);
        should(g.edgeNum() == nChainEdges);


        // remove edges from 1 to 5


        for(size_t e=1;e<=5;++e){



            should(g.edgeNum()==nChainEdges-(e-1));
            should(g.nodeNum()==nChainNodes-(e-1));
            // check that edge is there 
            should(g.hasEdgeId(e));
            should(g.edgeFromId(e)!=lemon::INVALID);


            // fist node is rep of e, second node still untouched e+1
            should(g.id(g.u(g.edgeFromId(e)))==g.reprNodeId(e));


            should(g.id(g.v(g.edgeFromId(e)))==e+1);

            // remove the edge

            g.mergeRegions(e);

            should(!g.hasEdgeId(e));
            should(g.edgeFromId(e)==lemon::INVALID);
            should(g.edgeNum()==nChainEdges-(e-1)-1);
            should(g.nodeNum()==nChainNodes-(e-1)-1);
        }
    }
};


 
struct RagMergeGraphAdaptorTestSuite
: public vigra::test_suite
{
    RagMergeGraphAdaptorTestSuite()
    : vigra::test_suite("RagMergeGraphTestSuite")
    {   
        
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragArcTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragFindEdgeTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragUVOrderTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragEdgeItTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragNodeItTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragArcItTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragIncEdgeItTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragInArcItTest));
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragOutArcItTest));


        // test which ragMergeChainTest
        add( testCase( &Rag2MergeGraphTest<vigra::UInt32>::ragMergeChainTest));

    }
};

int main(int argc, char ** argv)
{
    RagMergeGraphAdaptorTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

