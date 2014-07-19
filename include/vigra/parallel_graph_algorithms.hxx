/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
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

/**
 * This header provides definitions of graph-related algorithms
 */

#ifndef VIGRA_PARALLEL_GRAPH_ALGORITHMS_HXX
#define VIGRA_PARALLEL_GRAPH_ALGORITHMS_HXX

/*std*/
#include <algorithm>
#include <vector>
#include <functional>
#include <set>


/*vigra*/
#include "graphs.hxx"
#include "graph_generalization.hxx"
#include "multi_gridgraph.hxx"
#include "priority_queue.hxx"
#include "union_find.hxx"
#include "adjacency_list_graph.hxx"
#include "graph_maps.hxx"
#include "timing.hxx"
#include "thread_pool.hxx"



#include <boost/asio.hpp>
#include <boost/thread.hpp>

namespace vigra{


struct GraphEdgeLocking
{
    GraphEdgeLocking(const size_t maxEdgeId)
    : edgeLock_(maxEdgeId+1)
    {

    }

    std::vector< boost::mutex > edgeLock_;
};



template<class RAG, class GRAPH, class AFF, class LABELS>
struct FillAffEdges{
    typedef typename  GRAPH::EdgeIt EdgeIt;
    typedef typename  GRAPH::Edge Edge;
    FillAffEdges(
        const RAG & rag,
        const GRAPH & graph,
        AFF & affEdges,
        const LABELS & labels,
        const size_t rangeBegin, 
        const size_t rangeEnd,
        GraphEdgeLocking & edgeLocking
    )
    :
    rag_(&rag),
    graph_(&graph),
    affEdges_(&affEdges),
    labels_(&labels),
    rangeBegin_(rangeBegin),
    rangeEnd_(rangeEnd),
    edgeLocking_(&edgeLocking)
    {
        std::cout<<"my edge range "<<rangeBegin_<<" "<<rangeEnd_<<"\n";    
    }   

    void operator()(){


        std::cout<<"START edge range "<<rangeBegin_<<" "<<rangeEnd_<<"\n";    

        std::cerr<<"get first iter \n";
        EdgeIt start(*graph_,graph_->edgeFromId(rangeBegin_));
        std::cerr<<"get second iter \n";
        EdgeIt end(*graph_,graph_->edgeFromId(rangeEnd_));


        std::cerr<<"find start....\n";
        while(start==lemon::INVALID){
            ++rangeBegin_;
            start = EdgeIt(*graph_,graph_->edgeFromId(rangeBegin_));
        }

        if(rangeEnd_ == graph_->maxEdgeId()+1){
            end = GraphIteratorAccessor<GRAPH>::edgesEnd(*graph_);
        }

        std::cerr<<"find end.....\n";
        while(end==lemon::INVALID){
            ++rangeEnd_;
            end = EdgeIt(*graph_,graph_->edgeFromId(rangeEnd_));
        }

        std::cerr<<"REAL START edge range "<<rangeBegin_<<" "<<rangeEnd_<<"\n";    

        while(start!=end){

            const Edge edge(*start);
            const int lu = (*labels_)[graph_->u(edge)];
            const int lv = (*labels_)[graph_->v(edge)];
            //std::cout<<"edge between ?? "<<lu<<" "<<lv<<"\n";
            if(  lu!=lv ){
                //std::cout<<"foundEdge\n";
                typename RAG::Edge ragEdge= rag_->findEdge(rag_->nodeFromId(lu),rag_->nodeFromId(lv));

                edgeLocking_->edgeLock_[rag_->id(ragEdge)].lock();
                //std::cout<<"edge "<<rag_->id(ragEdge)<<" locked...\n";
                (*affEdges_)[ragEdge].push_back(edge);
                edgeLocking_->edgeLock_[rag_->id(ragEdge)].unlock();
            }

            ++start;
        }
    }

    const RAG * rag_;
    const GRAPH * graph_;
    AFF * affEdges_;
    const LABELS * labels_;
    size_t rangeBegin_;
    size_t rangeEnd_;
    GraphEdgeLocking * edgeLocking_;
};



    /// \brief make a region adjacency graph from a graph and labels w.r.t. that graph
    ///
    /// \param graphIn  : input graph
    /// \param labels   : labels w.r.t. graphIn
    /// \param[out] rag  : region adjacency graph 
    /// \param[out] affiliatedEdges : a vector of edges of graphIn for each edge in rag
    /// \param      ingoreLabel : optional label to ignore (default: -1 means no label will be ignored)
    ///
    template<
        class GRAPH_IN,
        class GRAPH_IN_NODE_LABEL_MAP
    >
    void makeRegionAdjacencyGraphFast(
        GRAPH_IN                   graphIn,
        GRAPH_IN_NODE_LABEL_MAP    labels,
        AdjacencyListGraph & rag,
        typename AdjacencyListGraph:: template EdgeMap< std::vector<typename GRAPH_IN::Edge> > & affiliatedEdges,
        const Int64   maxLabel,
        const Int64   reserveEdges
    ){
        rag=AdjacencyListGraph();
        typedef typename GraphMapTypeTraits<GRAPH_IN_NODE_LABEL_MAP>::Value LabelType;
        typedef GRAPH_IN GraphIn;
        typedef AdjacencyListGraph GraphOut;

        typedef typename GraphIn::Edge   EdgeGraphIn;
        //typedef typename GraphIn::NodeIt NodeItGraphIn;
        typedef typename GraphIn::EdgeIt EdgeItGraphIn;
        //typedef typename GraphOut::Edge   EdgeGraphOut; 

        std::cout<<"THREADED VERSION  lalalala\n";

        USETICTOC;
        std::cout<<"add nodes\n";
        TIC;
        rag.assignNodeRange(1,maxLabel+1);
        TOC;

        std::cout<<"add edges\n";
        TIC
        for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
            const EdgeGraphIn edge(*e);
            const LabelType lu = labels[graphIn.u(edge)];
            const LabelType lv = labels[graphIn.v(edge)];
            if(  lu!=lv  ){
                // if there is an edge between lu and lv no new edge will be added
                rag.addEdge( rag.nodeFromId(lu),rag.nodeFromId(lv));
            }
        }
        TOC


        std::cout<<"setup hyper edges\n";
        TIC




        // SET UP HYPEREDGES
        affiliatedEdges.assign(rag);


        

        typedef typename AdjacencyListGraph:: template EdgeMap< std::vector<typename GRAPH_IN::Edge> > AffType;
        typedef FillAffEdges<AdjacencyListGraph, GRAPH_IN, AffType,GRAPH_IN_NODE_LABEL_MAP> FillAffEdgesJob;



        {
            GraphEdgeLocking edgeLocking(rag.maxEdgeId());
            ThreadPool threadPool(20);
            const size_t perRun=graphIn.maxEdgeId()/100;
            size_t counting=0;
            EdgeItGraphIn rangeBegin(graphIn);
            EdgeItGraphIn rangeEnd = GraphIteratorAccessor<GRAPH_IN>::edgesEnd(graphIn);

            for(size_t id=0; id<=graphIn.maxEdgeId(); id+=perRun){

                size_t endId = id+perRun;
                endId = std::min(size_t(graphIn.maxEdgeId()+1),endId);


                threadPool.run_task(FillAffEdgesJob(rag, graphIn, affiliatedEdges, labels, id, endId, edgeLocking));
            }   
        }

        TOC;

       //// add edges
       //for(EdgeItGraphIn e(graphIn);e!=lemon::INVALID;++e){
       //    const EdgeGraphIn edge(*e);
       //    const LabelType lu = labels[graphIn.u(edge)];
       //    const LabelType lv = labels[graphIn.v(edge)];
       //    //std::cout<<"edge between ?? "<<lu<<" "<<lv<<"\n";
       //    if(  lu!=lv ){
       //        //std::cout<<"find edge between "<<lu<<" "<<lv<<"\n";
       //        EdgeGraphOut ragEdge= rag.findEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
       //        //std::cout<<"invalid?"<<bool(ragEdge==lemon::INVALID)<<" id "<<rag.id(ragEdge)<<"\n";
       //        affiliatedEdges[ragEdge].push_back(edge);
       //        //std::cout<<"write done\n";
       //    }
       //}

        //#pragma omp parallel for
        //for(size_t id=0; id<=graphIn.maxEdgeId(); ++id){
        //    const EdgeGraphIn edge(graphIn.edgeFromId(id));
        //    if(edge!=lemon::INVALID){
        //        const LabelType lu = labels[graphIn.u(edge)];
        //        const LabelType lv = labels[graphIn.v(edge)];
        //        //std::cout<<"edge between ?? "<<lu<<" "<<lv<<"\n";
        //        if(  lu!=lv ){
        //            
        //            #pragma omp critical
        //            {
        //            EdgeGraphOut ragEdge= rag.findEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
        //            affiliatedEdges[ragEdge].push_back(edge);
        //            }
        //            //std::cout<<"write done\n";
        //        }
        //    }
        //}
        //TOC
    }



} // namespace vigra

#endif // VIGRA_PARALLEL_GRAPH_ALGORITHMS_HXX
