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

#ifndef VIGRA_GRAPH_RAG_HXX
#define VIGRA_GRAPH_RAG_HXX

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

    template<class GRAPH>
    struct MakeEdgeRangeVector{

        typedef typename  GRAPH::Edge Edge;
        typedef std::pair<Edge, Edge> EdgePair;
        typedef std::vector<EdgePair> EdgeRangeVector;

    };

    template<class GRAPH>
    std::pair< typename GRAPH::Edge, typename GRAPH::index_type>   nexValidEdge(
        const GRAPH & g,
        typename GRAPH::index_type maybeValidId
    ){
        //std::cout<<"get next valid edge for "<<maybeValidId<<"\n";

        typedef typename GRAPH::Edge Edge;
        typedef typename GRAPH::index_type index_type;
        typedef std::pair< Edge, index_type> EdgeIndexPair;

        Edge edge = g.edgeFromId(maybeValidId);

        while(edge==lemon::INVALID && maybeValidId< g.maxEdgeId()){
            ++maybeValidId;
            //std::cout<<"   try"<<maybeValidId<<"\n";
            edge = g.edgeFromId(maybeValidId);
        }
        if(edge==lemon::INVALID || maybeValidId>g.maxEdgeId()){
            //std::cout<<"   no success"<<maybeValidId<<"\n";
            return EdgeIndexPair(Edge(lemon::INVALID), -1);
        }
        else{
            //std::cout<<"   success"<<maybeValidId<<"\n";
            return EdgeIndexPair(edge,maybeValidId);
        }
    }


    template<class GRAPH>
    void edgeRanges(
        const GRAPH & g,
        typename MakeEdgeRangeVector<GRAPH>::EdgeRangeVector & edgeRanges,
        const size_t numberOfRangesHint
    );

    // for grid graph
    template<unsigned int DIM, class DTAG>
    void edgeRanges(
        const GridGraph<DIM, DTAG> & g,
        typename MakeEdgeRangeVector<GridGraph<DIM, DTAG> >::EdgeRangeVector & edgeRanges,
        const size_t numberOfRangesHint
    ){
        typedef GridGraph<DIM, DTAG> Graph;
        typedef typename MakeEdgeRangeVector< Graph >::EdgePair EdgePair;
        //typedef typename MakeEdgeRangeVector<Graph >::EdgeRangeVector EdgeRangeVector;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::index_type index_type;
        typedef std::pair< Edge, index_type> EdgeIndexPair;

        edgeRanges.clear();
        edgeRanges.reserve(numberOfRangesHint);

        const index_type maxEdgeId = g.maxEdgeId();
        const index_type integralRange = maxEdgeId/numberOfRangesHint;

        index_type startId=0;


        while(startId<maxEdgeId){
            const index_type endId = startId+integralRange;
            const EdgeIndexPair edgeIndexPairStart = nexValidEdge(g, startId);
            const EdgeIndexPair edgeIndexPairEnd   = endId<g.maxEdgeId() ? 
                                                        nexValidEdge(g, endId) : 
                                                        EdgeIndexPair(Edge(lemon::INVALID),-1);
            startId+=integralRange;
            edgeRanges.push_back(EdgePair(edgeIndexPairStart.first, edgeIndexPairEnd.first));
        }
    } 




    struct RagOptions{

        RagOptions(
            const Int64 minLabel=-1,
            const Int64 maxLabel=-1,
            const Int64 isDense=-1,
            const Int64 ignoreLabel=-1,
            const Int64 reserveEdges=-1,
            const Int64 nThreads=-1,
            const bool parallel =false
        ) 
        :
        minLabel_(minLabel),
        maxLabel_(maxLabel),
        isDense_(isDense),
        ignoreLabel_(ignoreLabel),
        reserveEdges_(reserveEdges),
        nThreads_(nThreads),
        parallel_(parallel)
        {
        }

        Int64 minLabel_;
        Int64 maxLabel_;
        Int64 isDense_;
        Int64 ignoreLabel_;
        Int64 reserveEdges_;
        Int64 nThreads_;
        bool  parallel_;

    };


    template<class GRAPH>
    struct AffiliatedEdgesHelper{
        typedef typename GRAPH::Edge Edge;
        typedef std::vector<Edge> EdgeVec;
        typedef typename AdjacencyListGraph:: template EdgeMap< EdgeVec > AffiliatedEdgesType;
    };


    inline size_t ragDefaultNumThreads(){
        return 8;
    }



    namespace rag_detail{

    template<class GRAPH,class LABEL_MAP>
    void preprocess(
        const GRAPH & g, 
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & optionsIn,
        RagOptions & optionsOut
    ){
        // process the options
        std::cout<<"process the options \n";

        // check thread options
        if(optionsIn.nThreads_>1)
            optionsOut.parallel_=true;
        if(optionsIn.parallel_==true && optionsIn.nThreads_<1){
            // guess the number of threads
            optionsOut.nThreads_=static_cast<UInt64>(ragDefaultNumThreads());
        }


    }


    template<class GRAPH,class LABEL_MAP>
    void addNodesSerial(
        const GRAPH & g,  
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"addNodes Serial \n";

        typedef typename GRAPH::NodeIt NodeIt;
        typedef typename LABEL_MAP::Value LabelType;

        //typedef typename LABEL_MAP::Value LabelType;
        // dense, no ignore label, min and max are known
        if(opt.isDense_ && opt.minLabel_!=-1 && opt.maxLabel_ != -1 && opt.ignoreLabel_==-1){
            std::cout<<"addNodes batch \n";    
            rag.assignNodeRange(opt.minLabel_,opt.maxLabel_+1);
        }
        else{
            std::cout<<"addNodes explicit \n";    
            if(opt.maxLabel_>0){
                std::cout<<"reserve max node id \n"; 
                rag.reserveMaxNodeId(opt.maxLabel_);
            }

            if(opt.ignoreLabel_>=0){
                const LabelType ignoreLabel = static_cast<LabelType>(opt.ignoreLabel_);
                for(NodeIt iter(g);iter!=lemon::INVALID;++iter){
                    const LabelType l=labels[*iter];
                    if(static_cast<Int64>(l)!=ignoreLabel){
                        rag.addNode(l);
                    }
                }
            }
            else{
                for(NodeIt iter(g);iter!=lemon::INVALID;++iter){
                    const LabelType l=labels[*iter];
                    rag.addNode(l);
                }
            }
        }
    }


    template<class GRAPH,class LABEL_MAP>
    void addNodesParallel(
        const GRAPH & g,  
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"addNodes Parallel (fallback serial)\n";
        addNodesSerial(g, labels, rag, affiliatedEdges, opt);
    }





    template<class GRAPH,class LABEL_MAP>
    inline void addNodes(
        const GRAPH & g, 
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        if(opt.parallel_)
            addNodesParallel(g, labels, rag, affiliatedEdges, opt);
        else
            addNodesSerial(g, labels, rag, affiliatedEdges, opt);
    }
    




    template<class GRAPH,class LABEL_MAP>
    void addEdgesSerial(
        const GRAPH & g, 
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"addEdges Serial \n";
        typedef typename LABEL_MAP::Value LabelType;
        typedef typename GRAPH::Edge Edge;
        typedef typename GRAPH::EdgeIt EdgeIt;

        if(opt.ignoreLabel_>=0){
            const LabelType ignoreLabel = static_cast<LabelType>(opt.ignoreLabel_);
            for(EdgeIt e(g);e!=lemon::INVALID;++e){
                const Edge edge(*e);
                const LabelType lu = labels[g.u(edge)];
                const LabelType lv = labels[g.v(edge)];
                if(  lu!=lv &&  lu!=ignoreLabel  && lv!=ignoreLabel){
                    // if there is an edge between lu and lv no new edge will be added
                    rag.addEdge( rag.nodeFromId(lu),rag.nodeFromId(lv));
                }
            }
        }
        else{
            for(EdgeIt e(g);e!=lemon::INVALID;++e){
                const Edge edge(*e);
                const LabelType lu = labels[g.u(edge)];
                const LabelType lv = labels[g.v(edge)];
                if(  lu!=lv ){
                    // if there is an edge between lu and lv no new edge will be added
                    rag.addEdge( rag.nodeFromId(lu),rag.nodeFromId(lv));
                }
            }
        }
    }

    template<class GRAPH,class LABEL_MAP>
    void addEdgesParallel(
        const GRAPH & g,  
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"addEdges Parallel Serial Fallback\n";
        addEdgesSerial(g, labels, rag, affiliatedEdges, opt);
    }

    template<class GRAPH,class LABEL_MAP>
    void addEdges(
        const GRAPH & g, 
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"addEdges\n";
        if(opt.parallel_)
            addEdgesParallel(g, labels, rag, affiliatedEdges, opt);
        else
            addEdgesSerial(g, labels, rag, affiliatedEdges, opt);
    }









    template<class GRAPH,class LABEL_MAP>
    void setupAffiliatedEdgesSerial(
        const GRAPH & g, 
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"setupAffiliatedEdges Serial \n";

        typedef typename LABEL_MAP::Value LabelType;
        typedef typename GRAPH::Edge Edge;
        typedef typename GRAPH::EdgeIt EdgeIt;
        typedef typename AdjacencyListGraph::Edge RagEdge;
        // SET UP HYPEREDGES
        affiliatedEdges.assign(rag);

        if(opt.ignoreLabel_>=0){
            const LabelType ignoreLabel = static_cast<LabelType>(opt.ignoreLabel_);

            for(EdgeIt e(g);e!=lemon::INVALID;++e){
                const Edge edge(*e);
                const LabelType lu = labels[g.u(edge)];
                const LabelType lv = labels[g.v(edge)];
                //std::cout<<"edge between ?? "<<lu<<" "<<lv<<"\n";
                if(  lu!=lv && lu!=ignoreLabel  && lv!=ignoreLabel ){
                    RagEdge ragEdge= rag.findEdge(rag.nodeFromId(lu),rag.nodeFromId(lv));
                    affiliatedEdges[ragEdge].push_back(edge);
                }
            }
        }
    }

    template<class GRAPH,class LABEL_MAP>
    struct RagData{
        RagData(
            const GRAPH & g, 
            const LABEL_MAP & labels,
            AdjacencyListGraph & rag,
            typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
            const RagOptions & opt
        )
        :   g_(g),
            labels_(labels),
            rag_(rag),
            affiliatedEdges_(affiliatedEdges),
            opt_(opt){
        }

        const GRAPH & g_;
        const LABEL_MAP & labels_;
        AdjacencyListGraph & rag_;
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges_;
        const RagOptions & opt_;
    };


    template<class GRAPH, class LABEL_MAP>
    struct SetupAffiliatedEdgesWorker{
        typedef typename  GRAPH::Edge Edge;
        typedef typename  GRAPH::EdgeIt EdgeIt;
        typedef typename LABEL_MAP::Value LabelType;
        typedef typename AdjacencyListGraph::Edge RagEdge;
        typedef typename GRAPH::index_type index_type;
        SetupAffiliatedEdgesWorker(
            RagData<GRAPH, LABEL_MAP> & ragData,
            std::vector<boost::mutex > & edgeLock,
            const Edge & beginEdge,
            const Edge & endEdge
        )
        :   ragData_(&ragData),
            edgeLock_(&edgeLock),
            beginEdge_(beginEdge),
            endEdge_(endEdge)
        {

        }


        void operator()(){
            //std::cout<<"do stuff parallel_\n";
            //std::cout<<"get begin iter\n";
            vigra_precondition(beginEdge_!=lemon::INVALID,"invalid edge");
            EdgeIt beginIter(ragData_->g_, beginEdge_);



            //std::cout<<"get end iter\n";
            EdgeIt endIter  =  endEdge_==lemon::INVALID ? 
                GraphIteratorAccessor<GRAPH>::edgesEnd(ragData_->g_) :
                EdgeIt(ragData_->g_, endEdge_) ;

            //( endEdge_==lemon::INVALID )
           //      std::cout<<"end is REAL end  \n";    
           // std::cout<<"start loop  \n";    

            while(beginIter!=endIter && beginIter!=lemon::INVALID){
                //std::cout<<"----\n";
                const Edge edge(*beginIter);
                //std::cout<<"get uv labels from  "<<edge[0]<<" "<<edge[1]<<" "<<edge[2]<<" "<<edge[3]<<"\n";
                const LabelType lu = ragData_->labels_[ragData_->g_.u(edge)];
                const LabelType lv = ragData_->labels_[ragData_->g_.v(edge)];

                if(lu!=lv){
                    //std::cout<<"u=!v\n";
                    RagEdge ragEdge= ragData_->rag_.findEdge(
                        ragData_->rag_.nodeFromId(lu),
                        ragData_->rag_.nodeFromId(lv)
                    );
                    const index_type ragEdgeId = ragData_->rag_.id(ragEdge);
                    (*edgeLock_)[ragEdgeId].lock();
                    ragData_->affiliatedEdges_[ragEdge].push_back(edge);
                    (*edgeLock_)[ragEdgeId].unlock();
                }
                else{
                    //std::cout<<"u==v\n";
                }

                ++beginIter;
            }

        }


        RagData<GRAPH, LABEL_MAP> * ragData_;
        std::vector< boost::mutex> * edgeLock_;
        Edge beginEdge_;
        Edge endEdge_;
    };



    template<class GRAPH,class LABEL_MAP>
    void setupAffiliatedEdgesParallel(
        const GRAPH & g,  
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"setupAffiliatedEdges Parallel Serial Fallback \n";
        
        typedef typename MakeEdgeRangeVector<GRAPH>::EdgeRangeVector EdgeRangeVector;
        typedef typename GRAPH::Edge Edge;
        typedef SetupAffiliatedEdgesWorker<GRAPH, LABEL_MAP> WorkType;


        // SET UP HYPEREDGES
        affiliatedEdges.assign(rag);


        // get edge range 
        EdgeRangeVector edgeRangeVec;


        edgeRanges(g, edgeRangeVec, 40);



        // make rag data
        RagData<GRAPH, LABEL_MAP> ragData(g, labels, rag, affiliatedEdges, opt);

        //std::cout<<"rag edgeNum "<<rag.edgeNum()<<"\n";
        //std::cout<<"rag.maxEdgeId "<<rag.maxEdgeId()<<"\n";

        std::vector< boost::mutex>  edgeLock(rag.maxEdgeId()+1);

        {
            ThreadPool threadPool(opt.nThreads_);

            for(size_t i=0; i<edgeRangeVec.size(); ++i){
                //std::cout<<"i "<<i<<"\n";
                const Edge beginEdge = edgeRangeVec[i].first;
                const Edge endEdge = edgeRangeVec[i].second;

                //std::cout<<"start pwork \n";
                //std::cout<<"edge id range"<<g.id(beginEdge)<<" "<<g.id(endEdge)<<" max Id"<<g.maxEdgeId()<< "\n";

                //WorkType work(ragData, edgeLock, beginEdge, endEdge);
                //work();
                threadPool.run_task( WorkType(ragData, edgeLock, beginEdge, endEdge));
            }
        }

    }


    template<class GRAPH,class LABEL_MAP>
    inline void setupAffiliatedEdges(
        const GRAPH & g,  
        const LABEL_MAP & labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & opt
    ){
        std::cout<<"setupAffiliatedEdges \n";
        if(opt.parallel_)
            setupAffiliatedEdgesParallel(g, labels, rag, affiliatedEdges, opt);
        else
            setupAffiliatedEdgesSerial(g, labels, rag, affiliatedEdges, opt);
    }


    } // end namespace rag_detail



    template<class GRAPH,class LABEL_MAP>
    void makeRag(
        const GRAPH & g, 
        LABEL_MAP labels,
        AdjacencyListGraph & rag,
        typename AffiliatedEdgesHelper<GRAPH>::AffiliatedEdgesType & affiliatedEdges,
        const RagOptions & options = RagOptions()
    ){
        std::cout<<"make Rag \n";

        RagOptions opt=options;

        // check if input matches options
        // and get a new options object
        rag_detail::preprocess(g, labels, rag, affiliatedEdges, options, opt);

        // add nodes to the graph
        rag_detail::addNodes(g, labels, rag, affiliatedEdges, opt);

        // add edges to the graph
        rag_detail::addEdges(g, labels, rag, affiliatedEdges, opt);

        // setup affiliated edges
        rag_detail::setupAffiliatedEdges(g, labels, rag, affiliatedEdges, opt);
    }

}   

#endif /* VIGRA_GRAPH_RAG_HXX */

