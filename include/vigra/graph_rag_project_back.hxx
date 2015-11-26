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

#ifndef VIGRA_GRAPH_RAG_PROJECT_BACK_HXX
#define VIGRA_GRAPH_RAG_PROJECT_BACK_HXX

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



namespace vigra{

    /// \cond
    namespace detail_rag_project_back{

    template<
        class BASE_GRAPH,
        class BASE_GRAPH_LABELS,
        class RAG_FEATURES,
        class BASE_GRAPH_FEATURES
    >
    struct RagProjectBack{


        static void projectBack(
            const AdjacencyListGraph & rag,
            const BASE_GRAPH & bg,
            const Int64 ignoreLabel,
            const BASE_GRAPH_LABELS bgLabels,
            const RAG_FEATURES & ragFeatures,
            BASE_GRAPH_FEATURES & bgFeatures
        ){
            typedef BASE_GRAPH Bg;
            typedef typename Bg::NodeIt BgNodeIt;
            typedef typename Bg::Node BgNode;

            if(ignoreLabel==-1){
                for(BgNodeIt iter(bg); iter!=lemon::INVALID; ++iter){
                    const BgNode bgNode(*iter);
                    bgFeatures[bgNode] = ragFeatures[rag.nodeFromId(bgLabels[bgNode])];
                }
            }
            else{
                for(BgNodeIt iter(bg); iter!=lemon::INVALID; ++iter){
                    const BgNode bgNode(*iter);
                    if(static_cast<Int64>(bgLabels[bgNode])!=ignoreLabel)
                        bgFeatures[bgNode] = ragFeatures[rag.nodeFromId(bgLabels[bgNode])];
                }
            }
        }
    };


    template<
        class BASE_GRAPH_LABELS,
        class RAG_FEATURES,
        class BASE_GRAPH_FEATURES
    >
    struct RagProjectBack<
        vigra::GridGraph<3, boost::undirected_tag>,
        BASE_GRAPH_LABELS,
        RAG_FEATURES,
        BASE_GRAPH_FEATURES
    >{
        typedef vigra::GridGraph<3, boost::undirected_tag> BASE_GRAPH;

        static void projectBack(
            const AdjacencyListGraph & rag,
            const BASE_GRAPH & bg,
            const Int64 ignoreLabel,
            const BASE_GRAPH_LABELS bgLabels,
            const RAG_FEATURES & ragFeatures,
            BASE_GRAPH_FEATURES & bgFeatures
        ){
            typedef BASE_GRAPH Bg;
            typedef typename Bg::Node BgNode;


            vigra::TinyVector<Int64, 3> shape = bg.shape();


            if(ignoreLabel==-1){
                
// FIXME: replace with threadpool                #pragma omp parallel for
                for(Int64 z=0; z<shape[2]; ++z){    
                    BgNode node;
                    node[2]=z;
                    for(node[1]=0; node[1]<shape[1]; ++node[1])  
                    for(node[0]=0; node[0]<shape[0]; ++node[0]){
                        bgFeatures[node] = ragFeatures[rag.nodeFromId(bgLabels[node])];
                    }  
                }

            }
            else{
// FIXME: replace with threadpool                #pragma omp parallel for
                for(Int64 z=0; z<shape[2]; ++z){    
                    BgNode node;
                    node[2]=z;
                    for(node[1]=0; node[1]<shape[1]; ++node[1])  
                    for(node[0]=0; node[0]<shape[0]; ++node[0]){
                        if(static_cast<Int64>(bgLabels[node])!=ignoreLabel)
                            bgFeatures[node] = ragFeatures[rag.nodeFromId(bgLabels[node])];
                    }  
                }
            }
        }
    };



    }
    /// \endcond

    /// project node features of a region adjacency
    /// graph back to the base graph.
    ///
    /// This function can be used to show a segmentation
    /// or node features of RAG on pixel / voxel level
    template< class BASE_GRAPH,
                class BASE_GRAPH_LABELS,
                class RAG_FEATURES,
                class BASE_GRAPH_FEATURES 
    >
    inline void projectBack(
            const AdjacencyListGraph & rag,
            const BASE_GRAPH & bg,
            const Int64 ignoreLabel,
            const BASE_GRAPH_LABELS bgLabels,
            const RAG_FEATURES & ragFeatures,
            BASE_GRAPH_FEATURES & bgFeatures
    ){
        using namespace detail_rag_project_back;
        detail_rag_project_back::RagProjectBack< BASE_GRAPH,BASE_GRAPH_LABELS,RAG_FEATURES,BASE_GRAPH_FEATURES>::projectBack(rag,
            bg,ignoreLabel,bgLabels,ragFeatures,bgFeatures);
    }



}

#endif /* VIGRA_GRAPH_RAG_PROJECT_BACK_HXX */

