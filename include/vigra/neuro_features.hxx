/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
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



#ifndef NEURO_CLUSTER_OPERATOR
#define NEURO_CLUSTER_OPERATOR



/*std*/
#include <queue>          
#include <iomanip>

/*vigra*/
#include "priority_queue.hxx"
#include "priority_queue.hxx"
#include "metrics.hxx" 
#include "merge_graph_adaptor.hxx"
#include "adjacency_list_graph.hxx"

namespace vigra{      




    class NeuroDynamicFeatures{
    public:
        typedef NeuroDynamicFeatures SelfType;
        typedef AdjacencyListGraph Graph;
        typedef MergeGraphAdaptor<Graph> MergeGraph;



        typedef typename MergeGraph::MergeNodeCallBackType MergeNodeCallBackType;
        typedef typename MergeGraph::MergeEdgeCallBackType MergeEdgeCallBackType;
        typedef typename MergeGraph::EraseEdgeCallBackType EraseEdgeCallBackType;


        typedef MergeGraph::Edge MgEdge;
        typedef MergeGraph::Node MgNode;
        typedef Graph::Edge Edge;
        typedef Graph::Node Node;
        typedef MergeGraph::index_type index_type;


        NeuroDynamicFeatures(const Graph & graph, MergeGraph & mergeGraph)
        :   graph_(&graph),
            mergeGraph_(&mergeGraph){

        }



        void assignEdgeCues(const  MultiArrayView<2, float> & edgeCues ){
            edgeCues_= edgeCues;
        }
        void assignNodeCues(const  MultiArrayView<2, float> & nodeCues ){
            nodeCues_ = nodeCues;
        }
        void assignEdgeSizes(const MultiArrayView<1, float> & edgeSizes){
            edgeSizes_ = edgeSizes;
        }
        void assignNodeSizes(const MultiArrayView<1, float> & nodeSizes){
            nodeSizes_ = nodeSizes;
        }


        void mergeNodes(const MgNode & nodeA, const MgNode & nodeB){
            const index_type idA  = mergeGraph_->id(nodeA);
            const index_type idB  = mergeGraph_->id(nodeB);
            std::cout<<"merge nodes "<<idA<<" "<<idB<<"\n";

            // merge node sizes
            const float sizeA = nodeSizes_[idA];
            const float sizeB = nodeSizes_[idB];
            const float sizeAB = sizeA + sizeB;
            nodeSizes_[idA] = sizeAB;

            // merge node cues
            for(size_t f=0; f<nodeCues_.shape(1); ++f){
                nodeCues_(idA, f) = (sizeA*nodeCues_(idA, f) + sizeB*nodeCues_(idB, f))/sizeAB;
            }


        }   
        
        void mergeEdges(const MgEdge & edgeA, const MgEdge & edgeB){
            const index_type idA  = mergeGraph_->id(edgeA);
            const index_type idB  = mergeGraph_->id(edgeB);
            std::cout<<"merge edges "<<idA<<" "<<idB<<"\n";

        }

        void eraseEdge(const MgEdge & edge){
            const index_type id  = mergeGraph_->id(edge);
            std::cout<<"erase edge "<<id<<"\n";
        }



        void registerCallbacks(){
            MergeNodeCallBackType cbMn(MergeNodeCallBackType:: template from_method<SelfType,&SelfType::mergeNodes>(this));
            MergeEdgeCallBackType cbMe(MergeEdgeCallBackType:: template from_method<SelfType,&SelfType::mergeEdges>(this));
            EraseEdgeCallBackType cbEe(EraseEdgeCallBackType:: template from_method<SelfType,&SelfType::eraseEdge>(this));

            mergeGraph_->registerMergeNodeCallBack(cbMn);
            mergeGraph_->registerMergeEdgeCallBack(cbMe);
            mergeGraph_->registerEraseEdgeCallBack(cbEe);

        }

    private:

        const Graph * graph_;
        MergeGraph * mergeGraph_;


        MultiArray<2, float> edgeCues_;
        MultiArray<2, float> nodeCues_;
        MultiArray<1, float> edgeSizes_;
        MultiArray<1, float> nodeSizes_;

    };





}

#endif // NEURO_CLUSTER_OPERATOR
