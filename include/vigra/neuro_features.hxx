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

        void assignLabels(const MultiArrayView<1, UInt32> & edgeLabels){
            edgeLabelProb_.reshape(
                MultiArrayView<2, UInt32>::difference_type(edgeLabels.shape(0),4)
            );
            for(size_t i=0; i<edgeLabels.shape(0); ++i){
                // edge is inactive
                if(edgeLabels(i)==0){
                    edgeLabelProb_(i,0) = 1.0;
                    edgeLabelProb_(i,1) = 0.0;
                    edgeLabelProb_(i,2) = 0.0;
                }
                else if(edgeLabels(i)==1){
                    edgeLabelProb_(i,0) = 0.0;
                    edgeLabelProb_(i,1) = 1.0;
                    edgeLabelProb_(i,2) = 0.0;
                }
                else{
                    edgeLabelProb_(i,0) = 0.0;
                    edgeLabelProb_(i,1) = 0.0;
                    edgeLabelProb_(i,2) = 1.0;
                }
            }
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

            // merge edge sizes
            const float sizeA = edgeSizes_[idA];
            const float sizeB = edgeSizes_[idB];
            const float sizeAB = sizeA + sizeB;

            edgeSizes_[idA] = sizeAB;

            // merge edge cues
            for(size_t f=0; f<edgeCues_.shape(1); ++f){
                edgeCues_(idA, f) = (sizeA*edgeCues_(idA, f) + sizeB*edgeCues_(idB, f))/sizeAB;
            }

            if(edgeLabelProb_.size()>0){
                edgeLabelProb_(idA,0) = (sizeA*edgeLabelProb_(idA,0) + sizeB*edgeLabelProb_(idB,0))/sizeAB;
                edgeLabelProb_(idA,1) = (sizeA*edgeLabelProb_(idA,1) + sizeB*edgeLabelProb_(idB,1))/sizeAB;
                edgeLabelProb_(idA,2) = (sizeA*edgeLabelProb_(idA,2) + sizeB*edgeLabelProb_(idB,2))/sizeAB;
            }
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


        size_t numberOfFeatures()const{
            const size_t nMeanEdgeCues = edgeCues_.shape(1);
            const size_t nDegreeFeatures = 5;
            const size_t nSizeFeatures = 24;
            return  nMeanEdgeCues + nDegreeFeatures + nSizeFeatures;
        }


        void computeFeature(
            const MgEdge & edge, 
            MultiArrayView<1, float> & features
        ){
            size_t featureIndex = 0;

            meanEdgeCues(edge, features, featureIndex);
            degreeFeatures(edge, features, featureIndex);
            sizeFeatures(edge, features, featureIndex);
        }


        void meanEdgeCues(
            const MgEdge & edge, 
            MultiArrayView<1, float> & feature,
            size_t & featureIndex
        ){
            const int edgeId = mergeGraph_->id(edge);
            for (int i = 0; i < edgeCues_.shape(1); ++i){
                feature[featureIndex++] = edgeCues_(edgeId, i);
            }
        }

        void degreeFeatures(
            const MgEdge & edge, 
            MultiArrayView<1, float> & feature,
            size_t & featureIndex
        ){
            const MgNode u = mergeGraph_->u(edge);
            const MgNode v = mergeGraph_->v(edge);
            const size_t dU = mergeGraph_->degree(u);
            const size_t dV = mergeGraph_->degree(v);

            feature[featureIndex++] = std::min(dU, dV);
            feature[featureIndex++] = std::max(dU, dV);
            feature[featureIndex++] = float(std::max(dU, dV))/std::min(dU, dV);
            feature[featureIndex++] = std::max(dU, dV)-std::min(dU, dV);
            feature[featureIndex++] = dU + dV;
        }

        void sizeFeatures(
            const MgEdge & edge, 
            MultiArrayView<1, float> & feature,
            size_t & featureIndex
        ){
            const MgNode u = mergeGraph_->u(edge);
            const MgNode v = mergeGraph_->v(edge);

            const float edgeSize = edgeSizes_[mergeGraph_->id(edge)];

            const float sizeU = nodeSizes_[mergeGraph_->id(u)];
            const float sizeV = nodeSizes_[mergeGraph_->id(v)];
            const float sizeUV = sizeU+sizeV;
            const float maxSize = std::max(sizeU,sizeV);
            const float minSize = std::min(sizeU,sizeV);

            const float surfSizeU = computeSurfaceSize(u);
            const float surfSizeV = computeSurfaceSize(v);
            const float surfSizeUV = surfSizeU + surfSizeV;
            const float maxSurfSize = std::max(surfSizeU,surfSizeV);
            const float minSurfSize = std::min(surfSizeU,surfSizeV);

            const float sizeSurfRatioU = std::pow(sizeU, 1.0/3.0) / std::sqrt(surfSizeU);
            const float sizeSurfRatioV = std::pow(sizeV, 1.0/3.0) / std::sqrt(surfSizeV);
            const float sizeSurfRatioUV = std::pow(sizeUV, 1.0/3.0) / std::sqrt(surfSizeUV-2.0*edgeSize);

            // node size
            feature[featureIndex++] = sizeUV;
            feature[featureIndex++] = maxSize;
            feature[featureIndex++] = minSize;
            
            // surface size
            feature[featureIndex++] = surfSizeUV;
            feature[featureIndex++] = maxSurfSize;
            feature[featureIndex++] = minSurfSize;

            // size surface ratios
            feature[featureIndex++] = sizeSurfRatioU+sizeSurfRatioV;
            feature[featureIndex++] = std::min(sizeSurfRatioU, sizeSurfRatioV);
            feature[featureIndex++] = std::max(sizeSurfRatioU, sizeSurfRatioV);
            feature[featureIndex++] = sizeSurfRatioUV;
            feature[featureIndex++] = sizeSurfRatioUV-std::min(sizeSurfRatioU, sizeSurfRatioV);
            feature[featureIndex++] = sizeSurfRatioUV-std::max(sizeSurfRatioU, sizeSurfRatioV);
            feature[featureIndex++] = sizeSurfRatioUV/std::min(sizeSurfRatioU, sizeSurfRatioV);
            feature[featureIndex++] = sizeSurfRatioUV/std::max(sizeSurfRatioU, sizeSurfRatioV);

            // edge size
            feature[featureIndex++] = edgeSize;

            // edge Size / node size ratio
            feature[featureIndex++] = edgeSize / sizeUV;
            feature[featureIndex++] = edgeSize / maxSize;
            feature[featureIndex++] = edgeSize / minSize;

            feature[featureIndex++] = edgeSize / std::pow(sizeUV , 2.0/3.0);
            feature[featureIndex++] = edgeSize / std::pow(maxSize , 2.0/3.0);
            feature[featureIndex++] = edgeSize / std::pow(minSize , 2.0/3.0);


            // edge size surfaca fraction
            feature[featureIndex++] = edgeSize / surfSizeUV;
            feature[featureIndex++] = edgeSize / maxSurfSize;
            feature[featureIndex++] = edgeSize / minSurfSize;

        }




    private:

        float computeSurfaceSize(const MgNode node)const{
            float surfSize = 0.0;
            for(MergeGraph::IncEdgeIt iter(*mergeGraph_,node); iter!=lemon::INVALID; ++iter){
                surfSize+=edgeSizes_[mergeGraph_->id(*iter)];
            }
            return surfSize;
        }


        const Graph * graph_;
        MergeGraph * mergeGraph_;

        MultiArray<1, float> edgeSizes_;
        MultiArray<1, float> nodeSizes_;

        MultiArray<2, float> edgeCues_;
        MultiArray<2, float> nodeCues_;


        // raw histogram and pmap histogram

        std::vector < MultiArray<2, float>  > edgeHistograms_;
        std::vector < MultiArray<2, float>  > nodeHistograms_;







        MultiArray<2, float> edgeLabelProb_;


    };



    void mergeMeanAndVariance(
        float & meanA,
        float & varA,
        const float sizeA,
        const float meanB,
        const float varB,
        const float sizeB
    ){
        const float sizeAB = sizeA + sizeB;
        const float newMean = (sizeA*meanA + sizeB*meanB)/sizeAB;

        const float  fac =  1.0/(sizeAB - 1);

        varA = fac * (sizeA * std::pow(meanA - newMean, 2) + sizeB * std::pow(meanB - newMean, 2)   +  (sizeA-1)*varA + (sizeB-1)*varB );
        meanA = newMean;
    }  
 
}

#endif // NEURO_CLUSTER_OPERATOR
