#ifndef VIGRA_GRAPH_FEATURES_HXX
#define  VIGRA_GRAPH_FEATURES_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>





namespace vigra{




template<
    class LABELS_ARRAY,
    class DATA_ARRAY,
    class NODE_ACC_CHAIN,
    class EDGE_ACC_CHAIN
>
void nodeAndEdgeAccumlation(
    const AdjacencyListGraph &  graph,
    const LABELS_ARRAY &        labelsArray,
    const DATA_ARRAY &          dataArray,
    NODE_ACC_CHAIN &            nodeAccChain, 
    EDGE_ACC_CHAIN &            edgeAccChain
){
    typedef typename AdjacencyListGraph::Edge Edge;
    typedef typename AdjacencyListGraph::Node Node;
    typedef typename MultiArrayShape<LABELS_ARRAY::ActualDimension>::type Coord;
    typedef typename DATA_ARRAY::value_type DataType;
    typedef typename LABELS_ARRAY::value_type LabelType;

    const size_t nPassesNodeAcc = nodeAccChain.front().passesRequired();
    const size_t nPassesEdgeAcc = edgeAccChain.front().passesRequired();
    const size_t nPasses = std::max(nPassesNodeAcc, nPassesEdgeAcc);
    Coord shape = labelsArray.shape();


    for(size_t p=0; p<nPasses; ++p)
    for(size_t coord[2]=0; coord[2]<shape[2]; ++coord[2])
    for(size_t coord[1]=0; coord[1]<shape[1]; ++coord[1])
    for(size_t coord[0]=0; coord[0]<shape[0]; ++coord[0]){

        
        const LabelType uLabel = labelsArray[coord];
        const Node uNode = graph.nodeFromId(uLabel);
        const DataType uVal = dataArray[coord];

        if(p<nPassesNodeAcc){
            nodeAccChain[uNode].updatePassN(uNode, coord, p+1);
        }

        if(p<nPassesEdgeAcc){
            for(size_t d=0; d<LABELS_ARRAY::ActualDimension; ++d){
                Coord otherCoord = coord;
                otherCoord[d] += 1;
                if(otherCoord[d] < shape[d]){
                    const LabelType vLabel = labelsArray[otherCoord];

                    if(uLabel != vLabel){

                        const Node vNode = graph.nodeFromId(vLabel);
                        const DataType  vVal  = dataArray[coord];
                        const Edge e  = graph.findEdge(uNode, vNode);

                        edgeAccChain[e].updatePassN(uVal, coord, p+1); 
                        edgeAccChain[e].updatePassN(vVal, otherCoord, p+1); 

                    }
                }
            }
        }
    }
}












template<
    unsigned int DIM, 
    class LABEL_TYPE
>
class GridRagFeatureExtractor{



public:



    GridRagFeatureExtractor(
        const AdjacencyListGraph & graph,
        const MultiArrayView<DIM, LABEL_TYPE> & labels            
    )
    :   graph_(graph),
        labels_(labels){

    }

    template<
        class DATA_TYPE,
        class OUT_TYPE
    >
    void accumulateFeature(
        vigra::MultiArrayView<DIM, DATA_TYPE> & data,
        vigra::MultiArrayView<2, OUT_TYPE> & features 
    ){ 
        typedef acc::Select< 
            acc::DataArg<1>, acc::Variance, 
            acc::Mean, acc::StdDev, acc::Minimum, 
            acc::Maximum, acc::RootMeanSquares, 
            acc::Skewness, acc::Covariance, acc::RegionCenter
        >  SelectType;    

        typedef acc::StandAloneAccumulatorChain<3, DATA_TYPE, SelectType> FreeChain;

        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;

        NodeChainMap nodeAccChainMap(graph_);
        EdgeChainMap edgeAccChainMap(graph_);
        nodeAndEdgeAccumlation(graph_, labels_, data, nodeAccChainMap, edgeAccChainMap);
    }


private:
    const AdjacencyListGraph & graph_;
    const MultiArrayView<DIM, LABEL_TYPE> & labels_;        
};


}


#endif /* VIGRA_GRAPH_FEATURES_HXX */
