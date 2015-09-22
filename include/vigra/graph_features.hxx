#ifndef VIGRA_GRAPH_FEATURES_HXX
#define  VIGRA_GRAPH_FEATURES_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>





namespace vigra{




template<
    unsigned int DIM,
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
    typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
    typedef typename AdjacencyListGraph::NodeIt NodeIt;
    typedef typename MultiArrayShape<DIM>::type Coord;
    typedef typename DATA_ARRAY::value_type DataType;
    typedef typename LABELS_ARRAY::value_type LabelType;

    


    const size_t nPassesNodeAcc = nodeAccChain[*NodeIt(graph)].passesRequired();
    const size_t nPassesEdgeAcc = edgeAccChain[*EdgeIt(graph)].passesRequired();


    const size_t nPasses = std::max(nPassesNodeAcc, nPassesEdgeAcc);



    Coord shape(labelsArray.shape());
    Coord coord(0);

    if(DIM == 2){

        for(size_t p=0; p<nPasses; ++p)
        for(coord[1]=0; coord[1]<shape[1]; ++coord[1])
        for(coord[0]=0; coord[0]<shape[0]; ++coord[0]){

            
            const LabelType uLabel = labelsArray[coord];
            const Node uNode = graph.nodeFromId(uLabel);
            const DataType uVal = dataArray[coord];

            if(p<nPassesNodeAcc){
                nodeAccChain[uNode].updatePassN(uVal, coord, p+1);
            }

           
            if(p<nPassesEdgeAcc){


                for(size_t d=0; d<DIM; ++d){
                    Coord otherCoord = coord;
                    otherCoord[d] += 1;


                    if(otherCoord[d] < shape[d]){
                        const LabelType vLabel = labelsArray[otherCoord];
                        
                        if(uLabel != vLabel){
                            const Node vNode = graph.nodeFromId(vLabel);
                            const DataType  vVal  = dataArray[otherCoord];
                            const Edge e  = graph.findEdge(uNode, vNode);

    
                            edgeAccChain[e].updatePassN(uVal, coord, p+1); 
                            edgeAccChain[e].updatePassN(vVal, otherCoord, p+1); 

                        }
                    }
                }
            }
        }
    }
    if(DIM == 3){
        for(size_t p=0; p<nPasses; ++p)
        for(coord[2]=0; coord[2]<shape[2]; ++coord[2])
        for(coord[1]=0; coord[1]<shape[1]; ++coord[1])
        for(coord[0]=0; coord[0]<shape[0]; ++coord[0]){

            
            const LabelType uLabel = labelsArray[coord];
            const Node uNode = graph.nodeFromId(uLabel);
            const DataType uVal = dataArray[coord];

            if(p<nPassesNodeAcc){
                nodeAccChain[uNode].updatePassN(uVal, coord, p+1);
            }
            if(p<nPassesEdgeAcc){
                for(size_t d=0; d<DIM; ++d){
                    Coord otherCoord = coord;
                    otherCoord[d] += 1;
                    if(otherCoord[d] < shape[d]){
                        const LabelType vLabel = labelsArray[otherCoord];
                        if(uLabel != vLabel){

                            const Node vNode = graph.nodeFromId(vLabel);
                            const DataType  vVal  = dataArray[otherCoord];
                            const Edge e  = graph.findEdge(uNode, vNode);
                            edgeAccChain[e].updatePassN(uVal, coord, p+1); 
                            edgeAccChain[e].updatePassN(vVal, otherCoord, p+1); 

                        }
                    }
                }
            }
        }
    }
}


template<
    unsigned int DIM,
    class LABELS_ARRAY,
    class NODE_ACC_CHAIN,
    class EDGE_ACC_CHAIN
>
void nodeAndEdgeAccumlation(
    const AdjacencyListGraph &  graph,
    const LABELS_ARRAY &        labelsArray,
    NODE_ACC_CHAIN &            nodeAccChain, 
    EDGE_ACC_CHAIN &            edgeAccChain
){
    typedef typename AdjacencyListGraph::Edge Edge;
    typedef typename AdjacencyListGraph::Node Node;
    typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
    typedef typename AdjacencyListGraph::NodeIt NodeIt;
    typedef typename MultiArrayShape<DIM>::type Coord;
    typedef typename LABELS_ARRAY::value_type LabelType;

    


    const size_t nPassesNodeAcc = nodeAccChain[*NodeIt(graph)].passesRequired();
    const size_t nPassesEdgeAcc = edgeAccChain[*EdgeIt(graph)].passesRequired();


    const size_t nPasses = std::max(nPassesNodeAcc, nPassesEdgeAcc);



    Coord shape(labelsArray.shape());
    Coord coord(0);

    if(DIM == 2){

        for(size_t p=0; p<nPasses; ++p)
        for(coord[1]=0; coord[1]<shape[1]; ++coord[1])
        for(coord[0]=0; coord[0]<shape[0]; ++coord[0]){

            
            const LabelType uLabel = labelsArray[coord];
            const Node uNode = graph.nodeFromId(uLabel);


            if(p<nPassesNodeAcc){
                nodeAccChain[uNode].updatePassN(coord, p+1);
            }

           
            if(p<nPassesEdgeAcc){


                for(size_t d=0; d<DIM; ++d){
                    Coord otherCoord = coord;
                    otherCoord[d] += 1;


                    if(otherCoord[d] < shape[d]){
                        const LabelType vLabel = labelsArray[otherCoord];
                        
                        if(uLabel != vLabel){
                            const Node vNode = graph.nodeFromId(vLabel);
                            const Edge e  = graph.findEdge(uNode, vNode);
                            edgeAccChain[e].updatePassN(coord, p+1); 
                            edgeAccChain[e].updatePassN(otherCoord, p+1); 

                        }
                    }
                }
            }
        }
    }
    if(DIM == 3){
        for(size_t p=0; p<nPasses; ++p)
        for(coord[2]=0; coord[2]<shape[2]; ++coord[2])
        for(coord[1]=0; coord[1]<shape[1]; ++coord[1])
        for(coord[0]=0; coord[0]<shape[0]; ++coord[0]){

            
            const LabelType uLabel = labelsArray[coord];
            const Node uNode = graph.nodeFromId(uLabel);


            if(p<nPassesNodeAcc){
                nodeAccChain[uNode].updatePassN( coord, p+1);
            }
            if(p<nPassesEdgeAcc){
                for(size_t d=0; d<DIM; ++d){
                    Coord otherCoord = coord;
                    otherCoord[d] += 1;
                    if(otherCoord[d] < shape[d]){
                        const LabelType vLabel = labelsArray[otherCoord];
                        if(uLabel != vLabel){
                            const Node vNode = graph.nodeFromId(vLabel);
                            const Edge e  = graph.findEdge(uNode, vNode);
                            edgeAccChain[e].updatePassN(coord, p+1); 
                            edgeAccChain[e].updatePassN(otherCoord, p+1); 

                        }
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

        typedef acc::UserRangeHistogram<40> Hist;
        typedef acc::StandardQuantiles<Hist> Quants;
        typedef typename AdjacencyListGraph::Edge Edge;
        typedef typename AdjacencyListGraph::Node Node;
        typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
        typedef typename AdjacencyListGraph::NodeIt NodeIt;

public:



    GridRagFeatureExtractor(
        const AdjacencyListGraph & graph,
        const MultiArrayView<DIM, LABEL_TYPE> & labels            
    )
    :   graph_(graph),
        labels_(labels){

    }

    UInt64 edgeNum()const{
        return graph_.edgeNum();
    }

    template<
        class DATA_TYPE,
        class OUT_TYPE
    >
    void accumulatedFeatures(
        const vigra::MultiArrayView<DIM, DATA_TYPE> & data,
        const DATA_TYPE minVal,
        const DATA_TYPE maxVal,
        vigra::MultiArrayView<2, OUT_TYPE> & features 
    )const{ 


        typedef acc::Select< 
            acc::DataArg<1>, acc::Variance, acc::Count,
            acc::Mean, acc::StdDev, acc::Minimum, 
            acc::Maximum, acc::RootMeanSquares, 
            acc::Skewness, acc::Covariance, acc::RegionCenter,
            Hist,Quants
        >  SelectType;    

        typedef acc::StandAloneAccumulatorChain<DIM, DATA_TYPE, SelectType> FreeChain;
        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;


        NodeChainMap nodeAccChainMap(graph_);
        EdgeChainMap edgeAccChainMap(graph_);



        vigra::HistogramOptions histogram_opt;
        //histogram_opt = histogram_opt.setBinCount(50);
        histogram_opt = histogram_opt.setMinMax(minVal, maxVal); 

        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){
            edgeAccChainMap[*eIt].setHistogramOptions(histogram_opt); 
        }
        for(NodeIt nIt(graph_); nIt != lemon::INVALID; ++nIt){
            nodeAccChainMap[*nIt].setHistogramOptions(histogram_opt); 
        }

        nodeAndEdgeAccumlation<DIM>(graph_, labels_, data, nodeAccChainMap, edgeAccChainMap);
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){



            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);
            const FreeChain & eChain = edgeAccChainMap[edge];
            const float mean = acc::get<acc::Mean>(eChain);
            features(eid, 0) = acc::get<acc::Count>(eChain);
            features(eid, 1) = acc::get<acc::Mean>(eChain);
            features(eid, 2) = acc::get<Quants>(eChain)[3];
        }
    }


    template<
        class OUT_TYPE
    >
    void geometricFeatures(
        vigra::MultiArrayView<2, OUT_TYPE> & features 
    )const{ 


        typedef acc::Select< 
            acc::Count,acc::RegionCenter
        >  SelectType;    

        typedef acc::StandAloneDataFreeAccumulatorChain<DIM, SelectType> FreeChain;
        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;




        NodeChainMap nodeAccChainMap(graph_);
        EdgeChainMap edgeAccChainMap(graph_);




        nodeAndEdgeAccumlation<DIM>(graph_, labels_, nodeAccChainMap, edgeAccChainMap);
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){

            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);
            const Node u = graph_.u(edge);
            const Node v = graph_.v(edge);

            const FreeChain & eChain = edgeAccChainMap[edge];
            const FreeChain & uChain = nodeAccChainMap[u];
            const FreeChain & vChain = nodeAccChainMap[v];

            const TinyVector<Int32, DIM> eCenter = acc::get<acc::RegionCenter>(eChain);
            const TinyVector<Int32, DIM> uCenter = acc::get<acc::RegionCenter>(uChain);
            const TinyVector<Int32, DIM> vCenter = acc::get<acc::RegionCenter>(vChain);


            features(eid,0) = squaredNorm(uCenter-vCenter);
            features(eid,1) = std::min( squaredNorm(uCenter-eCenter), squaredNorm(vCenter-eCenter));
            features(eid,2) = std::max( squaredNorm(uCenter-eCenter), squaredNorm(vCenter-eCenter));
            features(eid,3) = squaredNorm(uCenter-eCenter) + squaredNorm(vCenter-eCenter);
            features(eid,4) = std::abs(squaredNorm(uCenter-eCenter) - squaredNorm(vCenter-eCenter));
        }
    }




private:
    const AdjacencyListGraph & graph_;
    MultiArrayView<DIM, LABEL_TYPE> labels_;        
};


}


#endif /* VIGRA_GRAPH_FEATURES_HXX */
