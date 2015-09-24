#ifndef VIGRA_GRAPH_FEATURES_HXX
#define  VIGRA_GRAPH_FEATURES_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>





namespace vigra{

template<class T, class KEY>
struct ArrayMap{

    typedef MultiArrayView<1, T>  View;
    typedef KEY    Key;
    typedef typename View::value_type  Value;
    typedef typename View::const_reference  ConstReference;
    typedef typename View::reference  Reference;

    ArrayMap(const AdjacencyListGraph & graph)
    :   graph_(graph),
        view_(){
    }

    ArrayMap(const AdjacencyListGraph & graph, const View & view)
    :   graph_(graph),
        view_(view){
    }

    void setArray(const MultiArrayView<1, T> & view){
        view_ = view;
    }

    Reference operator[](const Key & key){
       return view_(graph_.id(key));
    }

    ConstReference operator[](const Key & key)const{
        return view_(graph_.id(key));
    }
    MultiArrayView<1, T> view_;
    AdjacencyListGraph & graph_;
};


double squaredHellinger(
    const double mA, const double mB,
    const double sA, const double sB
){
    const double ssA = sA*sA;
    const double ssB = sB*sB;

    const double eArg = -0.25*std::pow(mA-mB,2)/(ssA+ssB);
    const double sArg = (2*sA*sB)/(ssA+ssB);
    return 1.0 - std::sqrt(sArg)*std::exp(eArg);
}



// TinyVector<double, 5> 
// fHellinger(
    // const double mE, const double mU, const double mV,
    // const double sE, const double sU, const double sV
// ){
// 
    // const double hUV = squaredHellinger(mU, mV, sU, sV);
    // const double hUE = squaredHellinger(mU, mE, sU, sE);
    // const double hVE = squaredHellinger(mV, mE, sV, sE);
// 
// }






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

        typedef acc::AutoRangeHistogram<40> AHist;
        typedef acc::UserRangeHistogram<40> UHist;
        typedef acc::StandardQuantiles<UHist> Quants;
        typedef acc::StandardQuantiles<AHist> AQuants;
        typedef typename AdjacencyListGraph::Edge Edge;
        typedef typename AdjacencyListGraph::Node Node;
        typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
        typedef typename AdjacencyListGraph::NodeIt NodeIt;
        typedef typename AdjacencyListGraph::OutArcIt OutArcIt;
        typedef std::set<Node> NodeSet;
        typedef typename NodeSet::const_iterator NodeSetIter;
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







    template<class DATA_TYPE>
    UInt64 nAccumulatedFeaturesSimple()const{
        return 2*11;
    }

    template<
        class DATA_TYPE,
        class OUT_TYPE
    >
    void accumulatedFeaturesSimple(
        const vigra::MultiArrayView<DIM, DATA_TYPE> & data,
        vigra::MultiArrayView<2, OUT_TYPE> & features 
    )const{ 


        typedef acc::Select< 
            acc::DataArg<1>,
            acc::Mean, acc::StdDev
        >  SelectType;    

        typedef acc::StandAloneAccumulatorChain<DIM, DATA_TYPE, SelectType> FreeChain;
        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;


        NodeChainMap nodeAccChainMap(graph_);
        EdgeChainMap edgeAccChainMap(graph_);



    
        nodeAndEdgeAccumlation<DIM>(graph_, labels_, data, nodeAccChainMap, edgeAccChainMap);
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){




            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);
            vigra::MultiArrayView<1, OUT_TYPE> edgeFeat = features.bindInner(eid);
            const Node u = graph_.u(edge);
            const Node v = graph_.v(edge);

            const FreeChain & eChain = edgeAccChainMap[edge];
            const FreeChain & uChain = nodeAccChainMap[u];
            const FreeChain & vChain = nodeAccChainMap[v];
            const float mean = acc::get<acc::Mean>(eChain);

            const float eM = acc::get<acc::Mean>(eChain);
            const float uM = acc::get<acc::Mean>(uChain);
            const float vM = acc::get<acc::Mean>(vChain);
            const float eS = acc::get<acc::StdDev>(eChain);
            const float uS = acc::get<acc::StdDev>(uChain);
            const float vS = acc::get<acc::StdDev>(vChain);
       


            UInt64 fIndex = 0;

            defaultFeat(fIndex,edgeFeat, eM,uM,vM);
            defaultFeat(fIndex,edgeFeat, eS,uS,vS);
            //for(size_t qi=0; qi<7;++qi)
            //    defaultFeat(fIndex,edgeFeat, eQnt[qi],uQnt[qi],vQnt[qi]);
        }
    }






    template<class DATA_TYPE>
    UInt64 nAccumulatedFeatures()const{
        return 9*11;
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
            acc::DataArg<1>,
            acc::Mean, acc::StdDev, Quants
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
            vigra::MultiArrayView<1, OUT_TYPE> edgeFeat = features.bindInner(eid);
            const Node u = graph_.u(edge);
            const Node v = graph_.v(edge);

            const FreeChain & eChain = edgeAccChainMap[edge];
            const FreeChain & uChain = nodeAccChainMap[u];
            const FreeChain & vChain = nodeAccChainMap[v];
            const float mean = acc::get<acc::Mean>(eChain);

            const float eM = acc::get<acc::Mean>(eChain);
            const float uM = acc::get<acc::Mean>(uChain);
            const float vM = acc::get<acc::Mean>(vChain);
            const float eS = acc::get<acc::StdDev>(eChain);
            const float uS = acc::get<acc::StdDev>(uChain);
            const float vS = acc::get<acc::StdDev>(vChain);
            const TinyVector<float, 7> eQnt = acc::get<Quants>(eChain);
            const TinyVector<float, 7> uQnt = acc::get<Quants>(uChain);
            const TinyVector<float, 7> vQnt = acc::get<Quants>(vChain);


            UInt64 fIndex = 0;

            defaultFeat(fIndex,edgeFeat, eM,uM,vM);
            defaultFeat(fIndex,edgeFeat, eS,uS,vS);
            for(size_t qi=0; qi<7;++qi)
                defaultFeat(fIndex,edgeFeat, eQnt[qi],uQnt[qi],vQnt[qi]);
        }
    }


    // 9 features
    template<class OUT_TYPE>
    void defaultFeat(
        UInt64 &fIndex,
        vigra::MultiArrayView<1, OUT_TYPE> & features,
        const float e,
        const float u, 
        const float v
    )const{
        const float dUV = std::abs(u-v);
        const float eUE = std::abs(u-e);
        const float eVE = std::abs(v-e);
        const float minUV = std::min(u,v);
        const float maxUV = std::max(u,v);
        const float meanUV = (u+v)/2.0;

        features(fIndex++) = e;
        features(fIndex++) = u+v;
        features(fIndex++) = dUV;
        features(fIndex++) = minUV;
        features(fIndex++) = maxUV;
        features(fIndex++) = meanUV;
        features(fIndex++) = std::abs(e-minUV);
        features(fIndex++) = std::abs(e-maxUV);
        features(fIndex++) = std::abs(e-meanUV);
        features(fIndex++) = std::max(eUE,eVE) - dUV;
        features(fIndex++) = std::min(eUE,eVE) - dUV;
    }



    UInt64 nGeometricFeatures()const{
        return 21;
    }

    template<class OUT_TYPE>
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



        UInt64 fi = 0;
        nodeAndEdgeAccumlation<DIM>(graph_, labels_, nodeAccChainMap, edgeAccChainMap);
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){

            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);
            const Node u = graph_.u(edge);
            const Node v = graph_.v(edge);

            const FreeChain & eChain = edgeAccChainMap[edge];
            const FreeChain & uChain = nodeAccChainMap[u];
            const FreeChain & vChain = nodeAccChainMap[v];

            const TinyVector<float, DIM>  eCenter = acc::get<acc::RegionCenter>(eChain);
            const TinyVector<float, DIM>  uCenter = acc::get<acc::RegionCenter>(uChain);
            const TinyVector<float, DIM>  vCenter = acc::get<acc::RegionCenter>(vChain);

            const float eCount = acc::get<acc::Count>(eChain);
            const float uCount = acc::get<acc::Count>(uChain);
            const float vCount = acc::get<acc::Count>(vChain);


            const float eNCount = std::sqrt(eCount);
            const float uNCount = std::pow(uCount,1.0/3.0);
            const float vNCount = std::pow(vCount,1.0/3.0);
            const float uvNCount = std::pow(uCount + vCount,1.0/3.0);

            const float dUV = squaredNorm(uCenter-vCenter);
            const float dEU = squaredNorm(eCenter-uCenter);
            const float dEV = squaredNorm(eCenter-vCenter);



            const float uRat = uNCount/dEU;
            const float vRat = vNCount/dEV;
            const float uvRat = uvNCount/dEV;
            
            
            // count based 
            features(eid, 0) = eCount;                            
            features(eid, 1) = std::min(uCount, vCount);
            features(eid, 2) = std::max(uCount, vCount);
            features(eid, 3) = uCount + vCount;
            features(eid, 4) = std::min(uCount, vCount) / std::max(uCount, vCount);
            features(eid, 5) = std::abs(uCount - vCount);
            features(eid, 6) = eNCount / std::min(uNCount, vNCount);
            features(eid, 7) = eNCount / std::max(uNCount, vNCount);
            features(eid, 8) = eNCount / uvNCount;

            // center of mass based
            features(eid, 9) = dUV;
            features(eid,10) = std::min(dEU, dEV);
            features(eid,11) = std::max(dEU, dEV);
            features(eid,12) = dEU + dEV;
            features(eid,13) = std::min(dEU, dEV)/std::max(dEU, dEV);
            features(eid,14) = std::abs(dEU - dEV);


            // count and center of mass based
            features(eid,15) = std::log(uvRat);
            features(eid,16) = std::log(std::min(uRat, vRat));
            features(eid,17) = std::log(std::max(uRat, vRat));
            features(eid,18) = std::log(uRat + vRat);
            features(eid,19) = std::log(std::abs(uRat - vRat)+0.5);
            features(eid,20) = std::log(std::max(uRat, vRat) / std::min(uRat, vRat));
            
        }
    }


    UInt64 nTopologicalFeatures()const{
        return 5;
    }

    template<class OUT_TYPE>
    void topologicalFeatures(vigra::MultiArrayView<2, OUT_TYPE> & features )const{

        UInt64 fi = 0;
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){

            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);
            const Node u = graph_.u(edge);
            const Node v = graph_.v(edge);

            const Int32 uDeg = graph_.degree(u);
            const Int32 vDeg = graph_.degree(v);
            
            NodeSet uSet,vSet;

            for(OutArcIt outArcIt(graph_,u); outArcIt!=lemon::INVALID;++outArcIt)
               uSet.insert(graph_.target(*outArcIt));

            for(OutArcIt outArcIt(graph_,v); outArcIt!=lemon::INVALID;++outArcIt)
               vSet.insert(graph_.target(*outArcIt));

            UInt64 intersectionSize = 0;
            for(NodeSetIter siter=uSet.begin(); siter!=uSet.end(); ++siter){
                const Node & node = *siter;
                if(vSet.find(node)!=vSet.end())
                    ++intersectionSize;
            }



            // degree based                            
            features(eid, 0) = std::min(uDeg, vDeg);
            features(eid, 1) = std::max(uDeg, vDeg);
            features(eid, 2) = uDeg + vDeg;
            features(eid, 3) = std::abs(uDeg - vDeg)*std::abs(uDeg - vDeg);


            // cycle based
            // find a nodes n which are connected to u and v
            


            features(eid, 4) = intersectionSize;
        }
    }

    /**
     * @brief      do ucm tranformation of an edge indicator
     *
     * @param      edgeIndicator  0 mean no edge, high value mean 'more edge'
     * @param      features       out
     */
    void ucmTransform(
        vigra::MultiArrayView<2, float> & edgeIndicators,
        vigra::MultiArrayView<1, float> & edgeSizes,   
        vigra::MultiArrayView<1, float> & nodeSizes,                     
        vigra::MultiArrayView<2, float> & features             
    )const{

        typedef MergeGraphAdaptor<AdjacencyListGraph> Mg;

        typedef ArrayMap<float, Edge> EdgeMapView;
        typedef ArrayMap<float, Node> NodeMapView;


        typedef cluster_operators::EdgeWeightedUcm<
            Mg,EdgeMapView,EdgeMapView,NodeMapView,EdgeMapView
        > ClusterOperator;


        typedef HierarchicalClustering<ClusterOperator> Hc;
        typedef typename Hc::Parameter HcParam;

        // buffers
        vigra::MultiArray<1, float > edgeIndicatorBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > ucmBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > edgeSizeBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > nodeSizeBuffer(nodeSizes.shape());

        EdgeMapView edgeIndicatorMap(graph_, edgeIndicatorBuffer);
        EdgeMapView edgeSizeMap(graph_, edgeSizeBuffer);
        NodeMapView nodeSizeMap(graph_, nodeSizeBuffer);
        EdgeMapView ucmMap(graph_, ucmBuffer);



        Mg mg(graph_);
        ClusterOperator cOp(mg, edgeIndicatorMap, edgeSizeMap,
                        nodeSizeMap, ucmMap, 1.0);

        HcParam hcParam(1,false,true);

        for(size_t i=0; i<edgeIndicators.shape(1); ++i){
            
            if(i>0)
                cOp.resetMgAndPq();
            // fill buffers with data
            edgeIndicatorBuffer = edgeIndicators.bindInner(i);
            ucmBuffer =  0.0f;
            edgeSizeBuffer = edgeSizes;
            nodeSizeBuffer = nodeSizes;

            Hc hc(cOp, hcParam);
            hc.cluster();


        }

    }


    UInt64 nCyclePropergationFeatures(

    )const{
        return 8;
    }

    void cyclePropergationFeatures(
        const vigra::MultiArrayView<1, float> & edgeFeatureIn,                    
        vigra::MultiArrayView<2, float> & features             
    )const{


        typedef acc::Select< 
            acc::DataArg<1>,
            acc::Mean, 
            AQuants
        >  SelectType;    

        typedef acc::StandAloneAccumulatorChain<DIM, float, SelectType> FreeChain;
        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        //typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;


        EdgeChainMap edgeAccChain(graph_);
        const size_t nPasses = edgeAccChain[*EdgeIt(graph_)].passesRequired();


        MultiArray<1, TinyVector<Int32, 3> >  cyclesArray;
        find3CyclesEdges(graph_, cyclesArray);


        TinyVector<Int32, DIM> fakeCoord;

        for(size_t p=0; p<nPasses; ++p){
            for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){
                const Edge edge = *eIt;
                const UInt32 eid = graph_.id(edge);
                edgeAccChain[eid].updatePassN(edgeFeatureIn[eid], fakeCoord, p+1); 
            }
            for(size_t ci=0; ci<cyclesArray.shape(0); ++ci){
                const TinyVector<Int32, 3> & ce = cyclesArray[ci];


                uint  indices[3]={0,1,2};
                float vals[3] ={
                    edgeFeatureIn[ce[0]],
                    edgeFeatureIn[ce[1]],
                    edgeFeatureIn[ce[2]]
                };

                indexSort(vals,vals+3, indices);

                for(size_t i=0; i<3; ++i){
                    const UInt64 ei = ce[i];
                    FreeChain & eChain = edgeAccChain[ei];


                    //for(size_t j=0; j<3; ++j){
                    //    if(i!=j){
                    //        eChain.updatePassN(edgeFeatureIn[ce[j]], fakeCoord, p+1);
                    //    }
                    //}
                    // smallest value
                    if(i==indices[0]){
                       // do nothing
                    }
                    else{
                       if(i==indices[1]){
                           UInt64 oe = ce[indices[2]];
                           eChain.updatePassN(edgeFeatureIn[oe], fakeCoord, p+1);
                       }
                       if(i==indices[2]){
                           UInt64 oe = ce[indices[1]];
                           eChain.updatePassN(edgeFeatureIn[oe], fakeCoord, p+1);
                       }
                    }
                }
            }
        }

        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){
            const Edge edge = *eIt;
            const UInt32 eid = graph_.id(edge);

            FreeChain & eChain = edgeAccChain[eid];
            const TinyVector<float, 7> eQnt = acc::get<AQuants>(eChain);
            features(eid,0) = acc::get<acc::Mean>(eChain);
            for(size_t i=0; i<7; ++i){
                features(eid,i+1) = eQnt[i];
            }
        }



    }


    // template<class OUT_TYPE>
    // void topologicalFeatures(
    //     vigra::MultiArrayView<2, OUT_TYPE> & features
    // )const{

    // }

private:




    const AdjacencyListGraph & graph_;
    MultiArrayView<DIM, LABEL_TYPE> labels_;        
};


}


#endif /* VIGRA_GRAPH_FEATURES_HXX */
