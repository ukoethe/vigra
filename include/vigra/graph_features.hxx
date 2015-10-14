#ifndef VIGRA_GRAPH_FEATURES_HXX
#define  VIGRA_GRAPH_FEATURES_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/for_each_coord.hxx>




namespace vigra{




template< 
    unsigned int DIM,class LABELS_ARRAY,
    class DATA_ARRAY,class NODE_ACC_CHAIN,class EDGE_ACC_CHAIN
>
void nodeAndEdgeAccumlation(
    const AdjacencyListGraph &  graph,
    const LABELS_ARRAY &        labelsArray,
    const DATA_ARRAY &          dataArray,
    NODE_ACC_CHAIN &            nodeAccChain, 
    EDGE_ACC_CHAIN &            edgeAccChain
){
    typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
    typedef typename AdjacencyListGraph::NodeIt NodeIt;
    typedef typename MultiArrayShape<DIM>::type Coord;
    const size_t nPassesNodeAcc = nodeAccChain[*NodeIt(graph)].passesRequired();
    const size_t nPassesEdgeAcc = edgeAccChain[*EdgeIt(graph)].passesRequired();
    for(size_t p=0; p<std::max(nPassesNodeAcc, nPassesEdgeAcc); ++p){
        ForEachCoord<DIM>::forEachCoord(labelsArray.shape(), [&](const Coord & coord){
            const auto uLabel = labelsArray[coord];
            const auto uNode = graph.nodeFromId(labelsArray[coord]);
            const auto uVal = dataArray[coord];
            if(p<nPassesNodeAcc){
                nodeAccChain[uNode].updatePassN(uVal, coord, p+1);
            }
            if(p<nPassesEdgeAcc){
                for(size_t d=0; d<DIM; ++d){
                    Coord otherCoord = coord;
                    otherCoord[d] += 1;
                    if(otherCoord[d] < labelsArray.shape(d)){
                        const auto vNode = graph.nodeFromId(labelsArray[otherCoord]);
                        if(uNode != vNode){
                            const auto vVal = dataArray[otherCoord];
                            const auto e = graph.findEdge(uNode, vNode);
                            edgeAccChain[e].updatePassN(uVal, coord, p+1); 
                            edgeAccChain[e].updatePassN(vVal, otherCoord, p+1); 
                        }
                    }
                }
            }
        });
    }
}

















template<unsigned int DIM, class LABEL_TYPE>
class GridRagFeatureExtractor{

    typedef acc::AutoRangeHistogram<40> AHist;
    typedef acc::StandardQuantiles<AHist> Quants;
    typedef typename AdjacencyListGraph::Edge Edge;
    typedef typename AdjacencyListGraph::Node Node;
    typedef typename AdjacencyListGraph::EdgeIt EdgeIt;
    typedef typename AdjacencyListGraph::NodeIt NodeIt;
    typedef typename AdjacencyListGraph::OutArcIt OutArcIt;
    typedef std::set<Node> NodeSet;
    typedef typename NodeSet::const_iterator NodeSetIter;


    struct ZeroVal{
        template<class T>
        int operator[](const T &)const{
            return int();
        }
    };

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
    UInt64 nAccumulatedFeatures()const{
        return 9*11;
    }

    template<
        class DATA_TYPE,
        class OUT_TYPE
    >
    std::vector<std::string> accumulatedFeatures(
        const vigra::MultiArrayView<DIM, DATA_TYPE> & data,
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


        std::vector<std::string> fnames;
        std::vector<std::string> baseNames = {
            "mean","StdDev","Minimum", "Q0.10", "Q0.25","Q0.50", "Q0.75", "Q0.9","Maximum"
        };
        std::vector<std::string> defNames = {
            "E","UVSum","absDiffUV","minUV","maxUV","meanUV",
            "absEminUVDiff","absEmaxUVDiff","absEmeanUVDiff",
            "min(eUE,eVE) - dUV","max(eUE,eVE) - dUV"
        };

        for(const auto & bn : baseNames){
            for(const auto & dn: defNames){
                fnames.push_back(bn +std::string("_")+ dn);
            }
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
        return fnames;
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
    std::vector<std::string> 
    geometricFeatures(
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



        std::vector<std::string> fnames{
            "eCount",
            "minUVCount",
            "maxUVCount",
            "sumUVCount",
            "minUVCount/maxUvCount",
            "absUVCountDiff",
            "eNCount/minNUVCount",
            "eNCount/maxNUVCount",
            "eNCount/NUVCount",

            "dUV",
            "mindEUdEV",
            "maxdEUdEV",
            "sumdEUdEV",
            "mindEUdEV/maxdEUdEV",
            "absdEUdEVDiff",

            "uvNCount/dEV",
            "min(uRat, vRat)",
            "max(uRat, vRat)",
            "uRat + vRat",
            "abs(uRat - vRat)+0.5",
            "max(uRat, vRat) / std::min(uRat, vRat)"
        };

        UInt64 fi = 0;
        nodeAndEdgeAccumlation<DIM>(graph_, labels_, ZeroVal(),nodeAccChainMap, edgeAccChainMap);
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



            const float uRat = uNCount/(dEU+1.0);
            const float vRat = vNCount/(dEV+1.0);
            const float uvRat = uvNCount/(dUV+1.0);
            
/*            std::cout<<"dEU"<<dEU<<"\n";
            std::cout<<"dEV"<<dEV<<"\n";
            std::cout<<"dUV"<<dUV<<"\n";
            std::cout<<"uNCount"<<uNCount<<"\n";
            std::cout<<"vNCount"<<vNCount<<"\n";
            std::cout<<"eNCount"<<eNCount<<"\n";*/
            
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
            features(eid,15) =uvRat;
            features(eid,16) =std::min(uRat, vRat);
            features(eid,17) =std::max(uRat, vRat);
            features(eid,18) =uRat + vRat;
            features(eid,19) =std::abs(uRat - vRat)+0.5;
            features(eid,20) =std::max(uRat, vRat) / std::min(uRat, vRat);
            
        }
        return fnames;
    }


    UInt64 nTopologicalFeatures()const{
        return 5;
    }

    template<class OUT_TYPE>
    std::vector<std::string> 
    topologicalFeatures(vigra::MultiArrayView<2, OUT_TYPE> & features )const{


        std::vector<std::string> fnames{
            "minUVDegree",
            "maxUVDegree",
            "sumUVDegree",
            "absUVDegreeDiff",
            "uvNhIntersectionSize"
        };
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
            features(eid, 3) = std::abs(uDeg - vDeg);


            // cycle based
            // find a nodes n which are connected to u and v
            


            features(eid, 4) = intersectionSize;
        }
        return fnames;
    }

     
    UInt64 nUcmTransformFeatures(
        const vigra::MultiArrayView<2, float> & edgeIndicators,
        const vigra::MultiArrayView<1, float> & wardness
    )const{
        return edgeIndicators.shape(1)*2*wardness.shape(0);
    }

    void ucmTransformFeatures(
        const vigra::MultiArrayView<2, float> & edgeIndicators,   
        const vigra::MultiArrayView<1, float> & wardness,                 
        vigra::MultiArrayView<2, float> & features             
    )const{

        typedef acc::Select< 
            acc::Count
        >  SelectType;    

        typedef acc::StandAloneDataFreeAccumulatorChain<DIM, SelectType> FreeChain;
        typedef typename AdjacencyListGraph:: template EdgeMap<FreeChain> EdgeChainMap;
        typedef typename AdjacencyListGraph:: template NodeMap<FreeChain> NodeChainMap;




    



        typedef MergeGraphAdaptor<AdjacencyListGraph> Mg;
        typedef ArrayMap<float,AdjacencyListGraph, Edge> EdgeMapView;
        typedef ArrayMap<float,AdjacencyListGraph, Node> NodeMapView;


        typedef cluster_operators::EdgeWeightedUcm<
            Mg,EdgeMapView,EdgeMapView,NodeMapView,EdgeMapView
        > ClusterOperator;


        typedef HierarchicalClustering<ClusterOperator> Hc;
        typedef typename Hc::Parameter HcParam;


        NodeChainMap nodeAccChainMap(graph_);
        EdgeChainMap edgeAccChainMap(graph_);

        nodeAndEdgeAccumlation<DIM>(graph_, labels_, ZeroVal(), nodeAccChainMap, edgeAccChainMap);


        //std::cout<<"a\n";

        vigra::MultiArray<1, float > edgeSizes(vigra::MultiArray<1, float >::difference_type(graph_.edgeNum()));
        vigra::MultiArray<1, float > nodeSizes(vigra::MultiArray<1, float >::difference_type(graph_.maxNodeId()+1));

        //std::cout<<"b\n";
        
        for(EdgeIt eIt(graph_); eIt != lemon::INVALID; ++eIt){
            const Edge edge  = *eIt;
            const FreeChain & eChain = edgeAccChainMap[edge];
            const UInt32 eid = graph_.id(edge);
            edgeSizes[eid] = acc::get<acc::Count>(eChain);
        }

        //std::cout<<"c\n";

        for(NodeIt nIt(graph_); nIt != lemon::INVALID; ++nIt){
            const Node node  = *nIt;
            const FreeChain & nChain = nodeAccChainMap[node];
            const UInt32 nid = graph_.id(node);
            nodeSizes[nid] = acc::get<acc::Count>(nChain);
        }

        //std::cout<<"d\n";

        // buffers
        vigra::MultiArray<1, float > edgeIndicatorBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > ucmBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > edgeSizeBuffer(edgeSizes.shape());
        vigra::MultiArray<1, float > nodeSizeBuffer(nodeSizes.shape());

        //std::cout<<"e\n";

        EdgeMapView edgeIndicatorMap(graph_, edgeIndicatorBuffer);
        EdgeMapView edgeSizeMap(graph_, edgeSizeBuffer);
        NodeMapView nodeSizeMap(graph_, nodeSizeBuffer);
        EdgeMapView ucmMap(graph_, ucmBuffer);


        //std::cout<<"f\n";

        edgeIndicatorBuffer = edgeIndicators.bindOuter(0);
        edgeSizeBuffer = edgeSizes;
        nodeSizeBuffer = nodeSizes;


        //std::cout<<"g\n";

        Mg mg(graph_);
        ClusterOperator cOp(mg, edgeIndicatorMap, edgeSizeMap,
                        nodeSizeMap, ucmMap, 0.0);

        HcParam hcParam(1,false,false);

        //std::cout<<"mg node num "<<mg.nodeNum()<<"\n";
        //std::cout<<"mg edge num "<<mg.edgeNum()<<"\n";


        size_t fi=0;

        for(auto w : wardness){
            for(size_t i=0; i<edgeIndicators.shape(1); ++i){

                edgeIndicatorBuffer = edgeIndicators.bindOuter(i);
                ucmBuffer =  0.0f;
                edgeSizeBuffer = edgeSizes;
                nodeSizeBuffer = nodeSizes;


    
                cOp.setWardness(w);
                cOp.resetMgAndPq();
                
                // fill buffers with data
                
                Hc hc(cOp, hcParam);
                hc.cluster();


                // ucm map itself
                hc.ucmTransform(ucmMap);
                features.bindOuter(fi++) = ucmBuffer;
                //  edge indicator
                hc.ucmTransform(edgeIndicatorMap);
                features.bindOuter(fi++) = edgeIndicatorBuffer;
            }
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
            Quants
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
            const TinyVector<float, 7> eQnt = acc::get<Quants>(eChain);
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
