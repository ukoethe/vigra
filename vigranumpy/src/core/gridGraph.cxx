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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpygraphs_PyArray_API
#define NO_IMPORT_ARRAY

#define WITH_BOOST_GRAPH

/*vigra*/
#include "export_graph_visitor.hxx"
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_algorithms.hxx>
#include <vigra/python_graph_generalization.hxx>
#include <vigra/graph_helper/on_the_fly_edge_map.hxx>
#include <vigra/merge_graph/distance.hxx>

/*boost*/
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>

namespace python = boost::python;

namespace vigra{


    template<unsigned int DIM,class DTAG,class LABEL_TYPE>
    void pyGetRag(
        const GridGraph<DIM,DTAG>                                                                & graph,
        NumpyArray<DIM,LABEL_TYPE>                                                                labelMap,
        AdjacencyListGraph                                                                       & rag,
        AdjacencyListGraph:: template EdgeMap< std::vector<typename GridGraph<DIM,DTAG>::Edge> > & hyperEdges,
        const Int64                                                                                ignoreLabel=-1
    ){
        makeRegionAdjacencyGraph(graph,labelMap,rag,hyperEdges,ignoreLabel);
    }



    template<unsigned int DIM,class DTAG,class DIST_FUNCTOR>
    NumpyAnyArray pyGridGraphFelzenszwalbSegmentationT(
        NumpyArray<DIM+1,Multiband<float>  >   image,
        NumpyArray<DIM  ,float >               edgeIndicatorImage,
        const float                           k, 
        NumpyArray<DIM,UInt32>                labeling = NumpyArray<DIM,UInt32>()
    ){
        typename NumpyArray<DIM,UInt32>::difference_type outShape;
        for(size_t i=0;i<DIM;++i){
            outShape[i]=image.shape(i);
        }
        //resize output? 
        labeling.reshapeIfEmpty(outShape);

        // make a grid graph
        typedef GridGraph< DIM,boost::undirected_tag > Graph;
        typedef typename Graph::Node Node;
        typedef NumpyMultibandNodeMap<Graph,  NumpyArray<DIM+1,Multiband<float> >  > NodeMultiFloatMap;
        typedef NumpyScalarNodeMap<Graph,  NumpyArray<DIM,float> > NodeFloatMap;
        typedef NumpyScalarNodeMap<Graph,  NumpyArray<DIM,UInt32> > NodeLabelMap;
        typedef distances::Manhattan<float> DistanceType;
        typedef OnTheFlyEdgeMap<Graph,NodeMultiFloatMap,DistanceType,float > ImplicitNodeDiffMap;
        typedef std::plus<float> AccType;
        typedef OnTheFlyEdgeMap<Graph,NodeFloatMap,AccType,float > ImplicitEdgeIndicatorMap;
        typedef std::plus<float> BinOp;
        typedef BinaryOpEdgeMap<Graph,ImplicitEdgeIndicatorMap,ImplicitNodeDiffMap,BinOp,float> CombinedEdgeWeightMap;

        // make grid graph
        const Graph g(outShape);

        // make a node map from image 
        // and make implicit edge map from that
        NodeMultiFloatMap imageMap(g,image);
        DistanceType distanceFunctor;
        ImplicitNodeDiffMap nodeDistEdgeMap(g,imageMap,distanceFunctor);

        // make a node map from edge edge indicator image
        // and make implicit edge map from that
        NodeFloatMap  edgeIndicatorImageNodeMap(g,edgeIndicatorImage);
        AccType accFunctor;
        ImplicitEdgeIndicatorMap edgeIndicatorEdgeMap(g,edgeIndicatorImageNodeMap,accFunctor);


        // combine both maps
        BinOp binOpFunctor;
        CombinedEdgeWeightMap combinedEdgeWeightMap(g,edgeIndicatorEdgeMap,nodeDistEdgeMap,binOpFunctor);

        // make the result label map from image
        NodeLabelMap nodeLabelingMap(g,labeling);

        // make sizemap filled with ones
        typename  Graph::template NodeMap<size_t>  nodeSizeMap(g,1.0);

        // run alg.
        felzenszwalbSegmentation(g,combinedEdgeWeightMap,nodeSizeMap,k,nodeLabelingMap);
        return labeling;
    }


    template<unsigned int DIM,class DTAG>
    NumpyAnyArray pyGridGraphFelzenszwalbSegmentation(
        NumpyArray<DIM+1,Multiband<float>  > image,
        NumpyArray<DIM  ,float >             edgeIndicatorImage,
        const float                          k, 
        const std::string &                  distanceType,
        NumpyArray<DIM,UInt32>               labeling = NumpyArray<DIM,UInt32>()
    ){
        if      (distanceType=="l1" || distanceType == "manhatten")
            return pyGridGraphFelzenszwalbSegmentationT<DIM,DTAG,distances::Manhattan<float> >(
                image,edgeIndicatorImage,k,labeling);
        else if (distanceType=="l2" || distanceType == "eucledian")
            return pyGridGraphFelzenszwalbSegmentationT<DIM,DTAG,distances::Norm<float> >(
                image,edgeIndicatorImage,k,labeling);
        else if (distanceType=="chiSquared")
            return pyGridGraphFelzenszwalbSegmentationT<DIM,DTAG,distances::ChiSquared<float> >(
                 image,edgeIndicatorImage,k,labeling);
        else{
            throw std::runtime_error("unknown distance type");
            return labeling;
        } 

            
    }

    template<unsigned int DIM,class DTAG>
    NumpyAnyArray pyGridGraphMinimumSpanningTreeSegmentation(
        NumpyArray<DIM+1,Multiband<float> >   image,
        const bool                            useWeightThreshold,
        const bool                            useNodeNumThreshold,
        const float                           weightThreshold,
        const size_t                          nodeNumThreshold, 
        NumpyArray<DIM,UInt32>                labeling = NumpyArray<DIM,UInt32>()
    ){
        typename NumpyArray<DIM,UInt32>::difference_type outShape;
        for(size_t i=0;i<DIM;++i){
            outShape[i]=image.shape(i);
        }
        //resize output? 
        labeling.reshapeIfEmpty(outShape);

        // make a grid graph
        typedef GridGraph< DIM,boost::undirected_tag > Graph;
        typedef typename Graph::Node Node;
        typedef NumpyMultibandNodeMap<Graph,  NumpyArray<DIM+1,Multiband<float> >  > NodeMultiFloatMap;
        typedef NumpyScalarNodeMap<Graph,  NumpyArray<DIM,UInt32> > NodeLabelMap;
        typedef distances::Manhattan<float> DistanceType;
        typedef OnTheFlyEdgeMap<Graph,NodeMultiFloatMap,DistanceType,float > ImplicitEdgeMap;

        // make grid graph
        const Graph g(outShape);

        // make a node map from edgeIndicator image
        // and make implicit edge map from that
        NodeMultiFloatMap imageMap(g,image);
        DistanceType distanceFunctor;
        ImplicitEdgeMap edgeIndicatorEdgeMap(g,imageMap,distanceFunctor);

        // make the result label map from image
        NodeLabelMap nodeLabelingMap(g,labeling);

        // run alg.
        minimumSpanningTreeSegmentation(g,edgeIndicatorEdgeMap,useWeightThreshold,useNodeNumThreshold,
            weightThreshold,nodeNumThreshold,nodeLabelingMap);
        return labeling;
    }

    template<unsigned int DIM,class DTAG>
    NumpyAnyArray pyGridGraphShortestPathDijkstra(
        NumpyArray<DIM,float>   edgeIndicatorImage,
        typename NumpyArray<DIM,float>::difference_type source,
        typename NumpyArray<DIM,float>::difference_type target
    ){
        // make a grid graph
        typedef GridGraph< DIM,boost::undirected_tag > Graph;
        typedef typename Graph::Node Node;
        typedef typename Graph::vertex_descriptor vertex_descriptor;
        typedef NumpyScalarNodeMap<Graph,  NumpyArray<DIM,float> > NodeFloatMap;
        typedef OnTheFlyEdgeMap<Graph,NodeFloatMap,std::plus<float> ,float> ImplicitEdgeMap;


        // make a node map from edgeIndicator image
        const Graph g(edgeIndicatorImage.shape());
        NodeFloatMap edgeIndicatorNodeMap(g,edgeIndicatorImage);
        std::plus<float> functor;
        ImplicitEdgeMap edgeIndicatorEdgeMap(g,edgeIndicatorNodeMap,functor);


        // make a distance map and predecessor map
        typedef typename Graph:: template NodeMap<vertex_descriptor> PredecessorMapType;
        typedef typename Graph:: template NodeMap<float> DistMapType;
        typedef typename Graph:: template NodeMap<int>   IndexMapType;

        //typedef NumpyArray<DIM,float>  b

        PredecessorMapType pmap(g);
        DistMapType        dmap(g);

        const Node s(source); 
        const Node t(target); 

        shortestPathDijkstra(g,s,edgeIndicatorEdgeMap,pmap,dmap,t);
        std::vector<Node> path;

        // has the target node been reached?
        if(pmap[t]==lemon::INVALID){
            NumpyArray<2,UInt32>  pathCooordinates(typename NumpyArray<2,UInt32>::difference_type(0,0));
            return pathCooordinates;
        }
        else{
            Node running = t;
            path.push_back(t);
            while(running!=s){
                running=pmap[running];
                path.push_back(running);
            }
            NumpyArray<2,UInt32>  pathCooordinates(typename NumpyArray<2,UInt32>::difference_type(path.size(),DIM));
            for(size_t i=0;i<path.size();++i){
                for(size_t d=0;d<DIM;++d)
                    pathCooordinates(path.size()-1-i,d)=path[i][d];
            }
            return pathCooordinates;
        }
    }



    template<unsigned int DIM,class DTAG>
    NumpyAnyArray pyGridGraphShortestPathDijkstraDistanceImage(
        NumpyArray<DIM,float>   edgeIndicatorImage,
        typename NumpyArray<DIM,float>::difference_type source,
        NumpyArray<DIM,float>   distanceImage = NumpyArray<DIM,float>()
    ){

        //resize output? 
        distanceImage.reshapeIfEmpty(edgeIndicatorImage.shape());

        // make a grid graph
        typedef GridGraph< DIM,boost::undirected_tag > Graph;
        typedef typename Graph::Node Node;
        typedef typename Graph::vertex_descriptor vertex_descriptor;
        typedef NumpyScalarNodeMap<Graph,  NumpyArray<DIM,float> > NodeFloatMap;
        typedef OnTheFlyEdgeMap<Graph,NodeFloatMap,std::plus<float> ,float> ImplicitEdgeMap;

        // make grid graph
        const Graph g(edgeIndicatorImage.shape());

        // make a node map from edgeIndicator image
        // and make implicit edge map from that
        NodeFloatMap edgeIndicatorNodeMap(g,edgeIndicatorImage);
        std::plus<float> functor;
        ImplicitEdgeMap edgeIndicatorEdgeMap(g,edgeIndicatorNodeMap,functor);


        // make a distance map from distance image
        NodeFloatMap dmap(g,distanceImage);

        // make a distance map and predecessor map
        typedef typename Graph:: template NodeMap<vertex_descriptor> PredecessorMapType;
        typedef typename Graph:: template NodeMap<float> DistMapType;
        typedef typename Graph:: template NodeMap<int>   IndexMapType;

        //typedef NumpyArray<DIM,float>  b

        PredecessorMapType pmap(g);
        //DistMapType        dmap(g);

        const Node s(source); 
        shortestPathDijkstra(g,s,edgeIndicatorEdgeMap,pmap,dmap);
        return distanceImage;
    }

    template<unsigned int DIM>
    void defineGridGraphT(const std::string & clsName){

        typedef GridGraph<DIM,boost::undirected_tag> Graph;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::Node Node;

        typedef typename AdjacencyListGraph:: template EdgeMap< std::vector<Edge> > HyperEdgeMap;
        typedef typename AdjacencyListGraph:: template NodeMap< std::vector<Node> > HyperNodeMap;


        const std::string hyperEdgeMapNamClsName = clsName + std::string("HyperEdgeMap");
        python::class_<HyperEdgeMap>(hyperEdgeMapNamClsName.c_str(),python::init<>())
        ;



        
        python::class_<Graph>(clsName.c_str(),python::init< TinyVector<Int64,DIM> >())
        .def(LemonUndirectedGraphCoreVisitor<Graph>(clsName))
        .def("getRegionAdjacencyGraph",registerConverters(&pyGetRag<DIM,boost::undirected_tag,UInt32> ),
           (
               python::arg("labels"),
               python::arg("rag"),
               python::arg("hyperEdges"),
               python::arg("ignoreLabel")=-1
           )
        )
        ;

        python::def("shortestImagePathDijkstra",registerConverters(&pyGridGraphShortestPathDijkstra<DIM,boost::undirected_tag>) );
       
      
        python::def("felzenszwalbSegmentation",registerConverters(&pyGridGraphFelzenszwalbSegmentation<DIM,boost::undirected_tag>),
            (
               python::arg("nodeFeatures"),
               python::arg("edgeWeights"),
               python::arg("k"),
               python::arg("distanceType")=std::string("manhatten"),
               python::arg("out")=python::object()
            )
        );
        python::def("minimumSpanningTreeSegmentation",registerConverters(&pyGridGraphMinimumSpanningTreeSegmentation<DIM,boost::undirected_tag>),
            (
               python::arg("weights"),
               python::arg("useWeightThreshold")=false,
               python::arg("useNodeNumThreshold")=true,
               python::arg("weightThreshold")=0.0,
               python::arg("nodeNumThreshold")=0,
               python::arg("out")=python::object()
            )
        );

        python::def("shortestImagePathDistanceImageDijkstra",registerConverters(&pyGridGraphShortestPathDijkstraDistanceImage<DIM,boost::undirected_tag>),
            (
               python::arg("weights"),
               python::arg("source"),
               python::arg("out")=python::object()
            )
        );
    }

    void defineGridGraph(){
        defineGridGraphT<2>("GridGraphUndirected2d");
        defineGridGraphT<3>("GridGraphUndirected3d");
    }

    void defineInvalid(){
        python::class_<lemon::Invalid>("Invalid",python::init<>())
        ;
    }

} 


