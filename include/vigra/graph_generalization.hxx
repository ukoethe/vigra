/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 Stefan Schmidt and Ullrich Koethe            */
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
 * This header provides definitions of graph-related types
 * and optionally provides a gateway to popular graph libraries
 * (for now, BGL is supported).
 */

#ifndef VIGRA_GRAPH_GENERALIZATION_HXX
#define VIGRA_GRAPH_GENERALIZATION_HXX


#include <vigra/graphs.hxx>
#include <vigra/multi_gridgraph.hxx>
namespace vigra{



    template<class MAP>
    struct GraphMapTypeTraits{
        typedef typename MAP::Value          Value;
        typedef typename MAP::Reference      Reference;
        typedef typename MAP::ConstReference ConstReference;
    };

    template<class GRAPH>
    struct GraphIteratorAccessor{
        typedef GRAPH Graph;
        typedef typename Graph::NodeIt NodeIt;
        typedef typename Graph::EdgeIt EdgeIt;
        typedef typename Graph::ArcIt  ArcIt;

        static NodeIt nodesBegin(const Graph & g){ return NodeIt(g);}
        static EdgeIt edgesBegin(const Graph & g){ return EdgeIt(g);}
        static ArcIt  arcsBegin( const Graph & g){ return ArcIt( g);}

        static NodeIt nodesEnd(const Graph & g){ return NodeIt(lemon::INVALID);}
        static EdgeIt edgesEnd(const Graph & g){ return EdgeIt(lemon::INVALID);}
        static ArcIt  arcsEnd( const Graph & g){ return ArcIt( lemon::INVALID);}
    };


    template<unsigned int DIM,class DIRECTED_TAG>
    struct GraphIteratorAccessor<GridGraph<DIM,DIRECTED_TAG> >{
        typedef GridGraph<DIM,DIRECTED_TAG> Graph;
        typedef typename Graph::NodeIt NodeIt;
        typedef typename Graph::EdgeIt EdgeIt;
        typedef typename Graph::ArcIt  ArcIt;

        static NodeIt nodesBegin(const Graph & g){ return NodeIt(g);}
        static EdgeIt edgesBegin(const Graph & g){ return g.get_edge_iterator();}
        static ArcIt  arcsBegin( const Graph & g){ return ArcIt( g);}

        static NodeIt nodesEnd(const Graph & g){ return g.get_vertex_end_iterator();}
        static EdgeIt edgesEnd(const Graph & g){ return g.get_edge_end_iterator();}
        static ArcIt  arcsEnd( const Graph & g){ return g.get_arc_end_iterator(); }
    };



    template<class GRAPH>
    class IntrinsicGraphShape{
    private:
        typedef GRAPH Graph;
    typedef typename Graph::index_type  index_type;
    public:
        typedef typename Graph::Node Node ;
        typedef typename Graph::Edge Edge ;
        typedef typename  Graph::Arc  Arc ;

        typedef TinyVector<Int64,1> IntrinsicNodeMapShape;
        typedef TinyVector<Int64,1> IntrinsicEdgeMapShape;
        typedef TinyVector<Int64,1>  IntrinsicArcMapShape;

        static IntrinsicNodeMapShape intrinsicNodeMapShape(const Graph & g){
            return IntrinsicNodeMapShape(g.maxNodeId()+1);
        }
        static IntrinsicEdgeMapShape intrinsicEdgeMapShape(const Graph & g){
            return IntrinsicEdgeMapShape(g.maxEdgeId()+1);
        }
        static IntrinsicArcMapShape intrinsicArcMapShape(const Graph & g){
            return  IntrinsicArcMapShape(g.maxArcId()+1);
        }


        static const unsigned int IntrinsicNodeMapDimension=1;
        static const unsigned int IntrinsicEdgeMapDimension=1;
        static const unsigned int IntrinsicArceMapDimension=1;
    };

    template<unsigned int DIM,class DIRECTED_TAG>
    class IntrinsicGraphShape<GridGraph<DIM,DIRECTED_TAG> >{
    private:
        typedef GridGraph<DIM,DIRECTED_TAG> Graph;
    typedef typename Graph::index_type  index_type;
    public:
        typedef typename Graph::Node Node ;
        typedef typename Graph::Edge Edge ;
        typedef typename  Graph::Arc  Arc ;

        typedef typename Graph::shape_type              IntrinsicNodeMapShape;
        typedef typename Graph::edge_propmap_shape_type IntrinsicEdgeMapShape;
        typedef typename Graph::edge_propmap_shape_type IntrinsicArcMapShape;

        static IntrinsicNodeMapShape intrinsicNodeMapShape(const Graph & g){
            return  g.shape();
        }
        static IntrinsicEdgeMapShape intrinsicEdgeMapShape(const Graph & g){
            return  g.edge_propmap_shape();
        }
        static IntrinsicArcMapShape intrinsicArcMapShape(const Graph & g){
            return  g.arc_propmap_shape();
        }
        static const unsigned int IntrinsicNodeMapDimension=DIM;
        static const unsigned int IntrinsicEdgeMapDimension=DIM+1;
        static const unsigned int IntrinsicArceMapDimension=DIM+1;
    };




    template<class GRAPH>
    class GraphDescriptorToMultiArrayIndex{
    private:
        typedef GRAPH Graph;
    typedef typename Graph::index_type  index_type;
    public:
        typedef typename Graph::Node Node ;
        typedef typename Graph::Edge Edge ;
        typedef typename  Graph::Arc  Arc ;

        typedef typename IntrinsicGraphShape<Graph>::IntrinsicNodeMapShape IntrinsicNodeMapShape;
        typedef typename IntrinsicGraphShape<Graph>::IntrinsicEdgeMapShape IntrinsicEdgeMapShape;
        typedef typename IntrinsicGraphShape<Graph>::IntrinsicArcMapShape  IntrinsicArcMapShape;

        static IntrinsicNodeMapShape intrinsicNodeCoordinate(const Graph & g,const Node & node){
            return IntrinsicNodeMapShape(g.id(node));
        }
        static IntrinsicEdgeMapShape intrinsicEdgeCoordinate(const Graph & g,const Edge & edge){
            return IntrinsicEdgeMapShape(g.id(edge));
        }
        static IntrinsicArcMapShape intrinsicArcCoordinate(const Graph & g,const Arc & arc){
            return  IntrinsicArcMapShape(g.id(arc));
        }

    };



    template<unsigned int DIM,class DIRECTED_TAG>
    class GraphDescriptorToMultiArrayIndex<GridGraph<DIM,DIRECTED_TAG> >{
    private:
        typedef GridGraph<DIM,DIRECTED_TAG> Graph;
        typedef typename Graph::index_type  index_type;
    public:
        typedef typename Graph::Node Node ;
        typedef typename Graph::Edge Edge ;
        typedef typename  Graph::Arc  Arc ;

        typedef typename IntrinsicGraphShape<Graph>::IntrinsicNodeMapShape IntrinsicNodeMapShape;
        typedef typename IntrinsicGraphShape<Graph>::IntrinsicEdgeMapShape IntrinsicEdgeMapShape;
        typedef typename IntrinsicGraphShape<Graph>::IntrinsicArcMapShape  IntrinsicArcMapShape;


        static Node intrinsicNodeCoordinate(const Graph & g,const Node & node){
            return node;
        }
        static Edge intrinsicEdgeCoordinate(const Graph & g,const Edge & edge){
            return edge;
        }
        static Arc  intrinsicArcCoordinate (const Graph & g,const Arc & arc){
            return arc;
        }

    };





} // namespace vigra

#endif // VIGRA_GRAPH_GENERALIZATION_HXX
