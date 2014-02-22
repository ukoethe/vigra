/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier Ullrich Koeth                   */
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

#ifndef VIGRA_GRAPH_MAP_ALGORITHMS_HXX
#define VIGRA_GRAPH_MAP_ALGORITHMS_HXX

/*std*/
#include <algorithm>
#include <vector>
#include <functional>


/*vigra*/
#include <vigra/graphs.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/priority_queue.hxx>
#include <vigra/union_find.hxx>
#include <vigra/adjacency_list_graph.hxx>
#include <vigra/graph_helper/on_the_fly_edge_map.hxx>



namespace vigra{
    namespace detail_graph_morphology{



        template<class T>
        struct MinFunctor{
            MinFunctor(const T & val)
            :   val_(val){
            }
            void neighbourValue(const T & val){
                val_=std::min(val,val_);
            }
            const T & resultValue()const{
                return val_;
            }
            T val_;
        };


        template<class MorphologyFunctor,class Graph, class NodeMapIn,class NodeMapOut>
        void graphMorphologyImpl(
            const Graph & g,
            const NodeMapIn & nodeMapIn,
            NodeMapOut &      nodeMapOut
        ){  
            typedef typename Graph::Node        Node;
            typedef typename Graph::NodeIt      NodeIt;
            typedef typename Graph::OutArcIt    OutArcIt;
            for(NodeIt n(g); n!=lemon::INVALID; ++n){
                MorphologyFunctor functor(nodeMapIn[*n]);
                for(OutArcIt a(g,n);a!=lemon::INVALID;++n){
                    functor.neighbourValue(nodeMapIn[ g.target(*a) ]);
                }
                nodeMapOut[*n]=functor.resultValue();
            }
        }
    }




    template<class Graph, class NodeMapIn,class NodeMapOut>
    void graphMorphologyImpl(
        const Graph & g,
        const NodeMapIn & nodeMapIn,
        NodeMapOut &      nodeMapOut
    ){  
        typedef detail_graph_morphology::MinFunctor<typename NodeMapIn::Value> Functor;
        graphMorphologyImpl<Functor>(g,nodeMapIn,nodeMapOut);
    }

}


#endif // VIGRA_GRAPH_MAP_ALGORITHMS_HXX
