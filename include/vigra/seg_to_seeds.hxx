/************************************************************************/
/*                                                                      */
/*               Copyright 2012-2013 by Ullrich Koethe                  */
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


#ifndef VIGRA_SEG_TO_SEED_HXX
#define VIGRA_SEG_TO_SEED_HXX

#include "multi_array.hxx"
#include "multi_gridgraph.hxx"
namespace vigra {


    template<unsigned int DIM, class LABEL_TYPE,class LABEL_TYPE_OUT>
    void shrinkLabels(
        MultiArrayView<DIM,LABEL_TYPE>     labels,
        const size_t shrinkNpixels,
        MultiArrayView<DIM,LABEL_TYPE_OUT> shrinkedLabels
    ){
        shrinkedLabels = labels;

        typedef GridGraph<DIM,boost::undirected_tag> Graph;
        typedef typename Graph::Node Node;
        typedef typename Graph::Edge Edge;
        typedef typename Graph::NodeIt graph_scanner;
        typedef typename Graph::OutArcIt neighbor_iterator;

        const Graph g(labels.shape());

        // INITAL LOOP
        for (graph_scanner n(g); n != lemon::INVALID; ++n){
            const Node node(*n);
            for (neighbor_iterator arc(g, node); arc != lemon::INVALID; ++arc){
                const Node otherNode = g.target(arc);

                if(labels[node]!=labels[otherNode]){
                    shrinkedLabels[node]=0;
                    shrinkedLabels[otherNode]=0;
                }
            }
        }

        MultiArray<DIM,bool> visited(labels.shape());
        for(size_t r=0;r<shrinkNpixels-1;++r){
            std::fill(visited.begin(),visited.end(),false);
            for (graph_scanner n(g); n != lemon::INVALID; ++n){
                const Node node(*n);
                if(!visited[n] && shrinkedLabels[node]==0){
                    for (neighbor_iterator arc(g, node); arc != lemon::INVALID; ++arc){
                        const Node otherNode = g.target(arc);
                        shrinkedLabels[otherNode]=0;
                        visited[otherNode]=true;
                    }
                }
            }
        }
    }


} // end namespace vigra

#endif // VIGRA_SEG_TO_SEED_HXX
