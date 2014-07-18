/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
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

#ifndef VIGRA_GRAPH_SERIALIZATION
#define VIGRA_GRAPH_SERIALIZATION

/*std*/
#include <algorithm>
#include <vector>
#include <functional>
#include <set>


/*vigra*/
#include "graphs.hxx"
#include "graph_generalization.hxx"
#include "multi_gridgraph.hxx"
#include "priority_queue.hxx"
#include "union_find.hxx"
#include "adjacency_list_graph.hxx"
#include "graph_maps.hxx"





namespace vigra{



    template<class GRAPH>
    class GraphSerialization{
    public:

        typedef GRAPH Graph;

        static UInt64 size(const Graph & g);

        template<class ITER>
        static void serialize(const Graph & g , const ITER begin, const ITER end);

        template<class ITER>
        static void deserialize(const ITER begin, const ITER end, Graph & g );

    };



    template<unsigned int DIM,class DTAG>
    class GraphSerialization< GridGraph<DIM, DTAG> >{
    public:

        typedef GridGraph<DIM, DTAG> Graph;

        static UInt64 size(const Graph & g){
            // shape (DIM)
            // directed tag (even if it is templated) (1)
            // 4-8 nh / 6-26 (1)
            return DIM+1;
        }

        template<class ITER>
        static void serialize(const Graph & g , const ITER begin, const ITER end){
            precondition(std::distance(begin,end)==size(),"wrong iterator begin-end distance");
            for(size_t d=0; d<DIM; ++d){
                *begin=g.shape()[d];
                ++begin;
            }
            precondition(false,"IMPL. ME");
            *begin= IsSameType<DTAG, directed_tag>::value > ? 1 : 0;
            ++begin; // FIXME DIRECTED TAG HERE

            *begin=0; 
            ++begin; // FIXME 4-8 nh 
        }

        template<class ITER>
        static void deserialize(const ITER begin, const ITER end, Graph & g ){
            precondition(std::distance(begin,end)==size(),"wrong iterator begin-end distance");
        }

    };

    template<unsigned int DIM,class DTAG>
    class GraphSerialization< GridGraph<DIM, DTAG> >{
    public:

        typedef AdjacencyListGraph Graph;

        static UInt64 size(const Graph & g){
            return g.serializationSize();
        }

        template<class ITER>
        static void serialize(const Graph & g , const ITER begin, const ITER end){
            g.serialize(begin);
        }

        template<class ITER>
        static void deserialize(const ITER begin, const ITER end, Graph & g ){
            g.deserialize(begin,end);
        }

    };



} // namespace vigra

#endif // VIGRA_GRAPH_SERIALIZATION
