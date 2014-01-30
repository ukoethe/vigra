/************************************************************************/
/*                                                                      */
/*     Copyright 2009-2010 by Ullrich Koethe and Hans Meine             */
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
#ifndef VIGRA_MERGE_GRAPH_PYTHON_CLUSTERING_OPERATOR_HXX
#define VIGRA_MERGE_GRAPH_PYTHON_CLUSTERING_OPERATOR_HXX

/*boost*/

#include <boost/python.hpp>

/*vigra*/
#include <vigra/priority_queue.hxx>
#include <vigra/merge_graph/distance.hxx>

namespace vigra{



    namespace cluster_operators{

    template<class MERGE_GRAPH>
    class PythonOperator{
        
        typedef PythonOperator<MERGE_GRAPH > SelfType;
    public:


        typedef float WeightType;
        typedef MERGE_GRAPH MergeGraph;
        typedef typename MergeGraph::Graph Graph;
        typedef typename Graph::Edge GraphEdge;
        typedef typename Graph::Node GraphNode;
        typedef typename MergeGraph::Edge Edge;
        typedef typename MergeGraph::Node Node;
        typedef typename MergeGraph::EdgeIt EdgeIt;
        typedef typename MergeGraph::NodeIt NodeIt;
        typedef typename MergeGraph::IncEdgeIt IncEdgeIt;
        typedef typename MergeGraph::index_type index_type;
        typedef MergeGraphItemHelper<MergeGraph,Edge> EdgeHelper;
        typedef MergeGraphItemHelper<MergeGraph,Node> NodeHelper;


        typedef NodeHolder<MERGE_GRAPH> NodeHolderType;
        typedef EdgeHolder<MERGE_GRAPH> EdgeHolderType;

        PythonOperator(
            MergeGraph & mergeGraph,
            boost::python::object object,
            const bool useMergeNodeCallback,
            const bool useMergeEdgesCallback,
            const bool useEraseEdgeCallback
        )
        :   mergeGraph_(mergeGraph),
            object_(object)
        {
            if(useMergeNodeCallback)
                mergeGraph_.registerMergeNodeCallBack(*this,& SelfType::mergeNodes);
            if(useMergeEdgesCallback)
                mergeGraph_.registerMergeEdgeCallBack(*this,& SelfType::mergeEdges);
            if(useEraseEdgeCallback)
                mergeGraph_.registerEraseEdgeCallBack(*this,& SelfType::eraseEdge);

        }

        void mergeEdges(const Edge & a,const Edge & b){
            const EdgeHolderType aa(mergeGraph_,a);
            const EdgeHolderType bb(mergeGraph_,b);
            object_.attr("mergeEdges")(aa,bb);
        }
        void mergeNodes(const Node & a,const Node & b){
            const NodeHolderType aa(mergeGraph_,a);
            const NodeHolderType bb(mergeGraph_,b);
            object_.attr("mergeNodes")(aa,bb);
        }
        void eraseEdge(const Edge & e){
            const EdgeHolderType ee(mergeGraph_,e);
            object_.attr("eraseEdge")(ee);
        }
        Edge contractionEdge(){
            EdgeHolderType eh = boost::python::extract<EdgeHolderType>(object_.attr("contractionEdge")());
            return eh;
        }
        WeightType contractionWeight()const{
            return boost::python::extract<WeightType>(object_.attr("contractionWeight")());
        }

        MergeGraph & mergeGraph(){
            return mergeGraph_;
        }
    private:
        MergeGraph & mergeGraph_;
        boost::python::object object_;
    };






    } // end namespace cluster_operators


} // end namespace vigra

#endif