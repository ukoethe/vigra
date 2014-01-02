    
#ifndef VIGRA_MERGE_GRAPH_ITERATORS_HXX
#define VIGRA_MERGE_GRAPH_ITERATORS_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* boost */
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>


/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graphs.hxx>

/* this project*/
#include "merge_graph_edge.hxx"
#include "merge_graph_node.hxx"
#include "merge_graph_item.hxx"
#include "merge_graph_callbacks.hxx"

/* this project (TO BE REFACTORED) */
#include "iterable_partition.hxx"
#include "macros.hxx"

namespace vigra {


// iterator functors
namespace merge_graph_detail{

    namespace filter{

        template<class ID_TYPE>
        struct SmallerThan{
            SmallerThan(const ID_TYPE id)
            : id_(id){
            }
            bool operator()(const ID_TYPE x) { 
                return x < id_; 
            }
            const ID_TYPE id_;
        };


        template<class MERGE_GRAPH>
        struct IsOutgoingArc{

            typedef MERGE_GRAPH MergeGraphType;
            typedef typename MergeGraphType::Edge Edge;
            typedef typename MergeGraphType::IdType IdType;

            IsOutgoingArc(const MergeGraphType & graph, const IdType nodeId)
            :   graph_(&graph),
                nodeId_(nodeId){
            }
            bool operator()(const Edge edge) { 
                return edge[0]==nodeId_;
            }
            const MergeGraphType * graph_;
            const IdType nodeId_;
        };
    }

    namespace transform{
        template<class MERGE_GRAPH,class GRAPH_ITEM>
        struct IdToGraphItem;


        template<class MERGE_GRAPH>
        struct IdToGraphItem<MERGE_GRAPH,typename MERGE_GRAPH::Edge>{
            IdToGraphItem()
            :graph_(NULL){
            }
            IdToGraphItem(const MERGE_GRAPH & graph)
            :  graph_(&graph){
            }
            typename MERGE_GRAPH::Edge operator()(typename MERGE_GRAPH::IdType id)const{
                return graph_->edgeFromId(id);
            }
            const MERGE_GRAPH * graph_;
        };

        template<class MERGE_GRAPH>
        struct IdToGraphItem<MERGE_GRAPH,typename MERGE_GRAPH::Node>{
            IdToGraphItem()
            :graph_(NULL){
            }
            IdToGraphItem(const MERGE_GRAPH & graph)
            :  graph_(&graph){
            }
            const typename MERGE_GRAPH::Node & operator()(typename MERGE_GRAPH::IdType id)const{
                return graph_->nodeFromId(id);
            }
            const MERGE_GRAPH * graph_;
        };

        // get other node id from edge id
        template<class MERGE_GRAPH>
        struct OtherNodeId{
            typedef MERGE_GRAPH MergeGraphType;
            typedef typename MergeGraphType::IdType IdType;

            OtherNodeId(const MergeGraphType & graph, const IdType ownNodeId)
            :   graph_(&graph),
                ownNodeId_(ownNodeId){
            }

            // return the other node's is
            IdType operator()(const IdType edgeId)const{
                return graph_->edgeFromId().otherNodeId(ownNodeId_);
            }
            const MergeGraphType graph_;
            IdType ownNodeId_;

        };
    }




    template<class ITER>
    class IsEndIter
    :   public ITER {
    public:
        IsEndIter()
        :   valid_(false){            
        }

        IsEndIter(ITER iter , ITER end)
        :   ITER(iter),
            end_(end),
            valid_(true){
        }

        bool isEnd()const{
            return !valid_ || *this==end_;
        }
    private:
        ITER end_;
        bool valid_;
    };

};


template<class MERGE_GRAPH,class GRAPH_ITEM_TYPE>
class  MergeGraphItemIterator
:   public boost::transform_iterator<
        merge_graph_detail::transform::IdToGraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>, 
        typename merge_graph_detail::GraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>::GraphItemIdIt,
        typename merge_graph_detail::GraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>::result
    > 
{
public: 
    MergeGraphItemIterator(const MERGE_GRAPH & graph,typename merge_graph_detail::GraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>::GraphItemIdIt iterator)
    : boost::transform_iterator<
        merge_graph_detail::transform::IdToGraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>, 
        typename merge_graph_detail::GraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>::GraphItemIdIt,
        typename merge_graph_detail::GraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>::result
    >(iterator,merge_graph_detail::transform::IdToGraphItem<MERGE_GRAPH,GRAPH_ITEM_TYPE>(graph)){
    }
    bool isBegin()const{
        return this->base().isBegin();
    }
    bool isEnd()const{
        return this->base().isEnd();
    }
};



template<class MERGE_GRAPH>
class  MergeGraphNeigbourhoodIdIterator
:   public boost::transform_iterator<
        merge_graph_detail::transform::OtherNodeId<MERGE_GRAPH>, 
        typename MERGE_GRAPH::Node::EdgeIdSet::const_iterator,const typename MERGE_GRAPH::Node & ,typename MERGE_GRAPH::Node
    > 
{
    typedef MERGE_GRAPH MergeGraphType;
    typedef typename MergeGraphType::IdType IdType;
    typedef typename MergeGraphType::Node::EdgeIdSet::const_iterator BaseIteratorType;
public:

    MergeGraphNeigbourhoodIdIterator(const MergeGraphType & graph,const IdType ownNodeId ,const  BaseIteratorType pos)
    :   boost::transform_iterator<
            merge_graph_detail::transform::OtherNodeId<MERGE_GRAPH>, 
            typename MERGE_GRAPH::Node::EdgeIdSet::const_iterator,const typename MERGE_GRAPH::Node & ,typename MERGE_GRAPH::Node
        > (pos,merge_graph_detail::transform::OtherNodeId<MERGE_GRAPH>(graph,ownNodeId))
    {
    }
    bool isEnd()const{
        //return this->base()==end_;
        return this->base().isEnd();
    }
private:
    //BaseIteratorType end_;
};





template<class MERGE_GRAPH>
class  MergeGraphNeigbourhoodIterator
:  public boost::iterator_facade
< 
    MergeGraphNeigbourhoodIterator<MERGE_GRAPH>,
    typename MERGE_GRAPH::Node const,
    boost::bidirectional_traversal_tag
>

{

private:
    typedef MERGE_GRAPH MergeGraphType;
    typedef typename MergeGraphType::IdType IdType;
    typedef typename MergeGraphType::Node Node;
    typedef typename MergeGraphType::Edge Edge;
    typedef typename Node::EdgeIdSet EdgeIdSet;
    typedef typename EdgeIdSet::const_iterator EdgeIdIter;
public: 
    MergeGraphNeigbourhoodIterator(){
        mergeGraph_=NULL;
        node_=NULL;
    }

    MergeGraphNeigbourhoodIterator(
        const MergeGraphType & mergeGraph , 
        const Node & node ,
        EdgeIdIter edgeIdIter 
    )
    :   mergeGraph_(&mergeGraph),
        node_(&node),
        edgeIdIter_(edgeIdIter)
    {

    }


    bool isBegin()const{
        return mergeGraph_!=NULL && edgeIdIter_== node_->edgeIdsBegin();
    }

    bool isEnd()const{
        return mergeGraph_==NULL || edgeIdIter_== node_->edgeIdsEnd();
    }

    bool equal(const MergeGraphNeigbourhoodIterator & other)const{
        return   (this->isEnd() && other.isEnd() )  || ((this->isEnd()==other.isEnd() ) && edgeIdIter_==other.edgeIdIter_);
    }

    void increment(){
        ++edgeIdIter_;
    }
    void decrement(){
        --edgeIdIter_;
    }

    const Node & dereference()const{
        const Edge edge =  mergeGraph_->edgeFromId(*edgeIdIter_);
        return mergeGraph_->nodeFromId(edge.otherNodeId(node_->id()));
    }


private:
    const MergeGraphType * mergeGraph_;
    Node                 * node_;
    EdgeIdIter             edgeIdIter_;
    IdType nodeId_;
    

};




namespace merge_graph_detail{

};  





}



#endif //VIGRA_MERGE_GRAPH_ITERATORS_HXX