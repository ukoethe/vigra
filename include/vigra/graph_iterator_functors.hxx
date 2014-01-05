    
#ifndef VIGRA_GRAPH_ITERATORS_HXX
#define VIGRA_GRAPH_ITERATORS_HXX

namespace vigra {


namespace detail{
    // iterator filters
    namespace filter{
        // functor to implement back neigh. iterator
        template<class ID_TYPE>
        struct SmallerThan{
            SmallerThan(){}
            SmallerThan(const ID_TYPE id): id_(id){}
            bool operator()(const ID_TYPE x) { return x < id_; }
            ID_TYPE id_;
        };
        // functor to implement back neigh.edge  iterator
        template<class GRAPH>
        struct BackEdgeIdFilter{
            typedef GRAPH GraphType;
            typedef typename GraphType::IdType IdType;
            BackEdgeIdFilter():graph_(NULL){}
            BackEdgeIdFilter(const GraphType & graph, const IdType ownNodeId)
            :   graph_(&graph),
                ownNodeId_(ownNodeId){
            }
            // return if the other nodes id is smaller than own id
            bool operator()(const IdType edgeId)const{
                const IdType otherNodeId =  graph_->id(
                    graph_->oppositeNode(   
                        graph_->nodeFromId(ownNodeId_) , graph_->edgeFromId(edgeId)   
                    ) 
                );
                return otherNodeId<ownNodeId_;
            }
            const GraphType * graph_;
            IdType ownNodeId_;
        };
    }

    namespace transform{

        template<class ID,class ITEM>
        struct IdToItem{
            ITEM operator()(const ID id)const{
                return ITEM(id);
            }
        };

        /*
        template<class GRAPH,class GRAPH_ITEM>
        struct IdToGraphItem;


        template<class GRAPH>
        struct IdToGraphItem<GRAPH,typename GRAPH::Arc>{
            IdToGraphItem():graph_(NULL){}
            IdToGraphItem(const GRAPH & graph):graph_(&graph){}
            typename GRAPH::Arc operator()(typename GRAPH::index_type id)const{
                return graph_->arcFromId(id);
            }
            const GRAPH * graph_;
        };

        template<class GRAPH>
        struct IdToGraphItem<GRAPH,typename GRAPH::Edge>{
            IdToGraphItem():graph_(NULL){}
            IdToGraphItem(const GRAPH & graph):graph_(&graph){}
            typename GRAPH::Edge operator()(typename GRAPH::index_type id)const{
                return graph_->edgeFromId(id);
            }
            const GRAPH * graph_;
        };


        template<class GRAPH>
        struct IdToGraphItem<GRAPH,typename GRAPH::Node>{
            IdToGraphItem():graph_(NULL){}
            IdToGraphItem(const GRAPH & graph):graph_(&graph){}
            typename GRAPH::Node  operator()(typename GRAPH::index_type id)const{
                return graph_->nodeFromId(id);
            }
            const GRAPH * graph_;
        };

        */

        // get other node id from edge id
        template<class GRAPH>
        struct OtherNodeId{
            typedef GRAPH GraphType;
            typedef typename GraphType::IdType IdType;
            OtherNodeId():graph_(NULL){}
            OtherNodeId(const GraphType & graph, const IdType ownNodeId)
            :   graph_(&graph),
                ownNodeId_(ownNodeId){
            }
            // return the other node's is
            IdType operator()(const IdType edgeId)const{
                //return graph_->edgeFromId(edgeId).otherNodeId(ownNodeId_);

                return graph_->id(
                    graph_->oppositeNode(   
                        graph_->nodeFromId(ownNodeId_) , graph_->edgeFromId(edgeId)   
                    ) 
                );
            }
            const GraphType * graph_;
            IdType ownNodeId_;
        };
        // get other node  from edge id
        template<class GRAPH>
        struct OtherNode{
            typedef GRAPH GraphType;
            typedef typename GraphType::IdType IdType;
            typedef typename GraphType::Node   Node;
            OtherNode():graph_(NULL){}
            OtherNode(const GraphType & graph, const IdType ownNodeId)
            :   graph_(&graph),
                ownNodeId_(ownNodeId){
            }
            // return the other node's is
            Node  operator()(IdType edgeId)const{
                //CGP_ASSERT_OP(graph_,!=,NULL);
                return graph_->oppositeNode(   graph_->nodeFromId(ownNodeId_) , graph_->edgeFromId(edgeId) ) ;
            }
            const GraphType * graph_;
            IdType ownNodeId_;
        };
    }


};


} // end namespace vigra



#endif //VIGRA_GRAPH_ITERATORS_HXX