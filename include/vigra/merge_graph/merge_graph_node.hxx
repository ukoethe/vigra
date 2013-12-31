    
#ifndef VIGRA_MERGE_GRAPH_NODE_HXX
#define VIGRA_MERGE_GRAPH_NODE_HXX

/* std library */
#include <set>


namespace vigra {


namespace merge_graph_detail{

template <class OUT_TYPE,class InputIterator1, class InputIterator2>
std::pair<OUT_TYPE,bool> first_set_intersection (
    InputIterator1 first1, 
    InputIterator1 last1,
    InputIterator2 first2, 
    InputIterator2 last2
)
{
    while (first1!=last1 && first2!=last2){
        if (*first1<*first2) 
            ++first1;
        else if (*first2<*first1) 
            ++first2;
        else {
            return std::pair<OUT_TYPE,bool>(static_cast<OUT_TYPE>(*first1),true);
        }
    }
    return std::pair<OUT_TYPE,bool>(OUT_TYPE( ),false); 
}

} // end namespace merge_graph_detail



template<class ID_TYPE>
class MergeGraphNode{

    public:
        typedef ID_TYPE IdType;


        typedef typename std::set<IdType>::const_iterator EdgeIdIt;


        MergeGraphNode(){}
        MergeGraphNode(const IdType id)
        :   id_(id){
        }
        // query
        size_t numberOfEdges()const{return edges_.size();}
        bool hasEdgeId(const IdType edge)const{return edges_.find(edge)!=edges_.end();}

        // modification
        void  mergeEdges(const MergeGraphNode & other){
            edges_.insert(other.edges_.begin(),other.edges_.end());
        }

        
        std::pair<IdType,bool> sharedEdge(const MergeGraphNode & other)const{
            return merge_graph_detail::first_set_intersection<IdType>(
                edgeIdsBegin(),edgeIdsEnd(),
                other.edgeIdsBegin(),other.edgeIdsEnd()
            );
        }
        

        bool eraseEdge(const size_t edgeIndex){
            return edges_.erase(edgeIndex)==1;
        }

        void eraseAndInsert(const IdType removeEdge,const IdType insertEdge){
            edges_.erase(removeEdge);
            edges_.insert(insertEdge);
        }

        EdgeIdIt edgeIdsBegin()const{
            return edges_.begin();
        }
        EdgeIdIt edgeIdsEnd()const{
            return edges_.end();
        }

        IdType id()const{
            return id_;
        }
        void clear(){
            edges_.clear();
        }

    std::set<IdType> edges_;
    IdType id_;
};

}



#endif //VIGRA_MERGE_GRAPH_NODE_HXX