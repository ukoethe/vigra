    
#ifndef VIGRA_MERGE_GRAPH_NODE_HXX
#define VIGRA_MERGE_GRAPH_NODE_HXX

/* std library */
#include <set>


namespace vigra {

template<class LABEL_TYPE>
class MergeGraphNode{

    public:
        typedef LABEL_TYPE LabelType;


        typedef typename std::set<LabelType>::const_iterator EdgeIterator;

        // query
        size_t numberOfEdges()const{return edges_.size();}
        bool hasEdge(const LabelType edge)const{return edges_.find(edge)!=edges_.end();}

        // modification
        void  mergeEdges(const MergeGraphNode & other){
            edges_.insert(other.edges_.begin(),other.edges_.end());
        }

        bool eraseEdge(const size_t edgeIndex){
            return edges_.erase(edgeIndex)==1;
        }

        void eraseAndInsert(const LabelType removeEdge,const LabelType insertEdge){
            edges_.erase(removeEdge);
            edges_.insert(insertEdge);
        }

        EdgeIterator edgesBegin()const{
            return edges_.begin();
        }
        EdgeIterator edgesEnd()const{
            return edges_.end();
        }
    std::set<LabelType> edges_;
};

}



#endif //VIGRA_MERGE_GRAPH_NODE_HXX