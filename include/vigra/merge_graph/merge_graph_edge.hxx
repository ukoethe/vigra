    
#ifndef VIGRA_MERGE_GRAPH_EDGE_HXX
#define VIGRA_MERGE_GRAPH_EDGE_HXX


/* this project (TO BE REFACTORED) */
#include "macros.hxx"

namespace vigra {


template<class LABEL_TYPE>
class MergeGraphEdge{

    public:
        typedef LABEL_TYPE LabelType;
        bool hasNode(const LabelType node)const{
            return node==first || node==second;
        }
        LabelType otherNode(const LabelType node)const{
            CGP_ASSERT_OP(hasNode(node),==,true);
            return (node==first ? second : first);
        }
        const LabelType & operator[](const LabelType i)const{
            return (i==0 ? first : second);
        }
        LabelType & operator[](const LabelType i){
            return (i==0 ? first : second);
        }
    //private:
        LabelType first;
        LabelType second;
};


} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_EDGE_HXX