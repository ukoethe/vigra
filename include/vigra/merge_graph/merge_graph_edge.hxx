    
#ifndef VIGRA_MERGE_GRAPH_EDGE_HXX
#define VIGRA_MERGE_GRAPH_EDGE_HXX


/* this project (TO BE REFACTORED) */
#include "macros.hxx"

namespace vigra {


template<class ID_TYPE>
class MergeGraphEdge{

    public:
        typedef ID_TYPE IdType;
        MergeGraphEdge(){}
        MergeGraphEdge(const IdType a,const IdType b)
        :   first(a),
            second(b){
        }
        bool hasNode(const IdType node)const{
            return node==first || node==second;
        }
        IdType otherNodeId(const IdType node)const{
            CGP_ASSERT_OP(hasNode(node),==,true);
            return (node==first ? second : first);
        }
        const IdType & operator[](const IdType i)const{
            return (i==0 ? first : second);
        }
        IdType & operator[](const IdType i){
            return (i==0 ? first : second);
        }
    //private:
        IdType first;
        IdType second;
};



template<class ID_TYPE>
class MergeGraphArc{

    public:
        typedef ID_TYPE IdType;
        MergeGraphArc(){}
        MergeGraphArc(const IdType a,const IdType b)
        :   first(a),
            second(b){
        }
        bool hasNode(const IdType node)const{
            return node==first || node==second;
        }
        IdType otherNodeId(const IdType node)const{
            CGP_ASSERT_OP(hasNode(node),==,true);
            return (node==first ? second : first);
        }
        const IdType & operator[](const IdType i)const{
            return (i==0 ? first : second);
        }
        IdType & operator[](const IdType i){
            return (i==0 ? first : second);
        }
    //private:
        IdType first;
        IdType second;
};



} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_EDGE_HXX