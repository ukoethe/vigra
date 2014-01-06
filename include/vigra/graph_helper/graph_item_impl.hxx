#ifndef VIGRA_NODE_IMPL_HXX
#define VIGRA_NODE_IMPL_HXX

#include <vigra/is_end_iterator.hxx>
#include <vigra/algorithm.hxx>
#include <vigra/tinyvector.hxx>

namespace vigra{
namespace detail{



template<class INDEX_TYPE,class SET_TYPE>
class GenericNodeImpl{

    public:
        typedef INDEX_TYPE index_type;
        typedef SET_TYPE EdgeIdSet;
        typedef IsEndIter< typename EdgeIdSet::const_iterator > EdgeIdIt;



    private:
        //GenericNodeImpl();                               // non empty-construction
        //GenericNodeImpl( const GenericNodeImpl& other );      // non construction-copyable
        //GenericNodeImpl & operator=( const GenericNodeImpl& ); // non assignable
    public:
        GenericNodeImpl(){

        }
        //GenericNodeImpl(const index_type id)
        //:   id_(id){
        // }
        // query
        size_t numberOfEdges()const{return edges_.size();}
        size_t edgeNum()const{return edges_.size();}
         size_t num_edges()const{return edges_.size();}

        bool hasEdgeId(const index_type edge)const{return edges_.find(edge)!=edges_.end();}

        // modification
        void  mergeEdges(const GenericNodeImpl & other){
            edges_.insert(other.edges_.begin(),other.edges_.end());
        }

        
        std::pair<index_type,bool> sharedEdge(const GenericNodeImpl & other)const{
            return detail::first_set_intersection<index_type>(
                edgeIdsBegin(),edgeIdsEnd(),
                other.edgeIdsBegin(),other.edgeIdsEnd()
            );
        }
        

        bool eraseEdge(const size_t edgeIndex){
            return edges_.erase(edgeIndex)==1;
        }

        void eraseAndInsert(const index_type removeEdge,const index_type insertEdge){
            edges_.erase(removeEdge);
            edges_.insert(insertEdge);
        }

        EdgeIdIt edgeIdsBegin()const{
            return EdgeIdIt(edges_.begin() ,edges_.end()) ;
        }
        EdgeIdIt edgeIdsEnd()const{
            return EdgeIdIt(edges_.end() ,edges_.end()) ;
        }

        index_type id()const{
            return id_;
        }
        void clear(){
            edges_.clear();
        }

        const EdgeIdSet & edgeIdSet()const{
            return edges_;
        }

        void insertEdgeId(const index_type id){
        	edges_.insert(id);
        }
    public:

        EdgeIdSet edges_;
        index_type id_;
};

template<class INDEX_TYPE>
class GenericEdgeImpl
:  public vigra::TinyVector<INDEX_TYPE,3> {
        // public typedefs
    public:
        typedef INDEX_TYPE index_type;

        GenericEdgeImpl(){}
    // public methods
    public:
        index_type u()const{return this->operator[](0);}
        index_type v()const{return this->operator[](1);}
        index_type id()const{return this->operator[](2);}

    // private 
    private:
};


template<class INDEX_TYPE>
class GenericEdge;

template<class INDEX_TYPE>
class GenericArc{
public:
    typedef INDEX_TYPE index_type;
    GenericArc(
        const index_type id = static_cast<index_type>(-1),
        const index_type edgeId = static_cast<index_type>(-1) 
    )
    :   id_(id),
        edgeId_(edgeId){

    }
    index_type id()const{return id_;}
    index_type edgeId()const{return id_;}

    operator GenericEdge<INDEX_TYPE> () const{
        return GenericEdge<INDEX_TYPE>(edgeId());
    }

private:
    index_type id_;
    index_type edgeId_;
};

template<class INDEX_TYPE>
class GenericEdge{
public:
    typedef INDEX_TYPE index_type;


    GenericEdge(const index_type id = static_cast<index_type>(-1) )
    : id_(id){

    }

    GenericEdge(const GenericArc<INDEX_TYPE> & arc)
    :   id_(arc.edgeId())
    {
    }




    index_type id()const{return id_;}
private:
    index_type id_;
};


template<class INDEX_TYPE>
class GenericNode{
public:
    typedef INDEX_TYPE index_type;
    GenericNode(const index_type id = static_cast<index_type>(-1) )
    : id_(id){
        
    }
    index_type id()const{return id_;}
private:
    index_type id_;
};

template<class INDEX_TYPE>
inline bool operator == (const GenericEdge<INDEX_TYPE> & edge,const lemon::Invalid & iv){
    return edge.id()==-1; 
}
template<class INDEX_TYPE>
inline bool operator == (const lemon::Invalid & iv , const GenericEdge<INDEX_TYPE> & edge){
    return edge.id()==-1;
}

template<class INDEX_TYPE>
inline bool operator != (const GenericEdge<INDEX_TYPE> & edge,const lemon::Invalid & iv){
    return edge.id()!=-1;
}
template<class INDEX_TYPE>
inline bool operator != (const lemon::Invalid & iv , const GenericEdge<INDEX_TYPE> & edge){
    return edge.id()!=-1;
}

template<class INDEX_TYPE>
inline bool operator == (const GenericNode<INDEX_TYPE> & node,const lemon::Invalid & iv){
    return node.id()==-1; 
}
template<class INDEX_TYPE>
inline bool operator == (const lemon::Invalid & iv , const GenericNode<INDEX_TYPE> & node){
    return node.id()==-1;
}

template<class INDEX_TYPE>
inline bool operator != (const GenericNode<INDEX_TYPE> & node,const lemon::Invalid & iv){
    return node.id()!=-1;
}
template<class INDEX_TYPE>
inline bool operator != (const lemon::Invalid & iv , const GenericNode<INDEX_TYPE> & node){
    return node.id()!=-1;
}


}
} // end namespace vigra



#endif // VIGRA_NODE_IMPL_HXX
