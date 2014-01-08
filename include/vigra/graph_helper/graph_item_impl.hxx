#ifndef VIGRA_NODE_IMPL_HXX
#define VIGRA_NODE_IMPL_HXX

#include <vigra/is_end_iterator.hxx>
#include <vigra/algorithm.hxx>
#include <vigra/tinyvector.hxx>
#include <boost/iterator/iterator_facade.hpp>

namespace vigra{
namespace detail{

/*
class MyIterator : public std::iterator<std::input_iterator_tag, int>
{
  int* p;
public:
  MyIterator(int* x) :p(x) {}
  MyIterator(const MyIterator& mit) : p(mit.p) {}
  MyIterator& operator++() {++p;return *this;}
  MyIterator operator++(int) {MyIterator tmp(*this); operator++(); return tmp;}
  bool operator==(const MyIterator& rhs) {return p==rhs.p;}
  bool operator!=(const MyIterator& rhs) {return p!=rhs.p;}
  int& operator*() {return *p;}
};
*/


template<class GRAPH>
struct NeighborNodeFilter{
    typedef typename GRAPH::Node ResultType;


    bool static valid(
    ResultType const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        return true;
    }


     static transform(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        return g.oppositeNode(g.nodeFromId(ownNodeId), g.edgeFromId(edegeId) );    
    }

    static bool IsFilter = false ; 
};

template<class GRAPH>
struct IncEdgeFilter{
    typedef typename GRAPH::Edge ResultType;

    template<GRAPH>
    bool static valid(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        return true;
    }

    template<GRAPH>
    ResultType static transform(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        return g.edgeFromId(edgeId);    
    }

    static bool IsFilter = false ; 
};

template<class GRAPH>
struct BackEdgeFilter{
    typedef typename GRAPH::Edge ResultType;

    bool static valid(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId
    ){
        return g.id(g.oppositeNode(g.nodeFromId(ownNodeId), g.edgeFromId(edegeId) )) < ownNodeId;
    } 

    ResultType static transform(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        return g.edgeFromId(edgeId);
    }

    static bool IsFilter = true ; 
};

template<class GRAPH>
struct IsOutFilter{
    typedef typename GRAPH::Arc ResultType;

    bool static valid(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId
    ){
        return  g.id(g.u(g.edgeId(edegeId)))==ownNodeId;
    } 
    ResultType static transform(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        const typename GRAPH::Edge edge g.edgeFromId(edgeId);
        return g.direct(edge);
        //direct (Edge const &e, Node const &n)
    }

    static bool IsFilter = true ; 
};

template<class GRAPH>
struct IsInFilter{
    typedef typename GRAPH::Arc ResultType;

    bool static valid(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId
    ){
        return  g.id(g.v(g.edgeId(edegeId)))==ownNodeId;
    } 
    ResultType static transform(
        const GRAPH & g,
        const typename GRAPH::index_type edegeId,
        const typename GRAPH::index_type ownNodeId, 
        const const typename GRAPH::index_type otherItemId
    ){
        const typename GRAPH::Edge edge g.edgeFromId(edgeId);
        return g.direct(edge,g.v(edge));
    }
    static bool IsFilter = true ; 
};


//  This class with the filter above can create 
//  the following lemon iterators :
//   - IncEdgeIt
//   - OutArcIt
//   - 
//
//
template<class GRAPH,class NODE_IMPL,class FILTER>    
class GenericIncEdgeIt

:  public boost::iterator_facade<
      GenericIncEdgeIt<GRAPH,NODE_IMPL,FILTER>,
      typename FILTER::ResultType const,
      boost::forward_traversal_tag
   >
{
public:

    typedef GRAPH Graph;
    typedef typename Graph::NodeIt NodeIt;
    typedef Graph::Edge Edge;
    typedef Graph::Node Node;
    typedef typename GraphItemHelper<typename FILTER::ResultType>  ResultItem

    // default constructor
    GenericIncEdgeIt(const lemon::Invalid & invalid = lemon::INVALID)
    :   nodeImpl_(NULL),
        graph_(NULL),
        ownNodeId_(-1),
        idIter_(),
        resultItem_(lemon::INVALID){
    }

    // from a given node
    GenericIncEdgeIt(const Graph & g , const Node & node)
    :   nodeImpl_(NULL),
        graph_(&g),
        idIter_(),
        ownNodeId_(g.id(node)),
        resultItem_(lemon::INVALID){
            // get first node
            const Node firstNode  * NodeIt(g);
            nodeImpl_(&g.nodeImpl(firstNode));
    }

private:
    friend class boost::iterator_core_access;

    typedef NODE_IMPL NodeImpl;
    typedef typename NodeImpl::EdgeIdIt EdgeIdIt;

    bool isEnd()const{
        return  (nodeImpl_==NULL  || idIter_==nodeImpl_.edgeIdsEnd());      
    }
    bool isBegin()const{
        return (nodeImpl_!=NULL &&  idIter_==nodeImpl_.edgeIdsBegin());
    }

    void increment()const{
        ++idIter_;
        if(FILTER::IsFilter){
            while(idIter_!=nodeImpl_->edgeIdsEnd() && !FILTER::valid(*graph_,*idIter_,ownNodeId_){
                ++idIter_;
            }
        }
    }

    // might no need to make this constant
    // therefore we would lose the "mutabe"
    const ResultItem & dereference()const{
        resultItem_ =  FILTER::transform(*graph_,*idIter_,ownNodeId_)
        return resultItem_;
    }


    const NODE_IMPL * nodeImpl_;
    const GRAPH     * graph_;
    const index_type  ownNodeId_;
    EdgeIdIt idIter_;
    mutable ResultItem resultItem_;
};



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

    GenericArc(const lemon::Invalid & iv)
    :   id_(-1),
        edgeId_(-1){

    }



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

    GenericEdge(const lemon::Invalid & iv)
    : id_(-1){

    }

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

    GenericNode(const lemon::Invalid & iv)
    : id_(-1){

    }

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
