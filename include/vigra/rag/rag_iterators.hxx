#ifndef VIGRA_RAG_ITERATORS_HXX
#define VIGRA_RAG_ITERATORS_HXX

/*boost*/
#include <boost/iterator/iterator_facade.hpp>

/*vigra*/
#include <vigra/is_end_enumeration_iterator.hxx>	

namespace vigra{
	namespace detail{
		    template<class GRAPH,class ITEM>
    struct RagItemIt 
    :   public EnumerationIterator< ITEM > {
    private:
        typedef EnumerationIterator< ITEM > BaseIterType;
    public:
        typedef GRAPH Graph;
        typedef typename Graph::index_type index_type;
        // make this private
        RagItemIt(const index_type pos,const index_type size)
        :   BaseIterType(pos,size){
        }
        // default constructor
        RagItemIt()
        :   BaseIterType(0,0,1){
        }
        RagItemIt(const RagItemIt & other)
        :   BaseIterType(other){
        }

        // Invalid constructor & conversion. 
        RagItemIt(const lemon::Invalid & invalid)
        :   BaseIterType(0,0,1){
        }
        RagItemIt(const Graph & g)
        :   BaseIterType(0, GraphItemHelper<GRAPH,ITEM>::itemNum(g) ,1){
        }
        RagItemIt(const Graph & g,const ITEM & item)
        :   BaseIterType(g.id(item)-1, GraphItemHelper<GRAPH,ITEM>::itemNum(g) , 1 ){
        }
    };

    // specialization for arc
    template<class GRAPH>
    class RagItemIt<GRAPH,typename GRAPH::Arc>
    : public boost::iterator_facade<
        RagItemIt<GRAPH,typename GRAPH::Arc>,
        const typename GRAPH::Arc,
        //boost::forward_traversal_tag 
        boost::forward_traversal_tag
    >
    {
    public:
        typedef GRAPH Graph;
        typedef typename Graph::index_type index_type;
        typedef typename Graph::Arc   Arc;

        RagItemIt(const lemon::Invalid &  invalid = lemon::INVALID)
        :   current_(0),
            maxEdgeId_(0),
            maxArcId_(0){
        }

        RagItemIt(const Graph & g)
        :   current_(1),
            maxEdgeId_(g.maxEdgeId()),
            maxArcId_(g.maxArcId()){
        }
        RagItemIt(const Graph & g,const Arc & arc)
        :   current_(g.id(arc)),
            maxEdgeId_(g.maxEdgeId()),
            maxArcId_(g.maxArcId()){
        }
    private:
    	friend class boost::iterator_core_access;
        bool isEnd()const{
            return  (maxArcId_==0 && current_==0 )||  current_>maxArcId_;
        }

        bool isBegin()const{
            return current_ == 0 && maxArcId_!=0;
        }



        

        void increment() {
            ++current_;
        }
        

        
        bool equal(const RagItemIt & other) const{
            return   (isEnd() && other.isEnd() ) || (current_ == other.current_);
        }

        const Arc & dereference() const { 
            if(current_<=maxEdgeId_){
                arc_ =  Arc(current_,current_);
            }
            else {
                arc_ =  Arc(current_,current_-maxEdgeId_);
            }
            return arc_;
        }




        ptrdiff_t current_;
        ptrdiff_t maxEdgeId_;
        ptrdiff_t maxArcId_;
        mutable Arc arc_;
    };


	}
}

	
#endif