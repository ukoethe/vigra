#ifndef VIGRA_GRAPH_CRTP_BASE_HXX
#define VIGRA_GRAPH_CRTP_BASE_HXX

/*boost*/
#include <boost/iterator/iterator_facade.hpp>

/*vigra*/
#include <vigra/graph_helper/graph_item_impl.hxx>
#include <vigra/is_end_filter_iterator.hxx>

namespace vigra{
namespace detail{



	template<class GRAPH,class ARC,class EDGE_IT>
	class FallbackArcIt
	: public boost::iterator_facade<
        FallbackArcIt<GRAPH,ARC,EDGE_IT>,
        const ARC,
        boost::forward_traversal_tag
    >

	{
	public:
		FallbackArcIt()
		:	graph_(NULL),
			pos_(),
			inFirstHalf_(false),
			veryEnd_(true){
		}

		FallbackArcIt(const GRAPH & g,EDGE_IT pos,const bool inFirstHalf,const bool veryEnd)
		:	graph_(&g),
			pos_(pos),
			inFirstHalf_(inFirstHalf),
			veryEnd_(veryEnd){
		}

		bool isEnd()const{
			return veryEnd_ || graph_==NULL;
		}

	private:
		bool isBegin()const{
			return veryEnd_==false && pos_ == graph_->edgesBegin();			
		}


		friend class boost::iterator_core_access;

        void increment() {




        	
        	//std::cout<<"increment \n";
        	if(inFirstHalf_){
        		//std::cout<<"first half very end"<<veryEnd_<<"\n";
        		++pos_;
        		if(pos_ == lemon::INVALID  ) {
        			//std::cout<<"invalid in first half\n";
        			pos_ =graph_->edgesBegin();
        			inFirstHalf_=false;
        		}
        		return;
        	}
        	else{
        		//std::cout<<"last half \n";
        		++pos_;
        		if(pos_ == lemon::INVALID){
        			//std::cout<<"invalid in secomd half\n";
        			veryEnd_=true;
        		}
        		return;
        	}
        
           
        }
        bool equal(FallbackArcIt const& other) const{
            return (
            	(
            		isEnd()==other.isEnd() 					&&
            		inFirstHalf_==other.inFirstHalf_ 
            	) &&
            	(isEnd() || graph_==NULL || pos_==other.pos_ )
            	);
            	
        }

        const ARC & dereference() const { 
        	std::cout<<graph_->id(*pos_)<<"\n";
        	if(inFirstHalf_)
            	arc_ = ARC(graph_->id(*pos_));
            else
            	arc_ = ARC(graph_->id(*pos_)+graph_->maxEdgeId()+1);
           	return arc_;
        }


		const GRAPH * graph_;
		EDGE_IT pos_;
		bool inFirstHalf_;
		bool veryEnd_;
		mutable ARC arc_;
	};


	template<class GRAPH>
	class OutEdgeOnlyFilter{
	public:
		typedef typename GRAPH::Edge 		Edge;
		typedef typename GRAPH::index_type 	index_type;
		OutEdgeOnlyFilter()
		: 	graph_(NULL),
			ownNodeId_(-1){
		}
		OutEdgeOnlyFilter(const GRAPH & graph, const index_type ownNodeId)
		: 	graph_(&graph),
			ownNodeId_(ownNodeId){
		}
		bool operator()(const Edge & edge)const{
			const index_type uid = graph_->id(graph_->u(edge));
			return uid==ownNodeId_;
		}
	private:
		const GRAPH * graph_;
		index_type ownNodeId_;
	};

	template<class GRAPH>
	class InEdgeOnlyFilter{
	public:
		typedef typename GRAPH::Edge 		Edge;
		typedef typename GRAPH::index_type 	index_type;
		InEdgeOnlyFilter()
		: 	graph_(NULL),
			ownNodeId_(-1){
		}
		InEdgeOnlyFilter(const GRAPH & graph, const index_type ownNodeId)
		: 	graph_(&graph),
			ownNodeId_(ownNodeId){
		}
		bool operator()(const Edge & edge)const{
			const index_type vid = graph_->id(graph_->v(edge));
			return vid==ownNodeId_;
		}
	private:
		const GRAPH * graph_;
		index_type ownNodeId_;
	};

	template<class GRAPH>
	class BackEdgeOnlyFilter{
	public:
		typedef typename GRAPH::Edge 		Edge;
		typedef typename GRAPH::index_type 	index_type;
		BackEdgeOnlyFilter()
		: 	graph_(NULL),
			ownNodeId_(-1){
		}
		BackEdgeOnlyFilter(const GRAPH & graph, const index_type ownNodeId)
		: 	graph_(&graph),
			ownNodeId_(ownNodeId){
		}
		bool operator()(const Edge & edge)const{
			const index_type otherNodeId = graph_->id(graph_->oppositeNode(graph_->nodeFromId(ownNodeId_),edge));
			return otherNodeId<ownNodeId_;
		}
	private:
		const GRAPH * graph_;
		index_type ownNodeId_;
	};




	template<	class GRAPH,	class INDEX_TYPE,
				class EDGE,		class NODE,
				class EDGE_IT,	class NODE_IT
	>
	class ArcHelper{

	public:


		typedef GenericArc<INDEX_TYPE> 			  Arc;
		typedef FallbackArcIt<GRAPH,Arc,EDGE_IT>  ArcIt;


		size_t arcNum()const{
			return graph().edgeNum()*2;
		}
		INDEX_TYPE maxArcId()const{
			return  graph().maxEdgeId()*2+1;
		}
		INDEX_TYPE arcId(const Arc & arc)const{
			return arc.id();
		}

		ArcIt arcsBegin()const{
			return ArcIt(graph(),graph().edgesBegin(),true,false);
		}
		ArcIt arcsEnd()const{
			return ArcIt(graph(),graph().edgesBegin(),false,true);
		}


        NODE oppositeNode(NODE const &n, const EDGE &e) const {
            const NODE uNode = graph().u(e);
            const NODE vNode = graph().v(e);
            if(graph().id(uNode)==graph().id(n)){
                return vNode;
            }
            else if(graph().id(vNode)==graph().id(n)){
                return uNode;
            }
            else{
                return NODE(-1);
            }
        }

		NODE source(const Arc & arc)const{
			const INDEX_TYPE arcIndex  = graph().id(arc);
			const INDEX_TYPE maxEdgeId = graph().maxEdgeId();

			if (arcIndex > maxEdgeId ){
				const INDEX_TYPE edgeIndex = arcIndex-maxEdgeId-1;
				const EDGE edge = graph().edgeFromId(edgeIndex);
				return graph().v(edge);
			}
			else{
				const INDEX_TYPE edgeIndex = arcIndex;
				const EDGE edge = graph().edgeFromId(edgeIndex);
				return graph().u(edge);
			}
		}

		NODE target(const Arc & arc)const{
			const INDEX_TYPE arcIndex  = graph().id(arc);
			const INDEX_TYPE maxEdgeId = graph().maxEdgeId();

			if (arcIndex > maxEdgeId ){
				const INDEX_TYPE edgeIndex = arcIndex-maxEdgeId-1;
				const EDGE edge = graph().edgeFromId(edgeIndex);
				return graph().u(edge);
			}
			else{
				const INDEX_TYPE edgeIndex = arcIndex;
				const EDGE edge = graph().edgeFromId(edgeIndex);
				return graph().v(edge);
			}
		}

	private:
		const GRAPH & graph()const{
			return * static_cast<const GRAPH*>(this);
		}
		GRAPH & graph(){
			return * static_cast<GRAPH*>(this);
		}

	};




}	
}


#endif // VIGRA_GRAPH_CRTP_BASE_HXX