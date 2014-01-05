#ifndef VIGRA_GRAPH_CRTP_BASE_HXX
#define VIGRA_GRAPH_CRTP_BASE_HXX

#include <vigra/graph_helper/graph_item_impl.hxx>
#include <boost/iterator/iterator_facade.hpp>

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





	template<class GRAPH,class INDEX,class EDGE,class NODE,class EDGE_IT,class NODE_IT>
	class ArcHelper{

	public:
		typedef INDEX index_type;
		typedef EDGE  Edge;
		typedef NODE  Node;
		typedef GenericArc<index_type> Arc;
		typedef FallbackArcIt<GRAPH,Arc,EDGE_IT>  ArcIt;


		size_t arcNum()const{
			return graph().edgeNum()*2;
		}
		index_type maxArcId()const{
			return  graph().arcNum()-1;
		}
		index_type arcId(const Arc & arc)const{
			return arc.id();
		}

		ArcIt arcsBegin()const{
			return ArcIt(graph(),graph().edgesBegin(),true,false);
		}
		ArcIt arcsEnd()const{
			return ArcIt(graph(),graph().edgesBegin(),false,true);
		}

		Node source(const Arc & arc)const{
			const index_type arcIndex  = graph().id(arc);
			const index_type maxEdgeId = graph().maxEdgeId();

			if (arcIndex > maxEdgeId ){
				const index_type edgeIndex = arcIndex-maxEdgeId-1;
				const Edge edge = graph().edgeFromId(edgeIndex);
				return graph().v(edge);
			}
			else{
				const index_type edgeIndex = arcIndex;
				const Edge edge = graph().edgeFromId(edgeIndex);
				return graph().u(edge);
			}
		}

		Node target(const Arc & arc)const{
			const index_type arcIndex  = graph().id(arc);
			const index_type maxEdgeId = graph().maxEdgeId();

			if (arcIndex > maxEdgeId ){
				const index_type edgeIndex = arcIndex-maxEdgeId-1;
				const Edge edge = graph().edgeFromId(edgeIndex);
				return graph().u(edge);
			}
			else{
				const index_type edgeIndex = arcIndex;
				const Edge edge = graph().edgeFromId(edgeIndex);
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