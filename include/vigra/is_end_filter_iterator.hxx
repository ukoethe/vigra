#ifndef VIGRA_IS_END_FILTER_ITERATOR_HXX
#define VIGRA_IS_END_FILTER_ITERATOR_HXX

#include <vigra/graphs.hxx>
#include <boost/iterator/filter_iterator.hpp>

namespace vigra{



    template<class FILTTER,class ITER>
    class FilterIter
    :   public boost::filter_iterator< FILTTER,ITER>
    {
    public:
        FilterIter()
        :   boost::filter_iterator< FILTTER,ITER>(){
        }
        FilterIter(const FILTTER & filter, ITER  pos, ITER end)
        :   boost::filter_iterator< FILTTER,ITER>(filter,pos,end){
        }
        bool isEnd()const{
            //return *this == this->end();
            return this->base().isEnd();
        }

    };


	template<class FILTER,class ITER>
	inline bool operator == (const FilterIter<FILTER,ITER> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator == (const lemon::Invalid & iv,const FilterIter<FILTER,ITER> & iter){
    	return iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator != (const FilterIter<FILTER,ITER> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class FILTER,class ITER>
	inline bool operator != (const lemon::Invalid & iv,const FilterIter<FILTER,ITER> & iter){
    	return !iter.isEnd();
	}



} // end namespace vigra






#endif //VIGRA_IS_END_FILTER_ITERATOR_HXX