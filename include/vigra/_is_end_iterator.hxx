#ifndef VIGRA_IS_END_ITERATOR_HXX
#define VIGRA_IS_END_ITERATOR_HXX

#include <vigra/graphs.hxx>

namespace vigra{

    template<class ITER>
    class IsEndIter
    :   public ITER {
    public:
        IsEndIter()
        :   ITER(),
            end_(),
            valid_(false){            
        }

        IsEndIter(ITER iter , ITER end)
        :   ITER(iter),
            end_(end),
            valid_(true){
        }

        bool isEnd()const{
            return !valid_ || *this==end_;
        }
    private:
        ITER end_;
        bool valid_;
    };





	template<class ITER>
	inline bool operator == (const IsEndIter<ITER> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class ITER>
	inline bool operator == (const lemon::Invalid & iv,const IsEndIter<ITER> & iter){
    	return iter.isEnd();
	}
	template<class ITER>
	inline bool operator != (const IsEndIter<ITER> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class ITER>
	inline bool operator != (const lemon::Invalid & iv,const IsEndIter<ITER> & iter){
    	return !iter.isEnd();
	}





} // end namespace vigra






#endif //VIGRA_IS_END_ITERATOR_HXX