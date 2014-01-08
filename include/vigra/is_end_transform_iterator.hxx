#ifndef VIGRA_IS_END_TRANSFORM_ITERATOR_HXX
#define VIGRA_IS_END_TRANSFORM_ITERATOR_HXX

#include <vigra/graphs.hxx>
#include <boost/iterator/transform_iterator.hpp>

namespace vigra{

    template<class TRANSFORM,class ITER,class REF,class VAL>
    class NewTransformIter
    :   public boost::transform_iterator< TRANSFORM,ITER,REF,VAL>
    {
    public:
        NewTransformIter()
        :boost::transform_iterator< TRANSFORM,ITER,REF,VAL>(){
        }
        NewTransformIter(ITER  pos ,const TRANSFORM & transform)
        :   boost::transform_iterator< TRANSFORM,ITER,REF,VAL>(pos,transform){
        }
    };


    template<class TRANSFORM,class ITER,class REF,class VAL>
    class TransformIter
    :   public boost::transform_iterator< TRANSFORM,ITER,REF,VAL>
    {
    public:
        TransformIter()
        :boost::transform_iterator< TRANSFORM,ITER,REF,VAL>(){
        }
        TransformIter(ITER  pos ,const TRANSFORM & transform)
        :   boost::transform_iterator< TRANSFORM,ITER,REF,VAL>(pos,transform){
        }
        bool isEnd()const{
            return this->base().isEnd();
        }

    };


	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator == (const TransformIter<TRANSFORM,ITER,REF,VAL> & iter,const lemon::Invalid & iv){
    	return iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator == (const lemon::Invalid & iv,const TransformIter<TRANSFORM,ITER,REF,VAL> & iter){
    	return iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator != (const TransformIter<TRANSFORM,ITER,REF,VAL> & iter,const lemon::Invalid & iv){
    	return !iter.isEnd();
	}
	template<class TRANSFORM,class ITER,class REF,class VAL>
	inline bool operator != (const lemon::Invalid & iv,const TransformIter<TRANSFORM,ITER,REF,VAL> & iter){
    	return !iter.isEnd();
	}


} // end namespace vigra






#endif //VIGRA_IS_END_TRANSFORM_ITERATOR_HXX