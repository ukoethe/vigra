#ifndef VIGRA_DISTANCES_HXX
#define VIGRA_DISTANCES_HXX

#include <vigra/numerictraits.hxx>
#include <vigra/multi_array.hxx>

namespace vigra{
namespace distances{



	
	template<bool IS_FUNDAMENTAL>
	struct TagSelctor;

	template<>
	struct TagSelector<true>{
		typename ScalarTag TagType;
	};

	template<>
	struct TagSelector<false>{
		typename ScalarTag NotScalarTag;
	};





	template<class T>
	struct ChiSquaredNew{
		ChiSquaredNew(){}

	public:


		// two scalars
	private:
		template<class A,class B>
		T op_impl(const A & a ,const B & b,const SclarTag){
			const T sum  = a + b;
			const T diff = a - b; 
			diff*diff/sum*static_cast<T>(0.5);
			return
		}

		// multiarray 1d
		template<class A,class B>
		T op_impl()(const vigra::MultiArrayView<1,A> & a ,const vigra::MultiArrayView<1,B> & b , const NotScalarTag)const{

			const size_t nA=a.shape(0);
			const size_t nB=b.shape(0);
			T res = 0.0;
			for(size_t i=0;i<nA;++i){
				const T aa=static_cast<T>(a(i));
				const T bb=static_cast<T>(b(i));
				const T sum  = aa + bb;
				const T diff = aa - bb; 
				if(sum<= static_cast<T>(0.0000001)){
				}
				else{
					res+=(diff*diff)/sum;
				}
			}
			return res*static_cast<T>(0.5);
		}

	public:





		template<class A,class B>
		T operator()(const A & a ,const B & b)const{
			return op_impl(a,b, 
				typename TagSelector<  std::numeric_limits<A>::is_specialized &&  std::numeric_limits<A>::is_specialized  >::TagType
			);
		}

	};


	/*
	template<class T>
	struct ChiSquared{
		ChiSquared(){}

		template<class A,class B>
		T operator()(const A & a ,const B & b)const{

			const size_t nA=a.shape(0);
			const size_t nB=b.shape(0);
			T res = 0.0;
			for(size_t i=0;i<nA;++i){
				const T aa=static_cast<T>(a(i));
				const T bb=static_cast<T>(b(i));
				const T sum  = aa + bb;
				const T diff = aa - bb; 
				if(sum<= static_cast<T>(0.0000001)){
				}
				else{
					res+=(diff*diff)/sum;
				}
			}
			return res*T(0.5);
		}

	};
	*/


	template<class T>
	struct SquaredNorm{
		SquaredNorm(){}

		template<class A,class B>
		T operator()(const A & a ,const B & b)const{

			const size_t nA=a.shape(0);
			const size_t nB=b.shape(0);
			T res = 0.0;
			for(size_t i=0;i<nA;++i){
				const T aa=static_cast<T>(a(i));
				const T bb=static_cast<T>(b(i));
				const T diff = aa-bb;
				res+=diff*diff;	
			}
			return res;
		}

	};

	template<class T>
	struct Norm{
		Norm(){}

		template<class A,class B>
		T operator()(const A & a ,const B & b)const{

			const size_t nA=a.shape(0);
			const size_t nB=b.shape(0);
			T res = 0.0;
			for(size_t i=0;i<nA;++i){
				const T aa=static_cast<T>(a(i));
				const T bb=static_cast<T>(b(i));
				const T diff = aa-bb;
				res+=diff*diff;	
			}
			return std::sqrt(res);
		}

	};

} // end namespace distances
} // end namepsace vigra


#endif //VIGRA_DISTANCES_HXX