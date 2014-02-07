#ifndef VIGRA_METRIC_HXX
#define VIGRA_METRIC_HXX

#include <vigra/numerictraits.hxx>
#include <vigra/multi_array.hxx>

namespace vigra{
namespace metrics{



	
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


	template<class T>
	struct Manhattan{
		Manhattan(){}

		template<class A,class B>
		T operator()(const A & a ,const B & b)const{

			const size_t nA=a.shape(0);
			const size_t nB=b.shape(0);
			T res = 0.0;
			for(size_t i=0;i<nA;++i){
				const T aa=static_cast<T>(a(i));
				const T bb=static_cast<T>(b(i));
				const T diff = aa-bb;
				res+=std::fabs(diff);	
			}
			return res;
		}

	};

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

} // end namespace metric
} // end namepsace vigra


#endif //VIGRA_METRIC_HXX