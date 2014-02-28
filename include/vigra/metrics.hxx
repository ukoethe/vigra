#ifndef VIGRA_METRIC_HXX
#define VIGRA_METRIC_HXX

#include <vigra/numerictraits.hxx>
#include <vigra/multi_array.hxx>

namespace vigra{
namespace metrics{





    template<class T>
    class ChiSquared{
    public:
        ChiSquared(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A>
        T operator()(const A & a,const A & b)const{
            return opImpl(a.begin(),a.end(),b.begin());
        } 
    private:
        template<class ITER_A,class ITER_B>
        T opImpl(
            ITER_A  iterA  ,ITER_A  endA   ,ITER_B  iterB 
        )const{   
            T res = 0.0;
            while(iterA!=endA){
                const T aa=static_cast<T>(*iterA);
                const T bb=static_cast<T>(*iterB);
                const T sum  = aa + bb;
                const T diff = aa - bb; 
                if(sum> static_cast<T>(0.0000001))
                    res+=(diff*diff)/sum;
                ++iterA;
                ++iterB;
            }
            return res*T(0.5);
        }
    };


    template<class T>
    class HellingerDistance{
    public:
        HellingerDistance(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A>
        T operator()(const A & a,const A & b)const{
            return opImpl(a.begin(),a.end(),b.begin());
        } 
    private:
        template<class ITER_A,class ITER_B>
        T opImpl(
            ITER_A  iterA  ,ITER_A  endA   ,ITER_B  iterB 
        )const{   
            T res = 0.0;
            while(iterA!=endA){
                const T aa=std::sqrt(static_cast<T>(*iterA));
                const T bb=std::sqrt(static_cast<T>(*iterB));
                const T diff = aa - bb; 
                sum+=diff*diff;
                ++iterA;
                ++iterB;
            }
            return std::sqrt(res)/std::sqrt(2.0);
        }
    };
    
    template<class T,unsigned int NORM,bool TAKE_ROOT=true>
    class PNorm{
    public:
        PNorm(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A>
        T operator()(const A & a,const A & b)const{
            return opImpl(a.begin(),a.end(),b.begin());
        } 
    private:
        template<class ITER_A,class ITER_B>
        T opImpl(
            ITER_A  iterA  ,ITER_A  endA   ,ITER_B  iterB 
        )const{
            T res = static_cast<T>(0.0);
            while(iterA!=endA){
                const T aa=static_cast<T>(*iterA);
                const T bb=static_cast<T>(*iterB);
                const T diff = aa-bb;
                res+= std::abs(std::pow(diff,NORM));  
                ++iterA;
                ++iterB;
            }
            return TAKE_ROOT ? std::pow(res,static_cast<T>(1)/static_cast<T>(NORM)) : res;
        }
    };

    template<class T>
    class SquaredNorm
    :  public PNorm<T,2,false>{
    public:
        SquaredNorm()
        :   PNorm<T,2,false>(){
        }
    };

    template<class T>
    class Norm
    :  public PNorm<T,2,true>{
    public:
        Norm()
        :   PNorm<T,2,true>(){
        }
    };

    template<class T>
    class Manhattan
    :  public PNorm<T,1,false>{
    public:
        Manhattan()
        :   PNorm<T,1,false>(){
        }
    };


} // end namespace metric
} // end namepsace vigra


#endif //VIGRA_METRIC_HXX