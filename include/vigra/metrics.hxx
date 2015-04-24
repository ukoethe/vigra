/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/
#ifndef VIGRA_METRIC_HXX
#define VIGRA_METRIC_HXX

#include <vigra/numerictraits.hxx>
#include <vigra/multi_array.hxx>

#include <cmath>

namespace vigra{
namespace metrics{


    template<class T>
    class ChiSquared{
    public:
        ChiSquared(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A, class B>
        T operator()(const A & a,const B & b)const{
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
        template<class A, class B>
        T operator()(const A & a,const B & b)const{
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
                res+=diff*diff;
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
        template<class A, class B>
        T operator()(const A & a,const B & b)const{
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
                res+= std::abs(std::pow((double)diff,(int)NORM));
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

    template<class T>
    class SymetricKlDivergenz{
    public:
        SymetricKlDivergenz(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A, class B>
        T operator()(const A & a,const B & b)const{
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
                const T val = std::log(aa/bb)*(aa - bb);
                if(!std::isinf(val) && !std::isnan(val))
                    res+=val;
                ++iterA;
                ++iterB;
            }
            return res/static_cast<T>(2.0);
        }
    };

    template<class T>
    class BhattacharyaDistance{
    public:
        BhattacharyaDistance(){}
        T operator()(const T & a,const T & b)const{
            return opImpl(&a,&a+1,&b);
        }
        template<class A, class B>
        T operator()(const A & a,const B & b)const{
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
                res+=std::sqrt(aa*bb);
                ++iterA;
                ++iterB;
            }
            return std::sqrt( static_cast<T>(1.0)-res);
        }
    };

    enum MetricType{
        ChiSquaredMetric=0,
        HellingerMetric=1,
        SquaredNormMetric=2,
        NormMetric=3,
        ManhattanMetric=4,
        SymetricKlMetric=5,
        BhattacharyaMetric=6
    };


    template<class T>
    class Metric{
    public:

        Metric(const MetricType metricType = ManhattanMetric)
        : metricType_(metricType){

        }

        template<class A, class B>
        T operator()(const A & a,const B & b)const{
            switch(static_cast<unsigned int>(metricType_)){
                case 0:
                    return chiSquared_(a,b);
                case 1:
                    return hellingerDistance_(a,b);
                case 2:
                    return squaredNorm_(a,b);
                case 3:
                    return norm_(a,b);
                case 4:
                    return manhattan_(a,b);
                case 5:
                    return symetricKlDivergenz_(a,b); 
                case 6 :
                    return bhattacharyaDistance_(a,b); 
                default :
                    return 0;
            }
        } 
    private:
        MetricType metricType_;
        ChiSquared<T> chiSquared_;
        HellingerDistance<T> hellingerDistance_;
        SquaredNorm<T> squaredNorm_;
        Norm<T> norm_;
        Manhattan<T> manhattan_;
        SymetricKlDivergenz<T> symetricKlDivergenz_;
        BhattacharyaDistance<T> bhattacharyaDistance_;
    };

} // end namespace metric
} // end namepsace vigra


#endif //VIGRA_METRIC_HXX
