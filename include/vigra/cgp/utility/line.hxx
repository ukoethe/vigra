#ifndef CGP2D_LINE_HXX
#define CGP2D_LINE_HXX


#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>


namespace vigra {

/****************************************************************************/
/* L i n e 2 d                                                              */
/****************************************************************************/

template<class T>
class Line2d{
public:
    typedef T ValueType;
    typedef vigra::TinyVector<ValueType,2> PointType;
    Line2d(){}
    template<class P>
    Line2d(const P & pA,const P & pB)
    :   pA_(pA),
        pB_(pB){

    }

    template<class ITER>
    ValueType accDistance( ITER begin,ITER end ){
        ValueType dist = 0.0 ;
        while(begin!=end){
            dist+=distance(*begin);
            ++begin;
        }
        return dist;
    }

    template<class V>
    ValueType distance(const vigra::TinyVector<V,2> & p){
        PointType n = pA_ - pB_;
        n/=vigra::norm(n);
        PointType ap = pA_-p;
        const float apn =  vigra::dot(ap,n);
        return vigra::norm(ap - (apn)*n);
    }


    vigra::TinyVector<float,2> angle()const{
        vigra::TinyVector<float,2> r=(pA_- pB_);
        r/=vigra::norm(r);
        return r;
    }
private:
    PointType pA_;
    PointType pB_;
};  

/****************************************************************************/
/* L i n e P a t h 2 d                                                      */
/****************************************************************************/

template<class T>
class LinePath2d{
public:
    typedef vigra::TinyVector<T,2>         PointType;
    typedef std::vector<PointType>  PointVector;
    typedef Line2d<float>           FLine2d;
    typedef std::vector<FLine2d>    FLine2dVec;

    LinePath2d(const PointVector & inPoints,std::vector<vigra::TinyVector<float,2> >  & angles)
    :   inPoints_(inPoints),
        lineVec_(),
        angles_(angles)
    {
        const size_t numPoints=inPoints_.size();

        if (numPoints==1){
            const PointType & p = inPoints_.front();
            // A|B
            // vertical  boundary point
            if(p[0]%2==1){
                FLine2d line(PointType(0,0),PointType(0,1));
                angles_[0]=line.angle();
            }
            else{
                FLine2d line(PointType(0,0),PointType(1,0));
                angles_[0]=line.angle();
            }
        }
        else if (numPoints<=3){  
            FLine2d line(inPoints_.front(),inPoints_.back());
            for(size_t i=0;i<numPoints;++i){
                angles_[i]=line.angle();
            }
            const float dist = line.accDistance(inPoints_.begin(),inPoints_.end());
            //std::cout<<"tiny ! dist "<<dist<<"\n";
        }
        else{
            //std::cout<<"to big\n";
            size_t end = 0;
            while(end!=inPoints_.size()){
                //std::cout<<"-----\n";
                size_t start = end;
                end = this->getAngle(start);
            }
            
        }



    }

    size_t getAngle(const size_t start){
        //std::cout<<"start "<<start<<"\n";
        size_t end = inPoints_.size();
        const float  threshold = 2.0;
        float distance = 9999999.0f;
        size_t nP  = end -start;
        while(true) {

            //std::cout<<" * end "<<end<<"\n";
            CGP_ASSERT_OP(nP,>=,2);
            CGP_ASSERT_OP(start,!=,end-1);
            FLine2d line(inPoints_[start],inPoints_[end-1]);
            distance = line.accDistance(inPoints_.begin()+start,inPoints_.begin()+end);
            if(nP<=3 || distance < threshold){
                for(size_t i=start;i<end;++i){
                    angles_[i]=line.angle();
                }
                break;
            }
            else{
                end = start+nP/2;
                nP  = end -start;
            }
        }
        //std::cout<<"at ending where end is "<<end<<"\n";
        return end; 
    }
private:
    const PointVector &        inPoints_;
    FLine2dVec                 lineVec_;
    std::vector<vigra::TinyVector<float,2> >  &       angles_;

};

}

#endif //CGP2D_LINE_HXX