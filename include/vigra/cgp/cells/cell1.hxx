#ifndef CGP2D_CELL_1_HXX
#define CGP2D_CELL_1_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>

/* vigra cgp */
#include <vigra/cgp/cells/cell_base.hxx>
#include <vigra/cgp/tgrid.hxx>



namespace vigra{


template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cell<COORDINATE_TYPE,LABEL_TYPE,1> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,1>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef vigra::TinyVector<float,2> FloatPointType;

    typedef LinePath2d<CoordinateType> LinePath2dType;

    void getAngles(){
        angles_.resize(this->size());
        LinePath2dType linePath(this->points_,angles_);            
    }

    float getAngles(const size_t start,const size_t center,const size_t end){
        //std::cout<<"s "<<start<<" c "<<center<<" e "<<end<<"\n";

        FloatPointType   sp   = this->points_[start];
        FloatPointType   akkP(0.0,0.0);
        size_t c=0;
        for(size_t i=start+1;i<end;++i){
            ++c;
            akkP+=this->points_[i];
        }
        CGP_ASSERT_OP(c,==,(end-start-1));
        akkP/=(end-start-1);
        akkP-=sp;

        const float result = std::atan2 (akkP[1],akkP[0]) * 180.0 / M_PI;
        return result;
    }


    void sortCells(){
        size_t finishedPoints=0;
        std::vector<bool>  finishedPoint(this->size(),false);
        std::deque<LabelType> sorted;
        PointType pFront = this->points_[0];
        PointType pBack  = this->points_[0];

        // insert first point 
        sorted.push_back(0);
        finishedPoint[0]=true;
        finishedPoints=1;

        while(finishedPoints < this->size()){
            const size_t oldSize=finishedPoints;
            for(size_t i=0;i<this->size();++i){
                if(finishedPoint[i]==false){
                    const PointType & point = this->points_[i];
                    if(adjacent2Cells(pFront,point)){
                        sorted.push_front(i);
                        pFront=point;
                        finishedPoint[i]=true;
                        ++finishedPoints;
                        break;
                    }
                    else if(adjacent2Cells(pBack,point)){
                        sorted.push_back(i);
                        pBack=point;
                        finishedPoint[i]=true;
                        ++finishedPoints;
                        break;
                    }
                }
            }
            if(oldSize+1!=finishedPoints){
                std::cout<<"\n\n\n\n size "<<this->size()<<"\n";
            }
            CGP_ASSERT_OP(oldSize+1,==,finishedPoints);
        }

        std::vector<PointType> points;
        points.reserve(this->size());
        while(sorted.empty()==false){
            points.push_back(this->points_[sorted.front()]);
            sorted.pop_front();
        }
        this->points_=points;

        this->getAngles();
    }
private:
    bool adjacent2Cells(const PointType & pa,const PointType & pb ){
        
        // A|B
        // vertical  boundary point
        if(pa[0]%2==1){
            // six possible neighbours:
            //
            //    |      case  1
            //  --*--    case 2,3
            //    |               <-self
            //  --*--    case 4,5
            //    |      case  6

            //    
            //    1
            //  2 * 3
            //    |
            //  4 * 5
            //    6

            //case 1 (with border check)
            if     (pa[1]!=0 && ( pa[0]  ==pb[0] && pa[1]-2==pb[1] ) ){return true;}
            //case 2 
            else if(pa[1]!=0 && ( pa[0]-1==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 3 
            else if(pa[1]!=0 && ( pa[0]+1==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 4 
            else if(            ( pa[0]-1==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 5 
            else if(            ( pa[0]+1==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 6 
            else if(            ( pa[0]  ==pb[0] && pa[1]+2==pb[1] ) ){return true;}
            
            return false;
        }
        // horizontal boundary
        else{
            // six possible neighbours:
            //
            //      |  |        
            //    --*--*--      
            //      |  | 
            //          
            //      2  4
            //    1 *--* 6
            //      3  5 

            //case 1 (with border check)
            if     (pa[0]!=0  && ( pa[0] -2 ==pb[0] && pa[1]  ==pb[1] ) ){return true;}
            //case 2 
            else if(pa[0]!=0  && ( pa[0] -1 ==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 3 
            else if(pa[0]!=0  && ( pa[0] -1 ==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            //case 4 
            else if(             ( pa[0] +1 ==pb[0] && pa[1]-1==pb[1] ) ){return true;}
            //case 5 
            else if(             ( pa[0] +1 ==pb[0] && pa[1]+1==pb[1] ) ){return true;}
            else if(             ( pa[0] +2 ==pb[0] && pa[1]  ==pb[1] ) ){return true;}

            return false;
        }
    }
public:
    std::vector<vigra::TinyVector<float,2> > angles_;
};

}

#endif // CGP2D_CELL_1_HXX