#ifndef CGP2D_CELL_0_HXX
#define CGP2D_CELL_0_HXX

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
class Cell<COORDINATE_TYPE,LABEL_TYPE,0> : public CellBase<COORDINATE_TYPE,LABEL_TYPE,0>{
public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    //typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;

    /*
    void getAngles(const TopologicalGridType & tgrid,const size_t radius,std::vector<float> & angles)const{
        const size_t numBoundaries = this->bounds_.size();
        const int r=static_cast<int>(radius);
        const CoordinateType tx=this->points_[0][0],ty=this->points_[0][1];
        const CoordinateType xmin=    static_cast<int>(tx)-r < 0 ? 0 : static_cast<CoordinateType>(static_cast<int>(tx)-r );
        const CoordinateType ymin=    static_cast<int>(ty)-r < 0 ? 0 : static_cast<CoordinateType>(static_cast<int>(ty)-r );
        const CoordinateType xmax=    static_cast<int>(tx)+r+1 > tgrid.shape(0) ?   tgrid.shape(0) : static_cast<CoordinateType>(static_cast<int>(tx)+r+1 );
        const CoordinateType ymax=    static_cast<int>(ty)+r+1 > tgrid.shape(1) ?   tgrid.shape(1) : static_cast<CoordinateType>(static_cast<int>(ty)+r+1 );

        //std::cout<<"min "<<xmin<<" , "<<ymin<<"\n";
        //std::cout<<"max "<<xmax<<" , "<<ymax<<"\n";

        typedef std::pair<PointType,LabelType>  MapItem;
        typedef std::map<LabelType,MapItem > AverageMapType;
        typedef typename AverageMapType::const_iterator AverageMapConstIter;
        typedef typename AverageMapType::iterator AverageMapIter;

        AverageMapType averageMap;
        for(size_t b=0;b<this->bounds_.size();++b){
            MapItem initItem= MapItem(PointType(0,0),0);
            averageMap[this->bounds_[b]]=initItem;
        }

        // collect
        for(CoordinateType tyy=ymin;tyy<ymax;++tyy)
        for(CoordinateType txx=xmin;txx<xmax;++txx){
            // if boundary
            if(  (txx%2==1 && tyy%2==0) || (txx%2==0 && tyy%2==1) ){

                LabelType cell1Label=tgrid(txx,tyy);
                if(cell1Label!=0){
                    AverageMapIter iter=averageMap.find(cell1Label);
                    if(iter!=averageMap.end()){
                        MapItem item=iter->second;
                        item.first+=PointType(txx,tyy);
                        ++item.second;
                        averageMap[cell1Label]=item;
                    }
                }
            }
        }

        angles.resize(numBoundaries);
        size_t index=0;
        for(AverageMapConstIter iter=averageMap.begin();iter!=averageMap.end();++iter,++index){\
            MapItem item=iter->second;
            PointType averagePoint = item.first;


            averagePoint/=item.second;

            const float x=static_cast<float>(tx);
            const float y=static_cast<float>(ty);
            const float ax=static_cast<float>(averagePoint[0]);
            const float ay=static_cast<float>(averagePoint[1]);


            const float rx=ax-x;
            const float ry=ay-y;


            //std::cout<<" num P "<<item.second<<"\n";
            //std::cout<<"point          "<< x << " , "<< y<<"\n";
            //std::cout<<"averge point   "<< ax<< " , "<<ay<<"\n";
            //std::cout<<"relative point "<< rx<< " , "<<ry<<"\n";
    
            const float result = std::atan2 (ry,rx) * 180.0 / M_PI;
            angles[index]=result;
        }
        //std::cout<<"\n";
    }
    */
};

}

#endif // CGP2D_CELL_0_HXX