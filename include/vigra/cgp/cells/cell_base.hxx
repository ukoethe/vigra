#ifndef CGP2D_CELL_BASE_HXX
#define CGP2D_CELL_BASE_HXX


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

namespace vigra {

struct CellType {
    enum Values {
        Junction =0,
        Boundary =1,
        Region   =2
    };
};



template<class COORDINATE_TYPE,class LABEL_TYPE>
class Cgp;

template<class COORDINATE_TYPE,class LABEL_TYPE,int CELLTYPE>
class Cell;


template<class COORDINATE_TYPE,class LABEL_TYPE,int CELLTYPE>
class CellBase {

    // friend classes
    friend class Cgp<COORDINATE_TYPE,LABEL_TYPE>;
public:

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef vigra::TinyVector<float,2> FloatPointType;
    //typedef BoundingBox<CoordinateType> BoundingBoxType;
    typedef LABEL_TYPE LabelType;
    size_t size()const{
        return points_.size();
    }
    const PointType & operator[](const size_t i) const{
        CGP_ASSERT_OP( i , < , points_.size() );
        return points_[i];
    }
    const PointType & operator[](const size_t i) {
        CGP_ASSERT_OP(i,<,points_.size());
        return points_[i];
    }
    LabelType label()const{
        return label_;
    }


    bool operator == (const CellBase & other){
        return label_==other.label_;
    }
    bool operator != (const CellBase & other){
        return label_!=other.label_;
    }



    const CgpType & cgp()const{
        return *cgp_;
    }
    
    size_t cellType()const{
        return static_cast<size_t>(CELLTYPE);
    }


    std::pair<PointType,PointType> boundingBox()const{
        PointType ul=points_[0];
        PointType lr=points_[0];
        for(size_t p=0;p<size();++p){
            ul[0] = points_[p][0] < ul[0] ? points_[p][0] : ul[0];
            ul[1] = points_[p][1] < ul[1] ? points_[p][1] : ul[1];
            lr[0] = points_[p][0] > lr[0] ? points_[p][0] : lr[0];
            lr[1] = points_[p][1] > lr[1] ? points_[p][1] : lr[1];
        }
        return std::pair<PointType,PointType>(ul,lr);
    }
    
    
    
    FloatPointType centerCoordinate()const{
        FloatPointType cp(0.0f,0.0f);
        for(size_t p=0;p<size();++p){
            cp+=points_[p];
        }
        cp/=size();
        return cp;
    }
    
    const std::vector<LabelType> & bounds()const{
        return bounds_;
    }

    std::vector<LabelType> & __bounds__(){
        return bounds_;
    }



    const std::vector<LabelType> & boundedBy()const{
        return boundedBy_;
    }

    const std::vector<PointType> & points()const{
        return points_;
    }
protected:
    void push_back_point(const PointType & p){
        points_.push_back(p);
    }

    void push_back_bound(const LabelType l){
        bounds_.push_back(l);
    }

    void push_back_bounded_by(const LabelType l){
        boundedBy_.push_back(l);
    }

    void sortAdjaceny(){
        std::sort(bounds_.begin(),bounds_.end());
        std::sort(boundedBy_.begin(),boundedBy_.end());
    }

    void setLabel(const LabelType l){
        label_=l;
    }

    LabelType label_;
    // coordinates
    std::vector<PointType> points_;

    // bounds
    std::vector<LabelType> bounds_;
    std::vector<LabelType> boundedBy_;
    std::vector<LabelType> adjaceny_;

    // cgp pointer
    CgpType * cgp_;
};

}

#endif //CGP2D_CELL_BASE_HXX



