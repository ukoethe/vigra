#ifndef CGP2D_HXX
#define CGP2D_HXX

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

/* opengm */
//#include "opengm/config.hxx"
//#include "opengm/utilities/metaprogramming.hxx"

/* this project */
#include "partition.hxx"
#include "macros.hxx"

#include "utility/line.hxx"
#include "cells/cell_base.hxx"
#include "cells/cell0.hxx"
#include "cells/cell1.hxx"
#include "cells/cell2.hxx"
#include "tgrid.hxx"


namespace vigra{


template<class COORDINATE_TYPE,class LABEL_TYPE>
class Geometry;




template<class COORDINATE_TYPE,class LABEL_TYPE,int CELLTYPE>
class Cell;



/****************************************************************************/
/* C g p                                                                    */
/****************************************************************************/




template<int CELL_TYPE,class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper{
};


template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<0,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,0>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry0(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry0(); }
};




template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<1,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,1>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry1(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry1(); }
};



template<class COORDINATE_TYPE, class LABEL_TYPE>
struct CgpHelper<2,COORDINATE_TYPE,LABEL_TYPE>{

    typedef Cgp<COORDINATE_TYPE,LABEL_TYPE> CgpType;

    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,2>  CellType;
    typedef std::vector< CellType >             CellsType;

    static CellsType &       getCells(CgpType       & cgp){ return cgp.geometry2(); }
    static const CellsType & getCells(const CgpType & cgp){ return cgp.geometry2(); }
};





template<class COORDINATE_TYPE, class LABEL_TYPE>
class Cgp {
    
    public:
    typedef LABEL_TYPE LabelType;
    typedef COORDINATE_TYPE CoordinateType;
    typedef TopologicalGrid<LabelType> TopologicalGridType;
    typedef typename TopologicalGridType::ShapeType ShapeType;
    typedef vigra::TinyVector<CoordinateType,2> PointType;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,0> Cell0;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,1> Cell1;
    typedef Cell<COORDINATE_TYPE,LABEL_TYPE,2> Cell2;
    typedef std::vector< Cell0 > Cells0;
    typedef std::vector< Cell1 > Cells1;
    typedef std::vector< Cell2 > Cells2;

    typedef std::vector<std::vector< LABEL_TYPE> > CellAdjacencyGraphVectorType;
    typedef std::vector<std::set< LABEL_TYPE> > CellAdjacencyGraphSetType;

    // Constructor
    Cgp(const TopologicalGridType & tgrid );
    // Query
    const Cells0  & geometry0()const;
    const Cells1  & geometry1()const;
    const Cells2  & geometry2()const;

    const TopologicalGridType & tgrid()const;

    size_t numCells(const size_t cellType)const{
        return tgrid_.numCells(cellType);
    }

    size_t shape(const size_t d)const{
        return tgrid_.shape(d);
    }

    const ShapeType & shapeTopologicalGrid()const{
        return tgrid_.shapeTopologicalGrid();
    }

    const ShapeType & shapeLabeling()const{
        return tgrid_.shapeLabeling();
    }


    LabelType operator()(const size_t x,const size_t y) const{
        return tgrid_(x,y);
    }

    LabelType operator()(const size_t x,const size_t y)     {
        return tgrid_(x,y);
    }
    


    template<int CELL_TYPE>
    const std::vector< Cell<COORDINATE_TYPE,LABEL_TYPE,CELL_TYPE> > & cells()const{
        return CgpHelper<CELL_TYPE,CoordinateType,LabelType>::getCells(*this);
    }

    template<int CELL_TYPE>
    std::vector< Cell<COORDINATE_TYPE,LABEL_TYPE,CELL_TYPE> > & cells(){
        return CgpHelper<CELL_TYPE,CoordinateType,LabelType>::getCells(*this);
    }



    template<int CELL_TYPE>
    size_t cellSizeT(const LabelType cellIndex)const{
        return this-> template cells<CELL_TYPE> ()[cellIndex].size();
    }


    size_t cellSize(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                return this-> template cells<0>()[cellIndex].size();
            case 1:
                return this-> template cells<1>()[cellIndex].size();
            case 2:
                return this-> template cells<2>()[cellIndex].size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }



    size_t nBounds(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                return this-> template cells<0>().bounds().size();
            case 1:
                return this-> template cells<1>().bounds().size();
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }
    
    template<int CELL_TYPE>
    size_t nBoundsT(const LabelType cellIndex)const{
        if(CELL_TYPE==2){
            CGP_ASSERT_OP(false,==,true);
        }
        return this-> template cells<CELL_TYPE>().bounds().size();

    }


    size_t nBoundedBy(const size_t cellType ,const LabelType cellIndex)const{
        switch(cellType){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return Cells1_[cellIndex].boundedBy().size();
            case 2:
                return Cells2_[cellIndex].boundedBy().size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }
    template<int CELL_TYPE>
    size_t nBoundedByT(const LabelType cellIndex)const{
        switch(CELL_TYPE){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return Cells1_[cellIndex].boundedBy().size();
            case 2:
                return Cells2_[cellIndex].boundedBy().size();
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    LabelType bound(const size_t cellType ,const LabelType cellIndex,const size_t boundNr)const{
        switch(cellType){
            case 0:
                return Cells0_[cellIndex].bounds()[boundNr];
            case 1:
                return Cells1_[cellIndex].bounds()[boundNr];
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    template<int CELL_TYPE>
    LabelType bound(const LabelType cellIndex,const size_t boundNr)const{
        switch(CELL_TYPE){
            case 0:
                return Cells0_[cellIndex].bounds()[boundNr];
            case 1:
                return Cells1_[cellIndex].bounds()[boundNr];
            case 2:
                CGP_ASSERT_OP(false,==,true);
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }


    LabelType boundedBy(const size_t cellType ,const LabelType cellIndex,const size_t boundedByNr)const{
        switch(cellType){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return Cells1_[cellIndex].boundedBy()[boundedByNr];
            case 2:
                return Cells2_[cellIndex].boundedBy()[boundedByNr];
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }

    template<int CELL_TYPE>
    LabelType boundedBy(const LabelType cellIndex,const size_t boundedByNr)const{
        switch(CELL_TYPE){
            case 0:
                CGP_ASSERT_OP(false,==,true);
            case 1:
                return Cells1_[cellIndex].boundedBy()[boundedByNr];
            case 2:
                return Cells2_[cellIndex].boundedBy()[boundedByNr];
            default :
                CGP_ASSERT_OP(false,==,true);
        }
    }



    void cellAdjacencyGraphVector(int cellType,CellAdjacencyGraphVectorType & graph)const{
        CGP_ASSERT_OP(cellType,>=,1);
        CGP_ASSERT_OP(cellType,<=,2);
        const size_t numCells=this->numCells(cellType);
        graph.clear();
        graph.resize(numCells);

        CellAdjacencyGraphSetType  graphSet;
        this->cellAdjacencySetGraph(cellType,graphSet);
        for(size_t ci=0;ci<numCells;++ci){
            graph[ci].assign(graphSet[ci].begin(),graphSet[ci].end());
        }
    }


    void cellAdjacencySetGraph(int cellType,CellAdjacencyGraphSetType & graph)const{
        CGP_ASSERT_OP(cellType,>=,1);
        CGP_ASSERT_OP(cellType,<=,2);

        const size_t numCells=this->numCells(cellType);
        graph.clear();
        graph.resize(numCells);

        if(cellType==1){

            const size_t nBoundaries = numCells;
            // iterate over all boundaries
            for(size_t bi=0;bi<nBoundaries;++bi){

                // iterate over all junctions which bound the boudary
                const size_t nBoundedBy = Cells1_[bi].boundedBy().size();
                CGP_ASSERT_OP(nBoundedBy,<=,2);

                for(size_t j=0;j<nBoundedBy;++j){
                    const LabelType  ji = Cells1_[bi].boundedBy()[j]-1.0;

                    // get the number boundaries for the junction
                    const size_t nBounds = Cells0_[ji].bounds().size();
                    CGP_ASSERT_OP(nBounds,>=,3);
                    CGP_ASSERT_OP(nBounds,<=,4);

                    for(size_t b=0;b<nBounds;++b){

                        // insert other region to to adj.
                        const LabelType biOther = Cells0_[ji].bounds()[b];
                        if(biOther!=bi){
                            graph[bi].insert(biOther);
                        }
                    }
                }
            }
        }

        if(cellType==2){
            const size_t nRegion = numCells;

            // iterate over all regions
            for(size_t ri=0;ri<nRegion;++ri){

                // iterate over all faces which bound the region
                const size_t nBoundedBy = Cells2_[ri].boundedBy().size();
                const bool assertTrue  = nBoundedBy>0 || nRegion==1;
                CGP_ASSERT_OP(assertTrue,==,true);

                for(size_t b=0;b<nBoundedBy;++b){
                    const LabelType  bi = Cells2_[ri].boundedBy()[b]-1.0;

                    // get the 2 regions of the boundary
                    const size_t nBounds = Cells1_[bi].bounds().size();
                    CGP_ASSERT_OP(nBounds,==,2);
                    for(size_t r=0;r<2;++r){

                        // insert other region to to adj.
                        const LabelType riOther = Cells1_[bi].bounds()[r];
                        if(riOther!=ri){
                            graph[ri].insert(riOther);
                        }
                    }
                }
            }
        }
    }


    /** given a consecutive list of topological points
     *  (that define a line in 2D),
     *  return a list of cartesian coordinates which can be used
     *  to draw the line.
     */
    std::vector<PointType> cartesianLine(const Cell1& c1) const {
        std::vector<PointType> res;
        for(int x=0; x<c1.size(); ++x) {
            const typename Cell1::PointType& pt = c1[x];

            PointType p1; //start
            PointType p2; //end
            if(pt[0] % 2 == 0) { //normal in y direction
                p1[0] = (pt[0])/2;
                p2[0] = (pt[0])/2+1;
                p1[1] = (pt[1]+1)/2;
                p2[1] = (pt[1]+1)/2;
            }
            else { //normal in x direction
                p1[0] = (pt[0]+1)/2;
                p2[0] = (pt[0]+1)/2;
                p1[1] = (pt[1])/2;
                p2[1] = (pt[1])/2+1;
            }
            res.push_back(p1);
            res.push_back(p2);
        }
        
        if(res.size() > 2) {
            //sort points such that they form a consecutive line.
            //to do this, look at consecutive pairs of points and reorder them.
            for(int x=0; x<res.size()-2; x+=2) {
                if(res[x+0] == res[x+3]) {
                    std::swap(res[x+0], res[x+1]);
                    std::swap(res[x+2], res[x+3]);
                }
                else if(res[x+1] == res[x+2]) {
                    continue;
                }
                else if(res[x+0] == res[x+2]) {
                    std::swap(res[x+0], res[x+1]);
                }
                else if(res[x+1] == res[x+3]) {
                    std::swap(res[x+2], res[x+3]);
                }
                else {
                    throw std::runtime_error("path is not connected");
                }
                if(res[x+1] != res[x+2]) {
                    throw std::runtime_error("err");
                }
            }
        }
        
        std::vector<PointType> res2;
        for(int x=0; x<res.size(); x+=2) {
            if(x==res.size()-2) {
                res2.push_back(res[x]);
                res2.push_back(res[x+1]);
            }
            else {
                res2.push_back(res[x]); 
            }
        }
        
        return res2;
    }
    
    std::vector<unsigned int> serialize() const {
        //do not wonder that this code looks odd,
        //it was copied and adapted from the 3D code
        
        typedef std::vector<const Cell1*> Lines;
        typedef std::map<LABEL_TYPE, Lines> Label2Lines;
        typedef unsigned int data_t;
        typedef std::vector<unsigned int> SerializeType; 
        
        std::map<LABEL_TYPE, std::vector<const Cell1*> > label2lines;
        for(typename Cells1::const_iterator it=Cells1_.begin(); it!=Cells1_.end(); ++it) {
            label2lines[it->label()].push_back(&(*it));
        }
        
        //determine size of serialized data
        int numPoints = 0;
        for(typename Label2Lines::const_iterator it=label2lines.begin(); it!=label2lines.end(); ++it) {
            numPoints += 2;               //<label><size>
            for(typename Lines::const_iterator l= it->second.begin();
                l!=it->second.end(); ++l)
            {
                numPoints += 2*((*l)->size() +1);
            }
            numPoints += 2*it->second.size()-2; //extra room for seperator tokens
        }

        //array in which to serialize into
        SerializeType M(numPoints); 

        data_t* m    = &M[0];
        data_t* mEnd = m + numPoints;

        const data_t sepToken = std::numeric_limits<data_t>::max();

        //now serialize the whole slice into the array
        data_t* k;
        for(typename Label2Lines::const_iterator it = label2lines.begin(); it!=label2lines.end(); ++it) {
            const LABEL_TYPE& label = it->first;
            const std::vector<const Cell1*>& lines = it->second;

            //write label
            *m = label; ++m;
            
            //indicate following length
            k = m; //save position
            ++m;

            for(typename Lines::const_iterator l=lines.begin(); l!=lines.end(); ++l) {
                const Cell1& c1 = **l;
                
                std::vector<PointType> line = cartesianLine(c1);
                
                for(int x=0; x<line.size(); ++x) {
                    const typename Cell1::PointType& pt = line[x];
                    
                    *m = pt[0]; ++m;
                    *m = pt[1]; ++m;
                }
                if(l!=(--lines.end())) {
                    //write seperator token
                    *m = sepToken; ++m;
                    *m = sepToken; ++m;
                }
            }
            if( (m-k-1) % 2 != 0) {throw std::runtime_error("not good"); }
            *(k)   = m-k-1;

        }
        if(m != mEnd) throw std::runtime_error("Oh no!");
        
        return M;
    }

    template<class CELL1_STATE_ITER>
    void merge2Cells (
        CELL1_STATE_ITER statesBegin,
        CELL1_STATE_ITER statesEnd,
        TopologicalGridType & newTGrid
    ) const {
        typedef partition::Partition<LabelType>    UfdType;

        const size_t numCell1=numCells(1);
        const size_t numCell2=numCells(2);

        UfdType cell2Ufd(numCell2);
        // merging
        for(size_t  cell1Index=0;cell1Index<numCell1;++cell1Index,++statesBegin){
            const size_t state = static_cast<size_t>(*statesBegin);
            const bool mergeCell2= state==0;
            if(mergeCell2){
                const LabelType cell2IndexA=Cells1_[cell1Index].bounds()[0]-1;
                const LabelType cell2IndexB=Cells1_[cell1Index].bounds()[1]-1;
                cell2Ufd.merge(cell2IndexA,cell2IndexB);
            }
        }
        std::map<LabelType,LabelType> denseRelabeling;
        cell2Ufd.representativeLabeling(denseRelabeling);
        //std::cout<<"number of sets "<<cell2Ufd.numberOfSets()<<"\n";
        typedef vigra::MultiArray<2,LabelType>  LabelImageType;
        typedef typename LabelImageType::difference_type ShapeType;

        ShapeType labelingShape( (shape(0)+1)/2,(shape(1)+1)/2 );
        LabelImageType newLabeling(labelingShape);

        for(size_t y=0;y<labelingShape[1];++y)
        for(size_t x=0;x<labelingShape[0];++x){
            // get old label from tgrid;
            const LabelType oldCell2Label=tgrid_(x*2,y*2);
            const LabelType oldCell1Index=oldCell2Label-1;
            // find old index in udf
            const LabelType foundIndex=cell2Ufd.find(oldCell1Index);
            // relabel to dense new INDEX
            const LabelType newCell2Index=denseRelabeling[foundIndex];
            // finaly the new label
            const LabelType newCell2Label=newCell2Index+1;
            newLabeling(x,y)=newCell2Label;
        }
        newTGrid=TopologicalGridType(newLabeling);
    }

    private:
    const TopologicalGridType & tgrid_;
    std::vector< Cell0 >  Cells0_;
    std::vector< Cell1 >  Cells1_;
    std::vector< Cell2 >  Cells2_;
};



// Implementation Cgp
template<class COORDINATE_TYPE,class LABEL_TYPE>
Cgp<COORDINATE_TYPE,LABEL_TYPE>::Cgp(const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::TopologicalGridType  & tgrid )
:   tgrid_(tgrid),
    Cells0_(tgrid.numCells(0)),
    Cells1_(tgrid.numCells(1)),
    Cells2_(tgrid.numCells(2))
{
    // set up geometry
    const typename TopologicalGridType::LabelImageType & grid=tgrid.tgrid();
    for(size_t ty=0;ty<tgrid.shape(1);++ty)
    for(size_t tx=0;tx<tgrid.shape(0);++tx){
        // Cell 2 (=Region)
        if(tx%2==0 && ty%2==0){
            int label = grid(tx,ty);
            CGP_ASSERT_OP(label,>,0)
            CGP_ASSERT_OP(label,<=,tgrid.numCells(2));
            Cells2_[label-1].push_back_point(PointType(tx,ty));
        }
        // Cell 0 (== Junction)
        else if(tx%2!=0 && ty%2!=0){
            int label = grid(tx,ty);
            if(label!=0){
                CGP_ASSERT_OP(label,>,0)
                CGP_ASSERT_OP(label,<=,tgrid.numCells(0));
                Cells0_[label-1].push_back_point(PointType(tx,ty));
            }
                
        }
        // Cell 1 (== Boundary)
        else{
            int label = grid(tx,ty);
            if(label!=0){
                CGP_ASSERT_OP(label,>,0)
                CGP_ASSERT_OP(label,<=,tgrid.numCells(1));
                Cells1_[label-1].push_back_point(PointType(tx,ty));
            }
        }
    }
    // check size of geometry
    CGP_ASSERT_OP(Cells0_.size(),==,tgrid.numCells(0));
    CGP_ASSERT_OP(Cells1_.size(),==,tgrid.numCells(1));
    CGP_ASSERT_OP(Cells2_.size(),==,tgrid.numCells(2));
    // set up bounds and bounded by

    // iterate over all 0-cells / junctions
    for(size_t cell0Index=0;cell0Index<tgrid.numCells(0);++cell0Index){
        const LabelType cell0Label=cell0Index+1;
        // set up label
        Cells0_[cell0Index].setLabel(cell0Label);
        // get coordinates
        const size_t tx=Cells0_[cell0Index][0][0];
        const size_t ty=Cells0_[cell0Index][0][1];
        // Loop over all possible Cell1's / boundaries of the Cell0 / Junction
        const int px[]={ 1, -1, 0, 0};
        const int py[]={ 0,  0, 1,-1};
        for(size_t b=0;b<4;++b){
            LabelType cell1Label=grid(int(tx)+px[b],int(ty)+py[b]);
            // check if Cell1 / boundary is active
            if(cell1Label!=0){

                CGP_ASSERT_OP(cell1Label,>,0)
                CGP_ASSERT_OP(cell1Label,<=,tgrid.numCells(1));

                LabelType cell1Index=cell1Label-1;
                // bounds ( boundaries of a juction)
                Cells0_[cell0Index].push_back_bound(cell1Label);
                // junctions of a boundaty
                Cells1_[cell1Index].push_back_bounded_by(cell0Label);
            }
        }
        CGP_ASSERT_OP(Cells0_[cell0Index].bounds().size(),>=,3);
        CGP_ASSERT_OP(Cells0_[cell0Index].bounds().size(),<=,4);
    }

    // iterate over all 1-cells / boundaries
    for(size_t cell1Index=0;cell1Index<tgrid.numCells(1);++cell1Index){
        const LabelType cell1Label=cell1Index+1;
    
        // set up label
        Cells1_[cell1Index].setLabel(cell1Label);
        // get tx and ty of SOME element of the boundary (the first in this case) 
        const size_t tx=Cells1_[cell1Index][0][0];
        const size_t ty=Cells1_[cell1Index][0][1];
        // bounds (region labels)
        LabelType cell2LabelA,cell2LabelB;
        // vertical boundary
        if(tx%2==1){
            cell2LabelA=static_cast<LabelType>(grid(tx-1,ty));
            cell2LabelB=static_cast<LabelType>(grid(tx+1,ty));

        }
        else{
            cell2LabelA=static_cast<LabelType>(grid(tx,ty-1));
            cell2LabelB=static_cast<LabelType>(grid(tx,ty+1));
        }
        CGP_ASSERT_OP(cell2LabelA,>,0)
        CGP_ASSERT_OP(cell2LabelA,<=,tgrid.numCells(2));
        CGP_ASSERT_OP(cell2LabelB,>,0)
        CGP_ASSERT_OP(cell2LabelB,<=,tgrid.numCells(2));
        const LabelType cell2IndexA=cell2LabelA-1;
        const LabelType cell2IndexB=cell2LabelB-1;

        // set up bounds (the 2 adj. regions to this boundary)
        Cells1_[cell1Index].push_back_bound(cell2LabelA);
        Cells1_[cell1Index].push_back_bound(cell2LabelB);
        // set up bounded by ( n adj. boundaries of a region)
        Cells2_[cell2IndexA].push_back_bounded_by(cell1Label);
        Cells2_[cell2IndexB].push_back_bounded_by(cell1Label);

        Cells1_[cell1Index].sortCells();
    }
    // sortAdjaceny

    // iterate over all 2-cells / regions 
    for(size_t cell2Index=0;cell2Index<tgrid.numCells(2);++cell2Index){
        // set up ptr
        Cells2_[cell2Index].cgp_=this;
        // set up label
        Cells2_[cell2Index].setLabel(cell2Index+1);
        // sortAdjaceny
        Cells2_[cell2Index].sortAdjaceny();
    }
    // iterate over all 1-cells / boundaries
    for(size_t cell1Index=0;cell1Index<tgrid.numCells(1);++cell1Index){
        // set up ptr
        Cells1_[cell1Index].cgp_=this;
        // sortAdjaceny// sortAdjaceny
        Cells1_[cell1Index].setLabel(cell1Index+1);
        Cells1_[cell1Index].sortAdjaceny();
    }
    // iterate over all 0-cells / junctions
    for(size_t cell0Index=0;cell0Index<tgrid.numCells(0);++cell0Index){
        // set up ptr
        Cells0_[cell0Index].cgp_=this;
        // sortAdjaceny
        Cells0_[cell0Index].setLabel(cell0Index+1);
        // sortAdjaceny
        Cells0_[cell0Index].sortAdjaceny();
    }

}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::TopologicalGridType & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::tgrid()const{
    return tgrid_;
}
template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::Cells0  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry0()const{
    return Cells0_;
}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::Cells1  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry1()const{
    return Cells1_;
}

template<class COORDINATE_TYPE,class LABEL_TYPE>
const typename Cgp<COORDINATE_TYPE,LABEL_TYPE>::Cells2  & 
Cgp<COORDINATE_TYPE,LABEL_TYPE>::geometry2()const{
    return Cells2_;
}

} /* namespace vigra */

#endif /* CGP2D_HXX */