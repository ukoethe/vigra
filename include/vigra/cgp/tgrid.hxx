#ifndef CGP2D_TGRID_HXX
#define CGP2D_TGRID_HXX

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

/* this project */
#include "macros.hxx"



namespace vigra{


/****************************************************************************/
/* T o p o l o g i c a l G r i d                                            */
/****************************************************************************/

template<class LABEL_TYPE>
class TopologicalGrid{
public:
    typedef partition::Partition<size_t>    UfdType;
    typedef LABEL_TYPE LabelType;
    typedef vigra::MultiArray<2,LabelType>  LabelImageType;
    typedef typename LabelImageType::difference_type ShapeType;

    TopologicalGrid(){
    
    }

    /*
    TopologicalGrid(const TopologicalGrid<LABEL_TYPE> & labelType){
        std::cout<<"bla bla 2\n";
    }

    template<class L>
    TopologicalGrid(const TopologicalGrid<L> & labelType){
        std::cout<<"bla bla 1 \n";
    }
    */



    // constructor 
    template<class INPUT_IMG>
    TopologicalGrid( const INPUT_IMG & seg);

    // query
    const LabelImageType & tgrid()const;
    size_t numCells(const size_t i)const;
    size_t shape(const size_t d)const;

    const ShapeType & shapeTopologicalGrid()const{
        return tShape_;
    }
    const ShapeType & shapeLabeling()const{
        return lShape_;
    }

    LabelType operator()(const size_t tx,const size_t ty)     {return tgrid_(tx,ty);}
    LabelType operator()(const size_t tx,const size_t ty)const{return tgrid_(tx,ty);}

private:

    size_t numCells_[3];
    LabelImageType tgrid_;
    ShapeType tShape_;
    ShapeType lShape_;

};


// Implementation ToopogicalGrid
template<class LABEL_TYPE>
template<class INPUT_IMG>
TopologicalGrid<LABEL_TYPE>::TopologicalGrid(const INPUT_IMG & seg)
: tgrid_(ShapeType(seg.shape(0)*2-1,seg.shape(1)*2-1)),
    tShape_(seg.shape(0)*2-1,seg.shape(1)*2-1),
    lShape_(seg.shape(0),seg.shape(1))
{
    const size_t dx=seg.shape(0);
    const size_t dy=seg.shape(1);
    const size_t tdx=seg.shape(0)*2-1;
    const size_t tdy=seg.shape(1)*2-1;

    //std::cout<<"dx  "<<dx<<"   dy  "<<dy<<" \n";
    //std::cout<<"tdx "<<tdx<<" tdy  "<<tdy<<" \n";


    size_t shape[] = { tdx,tdy};
    // counters
    size_t maxRegionLabel=0; //TODO TODO TODO
    size_t junctionIndex=0;
    size_t boundaryElementLabel=1;
    ////////////////
    // 1. PASS //
    ////////////////
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        //std::cout<<" tx "<<tx<<" ty "<<ty<<"\n";
        // if region
        if(tx%2==0 && ty%2==0){

            //std::cout<<"tx "<<tx<<"  ty "<<ty<<"  tx/2 "<<tx/2<<" ty/2 "<<ty/2<<" \n";
            //std::cout<< seg(tx/2,ty/2) <<"\n";
            size_t label=seg(tx/2,ty/2);
            CGP_ASSERT_OP(label,>,0);
            tgrid_(tx,ty)=label;
            maxRegionLabel=label>maxRegionLabel ? label : maxRegionLabel;
        }
        // if junction
        else if(tx%2!=0 && ty%2!=0){
            //  A|B
            //  _ _
            //  C|D
            // check labels of A,B,C and D
            std::set<LABEL_TYPE> lset;
            lset.insert( seg((tx-1)/2,(ty-1)/2));  // A
            lset.insert( seg((tx+1)/2,(ty-1)/2));  // B
            lset.insert( seg((tx-1)/2,(ty+1)/2));  // A
            lset.insert( seg((tx+1)/2,(ty+1)/2));  // A
            if(lset.size()>=3){
                tgrid_(tx,ty)=junctionIndex+1;
                ++junctionIndex;
            }
            else{
                tgrid_(tx,ty)=0;
            }
        }
        // boundary
        else{
            size_t l0,l1;
            // A|B
            // vertical  boundary 
            if(tx%2==1){
                l0=seg( (tx-1)/2, ty/2 );
                l1=seg( (tx+1)/2, ty/2 );
            }
            // horizontal boundary
            else{
                l0=seg( tx/2, (ty-1)/2);
                l1=seg( tx/2, (ty+1)/2);
            }
            // active boundary ?
            if(l0!=l1){
                //std::cout<<l0<<"/"<<l1<<"\n";
                tgrid_(tx,ty)=boundaryElementLabel;
                ++boundaryElementLabel;
            }
            else
                tgrid_(tx,ty)=0;
        }
    }
    /////////////////
    // 2. PASS //
    /////////////////
    UfdType boundaryUdf(boundaryElementLabel-1);
    const size_t num1Elements=boundaryElementLabel-1;
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        // boundary
        if((tx%2!=0 && ty%2!=1 ) || (tx%2!=1 && ty%2!=0 )) {
            if ( tgrid_(tx,ty)!=0){
                size_t ownIndex=tgrid_(tx,ty);
                // vertical boundary
                if(tx%2==1){
                    // each horizontal boundary has 6 candidate neighbours:
                    //  _|_
                    //  _|_  <= this is the boundary we are looking right now
                    //   |
                    // if junction down is inctive
                    LabelType other;
                    if (ty+1 < tdy  && tgrid_(tx,ty+1)==0){
                        // boundary up is active?
                        other=tgrid_(tx,ty+2);
                        if( other!=0){
                            CGP_ASSERT_OP(other-1, <, num1Elements);
                            boundaryUdf.merge(ownIndex-1,other-1);
                        }
                        // boundary left is active?
                        other=tgrid_(tx-1,ty+1);
                        if( other!=0){
                            CGP_ASSERT_OP(other-1 ,<, num1Elements);
                            boundaryUdf.merge(ownIndex-1,other-1 );
                        }
                        // boundary right is active?
                        if( tgrid_(tx+1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty+1)-1 );
                        }
                    }
                    // if junction up is inctive
                    if(ty > 0 && tgrid_(tx,ty-1)==0){
                        // boundary up is active?
                        if( tgrid_(tx,ty-2)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx,ty-2) -1);
                        }
                        // boundary left is active?
                        if( tgrid_(tx-1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty-1)-1 );
                        }
                        // boundary right is active?
                        if( tgrid_(tx+1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty-1)-1 );
                        }
                    }
                }
                // horizontal boundary 
                else{
                    //   each horizontal boundary has 6 candidate neighbours:
                    //   _|_|_
                    //    | |
                    // 
                    // if left junction inactive?     
                    if(tx >0 && tgrid_( tx-1,ty)==0){
                        // boundary left is active?
                        if( tgrid_(tx-2,ty)!=0){
                            //std::cout<<"merge left \n";
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-2,ty)-1 );
                        }
                        // boundary up is active?
                        if( tgrid_(tx-1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty-1)-1 );
                        }
                        // boundary down is active?
                        if( tgrid_(tx-1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx-1,ty+1)-1 );
                        }
                    }
                    // if right junction inactive?     
                    if(tx+1<tdx &&tgrid_( tx+1,ty)==0){
                        // boundary right is active?
                        if( tgrid_(tx+2,ty)!=0){
                            //std::cout<<"merge right \n";
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+2,ty)-1 );
                        }
                        // boundary up is active?
                        if( tgrid_(tx+1,ty-1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty-1)-1 );
                        }
                        // boundary down is active?
                        if( tgrid_(tx+1,ty+1)!=0){
                            boundaryUdf.merge(ownIndex-1,tgrid_(tx+1,ty+1)-1 );
                        }
                    }
                }
            }
        }
    }

    // dense relabeling
    std::map<size_t,size_t> relabel;
    boundaryUdf.representativeLabeling(relabel);
    /////////////////
    // 3. PASS //
    /////////////////
    for(size_t ty=0;ty<tdy;++ty)
    for(size_t tx=0;tx<tdx;++tx){
        // boundary
        if((tx%2!=0 && ty%2!=1 ) || (tx%2!=1 && ty%2!=0 )) {
            if(tgrid_(tx,ty)!=0){
                // relabel
                size_t notDenseIndex=boundaryUdf.find( tgrid_(tx,ty)-1 );
                size_t denseIndex=relabel[notDenseIndex];
                tgrid_(tx,ty)=denseIndex+1;
            }
        }
    }

    // update cell counters
    numCells_[2]=maxRegionLabel;
    CGP_ASSERT_OP(boundaryUdf.numberOfSets(),==,relabel.size());
    numCells_[1]=relabel.size();
    numCells_[0]=junctionIndex;
}

template<class LABEL_TYPE>
const typename TopologicalGrid<LABEL_TYPE>::LabelImageType & TopologicalGrid<LABEL_TYPE>::tgrid()const{
    return tgrid_;
}

template<class LABEL_TYPE>
size_t TopologicalGrid<LABEL_TYPE>::numCells(const size_t i)const{
    return numCells_[i];
}

template<class LABEL_TYPE>
size_t TopologicalGrid<LABEL_TYPE>::shape(const size_t d)const{
    return tgrid_.shape(d);
}


}

#endif //CGP2D_TGRID_HXX