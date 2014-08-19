#ifndef SKNEURO_UTILITIES_BLOCKING 
#define SKNEURO_UTILITIES_BLOCKING 

#include <vigra/box.hxx>
#include <iostream>


namespace vigra{

template<class COORDINATE, int DIMENSION>
class Block
: public Box<COORDINATE, DIMENSION>{

public:
    typedef Box<COORDINATE, DIMENSION> BaseType;
    typedef typename BaseType::Vector CoordType;


    Block()
    :   BaseType(){
        
    }

    Block(const CoordType & a, const CoordType & b)
    :   BaseType(a,b){

    }

    std::string str()const{
        std::stringstream ss;
        ss<<"(";
        for(int i=0; i<DIMENSION; ++i){
            ss<<this->begin()[i]<<" ";
        }
        ss<<") -- (";
        for(int i=0; i<DIMENSION; ++i){
            ss<<this->end()[i]<<" ";
        }
        ss<<" )";
        return ss.str();
    }

private:
};

template<class COORDINATE, int DIMENSION>
class BlockWithBorder{
public:
    typedef Block<COORDINATE, DIMENSION>  BlockType;
    typedef typename BlockType::CoordType CoordType;

    BlockWithBorder(){  
    }

    BlockWithBorder(const BlockType & core, const BlockType & coreWithBorder)
    :   block_(core),
        blockWithAddedBorder_(coreWithBorder){
    }

    std::string str()const{
        std::stringstream ss;
        ss<<"Core "<<block_<<" CoreWithBorder "<<blockWithAddedBorder_;
        return ss.str();
    }
    
    // core (coordinates with respect to the coreWithBorder)


    CoordType blockShape()const{
        return block_.size();
    }

    CoordType blockWithBorderShape()const{
        return blockWithAddedBorder_.size();
    }



    const BlockType & block()const{
        return block_;
    }

    const BlockType & blockWithBorder()const{
        return blockWithAddedBorder_;
    }

    BlockType blockLocalCoordinates()const{
        BlockType res=block_;
        res-=blockWithAddedBorder_.begin();
        return res;
    }


    BlockType blockWithBoarderLocalCoordinates()const{
        return blockWithAddedBorder_ - block_.begin();
    }


private:
    BlockType block_;
    BlockType blockWithAddedBorder_;
};


template<class COORDINATE,int DIMENSION>
std::ostream & operator<<(std::ostream & lhs, const Block<COORDINATE,DIMENSION > & block ){
    lhs<<block.str();
    return lhs;
}

template<class COORDINATE,int DIMENSION>
std::ostream & operator<<(std::ostream & lhs, const BlockWithBorder<COORDINATE,DIMENSION > & blockWithBorder ){
    lhs<<blockWithBorder.str();
    return lhs;
}


template<class COORDINATE, int DIMENSION>
class Blocking{

public:
    
    typedef Block<COORDINATE, DIMENSION>  BlockType;
    typedef BlockWithBorder<COORDINATE, DIMENSION>  BlockWithBorderType;
    typedef typename BlockType::CoordType CoordType;

    Blocking()
    :   shape_(),
        blockShape_(),
        totalBlock_(),
        blocking_(){

    }

    Blocking(const CoordType & shape, const CoordType & blockShape)
    :   shape_(shape),
        blockShape_(blockShape),
        totalBlock_(CoordType(0),shape),
        blocking_(){

        SKNEURO_CHECK_OP(DIMENSION,==,3,"currently only implemented for 3D");

        CoordType blockStart(0);
        //for(blockStart[0]=0; blockStart[0]<shape[0]; blockStart[0]+=blockShape[0])
        //for(blockStart[1]=0; blockStart[1]<shape[1]; blockStart[1]+=blockShape[1])
        //for(blockStart[2]=0; blockStart[2]<shape[2]; blockStart[2]+=blockShape[2]){

        for(blockStart[2]=0; blockStart[2]<shape[2]; blockStart[2]+=blockShape[2])
        for(blockStart[1]=0; blockStart[1]<shape[1]; blockStart[1]+=blockShape[1])
        for(blockStart[0]=0; blockStart[0]<shape[0]; blockStart[0]+=blockShape[0]){

            CoordType blockEnd = blockStart + blockShape;
            BlockType block(blockStart,blockEnd);
            // intersect
            block &= totalBlock_;
            blocking_.push_back(block);
        }   
    }

    size_t size()const{
        return blocking_.size();
    }

    const BlockType & operator[](const size_t i)const{
        SKNEURO_ASSERT_OP(i,<,blocking_.size());
        return blocking_[i];
    }

    BlockWithBorderType blockWithBorder(const size_t index , const size_t width)const{
        BlockType blockWithBorder = blocking_[index];
        blockWithBorder.addBorder(width);
        // intersect
        blockWithBorder &=  totalBlock_;
        const BlockWithBorderType bwb( blocking_[index], blockWithBorder);
        return bwb;
    }   

private:
    CoordType shape_;
    CoordType blockShape_;
    BlockType totalBlock_;
    std::vector<BlockType> blocking_;
};




template<class COORDINATE,class T_BLOCK, class T_TOTAL>
void extractBlock(
    const BlockWithBorder<COORDINATE,3> & blockWithBorder,
    const MultiArrayView<3,T_TOTAL> & totalData,
    MultiArrayView<3,T_BLOCK> & blockData
){
    typedef typename BlockWithBorder<COORDINATE,3>::CoordType CoordType;
    CoordType begin = blockWithBorder.blockWithBorder().begin();
    CoordType end   = blockWithBorder.blockWithBorder().end();


    CoordType totalCoord;
    CoordType blockCoord;
    //for(totalCoord[0]=begin[0],blockCoord[0]=0; totalCoord[0]<end[0]; ++totalCoord[0],++blockCoord[0])
    //for(totalCoord[1]=begin[1],blockCoord[1]=0; totalCoord[1]<end[1]; ++totalCoord[1],++blockCoord[1])
    //for(totalCoord[2]=begin[2],blockCoord[2]=0; totalCoord[2]<end[2]; ++totalCoord[2],++blockCoord[2]){

    for(totalCoord[2]=begin[2],blockCoord[2]=0; totalCoord[2]<end[2]; ++totalCoord[2],++blockCoord[2])
    for(totalCoord[1]=begin[1],blockCoord[1]=0; totalCoord[1]<end[1]; ++totalCoord[1],++blockCoord[1])
    for(totalCoord[0]=begin[0],blockCoord[0]=0; totalCoord[0]<end[0]; ++totalCoord[0],++blockCoord[0]){
        blockData[blockCoord]=totalData[totalCoord];
    }
}


template<class COORDINATE,class T_BLOCK, class T_TOTAL>
void writeFromBlock(
    const BlockWithBorder<COORDINATE,3> & blockWithBorder,
    const MultiArrayView<3,T_BLOCK> & blockData,
    MultiArrayView<3,T_TOTAL> & totalData
){
    typedef typename BlockWithBorder<COORDINATE,3>::CoordType CoordType;
    const CoordType tBegin  = blockWithBorder.block().begin();
    const CoordType tEnd = blockWithBorder.block().end();

    const CoordType bBegin = blockWithBorder.blockLocalCoordinates().begin();
    //const CoordType bEnd = blockWithBorder.blockLocalCoordinates().end();


    CoordType tCoord;
    CoordType bCoord;
    for(tCoord[2]=tBegin[2],bCoord[2]=bBegin[2]; tCoord[2]<tEnd[2]; ++tCoord[2],++bCoord[2])
    for(tCoord[1]=tBegin[1],bCoord[1]=bBegin[1]; tCoord[1]<tEnd[1]; ++tCoord[1],++bCoord[1])
    for(tCoord[0]=tBegin[0],bCoord[0]=bBegin[0]; tCoord[0]<tEnd[0]; ++tCoord[0],++bCoord[0]){
        totalData[tCoord]=blockData[bCoord];
    }
}


} // end namespace vigra

#endif /*SKNEURO_UTILITIES_BLOCKING */