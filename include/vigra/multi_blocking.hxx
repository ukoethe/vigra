/************************************************************************/
/*                                                                      */
/*               Copyright 2015 by Thorsten Beier                       */
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

#ifndef VIGRA_MULTI_BLOCKING_HXX
#define VIGRA_MULTI_BLOCKING_HXX

#include "vigra/tinyvector.hxx"
#include "vigra/box.hxx"
#include "vigra/multi_iterator.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/transform_iterator.hxx"

namespace vigra{

    // forward declaration
    template<unsigned int DIM, class C>
    class MultiBlocking;

    /// \cond
    namespace detail_multi_blocking{

        template<unsigned int DIM, class C>
        class BlockWithBorder{
        public:
            typedef C  PointValue;
            typedef TinyVector<PointValue, DIM> Point;
            typedef Point Shape;
            typedef Box<PointValue, DIM> Block;
            typedef MultiCoordinateIterator< DIM> MultiCoordIter;

            BlockWithBorder(const Block & core = Block(), const Block & border = Block())
            :   core_(core),
                border_(border){
            }

            /// get the core block
            const Block & core()const{
                return core_;
            }

            Block  localCore()const{
                return core_-border_.begin();
            }


            const Block & border()const{
                return border_;
            }

            bool operator == (const BlockWithBorder & rhs) const{
                return core_ == rhs.core_ && border_ == rhs.border_;
            }
            bool operator != (const BlockWithBorder & rhs) const{
                return core_ != rhs.core_ || border_ != rhs.border_;
            }

        private:
            Block core_;
            Block border_;
        };

        template<class VALUETYPE, unsigned int DIMENSION>
        std::ostream& operator<< (std::ostream& stream, const BlockWithBorder<DIMENSION,VALUETYPE> & bwb) {
            stream<<"["<<bwb.core()<<", "<<bwb.border()<<" ]";
            return stream;
        }


        template<class MB>
        class MultiCoordToBlockWithBoarder{
        public:
            typedef typename MB::Shape Shape;
            typedef typename MB::BlockDesc BlockDesc;
            typedef typename MB::BlockWithBorder result_type;
            MultiCoordToBlockWithBoarder()
            :   mb_(NULL),
                width_(){
            }
            MultiCoordToBlockWithBoarder(const MB & mb, const Shape & width)
            :   mb_(&mb),
                width_(width){
            }

            result_type operator()(const BlockDesc & blockDesc)const{
                return mb_->getBlockWithBorder(blockDesc, width_);
            }
        private:
            const MB * mb_;
            Shape width_;
        };

        template<class MB>
        class MultiCoordToBlock{
        public:
            typedef typename MB::Shape Shape;
            typedef typename MB::BlockDesc BlockDesc;
            typedef typename MB::Block result_type;
            MultiCoordToBlock()
            :   mb_(NULL){
            }
            MultiCoordToBlock(const MB & mb)
            :   mb_(&mb){
            }

            result_type operator()(const BlockDesc & blockDesc)const{
                return mb_->getBlock(blockDesc);
            }
        private:
            const MB * mb_;
        };
    }
    /// \endcond
    

    /**
        MultiBlocking is used to split a image / volume / multiarray
        into non-overlapping blocks.
        These non overlapping blocks are called cores.
        A border can be added to the core boxes.
        These 'core+border' blocks are just called border.
        The core block within the coordinate system 
        of the border block is called local core.
    */
    template<unsigned int DIM, class C = MultiArrayIndex>
    class MultiBlocking{
    public:
        typedef MultiBlocking<DIM, C> SelfType;
        friend class detail_multi_blocking::MultiCoordToBlock<SelfType>;
        friend class detail_multi_blocking::MultiCoordToBlockWithBoarder<SelfType>;
        typedef C  PointValue;
        typedef TinyVector<PointValue, DIM> Point;
        typedef Point Shape;
        typedef Point BlockDesc;
        typedef Box<PointValue, DIM> Block;
        typedef MultiCoordinateIterator< DIM> MultiCoordIter;
        typedef detail_multi_blocking::BlockWithBorder<DIM, PointValue> BlockWithBorder;


        // iterators 
        typedef detail_multi_blocking::MultiCoordToBlockWithBoarder<SelfType> CoordToBwb;
        typedef detail_multi_blocking::MultiCoordToBlock<SelfType> CoordToB;
        typedef EndAwareTransformIterator<CoordToBwb, MultiCoordIter> BlockWithBorderIter;
        typedef EndAwareTransformIterator<CoordToB,   MultiCoordIter> BlockIter;


        MultiBlocking(const Shape & shape,
                      const Shape & blockShape,
                      const Shape & roiBegin = Shape(0),
                      const Shape & roiEnd = Shape(0)
        )
        :   shape_(shape),
            roiBlock_(roiBegin,roiEnd == Shape(0) ? shape : roiEnd),
            blockShape_(blockShape),
            blocksPerAxis_(vigra::SkipInitialization),
            numBlocks_(1)
        {
            const Shape roiShape = roiBlock_.size();
            blocksPerAxis_ = roiShape / blockShape_;


            // blocking
            for(size_t d=0; d<DIM; ++d){
                if(blocksPerAxis_[d]*blockShape_[d] < roiShape[d] ){
                    ++blocksPerAxis_[d];
                }
                numBlocks_ *= blocksPerAxis_[d];
            }

            // total image border blocks
            Shape beginCA(0),endCB(shape);
            for(size_t d=0; d<DIM; ++d){
                {
                    // fix coordinate d to zero
                    Shape endCA(shape);
                    endCA[d] = 1;
                    volumeBorderBlocks_.push_back(Block(beginCA,endCA));
                }
                {
                    // fix coordinate d to shape[dim]-1
                    Shape beginCB(shape);
                    beginCB[d] -= 1;
                    volumeBorderBlocks_.push_back(Block(beginCB,endCB));
                }
            }

            insideVolBlock_.setBegin(Shape(1));
            Shape insideVolBlockShapeEnd(shape);
            insideVolBlockShapeEnd -= Shape(1);
            insideVolBlock_.setEnd(insideVolBlockShapeEnd);
        }

        /// total number of blocks
        size_t numBlocks()const{
            return numBlocks_;
        }

        BlockWithBorderIter blockWithBorderBegin(const Shape & width)const{
            return BlockWithBorderIter(MultiCoordIter(blocksPerAxis_),
                                       CoordToBwb(*this, width));
        }

        BlockWithBorderIter blockWithBorderEnd(const Shape & width)const{
            const MultiCoordIter beginIter(blocksPerAxis_);
            return BlockWithBorderIter(beginIter.getEndIterator(),
                                       CoordToBwb(*this, width));
        }

        BlockIter blockBegin()const{
            return BlockIter(MultiCoordIter(blocksPerAxis_),CoordToB(*this));
        }


        BlockIter blockEnd()const{
            const MultiCoordIter beginIter(blocksPerAxis_);
            return BlockIter(beginIter.getEndIterator(),CoordToB(*this));
        }


        Block blockDescToBlock(const BlockDesc & blockDesc)const{
            MultiCoordIter iter(blocksPerAxis_);
            iter+=blockDesc;
            return *BlockIter(iter,CoordToB(*this));
        }


        /// does this block intersect with the volume border
        bool containsVolumeBorder(const Block & block) const {
            return !insideVolBlock_.contains(block);
        }


        const Shape & roiBegin()const{
            return roiBlock_.begin();
        }

        const Shape & roiEnd()const{
            return roiBlock_.end();
        }

        const Shape & shape()const{
            return shape_;
        }

        const Shape & blockShape()const{
            return blockShape_;
        }

        const Shape & blocksPerAxis()const{
            return blocksPerAxis_;
        }

        const std::vector<Block> & volumeBorderBlocks()const{
            return volumeBorderBlocks_;
        }


        std::vector<UInt32> intersectingBlocks(
            const Shape roiBegin,
            const Shape roiEnd
        )const{
            size_t i=0;
            std::vector<UInt32> iBlocks;
            const Block testBlock(roiBegin, roiEnd);
            for(BlockIter iter=blockBegin(); iter!=blockEnd(); ++iter){
                if(testBlock.intersects(*iter)){
                    iBlocks.push_back(i);
                }
                ++i;
            }
            return std::move(iBlocks);
        }



    private:

        /// get a block with border
        BlockWithBorder getBlockWithBorder(const BlockDesc & blockDesc, const Shape & width )const{
            const Point blockStart(blockDesc * blockShape_ + roiBlock_.begin());
            const Point blockEnd(blockStart + blockShape_);
            const Block core = Block(blockStart, blockEnd) & roiBlock_ ;
            Block border = core;
            border.addBorder(width);
            border &= Block(shape_);
            return BlockWithBorder( core, border );
        }

        /// get a block (without any overlap)
        Block getBlock(const BlockDesc & blockDesc)const{
            const Point blockStart(blockDesc * blockShape_ + roiBlock_.begin());
            const Point blockEnd(blockStart + blockShape_);
            return Block(blockStart, blockEnd) & roiBlock_;
        }

        Shape shape_;           // total shape of the input volume
        Block roiBlock_;        // ROI in which to compute filters/algorithms
        Shape blockShape_;      // shape sub-block for each thread (without border pixels)
        Shape blocksPerAxis_;   // how many blocks are on each axis
        size_t numBlocks_;      // total number of blocks


        std::vector<Block> volumeBorderBlocks_;
        Block insideVolBlock_;
    };

}

#endif // VIGRA_MULTI_BLOCKING_HXX
