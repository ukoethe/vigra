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

    namespace detail_multi_blocking{
        template<unsigned int DIM, class C>
        class BlockWithBorder{
        public:
            typedef C  PointValue;
            typedef TinyVector<PointValue, DIM> Point;
            typedef Point Shape;
            typedef Box<PointValue, DIM> Block;
            typedef MultiCoordinateIterator< DIM> MultiCoordIter;

            BlockWithBorder(const Block & core, const Block & border)
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

        private:
            Block core_;
            Block border_;
        };

        
        
    }

    template<unsigned int DIM, class C>
    vigra::TinyVector<C, DIM> makePositivCoordinate(const vigra::TinyVector<C, DIM> & vec, const vigra::TinyVector<C, DIM> & shape){
        vigra::TinyVector<C, DIM> res(vigra::SkipInitialization);
        for(size_t d=0; d<DIM; ++d){
            res[d] = vec[d] < 0 ? vec[d] + shape[d] : vec[d];
        }
        return res;
    }


    template<class MB>
    class MultiCoordToBlockWithBoarder{
    public:
        typedef typename MB::Shape Shape;
        typedef typename MB::MultiCoordIter input_iter;
        typedef typename MB::BlockWithBorder result_type;
        MultiCoordToBlockWithBoarder(const MB & mb, const Shape & width)
        :   mb_(&mb),
            width_(width){
        }

        result_type operator()(const input_iter & iter)const{
            return mb_->getBlockWithBorder(*iter, width_);
        }
    private:
        const MB * mb_;
        Shape width_;
    };


    template<unsigned int DIM, class C = MultiArrayIndex>
    class MultiBlocking{
    public:
        typedef MultiBlocking<DIM, C> SelfType;
        typedef C  PointValue;
        typedef TinyVector<PointValue, DIM> Point;
        typedef Point Shape;
        typedef Point BlockDesc;
        typedef Box<PointValue, DIM> Block;
        typedef MultiCoordinateIterator< DIM> MultiCoordIter;
        typedef detail_multi_blocking::BlockWithBorder<DIM, PointValue> BlockWithBorder;


        // iterators 
        typedef MultiCoordToBlockWithBoarder<SelfType> CoordToBwb;
        typedef TransformIterator<CoordToBwb, MultiCoordIter> BlockWithBorderIter;



        MultiBlocking(const Shape & shape,
                      const Shape & blockShape,
                      const Shape & roiBegin = Shape(0),
                      const Shape & roiEnd = Shape(0)
        )
        :   shape_(shape),
            roiBlock_(makePositivCoordinate<DIM, C>(roiBegin, shape),
                      roiEnd == Shape(0) ? shape : makePositivCoordinate<DIM, C>(roiEnd, shape)),
            blockShape_(blockShape),
            blocksPerAxis_(vigra::SkipInitialization),
            blockDescIter_(),
            numBlocks_(1)
        {
            const Shape roiShape = roiBlock_.size();
            blocksPerAxis_ = roiShape / blockShape_;

            for(size_t d=0; d<DIM; ++d){
                if(blocksPerAxis_[d]*blockShape_[d] < roiShape[d] ){
                    ++blocksPerAxis_[d];
                }
                numBlocks_ *= blocksPerAxis_[d];
            }
            blockDescIter_ = MultiCoordIter(blocksPerAxis_);
        }

        /// total number of blocks
        size_t numBlocks()const{
            return numBlocks_;
        }

        BlockWithBorderIter BlockWithBorderBegin(const Shape & width)const{
            return BlockWithBorderIter(MultiCoordIter(blocksPerAxis_),
                                       CoordToBwb(*this, width));
        }


        BlockWithBorderIter BlockWithBorderEnd(const Shape & width)const{
            const MultiCoordIter beginIter(blocksPerAxis_);
            return BlockWithBorderIter(beginIter.getEndIterator(),CoordToBwb(*this, width));
        }

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


        /// get a block with border
        BlockWithBorder getBlockWithBorder(const size_t index, const Shape & width )const{
            return getBlockWithBorder(blockDescIter_[index], width);
        }


        /// get a block (without any overlap)
        Block getBlock(const BlockDesc & blockDesc)const{
            const Point blockStart(blockDesc * blockShape_ + roiBlock_.begin());
            const Point blockEnd(blockStart + blockShape_);
            return Block(blockStart, blockEnd) & roiBlock_;
        }

        /// get a block (without any overlap)
        Block getBlock(const size_t index)const{
            return getBlock(blockDescIter_[index]);
        }

        const Shape & roiBegin()const{
            return roiBlock_.begin();
        }

        const Shape & roiEnd()const{
            return roiBlock_.end();
        }

    private:
        Shape shape_;
        Block roiBlock_;
        Shape blockShape_;
        Shape blocksPerAxis_;
        MultiCoordIter blockDescIter_;
        size_t numBlocks_;
    };

}

#endif // VIGRA_MULTI_BLOCKING_HXX
