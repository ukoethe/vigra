#ifndef VIGRA_MULTI_BLOCKING_HXX
#define VIGRA_MULTI_BLOCKING_HXX

#include "vigra/tinyvector.hxx"
#include "vigra/box.hxx"
#include "vigra/multi_iterator.hxx"
#include "vigra/multi_convolution.hxx"

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


        template<unsigned int DIM, class C>
        class MultiBlockWithBorderIter
        : public MultiCoordinateIterator<DIM>
        {
            public:
                typedef C  PointValue;
                typedef TinyVector<PointValue, DIM> Point;
                typedef Point Shape;
                typedef Box<PointValue, DIM> Block;
                typedef MultiCoordinateIterator< DIM> MultiCoordIter;
                //typedef BlockWithBorder<DIM, PointValue> BlockWithBorder;
                typedef vigra::MultiBlocking<DIM, PointValue> MultiBlockingType;

                MultiBlockWithBorderIter(const MultiBlockingType & multiBlocking, const MultiCoordIter & baseIter)
                :   MultiCoordIter(baseIter),
                    multiBlocking_(&multiBlocking_)
                {

                }

            private:
                const MultiBlockingType * multiBlocking_;
        };
    }


    template<unsigned int DIM, class C = MultiArrayIndex>
    class MultiBlocking{
    public:
        typedef C  PointValue;
        typedef TinyVector<PointValue, DIM> Point;
        typedef Point Shape;
        typedef Point BlockDesc;
        typedef Box<PointValue, DIM> Block;
        typedef MultiCoordinateIterator< DIM> MultiCoordIter;
        typedef detail_multi_blocking::BlockWithBorder<DIM, PointValue> BlockWithBorder;

        MultiBlocking(const Shape & shape,
                      const Shape & blockShape)
        :   shape_(shape),
            blockShape_(blockShape),
            blocksPerAxis_(shape/blockShape),
            blockDescIter_(),
            numBlocks_(1)
        {

            for(size_t d=0; d<DIM; ++d){
                if(blocksPerAxis_[d]*blockShape_[d] < shape_[d] ){
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



        /// get a block with border
        BlockWithBorder getBlockWithBorder(const BlockDesc & blockDesc, const Shape & width )const{
            const Point blockStart(blockDesc * blockShape_);
            const Point blockEnd(blockStart + blockShape_);
            const Block core = Block(blockStart, blockEnd) & Block(shape_) ;
            Block border = core;
            border.addBorder(width);
            border &=  Block(shape_);
            return BlockWithBorder( core, border );
        }


        /// get a block with border
        BlockWithBorder getBlockWithBorder(const size_t index, const Shape & width )const{
            return getBlockWithBorder(blockDescIter_[index], width);
        }


        /// get a block (without any overlap)
        Block getBlock(const BlockDesc & blockDesc)const{
            const Point blockStart(blockDesc * blockShape_);
            const Point blockEnd(blockStart + blockShape_);
            return Block(blockStart, blockEnd) & Block(shape_);
        }

        /// get a block (without any overlap)
        Block getBlock(const size_t index)const{
            return getBlock(blockDescIter_[index]);
        }


    private:
        Shape shape_;
        Shape blockShape_;
        Shape blocksPerAxis_;
        MultiCoordIter blockDescIter_;
        size_t numBlocks_;
    };

}

#endif // VIGRA_MULTI_BLOCKING_HXX
