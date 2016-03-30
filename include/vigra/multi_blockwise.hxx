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


#ifndef VIGRA_MULTI_BLOCKWISE_HXX
#define VIGRA_MULTI_BLOCKWISE_HXX

#include <cmath>
#include <vector>
#include "multi_blocking.hxx"
#include "multi_convolution.hxx"
#include "multi_tensorutilities.hxx"
#include "threadpool.hxx"
#include "array_vector.hxx"

namespace vigra{

    /** Option base class for blockwise algorithms.

        Attaches blockshape to ParallelOptions.
    */
class BlockwiseOptions
: public ParallelOptions
{
public:
    typedef ArrayVector<MultiArrayIndex> Shape;

    BlockwiseOptions()
    :   ParallelOptions()
    ,   blockShape_()
    {}

        /** Retrieve block shape as a std::vector.

            If the returned vector is empty, a default block shape should be used.
            If the returned vector has length 1, square blocks of size
            <tt>getBlockShape()[0]</tt> should be used.
        */
    Shape const & getBlockShape() const
    {
        return blockShape_;
    }

        // for Python bindings
    Shape readBlockShape() const
    {
        return blockShape_;
    }

        /** Retrieve block shape as a fixed-size vector.

            Default shape specifications are appropriately expanded.
            An exception is raised if the stored shape's length is
            incompatible with dimension <tt>N</tt>.
        */
    template <int N>
    TinyVector<MultiArrayIndex, N> getBlockShapeN() const
    {
        if(blockShape_.size() > 1)
        {
            vigra_precondition(blockShape_.size() == (size_t)N,
                "BlockwiseOptions::getBlockShapeN(): dimension mismatch between N and stored block shape.");
            return TinyVector<MultiArrayIndex, N>(blockShape_.data());
        }
        else if(blockShape_.size() == 1)
        {
            return TinyVector<MultiArrayIndex, N>(blockShape_[0]);
        }
        else
        {
            return detail::ChunkShape<N>::defaultShape();
        }
    }

        /** Specify block shape as a std::vector of appropriate length.
            If <tt>blockShape.size() == 0</tt>, the default shape is used.
            If <tt>blockShape.size() == 1</tt>, square blocks of size
            <tt>blockShape[0]</tt> are used.

            Default: Use square blocks with side length <tt>VIGRA_DEFAULT_BLOCK_SHAPE</tt>.
        */
    BlockwiseOptions & blockShape(const Shape & blockShape){
        blockShape_ = blockShape;
        return *this;
    }

        // for Python bindings
    void setBlockShape(const Shape & blockShape){
        blockShape_ = blockShape;
    }

        /** Specify block shape by a fixed-size shape object.
        */
    template <class T, int N>
    BlockwiseOptions & blockShape(const TinyVector<T, N> & blockShape){
        Shape(blockShape.begin(), blockShape.end()).swap(blockShape_);
        return *this;
    }

        /** Specify square block shape by its side length.
        */
    BlockwiseOptions & blockShape(MultiArrayIndex blockShape){
        Shape(1, blockShape).swap(blockShape_);
        return *this;
    }

    BlockwiseOptions & numThreads(const int n)
    {
        ParallelOptions::numThreads(n);
        return *this;
    }

    void setNumThreads(const int n)
    {
        ParallelOptions::numThreads(n);
    }

private:
    Shape blockShape_;
};

    /** Option class for blockwise convolution algorithms.

        Simply derives from \ref vigra::BlockwiseOptions and
        \ref vigra::ConvolutionOptions to join their capabilities.
    */
template<unsigned int N>
class BlockwiseConvolutionOptions
:   public  BlockwiseOptions
,   public  ConvolutionOptions<N>{
public:
    BlockwiseConvolutionOptions()
    :   BlockwiseOptions(),
        ConvolutionOptions<N>()
    {}
};


namespace blockwise{

    /**
        helper function to create blockwise parallel filters.
        This implementation should be used if the filter functor
        does not support the ROI/sub array options.
    */
    template<
        unsigned int DIM,
        class T_IN, class ST_IN,
        class T_OUT, class ST_OUT,
        class FILTER_FUNCTOR,
        class C
    >
    void blockwiseCallerNoRoiApi(
        const vigra::MultiArrayView<DIM, T_IN,  ST_IN > & source,
        const vigra::MultiArrayView<DIM, T_OUT, ST_OUT> & dest,
        FILTER_FUNCTOR & functor,
        const vigra::MultiBlocking<DIM, C> & blocking,
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth,
        const BlockwiseConvolutionOptions<DIM>  & options
    ){

        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;

        auto beginIter  =  blocking.blockWithBorderBegin(borderWidth);
        auto endIter   =  blocking.blockWithBorderEnd(borderWidth);

        parallel_foreach(options.getNumThreads(),
            beginIter, endIter,
            [&](const int /*threadId*/, const BlockWithBorder bwb)
            {
                // get the input of the block as a view
                vigra::MultiArrayView<DIM, T_IN, ST_IN> sourceSub = source.subarray(bwb.border().begin(),
                                                                             bwb.border().end());
                // get the output as NEW allocated array
                vigra::MultiArray<DIM, T_OUT> destSub(sourceSub.shape());
                // call the functor
                functor(sourceSub, destSub);
                 // write the core global out
                vigra::MultiArrayView<DIM, T_OUT, ST_OUT> destSubCore = destSub.subarray(bwb.localCore().begin(),
                                                                                bwb.localCore().end());
                // write the core global out
                dest.subarray(bwb.core().begin()-blocking.roiBegin(),
                              bwb.core().end()  -blocking.roiBegin()  ) = destSubCore;
            },
            blocking.numBlocks()
        );

    }

    /**
        helper function to create blockwise parallel filters.
        This implementation should be used if the filter functor
        does support the ROI/sub array options.
    */
    template<
        unsigned int DIM,
        class T_IN, class ST_IN,
        class T_OUT, class ST_OUT,
        class FILTER_FUNCTOR,
        class C
    >
    void blockwiseCaller(
        const vigra::MultiArrayView<DIM, T_IN,  ST_IN > & source,
        const vigra::MultiArrayView<DIM, T_OUT, ST_OUT> & dest,
        FILTER_FUNCTOR & functor,
        const vigra::MultiBlocking<DIM, C> & blocking,
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth,
        const BlockwiseConvolutionOptions<DIM>  & options
    ){

        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;
        //typedef typename MultiBlocking<DIM, C>::BlockWithBorderIter BlockWithBorderIter;
        typedef typename MultiBlocking<DIM, C>::Block Block;


        auto beginIter  =  blocking.blockWithBorderBegin(borderWidth);
        auto endIter   =  blocking.blockWithBorderEnd(borderWidth);

        parallel_foreach(options.getNumThreads(),
            beginIter, endIter,
            [&](const int /*threadId*/, const BlockWithBorder bwb)
            {
                // get the input of the block as a view
                vigra::MultiArrayView<DIM, T_IN, ST_IN> sourceSub = source.subarray(bwb.border().begin(),
                                                                            bwb.border().end());
                // get the output of the blocks core as a view
                vigra::MultiArrayView<DIM, T_OUT, ST_OUT> destCore = dest.subarray(bwb.core().begin(),
                                                                            bwb.core().end());
                const Block localCore =  bwb.localCore();
                // call the functor
                functor(sourceSub, destCore, localCore.begin(), localCore.end());
            },
            blocking.numBlocks()
        );


    }

    #define CONVOLUTION_FUNCTOR(FUNCTOR_NAME, FUNCTION_NAME) \
    template<unsigned int DIM> \
    class FUNCTOR_NAME{ \
    public: \
        typedef ConvolutionOptions<DIM> ConvOpt; \
        FUNCTOR_NAME(const ConvOpt & convOpt) \
        : sharedOpt_(convOpt){} \
        template<class S, class D> \
        void operator()(const S & s, D & d)const{ \
            FUNCTION_NAME(s, d, sharedOpt_); \
        } \
        template<class S, class D,class SHAPE> \
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){ \
            ConvOpt localOpt(sharedOpt_); \
            localOpt.subarray(roiBegin, roiEnd); \
            FUNCTION_NAME(s, d, localOpt); \
        } \
    private: \
        ConvOpt  sharedOpt_; \
    };


    CONVOLUTION_FUNCTOR(GaussianSmoothFunctor,            vigra::gaussianSmoothMultiArray);
    CONVOLUTION_FUNCTOR(GaussianGradientFunctor,          vigra::gaussianGradientMultiArray);
    CONVOLUTION_FUNCTOR(SymmetricGradientFunctor,         vigra::symmetricGradientMultiArray);
    CONVOLUTION_FUNCTOR(GaussianDivergenceFunctor,        vigra::gaussianDivergenceMultiArray);
    CONVOLUTION_FUNCTOR(HessianOfGaussianFunctor,         vigra::hessianOfGaussianMultiArray);
    CONVOLUTION_FUNCTOR(LaplacianOfGaussianFunctor,       vigra::laplacianOfGaussianMultiArray);
    CONVOLUTION_FUNCTOR(GaussianGradientMagnitudeFunctor, vigra::gaussianGradientMagnitude);
    CONVOLUTION_FUNCTOR(StructureTensorFunctor,           vigra::structureTensorMultiArray);

    #undef CONVOLUTION_FUNCTOR

    template<unsigned int DIM>
    class HessianOfGaussianEigenvaluesFunctor{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;
        HessianOfGaussianEigenvaluesFunctor(const ConvOpt & convOpt)
        : sharedOpt_(convOpt){}
        template<class S, class D>
        void operator()(const S & s, D & d)const{
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(d.shape());
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, sharedOpt_);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, d);
        }
        template<class S, class D,class SHAPE>
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(roiEnd-roiBegin);
            ConvOpt localOpt(sharedOpt_);
            localOpt.subarray(roiBegin, roiEnd);
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, localOpt);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, d);
        }
    private:
        ConvOpt  sharedOpt_;
    };

    template<unsigned int DIM, unsigned int EV>
    class HessianOfGaussianSelectedEigenvalueFunctor{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;
        HessianOfGaussianSelectedEigenvalueFunctor(const ConvOpt & convOpt)
        : sharedOpt_(convOpt){}
        template<class S, class D>
        void operator()(const S & s, D & d)const{
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;

            // compute the hessian of gaussian and extract eigenvalue
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(s.shape());
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, sharedOpt_);

            vigra::MultiArray<DIM, TinyVector<RealType, DIM > >  allEigenvalues(s.shape());
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, allEigenvalues);

            d = allEigenvalues.bindElementChannel(EV);
        }
        template<class S, class D,class SHAPE>
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){

            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;

            // compute the hessian of gaussian and extract eigenvalue
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(roiEnd-roiBegin);
            ConvOpt localOpt(sharedOpt_);
            localOpt.subarray(roiBegin, roiEnd);
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, localOpt);

            vigra::MultiArray<DIM, TinyVector<RealType, DIM > >  allEigenvalues(roiEnd-roiBegin);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, allEigenvalues);

            d = allEigenvalues.bindElementChannel(EV);
        }
    private:
        ConvOpt  sharedOpt_;
    };


    template<unsigned int DIM>
    class HessianOfGaussianFirstEigenvalueFunctor
    : public HessianOfGaussianSelectedEigenvalueFunctor<DIM, 0>{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;
        HessianOfGaussianFirstEigenvalueFunctor(const ConvOpt & convOpt)
        : HessianOfGaussianSelectedEigenvalueFunctor<DIM, 0>(convOpt){}
    };

    template<unsigned int DIM>
    class HessianOfGaussianLastEigenvalueFunctor
    : public HessianOfGaussianSelectedEigenvalueFunctor<DIM, DIM-1>{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;
        HessianOfGaussianLastEigenvalueFunctor(const ConvOpt & convOpt)
        : HessianOfGaussianSelectedEigenvalueFunctor<DIM, DIM-1>(convOpt){}
    };






    /// \warning this functions is deprecated
    /// and should not be used from end users
    template<unsigned int N>
    vigra::TinyVector< vigra::MultiArrayIndex, N > getBorder(
        const BlockwiseConvolutionOptions<N> & opt,
        const size_t order,
        const bool usesOuterScale = false
    ){
        vigra::TinyVector< vigra::MultiArrayIndex, N > res(vigra::SkipInitialization);

        if(opt.getFilterWindowSize()<=0.00001){
            for(size_t d=0; d<N; ++d){
                double stdDev =  opt.getStdDev()[d];
                if(usesOuterScale)
                    stdDev += opt.getOuterScale()[d];
                res[d] = static_cast<MultiArrayIndex>(3.0 * stdDev  + 0.5*static_cast<double>(order)+0.5);
            }
        }
        else{
            throw std::runtime_error("blockwise filters do not allow a user defined FilterWindowSize");
        }
        return res;
    }

} // end namespace blockwise

#define VIGRA_BLOCKWISE(FUNCTOR, FUNCTION, ORDER, USES_OUTER_SCALE) \
template <unsigned int N, class T1, class S1, class T2, class S2> \
void FUNCTION( \
    MultiArrayView<N, T1, S1> const & source, \
    MultiArrayView<N, T2, S2> dest, \
    BlockwiseConvolutionOptions<N> const & options \
) \
{  \
    typedef  MultiBlocking<N, vigra::MultiArrayIndex> Blocking; \
    typedef typename Blocking::Shape Shape; \
    const Shape border = blockwise::getBorder(options, ORDER, USES_OUTER_SCALE); \
    BlockwiseConvolutionOptions<N> subOptions(options); \
    subOptions.subarray(Shape(0), Shape(0));  \
    const Blocking blocking(source.shape(), options.template getBlockShapeN<N>()); \
    blockwise::FUNCTOR<N> f(subOptions); \
    blockwise::blockwiseCaller(source, dest, f, blocking, border, options); \
}

VIGRA_BLOCKWISE(GaussianSmoothFunctor,                   gaussianSmoothMultiArray,                   0, false );
VIGRA_BLOCKWISE(GaussianGradientFunctor,                 gaussianGradientMultiArray,                 1, false );
VIGRA_BLOCKWISE(SymmetricGradientFunctor,                symmetricGradientMultiArray,                1, false );
VIGRA_BLOCKWISE(GaussianDivergenceFunctor,               gaussianDivergenceMultiArray,               1, false );
VIGRA_BLOCKWISE(HessianOfGaussianFunctor,                hessianOfGaussianMultiArray,                2, false );
VIGRA_BLOCKWISE(HessianOfGaussianEigenvaluesFunctor,     hessianOfGaussianEigenvaluesMultiArray,     2, false );
VIGRA_BLOCKWISE(HessianOfGaussianFirstEigenvalueFunctor, hessianOfGaussianFirstEigenvalueMultiArray, 2, false );
VIGRA_BLOCKWISE(HessianOfGaussianLastEigenvalueFunctor,  hessianOfGaussianLastEigenvalueMultiArray,  2, false );
VIGRA_BLOCKWISE(LaplacianOfGaussianFunctor,              laplacianOfGaussianMultiArray,              2, false );
VIGRA_BLOCKWISE(GaussianGradientMagnitudeFunctor,        gaussianGradientMagnitudeMultiArray,        1, false );
VIGRA_BLOCKWISE(StructureTensorFunctor,                  structureTensorMultiArray,                  1, true  );

#undef  VIGRA_BLOCKWISE

    // alternative name for backward compatibility
template <unsigned int N, class T1, class S1, class T2, class S2>
inline void
gaussianGradientMagnitude(
    MultiArrayView<N, T1, S1> const & source,
    MultiArrayView<N, T2, S2> dest,
    BlockwiseConvolutionOptions<N> const & options)
{
    gaussianGradientMagnitudeMultiArray(source, dest, options);
}


} // end namespace vigra

#endif // VIGRA_MULTI_BLOCKWISE_HXX
