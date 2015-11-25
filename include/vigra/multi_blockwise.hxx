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
#include "vigra/multi_blocking.hxx"
#include "vigra/multi_convolution.hxx"
#include "vigra/multi_tensorutilities.hxx"

#ifndef VIGRA_DEFAULT_BLOCK_SHAPE
    #define VIGRA_DEFAULT_BLOCK_SHAPE 64
#endif


namespace vigra{

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
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth
    ){

        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;
        typedef typename MultiBlocking<DIM, C>::BlockWithBorderIter BlockWithBorderIter;

// FIXME: replace with threadpool        #pragma omp parallel
        {
            BlockWithBorderIter iter  =  blocking.blockWithBorderBegin(borderWidth);
            //std::cout<<"blockshape "<<(*iter).core().size()<<"\n";

// FIXME: replace with threadpool            #pragma omp for
            for(int i=0 ; i<blocking.numBlocks(); ++i){

                const BlockWithBorder bwb = iter[i];

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
            }
        }
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
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth
    ){

        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;
        typedef typename MultiBlocking<DIM, C>::BlockWithBorderIter BlockWithBorderIter;
        typedef typename MultiBlocking<DIM, C>::Block Block;
// FIXME: replace with threadpool        #pragma omp parallel
        {
            BlockWithBorderIter iter  =  blocking.blockWithBorderBegin(borderWidth);
            //std::cout<<"blockshape "<<(*iter).core().size()<<"\n";

// FIXME: replace with threadpool            #pragma omp for
            for(int i=0 ; i<blocking.numBlocks(); ++i){

                const BlockWithBorder bwb = iter[i];

                // get the input of the block as a view
                vigra::MultiArrayView<DIM, T_IN, ST_IN> sourceSub = source.subarray(bwb.border().begin(),
                                                                            bwb.border().end());

                // get the output of the blocks core as a view
                vigra::MultiArrayView<DIM, T_OUT, ST_OUT> destCore = dest.subarray(bwb.core().begin(),
                                                                            bwb.core().end());


                const Block localCore =  bwb.localCore();


                // call the functor
                functor(sourceSub, destCore, localCore.begin(), localCore.end());

                // write the core global out
                //vigra::MultiArrayView<DIM, T_OUT, ST_OUT> destSubCore = destSub.subarray(bwb.localCore().begin(),
                //                                                                bwb.localCore().end());
                //dest.subarray(bwb.core().begin()-blocking.roiBegin(),
                //              bwb.core().end()  -blocking.roiBegin()  ) = destSubCore;
            }
        }
    }

    #define CONVOLUTION_FUNCTOR(FUNCTOR_NAME, FUNCTION_NAME) \
    template<unsigned int DIM> \
    class FUNCTOR_NAME{ \
    public: \
        typedef ConvolutionOptions<DIM> ConvOpt; \
        FUNCTOR_NAME(const ConvOpt & convOpt) \
        : convOpt_(convOpt){} \
        template<class S, class D> \
        void operator()(const S & s, D & d)const{ \
            FUNCTION_NAME(s, d, convOpt_); \
        } \
        template<class S, class D,class SHAPE> \
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){ \
            ConvOpt convOpt(convOpt_); \
            convOpt.subarray(roiBegin, roiEnd); \
            FUNCTION_NAME(s, d, convOpt); \
        } \
    private: \
        ConvOpt  convOpt_; \
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
        : convOpt_(convOpt){}
        template<class S, class D>
        void operator()(const S & s, D & d)const{
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(d.shape());
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, convOpt_);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, d);
        }
        template<class S, class D,class SHAPE>
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(roiEnd-roiBegin);
            convOpt_.subarray(roiBegin, roiEnd);
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, convOpt_);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, d);
        }
    private:
        ConvOpt  convOpt_;
    };

    template<unsigned int DIM, unsigned int EV>
    class HessianOfGaussianSelectedEigenvalueFunctor{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;
        HessianOfGaussianSelectedEigenvalueFunctor(const ConvOpt & convOpt)
        : convOpt_(convOpt){}
        template<class S, class D>
        void operator()(const S & s, D & d)const{
            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;

            // compute the hessian of gaussian and extract eigenvalue
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(s.shape());
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, convOpt_);

            vigra::MultiArray<DIM, TinyVector<RealType, DIM > >  allEigenvalues(s.shape());
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, allEigenvalues);

            d = allEigenvalues.bindElementChannel(EV);
        }
        template<class S, class D,class SHAPE>
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){

            typedef typename vigra::NumericTraits<typename S::value_type>::RealPromote RealType;

            // compute the hessian of gaussian and extract eigenvalue
            vigra::MultiArray<DIM, TinyVector<RealType, int(DIM*(DIM+1)/2)> >  hessianOfGaussianRes(roiEnd-roiBegin);
            convOpt_.subarray(roiBegin, roiEnd);
            vigra::hessianOfGaussianMultiArray(s, hessianOfGaussianRes, convOpt_);

            vigra::MultiArray<DIM, TinyVector<RealType, DIM > >  allEigenvalues(roiEnd-roiBegin);
            vigra::tensorEigenvaluesMultiArray(hessianOfGaussianRes, allEigenvalues);

            d = allEigenvalues.bindElementChannel(EV);
        }
    private:
        ConvOpt  convOpt_;
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

    /// concurrency type to use within parallel algorithms
    ///
    /// So far, only OpenMpConcurrency is implemented
    enum ConcurrencyType{
        DefaultConcurrency,
        OpenMpConcurrency,
        BoostThreadsConcurrency,
        Std11ThreadsConcurrency,
        NoConcurrency

    };

    /// base option class for parallel algorithms
    class ParallelOptions{
    public:
        ParallelOptions(const size_t numThreds = 0,
                        const ConcurrencyType  concurrencyType = OpenMpConcurrency)
        :   numThreads_(numThreds), // zero means AUTO
            concurrencyType_(concurrencyType){
                if(concurrencyType_!=OpenMpConcurrency){
                    throw std::runtime_error("currently only OpenMpConcurrency is implemented");
                }
        }
        size_t getNumThreads()const{
            return numThreads_;
        }
        void setNumThreads(const size_t numThreads){
            numThreads_ = numThreads;
        }
        ConcurrencyType getConcurrencyType()const{
            return concurrencyType_;
        }
        void setConcurencyType(const ConcurrencyType & concurrencyType){
            concurrencyType_ = concurrencyType;
        }
    private:
        size_t numThreads_;
        ConcurrencyType concurrencyType_;
    };

    /// base option class for blockwise algorithms
    /// attaches blockshape to ParallelOptions
    template<unsigned int N>
    class BlockwiseOptions
    : public ParallelOptions
    {
    public:
        typedef vigra::TinyVector< vigra::MultiArrayIndex, N> Shape;

        BlockwiseOptions(const Shape & blockShape = Shape(VIGRA_DEFAULT_BLOCK_SHAPE))
        :   ParallelOptions(),
            blockShape_(blockShape){
        }
        Shape getBlockShape()const{
            return blockShape_;
        }
        void setBlockShape(const Shape & blockShape){
            blockShape_ = blockShape;
        }
    private:
        Shape blockShape_;
    };

    template<unsigned int N>
    class BlockwiseConvolutionOptions
    :   public  BlockwiseOptions<N>, public vigra::ConvolutionOptions<N>{
    public:
        BlockwiseConvolutionOptions()
        :   BlockwiseOptions<N>(),
            vigra::ConvolutionOptions<N>(){
        }
    private:

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


    #define BLOCKWISE_FUNCTION_GEN(FUNCTOR, FUNCTION, ORDER, USES_OUTER_SCALE) \
    template <unsigned int N, class T1, class S1, class T2, class S2> \
    void FUNCTION( \
        MultiArrayView<N, T1, S1> const & source, \
        MultiArrayView<N, T2, S2> dest, \
        const BlockwiseConvolutionOptions<N> & options \
    ) \
    {  \
        typedef  MultiBlocking<N, vigra::MultiArrayIndex> Blocking; \
        typedef typename Blocking::Shape Shape; \
        const Shape border = getBorder(options, ORDER, USES_OUTER_SCALE); \
        BlockwiseConvolutionOptions<N> subOptions(options); \
        subOptions.subarray(Shape(0), Shape(0));  \
        const Blocking blocking(source.shape(), options.getBlockShape()); \
        FUNCTOR f(subOptions); \
        blockwiseCaller(source, dest, f, blocking, border); \
    }

    BLOCKWISE_FUNCTION_GEN(GaussianSmoothFunctor<N> ,                   gaussianSmoothMultiArray,                   0, false );
    BLOCKWISE_FUNCTION_GEN(GaussianGradientFunctor<N> ,                 gaussianGradientMultiArray,                 1, false );
    BLOCKWISE_FUNCTION_GEN(SymmetricGradientFunctor<N> ,                symmetricGradientMultiArray,                1, false );
    BLOCKWISE_FUNCTION_GEN(GaussianDivergenceFunctor<N> ,               gaussianDivergenceMultiArray,               1, false );
    BLOCKWISE_FUNCTION_GEN(HessianOfGaussianFunctor<N> ,                hessianOfGaussianMultiArray,                2, false );
    BLOCKWISE_FUNCTION_GEN(HessianOfGaussianEigenvaluesFunctor<N> ,     hessianOfGaussianEigenvaluesMultiArray,     2, false );
    BLOCKWISE_FUNCTION_GEN(HessianOfGaussianFirstEigenvalueFunctor<N> , hessianOfGaussianFirstEigenvalueMultiArray, 2, false );
    BLOCKWISE_FUNCTION_GEN(HessianOfGaussianLastEigenvalueFunctor<N> ,  hessianOfGaussianLastEigenvalueMultiArray,  2, false );
    BLOCKWISE_FUNCTION_GEN(LaplacianOfGaussianFunctor<N> ,              laplacianOfGaussianMultiArray,              2, false );
    BLOCKWISE_FUNCTION_GEN(GaussianGradientMagnitudeFunctor<N>,         gaussianGradientMagnitudeMultiArray,        1, false );
    BLOCKWISE_FUNCTION_GEN(StructureTensorFunctor<N> ,                  structureTensorMultiArray,                  1, true  );


    #undef  BLOCKWISE_FUNCTION_GEN

} // end namespace blockwise
} // end namespace vigra

#endif // VIGRA_MULTI_BLOCKWISE_HXX
