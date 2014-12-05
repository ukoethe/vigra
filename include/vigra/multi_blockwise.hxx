#ifndef VIGRA_MULTI_BLOCKWISE_HXX
#define VIGRA_MULTI_BLOCKWISE_HXX

#include "vigra/multi_blocking.hxx"
#include "vigra/multi_convolution.hxx"


#ifndef VIGRA_DEFAULT_BLOCK_SHAPE 
    #define VIGRA_DEFAULT_BLOCK_SHAPE 64
#endif 

namespace vigra{

/*
    Blockwise needs to implement:
    
    Simple Element Wise (with reduction):
        min
        max
        minmax

    Simple Convolution:
        gaussianSmoothMultiArray
        gaussianGradientMultiArray
        symmetricGradientMultiArray
        gaussianDivergenceMultiArray
        hessianOfGaussianMultiArray
        laplacianOfGaussianMultiArray
        gaussianGradientMagnitude
        structureTensorMultiArray

    Tensor Related
        hessianOfGaussianEigenvalues
        hessianOfGaussianTrace
        structureTensorEigenvalues
        structureTensoTrace
        tensorEigenvalue
        tensorTrace

    Distance Transform:
        truncatedDistanceTransform

    ParabolicFilters: 

*/
namespace blockwise{

    template<
        unsigned int DIM,
        class T_IN, class ST_IN,
        class T_OUT, class ST_OUT,
        class FILTER_FUCTOR,
        class C
    >
    void blockwiseCaller(
        const vigra::MultiArrayView<DIM, T_IN,  ST_IN > & source,
        const vigra::MultiArrayView<DIM, T_OUT, ST_OUT> & dest,
        FILTER_FUCTOR & functor,
        const vigra::MultiBlocking<DIM, C> & blocking,
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth
    ){



        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;
        typedef typename MultiBlocking<DIM, C>::BlockWithBorderIter BlockWithBorderIter;

        


        #pragma omp parallel
        {
            BlockWithBorderIter iter  =  blocking.blockWithBorderBegin(borderWidth);
            //std::cout<<"blockshape "<<(*iter).core().size()<<"\n";

            #pragma omp for
            for(size_t i=0 ; i<blocking.numBlocks(); ++i){

                const BlockWithBorder bwb = iter[i];

                // get the input of the block as a view
                vigra::MultiArrayView<DIM, T_IN, ST_IN> sourceSub = source.subarray(bwb.border().begin(),
                                                                             bwb.border().end());

                // get the output as NEW allocated array
                vigra::MultiArray<DIM, T_OUT> destSub(sourceSub.shape());

                // call the functor
                functor(sourceSub, destSub);

                 // write the core global out
                vigra::MultiArrayView<DIM, T_IN, ST_IN> destSubCore = destSub.subarray(bwb.localCore().begin(),
                                                                                bwb.localCore().end());
                // write the core global out
                dest.subarray(bwb.core().begin()-blocking.roiBegin(), 
                              bwb.core().end()  -blocking.roiBegin()  ) = destSubCore;

            }
        }
    }


    template<
        unsigned int DIM,
        class T_IN, class ST_IN,
        class T_OUT, class ST_OUT,
        class FILTER_FUCTOR,
        class C
    >
    void blockwiseCallerNew(
        const vigra::MultiArrayView<DIM, T_IN,  ST_IN > & source,
        const vigra::MultiArrayView<DIM, T_OUT, ST_OUT> & dest,
        FILTER_FUCTOR & functor,
        const vigra::MultiBlocking<DIM, C> & blocking,
        const typename vigra::MultiBlocking<DIM, C>::Shape & borderWidth
    ){



        typedef typename MultiBlocking<DIM, C>::BlockWithBorder BlockWithBorder;
        typedef typename MultiBlocking<DIM, C>::BlockWithBorderIter BlockWithBorderIter;
        typedef typename BlockWithBorder::Shape Shape;
        
        const Shape & totalRoiBegin = blocking.roiBegin();

        //#pragma omp parallel
        {
            BlockWithBorderIter iter  =  blocking.blockWithBorderBegin(borderWidth);
            //std::cout<<"blockshape "<<(*iter).core().size()<<"\n";

            #pragma omp for
            for(size_t i=0 ; i<blocking.numBlocks(); ++i){

                const BlockWithBorder bwb = iter[i];

                // get the input of the block as a view
                vigra::MultiArrayView<DIM, T_IN, ST_IN> sourceSub = source.subarray(bwb.border().begin()-totalRoiBegin,
                                                                                    bwb.border().end()-totalRoiBegin);


                    
                vigra::MultiArrayView<DIM, T_OUT, ST_OUT>  desSubCore = dest.subarray(bwb.core().begin()-totalRoiBegin, 
                                                                                      bwb.core().end()  -totalRoiBegin);


                // call the functor
                const Shape roiBegin = bwb.localCore().begin();
                const Shape roiEnd = bwb.localCore().end();
                functor(sourceSub, desSubCore, roiBegin, roiEnd);

            }
        }
    }



    template<unsigned int DIM>
    struct GaussianSmoothOld{
    public:
        typedef ConvolutionOptions<DIM> ConvOpt;

        GaussianSmoothOld(const double sigma)
        :   sigma_(sigma){
        }
        template<class S, class D>
        void operator()(const S & s, D & d)const{
            vigra::gaussianSmoothMultiArray(s, d, sigma_);
        }
    private:
        const double sigma_;
    };


    #define CONVOLUTION_FUNCTOR(FUCTOR_NAME, FUNCTION_NAME) \
    template<unsigned int DIM> \
    class FUCTOR_NAME{ \
    public: \
        typedef ConvolutionOptions<DIM> ConvOpt; \
        FUCTOR_NAME(const ConvOpt & convOpt) \
        : convOpt_(convOpt){} \
        template<class S, class D> \
        void operator()(const S & s, D & d)const{ \
            FUNCTION_NAME(s, d, convOpt_); \
        } \
        template<class S, class D,class SHAPE> \
        void operator()(const S & s, D & d, const SHAPE & roiBegin, const SHAPE & roiEnd){ \
            convOpt_.subarray(roiBegin, roiEnd); \
            FUNCTION_NAME(s, d, convOpt_); \
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

    enum ConcurrencyType{
        DefaultConcurrency,
        OpenMpConcurrency,
        BoostThreadsConcurrency,
        Std11ThreadsConcurrency,
        NoConcurrency

    };
        

    class ParallelOptions{
    public:
        ParallelOptions(const size_t numThreds = 0, 
                        const ConcurrencyType  concurrencyType = DefaultConcurrency)
        :   numThreads_(numThreds), // zero means AUTO
            concurrencyType_(concurrencyType){
        }
        size_t getNumThreads()const{
            return numThreads_;
        }
        ConcurrencyType getConcurrencyType()const{
            return concurrencyType_;
        }
    private:
        size_t numThreads_;
        ConcurrencyType concurrencyType_;
    };


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



    template<unsigned int N>
    vigra::TinyVector< vigra::MultiArrayIndex, N > getBorder(
        const BlockwiseConvolutionOptions<N> & opt,
        const size_t order,
        const bool usesOuterScale = false
    ){
        vigra::TinyVector< vigra::MultiArrayIndex, N > res(vigra::SkipInitialization);
        if(opt.getFilterWindowSize()<=0.00001){
            for(size_t d=0; d<N; ++d){
                const double stdDev =  opt.getStdDev()[d];
                res[d] = std::round(3.0 * stdDev  + 0.5*static_cast<double>(order));
            }
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


    BLOCKWISE_FUNCTION_GEN(GaussianSmoothFunctor<N> ,           gaussianSmoothMultiArray,       0, false );
    BLOCKWISE_FUNCTION_GEN(GaussianGradientFunctor<N> ,         gaussianGradientMultiArray,     1, false );
    BLOCKWISE_FUNCTION_GEN(SymmetricGradientFunctor<N> ,        symmetricGradientMultiArray,    1, false );
    BLOCKWISE_FUNCTION_GEN(GaussianDivergenceFunctor<N> ,       gaussianDivergenceMultiArray,   2, false );
    BLOCKWISE_FUNCTION_GEN(HessianOfGaussianFunctor<N> ,        hessianOfGaussianMultiArray,    2, false );
    BLOCKWISE_FUNCTION_GEN(LaplacianOfGaussianFunctor<N> ,      laplacianOfGaussianMultiArray,  2, false );
    BLOCKWISE_FUNCTION_GEN(GaussianGradientMagnitudeFunctor<N>, gaussianGradientMagnitude,      1, false );
    BLOCKWISE_FUNCTION_GEN(StructureTensorFunctor<N> ,          structureTensorMultiArray,      1, true  );


    #undef  BLOCKWISE_FUNCTION_GEN

} // end namespace blockwise
} // end namespace vigra

#endif // VIGRA_MULTI_BLOCKWISE_HXX
