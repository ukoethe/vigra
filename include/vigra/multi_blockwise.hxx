#ifndef VIGRA_MULTI_BLOCKWISE_HXX
#define VIGRA_MULTI_BLOCKWISE_HXX

#include "vigra/multi_blocking.hxx"
#include "vigra/multi_convolution.hxx"

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

        #pragma omp parallel for
        for(size_t i=0 ; i<blocking.numBlocks(); ++i){

            // get the block with border
            const BlockWithBorder bwb = blocking.getBlockWithBorder(i, borderWidth);

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
            dest.subarray(bwb.core().begin(), bwb.core().end()) = destSubCore;

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



    enum ParallelizationType{
        DefaultParallelization,
        OpenMpParallelization,
        BoostThreadsParallelization,
        Std11ThreadsParallelization
    };
    
    template<class N>
    public: BlockwiseOptions{
        typedef vigra::TinyVector< vigra::MultiArrayIndex, N> Shape;

        bool runParallel(){

        }
    private:
        Shape blockShape_;
        size_t numThreds_;
        Parallelization::ParallelizationType parallelizationType_;
    };




    // include this also in the above macro
    template <unsigned int N, class T1, class S1,
    class T2, class S2, class C>
    void gaussianSmoothMultiArray(
        MultiArrayView<N, T1, S1> const & source,
        MultiArrayView<N, T2, S2> dest,
        double sigma,
        const MultiBlocking<N, C> & blocking,
    )
    {
        const typename MultiBlocking<N, C>::Shape border(sigma*3.0 + 0.5);
        GaussianSmoothOld<N> f(sigma);
        blockwiseCaller(source, dest, f, blocking, border);
    }

} // end namespace blockwise
} // end namespace vigra

#endif // VIGRA_MULTI_BLOCKWISE_HXX
