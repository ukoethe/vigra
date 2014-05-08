/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon                    */

/***************************************************************************
*              3D Adaptive Multiresolution Non-Local Means Filter          *
* Pierrick Coupe a, Jose V. Manjon, Montserrat Robles , D. Louis Collins   *
***************************************************************************/


/*                          Details on ONLM filter                        */
/***************************************************************************
 *  The ONLM filter is described in:                                       *
 *                                                                         *
 *  P. Coup�, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
 *  Avril 2008                                                             *
 ***************************************************************************/


#ifndef VIGRA_NON_LOCAL_MEAN
#define VIGRA_NON_LOCAL_MEAN

#include "math.h"
//#include "mex.h"
#include <stdlib.h>
//#include "matrix.h"

/*std*/
#include <iomanip>

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

/*boost*/
#include <boost/thread/thread.hpp>


/*vigra*/
#include "multi_array.hxx"
#include "multi_convolution.hxx"
#include "error.hxx"


namespace vigra{

struct NonLocalMeanParameter{

    NonLocalMeanParameter(

        const double sigma = 1.0,
        const int searchRadius = 3,
        const int patchRadius = 1,
        const bool gaussNoise = true,
        const double sigmaMean = 1.0,
        const int nThreads = 8,
        const double epsilon = 0.00001,
        const double mu1 = 0.95,
        const double var1 = 0.5,
        const int stepSize = 2,
        const bool verbose = true
    ):
    sigma_(sigma),
    searchRadius_(searchRadius),
    patchRadius_(patchRadius),
    gaussNoise_(gaussNoise),
    sigmaMean_(sigmaMean),
    nThreads_(nThreads),
    epsilon_(epsilon),
    mu1_(mu1),
    var1_(var1),
    stepSize_(stepSize),
    verbose_(verbose)
    {
    }

    double sigma_;
    int searchRadius_;
    int patchRadius_;
    bool gaussNoise_;
    double sigmaMean_;
    int nThreads_;


    double epsilon_;
    double mu1_;
    double var1_;
    int stepSize_;
    bool verbose_;
};

template<int DIM, class C, class S>
inline void mirrorIfIsOutsidePoint(
    const vigra::TinyVector<C,DIM> & shape,
    vigra::TinyVector<S,DIM> &      coord
){
    for(int c=0;c<DIM;++c){
        if(coord[c]<0)
            coord[c]=-1*coord[c];
        else if(coord[c]>= shape[c])
            coord[c] = 2 * shape[c] - coord[c] - 1;
    }
}

template<int DIM, class C, class S>
inline bool isOutsidePoint(
    const vigra::TinyVector<C,DIM> & shape,
    const vigra::TinyVector<S,DIM> &      coord
){
    for(int c=0;c<DIM;++c){
        if(coord[c]<0 || coord[c]>= shape[c])
            return true;
    }
    return false;
}

template<int DIM, class C, class S>
inline bool isInsidePoint(
    const vigra::TinyVector<C,DIM> & shape,
    const vigra::TinyVector<S,DIM> &      coord
){
    return !isOutsidePoint(shape,coord);
}

template<int DIM, class C, class S>
inline bool isZeroPoint(
    const vigra::TinyVector<C,DIM> & shape,
    const vigra::TinyVector<S,DIM> &      coord
){
    for(int c=0;c<DIM;++c){
        if(coord[c]!=0)
            return false;
    }
    return true;
}



template<class T,int SIZE>
void makeNegtiveValuesZero(
    vigra::TinyVector<T,SIZE> & a
){
    for(int i=0;i<SIZE;++i)
        if(a[i]<static_cast<T>(0))
            a[i]=static_cast<T>(0);
}

template<class T>
void makeNegtiveValuesZero(
    T & a
){
    if(a<static_cast<T>(0))
        a=static_cast<T>(0);
}



template<class T,int SIZE>
typename vigra::NumericTraits<T>::Promote pixelSum(
    const vigra::TinyVector<T,SIZE> & a
){
    return vigra::sum(a);
}

template<class T>
typename vigra::NumericTraits<T>::Promote  pixelSum(
    const T & a
){
    return a;
}



template<class T>
struct PixelTypeLength{
    const static int Length=1;
};


template<class T,int SIZE>
struct PixelTypeLength<vigra::TinyVector<T,SIZE> >{
    const static int Length=SIZE;
};






template<int DIM, class PIXEL_TYPE_IN, class PIXEL_TYPE_OUT >
class BockWiseNonLocalMeanImpl{
private:
    typedef BockWiseNonLocalMeanImpl<DIM,PIXEL_TYPE_IN,PIXEL_TYPE_OUT> SelfType;
public:
    typedef PIXEL_TYPE_IN       PixelTypeIn;
    typedef typename vigra::NumericTraits<PixelTypeIn>::RealPromote        RealPromotePixelType;  
    typedef typename vigra::NumericTraits<RealPromotePixelType>::ValueType RealPromoteScalarType;

    typedef PIXEL_TYPE_OUT      PixelTypeOut;
    typedef UInt32              LabelType;
    typedef typename vigra::MultiArray<DIM,int>::difference_type Coordinate;

    struct ThreadArguments{
        // data
        int start;
        int end;
        NonLocalMeanParameter param;
        vigra::MultiArrayView<DIM,PixelTypeIn>            inImage;
        vigra::MultiArrayView<DIM,RealPromotePixelType>   meansImage,variancesImage,estimateImage;
        vigra::MultiArrayView<DIM,LabelType>    labelImage;
        // to give verbose information
        int threadIndex;
        vigra::MultiArrayView<1,double>   progress;
    };


    struct ThreadHelper{
        #ifdef _WIN32
            typedef HANDLE HandleType; /* Handles to the worker threads*/

            static HandleType * makeThreadArray(const size_t numThreads){
                return (HandleType *)malloc(numThreads * sizeof( HandleType ));
            }
            static ThreadArguments * makeThreadArguments(const size_t numThreads){
                return (ThreadArguments *) malloc( numThreads * sizeof(ThreadArguments));
            }

            static void close(HandleType & handle){
                CloseHandle( handle );
            }

        #else
            typedef pthread_t HandleType;

            static HandleType * makeThreadArray(const size_t numThreads){
                return (HandleType *)calloc(numThreads , sizeof( HandleType ));
            }
            static ThreadArguments * makeThreadArguments(const size_t numThreads){
                return (ThreadArguments *) calloc(numThreads , sizeof(ThreadArguments));
            }

            static void close(HandleType & handle){
                pthread_join(handle, NULL);
            }
        #endif
    };





    static void run(
        const vigra::MultiArrayView<DIM,PixelTypeIn> & image,
        const NonLocalMeanParameter & parameter,
        vigra::MultiArrayView<DIM,PixelTypeOut> & outImage
    ){
        SelfType::inspectParameter(parameter);
        vigra::MultiArray<DIM,RealPromotePixelType> meanImage(image.shape());
        vigra::MultiArray<DIM,RealPromotePixelType> variancesImage(image.shape());
        vigra::MultiArray<DIM,RealPromotePixelType> estimageImage(image.shape());
        vigra::MultiArray<DIM,LabelType> labelImage(image.shape());

        labelImage = 0;

        // compute mean and variance with gaussianSmoothing
        vigra::gaussianSmoothMultiArray(image, meanImage, parameter.sigmaMean_);
        // square raw data 
        // - use estimate Image to store temp. results
        //   since it is already allocated
        // - No need to allocate another (huge) array
        for(int scanOrderIndex=0;scanOrderIndex<image.size();++scanOrderIndex){
            //using namespace vigra::multi_math;
            estimageImage[scanOrderIndex]=vigra::pow(image[scanOrderIndex],2);
            labelImage[scanOrderIndex]=0;
        }
        // mean of squared image 
        vigra::gaussianSmoothMultiArray(estimageImage,variancesImage, parameter.sigmaMean_);
        for(int scanOrderIndex=0;scanOrderIndex<image.size();++scanOrderIndex){
            RealPromotePixelType var = variancesImage[scanOrderIndex] - vigra::pow(meanImage[scanOrderIndex],2);
            // callbyref
            makeNegtiveValuesZero(var);
            variancesImage[scanOrderIndex] = var;
            // make estimageImage clean  again!
            // THIS IS IMPORTANT ! (must be zero initalized)
            estimageImage[scanOrderIndex]=RealPromotePixelType(0.0);
        }

        // THREADING
        typename ThreadHelper::HandleType * threadArray =  ThreadHelper::makeThreadArray(parameter.nThreads_);
        ThreadArguments *threadArgArray = ThreadHelper::makeThreadArguments(parameter.nThreads_);

        // for verbose
        vigra::MultiArray<1,double> progress(typename vigra::MultiArray<1,double>::difference_type(parameter.nThreads_));



        // NEW THREADING
        typedef boost::thread ThreadType;
        std::vector<ThreadType *> threadVector(parameter.nThreads_);


        for (size_t i = 0; i < parameter.nThreads_; i++){

            threadArgArray[i].start = (i * image.shape(DIM-1)) / parameter.nThreads_;
            threadArgArray[i].end = ((i + 1) * image.shape(DIM-1)) / parameter.nThreads_;
            threadArgArray[i].param = parameter;
            threadArgArray[i].inImage = image;
            threadArgArray[i].meansImage = meanImage;
            threadArgArray[i].variancesImage = variancesImage;;
            threadArgArray[i].estimateImage = estimageImage;
            threadArgArray[i].labelImage = labelImage;

            threadArgArray[i].threadIndex = i;
            threadArgArray[i].progress = progress;

            threadVector[i]= new ThreadType(&SelfType::runThread,&threadArgArray[i]);
        }
        for (size_t i = 0; i < parameter.nThreads_; i++){
            threadVector[i]->join();
        }

        free(threadArgArray);
        free(threadArray);


        /* Aggregation of the estimators (i.e. means computation) */
        for(int scanOrderIndex=0; scanOrderIndex<labelImage.size(); ++scanOrderIndex){
            if (labelImage[scanOrderIndex] == 0)
                outImage[scanOrderIndex]=image[scanOrderIndex];
            else
                outImage[scanOrderIndex]=estimageImage[scanOrderIndex] / labelImage[scanOrderIndex];
        }

        for (size_t i = 0; i < parameter.nThreads_; i++){
            delete threadVector[i];
        }

    }

private:
    typedef vigra::MultiArrayView<DIM, PixelTypeIn> InArrayView;
    typedef vigra::MultiArrayView<DIM, PixelTypeIn> OutArrayView;
    typedef vigra::MultiArrayView<DIM, RealPromotePixelType>   RealPromotedArrayView;
    typedef vigra::MultiArrayView<DIM, LabelType>   LabelArrayView;

    typedef std::vector<RealPromotePixelType>  BlockAverageVector;


    static void inspectParameter(const NonLocalMeanParameter & param){
        vigra_precondition(param.stepSize_>=1,
            "NonLocalMean Parameter: \"stepSize>=1\" violated");
        vigra_precondition(param.searchRadius_>=1,
            "NonLocalMean Parameter: \"searchRadius >=1\" violated");
        vigra_precondition(param.patchRadius_>=1,
            "NonLocalMean Parameter: \"searchRadius >=1\" violated");
        vigra_precondition(param.stepSize_-1<=param.patchRadius_,
            "NonLocalMean Parameter: \"stepSize -1 <= patchRadius\"  violated");

    }


    static void *runThread( void *pArguments ){
        ThreadArguments arg = *(ThreadArguments *) pArguments;
        const int patchRadius  = arg.param.patchRadius_;
        const int Ndims = std::pow( (2 * patchRadius + 1),DIM);
        BlockAverageVector average(Ndims);

        Coordinate xyz;
        Coordinate shape = arg.inImage.shape();
        const int start = arg.start;
        const int end = arg.end;
        const int stepSize = arg.param.stepSize_;
        const int outerSize =  end-start;
        const int threadIndex = arg.threadIndex;
        const int numThreads = arg.progress.size();
        const bool verbose = arg.param.verbose_;
        size_t counter=0;

        size_t jobSize = (outerSize/stepSize+1);
        for(int d=0;d<DIM-1;++d){
            jobSize*=(shape[d]/stepSize+1);
        }

        if(verbose && threadIndex==numThreads-1){
            std::cout<<"progress"; 
        }

        if(DIM==2){
            for (xyz[1] = start; xyz[1]  < end;    xyz[1]  += stepSize)
            for (xyz[0] = 0;   xyz[0]  < shape[0]; xyz[0]  += stepSize){
                SelfType::processSinglePixel(arg,average,xyz);
                if(verbose){
                    ++counter;
                    SelfType::progressPrinter(arg,jobSize,counter);
                }
            }
        }
        if(DIM==3){
            for (xyz[2] = start; xyz[2]  < end;    xyz[2]  += stepSize)
            for (xyz[1] = 0;   xyz[1]  < shape[1]; xyz[1]  += stepSize)
            for (xyz[0] = 0;   xyz[0]  < shape[0]; xyz[0]  += stepSize){
                SelfType::processSinglePixel(arg,average,xyz);
                if(verbose){
                    ++counter;
                    SelfType::progressPrinter(arg,jobSize,counter);
                }
            }
        }
        if(DIM==4){
            for (xyz[3] = start; xyz[3]  < end;    xyz[3]  += stepSize)
            for (xyz[2] = 0;   xyz[2]  < shape[2]; xyz[2]  += stepSize)
            for (xyz[1] = 0;   xyz[1]  < shape[1]; xyz[1]  += stepSize)
            for (xyz[0] = 0;   xyz[0]  < shape[0]; xyz[0]  += stepSize){
                SelfType::processSinglePixel(arg,average,xyz);
                if(verbose){
                    ++counter;
                    SelfType::progressPrinter(arg,jobSize,counter);
                }
            }
        }
        if(verbose &&  threadIndex==numThreads-1){
            std::cout<<"\rprogress "<<std::setw(10)<<"100"<<" %%"<<"\n";
        }

        #ifdef _WIN32
        //_endthreadex(0);
        #else
        //pthread_exit(0);
        #endif

        

        return 0;
    }


    static void progressPrinter(ThreadArguments & arg,const size_t jobSize, const size_t counter){

        const int threadIndex = arg.threadIndex;
        const int numThreads = arg.progress.size();

        for(size_t ti=0;ti<numThreads;++ti)
            arg.progress[ti] = double(counter+1)/double(jobSize);
        if(threadIndex==numThreads-1){
            if(counter%100 == 0){
                double pr=0;
                for(size_t ti=0;ti<numThreads;++ti)
                    pr+=arg.progress[ti];
                pr/=numThreads;
                pr*=100.0;
                std::cout<<"\rprogress "<<std::setw(10)<<pr<<" %%"<<std::flush;
            }
        }
    }



    static void processSinglePixel(
        ThreadArguments & arg,
        BlockAverageVector & average,
        const Coordinate & xyz
    ){
            Coordinate nxyz,xxyyzz;
            const Coordinate & shape = arg.inImage.shape();
            const bool rician = !arg.param.gaussNoise_;
            const int searchRadius = arg.param.searchRadius_;
            const int patchRadius  = arg.param.patchRadius_;
            const double sigma = arg.param.sigma_;
            const double epsilon = arg.param.epsilon_;
            const double bias = 2 * sigma * sigma;

            std::fill(average.begin(),average.end(),RealPromotePixelType(0.0));
            RealPromoteScalarType totalweight(0.0);

            if( pixelSum(arg.meansImage[xyz]) > epsilon  && pixelSum(arg.variancesImage[xyz]) > epsilon){
                RealPromoteScalarType wmax = 0.0;


                if(DIM==2){
                    for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                    for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                        nxyz = xyz  + xxyyzz;
                        if(isZeroPoint(shape,xxyyzz)){
                            continue;
                        }
                        SelfType::processSinglePair(arg,average,xyz,nxyz,shape,wmax,totalweight);
                    }
                }
                else if(DIM==3){
                    for (xxyyzz[2] = -searchRadius; xxyyzz[2] <= searchRadius; xxyyzz[2]++)
                    for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                    for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                        nxyz = xyz  + xxyyzz;
                        if(isZeroPoint(shape,xxyyzz)){
                            continue;
                        }
                        SelfType::processSinglePair(arg,average,xyz,nxyz,shape,wmax,totalweight);
                    }
                }
                else if(DIM==4){
                    for (xxyyzz[3] = -searchRadius; xxyyzz[3] <= searchRadius; xxyyzz[3]++)
                    for (xxyyzz[2] = -searchRadius; xxyyzz[2] <= searchRadius; xxyyzz[2]++)
                    for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                    for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                        nxyz = xyz  + xxyyzz;
                        if(isZeroPoint(shape,xxyyzz)){
                            continue;
                        }
                        SelfType::processSinglePair(arg,average,xyz,nxyz,shape,wmax,totalweight);
                    }
                }


                if (wmax == 0.0){
                    wmax = 1.0;
                }
                // give pixel xyz  as much weight as
                // the maximum weighted other patch
                SelfType::extractPatchValues(arg.inImage, xyz, patchRadius, average, wmax,rician);
                totalweight += wmax;

                // this if seems total useless to me
                if (totalweight != 0.0){
                    SelfType::writePatchValues(arg.estimateImage, arg.labelImage, xyz, patchRadius, average, totalweight, bias,rician);
                }

            }
            else{
                const double wmax = 1.0;
                SelfType::extractPatchValues(arg.inImage, xyz, patchRadius, average, wmax,rician);
                totalweight += wmax;
                SelfType::writePatchValues(arg.estimateImage, arg.labelImage, xyz, patchRadius, average, totalweight, bias,rician);
            }
    }


    static void processSinglePair(
        ThreadArguments & arg,
        BlockAverageVector & average,
        const Coordinate & xyz,
        const Coordinate & nxyz,
        const Coordinate & shape,
        RealPromoteScalarType & wmax,
        RealPromoteScalarType & totalweight
    ){
            if(isInsidePoint(shape,nxyz)){
                const double epsilon = arg.param.epsilon_;
                if( pixelSum(arg.meansImage[nxyz]) > epsilon  && pixelSum(arg.variancesImage[nxyz]) > epsilon){

                    const bool rician = !arg.param.gaussNoise_;
                    const int patchRadius  = arg.param.patchRadius_;
                    const double sigma = arg.param.sigma_;
                    const double mu1 = arg.param.mu1_;
                    const double var1 = arg.param.var1_;
                    const RealPromoteScalarType hh = sigma * sigma;

                    // Compute Ratios
                    const double t1 = pixelSum(arg.meansImage[xyz]/arg.meansImage[nxyz])/ PixelTypeLength<RealPromotePixelType>::Length;
                    const double t2 = pixelSum(arg.variancesImage[xyz]/arg.variancesImage[nxyz])/ PixelTypeLength<RealPromotePixelType>::Length;

                    // here we check if to patches fit to each other
                    // one patch is around xyz 
                    // other patch is arround nxyz                
                    if (t1 > mu1 && t1 < (1.0 / mu1) && t2 > var1 && t2 < (1.0 / var1)){

                        const double d = SelfType::distanceFunction(arg.inImage, xyz, nxyz, patchRadius);
                        const RealPromoteScalarType w = exp(-d / (hh));
                        wmax = std::max(w,wmax);
                        SelfType::extractPatchValues(arg.inImage, nxyz, patchRadius, average, w,rician);
                        totalweight = totalweight + w;
                    }
                }
            }
    }


    static double distanceFunction( const InArrayView & ima, const Coordinate  & pA,
                                    const Coordinate  & pB,const int f
    ){
        Coordinate offset,nPa,nPb;
        double acu = 0;
        double distancetotal = 0;

        #define VIGRA_NLM_IN_LOOP_CODE                          \
            nPa = pA+offset;                                    \
            nPb = pB+offset;                                    \
            mirrorIfIsOutsidePoint(ima.shape(),nPa);            \
            mirrorIfIsOutsidePoint(ima.shape(),nPb);            \
            const RealPromotePixelType vA = ima[nPa];           \
            const RealPromotePixelType vB = ima[nPb];           \
            distancetotal += vigra::squaredNorm(vA-vB)/         \
                PixelTypeLength<RealPromotePixelType>::Length ; \
            acu = acu + 1

        if(DIM==2){
            for (offset[1] = -f; offset[1] <= f; offset[1]++)
            for (offset[0] = -f; offset[0] <= f; offset[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        else if(DIM==3){
            for (offset[2] = -f; offset[2] <= f; offset[2]++)
            for (offset[1] = -f; offset[1] <= f; offset[1]++)
            for (offset[0] = -f; offset[0] <= f; offset[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        else if(DIM==4){
            for (offset[3] = -f; offset[3] <= f; offset[3]++)
            for (offset[2] = -f; offset[2] <= f; offset[2]++)
            for (offset[1] = -f; offset[1] <= f; offset[1]++)
            for (offset[0] = -f; offset[0] <= f; offset[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        #undef VIGRA_NLM_IN_LOOP_CODE
        return distancetotal / acu;
    }


    /* Function which computes the value assigned to each voxel */
    static void writePatchValues( RealPromotedArrayView & estimateImage, LabelArrayView & labelImage,
                            const Coordinate xyz, int neighborhoodsize, 
                            BlockAverageVector & average, const double global_sum, const double bias,
                            const bool rician
    ){
        Coordinate abc,xyzPos,nhSize(neighborhoodsize);        
        int count = 0 ;
        const int ns = 2 * neighborhoodsize + 1;

        #define VIGRA_NLM_IN_LOOP_CODE                                                      \
                xyzPos = xyz + abc - nhSize;                                                \
                if (isInsidePoint(estimateImage.shape(),xyzPos)){                           \
                    RealPromotePixelType value = estimateImage[xyzPos];                     \
                    if(rician){                                                             \
                        RealPromotePixelType denoisedValue =                                \
                            (average[count] / global_sum) -  RealPromotePixelType(bias);    \
                        makeNegtiveValuesZero(denoisedValue);                               \
                        denoisedValue=vigra::sqrt(denoisedValue);                           \
                        value += denoisedValue;                                             \
                    }                                                                       \
                    else{                                                                   \
                        value += average[count] / global_sum;                               \
                    }                                                                       \
                    estimateImage[xyzPos] = value;                                          \
                    ++labelImage[xyzPos];                                                   \
                }                                                                           \
                count++

        if(DIM==2){
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        if(DIM==3){
            for (abc[2] = 0; abc[2] < ns; abc[2]++)
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        if(DIM==4){
            for (abc[3] = 0; abc[3] < ns; abc[3]++)
            for (abc[2] = 0; abc[2] < ns; abc[2]++)
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        #undef VIGRA_NLM_IN_LOOP_CODE
    }

    // extract values from a patch 
    // and write them into average vector 
    // weighted!
    static void extractPatchValues(const InArrayView & ima, const Coordinate  & xyz,
        const int neighborhoodsize, BlockAverageVector & average, const double weight,const bool rician
    ){
        Coordinate xyzPos,abc;
        Coordinate nhSize3(neighborhoodsize);
        const int ns = 2 * neighborhoodsize + 1;
        int count = 0;

        // todo: remove abc vector

        #define VIGRA_NLM_IN_LOOP_CODE                                  \
            xyzPos = xyz + abc - nhSize3;                               \
            if (rician){                                                \
                if (isOutsidePoint(ima.shape(),xyzPos))                 \
                    average[count] += ima[xyz]*ima[xyz]* weight;        \
                else                                                    \
                    average[count] += ima[xyzPos]*ima[xyzPos]* weight;  \
            }                                                           \
            else{                                                       \
                if (isOutsidePoint(ima.shape(),xyzPos))                 \
                    average[count] += ima[xyz]* weight;                 \
                else                                                    \
                    average[count] += ima[xyzPos]* weight;              \
            }                                                           \
            count++

        if(DIM==2){
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        else if(DIM==3){
            for (abc[2] = 0; abc[2] < ns; abc[2]++)
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }
        else if(DIM==4){
            for (abc[3] = 0; abc[3] < ns; abc[3]++)
            for (abc[2] = 0; abc[2] < ns; abc[2]++)
            for (abc[1] = 0; abc[1] < ns; abc[1]++)
            for (abc[0] = 0; abc[0] < ns; abc[0]++){
                VIGRA_NLM_IN_LOOP_CODE;
            }
        }

        #undef VIGRA_NLM_IN_LOOP_CODE
    }

};





template<int DIM, class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT>
void nonLocalMean(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & image,
    const NonLocalMeanParameter parameter,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> outImage
){
    BockWiseNonLocalMeanImpl<DIM,PIXEL_TYPE_IN,PIXEL_TYPE_OUT>::run(image,parameter,outImage);
}


} // end namespace vigra


#endif