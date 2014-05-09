/************************************************************************/
/*                                                                      */
/*     Copyright 2014 by Thorsten Beier and Ullrich Koethe              */
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
/* RE-IMPLEMENTAION OF THE WORK OF:                                     */
/* Pierrick Coupe - pierrick.coupe@gmail.com                            */
/* Jose V. Manjon - jmanjon@fis.upv.es                                  */
/* Brain Imaging Center, Montreal Neurological Institute.               */
/* Mc Gill University                                                   */
/*                                                                      */
/************************************************************************/
/* The ONLM filter is described in:                                     */
/*                                                                      */
/* P. Coup�, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.   */
/* An Optimized Blockwise Non Local Means Denoising Filter              */ 
/* for 3D Magnetic Resonance Images                                     */
/* . IEEE Transactions on Medical Imaging, 27(4):425-441,               */
/* Avril 2008                                                           */
/************************************************************************/


#ifndef VIGRA_NON_LOCAL_MEAN
#define VIGRA_NON_LOCAL_MEAN


/*std*/
#include <iomanip>

/*vigra*/
#include "multi_array.hxx"
#include "multi_convolution.hxx"
#include "error.hxx"
#include "threading.hxx"


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
        const double meanRatio = 0.95,
        const double varRatio = 0.5,
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
    meanRatio_(meanRatio),
    varRatio_(varRatio),
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
    double meanRatio_;
    double varRatio_;
    int stepSize_;
    bool verbose_;
};



template<class PIXEL_TYPE_IN>
class SelectByRatio{
public:
    typedef typename NumericTraits<PIXEL_TYPE_IN>::RealPromote PixelType;
    typedef typename NumericTraits<PIXEL_TYPE_IN>::ValueType   ValueType;


    SelectByRatio(const ValueType  meanRatio,  const ValueType  varRatio)
    :   meanRatio_(meanRatio),
        varRatio_(varRatio){

    }

    bool operator()(
        const PixelType & meanA, const PixelType & varA,
        const PixelType & meanB, const PixelType & varB
    )const{
        // Compute mean ratio of mean and variance
        const ValueType m = mean(meanA/meanB);
        const ValueType v = mean(varA/varB);          
        return (m > meanRatio_ && m < (1.0 / meanRatio_) && v > varRatio_ && v < (1.0 / varRatio_));
    }

private:
    ValueType meanRatio_;
    ValueType varRatio_;
};


template<class PIXEL_TYPE_IN>
class SelectByNorm{
public:
    typedef typename NumericTraits<PIXEL_TYPE_IN>::RealPromote PixelType;
    typedef typename NumericTraits<PIXEL_TYPE_IN>::ValueType   ValueType;


    SelectByNorm(const ValueType meanNorm,  const ValueType varNorm)
    :   meanNorm_(meanNorm),
        varNorm_(varNorm){

    }

    bool operator()(
        const PixelType & meanA, const PixelType & varA,
        const PixelType & meanB, const PixelType & varB
    )const{
        // Compute mean ratio of mean and variance
        const ValueType m = norm(meanA-meanB);
        const ValueType v = norm(varA-varB);
        return (m<meanNorm_ && v<varNorm_);
    }
    
private:
    ValueType meanNorm_;
    ValueType varNorm_;
};





template<int DIM, class PIXEL_TYPE_IN>
class BlockWiseNonLocalMeanThreadObject{
    typedef PIXEL_TYPE_IN       PixelTypeIn;
    typedef typename NumericTraits<PixelTypeIn>::RealPromote        RealPromotePixelType;  
    typedef typename NumericTraits<RealPromotePixelType>::ValueType RealPromoteScalarType;
    
    typedef typename MultiArray<DIM,int>::difference_type Coordinate;
    typedef NonLocalMeanParameter ParameterType;

public:
    typedef void result_type;
    typedef UInt32              LabelType;
    typedef MultiArrayView<DIM,PixelTypeIn>          InArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType> MeanArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType> VarArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType> EstimateArrayView;
    typedef MultiArrayView<DIM,UInt32>               LabelArrayView;
    typedef std::vector<RealPromotePixelType>               BlockAverageVectorType;
    // range type
    typedef TinyVector<int,2> RangeType;
    //typedef boost::mutex      MutexType;

    typedef threading::thread  ThreadType;
    typedef threading::mutex   MutexType;

    BlockWiseNonLocalMeanThreadObject(
        const InArrayView &         inImage,
        MeanArrayView &             meanImage,
        VarArrayView &              varImage,
        EstimateArrayView &         estimageImage,
        LabelArrayView &            labelImage,
        const ParameterType &       param,
        const size_t                nThreads,
        MutexType &                 estimateMutex,
        MultiArray<1,int> &         progress
    )
    : 
    inImage_(inImage),
    meanImage_(meanImage),
    varImage_(varImage),
    estimageImage_(estimageImage),
    labelImage_(labelImage),
    param_(param),
    lastAxisRange_(),
    threadIndex_(),
    nThreads_(nThreads),
    estimateMutexPtr_(&estimateMutex),
    progress_(progress),
    average_(std::pow(2*param.patchRadius_+1,DIM)),
    shape_(inImage.shape()),
    totalSize_()
    {
        totalSize_ = 1;
        for(int dim=0;dim<DIM;++dim)
            totalSize_*=(shape_[dim]/param.stepSize_);
    }

    void setRange(const RangeType & lastAxisRange){
        lastAxisRange_=lastAxisRange;
    }
    void setThreadIndex(const size_t threadIndex){
        threadIndex_=threadIndex;
    }


    void operator()();

private:
    void processSinglePixel(const Coordinate & xyz);
    void processSinglePair( const Coordinate & xyz,const Coordinate & nxyz,RealPromoteScalarType & wmax,RealPromoteScalarType & totalweight);
    RealPromoteScalarType patchDistance(const Coordinate & xyz,const Coordinate & nxyz);

    void patchExtractAndAcc(const Coordinate & xyz,const RealPromoteScalarType weight);
    void patchAccMeanToEstimate(const Coordinate & xyz,const RealPromoteScalarType globalSum);

    void progressPrinter(const int counter);

    void mirrorIfIsOutsidePoint(Coordinate & coord)const{
        for(int c=0;c<DIM;++c){
            if(coord[c]<0)
                coord[c]=-1*coord[c];
            else if(coord[c]>= inImage_.shape(c))
                coord[c] = 2 * inImage_.shape(c) - coord[c] - 1;
        }
    }


    // array views
    InArrayView         inImage_;
    MeanArrayView       meanImage_;
    VarArrayView        varImage_;
    EstimateArrayView   estimageImage_;
    LabelArrayView      labelImage_;

    // param obj.
    ParameterType param_;

    // thread related; 
    RangeType lastAxisRange_;
    size_t threadIndex_;
    size_t nThreads_;
    MutexType * estimateMutexPtr_;
    MultiArrayView<1,int>  progress_;


    // computations
    BlockAverageVectorType average_;
    Coordinate shape_;
    size_t totalSize_;
};

template<int DIM,class PIXEL_TYPE_IN>
inline void BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::progressPrinter(const int counter){
    progress_[threadIndex_] = counter;
    if(threadIndex_==nThreads_-1){
        if(counter%10 == 0){
            int c=0;
            for(size_t ti=0;ti<nThreads_;++ti)
                c+=progress_[ti];
            double pr=c;
            pr/=totalSize_;
            pr*=100.0;
            std::cout<<"\rprogress "<<std::setw(10)<<pr<<" %%"<<std::flush;
        }
    }
}

template<int DIM,class PIXEL_TYPE_IN>
void BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::operator()(){
    const int start = lastAxisRange_[0];
    const int end = lastAxisRange_[1];
    const int stepSize = param_.stepSize_;

    Coordinate xyz; 
    int c=0;
    if(param_.verbose_ && threadIndex_==nThreads_-1){
        std::cout<<"progress"; 
    }

    if(DIM==2){
        for (xyz[1] = start; xyz[1]  < end;    xyz[1]  += stepSize)
        for (xyz[0] = 0;      xyz[0]  < shape_[0]; xyz[0]  += stepSize){
            this->processSinglePixel(xyz);
            if(param_.verbose_)
                this->progressPrinter(c);
            ++c;
        }
    }
    if(DIM==3){
        for (xyz[2] = start; xyz[2]  < end;    xyz[2]  += stepSize)
        for (xyz[1] = 0;   xyz[1]  < shape_[1]; xyz[1]  += stepSize)
        for (xyz[0] = 0;   xyz[0]  < shape_[0]; xyz[0]  += stepSize){
            this->processSinglePixel(xyz);
            if(param_.verbose_)
                this->progressPrinter(c);
            ++c;
        }
    }
    if(DIM==4){
        for (xyz[3] = start; xyz[3]  < end;    xyz[3]  += stepSize)
        for (xyz[2] = 0;   xyz[2]  < shape_[2]; xyz[2]  += stepSize)
        for (xyz[1] = 0;   xyz[1]  < shape_[1]; xyz[1]  += stepSize)
        for (xyz[0] = 0;   xyz[0]  < shape_[0]; xyz[0]  += stepSize){
            this->processSinglePixel(xyz);
            if(param_.verbose_)
                this->progressPrinter(c);
            ++c;
        }
    }
    if(param_.verbose_ && threadIndex_==nThreads_-1){
        std::cout<<"\rprogress "<<std::setw(10)<<"100"<<" %%"<<"\n";
    }
}

template<int DIM,class PIXEL_TYPE_IN>
inline void BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::processSinglePixel(
    const Coordinate & xyz
){
        Coordinate nxyz,xxyyzz;
        const int searchRadius = param_.searchRadius_;
        std::fill(average_.begin(),average_.end(),RealPromotePixelType(0.0));
        RealPromoteScalarType totalweight = 0.0;

        if( sum(meanImage_[xyz]) > param_.epsilon_  && sum(varImage_[xyz]) > param_.epsilon_){
            RealPromoteScalarType wmax = 0.0;

            if(DIM==2){
                for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                    nxyz = xyz  + xxyyzz;
                    if(isZero(xxyyzz))
                        continue;
                    this->processSinglePair(xyz,nxyz,wmax,totalweight);
                }
            }
            else if(DIM==3){
                for (xxyyzz[2] = -searchRadius; xxyyzz[2] <= searchRadius; xxyyzz[2]++)
                for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                    nxyz = xyz  + xxyyzz;
                    if(isZero(xxyyzz))
                        continue;
                    this->processSinglePair(xyz,nxyz,wmax,totalweight);
                }
            }
            else if(DIM==4){
                for (xxyyzz[3] = -searchRadius; xxyyzz[3] <= searchRadius; xxyyzz[3]++)
                for (xxyyzz[2] = -searchRadius; xxyyzz[2] <= searchRadius; xxyyzz[2]++)
                for (xxyyzz[1] = -searchRadius; xxyyzz[1] <= searchRadius; xxyyzz[1]++)
                for (xxyyzz[0] = -searchRadius; xxyyzz[0] <= searchRadius; xxyyzz[0]++){
                    nxyz = xyz  + xxyyzz;
                    if(isZero(xxyyzz))
                        continue;
                    this->processSinglePair(xyz,nxyz,wmax,totalweight);
                }
            }


            if (wmax == 0.0){
                wmax = 1.0;
            }
            // give pixel xyz  as much weight as
            // the maximum weighted other patch
            this->patchExtractAndAcc(xyz,wmax);
            totalweight += wmax;

            // this if seems total useless to me
            if (totalweight != 0.0){
                this->patchAccMeanToEstimate(xyz,totalweight);
            }

        }
        else{
            const double wmax = 1.0;
            this->patchExtractAndAcc(xyz,wmax);
            totalweight += wmax;
            this->patchAccMeanToEstimate(xyz,totalweight);
        }
}

template<int DIM,class PIXEL_TYPE_IN>
inline void BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::processSinglePair(
    const Coordinate & xyz,
    const Coordinate & nxyz,
    RealPromoteScalarType & wmax,
    RealPromoteScalarType & totalweight
){
    if( inImage_.isInside(nxyz)){
        if( sum(meanImage_[nxyz]) > param_.epsilon_  && sum(varImage_[nxyz]) > param_.epsilon_){

            // Compute Ratios // TODO fix types
            const double mRatio = mean(meanImage_[xyz]/meanImage_[nxyz]);
            const double vRatio = mean(varImage_[xyz]/varImage_[nxyz]);

            // here we check if to patches fit to each other
            // one patch is around xyz 
            // other patch is arround nxyz                
            if (mRatio > param_.meanRatio_ && mRatio < (1.0 / param_.meanRatio_) && vRatio > param_.varRatio_ && vRatio < (1.0 / param_.varRatio_)){
                const double d =this->patchDistance(xyz,nxyz);
                const RealPromoteScalarType hh = param_.sigma_ * param_.sigma_; // store hh in param?
                const RealPromoteScalarType w = exp(-d / (hh));
                wmax = std::max(w,wmax);
                this->patchExtractAndAcc(nxyz,w);
                totalweight+=  w;
            }
        }
    }
}


template<int DIM,class PIXEL_TYPE_IN>
inline typename BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::RealPromoteScalarType 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::patchDistance(
    const Coordinate & pA,
    const Coordinate & pB
){

    // TODO : use a acculator like think to make this more beautiful ?
    const int f = param_.patchRadius_;
    Coordinate offset,nPa,nPb;
    RealPromoteScalarType acu = 0;
    RealPromoteScalarType distancetotal = 0;
    
    #define VIGRA_NLM_IN_LOOP_CODE                              \
        nPa = pA+offset;                                        \
        nPb = pB+offset;                                        \
        this->mirrorIfIsOutsidePoint(nPa);                      \
        this->mirrorIfIsOutsidePoint(nPb);                      \
        const RealPromotePixelType vA = inImage_[nPa];          \
        const RealPromotePixelType vB = inImage_[nPb];          \
        distancetotal += vigra::sizeDividedSquaredNorm(vA-vB);  \
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

template<int DIM,class PIXEL_TYPE_IN>
inline void 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::patchExtractAndAcc(
    const Coordinate & xyz,
    const RealPromoteScalarType weight
){
    Coordinate xyzPos,abc;
    Coordinate nhSize3(param_.patchRadius_);
    const int ns = 2 * param_.patchRadius_ + 1;
    int count = 0;

    // todo: remove abc vector

    #define VIGRA_NLM_IN_LOOP_CODE                                              \
        xyzPos = xyz + abc - nhSize3;                                           \
        if (!param_.gaussNoise_){                                               \
            if (inImage_.isOutside(xyzPos))                                     \
                average_[count] += inImage_[xyz]*inImage_[xyz]* weight;         \
            else                                                                \
                average_[count] += inImage_[xyzPos]*inImage_[xyzPos]* weight;   \
        }                                                                       \
        else{                                                                   \
            if (inImage_.isOutside(xyzPos))                                     \
                average_[count] += inImage_[xyz]* weight;                       \
            else                                                                \
                average_[count] += inImage_[xyzPos]* weight;                    \
        }                                                                       \
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

template<int DIM,class PIXEL_TYPE_IN>
inline void 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN>::patchAccMeanToEstimate(
    const Coordinate & xyz,
    const RealPromoteScalarType globalSum
){
    Coordinate abc,xyzPos,nhSize(param_.patchRadius_);        
    int count = 0 ;
    const int ns = 2 * param_.patchRadius_ + 1;

    #define VIGRA_NLM_IN_LOOP_CODE                                                      \
            xyzPos = xyz + abc - nhSize;                                                \
            if ( inImage_.isInside(xyzPos)){                                            \
                if(!param_.gaussNoise_){                                                \
                    const double bias = 2 * param_.sigma_ * param_.sigma_;              \
                    estimateMutexPtr_->lock();                                          \
                    RealPromotePixelType value = estimageImage_[xyzPos];                \
                    RealPromotePixelType denoisedValue =                                \
                        (average_[count] / globalSum) -  RealPromotePixelType(bias);    \
                    denoisedValue=clipLower(denoisedValue);                             \
                    denoisedValue=vigra::sqrt(denoisedValue);                           \
                    value += denoisedValue;                                             \
                    estimageImage_[xyzPos] = value;                                     \
                    ++labelImage_[xyzPos];                                              \
                    estimateMutexPtr_->unlock();                                        \
                }                                                                       \
                else{                                                                   \
                    estimateMutexPtr_->lock();                                          \
                    RealPromotePixelType value = estimageImage_[xyzPos];                \
                    value += average_[count] / globalSum;                               \
                    estimageImage_[xyzPos] = value;                                     \
                    ++labelImage_[xyzPos];                                              \
                    estimateMutexPtr_->unlock();                                        \
                }                                                                       \
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


template<int DIM,class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT>
inline void gaussianMeanAndVariance(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & inArray,
    const double sigma,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & meanArray,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & varArray,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & tmpArray
){

    // compute mean and variance 
    vigra::gaussianSmoothMultiArray(inArray, meanArray, sigma);
    // square raw data (use estimate Image to store temp. results)
    for(int scanOrderIndex=0;scanOrderIndex<inArray.size();++scanOrderIndex){
        tmpArray[scanOrderIndex]=vigra::pow(inArray[scanOrderIndex],2);
    }
    vigra::gaussianSmoothMultiArray(tmpArray,varArray, sigma);
    for(int scanOrderIndex=0;scanOrderIndex<inArray.size();++scanOrderIndex){
        PIXEL_TYPE_OUT var = varArray[scanOrderIndex] - vigra::pow(meanArray[scanOrderIndex],2);
        var  = clipLower(var);
        //makeNegtiveValuesZero(var);  // callbyref
        varArray[scanOrderIndex] = var;
    }      
}

template<int DIM,class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT>
inline void gaussianMeanAndVariance(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & inArray,
    const double sigma,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & meanArray,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & varArray
){
    vigra::MultiArray<DIM,PIXEL_TYPE_OUT>  tmpArray(inArray.shape());
    gaussianMeanAndVariance<DIM,PIXEL_TYPE_IN,PIXEL_TYPE_OUT>(inArray,sigma,meanArray,varArray,tmpArray);   
}




template<int DIM, class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT>
void nonLocalMean(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & image,
    const NonLocalMeanParameter param,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> outImage
){

    typedef PIXEL_TYPE_IN       PixelTypeIn;
    typedef typename vigra::NumericTraits<PixelTypeIn>::RealPromote        RealPromotePixelType;  
    //typedef typename vigra::NumericTraits<RealPromotePixelType>::ValueType RealPromoteScalarType;
    //typedef PIXEL_TYPE_OUT      PixelTypeOut;
    typedef BlockWiseNonLocalMeanThreadObject<DIM,PixelTypeIn> ThreadObjectType;
    typedef typename ThreadObjectType::LabelType LabelType;

    // inspect parameter
    vigra_precondition(param.stepSize_>=1,"NonLocalMean Parameter: \"stepSize>=1\" violated");
    vigra_precondition(param.searchRadius_>=1, "NonLocalMean Parameter: \"searchRadius >=1\" violated");
    vigra_precondition(param.patchRadius_>=1,"NonLocalMean Parameter: \"searchRadius >=1\" violated");
    vigra_precondition(param.stepSize_-1<=param.patchRadius_,"NonLocalMean Parameter: \"stepSize -1 <= patchRadius\"  violated");

    // allocate arrays
    vigra::MultiArray<DIM,RealPromotePixelType> meanImage(image.shape());
    vigra::MultiArray<DIM,RealPromotePixelType> varImage(image.shape());
    vigra::MultiArray<DIM,RealPromotePixelType> estimageImage(image.shape());
    vigra::MultiArray<DIM,LabelType> labelImage(image.shape());

    // compute mean and variance 
    // last argument is a "buffer" since within "gaussianMeanAndVariance" another array is needed
    // ==> to avoid an unnecessary allocation we use the estimageImage as a buffer
    //gaussianMeanAndVariance<DIM,PixelTypeIn,RealPromotePixelType>(image,param.sigmaMean_,meanImage,varImage,estimageImage);
    gaussianMeanAndVariance<DIM,PixelTypeIn,RealPromotePixelType>(image,param.sigmaMean_,meanImage,varImage);

    // initialize
    labelImage = 0;
    estimageImage = RealPromotePixelType(0.0);

    ///////////////////////////////////////////////////////////////
    {   // MULTI THREAD CODE STARTS HERE



        typedef threading::thread  ThreadType;
        typedef threading::mutex   MutexType;

        MutexType estimateMutex;
        //typedef boost::thread ThreadType;

        const size_t nThreads =  param.nThreads_;
        MultiArray<1,int> progress = MultiArray<1,int>(typename  MultiArray<1,int>::difference_type(nThreads));

        // allocate all thread objects
        // each thread object works on a portion of the data
        std::vector<ThreadObjectType> threadObjects(nThreads, 
            ThreadObjectType(image, meanImage, varImage, estimageImage, labelImage, 
                param, nThreads, estimateMutex,progress)
        );

        // thread ptr
        std::vector<ThreadType *> threadPtrs(nThreads);
        for(size_t i=0; i<nThreads; ++i){
            ThreadObjectType & threadObj = threadObjects[i];
            threadObj.setThreadIndex(i);
            typename ThreadObjectType::RangeType lastAxisRange;
            lastAxisRange[0]=(i * image.shape(DIM-1)) / nThreads;
            lastAxisRange[1]=((i+1) * image.shape(DIM-1)) / nThreads;
            threadObj.setRange(lastAxisRange);
            // this will start the threads and cal operator() 
            // of the threadObjects
            threadPtrs[i] = new ThreadType(threadObjects[i]);
        }
        for(size_t i=0; i<nThreads; ++i)
            threadPtrs[i]->join();
        for(size_t i=0; i<nThreads; ++i)
            delete threadPtrs[i];

    }   // MULTI THREAD CODE ENDS HERE
    ///////////////////////////////////////////////////////////////

    // normalize estimates by the number of labels
    // and write that in output
    for(int scanOrderIndex=0; scanOrderIndex<labelImage.size(); ++scanOrderIndex){
        if (labelImage[scanOrderIndex] == 0)
            outImage[scanOrderIndex]=image[scanOrderIndex];
        else
            outImage[scanOrderIndex]=estimageImage[scanOrderIndex] / labelImage[scanOrderIndex];
    }
}


} // end namespace vigra


#endif