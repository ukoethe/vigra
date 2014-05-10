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
#include "gaussians.hxx"

namespace vigra{

struct NonLocalMeanParameter{

    NonLocalMeanParameter(
        const double sigmaSpatial = 2.0,
        const int searchRadius = 3,
        const int patchRadius = 1,
        const double sigmaMean = 1.0,
        const int stepSize = 2,
        const int iterations=1,
        const int nThreads = 8,
        const bool verbose = true
    ):
    sigmaSpatial_(sigmaSpatial),
    searchRadius_(searchRadius),
    patchRadius_(patchRadius),
    sigmaMean_(sigmaMean),
    stepSize_(stepSize),
    iterations_(iterations),
    nThreads_(nThreads),
    verbose_(verbose){
    }
    double sigmaSpatial_;
    int searchRadius_;
    int patchRadius_;
    double sigmaMean_;
    int stepSize_;
    int iterations_;
    int nThreads_;
    bool verbose_;
};


// this has no template since this should
// be the same class for different dimensions 
// and pixel types
class RatioPolicyParameter{
public:
    RatioPolicyParameter(
        const double sigma     = 5.0,
        const double meanRatio = 0.95,
        const double varRatio  = 0.5,
        const double epsilon   = 0.00001
    ):
    sigma_(sigma),
    meanRatio_(meanRatio),
    varRatio_(varRatio),
    epsilon_(epsilon){
    }
    double sigma_;
    double meanRatio_;
    double varRatio_;
    double epsilon_;
};


template<class PIXEL_TYPE_IN>
class RatioPolicy{

    public:
        typedef RatioPolicyParameter ParameterType;
        typedef typename NumericTraits<PIXEL_TYPE_IN>::RealPromote PixelType;
        typedef typename NumericTraits<PIXEL_TYPE_IN>::ValueType   ValueType;


        RatioPolicy(const ParameterType & param)
        :   meanRatio_(static_cast<ValueType>(param.meanRatio_)),
            varRatio_(static_cast<ValueType>(param.varRatio_)),
            epsilon_(static_cast<ValueType>(param.epsilon_)),
            sigmaSquared_(param.sigma_*param.sigma_){

        }

        bool usePixel(const PixelType & meanA, const PixelType & varA)const{
            return sum(meanA) > epsilon_  &&  sum(varA) > epsilon_;
        }


        bool usePixelPair(
            const PixelType & meanA, const PixelType & varA,
            const PixelType & meanB, const PixelType & varB
        )const{
            // Compute mean ratio of mean and variance
            const ValueType m = mean(meanA/meanB);
            const ValueType v = mean(varA/varB);          
            return (m > meanRatio_ && m < (1.0 / meanRatio_) && v > varRatio_ && v < (1.0 / varRatio_));
        }

        ValueType distanceToWeight(const PixelType & meanA, const PixelType & varA, const ValueType distance){
            return  exp(-distance /sigmaSquared_);
        }

    private:
        ValueType meanRatio_;
        ValueType varRatio_;
        ValueType epsilon_;
        ValueType sigmaSquared_;
};



// this has no template since this should
// be the same class for different dimensions 
// and pixel types
class NormPolicyParameter{
public:
    NormPolicyParameter(
        const double sigma     = 5.0,
        const double meanDist  = 0.95,
        const double varRatio   = 0.5,
        const double epsilon   = 0.00001
    ):
    sigma_(sigma),
    meanDist_(meanDist),
    varRatio_(varRatio),
    epsilon_(epsilon){
    }
    double sigma_;
    double meanDist_;
    double varRatio_;
    double epsilon_;
};



template<class V,int SIZE>
inline bool equal(const TinyVector<V,SIZE> & a,const TinyVector<V,SIZE> b){
    for(int i=0;i<SIZE;++i)
        if(a[i]!=b[i])
            return false;
    return true;
}



template<int DIM, bool ALWAYS_INSIDE>
struct BorderHelper;

template< int DIM>
struct BorderHelper<DIM,true>{
    template<class COORD,class IMAGE>
    static bool isInside(const COORD & ,const IMAGE &  ){
        return true;
    }
    template<class COORD,class IMAGE>
    static bool isOutside(const COORD & ,const IMAGE &  ){
        return false;
    }

    template<class COORD,class IMAGE>
    static void mirrorIfIsOutsidePoint(COORD & ,IMAGE & ){
    }

};

template< int DIM>
struct BorderHelper<DIM,false>{
    template<class COORD,class IMAGE>
    static bool isInside(const COORD & c,const IMAGE &  img){
        return img.isInside(c);
    }
    template<class COORD,class IMAGE>
    static bool isOutside(const COORD & c,const IMAGE &  img){
        return img.isOutside(c);
    }

    template<class COORD,class IMAGE>
    static void mirrorIfIsOutsidePoint(COORD & coord,const IMAGE & img){
        for(int c=0;c<DIM;++c){
            if(coord[c]<0)
                coord[c]=-1*coord[c];
            else if(coord[c]>= img.shape(c))
                coord[c] = 2 * img.shape(c) - coord[c] - 1;
        }
    }
};






template<class PIXEL_TYPE_IN>
class NormPolicy{

    public:
        typedef NormPolicyParameter ParameterType;
        typedef typename NumericTraits<PIXEL_TYPE_IN>::RealPromote PixelType;
        typedef typename NumericTraits<PIXEL_TYPE_IN>::ValueType   ValueType;


        NormPolicy(const ParameterType & param)
        :   meanDist_(static_cast<ValueType>(param.meanDist_)),
            varRatio_(static_cast<ValueType>(param.varRatio_)),
            epsilon_(static_cast<ValueType>(param.epsilon_)),
            sigmaSquared_(param.sigma_*param.sigma_){

        }

        bool usePixel(const PixelType & meanA, const PixelType & varA)const{
            return sum(varA)>epsilon_;
        }


        bool usePixelPair(
            const PixelType & meanA, const PixelType & varA,
            const PixelType & meanB, const PixelType & varB
        )const{
            // Compute mean ratio of mean and variance
            const ValueType m = squaredNorm(meanA-meanB);
            const ValueType v = mean(varA/varB); 
            //std::cout<<"norms  "<<m<<" v "<<v<<"\n";
            return (m < meanDist_ && v > varRatio_ && v < (1.0 / varRatio_));
        }

        ValueType distanceToWeight(const PixelType & meanA, const PixelType & varA, const ValueType distance){
            return  exp(-distance /sigmaSquared_);
        }

    private:
        ValueType meanDist_;
        ValueType varRatio_;
        ValueType epsilon_;
        ValueType sigmaSquared_;
};


// UnstridedArrayTag

template<int DIM, class PIXEL_TYPE_IN, class SMOOTH_POLICY>
class BlockWiseNonLocalMeanThreadObject{
    typedef PIXEL_TYPE_IN       PixelTypeIn;
    typedef typename NumericTraits<PixelTypeIn>::RealPromote        RealPromotePixelType;  
    typedef typename NumericTraits<RealPromotePixelType>::ValueType RealPromoteScalarType;
    
    typedef typename MultiArray<DIM,int>::difference_type Coordinate;
    typedef NonLocalMeanParameter ParameterType;

public:
    typedef void result_type;
    typedef MultiArrayView<DIM,PixelTypeIn>           InArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType>  MeanArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType>  VarArrayView;
    typedef MultiArrayView<DIM,RealPromotePixelType>  EstimateArrayView;
    typedef MultiArrayView<DIM,RealPromoteScalarType> LabelArrayView;
    typedef std::vector<RealPromotePixelType>         BlockAverageVectorType;
    typedef std::vector<RealPromoteScalarType>        BlockGaussWeightVectorType;
    typedef SMOOTH_POLICY                             SmoothPolicyType;
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
        const SmoothPolicyType  &   smoothPolicy,
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
    smoothPolicy_(smoothPolicy),
    param_(param),
    lastAxisRange_(),
    threadIndex_(),
    nThreads_(nThreads),
    estimateMutexPtr_(&estimateMutex),
    progress_(progress),
    average_(std::pow(2*param.patchRadius_+1,DIM)),
    gaussWeight_(std::pow(2*param.patchRadius_+1,DIM)),
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

    template<bool ALWAYS_INSIDE>
    void processSinglePixel(const Coordinate & xyz);

    template<bool ALWAYS_INSIDE>
    void processSinglePair( const Coordinate & xyz,const Coordinate & nxyz,RealPromoteScalarType & wmax,RealPromoteScalarType & totalweight);

    template<bool ALWAYS_INSIDE>
    RealPromoteScalarType patchDistance(const Coordinate & xyz,const Coordinate & nxyz);

    template<bool ALWAYS_INSIDE>
    void patchExtractAndAcc(const Coordinate & xyz,const RealPromoteScalarType weight);

    template<bool ALWAYS_INSIDE>
    void patchAccMeanToEstimate(const Coordinate & xyz,const RealPromoteScalarType globalSum);


    bool isAlwaysInside(const Coordinate & coord)const{
        const Coordinate r = (Coordinate(param_.searchRadius_) + Coordinate(param_.patchRadius_) +1 );
        const Coordinate test1 = coord - r;
        const Coordinate test2 = coord + r;
        return inImage_.isInside(test1) && inImage_.isInside(test2);
    }

    void initalizeGauss();

    void progressPrinter(const int counter);


    // array views
    InArrayView         inImage_;
    MeanArrayView       meanImage_;
    VarArrayView        varImage_;
    EstimateArrayView   estimageImage_;
    LabelArrayView      labelImage_;

    // policy object 
    SmoothPolicyType smoothPolicy_;

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
    BlockGaussWeightVectorType gaussWeight_;
    Coordinate shape_;
    size_t totalSize_;
};


template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
inline void BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::initalizeGauss(){
    Coordinate xyz;
    const int pr = param_.patchRadius_;
    Gaussian<RealPromoteScalarType> gaussian(param_.sigmaSpatial_);
    int c=0;
    RealPromoteScalarType sum = RealPromoteScalarType(0.0);
    if(DIM==2){
        for (xyz[1] = -pr; xyz[1]  <=pr; ++xyz[1])
        for (xyz[0] = -pr; xyz[0]  <=pr; ++xyz[0],++c){
            const RealPromoteScalarType distance = norm(xyz);
            const RealPromoteScalarType w =gaussian(distance);
            sum+=w;
            gaussWeight_[c]=w;   
        }
    }
    if(DIM==3){
        for (xyz[2] = -pr; xyz[2]  <=pr; ++xyz[2])
        for (xyz[1] = -pr; xyz[1]  <=pr; ++xyz[1])
        for (xyz[0] = -pr; xyz[0]  <=pr; ++xyz[0],++c){
            const RealPromoteScalarType distance = norm(xyz);
            const RealPromoteScalarType w =gaussian(distance);
            sum+=w;
            gaussWeight_[c]=w;     
        }
    }
    if(DIM==4){
        for (xyz[3] = -pr; xyz[3]  <=pr; ++xyz[3])
        for (xyz[2] = -pr; xyz[2]  <=pr; ++xyz[2])
        for (xyz[1] = -pr; xyz[1]  <=pr; ++xyz[1])
        for (xyz[0] = -pr; xyz[0]  <=pr; ++xyz[0],++c){
            const RealPromoteScalarType distance = norm(xyz);
            const RealPromoteScalarType w =gaussian(distance);
            sum+=w;
            gaussWeight_[c]=w;   
        }
    }
    // normalize
    for(size_t i=0;i<gaussWeight_.size();++i){
        gaussWeight_[i]/=sum;
    }
}


template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
inline void BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::progressPrinter(const int counter){
    progress_[threadIndex_] = counter;
    if(threadIndex_==nThreads_-1){
        if(counter%100 == 0){
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

template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
void BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::operator()(){
    const int start = lastAxisRange_[0];
    const int end = lastAxisRange_[1];
    const int stepSize = param_.stepSize_;

    this->initalizeGauss();

    Coordinate xyz; 
    int c=0;
    if(param_.verbose_ && threadIndex_==nThreads_-1){
        std::cout<<"progress"; 
    }

    if(DIM==2){
        for (xyz[1] = start; xyz[1]  < end;    xyz[1]  += stepSize)
        for (xyz[0] = 0;      xyz[0]  < shape_[0]; xyz[0]  += stepSize){

            if(isAlwaysInside(xyz))
                this->processSinglePixel<true>(xyz);
            else
                this->processSinglePixel<false>(xyz);
            if(param_.verbose_)
                this->progressPrinter(c);
            ++c;
        }
    }
    if(DIM==3){
        for (xyz[2] = start; xyz[2]  < end;    xyz[2]  += stepSize)
        for (xyz[1] = 0;   xyz[1]  < shape_[1]; xyz[1]  += stepSize)
        for (xyz[0] = 0;   xyz[0]  < shape_[0]; xyz[0]  += stepSize){
            if(isAlwaysInside(xyz))
                this->processSinglePixel<true>(xyz);
            else
                this->processSinglePixel<false>(xyz);
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
            if(isAlwaysInside(xyz))
                this->processSinglePixel<true>(xyz);
            else
                this->processSinglePixel<false>(xyz);
            if(param_.verbose_)
                this->progressPrinter(c);
            ++c;
        }
    }
    if(param_.verbose_ && threadIndex_==nThreads_-1){
        std::cout<<"\rprogress "<<std::setw(10)<<"100"<<" %%"<<"\n";
    }
}



template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
template<bool ALWAYS_INSIDE>
inline void BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::processSinglePixel(
    const Coordinate & xyz
){
        Coordinate nxyz(SkipInitialization);
        const int searchRadius = param_.searchRadius_;
        std::fill(average_.begin(),average_.end(),RealPromotePixelType(0.0));
        RealPromoteScalarType totalweight = 0.0;

        if(smoothPolicy_.usePixel(meanImage_[xyz],varImage_[xyz])){

            RealPromoteScalarType wmax = 0.0;
            const Coordinate start = xyz- Coordinate(param_.searchRadius_);
            const Coordinate end   = xyz+ Coordinate(param_.searchRadius_);

            if(DIM==2){
                for (nxyz[1] = start[1]; nxyz[1] <= end[1]; nxyz[1]++)
                for (nxyz[0] = start[0]; nxyz[0] <= end[0]; nxyz[0]++){
                    //nxyz = xyz  + nxyz;
                    if(equal(nxyz,xyz))
                        continue;
                    this->processSinglePair<ALWAYS_INSIDE>(xyz,nxyz,wmax,totalweight);
                }
            }
            else if(DIM==3){
                for (nxyz[2] = start[2]; nxyz[2] <= end[2]; nxyz[2]++)
                for (nxyz[1] = start[1]; nxyz[1] <= end[1]; nxyz[1]++)
                for (nxyz[0] = start[0]; nxyz[0] <= end[0]; nxyz[0]++){
                    if(equal(nxyz,xyz))
                        continue;
                    this->processSinglePair<ALWAYS_INSIDE>(xyz,nxyz,wmax,totalweight);
                }
            }
            else if(DIM==4){
                for (nxyz[3] = start[3]; nxyz[3] <= end[3]; nxyz[3]++)
                for (nxyz[2] = start[2]; nxyz[2] <= end[2]; nxyz[2]++)
                for (nxyz[1] = start[1]; nxyz[1] <= end[1]; nxyz[1]++)
                for (nxyz[0] = start[0]; nxyz[0] <= end[0]; nxyz[0]++){
                    if(equal(nxyz,xyz))
                        continue;
                    this->processSinglePair<ALWAYS_INSIDE>(xyz,nxyz,wmax,totalweight);
                }
            }


            if (wmax == 0.0){
                wmax = 1.0;
            }
            // give pixel xyz  as much weight as
            // the maximum weighted other patch
            this->patchExtractAndAcc<ALWAYS_INSIDE>(xyz,wmax);
            totalweight += wmax;

            // this if seems total useless to me
            if (totalweight != 0.0){
                this->patchAccMeanToEstimate<ALWAYS_INSIDE>(xyz,totalweight);
            }

        }
        else{
            const double wmax = 1.0;
            this->patchExtractAndAcc<ALWAYS_INSIDE>(xyz,wmax);
            totalweight += wmax;
            this->patchAccMeanToEstimate<ALWAYS_INSIDE>(xyz,totalweight);
        }
}



template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
template<bool ALWAYS_INSIDE>
inline void BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::processSinglePair(
    const Coordinate & xyz,
    const Coordinate & nxyz,
    RealPromoteScalarType & wmax,
    RealPromoteScalarType & totalweight
){
    
    if(BorderHelper<DIM,ALWAYS_INSIDE>::isInside(nxyz,inImage_)){
    //if(ALWAYS_INSIDE || inImage_.isInside(nxyz)){
        if(smoothPolicy_.usePixel(meanImage_[nxyz],varImage_[nxyz])){
            // here we check if to patches fit to each other
            // one patch is around xyz 
            // other patch is arround nxyz
            if(smoothPolicy_.usePixelPair(meanImage_[xyz],varImage_[xyz],meanImage_[nxyz],varImage_[nxyz])){                
                const RealPromoteScalarType distance =this->patchDistance<ALWAYS_INSIDE>(xyz,nxyz);
                const RealPromoteScalarType w = smoothPolicy_.distanceToWeight(meanImage_[xyz],varImage_[xyz],distance);
                wmax = std::max(w,wmax);
                this->patchExtractAndAcc<ALWAYS_INSIDE>(nxyz,w);
                totalweight+=  w;
            }
        }
    }
}



template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
template<bool ALWAYS_INSIDE>
inline typename BlockWiseNonLocalMeanThreadObject<DIM, PIXEL_TYPE_IN, SMOOTH_POLICY>::RealPromoteScalarType 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN,SMOOTH_POLICY>::patchDistance(
    const Coordinate & pA,
    const Coordinate & pB
){

    // TODO : use a acculator like think to make this more beautiful ?
    const int f = param_.patchRadius_;
    Coordinate offset(SkipInitialization),nPa(SkipInitialization),nPb(SkipInitialization);
    int acu = 0;
    RealPromoteScalarType distancetotal = 0;
    int c =0 ;
    //this->mirrorIfIsOutsidePoint<ALWAYS_INSIDE>(nPa);       
   // this->mirrorIfIsOutsidePoint<ALWAYS_INSIDE>(nPb);       
    #define VIGRA_NLM_IN_LOOP_CODE                              \
        nPa = pA+offset;                                        \
        nPb = pB+offset;                                        \
        BorderHelper<DIM,ALWAYS_INSIDE>::mirrorIfIsOutsidePoint(nPa,inImage_); \
        BorderHelper<DIM,ALWAYS_INSIDE>::mirrorIfIsOutsidePoint(nPb,inImage_); \
        const RealPromoteScalarType gaussWeight = gaussWeight_[c]; \
        const RealPromotePixelType vA = inImage_[nPa];          \
        const RealPromotePixelType vB = inImage_[nPb];          \
        distancetotal += gaussWeight*vigra::sizeDividedSquaredNorm(vA-vB);  \
        ++acu;

    if(DIM==2){
        for (offset[1] = -f; offset[1] <= f; ++offset[1])
        for (offset[0] = -f; offset[0] <= f; ++offset[0],++c){
            VIGRA_NLM_IN_LOOP_CODE;
        }
    }
    else if(DIM==3){
        for (offset[2] = -f; offset[2] <= f; ++offset[2])
        for (offset[1] = -f; offset[1] <= f; ++offset[1])
        for (offset[0] = -f; offset[0] <= f; ++offset[0],++c){
            VIGRA_NLM_IN_LOOP_CODE;
        }
    }
    else if(DIM==4){
        for (offset[3] = -f; offset[3] <= f; ++offset[3])
        for (offset[2] = -f; offset[2] <= f; ++offset[2])
        for (offset[1] = -f; offset[1] <= f; ++offset[1])
        for (offset[0] = -f; offset[0] <= f; ++offset[0],++c){
            VIGRA_NLM_IN_LOOP_CODE;
        }
    }
    #undef VIGRA_NLM_IN_LOOP_CODE
    return distancetotal / acu;
}



template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
template<bool ALWAYS_INSIDE>
inline void 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN,SMOOTH_POLICY>::patchExtractAndAcc(
    const Coordinate & xyz,
    const RealPromoteScalarType weight
){
    Coordinate xyzPos(SkipInitialization),abc(SkipInitialization);
    Coordinate nhSize3(param_.patchRadius_);
    const int ns = 2 * param_.patchRadius_ + 1;
    int count = 0;

    // todo: remove abc vector

    #define VIGRA_NLM_IN_LOOP_CODE                                          \
        xyzPos = xyz + abc - nhSize3;                                       \
        if(BorderHelper<DIM,ALWAYS_INSIDE>::isOutside(xyzPos,inImage_))     \
            average_[count] += inImage_[xyz]* weight;                       \
        else                                                                \
            average_[count] += inImage_[xyzPos]* weight;                    \
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


template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY>
template<bool ALWAYS_INSIDE>
inline void 
BlockWiseNonLocalMeanThreadObject<DIM,PIXEL_TYPE_IN,SMOOTH_POLICY>::patchAccMeanToEstimate(
    const Coordinate & xyz,
    const RealPromoteScalarType globalSum
){
    Coordinate abc(SkipInitialization),xyzPos(SkipInitialization),nhSize(param_.patchRadius_);        
    int count = 0 ;
    const int ns = 2 * param_.patchRadius_ + 1;

    #define VIGRA_NLM_IN_LOOP_CODE                                              \
            xyzPos = xyz + abc - nhSize;                                        \
            if(BorderHelper<DIM,ALWAYS_INSIDE>::isInside(xyzPos,inImage_)){     \
                estimateMutexPtr_->lock();                                      \
                RealPromotePixelType value = estimageImage_[xyzPos];            \
                const RealPromoteScalarType gw = gaussWeight_[count];           \
                RealPromotePixelType tmp =(average_[count] / globalSum);        \
                tmp*=gw;                                                        \
                value +=tmp;                                                    \
                estimageImage_[xyzPos] = value;                                 \
                labelImage_[xyzPos]+=gw;                                        \
                estimateMutexPtr_->unlock();                                    \
            }                                                                   \
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



template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY,class PIXEL_TYPE_OUT>
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

template<int DIM,class PIXEL_TYPE_IN, class SMOOTH_POLICY,class PIXEL_TYPE_OUT>
inline void gaussianMeanAndVariance(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & inArray,
    const double sigma,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & meanArray,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> & varArray
){
    vigra::MultiArray<DIM,PIXEL_TYPE_OUT>  tmpArray(inArray.shape());
    gaussianMeanAndVariance<DIM,PIXEL_TYPE_IN,PIXEL_TYPE_OUT>(inArray,sigma,meanArray,varArray,tmpArray);   
}

namespace detail_non_local_means{

template<int DIM, class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT,class SMOOTH_POLICY>
void nonLocalMean1Run(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & image,
    const SMOOTH_POLICY & smoothPolicy,
    const NonLocalMeanParameter param,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> outImage
){

    typedef PIXEL_TYPE_IN       PixelTypeIn;
    typedef typename vigra::NumericTraits<PixelTypeIn>::RealPromote         RealPromotePixelType;  
     typedef typename vigra::NumericTraits<RealPromotePixelType>::ValueType RealPromoteScalarType;  
    typedef SMOOTH_POLICY SmoothPolicyType;

    typedef BlockWiseNonLocalMeanThreadObject<DIM,PixelTypeIn,SmoothPolicyType> ThreadObjectType;


    // make the policy object
    typedef typename SmoothPolicyType::ParameterType  PolicyParameterType;
    //const PolicyParameterType policyParam(param.sigma_, param.meanRatio_, param.varRatio_);
    //const SmoothPolicyType smoothPolicy(policyParam);

    // inspect parameter
    vigra_precondition(param.stepSize_>=1,"NonLocalMean Parameter: \"stepSize>=1\" violated");
    vigra_precondition(param.searchRadius_>=1, "NonLocalMean Parameter: \"searchRadius >=1\" violated");
    vigra_precondition(param.patchRadius_>=1,"NonLocalMean Parameter: \"searchRadius >=1\" violated");
    vigra_precondition(param.stepSize_-1<=param.patchRadius_,"NonLocalMean Parameter: \"stepSize -1 <= patchRadius\"  violated");

    // allocate arrays
    vigra::MultiArray<DIM,RealPromotePixelType> meanImage(image.shape());
    vigra::MultiArray<DIM,RealPromotePixelType> varImage(image.shape());
    vigra::MultiArray<DIM,RealPromotePixelType> estimageImage(image.shape());
    vigra::MultiArray<DIM,RealPromoteScalarType> labelImage(image.shape());

    // compute mean and variance 
    // last argument is a "buffer" since within "gaussianMeanAndVariance" another array is needed
    // ==> to avoid an unnecessary allocation we use the estimageImage as a buffer
    //gaussianMeanAndVariance<DIM,PixelTypeIn,RealPromotePixelType>(image,param.sigmaMean_,meanImage,varImage,estimageImage);
    gaussianMeanAndVariance<DIM,PixelTypeIn,RealPromotePixelType>(image,param.sigmaMean_,meanImage,varImage);

    // initialize
    labelImage = RealPromoteScalarType(0.0);
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
                smoothPolicy, param, nThreads, estimateMutex,progress)
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
        if (labelImage[scanOrderIndex] <= RealPromoteScalarType(0.00001))
            outImage[scanOrderIndex]=image[scanOrderIndex];
        else
            outImage[scanOrderIndex]=estimageImage[scanOrderIndex] / labelImage[scanOrderIndex];
    }
}

}

template<int DIM, class PIXEL_TYPE_IN,class PIXEL_TYPE_OUT,class SMOOTH_POLICY>
void nonLocalMean(
    const vigra::MultiArrayView<DIM,PIXEL_TYPE_IN> & image,
    const SMOOTH_POLICY & smoothPolicy,
    const NonLocalMeanParameter param,
    vigra::MultiArrayView<DIM,PIXEL_TYPE_OUT> outImage
){
    detail_non_local_means::nonLocalMean1Run<DIM,PIXEL_TYPE_IN,PIXEL_TYPE_OUT,SMOOTH_POLICY>(image,smoothPolicy,param,outImage);
    if(param.iterations_>1){

        vigra::MultiArray<DIM,PIXEL_TYPE_OUT> tmp(outImage.shape());        
        for(size_t i=0;i<param.iterations_-1;++i){
            tmp=outImage;
            detail_non_local_means::nonLocalMean1Run<DIM,PIXEL_TYPE_OUT,PIXEL_TYPE_OUT,SMOOTH_POLICY>(tmp,smoothPolicy,param,outImage);
        }
    }
}




} // end namespace vigra


#endif