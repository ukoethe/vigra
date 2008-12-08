/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>
#include <functional> 

//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraDistance

using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
    declScalarMinMax(int, backgroundMode, 1, 0, 2);
    declCharConstr(method, 4, MULT, MULT_SQUARED, IMAG_DIST_TRANS, INVERTEDCRAP, e);
    declOut(double);
    
    //Only supported with MULT
    TinyVector<double, 3> pitch3D;
    //Only supported with IMAG_DIST_TRANS
    declScalar(T, backgroundPixel, 0);
    declScalarMinMax(int, norm, 2, 0, 2);
    
    MultiArray<3, double> tmp3D;
    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    : base_data<T>(inputs), map(backgroundMode), map(method), map(backgroundPixel), map(norm),
      pitch3D(1.0,1.0,1.0), tmp3D(this->in3D.shape())
    {
        mapOut_SAME(double);
        if(inputs.size() == 2)
        {
            mxArray* pitchArr =mxGetField(inputs[1], 0, "pitch");
            if(pitchArr != NULL && mxIsNumeric(pitchArr))
            {
                if(this->numOfDim == IMAG){
                    TinyVectorView<double, 2> temp = matlab::getVector<double,2>(pitchArr);
                    pitch3D[0] = temp[0]>0? temp[0]:1.0;
                    pitch3D[1] = temp[1]>0? temp[1]:1.0;
                }else{
                    TinyVectorView<double, 3> temp = matlab::getVector<double,3>(pitchArr);
                    pitch3D[0] = temp[0]>0? temp[0]:1.0;
                    pitch3D[1] = temp[1]>0? temp[1]:1.0;
                    pitch3D[2] = temp[2]>0? temp[2]:1.0;
                }
                if(this->method != MULT){
                    mexWarnMsgTxt("pitch option is only supported with the method MULT");
                }
            }
        }
        if(this->method == IMAG_DIST_TRANS){
            backgroundMode = 0;
            if(this->numOfDim == VOLUME){
                this->method = MULT;
                mexWarnMsgTxt("IMAG_DIST_TRANS only works with 2D Images using default: MULT");
            }
        }

    }
};

    struct signChangeFunctor
    {
        template <class PixelType>
        PixelType operator()(PixelType const& v) const
        {
            return 0.0-v;
        }
    };

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP2_(a, b) cP<a, data<T>::b>::value
struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
        switch(cantorPair(o.backgroundMode == 2, o.method)){
            //In this case function pointers may have been more elegant.
            case cP2_(0, MULT):
                    separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode == 1, o.pitch3D);
                    break;
            case cP2_(1, MULT):
                separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 0, o.pitch3D);
                transformMultiArray(srcMultiArrayRange(o.out3D), destMultiArray(o.out3D), signChangeFunctor());
                separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.tmp3D), 1, o.pitch3D);
                combineTwoMultiArrays(
                            srcMultiArrayRange(o.tmp3D), 
                            srcMultiArray(o.out3D), 
                            destMultiArray(o.out3D),  
                            std::plus<double>());
                break;
            case cP2_(0, MULT_SQUARED):
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode == 1);
                break;
            case cP2_(1, MULT_SQUARED):
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 0);
                transformMultiArray(srcMultiArrayRange(o.out3D), destMultiArray(o.out3D), signChangeFunctor());
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.tmp3D), 1);
                combineTwoMultiArrays(
                            srcMultiArrayRange(o.tmp3D), 
                            srcMultiArray(o.out3D), 
                            destMultiArray(o.out3D),  
                            std::plus<double>());
                break;
            case cP2_(0, IMAG_DIST_TRANS):
                distanceTransform(srcImageRange(o.in), destImage(o.out),o.backgroundPixel, o.norm);
                break;
            default:
                mexErrMsgTxt("Precondition checking not complete - something went wrong");
        }
        
    }
};


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the 
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB 
function D = vigraDistance(inputArray)
function D = vigraDistance(inputArray, options);

D = vigraDistance(inputArray) computes the Distancetransform using the default options.
D = vigraDistance(inputImage, options)  does the same with user options.
options is a struct with possible fields: "method", "backgroundMode" and "backgroundPixel" and "norm"

"method":                 'MULT'(default), 'MULT_SQUARED', 'IMAG_DIST_TRANS'
                        MULT and MULT_SQUARED are the faster and newer distance transform methods defined with VIGRA-Multiarrays
                        (vigra::seperableMultiDist....).MULT_SQUARED returns the squared values of MULT. 
                        IMAG_DIST_TRANS is defined with BasicImage and is slower. Use it only if  you explicitely need the backgroundPixel
                        or norm option. (vigra::distanceTransform)
"backgroundMode":         0 ,  1(default) , 2:     
                        This option is only used with methods MULT and MULT_SQUARED
                        If the parameter background is 1, then the squared distance of all background pixels to the nearest object is calculated. 
                        Otherwise, the distance of all object pixels to the nearest background pixel is calculated.
                        If the parameter background is 2, then the distance if background pixels to the nearest object will be negative valued and
                        the distance of object pixels to the background pixels shall be positive
"backgroundPixel:         0(default) , arb. value in grayscale range:         
                        This option defines the background Pixel value. Only used with method = IMAG_DIST_TRANS
"norm":                     2(default), 1 , 0.  
                        Defines the norm used to calculate the distance (Only used with method = IMAG_DIST_TRANS)
"pitch"                [1.0, 1.0]  (2D)or [1.0, 1.0, 1.0](3D) , arb 2D or 3D array.
                        define a pitch if data has non-uniform resolution .

Usage:
    opt = struct('method' ,'IMAGE_DIST_TRANS' , 'backgroundPixel', 10 , 'norm' , 0);
    out = vigraDistance(in, opt);

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
    // 
    callMexFunctor<vigraFunctor>(outputs, inputs);
}
