/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>
#include <functional> 

//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data
: public base_data<T>
{
    declScalarMinMax(int, backgroundMode, 1, 0, 2);
    declCharConstr3(method, MULT, MULT_SQUARED, IMAG_DIST_TRANS);
    declOut(double);
    
    //Only supported with MULT
    TinyVector<double, 3> pitch3D;
    //Only supported with IMAG_DIST_TRANS
    declScalar(T, backgroundValue, 0);
    declScalarMinMax(int, norm, 2, 0, 2);
    
    MultiArray<3, double> tmp3D;
    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    : 	base_data<T>(inputs), 
		initOption(backgroundMode), 
		initOption(method), 
		initOption(backgroundValue), 
		initOption(norm),
      pitch3D(1.0,1.0,1.0), tmp3D(base_data<T>::in3D.shape())
    {
        initOut_SAME(double);
        if(this->options.isValid("pitch"))
        {
            if(this->method == IMAG_DIST_TRANS)
                mexErrMsgTxt("vigraDistance(): 'pitch' option not supported by method 'IMAG_DIST_TRANS'");
            if(this->numOfDim == IMAG)
            {
                TinyVectorView<double, 2> temp = matlab::getTinyVector<2,double>(this->options["pitch"]);
                pitch3D[0] = temp[0]>0? temp[0]:1.0;
                pitch3D[1] = temp[1]>0? temp[1]:1.0;
            }
            else
            {
                TinyVectorView<double, 3> temp = matlab::getTinyVector<3,double>(this->options["pitch"]);
                pitch3D[0] = temp[0]>0? temp[0]:1.0;
                pitch3D[1] = temp[1]>0? temp[1]:1.0;
                pitch3D[2] = temp[2]>0? temp[2]:1.0;
            }
        }
        if(this->method == IMAG_DIST_TRANS)
        {
            if(this->numOfDim == VOLUME)
                mexErrMsgTxt("vigraDistance(): method 'IMAG_DIST_TRANS' requires 2D data.");
            if(backgroundMode != 1)
                mexErrMsgTxt("vigraDistance(): method 'IMAG_DIST_TRANS' requires 'backgroundMode' = 1.");
        }
        else
        {
            if(backgroundValue != 0)
                mexErrMsgTxt("vigraDistance(): methods 'MULT' and 'MULT_SQUARED' require 'backgroundValue' = 0.");
            if(norm != 2)
                mexErrMsgTxt("vigraDistance(): methods 'MULT' and 'MULT_SQUARED' require 'norm' = 2.");
        }
    }
};

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP2_(a, b) cP<a, data<T>::b>::value
struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs)
    {
        using namespace vigra::functor;
        //Options
        data<T>  o(outputs, inputs);
        
        bool computeSignedDist = (o.backgroundMode == 2);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
        switch(cantorPair(computeSignedDist, o.method))
        {
            //In this case function pointers may have been more elegant.
            case cP2_(false, MULT):
                    separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode == 1, o.pitch3D);
                    break;
            case cP2_(true, MULT):
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 0, o.pitch3D);
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.tmp3D), 1, o.pitch3D);
                combineTwoMultiArrays(
                            srcMultiArrayRange(o.tmp3D), 
                            srcMultiArray(o.out3D), 
                            destMultiArray(o.out3D),  
                            ifThenElse(Arg1() > Param(0.0), sqrt(Arg1())-Param(0.5), Param(0.5)-sqrt(Arg2())));
                break;
            case cP2_(0, MULT_SQUARED):
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode == 1);
                break;
            case cP2_(1, MULT_SQUARED):
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 0);
                separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.tmp3D), 1);
                combineTwoMultiArrays(
                            srcMultiArrayRange(o.tmp3D), 
                            srcMultiArray(o.out3D), 
                            destMultiArray(o.out3D),  
                            ifThenElse(Arg1() > Param(0.0), sq(sqrt(Arg1())-Param(0.5)), -sq(sqrt(Arg2())-Param(0.5))));
                break;
            case cP2_(0, IMAG_DIST_TRANS):
                distanceTransform(srcImageRange(o.in), destImage(o.out),o.backgroundValue, o.norm);
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

D = vigraDistance(inputArray) computes the distance transform using the default options.
D = vigraDistance(inputImage, options)  does the same with user options.

inputArray  - a 2D or 3D array of numeric type
options     - a struct with the following possible fields (default will be used 
              if field is not present)
    'method':          'MULT'(default), 'MULT_SQUARED', 'IMAG_DIST_TRANS'
                       MULT and MULT_SQUARED are the faster and newer distance transform 
                       methods defined with VIGRA-Multiarrays (vigra::seperableMultiDist....).
                       MULT_SQUARED returns the squared values of MULT. 
                       IMAG_DIST_TRANS is defined with BasicImage (vigra::distanceTransform) and
                       less accurate. Use it only if you explicitely need the 'backgroundValue'
                       or 'norm' options.
    'backgroundValue': 0 (default) , arbitrary value (only supported by IMAG_DIST_TRANS)
                       This option defines the background value. In MULT and MULT_SQUARED, the
                       'backgroundValue' is always 0, but see option 'backgroundMode'.
    'backgroundMode':  0 , 1 (default) , 2:     
                       This option is only used with methods MULT and MULT_SQUARED.
                       In method IMAG_DIST_TRANS, the distance of background points
                         (according to 'backgroundValue' above) to the nearest 
                         non-background is computed.
                       If 'backgroundMode' is 1, then the (squared) distance of all background 
                         points to the nearest object is calculated. 
                       If 'backgroundMode' is 0, the (squared) distance of all object 
                         points to the nearest background is calculated. 
                       If 'backgroundMode' is 2, the signed (squared) distance of all points
                         to the contour will be calculated, such that negative values are
                         inside the objects, positive ones in the background. IMAG_DIST_TRANS
    'norm':            2 (default, Euclidean distance), 1 (L1 distance), 0 (L-infinity distance).  
                       Defines the norm used to calculate the distance.
                       Only supported by method IMAG_DIST_TRANS
    'pitch':           2D: [1.0, 1.0] (default), arbitrary int32-Array of length 2.
                       3D: [1.0, 1.0, 1.0] (default), arbitrary int32-Array of length 3.
                       Define the pixel distance if data has non-uniform resolution.
                       Only supported by methods MULT and MULT_SQUARED.

Usage:
    opt = struct('method' ,'IMAGE_DIST_TRANS' , 'backgroundValue', 10 , 'norm' , 0);
    out = vigraDistance(in, opt);

*/

