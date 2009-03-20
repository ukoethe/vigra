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




//#define RN_DEBUG
#define cP2_(a, b) cP<(int)a, b>::value
template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){

    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    //Load input Image
    MultiArrayView<3,T>         in3D        = inputs.getMultiArray<3,T>(0, v_required());
    BasicImageView<T>           in          = makeBasicImageView(in3D.bindOuter(0));
    int                         numOfDim    = inputs.getDimOfInput(0, v_required());

    //Load Method Option
    VIGRA_CREATE_ENUM_AND_STD_MAP3(Methods,MapName, MULT, MULT_SQUARED, IMAG_DIST_TRANS);
    Methods         method      = (Methods)inputs.getEnum("method", v_default((int)MULT), MapName );


    //Load backgroundValue/Mode/norm
    int                 norm                = inputs.getScalarMinMax<int>("norm", v_default(2), 0, 2);
    T                   backgroundValue     = inputs.getScalar<T>("backgroundValue", v_default(0));
    int                 backgroundMode      = inputs.getScalarMinMax<int>("backgroundMode", v_default(1), 0, 2);

    //Load Pitch
    TinyVector<double, 3>       defaultPitch3D(1,1, 1);
    TinyVector<double, 2>       defaultPitch(1,1);
    TinyVectorView<double, 2>   pitch       =  inputs.getTinyVector<double, 2> ( "pitch", v_default(defaultPitch));
    TinyVectorView<double, 3>   pitch3D     = (numOfDim == 3)?
                                                inputs.getTinyVector<double, 3> ( "pitch", v_default(defaultPitch3D))
                                            :   vigra::TinyVectorView<double, 3>(defaultPitch3D);


    //This is a cheap way of checking whether pitch option has been set - if not the pointers of pitch and defaultPitch
    //should be the same;
    //Some more errorchecking
    if(method == IMAG_DIST_TRANS)
    {
        if  (numOfDim == VOLUME)
                mexErrMsgTxt("vigraDistance(): method 'IMAG_DIST_TRANS' requires 2D data.");
        if  (pitch.data() != defaultPitch.data())
                mexErrMsgTxt("vigraDistance(): 'IMAG_DIST_TRANS' does not support 'pitch' Option");
        if  (backgroundMode != 1)
                mexErrMsgTxt("vigraDistance(): method 'IMAG_DIST_TRANS' requires 'backgroundMode' = 1.");
    }
    else
    {
        if  (backgroundValue != 0)
                mexErrMsgTxt("vigraDistance(): methods 'MULT' and 'MULT_SQUARED' require 'backgroundValue' = 0.");
        if  (norm != 2)
                mexErrMsgTxt("vigraDistance(): methods 'MULT' and 'MULT_SQUARED' require 'norm' = 2.");
    }

    //Allocate Memory for output
    typedef double outType;
    MultiArrayView<3,outType>   out3D       = outputs.createMultiArray      <3,outType>   (0, v_required(), in3D.shape());
    BasicImageView<outType>     out(out3D.data(), in3D.shape(0), in3D.shape(1));

    MultiArray<3, outType>      tmp3D(in3D.shape());
    bool                        computeSignedDist  = (backgroundMode == 2);

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    // catorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
    using namespace vigra::functor;

    switch(cantorPair(computeSignedDist, method))
    {
        //In this case function pointers may have been more elegant.
        case cP2_(false, MULT):
                separableMultiDistance(srcMultiArrayRange(in3D), destMultiArray(out3D), backgroundMode == 1, pitch3D);
                break;
        case cP2_(true, MULT):
            separableMultiDistSquared(srcMultiArrayRange(in3D), destMultiArray(out3D), 0, pitch3D);
            separableMultiDistSquared(srcMultiArrayRange(in3D), destMultiArray(tmp3D), 1, pitch3D);
            combineTwoMultiArrays(
                        srcMultiArrayRange(tmp3D),
                        srcMultiArray(out3D),
                        destMultiArray(out3D),
                        ifThenElse(Arg1() > Param(0.0), sqrt(Arg1())-Param(0.5), Param(0.5)-sqrt(Arg2())));
            break;
        case cP2_(0, MULT_SQUARED):
            separableMultiDistSquared(srcMultiArrayRange(in3D), destMultiArray(out3D), backgroundMode == 1);
            break;
        case cP2_(1, MULT_SQUARED):
            separableMultiDistSquared(srcMultiArrayRange(in3D), destMultiArray(out3D), 0);
            separableMultiDistSquared(srcMultiArrayRange(in3D), destMultiArray(tmp3D), 1);
            combineTwoMultiArrays(
                        srcMultiArrayRange(tmp3D),
                        srcMultiArray(out3D),
                        destMultiArray(out3D),
                        ifThenElse(Arg1() > Param(0.0), sq(sqrt(Arg1())-Param(0.5)), -sq(sqrt(Arg2())-Param(0.5))));
            break;
        case cP2_(0, IMAG_DIST_TRANS):
            distanceTransform(srcImageRange(in), destImage(out),backgroundValue, norm);
            break;
        default:
            mexErrMsgTxt("Precondition checking not complete - something went wrong");
    }

}

/***************************************************************************************************
**           VIGRA GATEWAY                                                                        **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */
    //Add classes as you feel
    switch(inputs.typeOf(0))
    {
        ALLOW_FD
	ALLOW_UINT_8_64
	ALLOW_INT_8_64
        default:
	    mexErrMsgTxt("Type of input 0 not supported");
    }
}

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

