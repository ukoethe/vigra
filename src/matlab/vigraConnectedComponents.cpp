/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/labelvolume.hxx>
#include <vigra/matlab_FLEXTYPE.hxx>


//this could be a typedef but if you want outType to be the same type as inType then you can just
//set outType to T


using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//#define RN_DEBUG
#define cP3_(a, b , c) cP3<a, b, c>::value
template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    //Options

    MultiArrayView<3,T>         in3D        = inputs.getMultiArray<3,T>(0, Required());
    BasicImageView<T>           in        = makeBasicImageView(in3D.bindOuter(0));
    int                         numOfDim = inputs.getDimOfInput(0, Required());
    int                         connectivity =  inputs.getScalar<int>("conn", Optional(8,26,numOfDim));
    {
        if(numOfDim == 2 && connectivity!= 8 && connectivity != 4)
            mexErrMsgTxt("Connectivity for 2D data must be 8 or 4");
        else if(numOfDim == 3 && connectivity!= 26 && connectivity != 6)
            mexErrMsgTxt("Connectivity for 3D data must be 26 or 6");
    }
    VarChecker  hasBackground;
    T                           backgroundValue = inputs.getScalar<T>("backgroundValue", Optional(hasBackground));
    MultiArrayView<3,T> out3D           = outputs.createMultiArray      <3,T>   (0, Required(), in3D.shape());
    BasicImageView<T> out(out3D.data(), in3D.shape(0), in3D.shape(1));

    int max_region_label = 0;


    switch(cantorPair(hasBackground.isSet , numOfDim, connectivity)){
        //cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
        //Code is basically the code on the VIGRA-reference page
        case cP3_(0, IMAG, 8):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, IMAG, 8");
            #endif
            max_region_label = labelImage(srcImageRange(in), destImage(out), true);
            break;
        case cP3_(0, IMAG, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, IMAG, 4");
            #endif
            max_region_label = labelImage(srcImageRange(in), destImage(out), false);
            break;
        case cP3_(0, VOLUME, 26):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, VOLUME, 26");
            #endif
            max_region_label = labelVolumeSix(srcMultiArrayRange(in3D), destMultiArray(out3D));
            break;
        case cP3_(0, VOLUME, 6):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, VOLUME, 6");
            #endif
            max_region_label = labelVolume(srcMultiArrayRange(in3D), destMultiArray(out3D), NeighborCode3DTwentySix());
            break;
        case cP3_(1, IMAG, 8):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, IMAG, 8");
            #endif
            max_region_label = labelImageWithBackground(srcImageRange(in), destImage(out), true, backgroundValue);
            break;
        case cP3_(1, IMAG, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, IMAG, 4");
            #endif
            max_region_label = labelImageWithBackground(srcImageRange(in), destImage(out), false, backgroundValue);
            break;
        case cP3_(1, VOLUME, 26):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, VOLUME, 26");
            #endif
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(in3D), destMultiArray(out3D),
                                                                                        NeighborCode3DSix(), backgroundValue);
            break;
        case cP3_(1, VOLUME, 6):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, VOLUME, 6");
            #endif
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(in3D), destMultiArray(out3D),
                                                                                        NeighborCode3DTwentySix(), backgroundValue);
            break;
        default:
            mexErrMsgTxt("Something went wrong");
    }

    outputs.createScalar<int> (1, Optional(), max_region_label);
}
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */

    mxClassID inClass;
    FLEX_TYPE(inClass, 0, in);
    switch(inClass)
    {
        ALLOW_D
        DEFAULT_ERROR;
    }
}


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB
function D = vigraConnectedComponents(inputArray)
function D = vigraConnectedComponents(inputArray, options);

D = vigraConnectedComponents(inputArray) performs connected components labeling
        using the default options.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.

inputArray - a 2D or 3D array of numeric type

options    - is a struct with the following possible fields:
    'conn':            The neighborhood to be used
                       2D: 4 (default) or  8
                       3D: 6 (default) or 26
    'backgroundValue': Specify the value of a background region not to be labeled (will be labeled 0
                       in the result array, even when it doesn't form a single connected component).
                       If this option in not present, the entire image/volume will be labeled.

Usage:
    opt = struct('conn', 8);
    out = vigraConnectedComponents(in, opt);

*/
