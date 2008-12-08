/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/labelimage.hxx>
#include <vigra/labelvolume.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraConnectedComponents

using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
    declScalar2D3D(int, conn, 4, 6);

    T backgroundValue;
    bool hasbackground;
    T get_backgroundValue()
    {
        if(!this->options.isValid("backgroundValue"))
        {
            hasbackground = false;
            return 0;
        }
        const mxArray* name = this->options["backgroundValue"];
        if(!mxIsNumeric(name))
            mexErrMsgTxt("option 'backgroundValue' must be a numeric value.");
        hasbackground = true;
        return matlab::getScalar<T>(name);
    }

    declOut(double);

    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                map(backgroundValue), map(conn)
    {
        
        if(this->numOfDim == IMAG && conn != 8 && conn !=4)
            mexErrMsgTxt("vigraConnectedComponents(): Connectivity for 2D images must be 4 or 8");
        if(this->numOfDim == VOLUME && conn != 26 && conn !=6)
             mexErrMsgTxt("vigraConnectedComponents(): Connectivity for 3D volumes must be 6 or 26");
        mapOut_SAME(double);
    }
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP3_(a, b , c) cP3<a, b, c>::value
struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
    int max_region_label = 0;
    switch(cantorPair(o.hasbackground,o.numOfDim, o.conn)){
        //cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
        //Code is basically the code on the VIGRA-reference page 
        case cP3_(0, IMAG, 8):  
            max_region_label = labelImage(srcImageRange(o.in), destImage(o.out), true);
            mexWarnMsgTxt("Breakpoint");
            break;
        case cP3_(0, IMAG, 4):
            max_region_label = labelImage(srcImageRange(o.in), destImage(o.out), false);
            break;
        case cP3_(0, VOLUME, 26):
            max_region_label = labelVolumeSix(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D));
            break;
        case cP3_(0, VOLUME, 6):
            max_region_label = labelVolume(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), NeighborCode3DTwentySix());
            break;
        case cP3_(1, IMAG, 8):  
            max_region_label = labelImageWithBackground(srcImageRange(o.in), destImage(o.out), true, o.backgroundValue);
            break;
        case cP3_(1, IMAG, 4):
            max_region_label = labelImageWithBackground(srcImageRange(o.in), destImage(o.out), false, o.backgroundValue);
            break;
        case cP3_(1, VOLUME, 26):
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 
                                                                                        NeighborCode3DSix(), o.backgroundValue);
            break;
        case cP3_(1, VOLUME, 6):
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), 
                                                                                        NeighborCode3DTwentySix(), o.backgroundValue);
            break;
        default:
            mexErrMsgTxt("Something went wrong");
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
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
    // 
    callMexFunctor<vigraFunctor>(outputs, inputs);
}

