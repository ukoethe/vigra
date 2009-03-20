/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/cornerdetection.hxx>


using namespace vigra;
using namespace matlab;





template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    BasicImageView<T>   in      =   inputs.getImage<T>(0, v_required());
    double              scale   =   inputs.getScalarMinMax<double>(1, v_default(1.0), 0.0, "inf");

    VIGRA_CREATE_ENUM_AND_STD_MAP4(Methods,MapName, Corner, Foerstner, Rohr, Beaudet);
    Methods             method  =   (Methods)inputs.getEnum(2,  v_default((int)Corner), MapName);

    BasicImageView<double> out  =   outputs.createImage<double>(0, v_required(), in.width(), in.height());

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    switch(method){
        case Corner:
            cornerResponseFunction (srcImageRange(in), destImage(out), scale);
            break;
        case Foerstner:
            foerstnerCornerDetector (srcImageRange(in), destImage(out), scale);
            break;
        case Rohr:
            rohrCornerDetector (srcImageRange(in), destImage(out), scale);
            break;
        case Beaudet:
            beaudetCornerDetector (srcImageRange(in), destImage(out), scale);
            break;
        default:
            mexErrMsgTxt("Some Error occured");
    }

}



/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
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


/** MATLAB
function D = vigraCorner(inputImage)
function D = vigraCorner(inputImage, options);

D = vigraCorner(inputArray) does Corner detection.
D = vigraCorner(inputImage, options)  does the same with user options.

inputImage - 2D input array
options    - struct with following possible fields:
   'method':  'Corner' (default, corenr response function according to Harris), 'Beaudet', 'Foerstner', 'Rohr'
              Use corresponding method to detect corners (see vigra reference for more details).
   'scale':   1.0 (default), any positive floating point value
              scale parameter for corner feature computation

Usage:
    opt = struct('method' ,value);
    out = vigraCorner(in, opt);

*/
