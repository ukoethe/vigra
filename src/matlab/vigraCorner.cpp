/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/cornerdetection.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
    declScalarMinMax(double, scale, 1.0, 0.0, "inf");
    declCharConstr4(method, Corner, Beaudet, Foerstner,Rohr);
    declOut(double);
    
    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :            base_data<T>(inputs),
                 initOption(scale), 
				 initOption(method)
    {
        if(this->numOfDim != 2)
            mexErrMsgTxt("vigraCorner only operates on 2D Images");
        initOut_SAME(double);
    }
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
        switch(o.method){
            case data<T>::Corner:
                cornerResponseFunction (srcImageRange(o.in), destImage(o.out), o.scale);
                break;
            case data<T>::Foerstner:
                foerstnerCornerDetector (srcImageRange(o.in), destImage(o.out), o.scale);
                break;
            case data<T>::Rohr:
                rohrCornerDetector (srcImageRange(o.in), destImage(o.out), o.scale);
                break;
            case data<T>::Beaudet:
                beaudetCornerDetector (srcImageRange(o.in), destImage(o.out), o.scale);
                break;
            default:
                mexErrMsgTxt("Some Error occured");
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
