/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/symmetry.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T


using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
    declScalarMinMax(double, scale, 1.0, 0.0, "inf");
    declOut(double);
    
    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                initOption(scale)
    {
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
        radialSymmetryTransform(srcImageRange(o.in), destImage(o.out), o.scale);
        
    }
};


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the 
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB 
function D = vigraRadialSymmetry(inputImage)
function D = vigraradialSymmetry(inputImage, options);

D = vigraRadialSymmetry(inputImage) computes the Fast Radial Symmetry Transform 
            using default options, see vigra::RadialSymmetryTransform for more information.
D = vigraRadialSymmetry(inputImage, options)  does the same with user options.

inputImage - 2D input array
options    - a struct with following possible fields:
    'scale':    1.0 (default), any positive floating point value
                scale parameter for the vigraRadialSymmetry


Usage:
    opt = struct('method' ,value);
    out = vigraRadialSymmetry(in, opt);

*/
