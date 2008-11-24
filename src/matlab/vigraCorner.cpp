/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/cornerdetection.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraCorner

using namespace vigra;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
	declScalarMinMax(double, scale, 1.0, 0.0, "inf");
	declCharConstr(method, 4, Corner, Beaudet,Foerstner,Rohr, a);
	declOut(double);
	
	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data(inputs),
				map(scale), map(method)
	{
	if(numOfDim != 2)
		mexErrMsgTxt("vigraCorner only operates on 2D Images");
	mapOut_SAME(double);
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
function D = vigraCorner(inputArray)
function D = vigraCorner(inputArray, options);

D = vigraCorner(inputArray) computes the Fast Radial Symmetry Transform using the default options. see vigra::RadialSymmetryTransform
for more information.
D = vigraCorner(inputImage, options)  does the same with user options.
options is a struct with possible fields: "method", "backgroundMode" and "backgroundPixel" and "norm"

"scale": 				1.0(default),any floating point value
						scale parameter for the vigraCornerdetector


Usage:
	opt = struct('method' ,value);
	out = vigraConnectedComponents(in, opt);

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}

