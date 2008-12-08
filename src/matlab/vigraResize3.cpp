/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/multi_resize.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraResize3

using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.
template <class T>
struct data{
    MultiArrayView<4, T> in; 
    MultiArrayView<4, T> out;
	
	declCharConstr(method, 2, BSpline, CatmullRom, o, a, e);
	declScalarMinMax(int, splineOrder, 3, 0, 5);

	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			map(method),
				map(splineOrder)
	{
		in = matlab::getMultiArray<4, T>(inputs[0]);
		if(inputs.size() == 2)
		{
			mxArray* shape =mxGetField(inputs[1], 0, "shape");
			if(shape != NULL && mxIsNumeric(shape))
			{
				MultiArrayShape<3>::type newShape = matlab::getShape<3>(shape);		
				MultiArrayShape<4>::type newShape4(newShape[0], newShape[1], newShape[2], in.shape(3));
				out = matlab::createMultiArray<4,T>(newShape4,outputs[0]);
			}
			else out = matlab::createMultiArray<4,T>(in.shape(), outputs[0]);
		}
		else out = matlab::createMultiArray<4,T>(in.shape(), outputs[0]);
		
		if(method == CatmullRom){
			splineOrder = 0;
			mexWarnMsgTxt("Ignoring splineOrder parameter.");
		}
	}
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP2_(a, b) cP<data<T>::a, b>::value
struct vigraFunctor
{
	template <class T>
	static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
		//Options
		data<T>  o(outputs, inputs);

		// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 

	    for(int k=0; k<o.in.shape(3); ++k)
	    {
	        MultiArrayView<3, T> ink = o.in.bindOuter(k);    
	        MultiArrayView<3, T> outk = o.out.bindOuter(k);
	        switch(cantorPair(o.method, o.splineOrder))
	        {
	        case cP2_(BSpline, 0):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<0>());
	            break;
	        case cP2_(BSpline, 1):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<1>());
	            break;
	        case cP2_(BSpline, 2):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<2>());
	            break;
	        case cP2_(BSpline, 3):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<3>());
	            break;
	        case cP2_(BSpline, 4):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<4>());
	            break;
	        case cP2_(BSpline, 5):
	            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<5>());
	            break;
			case cP2_(CatmullRom, 0):
				resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), CatmullRomSpline<double>());
	            break;
	         default:
	            mexErrMsgTxt("Something went wrong");
	        }
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
function resized = vigraResize3(original)
function resized = vigraResize3(original, options)

D = vigraresize3(inputArray) resizes original Volume data with default options.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.

	original    - 3D matrix with original volume data 
					(gray scale or multi-band/RGB, type 'uint8', 'single' or 'double')
	options
		shape    - int32 array of length 3 giving the new shape
			(default: same size as old volume)
		splineOrder - order of interpolation 
			(0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
		useCatmull - use CatmullRom Splines - (Ignores splineOrder value)
			derfault: 0;
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}
