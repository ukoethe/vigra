/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/resizeimage.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraResize3

using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.

template <class T>
struct data: public base_data<T>{

	
	declCharConstr(method, 3, BSpline, CatmullRom, Coscot, a, e);
	declScalarMinMax(int, splineOrder, 3, 0, 5);
	declOut(T);
	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data<T>(inputs),
				map(method),
				map(splineOrder)
	{
		if(inputs.size() == 2)
		{
			mxArray* shape =mxGetField(inputs[1], 0, "shape");
			if(shape != NULL && mxIsNumeric(shape))
			{
				MultiArrayShape<2>::type newShape = matlab::getShape<2>(shape);		
				mapOut_3D(T, newShape[0], newShape[1], this->in3D.shape(2));
			}
			else out3D = matlab::createMultiArray<3,T>(this->in3D.shape(), outputs[0]);
		}
		else out3D = matlab::createMultiArray<3,T>(this->in3D.shape(), outputs[0]);
		
		if(method == Coscot || method == CatmullRom){
			mexWarnMsgTxt("Ignoring splineOrder parameter.");
			splineOrder = 0;
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

	    for(int k=0; k<o.in3D.shape(3); ++k)
	    {
			BasicImageView<T> ink = makeBasicImageView(o.in3D.bindOuter(k));    
			BasicImageView<T> outk = makeBasicImageView(o.out3D.bindOuter(k));    
	        switch(cantorPair(o.method, o.splineOrder))
	        {
	        case cP2_(BSpline, 0):
	            resizeImageNoInterpolation(srcImageRange(ink), destImageRange(outk));
	            break;
	        case cP2_(BSpline, 1):
				resizeImageLinearInterpolation(srcImageRange(ink), destImageRange(outk));
	            break;
	        case cP2_(BSpline, 2):
				resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<2>());
	            break;
	        case cP2_(BSpline, 3):
				resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<3>());
	            break;
	        case cP2_(BSpline, 4):
				resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<4>());
	            break;
	        case cP2_(BSpline, 5):
				resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<5>());
	            break;
			case cP2_(CatmullRom, 0):
				resizeImageCatmullRomInterpolation(srcImageRange(ink), destImageRange(outk));
	            break;
			case cP2_(Coscot, 0):
				resizeImageCoscotInterpolation(srcImageRange(ink), destImageRange(outk));
	            break;			
	        default:
	            mexErrMsgTxt("Something went wrong");
	        }
	    }
		
	}
};

/** MATLAB 
function resized = vigraResize2(original)
function resized = vigraResize2(original, options)

D = vigraresize3(inputArray) resizes original Volume data with default options.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.

	original    - Rrray with original image data 
					(gray scale or multi-band/RGB, type)
	options
		shape    - int32 array of length 2 giving the new shape
			(default: same size as old Image)
		splineOrder - order of interpolation 
			(0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
		method- BSpline(default), Cascot or CatmullRom
			only BSpline uses the splineOrder method
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}
