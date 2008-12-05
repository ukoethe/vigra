/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>


//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraDistance

using namespace vigra;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
	declScalarMinMax(bool, backgroundMode, 0, 0, 1);
	declCharConstr(method, 4, MULT, MULT_SQUARED, IMAG_DIST_TRANS, INVERTEDCRAP, e);
	declOut(double);
	
	//Only supported with MULT
	TinyVector<double, 3> pitch3D;
	//Only supported with IMAG_DIST_TRANS
	declScalar(T, backgroundPixel, 0);
	declScalarMinMax(int, norm, 2, 0, 2);
	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data(inputs), map(backgroundMode), map(method), map(backgroundPixel), map(norm),
				pitch3D(1.0,1.0,1.0)
	{
		mapOut_SAME(double);
		if(inputs.size() == 2)
		{
			mxArray* pitchArr =mxGetField(inputs[1], 0, "pitch");
			if(pitchArr != NULL && mxIsNumeric(pitchArr))
			{
				if(numOfDim == IMAG){
					TinyVectorView<double, 2> temp = matlab::getVector<double,2>(pitchArr);
					pitch3D[0] = temp[0]>0? temp[0]:1.0;
					pitch3D[1] = temp[1]>0? temp[1]:1.0;
				}else{
					TinyVectorView<double, 3> temp = matlab::getVector<double,3>(pitchArr);
					pitch3D = temp;
				}
				if(this->method != MULT){
					mexWarnMsgTxt("pitch option is only supported with the method MULT");
				}
			}
		}
		else
		if(this->numOfDim == VOLUME && this->method == IMAG_DIST_TRANS){
			this->method = MULT;
			mexWarnMsgTxt("IMAG_DIST_TRANS only works with 2D Images using default: MULT");
		}
	}
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP2_(a, b) cP<data<T>::a, data<T>::b>::value
struct vigraFunctor
{
	template <class T>
	static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
		//Options
		data<T>  o(outputs, inputs);

		// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
		switch(cantorPair(o.numOfDim, o.method)){
			//In this case function pointers may have been more elegant.
			case cP2_(IMAG, MULT):
				separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode, o.pitch3D);
				break;
			case cP2_(VOLUME, MULT):
				separableMultiDistance(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode, o.pitch3D);
				break;
			case cP2_(IMAG, MULT_SQUARED):
				separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode);
				break;
			case cP2_(VOLUME, MULT_SQUARED):
				separableMultiDistSquared(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D), o.backgroundMode);
				break;
			case cP2_(IMAG, IMAG_DIST_TRANS):
				distanceTransform(srcImageRange(o.in), destImage(o.out),o.backgroundPixel, o.norm);
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

D = vigraDistance(inputArray) computes the Distancetransform using the default options.
D = vigraDistance(inputImage, options)  does the same with user options.
options is a struct with possible fields: "method", "backgroundMode" and "backgroundPixel" and "norm"

"method": 				'MULT'(default), 'MULT_SQUARED', 'IMAG_DIST_TRANS'
						MULT and MULT_SQUARED are the faster and newer distance transform methods defined with VIGRA-Multiarrays
						(vigra::seperableMultiDist....).MULT_SQUARED returns the squared values of MULT. 
						IMAG_DIST_TRANS is defined with BasicImage and is slower. Use it only if  you explicitely need the backgroundPixel
						or norm option. (vigra::distanceTransform)
"backgroundMode": 		0 (default) ,  1: 	
						This option is only used with methods MULT and MULT_SQUARED
						If the parameter background is 1, then the squared distance of all background pixels to the nearest object is calculated. 
						Otherwise, the distance of all object pixels to the nearest background pixel is calculated.
"backgroundPixel: 		0(default) , arb. value in grayscale range: 		
						This option defines the background Pixel value. Only used with method = IMAG_DIST_TRANS
"norm": 					2(default), 1 , 0.  
						Defines the norm used to calculate the distance (Only used with method = IMAG_DIST_TRANS)
"pitch"				[1.0, 1.0]  (2D)or [1.0, 1.0, 1.0](3D) , arb 2D or 3D array.
						define a pitch if data has non-uniform resolution .

Usage:
	opt = struct('method' ,'IMAGE_DIST_TRANS' , 'backgroundPixel', 10 , 'norm' , 0);
	out = vigraDistance(in, opt);

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}
