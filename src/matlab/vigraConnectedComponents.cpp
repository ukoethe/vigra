/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/labelimage.hxx>
#include <vigra/labelvolume.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraConnectedComponents

using namespace vigra;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
	declScalarMinMax(int, conn, 8, 4, 26);
    T backgroundValue;
	bool hasbackground;
	T get_backgroundValue(matlab::InputArray inputs)
	{
		if(inputs.size() == 1){
			hasbackground = false;
			return 0;
		}else{
			mxArray* backgroundValue =mxGetField(inputs[1], 0, "backgroundValue");
			if(backgroundValue!=NULL&&mxIsNumeric(backgroundValue)){
				hasbackground = true;
				return matlab::getScalar<T>(backgroundValue);
			}else{
				hasbackground = false;
				return 0;
			}
		}
	}
	declOut(double);

	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data(inputs),
				map(backgroundValue), map(conn)
	{
		if(numOfDim == 3)conn = 26;
		
		if(numOfDim == IMAG && conn != 4 && conn !=8){
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 8");
		}else if(numOfDim == VOLUME && conn != 6 && conn !=26){
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 26");
		}
		mapOut_SAME(double);
	}
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP3_(a, b , c) cP3<a, data<T>::b, c>::value
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
function D = vigraRadialSymmetry(inputArray)
function D = vigraradialSymmetry(inputArray, options);

D = vigraConnectedComponents(inputArray) computes the Fast Radial Symmetry Transform using the default options. see vigra::RadialSymmetryTransform
for more information.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.
options is a struct with possible fields: "method", "backgroundMode" and "backgroundPixel" and "norm"

"scale": 				1.0(default),any floating point value
						scale parameter for the vigraRadialSymmetry


Usage:
	opt = struct('method' ,value);
	out = vigraConnectedComponents(in, opt);

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}

