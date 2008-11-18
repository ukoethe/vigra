/*
This is a template file used to create all porting functions of VIGRA
*/

/*++++++++++++++++++++++++++INCLUDES+++++++++++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>
#include <string>

#define vigraFunc vigraDistance
/*++++++++++++++++++++++++++HELPERFUNC+++++++++++++++++++++++++++++++*/
/* This is used for better readibility of the test cases            .
/* Nothing to be done here.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
using namespace vigra;

int cantorPair(int x, int y){
		return (int)(((x+y)*(x+y+1))/2+y);
}
int cantorPair(int x, int y, int z){
		return cantorPair(cantorPair(x,y),z);
}

template <int x, int y>
struct cP{
	enum { value = (int)(((x+y)*(x+y+1))/2+y)};
};

template <int x, int y, int z>
struct cP3{
	enum { value = cP<cP<x, y>::value, z>::value};
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* The Optons struct contains all the necassary working data and 
/* options for the vigraFunc. This is the minimal struct
/* Add fields as necessary
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T >
struct options{
	int numOfDim;
	bool backgroundMode;
	int method;
	
	//Currently only supported in the 2 Version
	int backgroundPixel;
	int norm;
	
	
	BasicImageView<T>  in;
	MultiArrayView<3,T> in3D;

	BasicImageView<double>  out;
	MultiArrayView<3,double> out3D;
	
	options(int nofDim, bool back, int meth, int backp, int nrm){
		numOfDim = nofDim;
		backgroundMode = back;
		method = meth;
		backgroundPixel = backp;
		norm = nrm;
	}
};

#define fillOptNumericField(name); \
		mxArray* name =mxGetField(inputs[1], 0, #name);\
		opt.name = (name!=NULL&&mxIsNumeric(name))?\
							mxGetScalar(name) : opt.name;


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
void vigraFunc(matlab::OutputArray outputs, matlab::InputArray inputs){
	// Constant definition for readibility
	enum {IMAG = 2, VOLUME = 3}dim;
	enum {MULT = 1, MULT_SQUARED = 2, IMAG_DIST_TRANS =3}meth;
	
	//Default Options
	options<T> opt(mxGetNumberOfDimensions(inputs[0]), false, MULT, 0 , 2);
	{//Preconditions on default options
	if(opt.numOfDim > VOLUME)  
						mexErrMsgTxt("Currently InputArray may only have 2 or 3 dimensions");
	if(inputs.isEmpty(0)) 
						mexErrMsgTxt("Input Image is empty!");
					
	}	
	
	//Map data to option fields
	if(opt.numOfDim == IMAG){
		opt.in = matlab::getImage<T>(inputs[0]);
		opt.out = matlab::createImage<double>(opt.in.width(), opt.in.height(), outputs[0]);
		opt.in3D = matlab::getMultiArray<3, T>(inputs[0]);
		opt.out3D = matlab::createMultiArray<3,double>(opt.in3D.shape(), outputs[0]);
	}else{
		opt.in3D = matlab::getMultiArray<3, T>(inputs[0]);
		opt.out3D = matlab::createMultiArray<3,double>(opt.in3D.shape(), outputs[0]);
	}	
	
	//User supplied Options
	if(inputs.isValid(1)){	
		fillOptNumericField(backgroundMode);
		
		mxArray* method =mxGetField(inputs[1], 0, "method");
		if(method!=NULL&&mxIsChar(method)){
			std::string meth = matlab::getString(method);
			
			if(meth == "MULT_SQUARED") 
				opt.method = MULT_SQUARED;
			else if(meth == "IMAG_DIST_TRANS" && opt.numOfDim == IMAG){
				opt.method = IMAG_DIST_TRANS;
				fillOptNumericField(backgroundPixel);
				fillOptNumericField(norm);
			}else if(meth == "IMAG_DIST_TRANS")
				mexWarnMsgTxt("IMAG_DIST_TRANS only valid for 2D images using default: MULT");
			else if(meth != "MULT")
				mexWarnMsgTxt("User supplied backgroundmode not supported using default: MULT");
		}
	}





	// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
	switch(cantorPair(opt.numOfDim, opt.method)){
		//In this case function pointers may have been more elegant.
		case cP<IMAG, MULT>::value:
			separableMultiDistance(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), opt.backgroundMode);
			break;
		case cP<VOLUME, MULT>::value:
			separableMultiDistance(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), opt.backgroundMode);
			break;
		case cP<IMAG, MULT_SQUARED>::value:
			separableMultiDistSquared(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), opt.backgroundMode);
			break;
		case cP<VOLUME, MULT_SQUARED>::value:
			separableMultiDistSquared(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), opt.backgroundMode);
			break;
		case cP<IMAG, IMAG_DIST_TRANS>::value:
			distanceTransform(srcImageRange(opt.in), destImage(opt.out),opt.backgroundPixel, opt.norm);
			break;
		default:
			mexErrMsgTxt("Precondition checking not complete - something went wrong");
	}
	
	// Are there more than one output? nargout.
	
}



/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* DELETE LINES IF A CERTAIN CLASS IS NOT SUPPORTED
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

Usage:
	opt = struct('method' ,'IMAGE_DIST_TRANS' , 'backgroundPixel', 10 , 'norm' , 0);
	out = vigraDistance(in, opt);

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
	mxClassID inClass = mxGetClassID(inputs[0]);
	switch(inClass){
		case mxDOUBLE_CLASS:
			vigraFunc<double>(outputs, inputs);	break;
		case mxSINGLE_CLASS:
			vigraFunc<float>(outputs, inputs);		break;
        case mxINT8_CLASS:
			vigraFunc<Int8>(outputs, inputs);		break;
		case mxINT16_CLASS:
			vigraFunc<Int16>(outputs, inputs);		break;
		case mxINT32_CLASS:
			vigraFunc<Int32>(outputs, inputs);		break;
		case mxINT64_CLASS:
			vigraFunc<Int64>(outputs, inputs);		break;
        case mxUINT8_CLASS:
			vigraFunc<UInt8>(outputs, inputs);		break;
		case mxUINT16_CLASS:
			vigraFunc<UInt16>(outputs, inputs);	break;
		case mxUINT32_CLASS:
			vigraFunc<UInt32>(outputs, inputs);	break;
		case mxUINT64_CLASS:
			vigraFunc<UInt64>(outputs, inputs);	break;		
		default:
			mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");
	}
}