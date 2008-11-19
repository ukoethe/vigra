/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/labelimage.hxx>
#include <vigra/labelvolume.hxx>


//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T
#define vigraFunc vigraConnectedComponents
#define outType double
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
template <class T>
struct options{
	int numOfDim;
	int conn;
	int backgroundValue;
	

	BasicImageView<T>  in;
	MultiArrayView<3,T> in3D;

	BasicImageView<outType>  out;
	MultiArrayView<3,outType> out3D;
	
	options(int a, int b, int d){
		numOfDim = a;
		conn = b; 
		backgroundValue = d;
	}
};

// Quick and dirty macro for filling numerical options fields.
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
	enum {FourNeighbor = 4, SixNeighbor = 6, EightNeighbor = 8, TSixNeighbor = 26}nhood;
	
	//Default Options
	options<T> opt(mxGetNumberOfDimensions(inputs[0]), 
				mxGetNumberOfDimensions(inputs[0])==VOLUME? TSixNeighbor:EightNeighbor,
				-1);
	
	{//Preconditions on default options
	if(opt.numOfDim > VOLUME)  
						mexErrMsgTxt("Currently InputArray may only have 2 or 3 dimensions");
	if(inputs.isEmpty(0)) 
						mexErrMsgTxt("Input Image is empty!");
					
	}	
	
	//Map data to option fields
	if(opt.numOfDim == IMAG){
		opt.in = matlab::getImage<T>(inputs[0]);
		opt.out = matlab::createImage<outType>(opt.in.width(), opt.in.height(), outputs[0]);
		opt.in3D = matlab::getMultiArray<3, T>(inputs[0]);
		//Lets out3D View the same data as out.
		opt.out3D = MultiArrayView<3, outType>(opt.in3D.shape(), (outType*)opt.out.data());
	}else{
		opt.in3D = matlab::getMultiArray<3, T>(inputs[0]);
		opt.out3D = matlab::createMultiArray<3,outType>(opt.in3D.shape(), outputs[0]);
	}
	
	//User supplied Options
	if(inputs.isValid(1)){	
		fillOptNumericField(backgroundValue);
		fillOptNumericField(conn);					
		if(opt.numOfDim == IMAG && opt.conn != FourNeighbor && opt.conn !=EightNeighbor){
			opt.conn = 	EightNeighbor;
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 8");
		}else if(opt.numOfDim == VOLUME && opt.conn != SixNeighbor && opt.conn !=TSixNeighbor){
			opt.conn = TSixNeighbor;
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 26");
		}
	}

	

	// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
	int max_region_label = 0;
	switch(cantorPair(opt.backgroundValue != -1,opt.numOfDim, opt.conn)){
		//cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
		//Code is basically the code on the VIGRA-reference page 
		case cP3<0, IMAG, EightNeighbor>::value:	
			max_region_label = labelImage(srcImageRange(opt.in), destImage(opt.out), true);
			mexWarnMsgTxt("Breakpoint");
			break;
		case cP3<0, IMAG, FourNeighbor>::value:
			max_region_label = labelImage(srcImageRange(opt.in), destImage(opt.out), false);
			break;
		case cP3<0, VOLUME, TSixNeighbor>::value:
			max_region_label = labelVolumeSix(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D));
			break;
		case cP3<0, VOLUME, SixNeighbor>::value:
			max_region_label = labelVolume(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), NeighborCode3DTwentySix());
			break;
		case cP3<1, IMAG, EightNeighbor>::value:	
			max_region_label = labelImageWithBackground(srcImageRange(opt.in), destImage(opt.out), true, opt.backgroundValue);
			break;
		case cP3<1, IMAG, FourNeighbor>::value:
			max_region_label = labelImageWithBackground(srcImageRange(opt.in), destImage(opt.out), false, opt.backgroundValue);
			break;
		case cP3<1, VOLUME, TSixNeighbor>::value:
			max_region_label = labelVolumeWithBackground(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), 
																						NeighborCode3DSix(), opt.backgroundValue);
			break;
		case cP3<1, VOLUME, SixNeighbor>::value:
			max_region_label = labelVolumeWithBackground(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D), 
																						NeighborCode3DTwentySix(), opt.backgroundValue);
			break;
		default:
			mexErrMsgTxt("Something went wrong");
	}	
	
	// Are there more than one output? nargout.
	
}



/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* DELETE LINES IF A CERTAIN CLASS IS NOT SUPPORTED
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB 
function D = vigraConnectedComponents(inputArray)
function D = vigraConnectedComponents(inputArray, options);

D = vigraConnectedComponents(inputArray) computes the ConnectedComponents using the default options.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.
options is a struct with possible fields: "method", "backgroundMode" and "backgroundPixel" and "norm"

"conn": 				8(default),4 for 2D and 26(default) 6 for 3D
						Use selected neighborhood for checking connectedness.
"backgroundValue": 		0 (default) ,	
						Ignore Pixel with 'backgroundValue' in Connected Components labeling. 


Usage:
	opt = struct('method' ,value);
	out = vigraConnectedComponents(in, opt);

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
