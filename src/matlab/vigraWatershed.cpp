/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <string>
#include <vigra/watersheds.hxx>
#include <vigra/watersheds3D.hxx>
#include <vigra/matlab.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/seededregiongrowing3D.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraWatershed

using namespace vigra;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
	declScalar2D3D(int, conn, 8, 26);
	declCharConstr(method,2, UNION, SEED, a, e, i);
	declCharConstr(crack, 2, completeGrow, keepContours, o, u, n);
	declOut(Int32);
	
	SRGType SRGcrack;
	BasicImageView<T>  seed;
	MultiArrayView<3,T> seed3D;	
	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data(inputs),
				map(conn), map(crack)
	{
		if(numOfDim == IMAG && conn != 8 && conn !=4){
			conn = 8;
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 8");
		}else if(numOfDim == VOLUME && conn != 26 && conn !=6){
			conn = 26;
			mexWarnMsgTxt("Invalid User supplied connectivity - using Default: 26");
		}
		method = UNION;
		if(inputs.size()!= 1){
			mxArray* mxArrSeed =mxGetField(inputs[1], 0, "seed");
			if(mxArrSeed != NULL && mxIsNumeric(mxArrSeed)){
				if(numOfDim == IMAG){
					seed = matlab::getImage<T>(mxArrSeed);
					seed3D = matlab::getMultiArray<3, T>(mxArrSeed);
					conn = 4;
				}else{
					seed3D = matlab::getMultiArray<3, T>(mxArrSeed);
				}
				method = SEED;
				if(seed3D.shape() != in3D.shape()){
					conn = 6;
					method = UNION;
					mexWarnMsgTxt("Seed and Input dimension mismatch. Using UNION algorithm");
				}
			}
		}
		if(crack == completeGrow)SRGcrack = vigra::CompleteGrow;
		else SRGcrack = vigra::KeepContours;
		mapOut_SAME(Int32);
	}
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP3_(a, b , c) cP3<data<T>::a, data<T>::b, c>::value
struct vigraFunctor
{
	template <class T>
	static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
		//Options
		data<T>  o(outputs, inputs);

		// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
		int max_region_label = -1;
		// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
		switch(cantorPair(o.method, o.numOfDim, o.conn)){
			//cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
			//Code is basically the code on the VIGRA-reference page 
			case cP3_(UNION, IMAG, 8):	
				max_region_label = watersheds(srcImageRange(o.in), destImage(o.out));
				mexWarnMsgTxt("UNION IMAG 8");
				break;
			case cP3_(UNION, IMAG, 4):
				max_region_label = watersheds(srcImageRange(o.in), destImage(o.out), FourNeighborCode());
				mexWarnMsgTxt("UNION IMAG 4");
				break;
			case cP3_(UNION, VOLUME, 26):
				max_region_label = watersheds3DTwentySix(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D));
				mexWarnMsgTxt("UNION VOLUME 26");
				break;
			case cP3_(UNION, VOLUME, 6):
				max_region_label = watersheds3DSix(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D));
				mexWarnMsgTxt("UNION VOLUME 6");
				break;
			case cP3_(SEED, IMAG, 4):
			{				
				//find maximimum of seed Image
				mexWarnMsgTxt("SEED IMAG 4");
				FindMinMax<T> minmax;   // init functor
				inspectImage(srcImageRange(o.seed), minmax);
				max_region_label = minmax.max;
				
				ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
												gradstat(max_region_label);
												
				seededRegionGrowing(	srcImageRange(o.in), 
										srcImage(o.seed),
										destImage(o.out), 
										gradstat,o.SRGcrack );
				break;
			}
			case cP3_(SEED, VOLUME, 6):
			{
				mexWarnMsgTxt("SEED IMAG 4");
				FindMinMax<T> minmax;
				inspectMultiArray(srcMultiArrayRange(o.seed3D), minmax);
				max_region_label = minmax.max;
				
				ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
												gradstat(max_region_label);
												
				seededRegionGrowing3D(	srcMultiArrayRange(o.in3D),
										srcMultiArray(o.seed3D),
										destMultiArray(o.out3D),
										gradstat, o.SRGcrack);
				break;
			}		
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
function L = vigraWatershed(inputArray)
function L = vigraWatershed(inputArray, options);

L = vigraWatershed(inputArray) computes a label matrix identifying the watershed regions of the inputArray , which may be 2 or 3 dimensional. The elements of L are Int32 values  greater than 0. 
There are no watershedpixels. Do not use this  with Arrays of type other than single or double.
L = vigraWatershed(inputImage, options)  does the same with user options.
options is a struct with possible fields: "seed", "conn" and "crack

"seed": 				An Array of same size as inputArray. If supplied seeded Region growing shall be used. If not the union find algorithm will be used.
					As of now seed has to be of same type as the inputArray. - this will be changed in the next update.
"conn": 		 		2D 4(default for Seed.), 8 (default for Union) . 3D 6 (default for Seed), 26 (default for union);
					While using seeded region growing only the default values are used.
"crack"				'completeGrow' (Default), 'keepContours'
					Choose whether to keep watershed pixels or not.  Only used while using Seeded Region Growing.
Usage:
	opt = struct('fieldname' ,'value',....);
	out = vigraFunc(in, opt);
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
	// 
	callMexFunctor<vigraFunctor>(outputs, inputs);
}

