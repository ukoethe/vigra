#include <iostream>
#include <string>
#include <vigra/watersheds.hxx>
#include <vigra/watersheds3D.hxx>
#include <vigra/matlab.hxx>
#include <vigra/basicimageview.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/seededregiongrowing3D.hxx>

using namespace vigra;

int cantorPair(int x, int y){
		return (int)(((x+y)*(x+y+1))/2+y);
}

template <int x, int y>
struct cP{
	enum { value = (int)(((x+y)*(x+y+1))/2+y)};
};





template <class T>
void vigraWatershed(matlab::OutputArray outputs, matlab::InputArray inputs){
	//number of Dimensions
	int numOfDim = mxGetNumberOfDimensions(inputs[0]);
	if(numOfDim > 3)	mexErrMsgTxt("Currently InputArray may only have 2 or 3 dimensions");
	if(inputs.isEmpty(0)) mexErrMsgTxt("Input Image is empty!");
	
	//If seed Matrix exists use seeded region growing.
	std::string method = inputs.isValid(2)?   "seed":"union";
	
	//connectivity default connectivity for Watershed is 8/26, for seeded Region 4/6.
	int conn = (inputs.isValid(1)&&!inputs.isEmpty(1))
				?	matlab::getScalar<int>(inputs[1])
				:	method == "seed"
					? numOfDim == 2
						?	4:6
					:numOfDim == 2
						?	8:26;
	//connectivity check is done in the switch.
	
	
	int max_region_label = -1;
	// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
	int index = cantorPair(numOfDim, conn);
	
	
	if(method == "union"){
		switch(index){
			//cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
			case cP<2, 8>::value:
			{
				//initialisation of input and output Arrays
				BasicImageView<T> in = matlab::getImage<T>(inputs[0]);
				BasicImageView<Int32> labeling = matlab::createImage<Int32>(in.width(), in.height(), outputs[0]);	
				
				max_region_label = watersheds(srcImageRange(in), destImage(labeling));
				break;
			}	
			case cP<2, 4>::value:
			{
				BasicImageView<T> in = matlab::getImage<T>(inputs[0]);
				BasicImageView<Int32> labeling = matlab::createImage<Int32>(in.width(), in.height(), outputs[0]);
				
				max_region_label = watersheds(srcImageRange(in), destImage(labeling), FourNeighborCode());
				break;
			}
			case cP<3,26>::value:
			{
				MultiArrayView<3,T> in = matlab::getMultiArray<3, T>(inputs[0]);
				MultiArrayView<3,Int32> out = matlab::createMultiArray<3,Int32>(in.shape(), outputs[0]);	
				
				max_region_label = watersheds3DTwentySix(srcMultiArrayRange(in), destMultiArray(out));
				break;
			}	
			case cP<3, 6>::value:
			{
				MultiArrayView<3,T> in = matlab::getMultiArray<3,T>(inputs[0]);
				MultiArrayView<3,Int32> out = matlab::createMultiArray<3,Int32>(in.shape(), outputs[0]);	

				max_region_label = watersheds3DSix(srcMultiArrayRange(in), destMultiArray(out));
				break;
			}
			default:
				mexErrMsgTxt("Connectivity only be 4 or 8 for 2D and 6 or 26 for 3D Arrays");
		}
	}else{
		switch(index){
			case cP<2, 4>::value:
			{				
				BasicImageView<T> in = matlab::getImage<T>(inputs[0]);
				BasicImageView<Int32> labeling = matlab::createImage<Int32>(in.width(), in.height(), outputs[0]);	
				BasicImageView<T> seed = matlab::getImage<T>(inputs[2]);
				
				
				if(seed.width() != in.width() || seed.height() != in.height())
					mexErrMsgTxt("Seed Matrix and Heightmap mus have same size");
				
				//find maximimum of seed Image
				FindMinMax<T> minmax;   // init functor
				inspectImage(srcImageRange(seed), minmax);
				max_region_label = minmax.max;
				
				ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
												gradstat(max_region_label);
				
				seededRegionGrowing(	srcImageRange(in), 
										srcImage(seed),
										destImage(labeling), 
										gradstat,KeepContours );
				break;
			}
			case cP<3,6>::value:
			{
				MultiArrayView<3,T> in = matlab::getMultiArray<3,T>(inputs[0]);
				MultiArrayView<3,Int32> labeling = matlab::createMultiArray<3,Int32>(in.shape(), outputs[0]);	
				MultiArrayView<3,T> seed = matlab::getMultiArray<3,T>(inputs[2]);

				if(seed.shape() != in.shape())
					mexErrMsgTxt("Seed Matrix and Heightmap must have same size");
				
				FindMinMax<T> minmax;
				inspectMultiArray(srcMultiArrayRange(seed), minmax);
				max_region_label = minmax.max;
				
				ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
												gradstat(max_region_label);
												
				seededRegionGrowing3D(	srcMultiArrayRange(in),
										srcMultiArray(seed),
										destMultiArray(labeling),
										gradstat, KeepContours);
				break;
			}
			default:
				mexErrMsgTxt("Connectivity using seeded region growing may only be 4(2D) or 6(3D)");	
		}
		
	}
}

/** MATLAB 
function L = vigraWatershed(inputArray)
function L = vigraWatershed(inputArray, connectivity);
function L = vigraWatershed(inputArray, [], seedArray);

L = vigraWatershed(inputArray) computes a label matrix identifying the watershed regions of the inputArray , which may be 2 or 3 dimensional. The elements of L are Int32 values  greater than 0. 
There are no watershedpixels.
L = vigraWatershed(inputImage, connectivity)  does same as above using the connectivity given. connectivity may be 4 or 8 for 2dimensional and 6 or 26 for 3Dimensional data.
default values for connectivity are 8/26
L = vigraWatershed(inputArray, [], seedArray) uses seeded Region growing to identify the watershed regions. seedArray must be an array of same size and type as inputArray. The seeds must be consecutively labeled with whole numbers. 
do not give connectivity using this Algorithm as currently only 4 and 6 neighborhoods are supported.  Pixels labeled as 0 are watershedpixels and devide different watershed regions.

L is of type Int32

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    if(matlab::ValueType<UInt8>::check(inputs[0]))
        vigraWatershed<UInt8>(outputs, inputs);
    else if(matlab::ValueType<float>::check(inputs[0]))
        vigraWatershed<float>(outputs, inputs);
    else if(matlab::ValueType<double>::check(inputs[0]))
        vigraWatershed<double>(outputs, inputs);
    else if(matlab::ValueType<UInt16>::check(inputs[0]))
        vigraWatershed<UInt8>(outputs, inputs);
    else if(matlab::ValueType<UInt32>::check(inputs[0]))
        vigraWatershed<float>(outputs, inputs);
    else if(matlab::ValueType<UInt64>::check(inputs[0]))
        vigraWatershed<double>(outputs, inputs);
    else if(matlab::ValueType<Int8>::check(inputs[0]))
        vigraWatershed<double>(outputs, inputs);
    else if(matlab::ValueType<Int16>::check(inputs[0]))
        vigraWatershed<UInt8>(outputs, inputs);
    else if(matlab::ValueType<Int32>::check(inputs[0]))
        vigraWatershed<float>(outputs, inputs);
    else if(matlab::ValueType<Int64>::check(inputs[0]))
        vigraWatershed<double>(outputs, inputs);
	else
        mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");
}
