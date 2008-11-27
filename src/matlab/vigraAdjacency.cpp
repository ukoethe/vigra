/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/contourcirculator.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <vigra/diff2d.hxx>
#include <vigra/cellconfigurations.hxx>
#include <set>
#include <vigra/inspectimage.hxx>
#include <vigra/multi_pointoperators.hxx>


//this could be a typedef but if you want inType to be the same type as inType then you can just 
//set inType to T

#define vigraFunctor vigraAdjacency

using namespace vigra;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data: public base_data<T>{
	bool hasWatershedPixel;
	bool get_hasWatershedPixel(matlab::InputArray inputs)
	{
		if(inputs.size() == 1){
			mexWarnMsgTxt("It is recommended to set opt.hasWatershedPixel for speed");
			for(int ii = 0; ii < in3D.elementCount(); ii++){
				if(in3D[ii] == 0)return  1;
				
			}
			return 0;
		}else{
			bool in;
			mxArray* hasWatershedPixel =mxGetField(inputs[1], 0, "hasWatershedPixel");
			if(hasWatershedPixel!=NULL&&mxIsNumeric(hasWatershedPixel)){
				in = matlab::getScalar<bool>(hasWatershedPixel);
			}else{
				mexWarnMsgTxt("It is recommended to set opt.hasWatershedPixel for speed");
				for(int ii = 0; ii < in3D.elementCount(); ii++){
					if(in3D[ii] == 0)return  1;
				}
				return 0;
			}
			if(!is_in_range<bool>(in,0,1)){
				mexWarnMsgTxt("Invalid Value for has WatershedPixel != 0..1");
				for(int ii = 0; ii < in3D.elementCount(); ii++){
					if(in3D[ii] == 0)return  1;
					
				}
				return 0;
			}
			return in;
		}
	}
	SparseArray<int> adj_matrix;
	T max_region_label;
	
	data(matlab::OutputArray outputs, matlab::InputArray inputs)
	:			base_data(inputs),
				map(hasWatershedPixel)
	{
		FindMinMax<T> minmax;
		inspectMultiArray(srcMultiArrayRange(in3D), minmax);
		max_region_label = minmax.max;		
		adj_matrix.assign(max_region_label, max_region_label);
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
		//oions
		data<T>  o(outputs, inputs);
		
		// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
		switch(cantorPair(o.numOfDim, o.hasWatershedPixel)){
			//in case <XX, 0 (crackedges)> the Pixel/Voxel underneath, on the right (and behind) are considered for neighborhood.
			//in case<IMAG, 1>  the cell configuration around current pixel is compared to the ones in cellconfigurations.h (Configuration are Freeman coded (see PhD thesis, Chap 5.7)
			// if current configuration is of Type Line then adjacent pixels are considered as neighbors.
			// in case<VOLUME, 1> the number of adjacent regions are counted. the regions are considered adjacent if only two such regions exist
			case cP2_(IMAG, 0):
			{

				ImageIterator<T> upleft = o.in.upperLeft();
				ImageIterator<T> downright = o.in.lowerRight();
				Diff2D size = downright - upleft;
				for(int y = 0; y < size.y-1; ++y){
					for(int x = 0; x< size.x-1; ++x){
						o.adj_matrix(upleft(x,y)-1, upleft(x+1, y)-1)++;
						o.adj_matrix(upleft(x,y)-1, upleft(x, y+1)-1)++;
						o.adj_matrix(upleft(x+1, y)-1, upleft(x,y)-1)++;
						o.adj_matrix(upleft(x, y+1)-1, upleft(x,y)-1)++;
					}
				}
				
				for(int x = 0; x < size.x-1; ++x){
					o.adj_matrix(upleft(x, size.y-1)-1, upleft(x+1, size.y-1)-1)++;
					o.adj_matrix(upleft(x+1, size.y-1)-1, upleft(x, size.y-1)-1)++;
				}
				for(int y = 0; y < size.y-1; ++y){
					o.adj_matrix(upleft(size.x-1, y)-1, upleft(size.x-1, y+1)-1)++;
					o.adj_matrix(upleft(size.x-1, y+1)-1, upleft(size.x-1, y)-1)++;
				}
				break;
			}
			case cP2_(VOLUME, 0):
			{

				MultiArrayShape<3>::type sze = o.in3D.shape();
				for(int ii = 0; ii < sze[0]-1; ii++){
					for(int jj = 0; jj < sze[1]-1; jj++){
						for(int kk = 0; kk < sze[2]-1; kk++){
							o.adj_matrix(o.in3D(ii,jj,kk)-1, o.in3D(ii+1,jj,kk)-1)++;
							o.adj_matrix(o.in3D(ii,jj,kk)-1, o.in3D(ii,jj+1,kk)-1)++;
							o.adj_matrix(o.in3D(ii,jj,kk)-1, o.in3D(ii,jj,kk+1)-1)++;
							
							o.adj_matrix(o.in3D(ii+1,jj,kk)-1, o.in3D(ii,jj,kk)-1)++;
							o.adj_matrix(o.in3D(ii,jj+1,kk)-1, o.in3D(ii,jj,kk)-1)++;
							o.adj_matrix(o.in3D(ii,jj,kk+1)-1, o.in3D(ii,jj,kk)-1)++;
						}
					}
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					for(int kk = 0; kk < sze[2]-1; kk++){
						o.adj_matrix(o.in3D(sze[0]-1,jj,kk)-1, o.in3D(sze[0]-1,jj+1,kk)-1)++;
						o.adj_matrix(o.in3D(sze[0]-1,jj,kk)-1, o.in3D(sze[0]-1,jj,kk+1)-1)++;

						o.adj_matrix(o.in3D(sze[0]-1,jj+1,kk)-1, o.in3D(sze[0]-1,jj,kk)-1)++;
						o.adj_matrix(o.in3D(sze[0]-1,jj,kk+1)-1, o.in3D(sze[0]-1,jj,kk)-1)++;						
					}
				}
				for(int ii = 0; ii < sze[0]-1; ii++){
					for(int kk = 0; kk < sze[2]-1; kk++){
						o.adj_matrix(o.in3D(ii,sze[1]-1,kk)-1, o.in3D(ii+1,sze[1]-1,kk)-1)++;
						o.adj_matrix(o.in3D(ii,sze[1]-1,kk)-1, o.in3D(ii,sze[1]-1,kk+1)-1)++;	

						o.adj_matrix(o.in3D(ii+1,sze[1]-1,kk)-1, o.in3D(ii,sze[1]-1,kk)-1)++;
						o.adj_matrix(o.in3D(ii,sze[1]-1,kk+1)-1, o.in3D(ii,sze[1]-1,kk)-1)++;							
					}
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					for(int ii = 0; ii <sze[0] -1; ii++){
						o.adj_matrix(o.in3D(ii,jj,sze[2]-1)-1, o.in3D(ii,jj+1,sze[2]-1)-1)++;
						o.adj_matrix(o.in3D(ii,jj,sze[2]-1)-1, o.in3D(ii+1,jj,sze[2]-1)-1)++;
						
						o.adj_matrix(o.in3D(ii,jj+1,sze[2]-1)-1, o.in3D(ii,jj,sze[2]-1)-1)++;
						o.adj_matrix(o.in3D(ii+1,jj,sze[2]-1)-1, o.in3D(ii,jj,sze[2]-1)-1)++;
					}
				}
				
				for(int ii = 0; ii < sze[0]-1; ii++){
					o.adj_matrix(o.in3D(ii,sze[1]-1,sze[2]-1)-1, o.in3D(ii+1,sze[1]-1,sze[2]-1)-1)++;
					o.adj_matrix(o.in3D(ii+1,sze[1]-1,sze[2]-1)-1, o.in3D(ii,sze[1]-1,sze[2]-1)-1)++;
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					o.adj_matrix(o.in3D(sze[0]-1,jj,sze[2]-1)-1, o.in3D(sze[0]-1,jj+1,sze[2]-1)-1)++;
					o.adj_matrix(o.in3D(sze[0]-1,jj+1,sze[2]-1)-1, o.in3D(sze[0]-1,jj,sze[2]-1)-1)++;
				}
				for(int kk = 0; kk < sze[2]-1; kk++){
					o.adj_matrix(o.in3D(sze[0]-1,sze[1]-1,kk)-1, o.in3D(sze[0]-1,sze[1]-1,kk+1)-1)++;
					o.adj_matrix(o.in3D(sze[0]-1,sze[1]-1,kk+1)-1, o.in3D(sze[0]-1,sze[1]-1,kk)-1)++;
				}
				break;
			}
			case cP2_(IMAG, 1):
			{
				ImageIterator<T> upleft = o.in.upperLeft();
				ImageIterator<T> downright = o.in.lowerRight();
				Diff2D size = downright - upleft;
				for(int y = 1; y < size.y-1; ++y){
					for(int x = 1; x< size.x-1; ++x){
						if(upleft(x,y) == 0){
							NeighborhoodCirculator<BasicImageView<T>::Iterator, EightNeighborCode>
																circulator(upleft+ Diff2D(x, y));
							//circulator += 2;
							NeighborhoodCirculator<BasicImageView<T>::Iterator, EightNeighborCode>
																end(circulator);
							unsigned char BitField = 0;
							std::set<int> regions;
							
							do{
								if(*circulator == 0){
									BitField = BitField >> 1;
									BitField = BitField | 0x80;
								}else{
									regions.insert(*circulator);
									BitField = BitField >> 1;
									BitField = BitField | 0x00;
								}
							}while(++circulator != end);
							
							if(cellimage::cellConfigurations[(int)BitField] == cellimage::CellTypeLine){
								if(regions.size() == 2){
									std::set<int>::const_iterator iter = regions.begin();
									o.adj_matrix(*iter-1,*(iter++)-1)++;
								}
							}
						}
						
						std::ostringstream  oight;
					}
				}				
			}
				break;
			case cP2_(VOLUME, 1):
			{
				char ne[26*3] = 	{	1,	-1,	0,
				1, 	-1,	1,
				0,	-1,	1,
				-1,	-1,	1,
				-1,	-1,	0,
				-1,	-1,	-1,
				0,	-1,	-1,
				1,	-1,	-1,
				1,	0,	0,
				1, 	0,	1,
				0,	0,	1,
				-1,	0,	1,
				-1,	0,	0,
				-1,	0,	-1,
				0,	0,	-1,
				1,	0,	-1,
				1,	1,	0,
				1, 	1,	1,
				0,	1,	1,
				-1,	1,	1,
				-1,	1,	0,
				-1,	1,	-1,
				0,	1,	-1,
				1,	1,	-1,
				0,	1,	0,
				0,	-1,	0}	;
				std::set<int> regions;
				MultiArrayShape<3>::type sze = o.in3D.shape();
				for(int ii = 1; ii < sze[0]-1; ii++)
				{
					for(int jj = 1; jj < sze[1]-1; jj++)
					{
						for(int kk = 1; kk < sze[2]-1; kk++)
						{
							if(o.in3D(ii,jj,kk) == 0)
							{
								

								regions.clear();
								
								for(int ll = 0; ll < 26; ll++){
									if(o.in3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2]) != 0 ){
										regions.insert((o.in3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2])));
									}
								}
							
							if(regions.size() == 2)
								{
									std::set<int>::const_iterator iter = regions.begin();
									o.adj_matrix(*iter-1,*(iter++)-1)++;
								}
							}
						}
					}
				}
				
				break;
			}
		}
		for(int ii = 0; ii < o.max_region_label; ii++){
			for(int jj = ii; jj < o.max_region_label; jj++){
				if(ii == jj)o.adj_matrix(ii,ii) = 1;
				else if(o.adj_matrix.get(ii, jj) != 0)
					o.adj_matrix(jj,ii) = o.adj_matrix.get(ii,jj);
				else if(o.adj_matrix.get(jj, ii) != 0)
						o.adj_matrix(ii, jj) = o.adj_matrix.get(jj, ii);
			}
		}
		o.adj_matrix.mapToMxArray(outputs[0]);
	}

};


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the 
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB 
function D = vigraAdjacency(inputArray)
function D = vigraAdjacency(inputArray, options);

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
	callMexAllIntFunctor<vigraFunctor>(outputs, inputs);
}
