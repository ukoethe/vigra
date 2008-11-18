#include <iostream>
#include <string>
#include <vigra/watersheds.hxx>
#include <vigra/watersheds3D.hxx>
#include <vigra/matlab.hxx>
#include <vigra/basicimageview.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/seededregiongrowing3D.hxx>
#include <vigra/contourcirculator.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <vigra/diff2d.hxx>
#include <sstream>
#include <vigra/cellconfigurations.hxx>
#include <set>


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

template <class T>
struct options{
	int numOfDim;
	int conn;
	SRGType crack;
	int method;
	

	BasicImageView<T>  seed;
	MultiArrayView<3,T> seed3D;

	BasicImageView<T>  in;
	MultiArrayView<3,T> in3D;

	BasicImageView<Int32>  out;
	MultiArrayView<3,Int32> out3D;
	
	options(int a, int b,SRGType c, int d){
		numOfDim = a;
		conn = b; 
		crack = c;
		method = d;
	}
};

 //function not used (originally was planned to be used with 3D-Watershed pixels.)
void s3(unsigned char& a, unsigned char&b, unsigned char&c, int shifta, int shiftb, int shiftc){
	a = a|(0x01 << shifta);
	a = b|(0x01 << shiftb);
	a = c|(0x01 << shiftc);
}
void s2(int& a, int&b, int shifta, int shiftb){
	a = a|(0x01 << shifta);
	a = b|(0x01 << shiftb);
}


template <class T>
void vigraWatershed(matlab::OutputArray outputs, matlab::InputArray inputs){
	
	// Constant definition for readibility
	enum {UNION = 0, SEED = 1}meth;
	enum {IMAG = 2, VOLUME = 3}dim;
	enum {FourNeighbor = 4, EightNeighbor = 8, SixNeighbor = 6, TSixNeighbor = 26}nhood;

	//Default Options
	options<T> opt(mxGetNumberOfDimensions(inputs[0]), 
				mxGetNumberOfDimensions(inputs[0])==VOLUME? TSixNeighbor:EightNeighbor,
				CompleteGrow,
				UNION);
	
	//User supplied Options
	if(inputs.isValid(1)){	
		mxArray* seed =mxGetField(inputs[1], 0, "seed");
		mxArray* conn  =mxGetField(inputs[1], 0, "conn");
		mxArray* crack = mxGetField(inputs[1], 0, "crack");
		
		switch(cantorPair( !mxIsChar(seed), opt.numOfDim) ){
			case cP<SEED, IMAG>::value :
				opt.seed = matlab::getImage<T>(seed);
				opt.method = SEED;
				opt.conn = FourNeighbor;
				opt.crack = mxIsNumeric(crack)?
								(bool)mxGetScalar(crack)? 
									  CompleteGrow:KeepContours
								:KeepContours;
				break;
			case cP<SEED, VOLUME>::value:
				opt.seed3D = matlab::getMultiArray<3,T>(seed);
				opt.method = SEED;
				opt.conn = SixNeighbor;
				opt.crack = mxIsNumeric(crack)?
								(bool)mxGetScalar(crack)? 
									 CompleteGrow: KeepContours
								:KeepContours;
				break;
			case cP<UNION, IMAG>::value :{
				int connect = mxIsNumeric(conn)?
							mxGetScalar(conn) : opt.conn;
				if(connect == FourNeighbor|| connect == EightNeighbor)opt.conn = connect;
				else mexWarnMsgTxt("Warning: User supplied connectivity not supported, using default (8)");
				break;
			}
			case cP<UNION, VOLUME>::value:{
				int connect = mxIsNumeric(conn)?
							mxGetScalar(conn) : opt.conn;
				if(connect == SixNeighbor || connect == TSixNeighbor)opt.conn = connect;
				else mexWarnMsgTxt("Warning: User supplied connectivity not supported, using default (26)");
				break;
			}
			default:
				break;
		}
	
	}
	if(opt.numOfDim == IMAG){
		opt.in = matlab::getImage<T>(inputs[0]);
		opt.out = matlab::createImage<Int32>(opt.in.width(), opt.in.height(), outputs[0]);
	}else{
		opt.in3D = matlab::getMultiArray<3, T>(inputs[0]);
		opt.out3D = matlab::createMultiArray<3,Int32>(opt.in3D.shape(), outputs[0]);
	}

					
	{//Preconditions
	if(opt.numOfDim > VOLUME)  
						mexErrMsgTxt("Currently InputArray may only have 2 or 3 dimensions");
	if(inputs.isEmpty(0)) 
						mexErrMsgTxt("Input Image is empty!");
	
	if(		opt.method == SEED 
		&& 	opt.numOfDim == IMAG 
		&&	(	opt.seed.width() != opt.in.width() 
			|| 	opt.seed.height() != opt.in.height()	) )
						mexErrMsgTxt("Seed Matrix and Heightmap must have same size");
								
	if(		opt.method == SEED 
		&& 	opt.numOfDim == VOLUME 
		&&	(opt.seed3D.shape() != opt.in3D.shape() ) )
						mexErrMsgTxt("Seed Matrix and Heightmap must have same size");	
	}
	
	
	int max_region_label = -1;
	// contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 
	switch(cantorPair(opt.method, opt.numOfDim, opt.conn)){
		//cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
		//Code is basically the code on the VIGRA-reference page 
		case cP3<UNION, IMAG, EightNeighbor>::value:	
			max_region_label = watersheds(srcImageRange(opt.in), destImage(opt.out));
			break;
		case cP3<UNION, IMAG, FourNeighbor>::value:
			max_region_label = watersheds(srcImageRange(opt.in), destImage(opt.out), FourNeighborCode());
			break;
		case cP3<UNION, VOLUME, TSixNeighbor>::value:
			max_region_label = watersheds3DTwentySix(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D));
			break;
		case cP3<UNION, VOLUME, SixNeighbor>::value:
			max_region_label = watersheds3DSix(srcMultiArrayRange(opt.in3D), destMultiArray(opt.out3D));
			break;
		case cP3<SEED, IMAG, FourNeighbor>::value:
		{				
			//find maximimum of seed Image
			FindMinMax<T> minmax;   // init functor
			inspectImage(srcImageRange(opt.seed), minmax);
			max_region_label = minmax.max;
			
			ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
											gradstat(max_region_label);
											
			seededRegionGrowing(	srcImageRange(opt.in), 
									srcImage(opt.seed),
									destImage(opt.out), 
									gradstat,opt.crack );
			break;
		}
		case cP3<SEED, VOLUME, SixNeighbor>::value:
		{
			FindMinMax<T> minmax;
			inspectMultiArray(srcMultiArrayRange(opt.seed3D), minmax);
			max_region_label = minmax.max;
			
			ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
											gradstat(max_region_label);
											
			seededRegionGrowing3D(	srcMultiArrayRange(opt.in3D),
									srcMultiArray(opt.seed3D),
									destMultiArray(opt.out3D),
									gradstat, opt.crack);
			break;
		}		
		default:
			mexErrMsgTxt("Something went wrong");
	}	

	//adjacency matrix requested
	if(outputs.isValid(1)){
		BasicImageView<Int32>  adj_matrix = matlab::createImage<Int32>(max_region_label, max_region_label, outputs[1]);
		ImageIterator<Int32> adjupleft = adj_matrix.upperLeft();
		switch(cantorPair(opt.numOfDim, (int)(opt.crack == KeepContours))){
			//in case <XX, 0 (crackedges)> the Pixel/Voxel underneath, on the right (and behind) are considered for neighborhood.
			//in case<IMAG, 1>  the cell configuration around current pixel is compared to the ones in cellconfigurations.h (Configuration are Freeman coded (see PhD thesis, Chap 5.7)
			// if current configuration is of Type Line then adjacent pixels are considered as neighbors.
			// in case<VOLUME, 1> the number of adjacent regions are counted. the regions are considered adjacent if only two such regions exist
			case cP<IMAG, 0>::value:
			{

				ImageIterator<Int32> upleft = opt.out.upperLeft();
				ImageIterator<Int32> downright = opt.out.lowerRight();
				Diff2D size = downright - upleft;
				for(int y = 0; y < size.y-1; ++y){
					for(int x = 0; x< size.x-1; ++x){
						adj_matrix(upleft(x,y)-1, upleft(x+1, y)-1)++;
						adj_matrix(upleft(x,y)-1, upleft(x, y+1)-1)++;
						adj_matrix(upleft(x+1, y)-1, upleft(x,y)-1)++;
						adj_matrix(upleft(x, y+1)-1, upleft(x,y)-1)++;
					}
				}
				
				for(int x = 0; x < size.x-1; ++x){
					adjupleft(upleft(x, size.y-1)-1, upleft(x+1, size.y-1)-1)++;
					adjupleft(upleft(x+1, size.y-1)-1, upleft(x, size.y-1)-1)++;
				}
				for(int y = 0; y < size.y-1; ++y){
					adjupleft(upleft(size.x-1, y)-1, upleft(size.x-1, y+1)-1)++;
					adjupleft(upleft(size.x-1, y+1)-1, upleft(size.x-1, y)-1)++;
				}
				break;
			}
			case cP<VOLUME, 0>::value:
			{
				MultiArrayShape<3>::type sze = opt.out3D.shape();
				for(int ii = 0; ii < sze[0]-1; ii++){
					for(int jj = 0; jj < sze[1]-1; jj++){
						for(int kk = 0; kk < sze[2]-1; kk++){
							adj_matrix(opt.out3D(ii,jj,kk)-1, opt.out3D(ii+1,jj,kk)-1)++;
							adj_matrix(opt.out3D(ii,jj,kk)-1, opt.out3D(ii,jj+1,kk)-1)++;
							adj_matrix(opt.out3D(ii,jj,kk)-1, opt.out3D(ii,jj,kk+1)-1)++;
							
							adj_matrix(opt.out3D(ii+1,jj,kk)-1, opt.out3D(ii,jj,kk)-1)++;
							adj_matrix(opt.out3D(ii,jj+1,kk)-1, opt.out3D(ii,jj,kk)-1)++;
							adj_matrix(opt.out3D(ii,jj,kk+1)-1, opt.out3D(ii,jj,kk)-1)++;
						}
					}
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					for(int kk = 0; kk < sze[2]-1; kk++){
						adj_matrix(opt.out3D(sze[0]-1,jj,kk)-1, opt.out3D(sze[0]-1,jj+1,kk)-1)++;
						adj_matrix(opt.out3D(sze[0]-1,jj,kk)-1, opt.out3D(sze[0]-1,jj,kk+1)-1)++;

						adj_matrix(opt.out3D(sze[0]-1,jj+1,kk)-1, opt.out3D(sze[0]-1,jj,kk)-1)++;
						adj_matrix(opt.out3D(sze[0]-1,jj,kk+1)-1, opt.out3D(sze[0]-1,jj,kk)-1)++;						
					}
				}
				for(int ii = 0; ii < sze[0]-1; ii++){
					for(int kk = 0; kk < sze[2]-1; kk++){
						adj_matrix(opt.out3D(ii,sze[1]-1,kk)-1, opt.out3D(ii+1,sze[1]-1,kk)-1)++;
						adj_matrix(opt.out3D(ii,sze[1]-1,kk)-1, opt.out3D(ii,sze[1]-1,kk+1)-1)++;	

						adj_matrix(opt.out3D(ii+1,sze[1]-1,kk)-1, opt.out3D(ii,sze[1]-1,kk)-1)++;
						adj_matrix(opt.out3D(ii,sze[1]-1,kk+1)-1, opt.out3D(ii,sze[1]-1,kk)-1)++;							
					}
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					for(int ii = 0; ii <sze[0] -1; ii++){
						adj_matrix(opt.out3D(ii,jj,sze[2]-1)-1, opt.out3D(ii,jj+1,sze[2]-1)-1)++;
						adj_matrix(opt.out3D(ii,jj,sze[2]-1)-1, opt.out3D(ii+1,jj,sze[2]-1)-1)++;
						
						adj_matrix(opt.out3D(ii,jj+1,sze[2]-1)-1, opt.out3D(ii,jj,sze[2]-1)-1)++;
						adj_matrix(opt.out3D(ii+1,jj,sze[2]-1)-1, opt.out3D(ii,jj,sze[2]-1)-1)++;
					}
				}
				
				for(int ii = 0; ii < sze[0]-1; ii++){
					adj_matrix(opt.out3D(ii,sze[1]-1,sze[2]-1)-1, opt.out3D(ii+1,sze[1]-1,sze[2]-1)-1)++;
					adj_matrix(opt.out3D(ii+1,sze[1]-1,sze[2]-1)-1, opt.out3D(ii,sze[1]-1,sze[2]-1)-1)++;
				}
				for(int jj = 0; jj < sze[1]-1; jj++){
					adj_matrix(opt.out3D(sze[0]-1,jj,sze[2]-1)-1, opt.out3D(sze[0]-1,jj+1,sze[2]-1)-1)++;
					adj_matrix(opt.out3D(sze[0]-1,jj+1,sze[2]-1)-1, opt.out3D(sze[0]-1,jj,sze[2]-1)-1)++;
				}
				for(int kk = 0; kk < sze[2]-1; kk++){
					adj_matrix(opt.out3D(sze[0]-1,sze[1]-1,kk)-1, opt.out3D(sze[0]-1,sze[1]-1,kk+1)-1)++;
					adj_matrix(opt.out3D(sze[0]-1,sze[1]-1,kk+1)-1, opt.out3D(sze[0]-1,sze[1]-1,kk)-1)++;
				}
				break;
			}
			case cP<IMAG, 1>::value:
			{

				ImageIterator<Int32> upleft = opt.out.upperLeft();
				ImageIterator<Int32> downright = opt.out.lowerRight();
				Diff2D size = downright - upleft;
				for(int y = 1; y < size.y-1; ++y){
					for(int x = 1; x< size.x-1; ++x){

					if(upleft(x,y) == 0){
							NeighborhoodCirculator<BasicImageView<Int32>::Iterator, EightNeighborCode>
																circulator(upleft+ Diff2D(x, y));
							//circulator += 2;
							NeighborhoodCirculator<BasicImageView<Int32>::Iterator, EightNeighborCode>
																end(circulator);
							unsigned char BitField = 0;
							int regions[4] = {-1, -1, -1, -1};
							int temp = 0;
							int* visited = new int[max_region_label];
							for(int ii = 0; ii < max_region_label; ii++){
								visited[ii] = 0;
							}
							
							do{
								if(*circulator == 0){
									BitField = BitField >> 1;
									BitField = BitField | 0x80;
								}else{
									if(visited[*circulator] == 0){
										regions[temp] = *circulator;
										temp++;
										visited[*circulator] = 1;
									}
									BitField = BitField >> 1;
									BitField = BitField | 0x00;
								}
							}while(++circulator != end);
								std::stringstream oi;
								std::string au = "";
								std::string bl = "";
								//oi << regions[0] << regions[1] << regions[2] << regions[3];
								//mexWarnMsgTxt(oi.str().c_str());
							if(cellimage::cellConfigurations[(int)BitField] == cellimage::CellTypeLine){
								for(int ii = 0; ii < 3; ii ++){
									for(int jj = ii+1; jj < 4; jj++){
										if(regions[ii] != -1 && regions[jj] != -1){
											adjupleft(regions[ii]-1, regions[jj]-1)++;
										}
									}
								}
							}
							delete [] visited;
						}
					}
				}				
				break;
			}
			case cP<VOLUME, 1>::value:
			{
				MultiArrayShape<3>::type sze = opt.out3D.shape();
				for(int ii = 1; ii < sze[0]-1; ii++){
					for(int jj = 1; jj < sze[1]-1; jj++){
						for(int kk = 1; kk < sze[2]-1; kk++){
							if(opt.out3D(ii,jj,kk) == 0){
								
								unsigned char x[4] = {0, 0 ,0, 0}, y[4] = {0, 0, 0, 0}, z[4] = {0, 0, 0, 0};
								std::set<int> regions;
								
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
								
								int s[26*3] = 	{	0,	7,	6,
												1,	7,	7,
												2,	6,	7,
												3,	5,	7,
												4,	5,	6,
												5,	5,	5,
												6,	6,	5,
												7,	7,	5,
												0,	0,	8,
												1,	0,	0,
												2,	0,	8,
												3,	4,	0,
												4,	4,	8,
												5,	4,	4,
												6,	4,	8,
												7,	0,	4,
												0,	1,	2,
												1,	1,	1,
												2,	2,	1,
												3,	3,	1,
												4,	3,	2,
												4,	3,	3,
												6,	2,	3,
												7,	1,	3,
												2,	2,	8,
												6,	6,	8};
								
								int ind[26*3] = {	1,	2,	3,
												1,	3,	3,
												1,	3,	2,
												1,	3,	1,
												1,	2,	1,
												1,	1,	1,
												1,	1,	2,
												1,	1,	3,
												2,	2,	0,
												2,	3,	3,
												2,	2,	0,
												2,	3,	1,
												2,	2,	0,
												2,	1,	1,
												2,	2,	0,
												2,	1,	3,
												3,	2,	3,
												3,	3,	3,
												3,	3,	2,
												3,	3,	1,
												3,	2,	1,
												3,	1,	1,
												3,	1,	2,
												3,	1,	3};
								
								regions.clear();
								
								for(int ll = 0; ll < 26; ll++){
									if(opt.out3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2]) == 0 ){
												s3	(
														x[ind[ll*3]], 	y[ind[ll*3+1]], z[ind[ll*3+2]],
														s[ll*3],		s[ll*3+1],		s[ll*3+2]
														);
									}
									else{
										regions.insert((opt.out3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2])));
									}
								}
								//mexWarnMsgTxt("   ");
								if(cellimage::cellConfigurations[(int)x[1]] == cellimage::CellTypeLine)x[1] = 1;
								else x[1] = 0;
								if(cellimage::cellConfigurations[(int)x[2]] == cellimage::CellTypeLine)x[2] = 1;
								else x[2] = 0;
								if(cellimage::cellConfigurations[(int)x[3]] == cellimage::CellTypeLine)x[3] = 1;
								else x[3] = 0;
								if(cellimage::cellConfigurations[(int)y[1]] == cellimage::CellTypeLine)y[1] = 1;
								else y[1] = 0;
								if(cellimage::cellConfigurations[(int)y[2]] == cellimage::CellTypeLine)y[2] = 1;
								else y[2] = 0;
								if(cellimage::cellConfigurations[(int)y[3]] == cellimage::CellTypeLine)y[3] = 1;
								else y[3] = 0;
								if(cellimage::cellConfigurations[(int)z[1]] == cellimage::CellTypeLine)z[1] = 1;
								else z[1] = 0;
								if(cellimage::cellConfigurations[(int)z[2]] == cellimage::CellTypeLine)z[2] = 1;
								else z[2] = 0;
								if(cellimage::cellConfigurations[(int)z[3]] == cellimage::CellTypeLine)z[3] = 1;
								else z[3] = 0;
								
								if(x[1] == 1 && x[2] == 1 && x[3] == 1)x[1] = 1;
									else x[1] = 0;
								if(y[1] == 1 && y[2] == 1 && y[3] == 1)y[1] = 1;
									else y[1] = 0;
								if(z[1] == 1 && z[2] == 1 && z[3] == 1)z[1] = 1;
									else z[1] = 0;
									
								if(regions.size() == 2 /*true /*(x[1] == 1 && y[1] == 1)||(x[1] == 1 && z[1] == 1)||(y[1] == 1 && z[1] == 1)*/){
									std::set<int>::const_iterator iter = regions.begin();
									adjupleft(*iter-1,*(iter++)-1)++;
									/*
									for( 	std::set<int>::const_iterator iter = regions.begin();
											iter != regions.end();
											++iter ) {
										for( 	std::set<int>::const_iterator iter2 = regions.begin();
												iter2 != regions.end();
												++iter2 ) {
											if(*iter != *iter2){
												
											}
										}								

									}*/

								}
							}
						}
					}
				}
				
				break;
			}
		}
		
		//complete adj_matrix
		for(int ii = 0; ii < max_region_label; ii++){
			for(int jj = ii; jj < max_region_label; jj++){
				if(ii == jj)adjupleft(ii,ii) = 1;
				else if(adjupleft(ii, jj) != 0)
					adjupleft(jj,ii) = adjupleft(ii,jj);
				else if(adjupleft(jj, ii) != 0)
						adjupleft(ii, jj) = adjupleft(jj, ii);
			}
		}
	}
	

}

/** MATLAB 
function L = vigraWatershedb(inputArray)
function L = vigraWatershedb(inputArray, options);
function [L A] = vigraWatershedb(..);

L = vigraWatershed(inputArray) computes a label matrix identifying the watershed regions of the inputArray , which may be 2 or 3 dimensional. The elements of L are Int32 values  greater than 0. 
There are no watershedpixels. Do not use this  with Arrays of type other than single or double.
L = vigraWatershed(inputImage, options)  does the same with user options.
options is a struct with possible fields: "seed", "conn" and "crack

"seed": either a matrix of the size of input (seed) or "none", if "none" or field not supplied the union find algorithm will be used. in this case "conn" will change the connectivity settings
(4 or 8 in 2D and 6 ord 26 in 3D images - 8 and 26 are default for union find)
"conn": see above. only used if seed not supplied
"crack": only used if seeded region growing is  being used. crack = 1 means crack edge borders, crack = 0 means the watershed pixels are left. 
L A = vigraWatershed(...) 	additionally produces an adjacency matrix. If watershed pixels exist then the value in the matrix is the number of watershedpixels, that truely connect
					two regions. If watershed pixels do not exists then the value corresponds to the number of cracks between two regions.

L is of type Int32

*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
	mxClassID inClass = mxGetClassID(inputs[0]);
	switch(inClass){
		case mxDOUBLE_CLASS:
			vigraWatershed<double>(outputs, inputs);	break;
		case mxSINGLE_CLASS:
			vigraWatershed<float>(outputs, inputs);		break;
        case mxINT8_CLASS:
			vigraWatershed<Int8>(outputs, inputs);		break;
		case mxINT16_CLASS:
			vigraWatershed<Int16>(outputs, inputs);		break;
		case mxINT32_CLASS:
			vigraWatershed<Int32>(outputs, inputs);		break;
		case mxINT64_CLASS:
			vigraWatershed<Int64>(outputs, inputs);		break;
        case mxUINT8_CLASS:
			vigraWatershed<UInt8>(outputs, inputs);		break;
		case mxUINT16_CLASS:
			vigraWatershed<UInt16>(outputs, inputs);	break;
		case mxUINT32_CLASS:
			vigraWatershed<UInt32>(outputs, inputs);	break;
		case mxUINT64_CLASS:
			vigraWatershed<UInt64>(outputs, inputs);	break;		
		default:
			mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");
	}
/*
    if(matlab::ValueType<UInt8>::check(inputs[0]))
        vigraWatershed<UInt8>(outputs, inputs);
    else if(matlab::ValueType<float>::check(inputs[0]))
        vigraWatershed<float>(outputs, inputs);
    else if(matlab::ValueType<double>::check(inputs[0]))
        vigraWatershed<double>(outputs, inputs);
    else if(matlab::ValueType<UInt16>::check(inputs[0]))
        vigraWatershed<UInt16>(outputs, inputs);
    else if(matlab::ValueType<UInt32>::check(inputs[0]))
        vigraWatershed<UInt32>(outputs, inputs);
    else if(matlab::ValueType<UInt64>::check(inputs[0]))
        vigraWatershed<UInt64>(outputs, inputs);
    else if(matlab::ValueType<Int8>::check(inputs[0]))
        vigraWatershed<Int8>(outputs, inputs);
    else if(matlab::ValueType<Int16>::check(inputs[0]))
        vigraWatershed<Int16>(outputs, inputs);
    else if(matlab::ValueType<Int32>::check(inputs[0]))
        vigraWatershed<Int32>(outputs, inputs);
    else if(matlab::ValueType<Int64>::check(inputs[0]))
        vigraWatershed<Int64>(outputs, inputs);
	else
        mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");*/
}
