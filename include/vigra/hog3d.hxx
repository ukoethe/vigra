#include "multi_array.hxx"
#include "multi_math.hxx"
#include "multi_convolution.hxx"
#include "multi_fft.hxx"

#include <math.h>

/*  3D HOG Descriptor

This function computes a histogram of gradients on a three dimensional input volume. 

Gradients are aggregated in a spherical region around the center point (A), discretized (B) and normalized (C).

A) The function provides the possibility to collect and concatenate gradient histograms for different radii. 
- 1) In the simplest case (filterlevel=1) a histogram is constructed from gradients in the spherical region with a radius of 2^1 around the center point (i.e. 27-point neighborhood). 
- 2) For filterlevel=2, two gradient histograms are concatenated. The first histogram is equal to case (1). The second histogram is constructed from gradients in the spherical region with a radius of 2^2, 
	excluding the gradients that have been used for the first filter.
- 3) For filterlevel=3, three gradient histograms are concatenated. The first two histogram are equal to case (2). The third  ...   with a radius of 2^3, ...  excluding previously used gradients.
- 4) For filterlevel=4, four ....                                 . The first three ...                  case (3). The fourth ...   with a radius of 2^4, ...  excluding previously used gradients.
- 5) For filterlevel=5, five ...                                  . The first four ...                   case (4). The fifth  ...   with a radius of 2^5, ...  excluding previously used gradients.

B) The discretization of gradients (= number of bins within each histogram) can be set to 6, 12 or 20. They correspond to the faces of the polyhedrons: cube, dodecahedron and icosahedron. 
	The gradient magnitudes will be sorted into the bins by determining which surface normal is closest to the gradient vector (using the dot-product). Gradient vectors that lie on the border
	between bins will be split equally into the respective bins.
	The option of ignoring the gradient sign is given. If the gradient sign is ignored opposing faces of the polytope are mapped into the same bin. The number of bins in the histogram is reduced by half.
	
	
C) Two (independent) normalization options can be selected:
- Normalization per filterlevel:  The gradient histogram on every filter level can by normalized by the number pixels in its filter (i.e. the discretized volume of the spherical region used for this level)
- Normalization of the complete descriptor:  The concatenated histograms can be normalized to a L2-norm of 1.

Two function are provided:
1) hogdesc3d(): This computes the descriptors on every pixel in the volume by applying the convolution operation in fourier space. 
The runtime space complexity for a volume with N voxels is O(N * filterlevel * discretization). E.g. for a 200x200x200 double volume (~61 MByte), 12 bins, 3 filter levels, the function uses 2.2 GB of memory to store the output.
This function should only be used if enough memory is available and a dense set of descriptors is needed.
1) hogdesc3d_samp(): This computes the descriptors only on sample points provided by the user avoiding expensive convolutions on the whole volume. 
The function is faster than hogdesc3d() for fewer than ~ 5% sample points of the total voxel count.   The runtime space complexity is O(4*N) to store the gradient directions and magnitudes.

Both functions take input volume and output descriptors as float or double.

USAGE:
hogdesc3d	(input_vol, output_desc,  		int filterlevels, int discretization, bool ignore_sign, bool normalize_descriptor, bool normalize_filterlevel)
hogdesc3d_samp	(input_vol, output_desc, samplepoints, 	int filterlevels, int discretization, bool ignore_sign, bool normalize_descriptor, bool normalize_filterlevel)

PARAMETERS:
input_vol: 		A three dimensional multiarray of size A x B x C.
output_desc:		The resulting descriptor(s). 
				For hogdesc3d this is a 	five dimensional multiarry of size 	A x B x C x (discretization * 1/(1+ignore_sign) ) x filterlevels
				For hogdesc3d_samp this is a	three dimensional multiarray of size 	NoSamplePoints x (discretization * 1/(1+ignore_sign) ) x filterlevels
filterlevels:		Determines the number to filter	levels (1 <= filterlevels <= 5)
discretization:		Determines the binning in the histograms, (allowed values: 6, 12, 20)
ignore_sign:		Determines if the gradient sign is ignored.
normalize_descriptor:	Determines if the whole descriptor is normalized to give a L2-norm of 1.
normalize_filterlevel:	Determines if the gradient histogram on every filter level is normalized by the number pixels in its filter.
			
ONLY for hogdesc3d_samp():
samplepoints: 		This is a two dimensional integer multiarray of size NoSamplePoints x 3 to specify which descriptors have to be calculated. 


*/

namespace vigra
{
	

template< class T>
void getSurfaceNormals(int nbands, T snormals[][3])
{ // Precomputed surface normals on a cube (6 faces), dodecahedron (12 faces), icosahedron (20 faces) 
	if (nbands==6)
	{
snormals[ 0][0] = 0;   snormals[ 0][1] = 1;   snormals[ 0][2] = 0;
snormals[ 1][0] = 0;   snormals[ 1][1] = -1;  snormals[ 1][2] = 0;
snormals[ 2][0] = 1;   snormals[ 2][1] = 0;   snormals[ 2][2] = 0;
snormals[ 3][0] = -1;  snormals[ 3][1] = 0;   snormals[ 3][2] = 0;
snormals[ 4][0] = 0;   snormals[ 4][1] = 0;   snormals[ 4][2] = 1;
snormals[ 5][0] = 0;   snormals[ 5][1] = 0;   snormals[ 5][2] = -1;
	}
	if (nbands==12)
	{
snormals[ 0][0] =  0;                 snormals[ 0][1] = -0.850650808352040; snormals[ 0][2] =  0.525731112119134;
snormals[ 1][0] =  0;                 snormals[ 1][1] =  0.850650808352040; snormals[ 1][2] = -0.525731112119134;
snormals[ 2][0] = -0.525731112119134; snormals[ 2][1] =  0;                 snormals[ 2][2] =  0.850650808352040;
snormals[ 3][0] =  0.525731112119134; snormals[ 3][1] =  0;                 snormals[ 3][2] = -0.850650808352040;
snormals[ 4][0] = -0.850650808352040; snormals[ 4][1] = -0.525731112119134; snormals[ 4][2] =  0;
snormals[ 5][0] =  0.850650808352040; snormals[ 5][1] =  0.525731112119134; snormals[ 5][2] =  0;
snormals[ 6][0] =  0.525731112119134; snormals[ 6][1] =  0;                 snormals[ 6][2] =  0.850650808352040;
snormals[ 7][0] = -0.525731112119134; snormals[ 7][1] =  0;                 snormals[ 7][2] = -0.850650808352040;
snormals[ 8][0] =  0;                 snormals[ 8][1] =  0.850650808352040; snormals[ 8][2] =  0.525731112119134;
snormals[ 9][0] =  0;                 snormals[ 9][1] = -0.850650808352040; snormals[ 9][2] = -0.525731112119134;
snormals[10][0] =  0.850650808352040; snormals[10][1] = -0.525731112119134; snormals[10][2] =  0;
snormals[11][0] = -0.850650808352040; snormals[11][1] =  0.525731112119134; snormals[11][2] =  0;
	}
	if (nbands==20)
	{
snormals[ 0][0] = -0.356822089773090; snormals[ 0][1] = 0;                   snormals[ 0][2] =  0.934172358962716;
snormals[ 1][0] =  0.356822089773090; snormals[ 1][1] = 0;                   snormals[ 1][2] = -0.934172358962716;
snormals[ 2][0] = -0.356822089773090; snormals[ 2][1] = 0;                   snormals[ 2][2] = -0.934172358962716;
snormals[ 3][0] =  0.356822089773090; snormals[ 3][1] = 0;                   snormals[ 3][2] =  0.934172358962716;
snormals[ 4][0] =  0.934172358962716; snormals[ 4][1] = -0.356822089773090;  snormals[ 4][2] =  0;
snormals[ 5][0] = -0.934172358962716; snormals[ 5][1] =  0.356822089773090;  snormals[ 5][2] =  0;
snormals[ 6][0] = -0.934172358962716; snormals[ 6][1] = -0.356822089773090;  snormals[ 6][2] =  0;
snormals[ 7][0] =  0.934172358962716; snormals[ 7][1] =  0.356822089773090;  snormals[ 7][2] =  0;
snormals[ 8][0] = -0.577350269189626; snormals[ 8][1] = -0.577350269189626;  snormals[ 8][2] = -0.577350269189626;
snormals[ 9][0] =  0.577350269189626; snormals[ 9][1] =  0.577350269189626;  snormals[ 9][2] =  0.577350269189626;
snormals[10][0] =  0.577350269189626; snormals[10][1] = -0.577350269189626;  snormals[10][2] =  0.577350269189626;
snormals[11][0] = -0.577350269189626; snormals[11][1] =  0.577350269189626;  snormals[11][2] = -0.577350269189626;
snormals[12][0] =  0;                 snormals[12][1] = -0.934172358962716;  snormals[12][2] = -0.356822089773090;
snormals[13][0] =  0;                 snormals[13][1] =  0.934172358962716;  snormals[13][2] =  0.356822089773090;
snormals[14][0] =  0.577350269189626; snormals[14][1] = -0.577350269189626;  snormals[14][2] = -0.577350269189626;
snormals[15][0] = -0.577350269189626; snormals[15][1] =  0.577350269189626;  snormals[15][2] =  0.577350269189626;
snormals[16][0] =  0;                 snormals[16][1] = -0.934172358962716;  snormals[16][2] =  0.356822089773090;
snormals[17][0] =  0;                 snormals[17][1] =  0.934172358962716;  snormals[17][2] = -0.356822089773090;
snormals[18][0] =  0.577350269189626; snormals[18][1] =  0.577350269189626;  snormals[18][2] = -0.577350269189626;
snormals[19][0] = -0.577350269189626; snormals[19][1] = -0.577350269189626;  snormals[19][2] =  0.577350269189626;
	}
}

template< class T>
void getGradient(MultiArrayView<3,T> vorg, 
		 MultiArrayView<3,T>& voldx, 
		 MultiArrayView<3,T>& voldy, 
		 MultiArrayView<3,T>& voldz, 
		 MultiArrayView<3,T>& volmag)
{
	// get sobel derivative of volume
	Kernel1D<T> triangle, deriv;
	triangle.initExplicitly(-1, 1) = 1.0, 2.0, 1.0;
	deriv.initExplicitly(-1, 1) = 1.0, 0.0, -1.0;
	triangle.setBorderTreatment(BORDER_TREATMENT_REPEAT);
	deriv.setBorderTreatment(BORDER_TREATMENT_REPEAT);

	convolveMultiArrayOneDimension(srcMultiArrayRange(vorg), destMultiArray(voldx), 2, deriv);
	convolveMultiArrayOneDimension(srcMultiArrayRange(voldx), destMultiArray(volmag), 0, triangle);
	convolveMultiArrayOneDimension(srcMultiArrayRange(volmag), destMultiArray(voldx), 1, triangle);
	
	convolveMultiArrayOneDimension(srcMultiArrayRange(vorg), destMultiArray(voldy), 0, deriv);
	convolveMultiArrayOneDimension(srcMultiArrayRange(voldy), destMultiArray(volmag), 1, triangle);
	convolveMultiArrayOneDimension(srcMultiArrayRange(volmag), destMultiArray(voldy), 2, triangle);

	convolveMultiArrayOneDimension(srcMultiArrayRange(vorg), destMultiArray(voldz), 1, deriv);
	convolveMultiArrayOneDimension(srcMultiArrayRange(voldz), destMultiArray(volmag), 0, triangle);
	convolveMultiArrayOneDimension(srcMultiArrayRange(volmag), destMultiArray(voldz), 2, triangle);	
	
	// calculate gradient magnitude
	{ using namespace vigra::multi_math;
	volmag = sqrt(pow(voldx,2) + pow(voldy,2) + pow(voldz,2)); }
}

template< class T>
void hogdesc3d_internal(MultiArrayView<3,T> vorg, 
			MultiArrayView<5,T>& descr, 
			int no_sl, 
			int nbands, 
			bool ignore_sign, 
			bool normalize_descriptor, 
			bool normalize_sl)
{
	
	int sizex = vorg.shape(0);
	int sizey = vorg.shape(1);
	int sizez = vorg.shape(2);
	
	// precomputed filter sizes, normalization constants
	int outer_radius[5] = {2,4,8,16,32}; // 2^k
	int inner_radius[5] = {0,2,4,8,16}; // 2^(k-1), inner radius for the first filter is set to zero to include the center pixel
	int sizef[5] = {3, 7, 15, 31, 63};
	T numb_voxels_in_filter[5] = {27, 224, 1852, 14968, 119988}; // precomputed number of active voxels in filter, used for normalization of every filter level
	
	// check preconditions
	vigra_precondition(no_sl >=1  && no_sl <=5,  "hogdesc3d(): Filter levels have to be in the range of 1 <= x <= 5.");
	vigra_precondition(nbands == 6  || nbands == 12 || nbands == 20,  "hogdesc3d(): Only angular discretization of 6, 12 and 20 bins are allowed.");
	vigra_precondition(sizex > outer_radius[no_sl-1] 
			&& sizey > outer_radius[no_sl-1] 
			&& sizez > outer_radius[no_sl-1], "hogdesc3d(): Input volume is too small for largest filter.");

	// get surface normals
	T snormals[nbands][3];
	getSurfaceNormals(nbands,snormals);
	
	if (ignore_sign==1) 
		nbands = nbands / 2;

	vigra_precondition(descr.shape(0) == sizex 
			&& descr.shape(1) == sizey 
			&& descr.shape(2) == sizez 
			&& descr.shape(3) == nbands 
			&& descr.shape(4) == no_sl, "hogdesc3d(): Output volume has incorrect size.");
	
	// instead of new arrays use memory of descr array  
	MultiArrayView<4, T> descr_temp(descr.bindOuter(0));
	MultiArrayView<3, T> voldx(descr_temp.bindOuter(0)); // the following three assignments are secure because there are always >= 3 bins for the gradient
	MultiArrayView<3, T> voldy(descr_temp.bindOuter(1));	
	MultiArrayView<3, T> voldz(descr_temp.bindOuter(2));

	
	MultiArrayShape<3>::type shape(sizex, sizey, sizez); 
	MultiArray<3,T> volmag(shape);

	getGradient(vorg, voldx, voldy, voldz, volmag);

	MultiArrayShape<4>::type shape_with_bands(sizex, sizey, sizez, nbands);
	MultiArray<4,T> temp4D(shape_with_bands);
	
	// write magnitudes in corresponding bins, ambiguous points are equally split into multiple bins 
	int idx_k[nbands]; // array of indices of possible bins 
	int idx_no; //number of of border faces
	T dotproduct,temp;
	for (int zi = 0; zi < sizez; zi++)
		for (int yi = 0; yi < sizey; yi++)
			for (int xi = 0; xi < sizex; xi++)
			{

				idx_no=0;  
				idx_k[0] = 0;				
				dotproduct = 0;
				
				for (int k = 0; k<nbands; k++)
				{
					if (ignore_sign) 
						temp = std::abs(voldx(xi,yi,zi) * snormals[k*2][0]  + voldy(xi,yi,zi) * snormals[k*2][1] + voldz(xi,yi,zi) * snormals[k*2][2]);
					else
						temp = voldx(xi,yi,zi) * snormals[k][0]  + voldy(xi,yi,zi) * snormals[k][1] + voldz(xi,yi,zi) * snormals[k][2];


					
					if (temp > dotproduct + 1e-04) // if new strongest response found
					{  
						idx_k[0] = k;   // write first element
						idx_no=1;
						dotproduct = temp;
					}
					else if (temp > dotproduct - 1e-04 && temp < dotproduct + 1e-04) // new response equal to the strongest found -> ambiguous point, add bin to list
					{ 
						idx_k[idx_no] = k; // add other borders
						idx_no++;

					}
				}					
				for (int k = 0 ; k < idx_no ; k++)
					temp4D(xi,yi,zi,idx_k[k]) = volmag(xi,yi,zi) / (T)idx_no;  // divide magnitude into several bins if necessary
			} 


	// define convolution filters
	MultiArrayShape<3>::type shape_filter1(sizef[0], sizef[0], sizef[0]);
	MultiArrayShape<3>::type shape_filter2(sizef[1], sizef[1], sizef[1]);
	MultiArrayShape<3>::type shape_filter3(sizef[2], sizef[2], sizef[2]);
	MultiArrayShape<3>::type shape_filter4(sizef[3], sizef[3], sizef[3]);
	MultiArrayShape<3>::type shape_filter5(sizef[4], sizef[4], sizef[4]);

	MultiArray<3, T> filt[5] = {
		MultiArray<3, T>(shape_filter1),
		MultiArray<3, T>(shape_filter2),		
		MultiArray<3, T>(shape_filter3),
		MultiArray<3, T>(shape_filter4),
		MultiArray<3, T>(shape_filter5),
	};
	
	for (int k = 0; k < no_sl; k++)
	{
		int ce = outer_radius[k]-1; // center
		T dist;

		// iterate only over one quadrant -> sufficient for generating of complete filter
		for (int zi = 0; zi < outer_radius[k]; zi++)
			for (int yi = 0; yi < outer_radius[k]; yi++)
				for (int xi = 0; xi < outer_radius[k]; xi++)			
				{
					dist = std::sqrt(std::pow(xi,2) + std::pow(yi,2) + std::pow(zi,2));
					if (dist >= inner_radius[k] && dist < outer_radius[k])
					{
						filt[k](ce + xi, ce + yi, ce + zi) = 1;
						filt[k](ce + xi, ce + yi, ce - zi) = 1;
						filt[k](ce + xi, ce - yi, ce + zi) = 1;
						filt[k](ce + xi, ce - yi, ce - zi) = 1;
						filt[k](ce - xi, ce + yi, ce + zi) = 1;
						filt[k](ce - xi, ce + yi, ce - zi) = 1;
						filt[k](ce - xi, ce - yi, ce + zi) = 1;
						filt[k](ce - xi, ce - yi, ce - zi) = 1;
					}
				}
				

		
		{ using namespace vigra::multi_math;
		if (normalize_sl) // normalization by pixelcount in filter.
			filt[k] = filt[k] / numb_voxels_in_filter[k] ;	}
		
	}
	
	// perform filtering in fourier space
	for (int k = 0; k < no_sl; k++)
	{
		MultiArrayView<4, T> no_sl_view (descr.bindOuter(k));
		for (int m = 0; m < nbands; m++)
		{
			MultiArrayView<3, T> nband_view_dst (no_sl_view.bindOuter(m));
			MultiArrayView<3, T> nband_view_src (temp4D.bindOuter(m));	 
			convolveFFT(nband_view_src, filt[k], nband_view_dst);
		}
	}
	

	// normalize descriptor
	if (normalize_descriptor)
	{
		T norm;

		for (int zi = 0; zi < sizez; zi++)
			for (int yi = 0; yi < sizey; yi++)
				for (int xi = 0; xi < sizex; xi++)
				
				{
					norm=0;
					for(int m=0; m< nbands;m++) 
						for (int k = 0; k < no_sl; k++)
							norm += std::pow(descr(xi,yi,zi,m,k),2);
					norm = std::sqrt(norm);
					for(int m=0; m< nbands;m++) 
						for (int k = 0; k < no_sl; k++)
							descr(xi,yi,zi,m,k) /= norm + 1e-11;
				}	
	}
	
}


template< class T>
void hogdesc3d_samp_internal(MultiArrayView<3,T> vorg, 
			     MultiArrayView<3,T>& descr, 
			     MultiArrayView<2,int> spoints, 
			     int no_sl, 
			     int nbands, 
			     bool ignore_sign, 
			     bool normalize_descriptor, 
			     bool normalize_sl)
{
	int nopoints = spoints.shape(0);
	
	int sizex = vorg.shape(0);
	int sizey = vorg.shape(1);
	int sizez = vorg.shape(2);
		// precomputed filter sizes, normalization constants
	int outer_radius[5] = {2,4,8,16,32}; // 2^k
	int inner_radius[5] = {0,2,4,8,16}; // 2^(k-1), inner radius for the first filter is set to zero to include the center pixel
	int sizef[5] = {3, 7, 15, 31, 63};
	T numb_voxels_in_filter[5] = {27, 224, 1852, 14968, 119988}; // precomputed number of active voxels in filter, used for normalization of every filter level

	// check preconditions
	vigra_precondition(no_sl >=1  && no_sl <=5,  "hogdesc3d_samp(): Filter levels have to be in the range of 1 <= x <= 5.");
	vigra_precondition(nbands == 6  || nbands == 12 || nbands == 20,  "hogdesc3d_samp(): Only angular discretization of 6, 12 and 20 bins are allowed.");
	vigra_precondition(sizex > outer_radius[no_sl-1] 
			&& sizey > outer_radius[no_sl-1] 
			&& sizez > outer_radius[no_sl-1], "hogdesc3d_samp(): Input volume is too small for largest filter.");

	// get surface normals
	T snormals[nbands][3];
	getSurfaceNormals(nbands,snormals);
	
	if (ignore_sign==1) 
		nbands = nbands / 2;
	
	vigra_precondition(nopoints > 0 && nopoints < sizex*sizey*sizez, "hogdesc3d_samp(): Invalid number of sample points.");
	vigra_precondition(spoints.shape(1) == 3, "hogdesc3d_samp(): Wrong dimensionality of sample points.");
	
	for (int i = 0 ; i < nopoints; i++)	
		vigra_precondition(spoints(i,0) >= 0 && spoints(i,0) < sizex 
				&& spoints(i,1) >= 0 && spoints(i,1) < sizey 
				&& spoints(i,2) >= 0 && spoints(i,2) < sizez, "hogdesc3d_samp(): Sample points out of range.");
	vigra_precondition(descr.shape(0) == nopoints 
			&& descr.shape(1) == nbands 
			&& descr.shape(2) == no_sl, "hogdesc3d_samp(): Output volume has incorrect size.");

	
	MultiArrayShape<3>::type shape(sizex, sizey, sizez); 
	MultiArray<3,T> voldx(shape), voldy(shape), voldz(shape), volmag(shape);  // voldx is also later used as array for binning information, voldz is used to store to how many bins this grad beĺongs
	
	getGradient(vorg, voldx, voldy, voldz, volmag);
		
	// write magnitudes in corresponding bins, ambiguous points are equally split into multiple bins 
	int idx_k[nbands]; // array of indices of possible bins 
	int idx_no; //number of of border faces
	T dotproduct,temp;
	for (int zi = 0; zi < sizez; zi++)
		for (int yi = 0; yi < sizey; yi++)
			for (int xi = 0; xi < sizex; xi++)
			{

				idx_no=0;  
				idx_k[0] = 0;				
				dotproduct = 0;
				
				for (int k = 0; k<nbands; k++)
				{
					if (ignore_sign) 
						temp = std::abs(voldx(xi,yi,zi) * snormals[k*2][0]  + voldy(xi,yi,zi) * snormals[k*2][1] + voldz(xi,yi,zi) * snormals[k*2][2]);
					else
						temp = voldx(xi,yi,zi) * snormals[k][0]  + voldy(xi,yi,zi) * snormals[k][1] + voldz(xi,yi,zi) * snormals[k][2];

					if (temp > dotproduct + 1e-04) // if new strongest response found
					{  
						idx_k[0] = k;   // write first element
						idx_no=1;
						dotproduct = temp;
					}
					else if (temp > dotproduct - 1e-04 && temp < dotproduct + 1e-04) // new response equal to the strongest found -> ambiguous point, add bin to list
					{ 
						idx_k[idx_no] = k; // add other borders
						idx_no++;
					}
				}

				// use voldx and voldz as arrays for bin information:
				// (normal case:) if gradient belongs to only one bin, then store the bin number in voldx and set voldz to zero.
				// (border case:) if gradient belongs to more than one bin, store number of used bins in voldx and use voldz as a binary array with a 1 for every active bin
				voldz(xi,yi,zi) = 0.0;
				if (idx_no==1)
					voldx(xi,yi,zi) = idx_k[0];
				else
				{
					voldx(xi,yi,zi) = idx_no;
					for (int k = 0 ; k < idx_no ; k++)

						voldz(xi,yi,zi) += std::pow(2, idx_k[k]);
				}
			} 
		
	// precompute distances in filter
	MultiArrayShape<3>::type shape_filter(sizef[4], sizef[4], sizef[4]);
	MultiArray<3, int> distarr(shape_filter);
	int ce = outer_radius[4]-1;		
	T dist;
		
	// check only for one quadrant -> sufficient for generating of complete filter
	for (int zi = 0; zi < outer_radius[4]; zi++)
		for (int yi = 0; yi < outer_radius[4]; yi++)
			for (int xi = 0; xi < outer_radius[4]; xi++)			
			{
				dist = sqrt(pow(xi,2) + pow(yi,2) + pow(zi,2));
				char val = 0;
				if (dist < outer_radius[4]) val = 5;
				if (dist < outer_radius[3]) val = 4;
				if (dist < outer_radius[2]) val = 3;
				if (dist < outer_radius[1]) val = 2;
				if (dist < outer_radius[0]) val = 1;
				
				distarr(ce + xi, ce + yi, ce + zi) = val;
				distarr(ce + xi, ce + yi, ce - zi) = val;
				distarr(ce + xi, ce - yi, ce + zi) = val;
				distarr(ce + xi, ce - yi, ce - zi) = val;
				distarr(ce - xi, ce + yi, ce + zi) = val;
				distarr(ce - xi, ce + yi, ce - zi) = val;
				distarr(ce - xi, ce - yi, ce + zi) = val;
				distarr(ce - xi, ce - yi, ce - zi) = val;
			}
 
	int p[3];
	int bin=0;
	int shiftfilter = outer_radius[4]-1;
	int d=0;
		
	for(int ip = 0; ip < nopoints; ip++)  // do for every requested descriptor
	{

		// sum up all values of magnitueds within filterrange in corresponding bin. If the filter range exceeds the input array, flip voxels -> simulate BORDER_TREATMENT_REFLECT behaviour
		for(int pz = -outer_radius[no_sl-1]+1; pz < outer_radius[no_sl-1]; pz++) // go over all points that have to be considered for max. filter size
		{
			p[2] = spoints(ip,2) + pz;
			if (p[2] < 0) p[2] = -p[2];
			else if (p[2] >= sizez) p[2] = 2*sizez - p[2] - 2;

			for(int py = -outer_radius[no_sl-1]+1; py < outer_radius[no_sl-1]; py++)
			{
				p[1] = spoints(ip,1) + py;
				if (p[1] < 0) p[1] = -p[1];
				else if (p[1] >= sizey) p[1] = 2*sizey - p[1] - 2;
				
				for(int px = -outer_radius[no_sl-1]+1; px < outer_radius[no_sl-1]; px++)
				{
					p[0] = spoints(ip,0) + px;
					if (p[0] < 0) p[0] = -p[0];
					else if (p[0] >= sizex) p[0] = 2*sizex - p[0] - 2;
					
					d = distarr(px + shiftfilter, py + shiftfilter,pz + shiftfilter);	
					if (d > 0 && d <= no_sl) // part of the descriptor
					{
	
						int bininformation_a = (int)voldz(p[0],p[1],p[2]); // is zero if only one bin active, else holds the information which bins are active
						int bininformation_b = (int)voldx(p[0],p[1],p[2]); // holds the active bin index if only one is active, else holds the number of active bins
						T magnitude = volmag(p[0],p[1],p[2]);

						if (bininformation_a==0) // only one bin active, normal case
						{
							descr(ip, bininformation_b, d-1) +=  magnitude ; // Add magnitude to active bin.
						}
						else  // gradient lies on border between bins
						{
							int idx_b=0;
							magnitude /= (T)bininformation_b;
							while (bininformation_b > 0) // while (bininformation_a > 0)
							{
								if (bininformation_a % 2) // number is odd : rightmost bit is set -> active bin found
								{
									descr(ip, idx_b, d-1) +=  magnitude; // add fraction of magnitude to active bin.
									bininformation_a --;
									bininformation_b --;
								}
								bininformation_a = bininformation_a >> 1;
								idx_b++;
							}
						}
					}
				}
			}
		}
		
		//normalization on filter level
		if (normalize_sl)
		{
			for(int m=0; m< nbands;m++) 
				for (int k = 0; k < no_sl; k++)
					descr(ip, m, k) /= numb_voxels_in_filter[k];
		}
		
	}		
	//normalization of whole descriptor
	if (normalize_descriptor)
	{
		T norm_descr;
		for(int ip = 0; ip < nopoints; ip++)
		{
			norm_descr=0;
			for(int m=0; m< nbands;m++) 
				for (int k = 0; k < no_sl; k++)
					norm_descr += std::pow(descr(ip, m, k),2);
			norm_descr = sqrt(norm_descr);
			for(int m=0; m< nbands;m++) 
				for (int k = 0; k < no_sl; k++)
					descr(ip, m, k) /= norm_descr;
		}
	}
	
}

void hogdesc3d_samp(MultiArrayView<3,double> vorg, MultiArrayView<3,double>& descr, MultiArrayView<2,int> spoints, int no_sl, int nbands, bool ignore_sign, bool normalize_descriptor, bool normalize_sl)
{ hogdesc3d_samp_internal(vorg, descr, spoints, no_sl, nbands, ignore_sign, normalize_descriptor, normalize_sl); }
void hogdesc3d_samp(MultiArrayView<3,float> vorg, MultiArrayView<3,float>& descr, MultiArrayView<2,int> spoints, int no_sl, int nbands, bool ignore_sign, bool normalize_descriptor, bool normalize_sl)
{ hogdesc3d_samp_internal(vorg, descr, spoints, no_sl, nbands, ignore_sign, normalize_descriptor, normalize_sl); }

void hogdesc3d(MultiArrayView<3,double> vorg, MultiArrayView<5,double>& descr,  int no_sl, int nbands, bool ignore_sign, bool normalize_descriptor, bool normalize_sl)
{ hogdesc3d_internal(vorg, descr, no_sl, nbands, ignore_sign, normalize_descriptor, normalize_sl); }
void hogdesc3d(MultiArrayView<3,float> vorg, MultiArrayView<5,float>& descr,  int no_sl, int nbands, bool ignore_sign, bool normalize_descriptor, bool normalize_sl)
{ hogdesc3d_internal(vorg, descr, no_sl, nbands, ignore_sign, normalize_descriptor, normalize_sl); }

}//namespace vigra