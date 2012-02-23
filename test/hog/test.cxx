/************************************************************************/
/*                                                                      */
/*    Copyright 2008-2011 by Michael Hanselmann and Ullrich Koethe      */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/
#define HOG_TEST 1

#ifdef _MSC_VER
# pragma warning (disable : 4244)
#endif

#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <vigra/hog3d.hxx>
#include <vigra/multi_distance.hxx>
#include <unittest.hxx>
#include <vector>
#include <limits>
//#include "test_data.hxx"

#include <stdlib.h>


#include "boost/date_time/posix_time/posix_time.hpp"  // for time measurement; clock_t gives only CPU Clock cycles !
#include <boost/thread/thread.hpp> 



using namespace vigra;


class HogTest
{

public:

    HogTest() // load precomputed data
    {}

    void testHog3D()
    {

	/*
	GENERATE TEST DATA:
	
	***Test volume 1:
	Create a 30x30x30 Volume V and initialize with zero. Let Voxel (1,1,1) be an edge voxel in V.
	Apply the distance transform (using the manhatten distance) on volume V. For all non-border points the gradient direction should be the vector (1,1,1), and the squared gradient magnitude should 
	be 3072 everywhere. (using a 3D sobel filter). This can easily be verified with the following MatlabCode:
	
	t = zeros(3,3,3);
	t(1,1,1) = 1;
	t = bwdist(t,'cityblock');
	sobelx(:,:,1) = [1 2 1;2 4 2;1 2 1];      sobelx(:,:,2) = [0 0 0; 0 0 0; 0 0 0];    sobelx(:,:,3) = [-1 -2 -1;-2 -4 -2;-1 -2 -1];
	sobely(:,:,1) = [1 2 1; 0 0 0; -1 -2 -1]; sobely(:,:,2) = [2 4 2; 0 0 0; -2 -4 -2]; sobely(:,:,3) = [1 2 1; 0 0 0; -1 -2 -1];
	sobelz(:,:,1) = [1 0 -1; 2 0 -2; 1 0 -1]; sobelz(:,:,2) = [2 0 -2; 4 0 -4; 2 0 -2]; sobelz(:,:,3) = [1 0 -1;  2 0 -2 ; 1 0 -1];
	a1 = sobelx .* t; a2 = sobely .* t; a3 = sobelz .* t;
	sum(a1(:)).^2 + sum(a2(:)).^2 + sum(a3(:)).^2
	
	*** Test volume 2:
	Same as Test volume 1 except that the edge voxel is (30,30,30), the gradient vector thus becomes (-1,-1,-1)
	
	*** Test volume 3:
	Create a 30x30x30 Volume V and initialize with zero. Let Voxel (15,15,15) be an edge voxel in V.
	Apply the distance transform (using the euclidean distance) on volume V.
	

	PERFORM TESTS:
	
	The following test computes the 3DHOG descriptor on the voxel (15,15,15) with both methods (computation on whole volume, computation on sampling points), 
	using the parameters: 6,12 and 20 bands, 3 filter levels, no normalization, not ignoring the gradient sign.
	The following result is expected and will be tested:
	
	A) Both methods (computing descriptors on the whole volume or only on sampling points) always return the same results.
	
	B1) The first filter level aggregates the gradient magnitudes of 27 voxels (discretized volume of sphere 4/3*pi*r^3 with r = 2^1), and 
	    the second filter   ...                                     224 voxels (...                of sphere 4/3*pi*r^3 with r = 2^2 minus volume of first filter).
	    the third filter   ...                                     1852 voxels (...                of sphere 4/3*pi*r^3 with r = 2^3 minus volume of first and second filter).
	   
	      Angular discretization in 6 bins (faces of a cube) :
		The gradient vector (1,1,1) is equally distant from three of the six surface normals. We expect the magnitude of 27, 224 and 1852 voxels each with a grad.mag of sqrt(3072) to be
		equally distributed into these three bins. I.e. we expect a resulting descriptor of 
		d = (  27*x, 0,   27*x, 0,   27*x, 0, 
		      224*x, 0,  224*x, 0,  224*x, 0, 
		     1852*x, 0, 1852*x, 0, 1852*x, 0), 	with x = sqrt(3072)/3;     
		     
	B2)   Angular discretization in 12 bins (faces of a dodecahedron) :
		The gradient vector (1,1,1) is equally distant from three of the twelve surface normals. We expect the magnitude of 27, 224 and 1852 voxels each with a grad.mag of sqrt(3072) to be
		equally distributed into these three bins. I.e. we expect a resulting descriptor of 		
		d = (0,0,0,0,0,  27*x,   27*x, 0,   27*x, 0,0,0
		     0,0,0,0,0, 224*x,  224*x, 0,  224*x, 0,0,0
		     0,0,0,0,0,1852*x, 1852*x, 0, 1852*x, 0,0,0), 	with x = sqrt(3072)/3;  
		     
	B3)   Angular discretization in 20 bins (faces of a icosahedron) :
		The gradient vector (1,1,1) is parallel to one surface normal. We expect the magnitude of 27, 224 and 1852 voxels each with a grad.mag of sqrt(3072) to be aggregated in one bin.
		I.e. we expect a resulting descriptor of 
		d = (0, 0, 0, 0, 0, 0, 0, 0, 0,   27*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0, 0, 0, 0,  224*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                                                                                                                                                    
		     0, 0, 0, 0, 0, 0, 0, 0, 0, 1852*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 	with x = sqrt(3072);   
		     
	C)  Angular discretization in 6 bins (faces of a cube), ignoring the gradient sign, normalization of every filter level, normalization of the complete descriptor
		We expect the descriptor of B1 with the following changes
			- opposing bins have been combined (ignoring the gradient sign)
			- every filter level is normalized by the pixelcount in the filter
			- the complete descriptor is normalized to give a L2-Norm of 1.
		I.e. we expect the decriptor: 
		d = (27*x / 27,     27*x  / 27,    27*x  / 27, 
		   224*x / 224,    224*x / 224,   224*x / 224, 
		 1852*x / 1852,  1852*x / 1852, 1852*x / 1852)  / || d || ,    with x = sqrt(3072)/3;
		  = 	(1/3, 1/3, 1/3,
			 1/3, 1/3, 1/3,
			 1/3, 1/3, 1/3);

	D) To verify equal border treatment (BORDER_REFLECT_TREATMENT) we compare both methods for all eight corner descriptors (e.g. (1,1,1) and (30,30,30)) . 
		We expect the results of both methods to be equal. Use maximum filter level of 4. (otherwise same parameters as in B1)
		
	E) The descriptor on voxel (15,15,15) on the test volume 1 and test volume 2 should be equal if the gradient sign is ignored. (otherwise same parameters as in B2)
	
	F) The descriptor on voxel (15,15,15) on the test volume 3 should show that the gradient magnitudes are distributed equally into all bins, (same parameters as in B2).
	
	G) Test B1 is repeated for float and double. We expect double and float results to be equal.

	
	*/
	
	// Generate Test volumes
	int mxsize[3] = {200, 100, 100};	
	MultiArrayShape<3>::type shape(mxsize[0], mxsize[1],mxsize[2]); 
	MultiArray<3,double> testvol_1(shape), testvol_2(shape), testvol_3(shape);
	MultiArray<3,float> testvol_1_fl(shape);
	
	for (int xi=0;xi<mxsize[0];xi++)
		for (int yi=0;yi<mxsize[1];yi++)
			for (int zi=0;zi<mxsize[2];zi++)
			{
				testvol_1(xi,yi,zi) = xi+yi+zi;
				testvol_1_fl(xi,yi,zi)=xi+yi+zi;
				
				testvol_2(xi,yi,zi) = (29-xi)+(29-yi)+(29-zi);
			}
			
	testvol_3(15,15,15) = 1;
	separableMultiDistance(srcMultiArrayRange(testvol_3), destMultiArray(testvol_3), true);

	// Define expected results:
	double x = sqrt(3072)/3;
	double resB1[18] = {27*x, 0,   27*x, 0,   27*x, 0, 224*x, 0,  224*x, 0,  224*x, 0, 1852*x, 0, 1852*x, 0, 1852*x, 0};
	double resB2[36] = {0,0,0,0,0,27*x, 27*x, 0, 27*x, 0,0,0, 0,0,0,0,0, 224*x, 224*x, 0, 224*x, 0,0,0, 0,0,0,0,0,1852*x, 1852*x, 0, 1852*x, 0,0,0};
	x = sqrt(3072);
	double resB3[60] = {0, 0, 0, 0, 0, 0, 0, 0, 0,   27*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  224*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1852*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Define equality tolerance
	double eps = 1e-10;
	double eps_dpsp = 1e-4; // for double-single precision comparison
	
	// Test B1 / A
		int filterlevels = 3 ;
		int nbands = 6;
		bool ignore_sign=0;
		bool normalize_descriptor=0;
		bool normalize_sl=0;
		

	boost::posix_time::ptime cl_start,cl_step, cl_end;
	boost::posix_time::time_duration msdiff;
	
		MultiArrayShape<5>::type shape_test1_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res1_full(shape_test1_full);	

cl_start = boost::posix_time::microsec_clock::local_time();

		hogdesc3d(testvol_1,hog_res1_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);

cl_end = boost::posix_time::microsec_clock::local_time();
msdiff = cl_end-cl_start;
std::cout << " Time: " << double(msdiff.total_milliseconds())/1000 << std::endl;
		
		MultiArrayShape<2>::type shape_points(80000, 3); 
		MultiArray<2,int> spoints(shape_points);
		spoints(0,0) = 15; spoints(0,1) = 15; spoints(0,2) = 15;
		MultiArrayShape<3>::type shape_test1_samp(spoints.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res1_samp(shape_test1_samp);	

cl_start = boost::posix_time::microsec_clock::local_time();

		hogdesc3d_samp(testvol_1,hog_res1_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
cl_end = boost::posix_time::microsec_clock::local_time();
msdiff = cl_end-cl_start;
std::cout << " Time: " << double(msdiff.total_milliseconds())/1000 << std::endl;
		
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res1_full(15,15,15,m,k),hog_res1_samp(0,m,k), eps);	
				// TEST B1
				shouldEqualTolerance(hog_res1_full(15,15,15,m,k),resB1[k*nbands + m], eps);	
			}
			
	// Test B2 / A
		nbands = 12;
		MultiArrayShape<5>::type shape_test2_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res2_full(shape_test2_full);	
		hogdesc3d(testvol_1,hog_res2_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<3>::type shape_test2_samp(spoints.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res2_samp(shape_test2_samp);	
		hogdesc3d_samp	(testvol_1,hog_res2_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res2_full(15,15,15,m,k),hog_res2_samp(0,m,k), eps);	
				// TEST B2
				shouldEqualTolerance(hog_res2_full(15,15,15,m,k),resB2[k*nbands + m], eps);	
			}
			
	// Test B3 / A
		nbands = 20;
		MultiArrayShape<5>::type shape_test3_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res3_full(shape_test3_full);	
		hogdesc3d(testvol_1,hog_res3_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<3>::type shape_test3_samp(spoints.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res3_samp(shape_test3_samp);	
		hogdesc3d_samp	(testvol_1,hog_res3_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res3_full(15,15,15,m,k),hog_res3_samp(0,m,k), eps);	
				// TEST B3
				shouldEqualTolerance(hog_res3_full(15,15,15,m,k),resB3[k*nbands + m], eps);	
			}
			
	// Test C / A
	
		nbands = 6;
		int nbands_true=nbands/2; // every two opposing bands are combined
		ignore_sign=1;
		normalize_descriptor=1;
		normalize_sl=1;

		MultiArrayShape<5>::type shape_test4_full(mxsize[0], mxsize[1], mxsize[2],nbands_true,filterlevels); 
		MultiArray<5,double> hog_res4_full(shape_test4_full);	
		hogdesc3d(testvol_1,hog_res4_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<3>::type shape_test4_samp(spoints.shape(0), nbands_true, filterlevels); 
		MultiArray<3,double> hog_res4_samp(shape_test4_samp);	
		hogdesc3d_samp	(testvol_1,hog_res4_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
	
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands_true;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res4_full(15,15,15,m,k),hog_res4_samp(0,m,k), eps);	
				// TEST C
				shouldEqualTolerance(hog_res4_full(15,15,15,m,k),(double)1/3, eps);	
				
			}

	// Test D
		filterlevels = 4;
		nbands = 12;
		ignore_sign=0;
		normalize_descriptor=0;
		normalize_sl=0;
		
		MultiArrayShape<5>::type shape_test5_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res5_full(shape_test5_full);	
		hogdesc3d(testvol_1,hog_res5_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<2>::type shape_points2(8, 3); 
		MultiArray<2,int> spoints2(shape_points2);
		spoints2(0,0) = 0; spoints2(0,1) = 0; spoints2(0,2) = 0;
		spoints2(1,0) = 0; spoints2(1,1) = 0; spoints2(1,2) = 29;
		spoints2(2,0) = 0; spoints2(2,1) = 29; spoints2(2,2) = 0;
		spoints2(3,0) = 0; spoints2(3,1) = 29; spoints2(3,2) = 29;
		spoints2(4,0) = 29; spoints2(4,1) = 0; spoints2(4,2) = 0;
		spoints2(5,0) = 29; spoints2(5,1) = 0; spoints2(5,2) = 29;
		spoints2(6,0) = 29; spoints2(6,1) = 29; spoints2(6,2) = 0;
		spoints2(7,0) = 29; spoints2(7,1) = 29; spoints2(7,2) = 29;
		MultiArrayShape<3>::type shape_test5_samp(spoints2.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res5_samp(shape_test5_samp);	
		hogdesc3d_samp	(testvol_1,hog_res5_samp,spoints2,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);

		for (int p_idx=0;p_idx<spoints2.shape(0);p_idx++)
			for (int k=0;k<filterlevels;k++)
				for (int m=0;m<nbands;m++)
				{
					// TEST D
					shouldEqualTolerance(hog_res5_full(spoints2(p_idx,0),spoints2(p_idx,1),spoints2(p_idx,2),m,k),hog_res5_samp(p_idx,m,k), eps);						
				}

	// Test E
		filterlevels = 3;
		nbands = 12;
		ignore_sign=1;
		nbands_true=nbands/2; // every two opposing bands are combined
		normalize_descriptor=0;
		normalize_sl=0;

		MultiArrayShape<5>::type shape_test6_full(mxsize[0], mxsize[1], mxsize[2],nbands_true,filterlevels); 
		MultiArray<5,double> hog_res6a_full(shape_test6_full),hog_res6b_full(shape_test6_full);

		hogdesc3d(testvol_1,hog_res6a_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		hogdesc3d(testvol_2,hog_res6b_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		
		MultiArrayShape<3>::type shape_test6_samp(spoints.shape(0), nbands_true, filterlevels); 
		MultiArray<3,double> hog_res6a_samp(shape_test6_samp), hog_res6b_samp(shape_test6_samp);	
		hogdesc3d_samp	(testvol_1,hog_res6a_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		hogdesc3d_samp	(testvol_2,hog_res6b_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands_true;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res6a_full(15,15,15,m,k),hog_res6a_samp(0,m,k), eps);
				shouldEqualTolerance(hog_res6b_full(15,15,15,m,k),hog_res6b_samp(0,m,k), eps);	
				// TEST E
				shouldEqualTolerance(hog_res6a_full(15,15,15,m,k),hog_res6b_full(15,15,15,m,k), eps);	
			}
		
		
	// Test F	
		nbands = 12;
		ignore_sign=0;
		normalize_descriptor=0;
		normalize_sl=0;

		MultiArrayShape<5>::type shape_test7_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res7_full(shape_test7_full);	
		hogdesc3d(testvol_3,hog_res7_full,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<3>::type shape_test7_samp(spoints.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res7_samp(shape_test7_samp);	
		hogdesc3d_samp	(testvol_3,hog_res7_samp,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		for (int k=0;k<filterlevels;k++)
		{
			double tempval = hog_res7_full(15,15,15,0,k); // all values in this filter level should be equal;
			for (int m=0;m<nbands;m++)
			{
				// TEST A
				shouldEqualTolerance(hog_res7_full(15,15,15,m,k),hog_res7_samp(0,m,k), eps);	
				// TEST F
				shouldEqualTolerance(hog_res7_full(15,15,15,m,k),tempval, eps);	
			}
		}

	
	//Test G
		filterlevels = 3 ;
		nbands = 6;
		ignore_sign=0;
		normalize_descriptor=0;
		normalize_sl=0;
		
		MultiArrayShape<5>::type shape_test8_full(mxsize[0], mxsize[1], mxsize[2],nbands,filterlevels); 
		MultiArray<5,double> hog_res8_full_double(shape_test8_full);
		MultiArray<5,float>  hog_res8_full_float(shape_test8_full);	
		
		hogdesc3d(testvol_1,   hog_res8_full_double,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		hogdesc3d(testvol_1_fl,hog_res8_full_float ,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		MultiArrayShape<3>::type shape_test8_samp(spoints.shape(0), nbands, filterlevels); 
		MultiArray<3,double> hog_res8_samp_double(shape_test8_samp);
		MultiArray<3,float>  hog_res8_samp_float(shape_test8_samp);	
		
		hogdesc3d_samp(testvol_1,   hog_res8_samp_double,spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		hogdesc3d_samp(testvol_1_fl,hog_res8_samp_float, spoints,filterlevels,nbands, ignore_sign, normalize_descriptor, normalize_sl);
		
		for (int k=0;k<filterlevels;k++)
			for (int m=0;m<nbands;m++)
			{
				// TEST G
				shouldEqualTolerance(hog_res8_full_double(15,15,15,m,k),hog_res8_full_float(15,15,15,m,k), eps_dpsp);
				shouldEqualTolerance(hog_res8_samp_double(0,m,k),hog_res8_samp_float(0,m,k), eps_dpsp);	
			}
	
		
	//printf("%.15g ?= %.15g\n", hog_res4_full(15,15,15,m,k),hog_res4_samp(0,m,k));
    }

};



struct HogTestSuite : public vigra::test_suite
{
    HogTestSuite()
        : vigra::test_suite("HogTestSuite")
    {
        add(testCase(&HogTest::testHog3D));

    }
};


int main (int argc, char ** argv)
{
    HogTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
