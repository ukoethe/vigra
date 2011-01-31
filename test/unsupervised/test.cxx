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
#define UNSUPERVISED_TEST 1

#ifdef _MSC_VER
# pragma warning (disable : 4244)
#endif

#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <vigra/unsupervised_decomposition.hxx>
#include <unittest.hxx>
#include <vector>
#include <limits>
#include "vigra/hdf5impex.hxx"
#include "vigra/multi_array.hxx"


#include <stdlib.h>

using namespace vigra;


class PLSATest
{

public:

    PLSATest()
    {}

	void testPLSADecomposition()
	{
		char hdf5File[] = "example_data.h5";
		char hdf5group_1[] = "volume/data";

		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,double> features(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        readHDF5(infoHDF5_1, features);

		unsigned int numComponents = 3;
		unsigned int numFeatures = infoHDF5_1.shapeOfDimension(0);
		unsigned int numVoxels = infoHDF5_1.shapeOfDimension(1);

		MultiArray<2,double> FZ(MultiArrayShape<2>::type(numFeatures, numComponents));
		MultiArray<2,double> ZV(MultiArrayShape<2>::type(numComponents, numVoxels));

		PLSA plsa(vigra::PLSAOptions().numberOfComponents(numComponents));
		plsa.decompose(features, FZ, ZV);

		// TESTS
		// -----
		Matrix<double> fz = FZ;
		Matrix<double> zv = ZV;
		Matrix<double> feats = features;
		double eps = 1e-10;

		// test if result matrices have correct dimensionality
		should (rowCount(fz) == rowCount(feats));
		should (columnCount(fz) == numComponents);
		should (rowCount(zv) == numComponents);
		should (columnCount(zv) == columnCount(feats));

		// test if result matrices (approximately) satisfy normalization properties
		Matrix<double> colSumFZ = fz.sum(0);
		for(int i=0; i<columnCount(fz); ++i)
		{
			shouldEqualTolerance( colSumFZ(i, 0), 1, eps );
		}
		Matrix<double> colSumZV = zv.sum(0);
		for(int i=0; i<columnCount(zv); ++i)
		{
			shouldEqualTolerance( colSumZV(i, 0), 1, eps );
		}

		// all entries in FZ, ZV are >= 0
		for(int j=0; j<columnCount(zv); ++j)
		{
			for(int i=0; i<rowCount(zv); ++i)
			{
				should ( zv(i, j) >= 0 );
			}
		}
		for(int j=0; j<columnCount(fz); ++j)
		{
			for(int i=0; i<rowCount(fz); ++i)
			{
				should ( fz(i, j) >= 0 );
			}
		}

		// test if reconstruction is close to original
		// tricky - how to properly test that? it will never be identical!
        Matrix<double> fzv = fz*zv; 
	    Matrix<double> voxelSums = feats.sum(0);
		Matrix<double> ones(numFeatures, 1, 1);
        Matrix<double> model = (pmul((ones * voxelSums), fzv));
		double meanError = (feats - model).squaredNorm() / columnCount(feats);
		should ( meanError < 5000 );

		/*
		char hdf5File_2[] = "example_data_results.h5";
		char hdf5group_2[] = "FZ";
		char hdf5group_3[] = "ZV";
		writeHDF5(hdf5File_2, hdf5group_2, FZ);
		writeHDF5(hdf5File_2, hdf5group_3, ZV);
		*/

	}

};



struct PLSATestSuite : public vigra::test_suite
{
    PLSATestSuite()
        : vigra::test_suite("PLSATestSuite")
    {
		add(testCase(&PLSATest::testPLSADecomposition));
	}
};


int main (int argc, char ** argv)
{
    PLSATestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

	return failed != 0;
}
