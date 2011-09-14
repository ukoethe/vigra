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
#include "test_data.hxx"


#include <stdlib.h>

using namespace vigra;


class UnsupervisedDecompositionTest
{

public:

    UnsupervisedDecompositionTest()
    {}

    void testPCADecomposition()
    {
        unsigned int numComponents = 3;
        unsigned int numFeatures = 159;
        unsigned int numSamples = 1024;

        Matrix<double> features(numFeatures, numSamples, plsaData, ColumnMajor);
        Matrix<double> fz(Shape2(numFeatures, numComponents));
        Matrix<double> zv(Shape2(numComponents, numSamples));

        prepareRows(features, features, ZeroMean);

        principleComponents(features, fz, zv);

        Matrix<double> model = fz*zv;
        shouldEqualTolerance(squaredNorm(model-features), 1530214.34284834, 1e-10);

#if 0
        char hdf5File_2[] = "example_data_results.h5";
        char hdf5group_2[] = "FZ";
        char hdf5group_3[] = "ZV";
        writeHDF5(hdf5File_2, hdf5group_2, fz);
        writeHDF5(hdf5File_2, hdf5group_3, zv);
#endif    
    }

    void testPLSADecomposition()
    {
#if 0 // load data fro, HDF5 file
        char hdf5File[] = "example_data.h5";
        char hdf5group[] = "volume/data";

        HDF5ImportInfo infoHDF5(hdf5File, hdf5group);
        unsigned int numComponents = 3;
        unsigned int numFeatures = infoHDF5.shapeOfDimension(0);
        unsigned int numSamples = infoHDF5.shapeOfDimension(1);

        Matrix<double> features(Shape2(numFeatures, numSamples));
        readHDF5(infoHDF5, features);
#else // get data from header file
        unsigned int numComponents = 3;
        unsigned int numFeatures = 159;
        unsigned int numSamples = 1024;

        Matrix<double> features(numFeatures, numSamples, plsaData, ColumnMajor);
#endif
        Matrix<double> fz(Shape2(numFeatures, numComponents));
        Matrix<double> zv(Shape2(numComponents, numSamples));

        pLSA(features, fz, zv, PLSAOptions().normalizedComponentWeights(false));

        // TESTS
        // -----
        double eps = 1e-10;

        // test if result matrices (approximately) satisfy normalization properties
        Matrix<double> colSumFZ = fz.sum(0);
        for(int i=0; i<columnCount(fz); ++i)
        {
            shouldEqualTolerance( colSumFZ(0,i), 1, eps );
        }
        Matrix<double> colSumZV = zv.sum(0);
        Matrix<double> colSumFeat = features.sum(0);
        for(int i=0; i<columnCount(zv); ++i)
        {
            shouldEqualTolerance( colSumZV(0,i) / colSumFeat(0, i), 1, eps );
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
        Matrix<double> model = fz*zv; 
        double meanError = (features - model).squaredNorm() / columnCount(features);
        should ( meanError < 5000 );

#if 0
        char hdf5File_2[] = "example_data_results.h5";
        char hdf5group_2[] = "FZ";
        char hdf5group_3[] = "ZV";
        writeHDF5(hdf5File_2, hdf5group_2, fz);
        writeHDF5(hdf5File_2, hdf5group_3, zv);
#endif    
    }
};



struct UnsupervisedDecompositionTestSuite : public vigra::test_suite
{
    UnsupervisedDecompositionTestSuite()
        : vigra::test_suite("UnsupervisedDecompositionTestSuite")
    {
        add(testCase(&UnsupervisedDecompositionTest::testPCADecomposition));
        add(testCase(&UnsupervisedDecompositionTest::testPLSADecomposition));
    }
};


int main (int argc, char ** argv)
{
    UnsupervisedDecompositionTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
