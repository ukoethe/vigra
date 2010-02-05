/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Michael Hanselmann and Ullrich Koethe        */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "vigra/stdimage.hxx"
#include "unittest.hxx"
#include "vigra/hdf5impex.hxx"
#include "vigra/multi_array.hxx"

using namespace vigra;


class HDF5ExportImportTest
{

public:

    HDF5ExportImportTest()
    {}


	void testUnstridedHDF5ExportImport_rowMajor()
	{
		testUnstridedHDF5ExportImport(true);
	}

	void testUnstridedHDF5ExportImport_columnMajor()
	{
		testUnstridedHDF5ExportImport(false);
	}

	void testUnstridedHDF5ExportImport(bool rowMajorOrder)
	{
		// export and import data from and to unstrided array

		char hdf5File[] = "testfile1.hdf5";
        std::remove(hdf5File);

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
		{
            out_data_2.data () [i] = i;
		}

		// export
		// ...data 1
		char hdf5group_1[] = "group/subgroup/subsubgroup/data1";
		writeToHDF5File(hdf5File, hdf5group_1, out_data_1, rowMajorOrder);
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data2";
		writeToHDF5File(hdf5File, hdf5group_2, out_data_2, rowMajorOrder);

		// import
		// ...data 1
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,int> in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_1, in_data_1, rowMajorOrder);
		// ...data 2
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2, rowMajorOrder);

		// compare content
		// ...data 1
		should (in_data_1 == out_data_1);
		// ...data 2
		should (in_data_2 == out_data_2);
	}

	void testStridedHDF5ExportImport1_rowMajor()
	{
		testStridedHDF5ExportImport1(true);
	}

	void testStridedHDF5ExportImport1_columnMajor()
	{
		testStridedHDF5ExportImport1(false);
	}

	void testStridedHDF5ExportImport1(bool rowMajorOrder)
	{
		// export data from strided arrays and import to unstrided arrays

		char hdf5File[] = "testfile2.hdf5";
        std::remove(hdf5File);

		// int data in 2 dimensions (partly negative)
		MultiArray<3,int> out_data_3(MultiArrayShape<3>::type(2, 3, 4));
        // initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_3.data () [i] = i;
		//std::cout << "Test (0,0,0), (1,0,0), (0,1,0), (0,0,1): " << out_data_3(0,0,0) << " " << out_data_3(1,0,0) << " " << out_data_3(0,1,0) << " " << out_data_3(0,0,1) << std::endl;
		// bind inner dimension to test if strided data impex works
		MultiArrayView<2,int,StridedArrayTag> out_data_4(out_data_3.bindInner(1));
		// bind dimension in the middle to test if strided data impex works
		MultiArrayView<2,int,StridedArrayTag> out_data_5(out_data_3.bindAt(1, 1));

		// export the two sets
		char hdf5group_4[] = "group/subgroup/data4";
		writeToHDF5File(hdf5File, hdf5group_4, out_data_4, rowMajorOrder);
		char hdf5group_5[] = "group/subgroup/data5";
		writeToHDF5File(hdf5File, hdf5group_5, out_data_5, rowMajorOrder);

		// import test: copy data to unstrided array
		HDF5ImportInfo infoHDF5_4(hdf5File, hdf5group_4);
        MultiArray<2,int> in_data_4(MultiArrayShape<2>::type(infoHDF5_4.shapeOfDimension(0), infoHDF5_4.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_4, in_data_4, rowMajorOrder);
		HDF5ImportInfo infoHDF5_5(hdf5File, hdf5group_5);
        MultiArray<2,int> in_data_5(MultiArrayShape<2>::type(infoHDF5_5.shapeOfDimension(0), infoHDF5_5.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_5, in_data_5, rowMajorOrder);
		// compare content
		should (in_data_4 == out_data_4);
		should (in_data_5 == out_data_5);
	}

	void testStridedHDF5ExportImport2_rowMajor()
	{
		testStridedHDF5ExportImport2(true);
	}

	void testStridedHDF5ExportImport2_columnMajor()
	{
		testStridedHDF5ExportImport2(false);
	}

	void testStridedHDF5ExportImport2(bool rowMajorOrder)
	{
		// export data from unstrided arrays and import to strided arrays

		char hdf5File[] = "testfile3.hdf5";
        std::remove(hdf5File);

		// int data in 2 dimensions (partly negative)
		MultiArray<3,int> out_data_6(MultiArrayShape<3>::type(2, 3, 4));
        // initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_6.data () [i] = i;

		// export the two sets
		char hdf5group_6[] = "group/subgroup/data6";
		writeToHDF5File(hdf5File, hdf5group_6, out_data_6, rowMajorOrder);

		// import test: copy data to strided array
		HDF5ImportInfo infoHDF5_6(hdf5File, hdf5group_6);
        MultiArray<4,int> in_data_6a(MultiArrayShape<4>::type(3, infoHDF5_6.shapeOfDimension(0), infoHDF5_6.shapeOfDimension(1), infoHDF5_6.shapeOfDimension(2)));
		in_data_6a.init(42);
		MultiArrayView<3,int,StridedArrayTag> in_data_6b(in_data_6a.bindInner(0));
        loadFromHDF5File(infoHDF5_6, in_data_6b, rowMajorOrder);
		// compare content
		should (in_data_6b == out_data_6);
	}

	void testAppendNewDataToHDF5_rowMajor()
	{
		testAppendNewDataToHDF5(true);
	}

	void testAppendNewDataToHDF5_columnMajor()
	{
		testAppendNewDataToHDF5(false);
	}

	void testAppendNewDataToHDF5(bool rowMajorOrder)
	{
		// write data to file, close file, open file, write more data, compare

		char hdf5File[] = "testfile4.hdf5";
        std::remove(hdf5File);

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		// ...data 1
		char hdf5group_1[] = "group/subgroup/subsubgroup/data1";
		writeToHDF5File(hdf5File, hdf5group_1, out_data_1, rowMajorOrder);

		// append to existing file
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data2";
		writeToHDF5File(hdf5File, hdf5group_2, out_data_2, rowMajorOrder);

		// import
		// ...data 1
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,int> in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_1, in_data_1, rowMajorOrder);
		// ...data 2
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2, rowMajorOrder);

		// compare content
		// ...data 1
		should (in_data_1 == out_data_1);
		// ...data 2
		should (in_data_2 == out_data_2);
	}

	void testOverwriteExistingDataInHDF5_rowMajor()
	{
		testOverwriteExistingDataInHDF5(true);
	}

	void testOverwriteExistingDataInHDF5_columnMajor()
	{
		testOverwriteExistingDataInHDF5(false);
	}

	void testOverwriteExistingDataInHDF5(bool rowMajorOrder)
	{
		// write data to file, close file, open file, overwrite data, compare

		char hdf5File[] = "testfile5.hdf5";
        std::remove(hdf5File);

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		// ...data 1
		char hdf5group_1[] = "group/subgroup/data1";
		writeToHDF5File(hdf5File, hdf5group_1, out_data_1, rowMajorOrder);

		// overwrite in existing file
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data1";
		writeToHDF5File(hdf5File, hdf5group_1, out_data_2, rowMajorOrder);

		// import
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2, rowMajorOrder);

		// compare content
		should (in_data_2 == out_data_2);
	}

};

struct HDF5ImportExportTestSuite : public vigra::test_suite
{
    HDF5ImportExportTestSuite()
        : vigra::test_suite("HDF5ImportExportTestSuite")
    {
		// row major order (C)
        add(testCase(&HDF5ExportImportTest::testUnstridedHDF5ExportImport_rowMajor));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport1_rowMajor));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport2_rowMajor));
        add(testCase(&HDF5ExportImportTest::testAppendNewDataToHDF5_rowMajor));
        add(testCase(&HDF5ExportImportTest::testOverwriteExistingDataInHDF5_rowMajor));
		// column major order (VIGRA)
        add(testCase(&HDF5ExportImportTest::testUnstridedHDF5ExportImport_columnMajor));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport1_columnMajor));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport2_columnMajor));
        add(testCase(&HDF5ExportImportTest::testAppendNewDataToHDF5_columnMajor));
        add(testCase(&HDF5ExportImportTest::testOverwriteExistingDataInHDF5_columnMajor));
	}
};


int main (int argc, char ** argv)
{
    HDF5ImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

	return failed != 0;
}
