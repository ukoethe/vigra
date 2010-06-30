/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Michael Hanselmann and Ullrich Koethe        */
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

	void testScalarUnstridedHDF5ExportImport()
	{
		// export and import data from and to scalar unstrided array

		char hdf5File[] = "testfile_scalar_unstrided.hdf5";

		// data 1: 
		MultiArray<3,int> out_data_1(MultiArrayShape<3>::type(2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_1.data () [i] = i - 12;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// export
		// ...data 1: 
		char hdf5group_1[] = "group/subgroup/data_unstrided";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);

		// import
		// data 1: 
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<3,int> in_data_1(MultiArrayShape<3>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1), infoHDF5_1.shapeOfDimension(2)));
        readHDF5(infoHDF5_1, in_data_1);

		// compare content
		should (in_data_1 == out_data_1);
	}


	void testScalarStridedHDF5ExportImport()
	{
		// export and import data from and to scalar strided array

		char hdf5File[] = "testfile_scalar_strided.hdf5";

		// data 1: unstrided
		MultiArray<3,int> out_data_1(MultiArrayShape<3>::type(2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_1.data () [i] = i - 12;
		// data 2: strided
		MultiArrayView<2,int,StridedArrayTag> out_data_2(out_data_1.bindAt(1, 0));
		
		// export
		// ...data 2: strided
		char hdf5group_2[] = "group/subgroup/data_strided";
		writeHDF5(hdf5File, hdf5group_2, out_data_2);

		// import
		// data 2: strided
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
		// import to strided array
		MultiArray<3,int> in_data_1(MultiArrayShape<3>::type(3, infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1))); // strided
		in_data_1.init(42);
		MultiArrayView<2,int,StridedArrayTag> in_data_2(in_data_1.bindInner(0));
        readHDF5(infoHDF5_2, in_data_2);

		// compare content
		should (in_data_2 == out_data_2);
	}


	void testRGBValueUnstridedHDF5ExportImport()
	{
		// export and import data from and to scalar unstrided array

		char hdf5File[] = "testfile_rgbvalue_unstrided.hdf5";

		// data 1: 3D multi array
		// create RGB image from that
		MultiArray< 2, RGBValue<double> > out_data_1(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_1.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);

		// export
		// ...data 2: RGB image
		char hdf5group_1[] = "group/subgroup/data_rgb";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);

		// import
		// data 1: 
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray< 2, RGBValue<double> > in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(1), infoHDF5_1.shapeOfDimension(2)));
        readHDF5(infoHDF5_1, in_data_1);

		// compare content
		should (in_data_1 == out_data_1);
	}


	void testRGBValueStridedHDF5ExportImport()
	{
		// export and import data from and to scalar unstrided array

		char hdf5File[] = "testfile_rgbvalue_strided.hdf5";

		// data 1: 3D multi array
		// create RGB image from that
		MultiArray< 3, RGBValue<double> > out_data_1(MultiArrayShape<3>::type(5, 3, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 120; ++i)
            out_data_1.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);
		MultiArrayView< 2, RGBValue<double>, StridedArrayTag> out_data_2(out_data_1.bindAt(1, 0));

		// export
		// ...data 2: RGB image
		char hdf5group_2[] = "group/subgroup/data_rgb";
		writeHDF5(hdf5File, hdf5group_2, out_data_2);

		// import
		// data 2: strided
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
		// import to strided array
		MultiArray< 3, RGBValue<double> > in_data_1(MultiArrayShape<3>::type(3, infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2))); // strided
		//in_data_1.init(42);
		MultiArrayView< 2, RGBValue<double>, StridedArrayTag > in_data_2(in_data_1.bindInner(0));

        readHDF5(infoHDF5_2, in_data_2);

		// compare content
		should (in_data_2 == out_data_2);
	}

	void testTinyVectorUnstridedHDF5ExportImport()
	{
		// export and import data from and to scalar unstrided array

		char hdf5File[] = "testfile_tinyvector_unstrided.hdf5";

		// data 1: 3D multi array
		// create TinyVector image from that
		MultiArray< 2, TinyVector<double, 4> > out_data_1(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_1.data () [i] = TinyVector<double, 4>(i + 0.1, i + 0.2, i + 0.3, i + 0.4);

		// export
		// ...data 2: TinyVector image
		char hdf5group_1[] = "group/subgroup/data_tinyvector";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);

		// import
		// data 1: 
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray< 2, TinyVector<double, 4> > in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(1), infoHDF5_1.shapeOfDimension(2)));
        readHDF5(infoHDF5_1, in_data_1);

		// compare content
		should (in_data_1 == out_data_1);
	}


	void testTinyVectorStridedHDF5ExportImport()
	{
		// export and import data from and to scalar unstrided array

		char hdf5File[] = "testfile_tinyvector_strided.hdf5";

		// data 1: 3D multi array
		// create TinyVector image from that
		MultiArray< 3, TinyVector<double, 4> > out_data_1(MultiArrayShape<3>::type(5, 3, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 120; ++i)
            out_data_1.data () [i] = TinyVector<double, 4>(i + 0.1, i + 0.2, i + 0.3, i + 0.4);
		MultiArrayView< 2, TinyVector<double, 4>, StridedArrayTag> out_data_2(out_data_1.bindAt(1, 0));

		// export
		// ...data 2: TinyVector image
		char hdf5group_2[] = "group/subgroup/data_tinyvector";
		writeHDF5(hdf5File, hdf5group_2, out_data_2);

		// import
		// data 2: strided
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
		// import to strided array
		MultiArray< 3, TinyVector<double, 4> > in_data_1(MultiArrayShape<3>::type(3, infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2))); // strided
		MultiArrayView< 2, TinyVector<double, 4>, StridedArrayTag > in_data_2(in_data_1.bindInner(0));

        readHDF5(infoHDF5_2, in_data_2);

		// compare content
		should (in_data_2 == out_data_2);
	}


	void testScalarToRGBValueUnstridedHDF5ExportImport()
	{
		// export scalar 3D array and import as RGB image (unstrided)

		char hdf5File[] = "testfile_scalar2rgbvalue_unstrided.hdf5";

		// data 1: 3D multi array
		MultiArray< 3, double > out_data_1(MultiArrayShape<3>::type(3, 5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 120; ++i)
            out_data_1.data () [i] = i;

		// export
		// ...data 1: 
		char hdf5group_1[] = "group/subgroup/data";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);

		// import
		// data 1: 
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray< 2, RGBValue<double> > in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(1), infoHDF5_1.shapeOfDimension(2)));
        readHDF5(infoHDF5_1, in_data_1);

		// compare content
		int i = 0;
		for (int j=0; j<120; ++i, j+=3)
			should (in_data_1.data () [i] == RGBValue<double>((double)out_data_1.data () [j], (double)out_data_1.data () [j+1], (double)out_data_1.data () [j+2]) );
	}


	void testOverwriteExistingDataInHDF5()
	{
		// write data to file, close file, open file, overwrite data, compare

		char hdf5File[] = "testfile_overwrite_existing_data.hdf5";

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		// ...data 1: 
		char hdf5group_1[] = "group/subgroup/data1";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);
		// ...data 2, overwriting existing data
		writeHDF5(hdf5File, hdf5group_1, out_data_2);

		// import
		// data 1: 
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_1);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
        readHDF5(infoHDF5_2, in_data_2);

		// compare content
		should (in_data_2 == out_data_2);
	}

	void testAppendNewDataToHDF5()
	{
		// write data to file, close file, open file, write more data, compare

		char hdf5File[] = "testfile_append_new_data.hdf5";

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		// ...data 1
		char hdf5group_1[] = "group/data1";
		writeHDF5(hdf5File, hdf5group_1, out_data_1);

		// append to existing file
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data2";
		writeHDF5(hdf5File, hdf5group_2, out_data_2);

		// import
		// ...data 1
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,int> in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        readHDF5(infoHDF5_1, in_data_1);
		// ...data 2
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
        readHDF5(infoHDF5_2, in_data_2);

		// compare content
		// ...data 1
		should (in_data_1 == out_data_1);
		// ...data 2
		should (in_data_2 == out_data_2);
	}

    void testHDF5FileDataAccess()
    {
        //write some data and read it again. Only spot test general functionality.

        std::string file_name( "testfile_HDF5File_data_access.hdf5");

        // data 1: int data in 2 dimensions (partly negative)
        MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;

        // data 2: double data in 4 dimensions (partly negative)
        MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

        // data 3: 2+1D multi array
        MultiArray< 2, TinyVector<double, 4> > out_data_3(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_3.data () [i] = TinyVector<double, 4>(i + 0.1, i + 0.2, i + 0.3, i + 0.4);

        // data 4: RGB values
        MultiArray< 2, RGBValue<double> > out_data_4(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_4.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);

        //create a file
        HDF5File file (file_name, HDF5File::New);
        

        //write one dataset in each group level
        file.write("/dataset",out_data_1);
        file.cd_mk("/group/");
        file.write("/group/dataset",out_data_2);
        file.mkdir("subgroup1");
        file.write("subgroup1/dataset",out_data_3);
        file.cd("..");
        file.write("/dataset_rgb",out_data_4);
        file.writeAtomic("/atomicint", (int)-42);
        file.writeAtomic("/atomicuint", (unsigned int)42);
        file.writeAtomic("/atomicdouble", (double)3.1);


        //create a new dataset
        MultiArrayShape<3>::type shape (50,50,50);
        unsigned char init = 42;
        file.createDataset<3,unsigned char>("/newset", shape, init );


        // check if data is really written

        MultiArray<2,int> in_data_1 (MultiArrayShape<2>::type(10, 11));
        file.read("dataset",in_data_1);

        MultiArray<4,double> in_data_2 (MultiArrayShape<4>::type(10, 2, 3, 4));
        file.read("/group/dataset",in_data_2);

        MultiArray< 2, TinyVector<double, 4> > in_data_3 (MultiArrayShape<2>::type(5,8));
        file.read("/group/subgroup1/dataset",in_data_3);

        MultiArray< 2, RGBValue<double> > in_data_4(MultiArrayShape<2>::type(5,8));
        file.read("/dataset_rgb", in_data_4);

        MultiArray< 3, unsigned char > in_data_5 (shape);
        file.read("/newset",in_data_5);

        int atomicint;
        file.readAtomic("/atomicint",atomicint);
        int atomicuint;
        file.readAtomic("/atomicuint",atomicuint);
        double atomicdouble;
        file.readAtomic("/atomicdouble",atomicdouble);

        file.setAttribute("/newset", "attribute", "This is a string attribute!");
        file.flushToDisk();

        should(file.getAttribute("/newset",  "attribute") == "This is a string attribute!");

        // compare content
        // ...data 1
        should (in_data_1 == out_data_1);

        // ...data 2
        should (in_data_2 == out_data_2);
        // ...data 3
        should (in_data_3 == out_data_3);
        // ...data 4
        should (in_data_4 == out_data_4);
        // ...data 5
        should (in_data_5(1,2,3) == init);

        should (atomicint == -42);
        should (atomicuint == 42);
        should (atomicdouble == 3.1);

        // overwrite existing dataset
        file.write("/dataset",out_data_2);
        file.flushToDisk();

        MultiArray<4,double> in_data_overwrite (MultiArrayShape<4>::type(10, 2, 3, 4));
        file.read("/dataset",in_data_overwrite);

        should(in_data_overwrite == out_data_2);

    }


    void testHDF5FileBlockAccess()
    {
        // Create 3D dataset
        MultiArray< 3, double > out_data(MultiArrayShape<3>::type(10, 10, 10));
        // ...initialize the array to the test data
        for (int i = 0; i < 10; ++i){
            for (int j = 0; j < 10; ++j){
                for (int k = 0; k < 10; ++k){
                    out_data (i,j,k) = 100*i+10*j+k;
                }
            }
        }


        std::string file_name( "testfile_HDF5File_block_access.hdf5");

        //create a file and write data
        HDF5File file (file_name, HDF5File::New);

        file.write("/dataset",out_data);
        file.cd_mk("group");
        file.write("dataset2",out_data);
        file.cd("..");
        file.flushToDisk();


        int sz = 10;
        MultiArray< 3, double > in_data(MultiArrayShape<3>::type(sz, sz, sz));
        MultiArray< 3, double > in_data_group(MultiArrayShape<3>::type(sz, sz, sz));
        MultiArrayShape<3>::type block_offset (0,0,0);
        MultiArrayShape<3>::type block_shape (sz,sz,sz);

        file.readBlock("/dataset", block_offset, block_shape, in_data);
        file.readBlock("/group/dataset2", block_offset, block_shape, in_data_group);

        should(in_data == out_data);
        should(in_data_group == out_data);

        sz = 4;
        MultiArray< 3, double > in_data_2(MultiArrayShape<3>::type(sz, sz, sz));
        MultiArrayShape<3>::type block_offset_2 (2,2,2);
        MultiArrayShape<3>::type block_shape_2 (sz,sz,sz);

        file.readBlock("/dataset", block_offset_2, block_shape_2, in_data_2);

        shouldEqual(in_data_2(0,0,0), out_data(2,2,2));
        shouldEqual(in_data_2(1,0,0), out_data(3,2,2));
        shouldEqual(in_data_2(0,2,0), out_data(2,4,2));
        shouldEqual(in_data_2(0,0,3), out_data(2,2,5));

        // write the data to different position
        MultiArrayShape<3>::type block_offset_3 (6,2,2);
        file.writeBlock("/dataset", block_offset_3, in_data_2 );
        file.flushToDisk();

        // now read it again and compare
        MultiArray< 3, double > in_data_3(MultiArrayShape<3>::type(sz, sz, sz));
        file.readBlock("/dataset", block_offset_3, block_shape_2, in_data_3);

        should(in_data_2 == in_data_3);

    }


    void testHDF5FileChunks()
    {
        //write some data and read it again. Only spot test general functionality.

        std::string file_name( "testfile_HDF5File_chunks.hdf5");

        // data 1: int data in 2 dimensions (partly negative)
        MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;

        // data 2: double data in 4 dimensions (partly negative)
        MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

        // data 3: 2+1D multi array
        MultiArray< 2, TinyVector<double, 4> > out_data_3(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_3.data () [i] = TinyVector<double, 4>(i + 0.1, i + 0.2, i + 0.3, i + 0.4);

        // data 4: RGB values
        MultiArray< 2, RGBValue<double> > out_data_4(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_4.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);

        //create a file
        HDF5File file (file_name, HDF5File::New);

        //write one dataset in each group level
        file.write("dataset",out_data_1,5);
        file.write("group/dataset",out_data_2, MultiArrayShape<4>::type(5,2,3,2));
        file.write("group/subgroup1/dataset",out_data_3, MultiArrayShape<2>::type(5,4));
        file.write("dataset_rgb",out_data_4, MultiArrayShape<2>::type(5,4));


        //create a new dataset
        MultiArrayShape<3>::type shape (50,50,50);
        MultiArrayShape<3>::type chunks (10,10,10);
        unsigned char init = 42;
        file.createDataset<3,unsigned char>("newset", shape, init, chunks);


        // check if data is really written

        MultiArray<2,int> in_data_1 (MultiArrayShape<2>::type(10, 11));
        file.read("dataset",in_data_1);

        MultiArray<4,double> in_data_2 (MultiArrayShape<4>::type(10, 2, 3, 4));
        file.read("/group/dataset",in_data_2);

        MultiArray< 2, TinyVector<double, 4> > in_data_3 (MultiArrayShape<2>::type(5, 8));
        file.read("/group/subgroup1/dataset",in_data_3);

        MultiArray< 2, RGBValue<double> > in_data_4(MultiArrayShape<2>::type(5,8));
        file.read("/dataset_rgb", in_data_4);

        MultiArray< 3, unsigned char > in_data_5 (shape);
        file.read("/newset",in_data_5);

        // compare content
        // ...data 1
        should (in_data_1 == out_data_1);
        // ...data 2
        should (in_data_2 == out_data_2);
        // ...data 3
        should (in_data_3 == out_data_3);
        // ...data 4
        should (in_data_4 == out_data_4);
        // ...data 5
        should (in_data_5(1,2,3) == init);

    }


    void testHDF5FileCompression()
    {
        //write some data and read it again. Only spot test general functionality.

        std::string file_name( "testfile_HDF5File_compression.hdf5");

        // data 1: int data in 2 dimensions (partly negative)
        MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;

        // data 2: double data in 4 dimensions (partly negative)
        MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

        // data 3: 2+1D multi array
        MultiArray< 2, TinyVector<double, 4> > out_data_3(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_3.data () [i] = TinyVector<double, 4>(i + 0.1, i + 0.2, i + 0.3, i + 0.4);

        // data 4: RGB values
        MultiArray< 2, RGBValue<double> > out_data_4(MultiArrayShape<2>::type(5, 8));
        // ...initialize the array to the test data
        for (int i = 0; i < 40; ++i)
            out_data_4.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);

        //create a file
        HDF5File file (file_name, HDF5File::New);

        //write one dataset in each group level
        file.write("/dataset",out_data_1, 5, 3);
        file.write("/group/dataset",out_data_2, MultiArrayShape<4>::type(5,2,3,2), 5);
        file.write("/group/subgroup1/dataset",out_data_3, MultiArrayShape<2>::type(5,4), 7);
        file.write("/dataset_rgb",out_data_4, MultiArrayShape<2>::type(5,4), 9);


        //create a new dataset
        MultiArrayShape<3>::type shape (50,50,50);
        MultiArrayShape<3>::type chunks (10,10,10);
        unsigned char init = 42;
        file.createDataset<3,unsigned char>("/newset", shape, init, chunks, 5);

        // read block
        int sz = 10;
        MultiArray< 3, unsigned char > in_data_block(MultiArrayShape<3>::type(sz, sz, sz));
        MultiArrayShape<3>::type block_offset (5,5,5);
        MultiArrayShape<3>::type block_shape (sz,sz,sz);

        file.readBlock("/newset", block_offset, block_shape, in_data_block);


        // check if data is really written

        MultiArray<2,int> in_data_1 (MultiArrayShape<2>::type(10, 11));
        file.read("dataset",in_data_1);

        MultiArray<4,double> in_data_2 (MultiArrayShape<4>::type(10, 2, 3, 4));
        file.read("/group/dataset",in_data_2);

        MultiArray< 2, TinyVector<double, 4> > in_data_3 (MultiArrayShape<2>::type(5, 8));
        file.read("/group/subgroup1/dataset",in_data_3);

        MultiArray< 2, RGBValue<double> > in_data_4(MultiArrayShape<2>::type(5,8));
        file.read("/dataset_rgb", in_data_4);

        MultiArray< 3, unsigned char > in_data_5 (shape);
        file.read("/newset",in_data_5);

        // compare content
        // ...data 1
        should (in_data_1 == out_data_1);

        // ...data 2
        should (in_data_2 == out_data_2);
        // ...data 3
        should (in_data_3 == out_data_3);
        // ...data 4
        should (in_data_4 == out_data_4);
        // ...data 5
        should (in_data_5(1,2,3) == init);
        // ...data block
        should (in_data_block(1,2,3) == 42);
    }


    void testHDF5FileBrowsing()
    {
        //create groups, change current group, ...

        std::string file_name( "testfile_HDF5File_browsing.hdf5");

        //create a file
        HDF5File file (file_name, HDF5File::New);

        //we should be in root group in the beginning
        should(file.pwd() == "/" );

        //create group "group1"
        file.mkdir("group1");
        //we should still be in root (only created group)
        should(file.pwd() == "/" );

        //now change to new group (relative change)
        file.cd("group1");
        should(file.pwd() == "/group1" );

        //create a subgroup and change there
        file.cd_mk("subgroup");
        should(file.pwd() == "/group1/subgroup" );

        //create a group tree in root group
        file.mkdir("/group2/subgroup/subsubgroup");
        should(file.pwd() == "/group1/subgroup" );

        //go to parent group
        file.cd_up();
        should(file.pwd() == "/group1" );

        //change to recently created group tree (absolute change)
        file.cd("/group2/subgroup/subsubgroup");
        should(file.pwd() == "/group2/subgroup/subsubgroup" );

        //change to parent group
        file.cd("..");
        should(file.pwd() == "/group2/subgroup" );

        //change up 2 groups
        file.cd_up(2);
        should(file.pwd() == "/" );

        //try to change to parent of root group
        should( !file.cd_up() );

        std::vector<std::string> entries;
        entries = file.ls();
        should(entries.size() == 2);
        should(entries[0] == "group1/");
        should(entries[1] == "group2/");
    }




};



struct HDF5ImportExportTestSuite : public vigra::test_suite
{
    HDF5ImportExportTestSuite()
        : vigra::test_suite("HDF5ImportExportTestSuite")
    {
		// tests for scalar data
		add(testCase(&HDF5ExportImportTest::testScalarUnstridedHDF5ExportImport));
		add(testCase(&HDF5ExportImportTest::testScalarStridedHDF5ExportImport));

		// tests for non-scalar data
		// RGBValue
		add(testCase(&HDF5ExportImportTest::testRGBValueUnstridedHDF5ExportImport));
		add(testCase(&HDF5ExportImportTest::testRGBValueStridedHDF5ExportImport));
		// TinyVector
		add(testCase(&HDF5ExportImportTest::testTinyVectorUnstridedHDF5ExportImport));
		add(testCase(&HDF5ExportImportTest::testTinyVectorStridedHDF5ExportImport));

		// mixed forms
		add(testCase(&HDF5ExportImportTest::testScalarToRGBValueUnstridedHDF5ExportImport));

		// general tests
		add(testCase(&HDF5ExportImportTest::testOverwriteExistingDataInHDF5));
		add(testCase(&HDF5ExportImportTest::testAppendNewDataToHDF5));

        // HDF5File tests
        add(testCase(&HDF5ExportImportTest::testHDF5FileDataAccess));
        add(testCase(&HDF5ExportImportTest::testHDF5FileBlockAccess));
        add(testCase(&HDF5ExportImportTest::testHDF5FileChunks));
        add(testCase(&HDF5ExportImportTest::testHDF5FileCompression));
        add(testCase(&HDF5ExportImportTest::testHDF5FileBrowsing));
 
	}
};


int main (int argc, char ** argv)
{
    HDF5ImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

	return failed != 0;
}
