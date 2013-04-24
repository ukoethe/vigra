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
        file.write("/atomicint", (int)-42);
        file.write("/atomicuint", (unsigned int)42);
        file.write("/atomicdouble", (double)3.1);

        //create a new dataset
        MultiArrayShape<3>::type shape (50,50,50);
        unsigned char init = 42;
        file.createDataset<3,unsigned char>("/newset", shape, init );

        //test compressed
        MultiArray<3, double> out_data_6(Shape3(100,100,30), 42.0);
        file.write("/compressed", out_data_6, 0, 6);

        file.close();
        file.open(file_name, HDF5File::Open);
        
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

        MultiArray< 3, double > in_data_6 (out_data_6.shape());
        file.read("/compressed",in_data_6);

        int atomicint;
        file.read("/atomicint",atomicint);
        int atomicuint;
        file.read("/atomicuint",atomicuint);
        double atomicdouble;
        file.read("/atomicdouble",atomicdouble);

        file.flushToDisk();

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
        // ...data 6
        should (in_data_6 == out_data_6);

        should (atomicint == -42);
        should (atomicuint == 42);
        should (atomicdouble == 3.1);

        // readAndResize
        MultiArray<2,int> in_re_data_1;
        file.readAndResize("dataset",in_re_data_1);

        MultiArray<4,double> in_re_data_2;
        file.readAndResize("/group/dataset",in_re_data_2);

        MultiArray< 2, TinyVector<double, 4> > in_re_data_3;
        file.readAndResize("/group/subgroup1/dataset",in_re_data_3);

        MultiArray< 2, RGBValue<double> > in_re_data_4;
        file.readAndResize("/dataset_rgb", in_re_data_4);

        MultiArray< 3, unsigned char > in_re_data_5 ;
        file.readAndResize("/newset",in_re_data_5);

        // compare content
        // ...data 1
        should (in_re_data_1 == out_data_1);
        // ...data 2
        should (in_re_data_2 == out_data_2);
        // ...data 3
        should (in_re_data_3 == out_data_3);
        // ...data 4
        should (in_re_data_4 == out_data_4);
        // ...data 5
        should (in_re_data_5(1,2,3) == init);

        // overwrite existing dataset
        file.write("/dataset",out_data_2);
        file.flushToDisk();

        MultiArray<4,double> in_data_overwrite (MultiArrayShape<4>::type(10, 2, 3, 4));
        file.read("/dataset",in_data_overwrite);

        should(in_data_overwrite == out_data_2);
        
        file.cd("/");
        
        // writing and reading single values

        file.write("set_char",'A');
        file.write("set_int8",(Int8)8);
        file.write("set_int16",(Int16)16);
        file.write("set_int32",(Int32)32);
        file.write("set_int64",(Int64)64);
        file.write("set_uint8",(UInt8)8);
        file.write("set_uint16",(UInt16)16);
        file.write("set_uint32",(UInt32)32);
        file.write("set_uint64",(UInt64)64);
        file.write("set_float",(float)1.);
        file.write("set_double",(double)2.);
        file.write("set_longdouble",(long double)3.);
        file.write("set_string",std::string("abc").c_str());
        file.write("set_string2",std::string("abcdef"));
        
        char read_char = 0;
        file.read("set_char",read_char);
        should(read_char=='A');
        
        Int8 read_int8 = 0;
        file.read("set_int8",read_int8);
        should(read_int8==8);
        
        Int16 read_int16 = 0;
        file.read("set_int16",read_int16);
        should(read_int16 == 16);
        
        Int32 read_int32 = 0;
        file.read("set_int32",read_int32);
        should(read_int32 == 32);
        
        Int64 read_int64 = 0;
        file.read("set_int64",read_int64);
        should(read_int64 == 64);
        
        UInt8 read_uint8 = 0;
        file.read("set_uint8",read_uint8);
        should(read_uint8 == 8);
        
        UInt16 read_uint16 = 0;
        file.read("set_uint16",read_uint16);
        should(read_uint16 == 16);
        
        UInt32 read_uint32 = 0;
        file.read("set_uint32",read_uint32);
        should(read_uint32 == 32);
        
        UInt64 read_uint64 = 0;
        file.read("set_uint64",read_uint64);
        should(read_uint64 == 64);
        
        float read_float = 0;
        file.read("set_float",read_float);
        should(read_float == 1.);
        
        double read_double = 0;
        file.read("set_double",read_double);
        should(read_double == 2.);
        
        long double read_longdouble = 0;
        file.read("set_longdouble",read_longdouble);
        should(read_longdouble == 3.);
        
        std::string read_string = "";
        file.read("set_string",read_string);  
        should(read_string == "abc");

        std::string read_string2 = "";
        file.read("set_string2",read_string2);
        should(read_string2 == "abcdef");
       
        // test read-only opening
        file.close();
        file.open(file_name, HDF5File::OpenReadOnly);
        read_string2 = "";
        file.read("set_string2",read_string2);
        should(read_string2 == "abcdef");
        file.close();
    }

    // reading and writing attributes. get handles of groups, datasets and attributes
    void testHDF5FileAttributes()
    {
        std::string file_name( "testfile_HDF5File_data_attributes.hdf5");

        // double data in 4 dimensions (partly negative)
        MultiArray<4,double> out_data_1(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_1.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

        // variable length string MultiArrays
        MultiArray<2,const char* > out_data_2(MultiArrayShape<2>::type(2, 2));
        // ...initialize the array to the test data
        out_data_2(0,0) = "Brot";
        out_data_2(0,1) = "Huette";
        out_data_2(1,0) = "Kamuffel";
        out_data_2(1,1) = "Gorch Frog";

        //create a file
        HDF5File file (file_name, HDF5File::New);

        //write one dataset in each group level
        file.write("/double/dataset",out_data_1);
        file.write("/string/dataset",out_data_2);

        // check if data is really written
        MultiArray< 2, const char* > in_data_2(MultiArrayShape<2>::type(2, 2));
        file.read("/string/dataset",in_data_2);

        should(std::string(in_data_2(0,0)) == std::string("Brot"));
        should(std::string(in_data_2(0,1)) == std::string("Huette"));
        should(std::string(in_data_2(1,0)) == std::string("Kamuffel"));
        should(std::string(in_data_2(1,1)) == std::string("Gorch Frog"));

        // write Attributes

        // integer attribute
        MultiArray<2,int> out_attr_1(MultiArrayShape<2>::type(2,3));
        // ...initialize the array to the test data
        for (int i = 0; i < 6; ++i)
            out_attr_1.data() [i] = i;

        // variable length string attribute
        MultiArray<1,const char*> out_attr_2(MultiArrayShape<1>::type(4));
        // ...initialize the array to the test data
        out_attr_2(3) = "Brot";
        out_attr_2(2) = "Huette";
        out_attr_2(1) = "Kamuffel";
        out_attr_2(0) = "Gorch Frog";

        // tiny vector multi array
        MultiArray<2, TinyVector<double, 3> > out_attr_3(MultiArrayShape<2>::type(2,1));
        // ...initialize the array to the test data
        for (int i = 0; i < 2; ++i)
            out_attr_3.data () [i] = TinyVector<double, 3>(i + 0.1, i + 0.2, i + 0.3);

        // data 4: RGB values
        MultiArray< 2, RGBValue<double> > out_attr_4(MultiArrayShape<2>::type(2,3));
        // ...initialize the array to the test data
        for (int i = 0; i < 6; ++i)
            out_attr_4.data () [i] = RGBValue<double>(i + 0.1, i + 0.2, i + 0.3);

        file.writeAttribute("/double/dataset","int attribute", out_attr_1);
        file.writeAttribute("/double/dataset","string attribute", out_attr_2);
        file.writeAttribute("/string/dataset","tinyvector attribute", out_attr_3);
        file.writeAttribute("/string/dataset","rgb attribute", out_attr_4);


        // read attributes
        MultiArray<2,int> in_attr_1(MultiArrayShape<2>::type(2,3));
        MultiArray<1,const char*> in_attr_2(MultiArrayShape<1>::type(4));
        MultiArray<2, TinyVector<double, 3> > in_attr_3(MultiArrayShape<2>::type(2,1));
        MultiArray< 2, RGBValue<double> > in_attr_4(MultiArrayShape<2>::type(2,3));

        file.readAttribute("/double/dataset","int attribute", in_attr_1);
        file.readAttribute("/double/dataset","string attribute", in_attr_2);
        file.readAttribute("/string/dataset","tinyvector attribute", in_attr_3);
        file.readAttribute("/string/dataset","rgb attribute", in_attr_4);

        should(in_attr_1 == out_attr_1);
        should(std::string(in_attr_2(0)) == std::string(out_attr_2(0)));
        should(std::string(in_attr_2(1)) == std::string(out_attr_2(1)));
        should(std::string(in_attr_2(2)) == std::string(out_attr_2(2)));
        should(std::string(in_attr_2(3)) == std::string(out_attr_2(3)));
        should(in_attr_3 == out_attr_3);
        should(in_attr_4 == out_attr_4);



        // write and read attributes

        file.cd("/");
        file.write("attrset",std::string("Dataset with many attributes").c_str());

        file.writeAttribute("attrset","set_char",'A');
        file.writeAttribute("attrset","set_int8",(Int8)8);
        file.writeAttribute("attrset","set_int16",(Int16)16);
        file.writeAttribute("attrset","set_int32",(Int32)32);
        file.writeAttribute("attrset","set_int64",(Int64)64);
        file.writeAttribute("attrset","set_uint8",(UInt8)8);
        file.writeAttribute("attrset","set_uint16",(UInt16)16);
        file.writeAttribute("attrset","set_uint32",(UInt32)32);
        file.writeAttribute("attrset","set_uint64",(UInt64)64);
        file.writeAttribute("attrset","set_float",(float)1.);
        file.writeAttribute("attrset","set_double",(double)2.);
        file.writeAttribute("attrset","set_longdouble",(long double)3.);
        file.writeAttribute("attrset","set_string",std::string("abc").c_str());
        file.writeAttribute("attrset","set_string2",std::string("abcdef"));

        char read_char = 0;
        file.readAttribute("attrset","set_char",read_char);
        should(read_char=='A');

        Int8 read_int8 = 0;
        file.readAttribute("attrset","set_int8",read_int8);
        should(read_int8==8);

        Int16 read_int16 = 0;
        file.readAttribute("attrset","set_int16",read_int16);
        should(read_int16 == 16);

        Int32 read_int32 = 0;
        file.readAttribute("attrset","set_int32",read_int32);
        should(read_int32 == 32);

        Int64 read_int64 = 0;
        file.readAttribute("attrset","set_int64",read_int64);
        should(read_int64 == 64);

        UInt8 read_uint8 = 0;
        file.readAttribute("attrset","set_uint8",read_uint8);
        should(read_uint8 == 8);

        UInt16 read_uint16 = 0;
        file.readAttribute("attrset","set_uint16",read_uint16);
        should(read_uint16 == 16);

        UInt32 read_uint32 = 0;
        file.readAttribute("attrset","set_uint32",read_uint32);
        should(read_uint32 == 32);

        UInt64 read_uint64 = 0;
        file.readAttribute("attrset","set_uint64",read_uint64);
        should(read_uint64 == 64);

        float read_float = 0;
        file.readAttribute("attrset","set_float",read_float);
        should(read_float == 1.);

        double read_double = 0;
        file.readAttribute("attrset","set_double",read_double);
        should(read_double == 2.);

        long double read_longdouble = 0;
        file.readAttribute("attrset","set_longdouble",read_longdouble);
        should(read_longdouble == 3.);

        std::string read_string = "";
        file.readAttribute("attrset","set_string",read_string);
        should(read_string == "abc");

        std::string read_string2 = "";
        file.readAttribute("attrset","set_string2",read_string2);
        should(read_string2 == "abcdef");

        // get handles
        HDF5Handle group_handle = file.getGroupHandle("/string");
        HDF5Handle dataset_handle = file.getDatasetHandle("/string/dataset");
        HDF5Handle attribute_handle = file.getAttributeHandle("/string/dataset","rgb attribute");

        should(group_handle.get() > 0);
        should(dataset_handle.get() > 0);
        should(attribute_handle.get() > 0);

        should(group_handle.close() >= 0);
        should(dataset_handle.close() >= 0);
        should(attribute_handle.close() >= 0);
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

        //enhanced navigation with .. and .
        file.cd("/group1/");
        should(file.pwd() == "/group1");
        file.cd("../not/existing/../../group2/././subgroup/");
        should(file.pwd() == "/group2/subgroup");
    }




    void testHDF5FileTutorial()
    {
        // First create a new HDF5 file (thereby also testing default construction + open);
        HDF5File file;
        file.open("tutorial_HDF5File.h5", HDF5File::New);

        // we should be in root group in the beginning
        should(file.pwd() == "/" );

        // create the group "/group1/subgroup1" using absolute paths
        file.mkdir("/group1/subgroup1");

        // create the group "group2/subgroup2" using relative paths
        file.mkdir("group2");
        file.cd("group2");
        file.cd_mk("subgroup2"); //cd_mk first creates group and then opens it

        should(file.pwd() == "/group2/subgroup2" );


        // Writing Data

        // Create a new dataset with shape (10,10,10) of int, initialized with 42.
        MultiArrayShape<3>::type shape (10,10,10);
        file.createDataset<3,int>("new_dataset",shape,42);

        // Create a new large dataset (20,20,20) with anisotropic (5,10,10) chunks and compression level 5.
        shape = MultiArrayShape<3>::type(20,20,20);
        MultiArrayShape<3>::type chunks (5,10,10);
        file.createDataset<3,int>("new_dataset_chunks",shape, 42, chunks, 5);

        // Create a double MultiArray
        MultiArray<4,double> out_data(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array with some data
        for (int i = 0; i < 240; ++i)
            out_data.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

        // Write double MultiArray to "/group2". Use relative paths and ".." notation.
        file.write("../double_array",out_data);

        // Create a large MultiArray
        MultiArray<3,float> out_data_float (shape, 42.);

        // write with isotropic chunks (10,10,10) and compression 5.
        file.write("../float_array",out_data_float,10,5);

        // Write a (10,10,10) Block of data into the center of "/group2/float_array"
        MultiArray<3,float> float_block (MultiArrayShape<3>::type(10,10,10), 42.);
        MultiArrayShape<3>::type block_offset (5,5,5);

        file.writeBlock("/group2/float_array",block_offset,float_block);

        // Write a single long integer value in the root group
        file.write("/single_value",(long int)23);

        // Attach simple string attribute to dataset in root group
        file.writeAttribute("/single_value","some_attribute_name","This is a simple string attribute");

        // Attach a MultiArray attribute
        MultiArray<1,double> attr_data (MultiArrayShape<1>::type(3));
        attr_data(0) = 10; attr_data(1) = 100; attr_data(2) = 1000;

        file.writeAttribute("/group2/float_array","float_array_attribute",attr_data);

        file.flushToDisk();


        // Reading Data

        // it is possible to access one file with two different HDF5File instances.
        HDF5File file_open ("tutorial_HDF5File.h5", HDF5File::Open);

        file_open.cd("/");

        // read "/group2/subgroup2/new_dataset_chunks", prepare array with correct shape
        MultiArrayShape<3>::type read_shape (20,20,20);
        MultiArray<3,int> read_chunks (read_shape);

        file_open.read("/group2/subgroup2/new_dataset_chunks",read_chunks);

        // read "/group2/subgroup2/new_dataset", reshape MultiArray automatically
        MultiArray<3,int> read_reshape;

        file_open.readAndResize("/group2/subgroup2/new_dataset",read_reshape);

        // read a block of data from "/group2/float_array"
        block_offset = MultiArrayShape<3>::type(0,0,0);
        MultiArrayShape<3>::type block_shape (15,15,15);
        MultiArray<3,float> read_block (block_shape);

        file_open.readBlock("/group2/float_array",block_offset,block_shape,read_block);

        // read the long int value from "/single_value"
        long int single_value;
        file_open.read("/single_value",single_value);

        // read the attribute of "/single_value"
        std::string attribute_string;
        file_open.readAttribute("/single_value","some_attribute_name", attribute_string);

        // read the MultiArray Attribute of "/group2/float_array"
        MultiArray<1,double> read_attr (MultiArrayShape<1>::type(3));
        file_open.readAttribute("/group2/float_array","float_array_attribute",read_attr);
    }

    struct HDF5File_close_test : public HDF5File
    {
        HDF5File_close_test(const std::string & name, HDF5File::OpenMode mode = HDF5File::New)
            : HDF5File(name, mode) {}

        hid_t get_file_id() const
        {
            return fileHandle_;
        }

        void closeCurrentGroup()
        {
            cGroupHandle_.close();
        }

            // this function is for debugging the test
            // (the leak test should be strong enough to detect the leak) 
        hid_t forceGroupLeak(std::string name)
        {
            return openCreateGroup_(name);
        }
    };
    
    void test_file_closing()
    {
        hid_t file_id = 0;
        hid_t leak_id = 0;
        {   // open a new block on purpose.
            HDF5File_close_test test_file("open_file_test.hdf5");
            file_id = test_file.get_file_id();

            // mess around with the file in order to maybe trigger
            // leaking hdf5 object descriptors that would block
            // closing the file
            test_file.cd_mk("subgroup_a"); // this at least used to leak.
            test_file.cd("/");
            test_file.cd_mk("subgroup_c");
            test_file.mkdir("group1");
            test_file.cd_mk("subgroup_b");
            test_file.mkdir("group1");
            test_file.cd("/");
            test_file.mkdir("group1");
            test_file.cd("/");
            test_file.cd("group1");
            test_file.mkdir("/group2/subgroup/subsubgroup");
            test_file.cd_up();
            test_file.cd("/group2/subgroup/subsubgroup");
            test_file.cd("..");
            test_file.cd_up(2);
            test_file.cd_up();
            test_file.ls();
            test_file.cd("/group1/");
            test_file.cd("../not/existing/../../group2/././subgroup/");
            MultiArrayShape<3>::type shape (10, 10, 10);
            test_file.createDataset<3, float>("new_dataset", shape, 1.23f);
            test_file.closeCurrentGroup();
            
            // the file handle must be the only remaining open object
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 1);

            hid_t new_file_id = H5Freopen(file_id);
            should(new_file_id >= 0); // the file must still be open
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 2);

            H5Fclose(new_file_id);
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 1);
        } // this calls ~HDF5File() and therefore absolutely must close the file
        {
            HDF5File_close_test test_file("open_file_test.hdf5", HDF5File::Open);
            file_id = test_file.get_file_id();
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 2);
        }
        {
            HDF5File_close_test test_file("open_file_test.hdf5");
            file_id = test_file.get_file_id();

            // create an intentional leak
            leak_id = test_file.forceGroupLeak("/group1/");

            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 3);
        }
        {
            HDF5File_close_test test_file("open_file_test.hdf5", HDF5File::Open);
            file_id = test_file.get_file_id();
            // check that the file has one more open object (the leak) than it should
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 3);
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL | H5F_OBJ_LOCAL), 2);
            // check that the 'leaked' object is still usable
            H5Gclose(leak_id);
            shouldEqual(H5Fget_obj_count(file_id, H5F_OBJ_ALL), 2);
        }
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
        add(testCase(&HDF5ExportImportTest::testHDF5FileAttributes));
        add(testCase(&HDF5ExportImportTest::testHDF5FileTutorial));
        add(testCase(&HDF5ExportImportTest::test_file_closing));
    }
};


int main (int argc, char ** argv)
{
    HDF5ImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
