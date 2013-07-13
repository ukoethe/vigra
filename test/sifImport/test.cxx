/************************************************************************/
/*                                                                    */
/*       Copyright 2010 by Joachim Schleicher and Ullrich Koethe        */
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

#include "unittest.hxx"
#include "vigra/sifImport.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_impex.hxx"

#ifdef HasHDF5
# include "vigra/hdf5impex.hxx"
#endif

#include "testSif_ref_4_16_30001.hxx"

using namespace vigra;


class SifImportTest
{

public:

    SifImportTest()
    {}

    // Test Sif from Andor software version 4.16.30001
    void testSifImport_4_16() {
        MultiArray<3,float> reference_data(MultiArrayShape<3>::type(128, 128, 1), ref_4_16);
        
        // read SIF data to MultiArrayView
        char sifFile[] = "testSif_4_16_30001.sif";
        {
            SIFImportInfo infoSIF(sifFile);
            MultiArray<3,float> in_data(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
            readSIF(infoSIF, in_data);

            // compare content
            should (infoSIF.shape()[0] == 128);
            should (infoSIF.shape()[1] == 128);
            should (infoSIF.shape()[2] == 1);
            should (in_data == reference_data);
        }
        {
            // read via importVolume()
            VolumeImportInfo info(sifFile);

            // compare content
            should (info.getFileType() == std::string("SIF"));
            should (info.shape() == Shape3(128, 128, 1));
            should (info.getPixelType() == std::string("FLOAT"));
            should (info.numBands() == 1);

            MultiArray<3,float> in_data(info.shape());
            importVolume(info, in_data);
            should (in_data == reference_data);
        }
    }

    // SIF from Andor software version 4.13.30000
    void testSifImport_4_13() {
        char sifFile[] = "testSif_4_13_30000.sif";

        const float reference_array[] = 
            {0.f, 27.f, 61.f, 96.f, 132.f, 166.f, 200.f, 235.f, 
             0.f, 27.f, 62.f, 97.f, 131.f, 166.f, 200.f, 235.f};
        MultiArray<3, float> reference_data(MultiArrayShape<3>::type(8, 2, 1), reference_array);
        

        // read SIF
        SIFImportInfo infoSIF(sifFile);
        // import to MultiArrayView
        MultiArray<3,float> in_data(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
        readSIF(infoSIF, in_data);
        
        
        // compare
        should (infoSIF.shape()[0] == 8);
        should (infoSIF.shape()[1] == 2);
        should (infoSIF.shape()[2] == 1);
        should (in_data == reference_data);
    }
    
    // SIF from Andor software version 4.6.3
    void testSifImport_4_6() {
        char sifFile[] = "testSif_4_6_30000.sif";

        const float reference_array[] = 
            {124.f, 176.114f, 228.f,
             68.f, 120.f, 172.f,
             12.299f, 64.f, 116.f,
             0.f, 8.299f, 60.f};
        MultiArray<3, float> reference_data(MultiArrayShape<3>::type(3, 4, 1), reference_array);
        

        // read SIF
        SIFImportInfo infoSIF(sifFile);
        // import to MultiArrayView
        MultiArray<3,float> in_data(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
        readSIF(infoSIF, in_data);
        
        // compare
        should (infoSIF.shape()[0] == 3);
        should (infoSIF.shape()[1] == 4);
        should (infoSIF.shape()[2] == 1);
        should (in_data == reference_data);
    }

    // check for consistency of shape() vs shapeOfDimension(i)
    void testShapeOfDimension() {
        char sifFile[] = "testSif_4_6_30000.sif";

        SIFImportInfo infoSIF(sifFile);

        // compare
        should (infoSIF.shape()[0] == infoSIF.width());
        should (infoSIF.shape()[1] == infoSIF.height());
        should (infoSIF.shape()[2] == infoSIF.stacksize());
        for (int i = 0; i < 3; ++i) {
            should (infoSIF.shape()[i] == infoSIF.shapeOfDimension(i));
        }
    }

    void testReadBlock() {
        char sifFile[] = "testSif_forBlocks_4_16_30001.sif";

        SIFImportInfo infoSIF(sifFile);
        should (infoSIF.stacksize() == 3);
        should (infoSIF.width() == 4);
        should (infoSIF.height() == 5);
        MultiArray<3,float> in_data_block(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), 1));

        // import whole volume to MultiArray
        MultiArray<3,float> in_data_volume(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
        readSIF(infoSIF, in_data_volume);

        // compare with readBlock() function
        for (int i=0; i<infoSIF.stacksize(); ++i) {
            readSIFBlock(infoSIF, Shape3(0,0,i), Shape3(4,5,1), in_data_block); // read one frame at a time
            should (in_data_block == in_data_volume.subarray(Shape3(0,0,i), Shape3(4,5,i+1)));

            // additional check: compare values explicitly
            for (int xx=0; xx<infoSIF.width(); ++xx) {
                for(int yy=0; yy<infoSIF.height(); ++yy) {
                    should ( in_data_block(xx,infoSIF.height()-1-yy) == (float)(xx | (yy<<4) | ((i+1)<<8)) );
                }
            }
        }
    }

};



struct SifImportTestSuite : public vigra::test_suite
{
    SifImportTestSuite()
        : vigra::test_suite("SifImportTestSuite")
    {
        // tests for import
        add(testCase(&SifImportTest::testSifImport_4_16));
        add(testCase(&SifImportTest::testSifImport_4_13));
        add(testCase(&SifImportTest::testSifImport_4_6));
        add(testCase(&SifImportTest::testShapeOfDimension));
        add(testCase(&SifImportTest::testReadBlock));
 
    }
};


int main (int argc, char ** argv)
{
    SifImportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
