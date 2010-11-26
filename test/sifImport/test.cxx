/************************************************************************/
/*                                                                      */
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
		SIFImportInfo infoSIF(sifFile);
		MultiArray<3,float> in_data(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
        readSIF(infoSIF, in_data);

		// compare content
		should (in_data == reference_data);
	}

	// SIF from Andor software version 4.13.30000
	void testSifImport_4_13() {
		char sifFile[] = "testSif_4_13_30000.sif";

		const float reference_array[] = 
			{0, 27, 61, 96, 132, 166, 200, 235, 
			 0, 27, 62, 97, 131, 166, 200, 235};
		MultiArray<3, float> reference_data(MultiArrayShape<3>::type(8, 2, 1), reference_array);
		

		// read SIF
		SIFImportInfo infoSIF(sifFile);
		// import to MultiArrayView
		MultiArray<3,float> in_data(MultiArrayShape<3>::type(infoSIF.width(), infoSIF.height(), infoSIF.stacksize()));
        readSIF(infoSIF, in_data);
		
		
		// compare
		should (in_data == reference_data);
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
 
	}
};


int main (int argc, char ** argv)
{
    SifImportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

	return failed != 0;
}
