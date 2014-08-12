/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
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

#include <functional>
#include "vigra/unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/error.hxx"
#include "vigra/openmp_vigra.h"

using namespace std;
using namespace vigra;

struct MultiArrayPointoperatorsTest {
	typedef float PixelType;
	typedef MultiArray<5, PixelType> Image5D;
	typedef MultiArrayView<5, PixelType> View5D;
	typedef Image5D::difference_type Size5;
	Image5D img;
	MultiArrayPointoperatorsTest() :
			img(Size5( 60, 60, 60, 60, 60)) {
		int i;
		PixelType c = 0.1f;
		for (i = 0; i < img.elementCount(); ++i, ++c)
			img.data()[i] = c;
	}
	void testCopy() {
		Image5D res(img.shape(), 1.0);
		double start_time = omp_get_wtime();
		vigra::copyMultiArray(srcMultiArrayRange(img), destMultiArray(res));
		double end_time = omp_get_wtime();
		std::cout << "Time " << (end_time-start_time)*1000 << "ms" << std::endl;
		//vigra::copyMultiArray(img, res1);

//        int x,y,z,k;
//        std::cout << "res1=" << std::endl;
//        for(k=0; k<res1.shape(3); ++k)
//        	for(z=0; z<res1.shape(2); ++z)
//        		for(y=0; y<res1.shape(1); ++y)
//        			for(x=0; x<res1.shape(0); ++x)
//        				std::cout << res1(x,y,z,k) << " ";
//        std::cout << std::endl;
//        std::cout << "img=" << std::endl;
//        for(k=0; k<img.shape(3); ++k)
//        	for(z=0; z<img.shape(2); ++z)
//                for(y=0; y<img.shape(1); ++y)
//                    for(x=0; x<img.shape(0); ++x)
//                        std::cout << img(x,y,z,k) << " ";

		should(img == res);
		//should(img == res1);
	}
};
struct MultiArrayPointOperatorsTestSuite: public vigra::test_suite {
	MultiArrayPointOperatorsTestSuite() :
			vigra::test_suite("MultiArrayPointOperatorsTestSuite") {
		add( testCase( &MultiArrayPointoperatorsTest::testCopy ));
	}
};

int main(int argc, char ** argv) {

	MultiArrayPointOperatorsTestSuite test2;
	int failed =  test2.run(vigra::testsToBeExecuted(argc, argv));
	std::cout << test2.report() << std::endl;

	return (failed != 0);
}
