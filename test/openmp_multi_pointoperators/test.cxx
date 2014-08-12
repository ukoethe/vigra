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
#include "vigra/timing.hxx"

using namespace std;
using namespace vigra;

struct MultiArrayPointoperatorsTest {
    typedef float PixelType;
    typedef MultiArray<5, PixelType> Image5D;
    typedef Image5D::difference_type Size5;
    Image5D img;
    int max_threads;

    MultiArrayPointoperatorsTest() :
        img(Size5( 60, 60, 60, 60, 60))
#ifdef OPENMP
        , max_threads(std::max(omp_get_num_procs() / 2, 2))
#else
        , max_threads(1)
#endif
    {
#ifdef OPENMP
        omp_set_num_threads(max_threads);
#endif
        int i;
        PixelType c = 0.1f;
        for (i = 0; i < img.elementCount(); ++i, ++c)
            img.data()[i] = c;
    }

    void testCopy() {
        USETICTOC
        Image5D res(img.shape(), 1.0);
        TIC
        vigra::copyMultiArray(img, res);
        TOC
        should(img == res);
    }
};

struct MultiArrayPointOperatorsTestSuite: public vigra::test_suite {
    MultiArrayPointOperatorsTestSuite() :
        vigra::test_suite("MultiArrayPointOperatorsTestSuite") {
        add( testCase( &MultiArrayPointoperatorsTest::testCopy ));
    }
};

int main(int argc, char ** argv) {
    MultiArrayPointOperatorsTestSuite test;
    int failed =  test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return (failed != 0);
}
