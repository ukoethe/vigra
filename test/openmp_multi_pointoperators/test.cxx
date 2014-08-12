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
#include "vigra/error.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/multi_math.hxx"
#include "vigra/openmp_vigra.h"
#include "vigra/timing.hxx"

using namespace std;
using namespace vigra;

#define DIMENSION_LENGTH 60

struct OpenMPMultiArrayPointOperatorsTest {
    typedef float PixelType;
    typedef MultiArray<5, PixelType> Image5D;
    typedef Image5D::difference_type Size5;
    Image5D img;
    int max_threads;

    OpenMPMultiArrayPointOperatorsTest() :
        img(Size5( DIMENSION_LENGTH, DIMENSION_LENGTH,
                   DIMENSION_LENGTH, DIMENSION_LENGTH, DIMENSION_LENGTH))
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


//TODO: test function returns number of threads involved in parallel region
//TODO: Meaningful runtime.
    void testInit()
        {
            Image5D res(img.shape());
            const Image5D::value_type ini = 1.1f;

            should(res.shape() == Size5( DIMENSION_LENGTH, DIMENSION_LENGTH,
                                         DIMENSION_LENGTH, DIMENSION_LENGTH,
                                         DIMENSION_LENGTH ));

            USETICTOC
            TIC
            initMultiArray(res, ini);
            TOC

            using namespace multi_math;
            should(all(res == ini));

            initMultiArray(res, 2.2f);
            should(all(res == 2.2f));

            res = 3.3f;
            should(all(res == 3.3f));

            res.init(4.4f);
            should(all(res == 4.4f));
        }


    void testCopy() {

        Image5D res(img.shape(), 1.0);

        USETICTOC
        TIC
        vigra::copyMultiArray(img, res);
        TOC

        should(img == res);
    }


    void testTransform()
    {
        using namespace vigra::functor;

        Image5D res(img.shape());

        USETICTOC
        TIC
        transformMultiArray(img, res, Arg1() + Arg1());
        TOC

        using namespace multi_math;
        should(all(2.0*img == res));
    }


    void testTransformOuterExpand()
    {
        using namespace functor;

        Image5D res(img.shape());

        USETICTOC
        TIC
        transformMultiArray(img.subarray(Size5(0,0,0,0,0), Size5(DIMENSION_LENGTH,1,1,1,1)), res,
                                Arg1() + Arg1());
        TOC

        int x,y,z,t,k;
        for(k=0; k<img.shape(4); ++k)
            for(t=0; t<img.shape(3); ++t)
                for(z=0; z<img.shape(2); ++z)
                    for(y=0; y<img.shape(1); ++y)
                        for(x=0; x<img.shape(0); ++x)
                            shouldEqual(res(x,y,z,t,k), 2.0*img(x,0,0,0,0));
    }
};

struct MultiArrayPointOperatorsTestSuite: public vigra::test_suite {
    MultiArrayPointOperatorsTestSuite() :
        vigra::test_suite("MultiArrayPointOperatorsTestSuite") {
        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testInit ) );
        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testCopy ));
        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testTransform ) );
        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testTransformOuterExpand ) );
//        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testCombine2 ) );
//        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testCombine3 ) );
//        add( testCase( &OpenMPMultiArrayPointOperatorsTest::testInspect ) );
    }
};

int main(int argc, char ** argv) {
    MultiArrayPointOperatorsTestSuite test;
    int failed =  test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return (failed != 0);
}
