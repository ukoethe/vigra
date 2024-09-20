/************************************************************************/
/*                                                                      */
/*    Copyright 2012-2013 by Ullrich Koethe and Thorsten Beier          */
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
#include <functional>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "vigra/unittest.hxx"

#include <vigra/slic.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/multi_math.hxx>
#include <vigra/colorconversions.hxx>


using namespace vigra;

template <unsigned int N>
struct SlicTest
{
    typedef MultiArray<N, float>                  FArray;
    typedef MultiArray<N, RGBValue<float> >       FRGBArray;
    typedef MultiArray<N, unsigned int>           IArray;
    typedef typename MultiArrayShape<N>::type     Shape;

    ImageImportInfo info;
    FRGBArray astronautImage;

    SlicTest()
    :   info("astronaut.xv"),
        astronautImage(info.shape())
    {
        importImage(info, destImage(astronautImage));
        transformMultiArray(srcMultiArrayRange(astronautImage), destMultiArray(astronautImage), RGBPrime2LabFunctor<float>());
    }

    void test_seeding()
    {
        Shape seeds_ref[] = {
               Shape(26, 22),
               Shape(65, 22),
               Shape(103, 20),
               Shape(26, 59),
               Shape(65, 59),
               Shape(104, 61),
               Shape(26, 98),
               Shape(65, 100),
               Shape(102, 99)
        };

        // get grad mag image
        FArray gradMag(astronautImage.shape());
        gaussianGradientMagnitude(astronautImage, gradMag, 3.0);

        IArray labels(astronautImage.shape());
        int maxSeedlabel = generateSlicSeeds(gradMag, labels, 39, 1);
        shouldEqual(maxSeedlabel, 9);

        typename IArray::iterator iter = labels.begin(),
                                  end  = iter.getEndIterator();
        int count = 0;
        for(; iter != end; ++iter)
        {
            if(*iter == 0)
                continue;
            should(*iter <= 9 && *iter > 0);
            shouldEqual(iter.point(), seeds_ref[*iter-1]);
            ++count;
        }
        shouldEqual(count, 9);
    }

    void test_slic()
    {
        IArray labels(astronautImage.shape()), labels_ref(astronautImage.shape());

        int seedDistance = 9;
        // compute seeds automatically
        int maxlabel = slicSuperpixels(astronautImage, labels, 20.0, seedDistance, SlicOptions().minSize(0).iterations(40));

        shouldEqual(maxlabel, 201);

        // exportImage(srcImageRange(labels), ImageExportInfo("slic.xv").setForcedRangeMapping(0.0, 1.0, 0.0, 1.0));
        importImage(ImageImportInfo("slic.xv"), destImage(labels_ref));

        should(labels == labels_ref);
    }
};


struct Slic2dTestSuite
: public test_suite
{
    Slic2dTestSuite()
    : test_suite("Slic2dTestSuite")
    {
        add( testCase( &SlicTest<2>::test_seeding));
        add( testCase( &SlicTest<2>::test_slic));
    }
};

int main(int argc, char ** argv)
{
    Slic2dTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

