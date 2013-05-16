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

#include "unittest.hxx"

#include <vigra/slic.hxx>
#include <vigra/impex.hxx>
#include <vigra/convolution.hxx>
#include <vigra/multi_math.hxx>


using namespace vigra;

template <unsigned int N>
struct SlicTest
{
    typedef MultiArray<N, float>              FArray;
    typedef MultiArray<N, int>                IArray;
    typedef typename MultiArrayShape<N>::type Shape;

    ImageImportInfo info;
    FArray lennaImage;

    SlicTest()
    :   info("lenna.xv"),
        lennaImage(info.shape())
    {
        importImage(info, destImage(lennaImage));
    }

    void test_seeding()
    {
        Shape seeds_ref[] = {
               Shape(20, 20),
               Shape(59, 18),
               Shape(98, 18),
               Shape(20, 59),
               Shape(57, 57),
               Shape(97, 59),
               Shape(20, 98),
               Shape(57, 97),
               Shape(98, 98)
        };

        // get grad mag image
        FArray gradMag(lennaImage.shape());
        gaussianGradientMagnitude(srcImageRange(lennaImage),destImage(gradMag),3.0);

        IArray labels(lennaImage.shape());
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
        // get grad mag image
        FArray gradMag(lennaImage.shape());
        gaussianGradientMagnitude(srcImageRange(lennaImage), destImage(gradMag), 3.0);

        IArray labels(lennaImage.shape()), labels_ref(lennaImage.shape());
        int seedDistance = 39;
        generateSlicSeeds(gradMag, labels, seedDistance, 1);

        slicSuperpixels(lennaImage, labels, 10.0, seedDistance);

        //exportImage(srcImageRange(labels), ImageExportInfo("slic.xv"));
        importImage(ImageImportInfo("slic.xv"), destImage(labels_ref));

        using namespace vigra::multi_math;
        should(all(labels == labels_ref));
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

