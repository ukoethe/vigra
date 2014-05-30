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
#include "vigra/stdimage.hxx"
#include "vigra/openmp_vigra.h"
#include "vigra/multi_array.hxx"

using namespace vigra;

struct OpenMPWrapperTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;

    OpenMPWrapperTest()
    : img(3,3)
    {
        Image::Accessor acc = img.accessor();
        Image::ScanOrderIterator i = img.begin();

        acc.set(1.1, i);
        ++i;
        acc.set(2.2, i);
        ++i;
        acc.set(3.3, i);
        ++i;
        acc.set(4.4, i);
        ++i;
        acc.set(5.5, i);
        ++i;
        acc.set(6.6, i);
        ++i;
        acc.set(7.7, i);
        ++i;
        acc.set(8.8, i);
        ++i;
        acc.set(9.9, i);
        ++i;
        should(i == img.end());

    }

    void additionTest()
    {
        Image img1(3,3), img2(3,3), img3(3,3);

        std::plus<Image::value_type> add;

        combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), add);
        combineTwoImages(View(img), View(img), View(img2), add);

        using namespace functor;
        combineThreeImages(View(img), View(img), View(img1), View(img3), Arg1() + Arg2() + Arg3());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i2 = img2.begin();
        Image::ScanOrderIterator i3 = img3.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1, ++i2, ++i3)
        {
            should(2.0*acc(i) == acc(i1));
            should(2.0*acc(i) == acc(i2));
            should(4.0*acc(i) == acc(i3));
        }
    }

    Image img;
};

struct OpenMPWrapperTestSuite
: public vigra::test_suite
{
    OpenMPWrapperTestSuite()
    : vigra::test_suite("OpenMPWrapperTestSuite")
    {
        add( testCase( &OpenMPWrapperTest::additionTest));
       
    }
};

int main(int argc, char ** argv)
{
    OpenMPWrapperTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
