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

#include <iostream>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/flatmorphology.hxx"
#include "vigra/multi_array.hxx"

using namespace vigra;

struct FlatMorphologyTest
{
    typedef vigra::BImage Image;
    typedef vigra::MultiArrayView<2, unsigned char> View;

    FlatMorphologyTest()
    : img(7,7), mask(7,7)
    {
        static const unsigned char in[] = {
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6,
            0, 1, 2, 3, 4, 5, 6};
        
        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const unsigned char * p = in;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }

        static const unsigned char in1[] = {
            1, 1, 1, 1, 1, 1, 1,
            0, 1, 1, 1, 1, 1, 0,
            0, 1, 1, 1, 1, 1, 0,
            0, 1, 1, 1, 1, 1, 0,
            0, 1, 1, 1, 1, 1, 0,
            0, 1, 1, 1, 1, 1, 0,
            1, 1, 1, 1, 1, 1, 1};
        
        i = mask.begin();
        end = mask.end();
        p = in1;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }
    
    void erosionTest()
    {
        Image res(img.size()), res1(img.size());
        
        discErosion(srcImageRange(img), destImage(res), 2);
        
        static const unsigned char desired[] = {
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discErosion(View(img), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    void erosionWithMaskTest()
    {
        Image res(img.size(), 9), res1(img.size(), 9);
        
        discErosionWithMask(srcImageRange(img), maskImage(mask),
                            destImage(res), 2);
        
        static const unsigned char desired[] = {
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 1, 1, 2, 3, 4,
                   1, 1, 1, 1, 2, 3, 4,
                   0, 0, 1, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4,
                   0, 0, 0, 1, 2, 3, 4};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discErosionWithMask(View(img), View(mask), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    void dilationTest()
    {
        Image res(img.size()), res1(img.size());
        
        discDilation(srcImageRange(img), destImage(res), 2);
        
        static const unsigned char desired[] = {
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discDilation(View(img), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    void dilationWithMaskTest()
    {
        Image res(img.size(), 9), res1(img.size(), 9);
        
        discDilationWithMask(srcImageRange(img), maskImage(mask),
                            destImage(res), 2);
        
        static const unsigned char desired[] = {
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 5, 6, 6,
                   2, 3, 4, 5, 5, 5, 5,
                   2, 3, 4, 5, 5, 6, 6,
                   2, 3, 4, 5, 6, 6, 6,
                   2, 3, 4, 5, 6, 6, 6};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discDilationWithMask(View(img), View(mask), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    void medianTest()
    {
        Image res(img.size()), res1(img.size());
        
        discMedian(srcImageRange(img), destImage(res), 2);
        
        static const unsigned char desired[] = {
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5,
            1, 1, 2, 3, 4, 5, 5};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discMedian(View(img), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    void medianWithMaskTest()
    {
        Image res(img.size(), 9), res1(img.size(), 9);
        
        discMedianWithMask(srcImageRange(img), maskImage(mask),
                            destImage(res), 2);
        
        static const unsigned char desired[] = {
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5,
                   1, 2, 2, 3, 4, 4, 5};
                                         
        const unsigned char * i1 = desired;
        const unsigned char * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }

        discMedianWithMask(View(img), View(mask), View(res1), 2);
        should(View(res) == View(res1));
    }
    
    Image img, mask;
};

        
struct MorphologyTestSuite
: public vigra::test_suite
{
    MorphologyTestSuite()
    : vigra::test_suite("MorphologyTestSuite")
    {
        add( testCase( &FlatMorphologyTest::erosionTest));
        add( testCase( &FlatMorphologyTest::erosionWithMaskTest));
        add( testCase( &FlatMorphologyTest::dilationTest));
        add( testCase( &FlatMorphologyTest::dilationWithMaskTest));
        add( testCase( &FlatMorphologyTest::medianTest));
        add( testCase( &FlatMorphologyTest::medianWithMaskTest));
    }
};

int main(int argc, char ** argv)
{
    MorphologyTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

