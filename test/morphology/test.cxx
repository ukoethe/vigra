#include <iostream>
#include "unittest.h"
#include "vigra/stdimage.hxx"
#include "vigra/flatmorphology.hxx"

struct FlatMorphologyTest
{
    typedef BImage Image;

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
        Image res(img);
        
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
    }
    
    void erosionWithMaskTest()
    {
        Image res(img);
        res = 9;
        
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
    }
    
    void dilationTest()
    {
        Image res(img);
        
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
    }
    
    void dilationWithMaskTest()
    {
        Image res(img);
        res = 9;
        
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
    }
    
    void medianTest()
    {
        Image res(img);
        res = 0;
        
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
    }
    
    void medianWithMaskTest()
    {
        Image res(img);
        res = 9;
        
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
    }
    
    Image img, mask;
};

        
struct MorphologyTestSuite
: public TestSuite
{
    MorphologyTestSuite()
    : TestSuite("MorphologyTestSuite")
    {
        add( testCase( &FlatMorphologyTest::erosionTest));
        add( testCase( &FlatMorphologyTest::erosionWithMaskTest));
        add( testCase( &FlatMorphologyTest::dilationTest));
        add( testCase( &FlatMorphologyTest::dilationWithMaskTest));
        add( testCase( &FlatMorphologyTest::medianTest));
        add( testCase( &FlatMorphologyTest::medianWithMaskTest));
    }
};

int main()
{
    MorphologyTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

