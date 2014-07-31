/************************************************************************/
/*                                                                      */
/*                 Copyright 2013 by Benjamin Seppke                    */      
/*                                                                      */
/************************************************************************/

#include <iostream>

#include "vigra/unittest.hxx"

#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"

#include "vigra/correlation.hxx"

using namespace vigra;

static double test_epsilon = 1.0e-5;
static double test_vs_epsilon = 1.0e-3;

template <class T>
void printBasicImage(const BasicImage<T> & img)
{    
    for (int y=0; y<img.height(); ++y)
    {
        for (int x=0; x<img.width(); ++x)
        {
            printf("img(%d,%d) = %1.10f; ", x , y , img(x,y));
        }      
        printf("\n");  
    }
}


struct FastVsSlowCorrelationTest
{
    FImage img;
    FImage mask;
    
    FastVsSlowCorrelationTest()
    : img(10,10),
      mask(5,5)
    {
        FImage::iterator iter = img.begin();
        
        for ( ; iter != img.end(); ++iter)
        {
            *iter =  rand() % 10 + 1;
        }
        
        for (iter = mask.begin(); iter != mask.end(); ++iter)
        {
            *iter =  rand() % 10 + 1;
        }
    }
    
    
    void testCorrelation()
    {
        FImage result_slow(10,10);
        FImage result_fast(10,10);
        
        for(int mask_size=3; mask_size<10; mask_size+=2)
        {
            fastCrossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(mask_size,mask_size), img.accessor()),
                                            srcImageRange(img),
                                            destImage(result_fast));
            
            crossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(mask_size,mask_size), img.accessor()),
                                        srcImageRange(img),
                                        destImage(result_slow));
            
            shouldEqualSequenceTolerance(result_fast.begin(), result_fast.end(), result_slow.begin(), test_vs_epsilon);
        }
    }
    
    
    void testNormalizedCorrelation()
    {
        FImage result_slow(10,10);
        FImage result_fast(10,10);
        
        for(int mask_size=3; mask_size<10; mask_size+=2)
        {
            fastNormalizedCrossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(mask_size,mask_size), img.accessor()),
                                                      srcImageRange(img),
                                                      destImage(result_fast));
            
            normalizedCrossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(mask_size,mask_size), img.accessor()),
                                                  srcImageRange(img),
                                                  destImage(result_slow));
            
            shouldEqualSequenceTolerance(result_fast.begin(), result_fast.end(), result_slow.begin(), test_vs_epsilon);
        }
    }
};


struct FastNormalizedCrossCorrelationEssentialTest
{
    FImage img;
    FImage mask;
    
    FastNormalizedCrossCorrelationEssentialTest()
    : img(10,10),
      mask(5,5)
    {
        FImage::iterator iter = img.begin();
        
        for ( ; iter != img.end(); ++iter)
        {
            *iter =  rand() % 10 + 1;
        }
        
        for (iter = mask.begin(); iter != mask.end(); ++iter)
        {
            *iter =  rand() % 10 + 1;
        }
    }
    
    void testImagePatch()
    {
        FImage result(10,10);
        
        for(int mask_size=3; mask_size<10; mask_size+=2)
        {
            fastNormalizedCrossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(mask_size,mask_size), img.accessor()),
                                                      srcImageRange(img),
                                                      destImage(result));
            
            FindMinMax<float> minmax;
            inspectImage(srcImageRange(result),minmax);
            should(minmax.min >= -1.0); //correlation coeff. minimum should always be >= 0.0
            should(minmax.max <=  1.0); //correlation coeff. maximum should always be <= 1.0
            
            shouldEqualTolerance(result(mask_size/2,mask_size/2), 1.0, test_epsilon);
        }
    }
    void testRandomPatch()
    {
        FImage result(10,10);
        
        fastNormalizedCrossCorrelation(srcIterRange(img.upperLeft(), img.upperLeft() + Diff2D(3,3), img.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        FindMinMax<float> minmax;
        inspectImage(srcImageRange(result),minmax);
        should(minmax.min >= -1.0); //correlation coeff. minimum should always be >= 0.0
        should(minmax.max <=  1.0); //correlation coeff. maximum should always be <= 1.0
    }
};
 
struct FastNormalizedCrossCorrelationExactTest
{
    FImage img;
    FImage mask;
            
    FastNormalizedCrossCorrelationExactTest()
    : img(7,7),
      mask(5,5)
    {
        img(0,0) = 1;  img(1,0) = 1;  img(2,0) = 1;  img(3,0) = 1;  img(4,0) = 1;  img(5,0) = 1;  img(6,0) = 1;
        img(0,1) = 1;  img(1,1) = 2;  img(2,1) = 2;  img(3,1) = 2;  img(4,1) = 2;  img(5,1) = 2;  img(6,1) = 1;
        img(0,2) = 1;  img(1,2) = 2;  img(2,2) = 3;  img(3,2) = 3;  img(4,2) = 3;  img(5,2) = 2;  img(6,2) = 1;
        img(0,3) = 1;  img(1,3) = 2;  img(2,3) = 3;  img(3,3) = 4;  img(4,3) = 3;  img(5,3) = 2;  img(6,3) = 1;
        img(0,4) = 1;  img(1,4) = 2;  img(2,4) = 3;  img(3,4) = 3;  img(4,4) = 3;  img(5,4) = 2;  img(6,4) = 1;
        img(0,5) = 1;  img(1,5) = 2;  img(2,5) = 2;  img(3,5) = 2;  img(4,5) = 2;  img(5,5) = 2;  img(6,5) = 1;
        img(0,6) = 1;  img(1,6) = 1;  img(2,6) = 1;  img(3,6) = 1;  img(4,6) = 1;  img(5,6) = 1;  img(6,6) = 1;
        
        mask(0,0) = 3;  mask(1,0) = 3;  mask(2,0) = 3;  mask(3,0) = 3;   mask(4,0) = 3;  
        mask(0,1) = 3;  mask(1,1) = 4;  mask(2,1) = 3;  mask(3,1) = 3;   mask(4,1) = 3;  
        mask(0,2) = 3;  mask(1,2) = 3;  mask(2,2) = 3;  mask(3,2) = 3;   mask(4,2) = 3;  
        mask(0,3) = 3;  mask(1,3) = 3;  mask(2,3) = 3;  mask(3,3) = 3;   mask(4,3) = 3;  
        mask(0,4) = 3;  mask(1,4) = 3;  mask(2,4) = 3;  mask(3,4) = 3;   mask(4,4) = 3;  
    }
    
    void testMask1x1()
    {
        FImage result(7,7), ref_img(7,7);
        result = 0.0;
        
        //Nothing is more or less equal to a one element mask
        ref_img = 0.0;
                
        fastNormalizedCrossCorrelation(srcIterRange(mask.upperLeft(), mask.upperLeft() + Diff2D(1,1), mask.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        shouldEqualSequenceTolerance(result.begin(), result.end(), ref_img.begin(), test_epsilon);
    }
    
    void testMask3x3()
    {
        FImage result(7,7), ref_img(7,7);
        result = 0.0;
        ref_img = 0.0;
        
        ::fastNormalizedCrossCorrelation(srcIterRange(mask.upperLeft(), mask.upperLeft() + Diff2D(3,3), mask.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        ref_img(1,1) = 0.2294157356f; ref_img(2,1) = 0.0533001795f; ref_img(3,1) = 0.0000000000f; ref_img(4,1) = 0.0533001795f; ref_img(5,1) = 0.2294157356f;
        ref_img(1,2) = 0.0533001795f; ref_img(2,2) = 0.2294157356f; ref_img(3,2) = 0.1250000009f; ref_img(4,2) = 0.2294157356f; ref_img(5,2) = 0.0533001795f;
        ref_img(1,3) = 0.0000000000f; ref_img(2,3) = 0.1250000009f; ref_img(3,3) = 1.0000000075f; ref_img(4,3) = 0.1250000009f; ref_img(5,3) = 0.0000000000f;
        ref_img(1,4) = 0.0533001795f; ref_img(2,4) = 0.2294157356f; ref_img(3,4) = 0.1250000009f; ref_img(4,4) = 0.2294157356f; ref_img(5,4) = 0.0533001795f;
        ref_img(1,5) = 0.2294157356f; ref_img(2,5) = 0.0533001795f; ref_img(3,5) = 0.0000000000f; ref_img(4,5) = 0.0533001795f; ref_img(5,5) = 0.2294157356f;
        
        shouldEqualSequenceTolerance(result.begin(), result.end(), ref_img.begin(), test_epsilon);
    }
    
    
    void testMask3x5()
    {
        FImage result(7,7), ref_img(7,7);
        result = 0.0;
        ref_img = 0.0;
        
        fastNormalizedCrossCorrelation(srcIterRange(mask.upperLeft() + Diff2D(1,0), mask.upperLeft() + Diff2D(4,5), mask.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        ref_img(1,2) = -0.2539664209f; ref_img(2,2) = -0.0834784210f; ref_img(3,2) = -0.1410190463f; ref_img(4,2) = -0.0834784210f; ref_img(5,2) = 0.0923513919f; 
        ref_img(1,3) = -0.3225896060f; ref_img(2,3) = -0.2017366886f; ref_img(3,3) =  0.1494035274f; ref_img(4,3) =  0.2305561155f; ref_img(5,3) = 0.4218478799f;
        ref_img(1,4) = -0.2539664209f; ref_img(2,4) = -0.0834784210f; ref_img(3,4) =  0.1611645669f; ref_img(4,4) =  0.5426095128f; ref_img(5,4) = 0.4386692047f; 
        
        shouldEqualSequenceTolerance(result.begin(), result.end(), ref_img.begin(), test_epsilon);
    }
    
    
    void testMask5x3()
    {
        FImage result(7,7), ref_img(7,7);
        result = 0.0;
        ref_img = 0.0;
        
        fastNormalizedCrossCorrelation(srcIterRange(mask.upperLeft() + Diff2D(0,1) , mask.upperLeft() + Diff2D(5,4), mask.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        ref_img(2,1) = -0.2539664209f; ref_img(3,1) = -0.3225896060f; ref_img(4,1) = -0.2539664209f; 
        ref_img(2,2) = -0.0834784210f; ref_img(3,2) = -0.2017366886f; ref_img(4,2) = -0.0834784210f; 
        ref_img(2,3) = -0.1410190463f; ref_img(3,3) =  0.1494035274f; ref_img(4,3) =  0.1611645669f;
        ref_img(2,4) = -0.0834784210f; ref_img(3,4) =  0.2305561155f; ref_img(4,4) =  0.5426095128f;
        ref_img(2,5) =  0.0923513919f; ref_img(3,5) =  0.4218478799f; ref_img(4,5) =  0.4386692047f;
        
        shouldEqualSequenceTolerance(result.begin(), result.end(), ref_img.begin(), test_epsilon);
    }
    
    void testMask5x5()
    {
        FImage result(7,7), ref_img(7,7);
        result = 0.0;
        ref_img = 0.0;
        
        fastNormalizedCrossCorrelation(srcIterRange(mask.upperLeft(), mask.upperLeft() + Diff2D(5,5), mask.accessor()),
                                                 srcImageRange(img),
                                                 destImage(result));
        
        
        ref_img(2,2) = -0.0089172045f; ref_img(3,2) = -0.0510310352f; ref_img(4,2) = -0.0089172045f; 
        ref_img(2,3) = -0.0510310352f; ref_img(3,3) =  0.2165063461f; ref_img(4,3) =  0.2041241407f;
        ref_img(2,4) = -0.0089172045f; ref_img(3,4) =  0.2041241407f; ref_img(4,4) =  0.4369430199f;
        
        shouldEqualSequenceTolerance(result.begin(), result.end(), ref_img.begin(), test_epsilon);
    }
};


struct FastNormalizedCrossCorrelationTestSuite
: public test_suite
{
    FastNormalizedCrossCorrelationTestSuite()
    : test_suite("FastNormalizedCrossCorrelationTestSuite")
    { 
        add( testCase( &FastVsSlowCorrelationTest::testCorrelation));
        add( testCase( &FastVsSlowCorrelationTest::testNormalizedCorrelation));
        
        add( testCase( &FastNormalizedCrossCorrelationEssentialTest::testImagePatch));
        add( testCase( &FastNormalizedCrossCorrelationEssentialTest::testRandomPatch));
        
        add( testCase( &FastNormalizedCrossCorrelationExactTest::testMask1x1));
        add( testCase( &FastNormalizedCrossCorrelationExactTest::testMask3x3));
        add( testCase( &FastNormalizedCrossCorrelationExactTest::testMask3x5));
        add( testCase( &FastNormalizedCrossCorrelationExactTest::testMask5x3));
        add( testCase( &FastNormalizedCrossCorrelationExactTest::testMask5x5));
   }
};


struct CorrelationTestCollection
: public test_suite
{
    CorrelationTestCollection()
    : test_suite("CorrelationTestCollection")
    {
        add( new FastNormalizedCrossCorrelationTestSuite);
   }
};

int main(int argc, char ** argv)
{
    CorrelationTestCollection test;
 
    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

