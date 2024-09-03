/************************************************************************/
/*                                                                      */
/*                 Copyright 2013 by Benjamin Seppke                    */      
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <cstdio>

#include "vigra/unittest.hxx"

#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"

#include "vigra/medianfilter.hxx"
#include "vigra/shockfilter.hxx"
#include "vigra/specklefilters.hxx"

using namespace vigra;

template <class T>
void printBasicImage(const BasicImage<T> & img)
{    
    for (int y=0; y<img.height(); ++y)
    {
        for (int x=0; x<img.width(); ++x)
        {
            std::printf("img(%d,%d) = %10.3f; ", x , y , img(x,y));
        }      
        std::printf("\n");  
    }
}


struct MedianFilterEssentialTest
{
    FImage img;
    Diff2D filterShape;
    
    MedianFilterEssentialTest()
    : img(10,10),
      filterShape(3,3)
    {
    }
    
    void testZeros()
    {
        img = 0;
        FImage result(img.size());
        
        //testing all reflect modes, none should produce a difference for an
        //image consisting only of zeros..
        
        //1. BORDER_TREATMENT_AVOID
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_AVOID);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //2. BORDER_TREATMENT_REPEAT
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REPEAT);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //3. BORDER_TREATMENT_REFLECT
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REFLECT);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //4. BORDER_TREATMENT_WRAP
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_WRAP);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
    
        //5. BORDER_TREATMENT_ZEROPAD
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_ZEROPAD);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
    }
    
    void testOnes()
    {
        float value=1;
        img = value;
        FImage result(img.size());
                
        //testing all reflect modes, none should produce a difference for an
        //image consisting only of zeros..
        
        //1. BORDER_TREATMENT_AVOID
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_AVOID);        
        border_zero_checker(result, value,false);
        
        //2. BORDER_TREATMENT_REPEAT
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REPEAT);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //3. BORDER_TREATMENT_REFLECT
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REFLECT);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //4. BORDER_TREATMENT_WRAP
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_WRAP);        
        shouldEqualSequence(img.begin(),img.end(), result.begin());
        
        //5. BORDER_TREATMENT_ZEROPAD
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_ZEROPAD);  
        border_zero_checker(result, value, true);
    }
    
    
    void border_zero_checker(const FImage & result, float value, bool corners_only=false)
    {
        for (int y=0; y!=result.height(); ++y)
        {
            bool top    = y==0, bottom = y==result.height()-1;
            
            for (int x=0; x!=result.width(); ++x)
            {
                bool left  = x==0, right = x==result.width()-1;
                
                //std::cerr << "(" << x << ", " << y << ") -> "<< result(x,y) << std::endl;
                if (left || right || top || bottom)
                {
                    if (corners_only && (left || right) && (top || bottom))
                    {
                        shouldEqual(result(x,y), 0);                        
                    }
                    else if (!corners_only) 
                    {
                        shouldEqual(result(x,y), 0);
                    }
                }
                else
                {
                    shouldEqual(result(x,y), value);
                }
            }
        }
    }
    
};


struct MedianFilterExactTest
{
    FImage img;
    FImage ref_img;
    Diff2D filterShape;
    
    MedianFilterExactTest()
    : img(5,5),
      ref_img(5,5),
      filterShape(3,3)
    {
        img(0,0) = 1; img(0,1) = 2; img(0,2) = 3; img(0,3) = 2; img(0,4) = 1; 
        img(1,0) = 2; img(1,1) = 3; img(1,2) = 4; img(1,3) = 3; img(1,4) = 2; 
        img(2,0) = 3; img(2,1) = 4; img(2,2) = 5; img(2,3) = 4; img(2,4) = 3; 
        img(3,0) = 2; img(3,1) = 3; img(3,2) = 4; img(3,3) = 3; img(3,4) = 2; 
        img(4,0) = 1; img(4,1) = 2; img(4,2) = 3; img(4,3) = 2; img(4,4) = 1; 
    }
    
    void testAVOID()
    {
        FImage result(img.size());
        
        ref_img(0,0) =      0.000; ref_img(1,0) =      0.000; ref_img(2,0) =      0.000; ref_img(3,0) =      0.000; ref_img(4,0) =      0.000; 
        ref_img(0,1) =      0.000; ref_img(1,1) =      3.000; ref_img(2,1) =      3.000; ref_img(3,1) =      3.000; ref_img(4,1) =      0.000; 
        ref_img(0,2) =      0.000; ref_img(1,2) =      3.000; ref_img(2,2) =      4.000; ref_img(3,2) =      3.000; ref_img(4,2) =      0.000; 
        ref_img(0,3) =      0.000; ref_img(1,3) =      3.000; ref_img(2,3) =      3.000; ref_img(3,3) =      3.000; ref_img(4,3) =      0.000; 
        ref_img(0,4) =      0.000; ref_img(1,4) =      0.000; ref_img(2,4) =      0.000; ref_img(3,4) =      0.000; ref_img(4,4) =      0.000; 
        
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_AVOID);        
        shouldEqualSequence(result.begin(), result.end(), ref_img.begin());
    }
    
    void testREPEAT()
    {
        FImage result(img.size());
        
        ref_img(0,0) =      2.000; ref_img(1,0) =      2.000; ref_img(2,0) =      3.000; ref_img(3,0) =      2.000; ref_img(4,0) =      2.000; 
        ref_img(0,1) =      2.000; ref_img(1,1) =      3.000; ref_img(2,1) =      3.000; ref_img(3,1) =      3.000; ref_img(4,1) =      2.000; 
        ref_img(0,2) =      3.000; ref_img(1,2) =      3.000; ref_img(2,2) =      4.000; ref_img(3,2) =      3.000; ref_img(4,2) =      3.000; 
        ref_img(0,3) =      2.000; ref_img(1,3) =      3.000; ref_img(2,3) =      3.000; ref_img(3,3) =      3.000; ref_img(4,3) =      2.000; 
        ref_img(0,4) =      2.000; ref_img(1,4) =      2.000; ref_img(2,4) =      3.000; ref_img(3,4) =      2.000; ref_img(4,4) =      2.000; 
    
        
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REPEAT);        
        shouldEqualSequence(result.begin(), result.end(), ref_img.begin());
    }
    
    void testREFLECT()
    {
        FImage result(img.size());
        
        ref_img(0,0) =      2.000; ref_img(1,0) =      2.000; ref_img(2,0) =      3.000; ref_img(3,0) =      2.000; ref_img(4,0) =      2.000; 
        ref_img(0,1) =      2.000; ref_img(1,1) =      3.000; ref_img(2,1) =      3.000; ref_img(3,1) =      3.000; ref_img(4,1) =      2.000; 
        ref_img(0,2) =      3.000; ref_img(1,2) =      3.000; ref_img(2,2) =      4.000; ref_img(3,2) =      3.000; ref_img(4,2) =      3.000; 
        ref_img(0,3) =      2.000; ref_img(1,3) =      3.000; ref_img(2,3) =      3.000; ref_img(3,3) =      3.000; ref_img(4,3) =      2.000; 
        ref_img(0,4) =      2.000; ref_img(1,4) =      2.000; ref_img(2,4) =      3.000; ref_img(3,4) =      2.000; ref_img(4,4) =      2.000; 
        
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_REFLECT);    
        shouldEqualSequence(result.begin(), result.end(), ref_img.begin());
    }
    void testWRAP()
    {
        FImage result(img.size());
        
        ref_img(0,0) =      2.000; ref_img(1,0) =      2.000; ref_img(2,0) =      3.000; ref_img(3,0) =      2.000; ref_img(4,0) =      2.000; 
        ref_img(0,1) =      2.000; ref_img(1,1) =      3.000; ref_img(2,1) =      3.000; ref_img(3,1) =      3.000; ref_img(4,1) =      2.000; 
        ref_img(0,2) =      3.000; ref_img(1,2) =      3.000; ref_img(2,2) =      4.000; ref_img(3,2) =      3.000; ref_img(4,2) =      3.000; 
        ref_img(0,3) =      2.000; ref_img(1,3) =      3.000; ref_img(2,3) =      3.000; ref_img(3,3) =      3.000; ref_img(4,3) =      2.000; 
        ref_img(0,4) =      2.000; ref_img(1,4) =      2.000; ref_img(2,4) =      3.000; ref_img(3,4) =      2.000; ref_img(4,4) =      2.000; 
        
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_WRAP);        
        shouldEqualSequence(result.begin(), result.end(), ref_img.begin());
    }
    void testZEROPAD()
    {
        FImage result(img.size());
        
        ref_img(0,0) =      0.000; ref_img(1,0) =      2.000; ref_img(2,0) =      2.000; ref_img(3,0) =      2.000; ref_img(4,0) =      0.000; 
        ref_img(0,1) =      2.000; ref_img(1,1) =      3.000; ref_img(2,1) =      3.000; ref_img(3,1) =      3.000; ref_img(4,1) =      2.000; 
        ref_img(0,2) =      2.000; ref_img(1,2) =      3.000; ref_img(2,2) =      4.000; ref_img(3,2) =      3.000; ref_img(4,2) =      2.000; 
        ref_img(0,3) =      2.000; ref_img(1,3) =      3.000; ref_img(2,3) =      3.000; ref_img(3,3) =      3.000; ref_img(4,3) =      2.000; 
        ref_img(0,4) =      0.000; ref_img(1,4) =      2.000; ref_img(2,4) =      2.000; ref_img(3,4) =      2.000; ref_img(4,4) =      0.000; 
        
        medianFilter(srcImageRange(img), destImage(result), filterShape, BORDER_TREATMENT_ZEROPAD);   
        shouldEqualSequence(result.begin(), result.end(), ref_img.begin());
    }
    
};

struct MedianFilterTestSuite
: public vigra::test_suite
{
    MedianFilterTestSuite()
    : vigra::test_suite("MedianFilterTestSuite")
    {
        add( testCase( &MedianFilterEssentialTest::testZeros));
        add( testCase( &MedianFilterEssentialTest::testOnes));
        add( testCase( &MedianFilterExactTest::testAVOID));
        add( testCase( &MedianFilterExactTest::testREPEAT));
        add( testCase( &MedianFilterExactTest::testREFLECT));
        add( testCase( &MedianFilterExactTest::testWRAP));
        add( testCase( &MedianFilterExactTest::testZEROPAD));
   }
};


struct ShockFilterTest
{
    FImage img;
    
    ShockFilterTest()
    : img(100,100)
    {
    }
    
    void test()
    {
        img = 0;
        FImage result(img.size());
        
        float sigma = 0.7f;
        float rho   = 3.0f;
        float upwind_factor_h = 0.3f;
        unsigned int iterations = 10;
        
        //Just test, if it's running with properly set parameters
        shockFilter(srcImageRange(img), destImage(result), sigma, rho, upwind_factor_h, iterations);
    }
};


struct ShockFilterTestSuite
: public vigra::test_suite
{
    ShockFilterTestSuite()
    : vigra::test_suite("ShockFilterTestSuite")
    {
        add( testCase( &ShockFilterTest::test));
    }
};


struct SpeckleFilterTest
{
    FImage img;
    Diff2D filterShape;
    
    SpeckleFilterTest()
    : img(100,100),
      filterShape(3,3)
    {
    }
    
    void testFrostFilter()
    {
        img = 0;
        FImage result(img.size());
        
        float k = 0.7f;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        frostFilter(srcImageRange(img), destImage(result), filterShape, k, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        frostFilter(srcImageRange(img), destImage(result), filterShape, k, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        frostFilter(srcImageRange(img), destImage(result), filterShape, k, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        frostFilter(srcImageRange(img), destImage(result), filterShape, k, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        frostFilter(srcImageRange(img), destImage(result), filterShape, k, BORDER_TREATMENT_ZEROPAD);
    }
    void testEnhancedFrostFilter()
    {
        img = 0;
        FImage result(img.size());
        
        float k = 0.5;
        unsigned int enl = 10;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        enhancedFrostFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        enhancedFrostFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        enhancedFrostFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        enhancedFrostFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        enhancedFrostFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_ZEROPAD);
    }
    
    void testGammaMAPFilter()
    {
        img = 0;
        FImage result(img.size());
        
        unsigned int enl = 10;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        gammaMAPFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        gammaMAPFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        gammaMAPFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        gammaMAPFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        gammaMAPFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_ZEROPAD);
    }
    
    
    void testKuanFilter()
    {
        img = 0;
        FImage result(img.size());
        
        unsigned int enl = 10;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        kuanFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        kuanFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        kuanFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        kuanFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        kuanFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_ZEROPAD);
    }
    
    void testLeeFilter()
    {
        img = 0;
        FImage result(img.size());
        
        unsigned int enl = 10;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        leeFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        leeFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        leeFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        leeFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        leeFilter(srcImageRange(img), destImage(result), filterShape, enl, BORDER_TREATMENT_ZEROPAD);
    }
    
    void testEnhancedLeeFilter()
    {
        img = 0;
        FImage result(img.size());
        
        float k = 0.5;
        unsigned int enl = 10;
        
        //Just test, if it's running with properly set parameters
        //1. BORDER_TREATMENT_AVOID
        enhancedLeeFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_AVOID);
        
        //2. BORDER_TREATMENT_REPEAT
        enhancedLeeFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_REPEAT);
        
        //3. BORDER_TREATMENT_REFLECT
        enhancedLeeFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_REFLECT);
        
        //4. BORDER_TREATMENT_WRAP
        enhancedLeeFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_WRAP);
        
        //5. BORDER_TREATMENT_ZEROPAD
        enhancedLeeFilter(srcImageRange(img), destImage(result), filterShape, k, enl, BORDER_TREATMENT_ZEROPAD);
    }
};


struct SpeckleFilterTestSuite
: public vigra::test_suite
{
    SpeckleFilterTestSuite()
    : vigra::test_suite("SpeckleFilterTestSuite")
    {
        add( testCase( &SpeckleFilterTest::testFrostFilter));
        add( testCase( &SpeckleFilterTest::testEnhancedFrostFilter));
        add( testCase( &SpeckleFilterTest::testGammaMAPFilter));
        add( testCase( &SpeckleFilterTest::testKuanFilter));
        add( testCase( &SpeckleFilterTest::testLeeFilter));
        add( testCase( &SpeckleFilterTest::testEnhancedLeeFilter));
    }
};

struct FilterTestCollection
: public vigra::test_suite
{
    FilterTestCollection()
    : vigra::test_suite("FilterTestCollection")
    {
        add( new MedianFilterTestSuite);
        add( new ShockFilterTestSuite);
        add( new SpeckleFilterTestSuite);
   }
};

int main(int argc, char ** argv)
{
    FilterTestCollection test;
 
    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

