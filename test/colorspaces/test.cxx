#include <algorithm>
#include "unittest.h"
#include "vigra/colorconversions.hxx"

struct ColorConversionsTest
{
    typedef vigra::RGBValue<double> RGB;
    typedef vigra::TinyVector<double, 3> Color;
    typedef vigra::TinyVector<float, 3> (*PolarFct)(vigra::TinyVector<float, 3> const &);
    enum { count = 15 };
    RGB black, red, green, blue, cyan, magenta, 
         yellow, white, gray, orange, violet, bordeaux, brown, bluegreen, darkgreen;
    Color original[count], transformed[count], back[count];

#if 0
    tblack, tred, tgreen, tblue, tcyan, tmagenta, 
         tyellow, twhite, tgray, torange, tviolet, tbordeaux, tbrown, tbluegreen;
    FRGB bblack, bred, bgreen, bblue, bcyan, bmagenta, 
         byellow, bwhite, bgray, borange, bviolet, bbordeaux, bbrown, bbluegreen;
#endif /* #if 0 */

    
    ColorConversionsTest()
    : black(0.0,0.0,0.0), 
      red(255.0, 0.0, 0.0), 
      green(0.0, 255.0, 0.0), 
      blue(0.0, 0.0, 255.0), 
      cyan(0.0, 255.0, 255.0), 
      magenta(255.0, 0.0, 255.0), 
      yellow(255.0, 255.0, 0.0), 
      white(255.0, 255.0, 255.0), 
      gray(128.0, 128.0, 128.0),
      orange(255.0, 150.0, 0.0),
      violet(200.0, 0.0, 255.0),
      bordeaux(255.0, 0.0, 140.0),
      brown(180.0, 90.0, 0.0),
      bluegreen(0.0, 160.0, 178.0),
      darkgreen(90.0, 170.0, 100.0)
    {
        original[0] = black;
        original[1] = red;
        original[2] = green;
        original[3] = blue;
        original[4] = cyan;
        original[5] = magenta;
        original[6] = yellow;
        original[7] = white;
        original[8] = gray;
        original[9] = orange;
        original[10] = violet;
        original[11] = bordeaux;
        original[12] = brown;
        original[13] = bluegreen;
        original[14] = darkgreen;
    }
   
    static bool equalColors(Color const & c1, Color const & c2, double epsilon = 0.001)
    {
        return (c1 - c2).magnitude() < epsilon;
    }
    
    static bool equalColors(Color * c1, Color * c1end, Color * c2, double epsilon = 0.001)
    {
        for(; c1 != c1end; ++c1, ++c2)
        {
            if(!equalColors(*c1, *c2, epsilon))
                return false;
        }
        return true;
    }
    
    static bool checkLimits(Color * c, Color * cend, 
               Color const & min, Color const & max, double epsilon = 0.001)
    {
        for(; c != cend; ++c)
        {
            for(int i=0; i<3; ++i)
                if(min[i] - (*c)[i] > epsilon || (*c)[i] - max[i] > epsilon)
                    return false;
        }
        return true;
    }
    
    static bool checkGray(Color const & c, double epsilon = 0.001)
    {
        return fabs(c[1]) < epsilon && fabs(c[2]) < epsilon;
    }
    
    static void printColor(Color const & c)
    {
        std::cerr << "[" << c[0] << ", " << c[1] << ", " << c[2] << "]\n";
    }
    
    static void printColors(Color * c, Color * cend)
    {
        for(; c != cend; ++c)
        {
            printColor(*c);
        }
    }
    
    void testRGBPrime2RGB()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2RGBFunctor<double, double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::RGB2RGBPrimeFunctor<double, double>());
        
        should(checkLimits(transformed, transformed+count, RGB(0.0,0.0,0.0), RGB(255.0,255.0,255.0)));
        should(equalColors(original, original+count, back));
        should(equalColors(transformed[count-1], RGB(25.202, 103.568, 31.8506)));
    }
    
    void testRGB2XYZ()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGB2XYZFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::XYZ2RGBFunctor<double>());
        
        should(equalColors(original, original+count, back));
        should(equalColors(transformed[count-1], RGB(0.454712, 0.580135, 0.458924)));
    }
    
    void testRGBPrime2XYZ()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2XYZFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::XYZ2RGBPrimeFunctor<double>());
        
        should(equalColors(original, original+count, back, 0.01));
        should(equalColors(transformed[count-1], RGB(0.20853, 0.320495, 0.169008)));
    }
    
    void testRGB2Lab()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGB2LabFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::Lab2RGBFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-86.1813,-107.862), RGB(100.0,98.2352,94.4758)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(80.7463, -25.9546, 16.8469)));
    }    
    
    void testRGBPrime2Lab()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2LabFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::Lab2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-86.1813,-107.862), RGB(100.0,98.2352,94.4758)));
        should(equalColors(original, original+count, back, 0.01));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(63.3838, -40.6056, 29.3815)));
    }    
    
    void testRGB2Luv()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGB2LuvFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::Luv2RGBFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-83.077,-134.101), RGB(100.0,175.015,107.393)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(80.7463, -26.4172, 28.6933)));
    }
    
    void testRGBPrime2Luv()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2LuvFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::Luv2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-83.077,-134.101), RGB(100.0,175.015,107.393)));
        should(equalColors(original, original+count, back, 0.01));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(63.3838, -38.5724, 44.4313)));
    }
    
    void testRGBPrime2YPrimePbPr()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimePbPrFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::YPrimePbPr2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-0.5,-0.5), RGB(1.0,0.5,0.5)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(0.541569, -0.0843182, -0.134542)));
    }
    
    void testRGBPrime2YPrimeCbCr()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeCbCrFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::YPrimeCbCr2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(16.0,16.0,16.0), RGB(235.0,240.0,240.0)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]-RGB(0.0,128.0,128.0)));
        should(checkGray(transformed[7]-RGB(0.0,128.0,128.0)));
        should(checkGray(transformed[8]-RGB(0.0,128.0,128.0)));
        should(equalColors(transformed[count-1], RGB(134.604, 109.113, 97.8627)));
    }
    
    void testRGBPrime2YPrimeIQ()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeIQFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::YPrimeIQ2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-0.596,-0.523), RGB(1.0,0.596,0.523)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(0.541569, -0.0985882, -0.151882)));
    }
    
    void testRGBPrime2YPrimeUV()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeUVFunctor<double>());
        std::transform(transformed, transformed+count, back, 
                        vigra::YPrimeUV2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,-0.436,-0.615), RGB(1.0,0.436,0.615)));
        should(equalColors(original, original+count, back));
        should(checkGray(transformed[0]));
        should(checkGray(transformed[7]));
        should(checkGray(transformed[8]));
        should(equalColors(transformed[count-1], RGB(0.541569, -0.0735254, -0.165463)));
    }
    
    void testLabPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2LabFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::lab2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2Lab);
        std::transform(back, back+count, back,
                        vigra::Lab2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));
        
        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(104.113, 0.633838, 0.374569)));
    }
    
    void testLuvPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2LuvFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::luv2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2Luv);
        std::transform(back, back+count, back,
                        vigra::Luv2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));
        
        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(118.79, 0.633838, 0.328633)));
    }
    
    void testYPrimePbPrPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimePbPrFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::yPrimePbPr2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2YPrimePbPr);
        std::transform(back, back+count, back,
                        vigra::YPrimePbPr2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));
        
        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(129.276, 0.541569, 0.297403)));
    }
    
    void testYPrimeCbCrPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeCbCrFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::yPrimeCbCr2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2YPrimeCbCr);
        std::transform(back, back+count, back,
                        vigra::YPrimeCbCr2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));
        
        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(129.276, 0.541569, 0.297403)));
    }
    
    void testYPrimeIQPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeIQFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::yPrimeIQ2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2YPrimeIQ);
        std::transform(back, back+count, back,
                        vigra::YPrimeIQ2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));
        
        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(142.569, 0.541569, 0.286246)));
    }
    
    void testYPrimeUVPolar()
    {
        std::transform(original, original+count, transformed, 
                        vigra::RGBPrime2YPrimeUVFunctor<double>());
        std::transform(transformed, transformed+count, transformed, &vigra::yPrimeUV2Polar);
        std::transform(transformed, transformed+count, back, (PolarFct)&vigra::polar2YPrimeUV);
        std::transform(back, back+count, back,
                        vigra::YPrimeUV2RGBPrimeFunctor<double>());
        
        should(checkLimits(transformed, transformed+count, 
               RGB(0.0,0.0,0.0), RGB(360.0,1.0, 1.0)));
        should(equalColors(original, original+count, back, 0.3));

        should(transformed[0][2] < 0.001);
        should(transformed[7][2] < 0.001);
        should(transformed[8][2] < 0.001);
        should(transformed[1][0] < 0.001 || transformed[1][0] > 359.999);
        
        should(equalColors(transformed[count-1], RGB(142.585, 0.541569, 0.286346)));
    }
};


struct ColorConversionsTestSuite
: public vigra::test_suite
{
                           
    ColorConversionsTestSuite()
    : vigra::test_suite("ColorConversionsTest")
    {
        add( testCase(&ColorConversionsTest::testRGBPrime2RGB));
        add( testCase(&ColorConversionsTest::testRGB2XYZ));
        add( testCase(&ColorConversionsTest::testRGBPrime2XYZ));
        add( testCase(&ColorConversionsTest::testRGB2Lab));
        add( testCase(&ColorConversionsTest::testRGBPrime2Lab));
        add( testCase(&ColorConversionsTest::testRGB2Luv));
        add( testCase(&ColorConversionsTest::testRGBPrime2Luv));
        add( testCase(&ColorConversionsTest::testRGBPrime2YPrimePbPr));
        add( testCase(&ColorConversionsTest::testRGBPrime2YPrimeCbCr));
        add( testCase(&ColorConversionsTest::testRGBPrime2YPrimeIQ));
        add( testCase(&ColorConversionsTest::testRGBPrime2YPrimeUV));
        add( testCase(&ColorConversionsTest::testLabPolar));
        add( testCase(&ColorConversionsTest::testLuvPolar));
        add( testCase(&ColorConversionsTest::testYPrimePbPrPolar));
        add( testCase(&ColorConversionsTest::testYPrimeCbCrPolar));
        add( testCase(&ColorConversionsTest::testYPrimeIQPolar));
        add( testCase(&ColorConversionsTest::testYPrimeUVPolar));
    }
};

int main()
{
    ColorConversionsTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
