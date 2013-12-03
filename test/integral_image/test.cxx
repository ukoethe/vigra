/************************************************************************/
/*                                                                      */
/*    Copyright 2012-2013 by Ullrich Koethe and Anna Kreshuk            */
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

#include <vigra/integral_image.hxx>

using namespace vigra;

struct IntegralImageTest
{
    typedef MultiArray<2, int>                  Image2;
    typedef MultiArray<2, TinyVector<int, 2> >  Vector2Image2;
    typedef MultiArray<3, int>                  Image3;
    typedef MultiArray<4, int>                  Image4;

    void test_2d()
    {        
        int desired[] = {
            2,  4,  6,  8,
            4,  8, 12, 16,
            6, 12, 18, 24,
            
            2,  4,  6,  8,
            4,  8, 12, 16,
            6, 12, 18, 24
        };
        
        {
            Image2 in(Shape2(4,3), 2);
            Image2 result(in.shape());
            
            integralMultiArray(in, result);
            shouldEqualSequence(result.begin(), result.end(), desired);
        }
        
        {
            Image3 in(Shape3(4,3,2), 2);
            Image3 result(in.shape());
            
            integralMultiArray(in.multiband(), result.multiband());
            shouldEqualSequence(result.begin(), result.end(), desired);
        }
        
        {
            Image2 in(Shape2(4,3), 2);
            Image2 result(in.shape());
            Image2 desired_squared(in.shape(), desired);
            desired_squared *= 2;
            
            integralMultiArraySquared(in, result);
            shouldEqualSequence(result.begin(), result.end(), desired_squared.begin());
        }
        
        {
            Image3 in(Shape3(4,3,2), 2);
            Image3 result(in.shape());
            Image3 desired_squared(in.shape(), desired);
            desired_squared *= 2;
            
            integralMultiArraySquared(in.multiband(), result.multiband());
            shouldEqualSequence(result.begin(), result.end(), desired_squared.begin());
        }
    }

    void test_3d()
    {
        int desired[] = { 
             2,  4,  6,
             4,  8, 12,
             6, 12, 18,
             
             4,  8, 12,
             8, 16, 24,
            12, 24, 36,
             
             6, 12, 18,
            12, 24, 36,
            18, 36, 54,
            
             2,  4,  6,
             4,  8, 12,
             6, 12, 18,
             
             4,  8, 12,
             8, 16, 24,
            12, 24, 36,
             
             6, 12, 18,
            12, 24, 36,
            18, 36, 54
        };
        
        {
            Image3 in(Shape3(3), 2);
            Image3 result(in.shape());
            
            integralMultiArray(in, result);        
            shouldEqualSequence(result.begin(), result.end(), desired);
        }
        
        {
            Image4 in(Shape4(3,3,3,2), 2);
            Image4 result(in.shape());
            
            integralMultiArray(in.multiband(), result.multiband());        
            shouldEqualSequence(result.begin(), result.end(), desired);
        }
        
        {
            Image3 in(Shape3(3), 2);
            Image3 result(in.shape());
            Image3 desired_squared(in.shape(), desired);
            desired_squared *= 2;
                        
            integralMultiArraySquared(in, result);        
            shouldEqualSequence(result.begin(), result.end(), desired_squared.begin());
        }
        
        {
            Image4 in(Shape4(3,3,3,2), 2);
            Image4 result(in.shape());
            Image4 desired_squared(in.shape(), desired);
            desired_squared *= 2;
            
            integralMultiArraySquared(in.multiband(), result.multiband());        
            shouldEqualSequence(result.begin(), result.end(), desired_squared.begin());
        }
    }

    void test_4d()
    {
        int desired[] = { 
             2,  4,  6,
             4,  8, 12,
             6, 12, 18,
             
             4,  8, 12,
             8, 16, 24,
            12, 24, 36,
             
             6, 12, 18,
            12, 24, 36,
            18, 36, 54,
            
             4,  8, 12,
             8, 16, 24,
            12, 24, 36,
             
             8, 16, 24,
            16, 32, 48,
            24, 48, 72,
             
            12, 24, 36,
            24, 48, 72,
            36, 72,108
        };
        
        {
            Image4 in(Shape4(3,3,3,2), 2);
            Image4 result(in.shape());
            
            integralMultiArray(in, result);        
            shouldEqualSequence(result.begin(), result.end(), desired);
        }
    }

    void test_vector()
    {        
        int desired[] = {
            2,  4,  6,  8,
            4,  8, 12, 16,
            6, 12, 18, 24,
            
            2,  4,  6,  8,
            4,  8, 12, 16,
            6, 12, 18, 24
        };
        
        {
            Vector2Image2 in(Shape2(4,3), TinyVector<int, 2>(2));
            Vector2Image2 result(in.shape());
            
            integralMultiArray(in, result);
            for(int c=0; c<2; ++c)
            {
                MultiArrayView <2, int> band = result.bindElementChannel(c);
                shouldEqualSequence(band.begin(), band.end(), desired);
            }
        }
    }
};


struct IntegralImageTestSuite
: public test_suite
{
    IntegralImageTestSuite()
    : test_suite("IntegralImageTestSuite")
    {
        add( testCase( &IntegralImageTest::test_2d));
        add( testCase( &IntegralImageTest::test_3d));
        add( testCase( &IntegralImageTest::test_4d));
        add( testCase( &IntegralImageTest::test_vector));
    }
};

int main(int argc, char ** argv)
{
    IntegralImageTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

