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
    typedef MultiArray<3, int>                  Image3;
    typedef MultiArray<4, int>                  Image4;
    typedef MultiArrayView<2, int>              View2;
    typedef MultiArrayView<3, int>              View3;
    typedef MultiArrayView<3, Multiband<int> >  ChannelView2;
    typedef MultiArrayView<4, Multiband<int> >  ChannelView3;

    Image2  img2;
    Image3  channel_img2;
    Image3  img3;
    Image4  channel_img3;
    
    IntegralImageTest()
    :   img2(Shape2(3), 2),
        channel_img2(Shape3(3), 2),
        img3(Shape3(3), 2),
        channel_img3(Shape4(3), 2)
    {
    }

    void test_2d()
    {
        int reference_data[] = { 
             2,  4,  6,
             4,  8, 12,
             6, 12, 18
        };
        
        Image2 result(img2.shape());
        
        integralImage(img2, result);        
        shouldEqualSequence(result.begin(), result.end(), reference_data);
        
        result = 0;
        integralImageNew(img2, result);        
        shouldEqualSequence(result.begin(), result.end(), reference_data);
        
        Image3 channel_result(channel_img2.shape());        
        integralImage(ChannelView2(channel_img2), ChannelView2(channel_result));
        
        for(int c=0; c<channel_result.shape(2); ++c)
            shouldEqualSequence(channel_result.bindOuter(c).begin(), channel_result.bindOuter(c).end(), 
                                reference_data);
        
        int reference_data_squared[] = { 
             4,  8, 12,
             8, 16, 24,
            12, 24, 36
        };
        
        result = 0;
        integralImage2(img2, result);        
        shouldEqualSequence(result.begin(), result.end(), reference_data_squared);
        
        channel_result = 0;
        integralImage2(ChannelView2(channel_img2), ChannelView2(channel_result));
        for(int c=0; c<channel_result.shape(2); ++c)
            shouldEqualSequence(channel_result.bindOuter(c).begin(), channel_result.bindOuter(c).end(), 
                                reference_data_squared);
    }

    void test_3d()
    {
        int reference_data[] = { 
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
        
        Image3 result(img3.shape());
        
        integralVolume(img3, result);        
        shouldEqualSequence(result.begin(), result.end(), reference_data);
        
        Image4 channel_result(channel_img3.shape());        
        integralVolume(ChannelView3(channel_img3), ChannelView3(channel_result));
        
        for(int c=0; c<channel_result.shape(2); ++c)
        {
            shouldEqualSequence(channel_result.bindOuter(c).begin(), channel_result.bindOuter(c).end(), 
                                reference_data);
        }
        
        int reference_data_squared[] = { 
             4,  8, 12,
             8, 16, 24,
            12, 24, 36,
             
             8, 16, 24,
            16, 32, 48,
            24, 48, 72,
             
            12, 24, 36,
            24, 48, 72,
            36, 72, 108
        };
        
        result = 0;
        integralVolume2(img3, result);        
        shouldEqualSequence(result.begin(), result.end(), reference_data_squared);
        
        channel_result = 0;
        integralVolume2(ChannelView3(channel_img3), ChannelView3(channel_result));
        
        for(int c=0; c<channel_result.shape(2); ++c)
        {
            shouldEqualSequence(channel_result.bindOuter(c).begin(), channel_result.bindOuter(c).end(), 
                                reference_data_squared);
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
    }
};

int main(int argc, char ** argv)
{
    IntegralImageTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

