/************************************************************************/
/*                                                                      */
/*             Copyright 2011-2012 by Ullrich Koethe                    */
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
#include <sstream>
#include <map>
#include <set>

#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/accumulator.hxx>


using namespace vigra;

// mask cl.exe shortcomings
#if defined(_MSC_VER)
#pragma warning( disable : 4503 )
#endif


struct LeEnd{};

struct GraphFeatureTest
{
    GraphFeatureTest()
    {}
    
    void testStandardizeTag()
    {
        using namespace vigra::acc;

        typedef double DataType;
        typedef typename MultiArrayShape<3>::type CoordType;
        typedef typename MultiArrayShape<3>::type   Point;


        typedef Select< 
            DataArg<1>,
            Variance, 
            Mean, 
            StdDev, 
            Minimum, 
            Maximum, 
            RootMeanSquares, 
            Skewness,
            Covariance,
            RegionCenter
        >  SelectType;



        typedef StandAloneAccumulatorChain<3, DataType, SelectType> FreeChain;
        FreeChain a;

        a.updatePassN(1.0, Point(3,0,0), 1);
        a.updatePassN(2.0, Point(0,3,0), 1);
        a.updatePassN(3.0, Point(0,0,3), 1);
        


        a.updatePassN(1.0, Point(3,0,0), 2);
        a.updatePassN(2.0, Point(0,3,0), 2);
        a.updatePassN(3.0, Point(0,0,3), 2);

        shouldEqualTolerance( get<Mean>(a), 2.0, 0.000001);

        CoordType rCenter = get<RegionCenter>(a);
        CoordType trueCenter(1,1,1);
        shouldEqualSequence(rCenter.begin(),rCenter.end(), trueCenter.begin());

    }

};


struct GraphFeaturesTestSuite : public vigra::test_suite
{
    GraphFeaturesTestSuite()
        : vigra::test_suite("GraphFeaturesTestSuite")
    {
        add(testCase(&GraphFeatureTest::testStandardizeTag));
    }
};

int main(int argc, char** argv)
{
    GraphFeaturesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
