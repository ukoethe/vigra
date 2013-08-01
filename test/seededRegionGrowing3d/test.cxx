/************************************************************************/
/*                                                                      */
/*       Copyright 2007 by F. Heinrich, B. Seppke, Ullrich Koethe       */
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
#include "vigra/unittest.hxx"

#include "vigra/seededregiongrowing3d.hxx"

using namespace vigra;

struct SeededRegionGrowing3DTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DoubleVolume;

    SeededRegionGrowing3DTest() :
            vol1(IntVolume::difference_type(5,5,5)),
            vol2(DoubleVolume::difference_type(4,4,4)),
            vol3(IntVolume::difference_type(5,5,5)),
            distvol1(DoubleVolume::difference_type(5,5,5)),
            distvol2(DoubleVolume::difference_type(4,4,4))
    {
        static const int in1[] = { 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0,  
                                   0, 0, 1, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,

                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,

                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,

                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,

                                   0, 0, 0, 0, 0,  
                                   0, 0, 0, 0, 0,  
                                   0, 0, 2, 0, 0,  
                                   0, 0, 0, 0, 0,   
                                   0, 0, 0, 0, 0};

        IntVolume::iterator i = vol1.begin();
        IntVolume::iterator end = vol1.end();
        const int * p = in1;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

        for(int z=0; z<5; ++z)
            for(int y=0; y<5; ++y)
                for(int x=0; x<5; ++x){
                    distvol1(x,y,z)=std::min( ((x-2.0)*(x-2.0)+(y-2.0)*(y-2.0)+(z-0.0)*(z-0.0)), 
                                              ((x-2.0)*(x-2.0)+(y-2.0)*(y-2.0)+(z-4.0)*(z-4.0)) );
                }

        static const double in2[] = { 0.0, 0.0, 0.0, 0.0, 
                                      0.0, 1.0, 0.0, 0.0, 
                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,

                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,

                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0,

                                      0.0, 0.0, 0.0, 0.0, 
                                      0.0, 0.0, 0.0, 0.0, 
                                      0.0, 0.0, 2.0, 0.0,  
                                      0.0, 0.0, 0.0, 0.0};

        DoubleVolume::iterator id = vol2.begin();
        DoubleVolume::iterator endd = vol2.end();
        const double * pd = in2;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

        for(int z=0; z<4; ++z)
            for(int y=0; y<4; ++y)
                for(int x=0; x<4; ++x){
                    distvol2(x,y,z)=std::min( ((x-1.0)*(x-1.0)+(y-1.0)*(y-1.0)+(z-0.0)*(z-0.0)), 
                                              ((x-2.0)*(x-2.0)+(y-2.0)*(y-2.0)+(z-3.0)*(z-3.0)) );
                }

        static const int in3[] = { 1, 1, 1, 1, 1, 
                                   1, 1, 1, 1, 1,  
                                   1, 1, 1, 1, 1,  
                                   1, 1, 1, 1, 1,  
                                   1, 1, 1, 1, 1,

                                   1, 1, 1, 1, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 1, 1, 1, 1,

                                   1, 1, 1, 1, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 2, 3, 2, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 1, 1, 1, 1,

                                   1, 1, 1, 1, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 2, 2, 2, 1,  
                                   1, 1, 1, 1, 1,

                                   4, 4, 4, 4, 4,  
                                   4, 4, 4, 4, 4,  
                                   4, 4, 4, 4, 4,  
                                   4, 4, 4, 4, 4,   
                                   4, 4, 4, 4, 4};

        i = vol3.begin();
        end = vol3.end();
        p = in3;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

    }

    struct DirectCostFunctor
    {
        typedef double argument_type;
        typedef double result_type;
        typedef double cost_type;

        void operator()(double const &) {}

        double const & cost(double const & v) const
        {
            return v;
        }
    };

    void voronoiTest()
    {
        DoubleVolume res(vol2);

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing3D(distvol2, vol2, res, cost, CompleteGrow);

        DoubleVolume::iterator i = res.begin();

        int x,y,z;

        for(z=0; z<4; ++z)
        {
            for(y=0; y<4; ++y)
            {
                for(x=0; x<4; ++x)
                {
                    double dist = *i++;
                    double dist1 = VIGRA_CSTD::sqrt((1.0 - x)*(1.0 - x) +
                                                    (1.0 - y)*(1.0 - y) +
                                                    (0.0 - z)*(0.0 - z)  );
                                                    
                    double dist2 = VIGRA_CSTD::sqrt((2.0 - x)*(2.0 - x) +
                                                    (2.0 - y)*(2.0 - y) +
                                                    (3.0 - z)*(3.0 - z)  );
                                                    
                    double desired = (dist1 <= dist2) ? 1 : 2;

                    if(VIGRA_CSTD::fabs(dist1 - dist2) > 1e-10)
                        shouldEqual(dist, desired);
                }
            }
        }
    }

    void voronoiTestWithBorder()
    {

        static const int desired[] = {  1, 1, 1, 1, 1, 
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,

                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,  
                                        1, 1, 1, 1, 1,

                                        0, 0, 0, 0, 0,  
                                        0, 0, 0, 0, 0,  
                                        0, 0, 0, 0, 0,  
                                        0, 0, 0, 0, 0,  
                                        0, 0, 0, 0, 0,

                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,

                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,  
                                        2, 2, 2, 2, 2,   
                                        2, 2, 2, 2, 2};

        IntVolume res(vol1);

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing3D(srcMultiArrayRange(distvol1), srcMultiArray(vol1),
                              destMultiArray(res), cost, KeepContours);

        //int c=1;
        //for(IntVolume::iterator iter=res.begin(); iter!=res.end(); ++iter, ++c){
        //    std::cerr << *iter << ", ";
        //    if(c%5==0) std::cerr << std::endl;
        //    if(c%25==0) std::cerr << std::endl;
        //}

        shouldEqualSequence(res.begin(), res.end(), desired);

    }

    void simpleTest()
    {
        IntVolume res(vol3);

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(4);
        seededRegionGrowing3D(srcMultiArrayRange(vol3), srcMultiArray(vol3),
                              destMultiArray(res), cost, CompleteGrow);

        shouldEqualSequence(res.begin(), res.end(), vol3.begin());
    }
    
    IntVolume    vol1;
    DoubleVolume vol2;
    IntVolume    vol3;
    DoubleVolume distvol1;
    DoubleVolume distvol2;
};



struct SeededRegionGrowing3DTestSuite
: public vigra::test_suite
{
    SeededRegionGrowing3DTestSuite()
    : vigra::test_suite("SeededRegionGrowing3DTestSuite")
    {
        add( testCase( &SeededRegionGrowing3DTest::voronoiTest));
        add( testCase( &SeededRegionGrowing3DTest::voronoiTestWithBorder));
        add( testCase( &SeededRegionGrowing3DTest::simpleTest));
    }
};

int main(int argc, char ** argv)
{
    SeededRegionGrowing3DTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

