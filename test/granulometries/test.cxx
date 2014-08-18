/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2014 by                                 */
/*               Ullrich Koethe,                                        */
/*               Esteban Pardo                                          */
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

#include "vigra/unittest.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/granulometries.hxx"


using namespace vigra;
using namespace vigra::morpho;

struct GranulometiresTest
{

    MultiArray<2, int> mask;
    
    GranulometiresTest()
    {
    
        int size = 6;
        mask = MultiArray<2, int> (Shape2(size, size));

        mask(1, 1) = 1;
        mask(2, 1) = 1;
        mask(3, 1) = 1;
        mask(1, 2) = 1;
        mask(2, 2) = 1;
        mask(3, 2) = 1;
        mask(1, 3) = 1;
        mask(2, 3) = 1;
        mask(3, 3) = 1;
    
    }
    
    void testDiscGranulometries()
    {

        using namespace std;
        
        vector<double> granulometricCurve;
        vector<int> axisAreas;
        
        int realAxisAreas[] = {
        5, 13, 29
        };
        
        double realGranulometricCurve[] = {
        0.44, 0.55, 0
        };
        
        granulometries(mask, 1, 3, granulometricCurve, axisAreas);


        for (unsigned i = 0; i < granulometricCurve.size(); ++i) {
            shouldEqualTolerance(granulometricCurve[i], realGranulometricCurve[i], 1e-1);
        }

        for (unsigned i = 0; i < axisAreas.size(); ++i) {
            shouldEqual(axisAreas[i], realAxisAreas[i]);
        }

    }
    
    void testStrElGranulometries()
    {

        using namespace std;
        
        vector<double> granulometricCurve;
        vector<int> axisAreas;
        
        int realAxisAreas[] = {
        5, 13, 25
        };
        
        double realGranulometricCurve[] = {
        0.44, 0.55, 0
        };
        
        std::vector<Diff2D> strElCoordinates;

        strElCoordinates.push_back(Diff2D(0, 0));
        strElCoordinates.push_back(Diff2D(1, 0));
        strElCoordinates.push_back(Diff2D(-1, 0));
        strElCoordinates.push_back(Diff2D(0, 1));
        strElCoordinates.push_back(Diff2D(0, -1));

        StructuringElement2D crossSE(strElCoordinates);
        
        granulometries(mask, crossSE, 3, granulometricCurve, axisAreas);


        for (unsigned i = 0; i < granulometricCurve.size(); ++i) {
            shouldEqualTolerance(granulometricCurve[i], realGranulometricCurve[i], 1e-1);
        }

        for (unsigned i = 0; i < axisAreas.size(); ++i) {
            shouldEqual(axisAreas[i], realAxisAreas[i]);
        }

    }
    



};

struct GranulometriesTestSuite : public vigra::test_suite
{
    GranulometriesTestSuite()
        : vigra::test_suite("GranulometriesTestSuite")
    {
        add(testCase(&GranulometiresTest::testDiscGranulometries));
        add(testCase(&GranulometiresTest::testStrElGranulometries));
    }
};

int main(int argc, char** argv)
{
    GranulometriesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
