

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

        structuringElement2D crossSE(strElCoordinates);
        
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
