

#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/convexhull.hxx>


using namespace vigra;
using namespace vigra::detail;

struct ConvexHullTest
{

    MultiArray<2, int> mask;
    
    ConvexHullTest()
    {
    
        int size = 6;
        mask = MultiArray<2, int> (Shape2(size, size));

        mask(1, 1) = 1;
        mask(1, 2) = 1;
        mask(2, 2) = 1;
        mask(3, 2) = 1;
        mask(3, 1) = 1;
        mask(1, 3) = 1;
        mask(3, 3) = 1;
        mask(1, 4) = 1;
        mask(3, 4) = 1;
    
    }
    
    void testPushLinePoints()
    {
        using namespace std;

        std::vector<TinyVector<float, 2> > middle_points;
        TinyVector<float, 2> a(0, 0);
        TinyVector<float, 2> b(2, 0.5);
        pushLinePoints(a, b, middle_points);
        
        TinyVector<float, 2> real_middle_points[] = {
        TinyVector<float, 2>(1, 0.25)
        };
        
        for (int i = 0; i < 1; ++i) {
            shouldEqualTolerance(middle_points[i][0], real_middle_points[i][0], 1e-2);
            shouldEqualTolerance(middle_points[i][1], real_middle_points[i][1], 1e-2);
        }

    }
    
    void testExtractContour()
    {
        using namespace std;

        TinyVector<int, 2> anchor_point;
        findAnchorPoint(mask, anchor_point);
        vector<TinyVector<float, 2> > contour_points;

        extractContour(mask, anchor_point, contour_points);
        
        TinyVector<float, 2> real_contour_points[] = {
        TinyVector<float, 2>(1, 0.5),
        TinyVector<float, 2>(0.5, 1),
        TinyVector<float, 2>(0.5, 2),
        TinyVector<float, 2>(0.5, 3),
        TinyVector<float, 2>(0.5, 4),
        TinyVector<float, 2>(1, 4.5),
        TinyVector<float, 2>(1.5, 4),
        TinyVector<float, 2>(1.5, 3),
        TinyVector<float, 2>(2, 2.5),
        TinyVector<float, 2>(2.5, 3),
        TinyVector<float, 2>(2.5, 4),
        TinyVector<float, 2>(3, 4.5),
        TinyVector<float, 2>(3.5, 4),
        TinyVector<float, 2>(3.5, 3),
        TinyVector<float, 2>(3.5, 2),
        TinyVector<float, 2>(3.5, 1),
        TinyVector<float, 2>(3, 0.5),
        TinyVector<float, 2>(2.5, 1),
        TinyVector<float, 2>(2, 1.5),
        TinyVector<float, 2>(1.5, 1),
        };
        
        for (int i = 0; i < 20; ++i) {
            shouldEqualTolerance(contour_points[i][0], real_contour_points[i][0], 1e-2);
            shouldEqualTolerance(contour_points[i][1], real_contour_points[i][1], 1e-2);
        }

    }
    
    void testConvexHullFeatures()
    {
        using namespace vigra::acc;
        using namespace std;

        int size = 6;
        MultiArray<2, int> mask(vigra::Shape2(size, size));

        mask(1, 1) = 1;
        mask(2, 1) = 1;
        mask(2, 2) = 1;
        mask(2, 3) = 1;
        mask(1, 3) = 1;
        mask(3, 1) = 1;
        mask(3, 3) = 1;
        mask(4, 1) = 1;
        mask(4, 3) = 1;

        ConvexHullFeatures<MultiArray<2, int> > chf(mask);

        shouldEqual(chf.getInputArea(), 9);

        shouldEqual(chf.getConvexHullArea(), 12);


        shouldEqualTolerance(chf.getConvexity(), 0.75, 1e-2);

        shouldEqual(chf.getConvexityDefectCount(), 2);

        shouldEqualTolerance(chf.getConvexityDefectAreaMean(), 1.5, 1e-1);

        shouldEqualTolerance(chf.getConvexityDefectAreaVariance(), 0.25, 1e-2);

        shouldEqualTolerance(chf.getConvexityDefectAreaSkewness(), 0.0, 1e-1);

        shouldEqualTolerance(chf.getConvexityDefectAreaKurtosis(), -2.0, 1e-1);

        shouldEqual(chf.getInputPerimeter(), 20);

        shouldEqual(chf.getConvexHullPerimeter(), 14);

        shouldEqualTolerance(chf.getRugosity(), 1.428571, 1e-6);
    }

};

struct ConvexHullTestSuite : public vigra::test_suite
{
    ConvexHullTestSuite()
        : vigra::test_suite("ConvexHullTestSuite")
    {
        add(testCase(&ConvexHullTest::testConvexHullFeatures));
        
        add(testCase(&ConvexHullTest::testExtractContour));
        
        add(testCase(&ConvexHullTest::testPushLinePoints));
    }
};

int main(int argc, char** argv)
{
    ConvexHullTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
