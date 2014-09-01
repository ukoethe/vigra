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

#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/convexhull.hxx>


using namespace vigra;
using namespace vigra::detail;

struct ConvexHullTest
{

    MultiArray<2, int> mask;
    
    ConvexHullTest()
    : mask(Shape2(6, 6))
    {
        mask(1, 1) = 1;
        mask(2, 1) = 1;
        mask(2, 2) = 1;
        mask(2, 3) = 1;
        mask(1, 3) = 1;
        mask(3, 1) = 1;
        mask(3, 3) = 1;
        mask(4, 1) = 1;
        mask(4, 3) = 1;
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
        typedef TinyVector<float, 2> Point;

        mask(3,0) = 1; // have the object touch the image border

        TinyVector<int, 2> anchor_point;
        findAnchorPoint(mask, anchor_point);
        vigra::Polygon<Point> contour_points;

        extractContour(mask, anchor_point, contour_points);
        
        Point real_contour_points[] = {
            Point(1, 0.5),
            Point(0.5, 1),
            Point(1, 1.5),
            Point(1.5, 2),
            Point(1, 2.5),
            Point(0.5, 3),
            Point(1, 3.5),
            Point(2, 3.5),
            Point(3, 3.5),
            Point(4, 3.5),
            Point(4.5, 3),
            Point(4, 2.5),
            Point(3, 2.5),
            Point(2.5, 2),
            Point(3, 1.5),
            Point(4, 1.5),
            Point(4.5, 1),
            Point(4, 0.5),
            Point(3.5, 0),
            Point(3, -0.5),
            Point(2.5, 0),
            Point(2, 0.5),
            Point(1, 0.5),
        };

        //Point real_contour_points[] = {
        //    Point(1, 0.5),
        //    Point(0.5, 1),
        //    Point(1, 1.5),
        //    Point(1.5, 2),
        //    Point(1, 2.5),
        //    Point(0.5, 3),
        //    Point(1, 3.5),
        //    Point(2, 3.5),
        //    Point(3, 3.5),
        //    Point(4, 3.5),
        //    Point(4.5, 3),
        //    Point(4, 2.5),
        //    Point(3, 2.5),
        //    Point(2.5, 2),
        //    Point(3, 1.5),
        //    Point(4, 1.5),
        //    Point(4.5, 1),
        //    Point(4, 0.5),
        //    Point(3, 0.5),
        //    Point(2, 0.5),
        //    Point(1, 0.5),
        //};

        shouldEqual(23, contour_points.size());
        //shouldEqual(21, contour_points.size());
        shouldEqualSequence(contour_points.begin(), contour_points.end(), real_contour_points);
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

        //for(auto i = mask.begin(); i != mask.end(); ++i)
        //{
        //    std::cerr << (*i ? "*" : " ");
        //    if(i.point()[0] == mask.shape(0)-1)
        //        std::cerr << "\n";
        //}

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

