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
#include <vigra/polygon.hxx>
#include <vigra/convexhull.hxx>
#include "convex_hull_test.hxx"


using namespace vigra;
using namespace vigra::detail;

struct PolygonTest
{

    MultiArray<2, int> mask;
    
    PolygonTest()
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

    void print(MultiArray<2, int> const & a)
    {
        for(auto i = a.begin(); i != a.end(); ++i)
        {
            std::cerr << (char)*i;
            if(i.point()[0] == a.shape(0)-1)
                std::cerr << "\n";
        }
    }

    void testFillPolygon()
    {
        typedef TinyVector<double, 2> P;
        typedef vigra::Polygon<P> Poly;

        {
            // triangle in general position
            P p[] = { P(0.5, 0.5), P(1.5, 4.1), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+4);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
            for(auto i=mask.begin(); i != mask.end(); ++i)
                if(*i == ' ')
                    should(!poly.contains(i.point()));
                else
                    should(poly.contains(i.point()));
        }
        {
            // triangle, lower knot on scanline
            P p[] = { P(0.5, 0.5), P(1.5, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+4);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // triangle, lower knot on pixel
            P p[] = { P(0.5, 0.5), P(1.0, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+4);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " *        "
                " *        "
                " *        "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // triangle, all knots on scanlines
            P p[] = { P(0.5, 1.0), P(1.5, 4.0), P(2.5, 1.0), P(0.5, 1.0) };
            Poly poly(p, p+4);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // triangle, all knots on scanlines
            P p[] = { P(0.5, 1.0), P(1.5, 4.0), P(2.5, 1.0), P(2.2, 1.0), P(0.5, 1.0) };
            Poly poly(p, p+5);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // quadrilateral in general position
            P p[] = { P(0.5, 0.5), P(0.9, 2.5), P(1.5, 4.1), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+5);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
            for(auto i=mask.begin(); i != mask.end(); ++i)
                if(*i == ' ')
                    should(!poly.contains(i.point()));
                else
                    should(poly.contains(i.point()));
        }
        {
            // quadrilateral, one knot on scanline
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.5, 4.1), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+5);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // quadrilateral, two knots on scanline
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.5, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+5);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // quadrilateral, one knot on scanline, one knot on pixel
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.0, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+5);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " *        "
                " *        "
                " *        "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon in general position
            P p[] = { P(0.5, 0.5), P(0.9, 2.5), P(1.2, 2.5), P(1.5, 4.1), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
            for(auto i=mask.begin(); i != mask.end(); ++i)
                if(*i == ' ')
                    should(!poly.contains(i.point()));
                else
                    should(poly.contains(i.point()));
        }
        {
            // concave pentagon, two knots on scanline
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.2, 2.0), P(1.5, 4.1), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon, three knots on scanline
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.2, 2.0), P(1.5, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "          "
                "          "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon, two knots on scanline, one knot on pixel
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.2, 2.0), P(2.0, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "  *       "
                "  *       "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon, three knots on scanline, one knot on pixel
            P p[] = { P(0.5, 0.5), P(0.9, 2.0), P(1.1, 2.0), P(1.3, 2.0), P(2.0, 4.0), P(2.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " **       "
                " **       "
                "  *       "
                "  *       "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon in general position
            P p[] = { P(0.5, 0.5), P(1.0, 4.1), P(2.0, 2.9), P(3.0, 4.1), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
            for(auto i=mask.begin(); i != mask.end(); ++i)
                if(*i == ' ')
                    should(!poly.contains(i.point()));
                else
                    should(poly.contains(i.point()));
        }
        {
            // concave pentagon, one knot on scanline
            P p[] = { P(0.5, 0.5), P(1.0, 4.1), P(1.9, 2.0), P(3.0, 4.1), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon, one knot on pixel
            P p[] = { P(0.5, 0.5), P(1.0, 4.1), P(2.0, 2.0), P(3.0, 4.1), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon, three knots on pixel
            P p[] = { P(0.5, 0.5), P(1.0, 4.0), P(2.0, 2.0), P(3.0, 4.0), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave pentagon, all knots on pixel
            P p[] = { P(1.0, 1.0), P(1.0, 4.0), P(2.0, 2.0), P(3.0, 4.0), P(3.0, 1.0), P(1.0, 1.0) };
            Poly poly(p, p+6);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon in general position
            P p[] = { P(0.5, 0.5), P(1.0, 4.1), P(1.9, 2.9), P(2.1, 2.9), P(3.0, 4.1), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
            for(auto i=mask.begin(); i != mask.end(); ++i)
                if(*i == ' ')
                    should(!poly.contains(i.point()));
                else
                    should(poly.contains(i.point()));
        }
        {
            // concave hexagon, two knots on scanline
            P p[] = { P(0.5, 0.5), P(1.0, 4.1), P(1.9, 2.0), P(2.1, 2.0), P(3.0, 4.1), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon, four knots on scanline
            P p[] = { P(0.5, 0.5), P(1.0, 4.0), P(1.9, 2.0), P(2.1, 2.0), P(3.0, 4.0), P(3.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ***      "
                " ***      "
                " * *      "
                " * *      "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon, four knots on scanline
            P p[] = { P(0.5, 0.5), P(1.0, 4.0), P(1.9, 2.0), P(3.1, 2.0), P(4.0, 4.0), P(4.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ****     "
                " ****     "
                " *  *     "
                " *  *     "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon, four knots on scanline
            P p[] = { P(0.5, 0.5), P(1.0, 4.0), P(1.9, 2.0), P(3.1, 2.0), P(4.0, 2.0), P(4.5, 0.5), P(0.5, 0.5) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ****     "
                " ****     "
                " *        "
                " *        "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave hexagon, all knots on scanline
            P p[] = { P(0.5, 1.0), P(1.0, 4.0), P(1.9, 2.0), P(3.1, 2.0), P(4.0, 2.0), P(4.5, 1.0), P(0.5, 1.0) };
            Poly poly(p, p+7);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ****     "
                " ****     "
                " *        "
                " *        "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
        }
        {
            // concave heptagon, all knots on scanline
            P p[] = { P(0.5, 1.0), P(1.0, 4.0), P(1.9, 2.0), P(3.1, 2.0), P(4.0, 2.0), P(4.5, 1.0), P(3.2, 1.0), P(0.5, 1.0) };
            Poly poly(p, p+8);
            MultiArray<2, int> mask(Shape2(10,6), ' ');
            fillPolygon(poly, mask, int('*'));

            std::string ref =
                "          "
                " ****     "
                " ****     "
                " *        "
                " *        "
                "          "
            ;
            shouldEqualSequence(mask.begin(), mask.end(), ref.begin());
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

        shouldEqual(23, contour_points.size());
        shouldEqualSequence(contour_points.begin(), contour_points.end(), real_contour_points);

        MultiArray<2, int> render(mask.shape());
        fillPolygon(contour_points, render, 1);

        shouldEqualSequence(mask.begin(), mask.end(), render.begin());
    }
    
   void testConvexHull()
    {
        typedef TinyVector<double, 2> Point;
        ArrayVector<Point> points, reference, hull;
        points.push_back(Point(0.0, 0.0));
        points.push_back(Point(2.0, 1.0));
        points.push_back(Point(2.0, -1.0));
        points.push_back(Point(0.0, 2.0));
        points.push_back(Point(-2.0, 1.0));
        points.push_back(Point(-2.0, -1.0));
        points.push_back(Point(0.0, -2.0));

        reference.push_back(Point(0.0, -2.0));
        reference.push_back(Point(-2.0, -1.0));
        reference.push_back(Point(-2.0, 1.0));
        reference.push_back(Point(0.0, 2.0));
        reference.push_back(Point(2.0, 1.0));
        reference.push_back(Point(2.0, -1.0));
        reference.push_back(Point(0.0, -2.0));
        
        vigra::convexHull(points, hull);

        shouldEqual(7u, hull.size());
        shouldEqualSequence(reference.begin(), reference.end(), hull.begin());

        hull.clear();

        vigra::convexHull(reference, hull);
        
        shouldEqual(7u, hull.size());
        shouldEqualSequence(reference.begin(), reference.end(), hull.begin());
        
        typedef Point P;
        P p[200] = { P(0.0, 0.0), P(42.0, 468.0), P(335.0, 501.0), P(170.0, 725.0), P(479.0, 359.0), 
                  P(963.0, 465.0), P(706.0, 146.0), P(282.0, 828.0), P(962.0, 492.0), 
                  P(996.0, 943.0), P(828.0, 437.0), P(392.0, 605.0), P(903.0, 154.0), 
                  P(293.0, 383.0), P(422.0, 717.0), P(719.0, 896.0), P(448.0, 727.0), 
                  P(772.0, 539.0), P(870.0, 913.0), P(668.0, 300.0), P(36.0, 895.0), 
                  P(704.0, 812.0), P(323.0, 334.0), P(674.0, 665.0), P(142.0, 712.0), 
                  P(254.0, 869.0), P(548.0, 645.0), P(663.0, 758.0), P(38.0, 860.0), 
                  P(724.0, 742.0), P(530.0, 779.0), P(317.0, 36.0), P(191.0, 843.0), 
                  P(289.0, 107.0), P(41.0, 943.0), P(265.0, 649.0), P(447.0, 806.0), 
                  P(891.0, 730.0), P(371.0, 351.0), P(7.0, 102.0), P(394.0, 549.0), 
                  P(630.0, 624.0), P(85.0, 955.0), P(757.0, 841.0), P(967.0, 377.0), 
                  P(932.0, 309.0), P(945.0, 440.0), P(627.0, 324.0), P(538.0, 539.0), 
                  P(119.0, 83.0), P(930.0, 542.0), P(834.0, 116.0), P(640.0, 659.0), 
                  P(705.0, 931.0), P(978.0, 307.0), P(674.0, 387.0), P(22.0, 746.0), 
                  P(925.0, 73.0), P(271.0, 830.0), P(778.0, 574.0), P(98.0, 513.0), 
                  P(987.0, 291.0), P(162.0, 637.0), P(356.0, 768.0), P(656.0, 575.0), 
                  P(32.0, 53.0), P(351.0, 151.0), P(942.0, 725.0), P(967.0, 431.0), 
                  P(108.0, 192.0), P(8.0, 338.0), P(458.0, 288.0), P(754.0, 384.0), 
                  P(946.0, 910.0), P(210.0, 759.0), P(222.0, 589.0), P(423.0, 947.0), 
                  P(507.0, 31.0), P(414.0, 169.0), P(901.0, 592.0), P(763.0, 656.0), 
                  P(411.0, 360.0), P(625.0, 538.0), P(549.0, 484.0), P(596.0, 42.0), 
                  P(603.0, 351.0), P(292.0, 837.0), P(375.0, 21.0), P(597.0, 22.0), 
                  P(349.0, 200.0), P(669.0, 485.0), P(282.0, 735.0), P(54.0, 1000.0), 
                  P(419.0, 939.0), P(901.0, 789.0), P(128.0, 468.0), P(729.0, 894.0), 
                  P(649.0, 484.0), P(808.0, 422.0), P(311.0, 618.0), P(814.0, 515.0), 
                  P(310.0, 617.0), P(936.0, 452.0), P(601.0, 250.0), P(520.0, 557.0), 
                  P(799.0, 304.0), P(225.0, 9.0), P(845.0, 610.0), P(990.0, 703.0), 
                  P(196.0, 486.0), P(94.0, 344.0), P(524.0, 588.0), P(315.0, 504.0), 
                  P(449.0, 201.0), P(459.0, 619.0), P(581.0, 797.0), P(799.0, 282.0), 
                  P(590.0, 799.0), P(10.0, 158.0), P(473.0, 623.0), P(539.0, 293.0), 
                  P(39.0, 180.0), P(191.0, 658.0), P(959.0, 192.0), P(816.0, 889.0), 
                  P(157.0, 512.0), P(203.0, 635.0), P(273.0, 56.0), P(329.0, 647.0), 
                  P(363.0, 887.0), P(876.0, 434.0), P(870.0, 143.0), P(845.0, 417.0), 
                  P(882.0, 999.0), P(323.0, 652.0), P(22.0, 700.0), P(558.0, 477.0), 
                  P(893.0, 390.0), P(76.0, 713.0), P(601.0, 511.0), P(4.0, 870.0), 
                  P(862.0, 689.0), P(402.0, 790.0), P(256.0, 424.0), P(3.0, 586.0), 
                  P(183.0, 286.0), P(89.0, 427.0), P(618.0, 758.0), P(833.0, 933.0), 
                  P(170.0, 155.0), P(722.0, 190.0), P(977.0, 330.0), P(369.0, 693.0), 
                  P(426.0, 556.0), P(435.0, 550.0), P(442.0, 513.0), P(146.0, 61.0), 
                  P(719.0, 754.0), P(140.0, 424.0), P(280.0, 997.0), P(688.0, 530.0), 
                  P(550.0, 438.0), P(867.0, 950.0), P(194.0, 196.0), P(298.0, 417.0), 
                  P(287.0, 106.0), P(489.0, 283.0), P(456.0, 735.0), P(115.0, 702.0), 
                  P(317.0, 672.0), P(787.0, 264.0), P(314.0, 356.0), P(186.0, 54.0), 
                  P(913.0, 809.0), P(833.0, 946.0), P(314.0, 757.0), P(322.0, 559.0), 
                  P(647.0, 983.0), P(482.0, 145.0), P(197.0, 223.0), P(130.0, 162.0), 
                  P(536.0, 451.0), P(174.0, 467.0), P(45.0, 660.0), P(293.0, 440.0), 
                  P(254.0, 25.0), P(155.0, 511.0), P(746.0, 650.0), P(187.0, 314.0), 
                  P(475.0, 23.0), P(169.0, 19.0), P(788.0, 906.0), P(959.0, 392.0), 
                  P(203.0, 626.0), P(478.0, 415.0), P(315.0, 825.0), P(335.0, 875.0), 
                  P(373.0, 160.0), P(834.0, 71.0), P(488.0, 298.0) };
                  
        P ref[10] = { P(0.0, 0.0), 
                      P(4.0, 870.0),
                      P(54.0, 1000.0),
                      P(882.0, 999.0),
                      P(996.0, 943.0),
                      P(987.0, 291.0),
                      P(959.0, 192.0), 
                      P(925.0, 73.0),
                      P(597.0, 22.0), 
                      P(0.0, 0.0) };

        points = vigra::ArrayVector<Point>(p, p+200);
        hull.clear();
        
        vigra::convexHull(points, hull);
        
        shouldEqual(10u, hull.size());
        shouldEqualSequence(ref, ref+10, hull.begin());

        int size = sizeof(convexHullInputs) / sizeof(Point);
        points = vigra::ArrayVector<Point>(convexHullInputs, convexHullInputs+size);
        hull.clear();
        
        vigra::convexHull(points, hull);
        
        shouldEqual(17u, hull.size());
        shouldEqualSequence(convexHullReference, convexHullReference+17, hull.begin());
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

struct PolygonTestSuite : public vigra::test_suite
{
    PolygonTestSuite()
        : vigra::test_suite("PolygonTestSuite")
    {
        add(testCase(&PolygonTest::testFillPolygon));
        add(testCase(&PolygonTest::testExtractContour));
        add(testCase(&PolygonTest::testConvexHull));
        add(testCase(&PolygonTest::testConvexHullFeatures));
    }
};

int main(int argc, char** argv)
{
    PolygonTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cerr << test.report() << std::endl;

    return failed != 0;
}

