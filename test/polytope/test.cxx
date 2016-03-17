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

#define VIGRA_CHECK_BOUNDS

#include <limits>
#include <algorithm>
#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/polytope.hxx>


namespace vigra
{

struct StarPolytopeTest
{
    typedef TinyVector<double, 2> Vector2;
    typedef TinyVector<double, 3> Vector3;
    typedef StarPolytope<2, double> Polytope2;
    typedef StarPolytope<3, double> Polytope3;

    StarPolytopeTest()
    : eps_(std::numeric_limits<double>::epsilon() * 2)
    {}

    void testContains2D()
    {
        Polytope2 poly(
                Vector2(0.  , 0.  ),
                Vector2(1.  , 0.  ),
                Vector2(0.  , 1.  ),
                Vector2(0.25, 0.25));
        // Internal
        shouldEqual(poly.contains(Vector2(0.25, 0.25)), true);
        // External
        shouldEqual(poly.contains(Vector2(1.  , 1.  )), false);
        // Edge
        shouldEqual(poly.contains(Vector2(0.5 , 0.  )), true);
        shouldEqual(poly.contains(Vector2(0.  , 0.5 )), true);
        shouldEqual(poly.contains(Vector2(0.5 , 0.5 )), true);
        // Vertex
        shouldEqual(poly.contains(Vector2(0.  , 0.  )), true);
        shouldEqual(poly.contains(Vector2(1.  , 0.  )), true);
        shouldEqual(poly.contains(Vector2(0.  , 1.  )), true);
    }

    void testContains3D()
    {
        Polytope3 poly(
                Vector3(0.  , 0.  , 0.  ),
                Vector3(1.  , 0.  , 0.  ),
                Vector3(0.  , 1.  , 0.  ),
                Vector3(0.  , 0.  , 1.  ),
                Vector3(0.25, 0.25, 0.25));
        // Internal
        shouldEqual(poly.contains(Vector3(0.25, 0.25, 0.25)), true);
        // External
        shouldEqual(poly.contains(Vector3(1.  , 1.  , 1.  )), false);
        // Facet
        shouldEqual(poly.contains(Vector3(0.  , 0.2 , 0.2)), true);
        shouldEqual(poly.contains(Vector3(0.2 , 0.  , 0.2)), true);
        shouldEqual(poly.contains(Vector3(0.2 , 0.2 , 0. )), true);
        shouldEqual(poly.contains(Vector3(0.25, 0.25, 0.5)), true);
        // Edge
        shouldEqual(poly.contains(Vector3(0.5 , 0.  , 0. )), true);
        shouldEqual(poly.contains(Vector3(0.  , 0.5 , 0. )), true);
        shouldEqual(poly.contains(Vector3(0.  , 0.  , 0.5)), true);
        shouldEqual(poly.contains(Vector3(0.  , 0.5 , 0.5)), true);
        shouldEqual(poly.contains(Vector3(0.5 , 0.  , 0.5)), true);
        shouldEqual(poly.contains(Vector3(0.5 , 0.5 , 0.0)), true);
        // Vertex
        shouldEqual(poly.contains(Vector3(0.  , 0.  , 0. )), true);
        shouldEqual(poly.contains(Vector3(1.  , 0.  , 0. )), true);
        shouldEqual(poly.contains(Vector3(0.  , 1.  , 0. )), true);
        shouldEqual(poly.contains(Vector3(0.  , 0.  , 1. )), true);
    }

    void testNVolume2D()
    {
        {
            Polytope2 poly(
                    Vector2(0.  , 0.  ),
                    Vector2(1.  , 0.  ),
                    Vector2(0.  , 1.  ),
                    Vector2(0.25, 0.25));
            shouldEqual(abs(poly.nVolume() - .5) < eps_, true);
        }
        {
            Polytope2 poly(
                    Vector2(0.5 , 0.5 ),
                    Vector2(0.5 , 1.  ),
                    Vector2(1.  , 0.5 ),
                    Vector2(0.6 , 0.6 ));
            shouldEqual(abs(poly.nVolume() - .125) < eps_, true);
        }
    }

    void testNVolume3D()
    {
        {
            Polytope3 poly(
                    Vector3(0.  , 0.  , 0.  ),
                    Vector3(1.  , 0.  , 0.  ),
                    Vector3(0.  , 1.  , 0.  ),
                    Vector3(0.  , 0.  , 1.  ),
                    Vector3(0.25, 0.25, 0.25));
            shouldEqual(abs(poly.nVolume() - 1./6.) < eps_, true);
        }
        {
            Polytope3 poly(
                    Vector3(0.5 , 0.5 , 0.5 ),
                    Vector3(1.  , 0.5 , 0.5 ),
                    Vector3(0.5 , 1.  , 0.5 ),
                    Vector3(0.5 , 0.5 , 1.  ),
                    Vector3(0.6 , 0.6 , 0.6 ));
            shouldEqual(abs(poly.nVolume() - 1./(6.*8.)) < eps_, true);
        }
    }

    void testNSurface2D()
    {
        Polytope2 poly(
                Vector2(0.  , 0.  ),
                Vector2(1.  , 0.  ),
                Vector2(0.  , 1.  ),
                Vector2(0.25, 0.25));
        shouldEqual(abs(poly.nSurface() - (2. + sqrt(2.))) < eps_, true);
    }

    void testNSurface3D()
    {
        Polytope3 poly(
                Vector3(0.  , 0.  , 0.  ),
                Vector3(1.  , 0.  , 0.  ),
                Vector3(0.  , 1.  , 0.  ),
                Vector3(0.  , 0.  , 1.  ),
                Vector3(0.25, 0.25, 0.25));
        const double surf = (3. + sqrt(3.)) / 2.;
        shouldEqual(abs(poly.nSurface() - surf) < eps_, true);
    }

    void testClosed2D()
    {
        Polytope2 poly_closed(
                Vector2(0.  , 0.  ),
                Vector2(1.  , 0.  ),
                Vector2(0.  , 1.  ),
                Vector2(0.25, 0.25));
        shouldEqual(poly_closed.closed(), true);

        Polytope2 poly_open(Vector2(0.25, 0.25));
        Polytope2::node_type n1 =  poly_open.addVertex(Vector2(0, 0));
        Polytope2::node_type n2 =  poly_open.addVertex(Vector2(1, 0));
        Polytope2::node_type n3 =  poly_open.addVertex(Vector2(0, 1));
        Polytope2::node_type f1 = poly_open.addFacet(n1, n2);
        Polytope2::node_type f2 = poly_open.addFacet(n1, n3);
        shouldEqual(poly_open.closed(), false);
        shouldEqual(poly_open.closed(f1), false);
        shouldEqual(poly_open.closed(f2), false);
    }

    void testClosed3D()
    {
        Polytope3 poly(Vector3(0.1, 0.1, 0.1));
        Polytope3::node_type n1 =  poly.addVertex(Vector3(0, 0, 0));
        Polytope3::node_type n2 =  poly.addVertex(Vector3(1, 0, 0));
        Polytope3::node_type n3 =  poly.addVertex(Vector3(0, 1, 0));
        Polytope3::node_type n4 =  poly.addVertex(Vector3(0, 0, 1));
        Polytope3::node_type f1 = poly.addFacet(n2, n3, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f2 = poly.addFacet(n1, n3, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f3 = poly.addFacet(n1, n2, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f4 = poly.addFacet(n1, n2, n3);
        shouldEqual(poly.closed(), true);
    }

    void testLitFacets2D()
    {
        Polytope2 poly(Vector2(0.25, 0.25));
        Polytope2::node_type n1 =  poly.addVertex(Vector2(0, 0));
        Polytope2::node_type n2 =  poly.addVertex(Vector2(1, 0));
        Polytope2::node_type n3 =  poly.addVertex(Vector2(0, 1));
        Polytope2::node_type f1 = poly.addFacet(n2, n3);
        Polytope2::node_type f2 = poly.addFacet(n1, n3);
        Polytope2::node_type f3 = poly.addFacet(n1, n2);
        auto lit_v1 = poly.litFacets(Vector2(-1. , -1. ));
        auto lit_v2 = poly.litFacets(Vector2( 2. , -0.5));
        auto lit_v3 = poly.litFacets(Vector2(-0.5,  2. ));
        auto lit_e1 = poly.litFacets(Vector2( 1. ,  1. ));
        auto lit_e2 = poly.litFacets(Vector2(-2. ,  0.5));
        auto lit_e3 = poly.litFacets(Vector2( 0.5, -2. ));
        auto lit0 = poly.litFacets(Vector2( 0.2,  0.2));
        shouldEqual(lit0.size(), 0);
        shouldEqual(lit_v1.size(), 2);
        shouldEqual(lit_v2.size(), 2);
        shouldEqual(lit_v3.size(), 2);
        shouldEqual(lit_e1.size(), 1);
        shouldEqual(lit_e2.size(), 1);
        shouldEqual(lit_e3.size(), 1);
        shouldEqual(std::count(lit_v1.begin(), lit_v1.end(), f2), 1);
        shouldEqual(std::count(lit_v1.begin(), lit_v1.end(), f3), 1);
        shouldEqual(std::count(lit_v2.begin(), lit_v2.end(), f1), 1);
        shouldEqual(std::count(lit_v2.begin(), lit_v2.end(), f3), 1);
        shouldEqual(std::count(lit_v3.begin(), lit_v3.end(), f1), 1);
        shouldEqual(std::count(lit_v3.begin(), lit_v3.end(), f2), 1);
        shouldEqual(std::count(lit_e1.begin(), lit_e1.end(), f1), 1);
        shouldEqual(std::count(lit_e2.begin(), lit_e2.end(), f2), 1);
        shouldEqual(std::count(lit_e3.begin(), lit_e3.end(), f3), 1);
    }

    void testFindNeighbor2D()
    {
        Polytope2 poly(Vector2(0.25, 0.25));
        Polytope2::node_type n1 =  poly.addVertex(Vector2(0, 0));
        Polytope2::node_type n2 =  poly.addVertex(Vector2(1, 0));
        Polytope2::node_type n3 =  poly.addVertex(Vector2(0, 1));
        Polytope2::node_type f1 = poly.addFacet(n2, n3);
        {
            auto aligns1 = poly.aligns_map_[f1];
            shouldEqual(aligns1.size(), 0);
        }
        Polytope2::node_type f2 = poly.addFacet(n1, n3);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            shouldEqual(aligns1.size(), 1);
            shouldEqual(aligns1.count(f2), 1);
            shouldEqual(aligns2.size(), 1);
            shouldEqual(aligns2.count(f1), 1);
        }
        Polytope2::node_type f3 = poly.addFacet(n1, n2);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];
            shouldEqual(aligns1.size(), 2);
            shouldEqual(aligns1.count(f2), 1);
            shouldEqual(aligns1.count(f3), 1);
            shouldEqual(aligns2.size(), 2);
            shouldEqual(aligns2.count(f1), 1);
            shouldEqual(aligns2.count(f3), 1);
            shouldEqual(aligns3.size(), 2);
            shouldEqual(aligns3.count(f2), 1);
            shouldEqual(aligns3.count(f1), 1);
        }
    }

    void testFindNeighbor3D()
    {
        Polytope3 poly(Vector3(0.1, 0.1, 0.1));
        Polytope3::node_type n1 =  poly.addVertex(Vector3(0, 0, 0));
        Polytope3::node_type n2 =  poly.addVertex(Vector3(1, 0, 0));
        Polytope3::node_type n3 =  poly.addVertex(Vector3(0, 1, 0));
        Polytope3::node_type n4 =  poly.addVertex(Vector3(0, 0, 1));
        Polytope3::node_type f1 = poly.addFacet(n2, n3, n4);
        {
            auto aligns1 = poly.aligns_map_[f1];
            shouldEqual(aligns1.size(), 0);
        }
        Polytope3::node_type f2 = poly.addFacet(n1, n3, n4);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            shouldEqual(aligns1.size(), 1);
            shouldEqual(aligns1.count(f2), 1);
            shouldEqual(aligns2.size(), 1);
            shouldEqual(aligns2.count(f1), 1);
        }
        Polytope3::node_type f3 = poly.addFacet(n1, n2, n4);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];

            shouldEqual(aligns1.size(), 2);
            shouldEqual(aligns1.count(f2), 1);
            shouldEqual(aligns1.count(f3), 1);
            shouldEqual(aligns2.size(), 2);
            shouldEqual(aligns2.count(f1), 1);
            shouldEqual(aligns2.count(f3), 1);
            shouldEqual(aligns3.size(), 2);
            shouldEqual(aligns3.count(f1), 1);
            shouldEqual(aligns3.count(f2), 1);
        }
        Polytope3::node_type f4 = poly.addFacet(n1, n2, n3);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];
            auto aligns4 = poly.aligns_map_[f4];
            shouldEqual(aligns1.size(), 3);
            shouldEqual(aligns1.count(f2), 1);
            shouldEqual(aligns1.count(f3), 1);
            shouldEqual(aligns1.count(f4), 1);
            shouldEqual(aligns2.size(), 3);
            shouldEqual(aligns2.count(f1), 1);
            shouldEqual(aligns2.count(f3), 1);
            shouldEqual(aligns2.count(f4), 1);
            shouldEqual(aligns3.size(), 3);
            shouldEqual(aligns3.count(f1), 1);
            shouldEqual(aligns3.count(f2), 1);
            shouldEqual(aligns3.count(f4), 1);
            shouldEqual(aligns4.size(), 3);
            shouldEqual(aligns4.count(f1), 1);
            shouldEqual(aligns4.count(f2), 1);
            shouldEqual(aligns4.count(f3), 1);
        }
    }

    double eps_;
};

struct ConvexPolytopeTest
{
    typedef TinyVector<double, 2> Vector2;
    typedef TinyVector<double, 3> Vector3;
    typedef ConvexPolytope<2, double> Polytope2;
    typedef ConvexPolytope<3, double> Polytope3;

    ConvexPolytopeTest()
    : eps_(std::numeric_limits<double>::epsilon() * 2)
    {}

    void testAddExtremeVertex2D()
    {
        const int N = 100;
        const double eps = 1. / sqrt(N);
        Polytope2 poly(
                Vector2( 1.,  0.),
                Vector2(-1.,  0.),
                Vector2( 0.,  1.));
        poly.addExtremeVertex(Vector2(0., -1.));
        shouldEqual(abs(poly.nVolume() - 2.) < eps_, true);
        shouldEqual(abs(poly.nSurface() - 4. * sqrt(2.)) < eps, true);
        for (int n = 0; n < N; n++)
        {
            Vector2 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX));
            vec /= vec.magnitude();
            shouldEqual(poly.contains(vec), false);
            shouldEqual(poly.litFacets(vec).size(), 1);
            poly.addExtremeVertex(vec);
            shouldEqual(poly.closed(), true);
        }
        const double sur_err = (2*M_PI - poly.nSurface()) / (2.*M_PI);
        shouldEqual(sur_err < eps, true);
        shouldEqual(sur_err > 0, true);
        const double vol_err = (M_PI - poly.nVolume()) / (M_PI);
        shouldEqual(vol_err < eps, true);
        shouldEqual(vol_err > 0, true);
        for (int n = 0; n < 100; n++)
        {
            Vector2 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX));
            if (abs(vec.magnitude() - 1) > eps)
            {
                shouldEqual(poly.contains(vec), vec.magnitude() < 1.);
            }
        }
    }

    void testAddExtremeVertex3D()
    {
        const int N = 100;
        const double eps = 3. / sqrt(N);
        Polytope3 poly(
                Vector3( 1.,  0.,  0.),
                Vector3(-1.,  0.,  0.),
                Vector3( 0.,  1.,  0.),
                Vector3( 0.,  0.,  1.));
        for (int n = 0; n < N; n++)
        {
            Vector3 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX));
            vec /= vec.magnitude();
            shouldEqual(poly.contains(vec), false);
            poly.addExtremeVertex(vec);
            shouldEqual(poly.closed(), true);
        }
        const double sur_err = (4.*M_PI - poly.nSurface()) / (4.*M_PI);
        shouldEqual(sur_err < eps, true);
        shouldEqual(sur_err > 0, true);
        const double vol_err = (4./3.*M_PI - poly.nVolume()) / (4./3.*M_PI);
        shouldEqual(vol_err < eps, true);
        shouldEqual(vol_err > 0, true);
        for (int n = 0; n < 100; n++)
        {
            Vector3 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX));
            if (abs(vec.magnitude() - 1) > eps)
            {
                shouldEqual(poly.contains(vec), vec.magnitude() < 1.);
            }
        }
    }

    double eps_;
};

struct PolytopeTestSuite : public vigra::test_suite
{
    PolytopeTestSuite() : vigra::test_suite("PolytopeTestSuite")
    {
        add(testCase(&StarPolytopeTest::testClosed2D));
        add(testCase(&StarPolytopeTest::testClosed3D));
        add(testCase(&StarPolytopeTest::testContains2D));
        add(testCase(&StarPolytopeTest::testContains3D));
        add(testCase(&StarPolytopeTest::testFindNeighbor2D));
        add(testCase(&StarPolytopeTest::testFindNeighbor3D));
        add(testCase(&StarPolytopeTest::testLitFacets2D));
        add(testCase(&StarPolytopeTest::testNSurface2D));
        add(testCase(&StarPolytopeTest::testNSurface3D));
        add(testCase(&StarPolytopeTest::testNVolume2D));
        add(testCase(&StarPolytopeTest::testNVolume3D));
        add(testCase(&ConvexPolytopeTest::testAddExtremeVertex2D));
        add(testCase(&ConvexPolytopeTest::testAddExtremeVertex3D));
    }
};

} // namespace vigra

int main(int argc, char** argv)
{
    vigra::PolytopeTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cerr << test.report() << std::endl;

    return failed != 0;
}

