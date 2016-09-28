#define VIGRA_CHECK_BOUNDS

#include <limits>
#include <algorithm>
#include <cmath>
#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/polytope.hxx>


namespace vigra
{

template <class Iterable, class T>
unsigned int count(const Iterable & vec, const T & value)
{
    return std::count(
            vec.begin(),
            vec.end(),
            value);
}

struct FloatStarPolytopeTest
{
    typedef TinyVector<double, 2> Vector2;
    typedef TinyVector<double, 3> Vector3;
    typedef StarPolytope<2, double> Polytope2;
    typedef StarPolytope<3, double> Polytope3;

    FloatStarPolytopeTest()
    : eps_(std::numeric_limits<double>::epsilon() * 3)
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

    void testFill2D()
    {
        Polytope2 poly(
                Vector2(0.  , 0.  ),
                Vector2(1.  , 0.  ),
                Vector2(0.  , 1.  ),
                Vector2(0.25, 0.25));
        MultiArray<2, unsigned int> label_image(vigra::Shape2(5, 5));
        for (auto it = label_image.begin(); it != label_image.end(); it++)
        {
            *it = 0;
        }
        unsigned int ref[25] = {
                0, 0, 0, 0, 0,
                0, 0, 0, 0, 0,
                0, 0, 1, 1, 1,
                0, 0, 1, 1, 0,
                0, 0, 1, 0, 0};
        Vector2 offset(-1., -1.);
        Vector2 scale(0.5, 0.5);
        poly.fill(label_image, 1, offset, scale);
        shouldEqualSequence(label_image.begin(), label_image.end(), ref);
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
        shouldEqual(count(lit_v1, f2), 1);
        shouldEqual(count(lit_v1, f3), 1);
        shouldEqual(count(lit_v2, f1), 1);
        shouldEqual(count(lit_v2, f3), 1);
        shouldEqual(count(lit_v3, f1), 1);
        shouldEqual(count(lit_v3, f2), 1);
        shouldEqual(count(lit_e1, f1), 1);
        shouldEqual(count(lit_e2, f2), 1);
        shouldEqual(count(lit_e3, f3), 1);
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
            shouldEqual(count(aligns1, lemon::INVALID), 2);
        }
        Polytope2::node_type f2 = poly.addFacet(n1, n3);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            shouldEqual(count(aligns1, f2), 1);
            shouldEqual(count(aligns1, lemon::INVALID), 1);
            shouldEqual(count(aligns2, f1), 1);
            shouldEqual(count(aligns2, lemon::INVALID), 1);
        }
        Polytope2::node_type f3 = poly.addFacet(n1, n2);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];
            shouldEqual(count(aligns1, f2), 1);
            shouldEqual(count(aligns1, f3), 1);
            shouldEqual(count(aligns2, f1), 1);
            shouldEqual(count(aligns2, f3), 1);
            shouldEqual(count(aligns3, f2), 1);
            shouldEqual(count(aligns3, f1), 1);
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
            shouldEqual(count(aligns1, lemon::INVALID), 3);
        }
        Polytope3::node_type f2 = poly.addFacet(n1, n3, n4);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            shouldEqual(count(aligns1, f2), 1);
            shouldEqual(count(aligns1, lemon::INVALID), 2);
            shouldEqual(count(aligns2, f1), 1);
            shouldEqual(count(aligns2, lemon::INVALID), 2);
        }
        Polytope3::node_type f3 = poly.addFacet(n1, n2, n4);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];

            shouldEqual(count(aligns1, f2), 1);
            shouldEqual(count(aligns1, f3), 1);
            shouldEqual(count(aligns1, lemon::INVALID), 1);
            shouldEqual(count(aligns2, f1), 1);
            shouldEqual(count(aligns2, f3), 1);
            shouldEqual(count(aligns2, lemon::INVALID), 1);
            shouldEqual(count(aligns3, f1), 1);
            shouldEqual(count(aligns3, f2), 1);
            shouldEqual(count(aligns3, lemon::INVALID), 1);
        }
        Polytope3::node_type f4 = poly.addFacet(n1, n2, n3);
        {
            auto aligns1 = poly.aligns_map_[f1];
            auto aligns2 = poly.aligns_map_[f2];
            auto aligns3 = poly.aligns_map_[f3];
            auto aligns4 = poly.aligns_map_[f4];
            shouldEqual(count(aligns1, f2), 1);
            shouldEqual(count(aligns1, f3), 1);
            shouldEqual(count(aligns1, f4), 1);
            shouldEqual(count(aligns2, f1), 1);
            shouldEqual(count(aligns2, f3), 1);
            shouldEqual(count(aligns2, f4), 1);
            shouldEqual(count(aligns3, f1), 1);
            shouldEqual(count(aligns3, f2), 1);
            shouldEqual(count(aligns3, f4), 1);
            shouldEqual(count(aligns4, f1), 1);
            shouldEqual(count(aligns4, f2), 1);
            shouldEqual(count(aligns4, f3), 1);
        }
    }

    double eps_;
};

struct IntStarPolytopeTest
{
    typedef TinyVector<int, 2> Vector2;
    typedef TinyVector<int, 3> Vector3;
    typedef StarPolytope<2, int> Polytope2;
    typedef StarPolytope<3, int> Polytope3;
    typedef NumericTraits<int>::RealPromote RealPromote;

    IntStarPolytopeTest()
    : eps_(std::numeric_limits<RealPromote>::epsilon() * 2)
    {}

    void testContains2D()
    {
        Polytope2 poly(
                Vector2(0, 0),
                Vector2(4, 0),
                Vector2(0, 2),
                Vector2(1, 1));
        // Internal
        shouldEqual(poly.contains(Vector2( 1,  1)), true);
        // External
        shouldEqual(poly.contains(Vector2( 3,  1)), false);
        shouldEqual(poly.contains(Vector2(-1,  1)), false);
        shouldEqual(poly.contains(Vector2( 2, -1)), false);
        // Edge
        shouldEqual(poly.contains(Vector2( 2,  1)), true);
        shouldEqual(poly.contains(Vector2( 0,  1)), true);
        shouldEqual(poly.contains(Vector2( 2,  0)), true);
        // Vertex
        shouldEqual(poly.contains(Vector2( 0,  0)), true);
        shouldEqual(poly.contains(Vector2( 4,  0)), true);
        shouldEqual(poly.contains(Vector2( 0,  2)), true);
    }

    void testContains3D()
    {
        Polytope3 poly(
                Vector3( 0,  0,  0),
                Vector3( 6,  0,  0),
                Vector3( 0,  6,  0),
                Vector3( 0,  0,  6),
                Vector3( 1,  1,  1));
        // Internal
        shouldEqual(poly.contains(Vector3( 1,  1,  1)), true);
        // External
        shouldEqual(poly.contains(Vector3( 6,  6,  6)), false);
        // Facet
        shouldEqual(poly.contains(Vector3( 2,  2,  2)), true);
        shouldEqual(poly.contains(Vector3( 0,  1,  1)), true);
        shouldEqual(poly.contains(Vector3( 1,  0,  1)), true);
        shouldEqual(poly.contains(Vector3( 1,  1,  0)), true);
        // Edge
        shouldEqual(poly.contains(Vector3( 1,  0,  0)), true);
        shouldEqual(poly.contains(Vector3( 0,  1,  0)), true);
        shouldEqual(poly.contains(Vector3( 0,  0,  1)), true);
        shouldEqual(poly.contains(Vector3( 0,  3,  3)), true);
        shouldEqual(poly.contains(Vector3( 3,  0,  3)), true);
        shouldEqual(poly.contains(Vector3( 3,  3,  0)), true);
        // Vertex
        shouldEqual(poly.contains(Vector3( 0,  0,  0)), true);
        shouldEqual(poly.contains(Vector3( 6,  0,  0)), true);
        shouldEqual(poly.contains(Vector3( 0,  6,  0)), true);
        shouldEqual(poly.contains(Vector3( 0,  0,  6)), true);
    }

    void testNVolume2D()
    {
        {
            Polytope2 poly(
                    Vector2( 0,  0),
                    Vector2( 3,  0),
                    Vector2( 0,  3),
                    Vector2( 1,  1));
            RealPromote n_volume = poly.nVolume();
            shouldEqualTolerance(n_volume, 4.5, eps_);
        }
        {
            Polytope2 poly(
                    Vector2( 1,  1),
                    Vector2( 4,  1),
                    Vector2( 1,  4),
                    Vector2( 2,  2));
            shouldEqualTolerance(poly.nVolume(),  4.5, eps_);
        }
    }

    void testNVolume3D()
    {
        {
            Polytope3 poly(
                    Vector3(0,  0,  0),
                    Vector3(6,  0,  0),
                    Vector3(0,  6,  0),
                    Vector3(0,  0,  6),
                    Vector3(1,  1,  1));
            shouldEqualTolerance(poly.nVolume(), 36., eps_);
        }
        {
            Polytope3 poly(
                    Vector3(1,  1,  1),
                    Vector3(7,  1,  1),
                    Vector3(1,  7,  1),
                    Vector3(1,  1,  7),
                    Vector3(2,  2,  2));
            shouldEqualTolerance(poly.nVolume(), 36., eps_);
        }
    }

    void testNSurface2D()
    {
        Polytope2 poly(
                Vector2(0,  0),
                Vector2(3,  0),
                Vector2(0,  3),
                Vector2(1,  1));
        const RealPromote surf = 2 * 3 + 3 * sqrt(2);
        shouldEqualTolerance(poly.nSurface(), surf, eps_);
    }

    void testNSurface3D()
    {
        Polytope3 poly(
                Vector3(0,  0,  0),
                Vector3(6,  0,  0),
                Vector3(0,  6,  0),
                Vector3(0,  0,  6),
                Vector3(1,  1,  1));
        const RealPromote surf = 6 * 6 * (3. + sqrt(3.)) / 2.;
        shouldEqualTolerance(poly.nSurface(), surf, eps_);
    }

    void testClosed2D()
    {
        Polytope2 poly_closed(
                Vector2(0,  0),
                Vector2(3,  0),
                Vector2(0,  3),
                Vector2(1,  1));
        shouldEqual(poly_closed.closed(), true);

        Polytope2 poly_open(Vector2(1,  1));
        Polytope2::node_type n1 =  poly_open.addVertex(Vector2(0, 0));
        Polytope2::node_type n2 =  poly_open.addVertex(Vector2(3, 0));
        Polytope2::node_type n3 =  poly_open.addVertex(Vector2(0, 3));
        Polytope2::node_type f1 = poly_open.addFacet(n1, n2);
        Polytope2::node_type f2 = poly_open.addFacet(n1, n3);
        shouldEqual(poly_open.closed(), false);
        shouldEqual(poly_open.closed(f1), false);
        shouldEqual(poly_open.closed(f2), false);
    }

    void testClosed3D()
    {
        Polytope3 poly(Vector3(1, 1, 1));
        Polytope3::node_type n1 =  poly.addVertex(Vector3(0, 0, 0));
        Polytope3::node_type n2 =  poly.addVertex(Vector3(6, 0, 0));
        Polytope3::node_type n3 =  poly.addVertex(Vector3(0, 6, 0));
        Polytope3::node_type n4 =  poly.addVertex(Vector3(0, 0, 6));
        Polytope3::node_type f1 = poly.addFacet(n2, n3, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f2 = poly.addFacet(n1, n3, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f3 = poly.addFacet(n1, n2, n4);
        shouldEqual(poly.closed(), false);
        Polytope3::node_type f4 = poly.addFacet(n1, n2, n3);
        shouldEqual(poly.closed(), true);
    }

    // void testLitFacets2D()
    // {
    //     Polytope2 poly(Vector2(0.25, 0.25));
    //     Polytope2::node_type n1 =  poly.addVertex(Vector2(0, 0));
    //     Polytope2::node_type n2 =  poly.addVertex(Vector2(1, 0));
    //     Polytope2::node_type n3 =  poly.addVertex(Vector2(0, 1));
    //     Polytope2::node_type f1 = poly.addFacet(n2, n3);
    //     Polytope2::node_type f2 = poly.addFacet(n1, n3);
    //     Polytope2::node_type f3 = poly.addFacet(n1, n2);
    //     auto lit_v1 = poly.litFacets(Vector2(-1. , -1. ));
    //     auto lit_v2 = poly.litFacets(Vector2( 2. , -0.5));
    //     auto lit_v3 = poly.litFacets(Vector2(-0.5,  2. ));
    //     auto lit_e1 = poly.litFacets(Vector2( 1. ,  1. ));
    //     auto lit_e2 = poly.litFacets(Vector2(-2. ,  0.5));
    //     auto lit_e3 = poly.litFacets(Vector2( 0.5, -2. ));
    //     auto lit0 = poly.litFacets(Vector2( 0.2,  0.2));
    //     shouldEqual(lit0.size(), 0);
    //     shouldEqual(lit_v1.size(), 2);
    //     shouldEqual(lit_v2.size(), 2);
    //     shouldEqual(lit_v3.size(), 2);
    //     shouldEqual(lit_e1.size(), 1);
    //     shouldEqual(lit_e2.size(), 1);
    //     shouldEqual(lit_e3.size(), 1);
    //     shouldEqual(std::count(lit_v1.begin(), lit_v1.end(), f2), 1);
    //     shouldEqual(std::count(lit_v1.begin(), lit_v1.end(), f3), 1);
    //     shouldEqual(std::count(lit_v2.begin(), lit_v2.end(), f1), 1);
    //     shouldEqual(std::count(lit_v2.begin(), lit_v2.end(), f3), 1);
    //     shouldEqual(std::count(lit_v3.begin(), lit_v3.end(), f1), 1);
    //     shouldEqual(std::count(lit_v3.begin(), lit_v3.end(), f2), 1);
    //     shouldEqual(std::count(lit_e1.begin(), lit_e1.end(), f1), 1);
    //     shouldEqual(std::count(lit_e2.begin(), lit_e2.end(), f2), 1);
    //     shouldEqual(std::count(lit_e3.begin(), lit_e3.end(), f3), 1);
    // }

    // void testFindNeighbor2D()
    // {
    //     Polytope2 poly(Vector2(0.25, 0.25));
    //     Polytope2::node_type n1 =  poly.addVertex(Vector2(0, 0));
    //     Polytope2::node_type n2 =  poly.addVertex(Vector2(1, 0));
    //     Polytope2::node_type n3 =  poly.addVertex(Vector2(0, 1));
    //     Polytope2::node_type f1 = poly.addFacet(n2, n3);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         shouldEqual(aligns1.size(), 0);
    //     }
    //     Polytope2::node_type f2 = poly.addFacet(n1, n3);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         auto aligns2 = poly.aligns_map_[f2];
    //         shouldEqual(aligns1.size(), 1);
    //         shouldEqual(aligns1.count(f2), 1);
    //         shouldEqual(aligns2.size(), 1);
    //         shouldEqual(aligns2.count(f1), 1);
    //     }
    //     Polytope2::node_type f3 = poly.addFacet(n1, n2);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         auto aligns2 = poly.aligns_map_[f2];
    //         auto aligns3 = poly.aligns_map_[f3];
    //         shouldEqual(aligns1.size(), 2);
    //         shouldEqual(aligns1.count(f2), 1);
    //         shouldEqual(aligns1.count(f3), 1);
    //         shouldEqual(aligns2.size(), 2);
    //         shouldEqual(aligns2.count(f1), 1);
    //         shouldEqual(aligns2.count(f3), 1);
    //         shouldEqual(aligns3.size(), 2);
    //         shouldEqual(aligns3.count(f2), 1);
    //         shouldEqual(aligns3.count(f1), 1);
    //     }
    // }

    // void testFindNeighbor3D()
    // {
    //     Polytope3 poly(Vector3(0.1, 0.1, 0.1));
    //     Polytope3::node_type n1 =  poly.addVertex(Vector3(0, 0, 0));
    //     Polytope3::node_type n2 =  poly.addVertex(Vector3(1, 0, 0));
    //     Polytope3::node_type n3 =  poly.addVertex(Vector3(0, 1, 0));
    //     Polytope3::node_type n4 =  poly.addVertex(Vector3(0, 0, 1));
    //     Polytope3::node_type f1 = poly.addFacet(n2, n3, n4);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         shouldEqual(aligns1.size(), 0);
    //     }
    //     Polytope3::node_type f2 = poly.addFacet(n1, n3, n4);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         auto aligns2 = poly.aligns_map_[f2];
    //         shouldEqual(aligns1.size(), 1);
    //         shouldEqual(aligns1.count(f2), 1);
    //         shouldEqual(aligns2.size(), 1);
    //         shouldEqual(aligns2.count(f1), 1);
    //     }
    //     Polytope3::node_type f3 = poly.addFacet(n1, n2, n4);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         auto aligns2 = poly.aligns_map_[f2];
    //         auto aligns3 = poly.aligns_map_[f3];

    //         shouldEqual(aligns1.size(), 2);
    //         shouldEqual(aligns1.count(f2), 1);
    //         shouldEqual(aligns1.count(f3), 1);
    //         shouldEqual(aligns2.size(), 2);
    //         shouldEqual(aligns2.count(f1), 1);
    //         shouldEqual(aligns2.count(f3), 1);
    //         shouldEqual(aligns3.size(), 2);
    //         shouldEqual(aligns3.count(f1), 1);
    //         shouldEqual(aligns3.count(f2), 1);
    //     }
    //     Polytope3::node_type f4 = poly.addFacet(n1, n2, n3);
    //     {
    //         auto aligns1 = poly.aligns_map_[f1];
    //         auto aligns2 = poly.aligns_map_[f2];
    //         auto aligns3 = poly.aligns_map_[f3];
    //         auto aligns4 = poly.aligns_map_[f4];
    //         shouldEqual(aligns1.size(), 3);
    //         shouldEqual(aligns1.count(f2), 1);
    //         shouldEqual(aligns1.count(f3), 1);
    //         shouldEqual(aligns1.count(f4), 1);
    //         shouldEqual(aligns2.size(), 3);
    //         shouldEqual(aligns2.count(f1), 1);
    //         shouldEqual(aligns2.count(f3), 1);
    //         shouldEqual(aligns2.count(f4), 1);
    //         shouldEqual(aligns3.size(), 3);
    //         shouldEqual(aligns3.count(f1), 1);
    //         shouldEqual(aligns3.count(f2), 1);
    //         shouldEqual(aligns3.count(f4), 1);
    //         shouldEqual(aligns4.size(), 3);
    //         shouldEqual(aligns4.count(f1), 1);
    //         shouldEqual(aligns4.count(f2), 1);
    //         shouldEqual(aligns4.count(f3), 1);
    //     }
    // }

    RealPromote eps_;
};

struct FloatConvexPolytopeTest
{
    typedef TinyVector<double, 2> Vector2;
    typedef TinyVector<double, 3> Vector3;
    typedef ConvexPolytope<2, double> Polytope2;
    typedef ConvexPolytope<3, double> Polytope3;

    FloatConvexPolytopeTest()
    : eps_(std::numeric_limits<double>::epsilon() * 3)
    {}

    void testAddExtremeVertex2D()
    {
        const int N = 100;
        Polytope2 poly(
                Vector2( 1.,  0.),
                Vector2(-1.,  0.),
                Vector2( 0.,  1.));
        poly.addExtremeVertex(Vector2( 0., -1.));
        shouldEqualTolerance(poly.nVolume(), 2., eps_);
        shouldEqualTolerance(poly.nSurface(), 4. * std::sqrt(2.), eps_);
        for (int n = 0; n < N; n++)
        {
            Vector2 vec(
                    std::cos(2*M_PI*n/N),
                    std::sin(2*M_PI*n/N));
            shouldEqualTolerance(vec.magnitude(), 1., eps_);
            poly.addExtremeVertex(vec);
            shouldEqual(poly.closed(), true);
        }
        const double sur_tgt = 2.*N*std::sin(M_PI / N);
        shouldEqualTolerance(poly.nSurface(), sur_tgt, eps_);
        const double vol_tgt = 1.*N*std::sin(M_PI / N) * std::cos(M_PI / N);
        shouldEqualTolerance(poly.nVolume(), vol_tgt, eps_);
        for (int n = 0; n < 100; n++)
        {
            Vector2 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX));
            if (vec.magnitude() > 1)
            {
                shouldEqual(poly.contains(vec), false);
            }
            else if (vec.magnitude() < std::cos(M_PI / N))
            {
                shouldEqual(poly.contains(vec), true);
            }
        }
    }
    
    void testAddExtremeVertex3D()
    {
        const int N = 100;
        Polytope3 poly(
                Vector3( 1.,  0.,  1.),
                Vector3(-1.,  0.,  1.),
                Vector3( 0.,  1.,  1.),
                Vector3( 0.,  0.,  0.));
        for (int n = 0; n < N; n++)
        {
            Vector3 vec(
                    std::cos(2*M_PI*n/N),
                    std::sin(2*M_PI*n/N),
                    1.);
            poly.addExtremeVertex(vec);
            shouldEqual(poly.closed(), true);
        }
        const double sur_tgt = N * std::sin(M_PI / N) * (
                std::cos(M_PI / N) +  std::sqrt(
                        std::cos(M_PI / N) * std::cos(M_PI / N) + 1.));
        shouldEqualTolerance(poly.nSurface(), sur_tgt, eps_ * N);
        const double vol_tgt = N * std::sin(M_PI / N) * std::cos(M_PI / N) / 3.;
        shouldEqualTolerance(poly.nVolume(), vol_tgt, eps_ * N);
        for (int n = 0; n < 100; n++)
        {
            Vector3 vec(
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    (2*rand() - 1)/static_cast<double>(RAND_MAX),
                    rand()/static_cast<double>(RAND_MAX));
            double dist = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
            if (dist > vec[2])
            {
                shouldEqual(poly.contains(vec), false);
            }
            else if (dist < (vec[2] * std::cos(M_PI / N)))
            {
                shouldEqual(poly.contains(vec), true);
            }
        }
    }

    void testAddNonExtremeVertex2D()
    {
        const int N = 1000;
        const double eps = 4. / sqrt(N);
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
            if (vec.magnitude() <= 1.)
            {
                poly.addExtremeVertex(vec);
                shouldEqual(poly.closed(), true);
            }
        }
        const double sur_err = (2*M_PI - poly.nSurface()) / (2.*M_PI);
        shouldEqual(sur_err < eps, true);
        shouldEqual(sur_err > 0, true);
        const double vol_err = (M_PI - poly.nVolume()) / (M_PI);
        shouldEqualTolerance(vol_err, 0, eps);
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

    void testAddNonExtremeVertex3D()
    {
        const int N = 1000;
        const double eps = 9. / sqrt(N);
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
            if (vec.magnitude() <= 1.)
            {
                poly.addExtremeVertex(vec);
                shouldEqual(poly.closed(), true);
            }
        }
        const double sur_err = (4.*M_PI - poly.nSurface()) / (4.*M_PI);
        shouldEqualTolerance(sur_err, 0, eps);
        shouldEqual(sur_err > 0, true);
        const double vol_err = (4./3.*M_PI - poly.nVolume()) / (4./3.*M_PI);
        shouldEqualTolerance(vol_err, 0, eps);
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
        add(testCase(&FloatStarPolytopeTest::testClosed2D));
        add(testCase(&FloatStarPolytopeTest::testClosed3D));
        add(testCase(&FloatStarPolytopeTest::testContains2D));
        add(testCase(&FloatStarPolytopeTest::testContains3D));
        add(testCase(&FloatStarPolytopeTest::testFill2D));
        add(testCase(&FloatStarPolytopeTest::testFindNeighbor2D));
        add(testCase(&FloatStarPolytopeTest::testFindNeighbor3D));
        add(testCase(&FloatStarPolytopeTest::testLitFacets2D));
        add(testCase(&FloatStarPolytopeTest::testNSurface2D));
        add(testCase(&FloatStarPolytopeTest::testNSurface3D));
        add(testCase(&FloatStarPolytopeTest::testNVolume2D));
        add(testCase(&FloatStarPolytopeTest::testNVolume3D));
        add(testCase(&FloatConvexPolytopeTest::testAddExtremeVertex2D));
        add(testCase(&FloatConvexPolytopeTest::testAddExtremeVertex3D));
        /*
        add(testCase(&FloatConvexPolytopeTest::testAddNonExtremeVertex2D));
        add(testCase(&FloatConvexPolytopeTest::testAddNonExtremeVertex3D));
        */
        add(testCase(&IntStarPolytopeTest::testClosed2D));
        add(testCase(&IntStarPolytopeTest::testClosed3D));
        add(testCase(&IntStarPolytopeTest::testContains2D));
        add(testCase(&IntStarPolytopeTest::testContains3D));
        add(testCase(&IntStarPolytopeTest::testNSurface2D));
        add(testCase(&IntStarPolytopeTest::testNSurface3D));
        add(testCase(&IntStarPolytopeTest::testNVolume2D));
        add(testCase(&IntStarPolytopeTest::testNVolume3D));
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

