#include <iostream>
#include <algorithm>
#include "unittest.hxx"
#include "vigra/polynomial.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/splines.hxx"
#include "vigra/gaussians.hxx"
#include "vigra/rational.hxx"

static double coefficients[][12] = 
{
    { 5.0, -416.0, 720.0, -464.0, 136.0, -18.0, 1.0 },
    { 8.0, 40320.0, -109584.0, 118124.0, -67284.0, 22449.0, -4536.0, 546.0, -36.0, 1.0},
    { 3.0, 1e10, -1e10, -1e-10, 1e-10},
    { 3.0, 1e-10, -1e-10, -1e10, 1e10},
    { 10.0, 2.88e-8, -1.848e-6, 0.00005204, -0.0008458, 0.008777,
           -0.06072, 0.2835, -0.882, 1.75, -2.0, 1.0},
    { 5.0, 0.3411268890719874, 0.48265610836623374, 0.29941395284477745, 
           0.13065520631476124, 0.68342489290545338, 0.0017437185812028133 },
    { 3.0, -1.0, 1000001000001.0 / 1e6, -1000001000001.0 / 1e6, 1.0},
    { 8.0, 36.0, 0.0, 85.0, 0.0, 63.0, 0.0, 15.0, 0.0, 1.0 }
};

typedef std::complex<double> C;

static C reference[][12] = 
{
    { C(1e-12), C(2.0), C(2.0), C(2.0), C(6.0, -4.0), C(6.0, 4.0) },
    { C(1e-11), C(6.0), C(5.0), C(4.0), C(7.0), C(3.0), C(8.0), C(2.0), C(1.0)  },
    { C(1e-12), C(1.0), C(-1e10), C(1e10) },
    { C(1e-12), C(1.0), C(-1e-10), C(1e-10) },
    { C(1e-5), C(0.2), C(0.2), C(0.2), C(0.3), C(0.3), 
               C(0.1), C(0.1), C(0.1), C(0.4), C(0.1) },
    { C(1e-12), C(0.47331479192572767, 0.89542786425410759), C(0.47331479192572767, -0.89542786425410759), 
                C(-0.56839260551055271, 0.4046562986541693), C(-391.74516023901123), 
                C(-0.56839260551055271, -0.4046562986541693) },
    { C(1e-12), C(1e6), C(1.0), C(1e-6) },
    { C(1e-12), C(0.0, 1.0), C(0.0, 1.0), C(0.0, -1.0), C(0.0, -2.0), 
                C(0.0, 2.0), C(0.0, -1.0), C(0.0, 3.0), C(0.0, -3.0) }
};

#if 0
#undef should
#define should(v) (std::cerr << #v << ": " << (v) << std::endl)
#undef shouldEqual
#define shouldEqual(v1, v2) \
(std::cerr << #v1 << " == " << #v2 << ": " << (v1) << " " << (v2) << std::endl)
#undef shouldEqualTolerance
#define shouldEqualTolerance(v1, v2, e) \
(std::cerr << #v1 << " == " << #v2 << ": " << (v1) << " " << (v2) << std::endl)
#endif

template <unsigned int N, class POLYNOMIAL>
struct PolynomialTest
{
    void testPolynomial()
    {
        double epsilon = reference[N][0].real();
        unsigned int order = (unsigned int)(coefficients[N][0] + 0.5);
        POLYNOMIAL p(coefficients[N]+1, order);
        
        vigra::ArrayVector<std::complex<double> > roots;
        
        should(polynomialRoots(p, roots));
        shouldEqual(roots.size(), order);
        for(unsigned int i = 0; i<roots.size(); ++i)
        {
            shouldEqualTolerance(roots[i].real(), reference[N][i+1].real(), epsilon);
            shouldEqualTolerance(roots[i].imag(), reference[N][i+1].imag(), epsilon);
        }
    }
};

struct HighOrderPolynomialTest
{
    void testPolynomial()
    {
        unsigned int order = 80;
        double epsilon = 1e-12;
        vigra::ArrayVector<double> coeffs(order+1, 0.0);
        coeffs[0] = -1.0;
        coeffs[order] = 1.0;
        vigra::Polynomial<double> p(coeffs.begin(), order);
        
        vigra::ArrayVector<std::complex<double> > roots;
        
        should(vigra::polynomialRoots(p, roots));        
        shouldEqual(roots.size(), order);
        for(unsigned int i = 0; i<roots.size(); ++i)
        {
            shouldEqualTolerance(std::abs(roots[i]), 1.0, epsilon);
            C r = p(roots[i]);
            shouldEqualTolerance(r.real(), 0.0, epsilon);
            shouldEqualTolerance(r.imag(), 0.0, epsilon);
        }
        vigra::ArrayVector<double> rroots;
        should(polynomialRealRoots(p, rroots));
        shouldEqual(rroots.size(), 2);
        shouldEqualTolerance(rroots[0], 1.0, epsilon);
        shouldEqualTolerance(rroots[1], -1.0, epsilon);
    }
};

template <int ORDER>
struct SplineTest
{
    typedef vigra::BSpline<ORDER, double> BS;
    typedef vigra::BSplineBase<ORDER, double> BSB;
    BS spline;
    BSB splineBase;
    
    void testValues()
    {
        double r = spline.radius();
        shouldEqual(r, splineBase.radius());
        
        for(int d = 0; d <= ORDER+1; ++d)
        {
            for(double x = -r-0.5; x <= r+0.5; x += 0.5)
                shouldEqualTolerance(spline(x, d), splineBase(x, d), 1e-15);
        }
    }
    
    void testPrefilterCoefficients()
    {
        int n = ORDER / 2;
        double const * ps = BS::prefilterCoefficients();
        double const * psb = BSB::prefilterCoefficients();
        
        if(n == 0)
        {
            shouldEqual(ps[0], 0.0);
            shouldEqual(psb[0], 0.0);
        }
        else
        {
            double * psb1 = const_cast<double *>(psb);
            std::sort(psb1, psb1 + n);
            for(int i = 0; i < n; ++i)
                shouldEqualTolerance(ps[i], psb[i], 1e-14);
        }
    }
    
    void testWeightMatrix()
    {
        int n = ORDER + 1;
        typename BS::WeightMatrix & ws = BS::weights();
        typename BSB::WeightMatrix & wsb = BSB::weights();
        
        for(int d = 0; d < n; ++d)
            for(int i = 0; i < n; ++i)
                shouldEqualTolerance(ws[d][i], wsb[d][i], 1e-14);
    }
};

struct FunctionsTest
{
    void testGaussians()
    {
        vigra::Gaussian<double> g,
                          g1(2.0, 1),
                          g2(1.0, 2),
                          g3(2.0, 3),
                          g4(2.0, 4),
                          g5(2.0, 5);
        
        double epsilon = 1e-15;
        shouldEqual(g.derivativeOrder(), 0);
        shouldEqual(g.sigma(), 1.0);
        shouldEqualTolerance(g(0.0), 0.3989422804014327, epsilon);
        shouldEqualTolerance(g(0.5), 0.35206532676429952, epsilon);
        shouldEqualTolerance(g(1.0), 0.24197072451914337, epsilon);
        shouldEqualTolerance(g(-1.0), 0.24197072451914337, epsilon);
        
        shouldEqual(g1.derivativeOrder(), 1);
        shouldEqual(g1.sigma(), 2.0);
        shouldEqualTolerance(g1(0.0), 0, epsilon);
        shouldEqualTolerance(g1(0.5), -0.024166757300178077, epsilon);
        shouldEqualTolerance(g1(1.0), -0.044008165845537441, epsilon);
        shouldEqualTolerance(g1(-1.0), 0.044008165845537441, epsilon);
        
        shouldEqual(g2.derivativeOrder(), 2);
        shouldEqual(g2.sigma(), 1.0);
        shouldEqualTolerance(g2(0.0), -0.3989422804014327, epsilon);
        shouldEqualTolerance(g2(0.5), -0.26404899507322466, epsilon);
        shouldEqualTolerance(g2(1.0), 0, epsilon);
        shouldEqualTolerance(g2(-1.0), 0, epsilon);
        shouldEqualTolerance(g2(1.5), 0.16189699458236467, epsilon);
        shouldEqualTolerance(g2(-1.5), 0.16189699458236467, epsilon);

        shouldEqual(g3.derivativeOrder(), 3);
        shouldEqual(g3.sigma(), 2.0);
        shouldEqualTolerance(g3(0.0), 0, epsilon);
        shouldEqualTolerance(g3(0.5), 0.017747462392318277, epsilon);
        shouldEqualTolerance(g3(1.0), 0.030255614018806987, epsilon);
        shouldEqualTolerance(g3(-1.0), -0.030255614018806987, epsilon);
        shouldEqualTolerance(g3(2.0*VIGRA_CSTD::sqrt(3.0)), 0, epsilon);
        shouldEqualTolerance(g3(-2.0*VIGRA_CSTD::sqrt(3.0)), 0, epsilon);
    
        shouldEqualTolerance(g4(0.0), 0.037400838787634318, epsilon);
        shouldEqualTolerance(g4(1.0), 0.017190689783413062, epsilon);
        shouldEqualTolerance(g4(-1.0), 0.017190689783413062, epsilon);
        shouldEqualTolerance(g4(1.483927568605452), 0, epsilon);
        shouldEqualTolerance(g4(4.668828436677955), 0, epsilon);
        shouldEqualTolerance(g5(0.0), 0, epsilon);
        shouldEqualTolerance(g5(1.0), -0.034553286464660257, epsilon);
        shouldEqualTolerance(g5(-1.0), 0.034553286464660257, epsilon);
        shouldEqualTolerance(g5(2.711252359948531), 0, epsilon);
        shouldEqualTolerance(g5(5.713940027745611), 0, epsilon);
    }
};

struct RationalTest
{
    void testGcdLcm()
    {
        shouldEqual(vigra::gcd(24, 18), 6);
        shouldEqual(vigra::lcm(6, 4), 12);
    }
    
    void testOperators()
    {
        typedef vigra::Rational<int> R;
        shouldEqual(R(3,4) + R(12,6), R(11,4));
        shouldEqual(R(3,4) - R(12,6), R(-5,4));
        shouldEqual(R(3,4) * R(12,6), R(3,2));
        shouldEqual(R(3,4) / R(12,6), R(3,8));
        shouldEqual(abs(R(-3,4)), R(3,4));
        should(R(3,4) == R(9,12));
        should(R(3,4) != R(12,6));
        should(R(3,4) < R(12,6));
        should(R(19,4) > R(12,6));
        should(R(3,4) <= R(12,6));
        should(R(19,4) >= R(12,6));

        shouldEqual(R(3,4) + 2, R(11,4));
        shouldEqual(R(3,4) - 2, R(-5,4));
        shouldEqual(R(3,4) * 2, R(3,2));
        shouldEqual(R(3,4) / 2, R(3,8));
        should(!(R(3,4) == 2));
        should(R(3,4) != 2);
        should(R(3,4) < 2);
        should(R(19,4) > 2);
        should(R(3,4) <= 2);
        should(R(19,4) >= 2);

        shouldEqual(2 + R(3,4), R(11,4));
        shouldEqual(2 - R(3,4), R(5,4));
        shouldEqual(2 * R(3,4), R(3,2));
        shouldEqual(2 / R(3,4), R(8, 3));
        should(!(2 == R(3,4)));
        should(2 != R(3,4));
        should(2 > R(3,4));
        should(2 < R(19,4));
        should(2 >= R(3,4));
        should(2 <= R(19,4));
    }
};

struct MathTestSuite
: public vigra::test_suite
{
    MathTestSuite()
    : vigra::test_suite("MathTest")
    {
        typedef vigra::Polynomial<double> P1;
        typedef vigra::StaticPolynomial<10, double> P2;
        add( testCase((&PolynomialTest<0, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<1, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<2, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<3, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<4, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<5, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<6, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<7, P1>::testPolynomial)));
        add( testCase((&PolynomialTest<0, P2>::testPolynomial)));
        add( testCase((&PolynomialTest<1, P2>::testPolynomial)));
        add( testCase((&PolynomialTest<2, P2>::testPolynomial)));
        add( testCase(&HighOrderPolynomialTest::testPolynomial));
        add( testCase(&SplineTest<0>::testValues));
        add( testCase(&SplineTest<1>::testValues));
        add( testCase(&SplineTest<2>::testValues));
        add( testCase(&SplineTest<3>::testValues));
        add( testCase(&SplineTest<5>::testValues));
        add( testCase(&SplineTest<0>::testPrefilterCoefficients));
        add( testCase(&SplineTest<1>::testPrefilterCoefficients));
        add( testCase(&SplineTest<2>::testPrefilterCoefficients));
        add( testCase(&SplineTest<3>::testPrefilterCoefficients));
        add( testCase(&SplineTest<5>::testPrefilterCoefficients));
        add( testCase(&SplineTest<0>::testWeightMatrix));
        add( testCase(&SplineTest<1>::testWeightMatrix));
        add( testCase(&SplineTest<2>::testWeightMatrix));
        add( testCase(&SplineTest<3>::testWeightMatrix));
        add( testCase(&SplineTest<5>::testWeightMatrix));
        add( testCase(&FunctionsTest::testGaussians));
        add( testCase(&RationalTest::testGcdLcm));
        add( testCase(&RationalTest::testOperators));
    }
};

int main()
{
    MathTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
