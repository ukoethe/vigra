#include <iostream>
#include "unittest.hxx"
#include "vigra/polynomial.hxx"
#include "vigra/array_vector.hxx"

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
        
        should(polynomialRoots(p, roots));        
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
    }
};

int main()
{
    MathTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
