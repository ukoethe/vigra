#include <iostream>
#include "unittest.hxx"
#include "vigra/polynomial.hxx"

static double coefficients[][12] = 
{
    { 5.0, -416.0, 720.0, -464.0, 136.0, -18.0, 1.0 },
    { 8.0, 40320.0, -109584.0, 118124.0, -67284.0, 22449.0, -4536.0, 546.0, -36.0, 1.0},
    { 3.0, 1e10, -1e10, -1e-10, 1e-10},
    { 3.0, 1e-10, -1e-10, -1e10, 1e10},
    { 10.0, 2.88e-8, -1.848e-6, 0.00005204, -0.0008458, 0.008777,
           -0.06072, 0.2835, -0.882, 1.75, -2.0, 1.0},
    { 5.0, 0.3411268890719874, 0.48265610836623374, 0.29941395284477745, 
           0.13065520631476124, 0.68342489290545338, 0.0017437185812028133 }
};

typedef std::complex<double> C;

static C reference[][12] = 
{
    { C(1e-12), C(2.0), C(2.0), C(2.0), C(6.0, 4.0), C(6.0, -4.0) },
    { C(1e-6), C(6.0), C(5.0), C(4.0), C(3.0), C(7.0), C(8.0), C(2.0), C(1.0) },
    { C(1e-12), C(1.0), C(-1e10), C(1e10) },
    { C(1e-12), C(1.0), C(-1e-10), C(1e-10) },
    { C(1e-5), C(0.2), C(0.2), C(0.2), C(0.3), C(0.3), 
               C(0.1), C(0.1), C(0.1), C(0.4), C(0.1) },
    { C(1e-6), C(0.47331479192572767, 0.89542786425410759), C(0.47331479192572767, -0.89542786425410759), 
               C(-0.56839260551055271, 0.4046562986541693), C(-391.74516023901123), 
               C(-0.56839260551055271, -0.4046562986541693) }
};

template <unsigned int N>
struct MathTest
{
    void testPolynomial()
    {
        double epsilon = reference[N][0].real();
        unsigned int order = (unsigned int)(coefficients[N][0] + 0.5);
        vigra::Polynomial<double> p(coefficients[N]+1, coefficients[N]+2+order);
        std::vector<std::complex<double> > roots;
        polynomialRoots(p, roots);
        should(roots.size() == order);
        for(unsigned int i = 0; i<roots.size(); ++i)
        {
            shouldEqualTolerance(roots[i].real(), reference[N][i+1].real(), epsilon);
            shouldEqualTolerance(roots[i].imag(), reference[N][i+1].imag(), epsilon);
        }
    }
};


struct MathTestSuite
: public vigra::test_suite
{
    MathTestSuite()
    : vigra::test_suite("MathTest")
    {
        add( testCase(&MathTest<0>::testPolynomial));
        add( testCase(&MathTest<1>::testPolynomial));
        add( testCase(&MathTest<2>::testPolynomial));
        add( testCase(&MathTest<3>::testPolynomial));
        add( testCase(&MathTest<4>::testPolynomial));
        add( testCase(&MathTest<5>::testPolynomial));
    }
};

int main()
{
    MathTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
