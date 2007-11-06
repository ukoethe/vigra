/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#if defined(__GNUC__) && __GNUC__ == 2 && __GNUC_MINOR__ == 95
// deactivate broken std::relops
#  define __SGI_STL_INTERNAL_RELOPS
#endif  // __GNUC__

#include <typeinfo>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include "unittest.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/polynomial.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/splines.hxx"
#include "vigra/gaussians.hxx"
#include "vigra/rational.hxx"
#include "vigra/fixedpoint.hxx"
#include "vigra/linear_algebra.hxx"
#include "vigra/singular_value_decomposition.hxx"
#include "vigra/regression.hxx"

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
    { C(1e-11), C(1.0), C(2.0), C(3.0), C(4.0), C(5.0), C(6.0), C(7.0), C(8.0) },
    { C(1e-12), C(-1e10), C(1.0), C(1e10) },
    { C(1e-12), C(-1e-10), C(1e-10), C(1.0) },
    { C(1e-5), C(0.1), C(0.1), C(0.1), C(0.1), C(0.2), C(0.2), C(0.2),
               C(0.3), C(0.3), C(0.4) },
    { C(1e-12), C(-391.74516023901123),
                C(-0.56839260551055271, -0.4046562986541693), C(-0.56839260551055271, 0.4046562986541693),
                C(0.47331479192572767, -0.89542786425410759), C(0.47331479192572767, 0.89542786425410759) },
    { C(1e-12), C(1e-6), C(1.0), C(1e6) },
    { C(1e-12), C(0.0, -3.0), C(0.0, -2.0), C(0.0, -1.0), C(0.0, -1.0),
                C(0.0, 1.0), C(0.0, 1.0), C(0.0, 2.0), C(0.0, 3.0) }
};

#if 0
#undef should
#define should(v) (std::cerr << #v << ": " << (v) << std::endl)
#undef shouldEqual
#define shouldEqual(v1, v2) \
(std::cerr << #v1 << " == " << #v2 << ": " << (v1) << " " << (v2) << std::endl)
#undef shouldEqualTolerance
#define shouldEqualTolerance(v1, v2, e) \
(std::cerr << #v1 << " == " << #v2 << ": " << (v1) << " " << (v2) << " " << (v1 - v2) << std::endl)
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

    void testPolynomialEigenvalueMethod()
    {
        double epsilon = 1e-7;
        unsigned int order = (unsigned int)(coefficients[N][0] + 0.5);
        POLYNOMIAL p(coefficients[N]+1, order);
        p.setEpsilon(epsilon);

        vigra::ArrayVector<std::complex<double> > roots;

        should(polynomialRootsEigenvalueMethod(p, roots));
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
        shouldEqualTolerance(rroots[0], -1.0, epsilon);
        shouldEqualTolerance(rroots[1], 1.0, epsilon);
    }

    void testPolynomialEigenvalueMethod()
    {
        unsigned int order = 80;
        double epsilon = 1e-12;
        vigra::ArrayVector<double> coeffs(order+1, 0.0);
        coeffs[0] = -1.0;
        coeffs[order] = 1.0;
        vigra::Polynomial<double> p(coeffs.begin(), order);

        vigra::ArrayVector<std::complex<double> > roots;

        should(vigra::polynomialRootsEigenvalueMethod(p, roots));
        shouldEqual(roots.size(), order);
        for(unsigned int i = 0; i<roots.size(); ++i)
        {
            shouldEqualTolerance(std::abs(roots[i]), 1.0, epsilon);
            C r = p(roots[i]);
            shouldEqualTolerance(r.real(), 0.0, epsilon);
            shouldEqualTolerance(r.imag(), 0.0, epsilon);
        }
        vigra::ArrayVector<double> rroots;
        should(polynomialRealRootsEigenvalueMethod(p, rroots));
        shouldEqual(rroots.size(), 2);
        shouldEqualTolerance(rroots[0], -1.0, epsilon);
        shouldEqualTolerance(rroots[1], 1.0, epsilon);
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

    void testFixedPointValues()
    {
        double r = spline.radius();
        shouldEqual(r, splineBase.radius());

        for(double x = -r-0.5; x <= r+0.5; x += 0.5)
        {
            vigra::FixedPoint<11,20> fpx20(x);
            vigra::FixedPoint<11,15> fpx15(x);
            should(vigra::abs(vigra::fixed_point_cast<double>(spline(fpx20)) - spline(x)) < 1e-6);
            should(vigra::abs(vigra::fixed_point_cast<double>(spline(fpx15)) - spline(x)) < 4e-5);

        }
    }

    void testPrefilterCoefficients()
    {
        int n = ORDER / 2;
        vigra::ArrayVector<double> const & ps = spline.prefilterCoefficients();
        vigra::ArrayVector<double> const & psb = splineBase.prefilterCoefficients();

        if(n == 0)
        {
            shouldEqual(ps.size(), 0);
            shouldEqual(psb.size(), 0);
        }
        else
        {
            vigra::ArrayVector<double> & psb1 =
                const_cast<vigra::ArrayVector<double> &>(psb);
            std::sort(psb1.begin(), psb1.end());

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
    
    void intSquareRootTest()
    {
        for(int i = 0; i < 1024; ++i)
        {
            shouldEqual(vigra::sqrti(i), (vigra::Int32)vigra::floor(vigra::sqrt((double)i)));
        }
    }

    void testSpecialFunctions()
    {
        shouldEqualTolerance(vigra::ellipticIntegralE(M_PI / 2.0, 0.0), M_PI / 2.0, 1e-14);
        shouldEqualTolerance(vigra::ellipticIntegralF(0.3, 0.3), 0.30039919311549118, 1e-14);
        shouldEqualTolerance(vigra::ellipticIntegralE(0.3, 0.3), 0.29960175507025716, 1e-14);

        shouldEqualTolerance(vigra::erf(0.3), 0.32862675945912745, 1e-7);
        
        should(vigra::noncentralChi2CDFApprox(200, 0.0, 200.0) > 0.5);
        should(vigra::noncentralChi2CDFApprox(200, 0.0, 199.0) < 0.5);
        should(vigra::noncentralChi2CDF(200, 0.0, 200.0) > 0.5);
        should(vigra::noncentralChi2CDF(200, 0.0, 199.0) < 0.5);
        
        shouldEqualTolerance(vigra::noncentralChi2CDF(2, 2.0, 2.0), 0.34574583872316456, 1e-7);
        shouldEqualTolerance(vigra::noncentralChi2(2, 2.0, 2.0), 0.154254161276835, 1e-7);
        shouldEqualTolerance(vigra::noncentralChi2CDF(3, 2.0, 2.0), 0.22073308707450343, 1e-7);
        shouldEqualTolerance(vigra::noncentralChi2(3, 2.0, 2.0), 0.13846402271767755, 1e-7);
        shouldEqualTolerance(vigra::noncentralChi2CDFApprox(2, 2.0, 2.0), 0.34574583872316456, 1e-1);
        shouldEqualTolerance(vigra::noncentralChi2CDFApprox(3, 2.0, 2.0), 0.22073308707450343, 1e-1);
    }

    void closeAtToleranceTest()
    {
        double a = 0.0, b = vigra::NumericTraits<double>::epsilon(), c = 1000.0, d = 1000.1;

        using vigra::closeAtTolerance;
        should(closeAtTolerance(a, b));
        should(closeAtTolerance(c, c + b));
        should(closeAtTolerance(c, d, 1.0));
        should(closeAtTolerance(-a, -b));
        should(closeAtTolerance(-c, -c + b));
        should(closeAtTolerance(-c, -d, 1.0));
        should(!closeAtTolerance(c, -c));
        should(!closeAtTolerance(a, c));
        should(!closeAtTolerance(c, d));
        should(!closeAtTolerance(-a, -c));
        should(!closeAtTolerance(-c, -d));
    }

    void testArgMinMax()
    {
        using namespace vigra;
        using namespace vigra::functor;
        
        double data[] = {1.0, 5.0,
                         3.0, 2.0,
                        -2.0, 4.0};
        double *end = data + 6;
        
        shouldEqual(argMin(data, end), data+4);
        shouldEqual(argMax(data, end), data+1);
        shouldEqual(argMinIf(data, end, Arg1() > Param(0.0)), data);
        shouldEqual(argMinIf(data, end, Arg1() > Param(5.0)), end);
        shouldEqual(argMaxIf(data, end, Arg1() < Param(5.0)), data+5);
        shouldEqual(argMaxIf(data, end, Arg1() < Param(-2.0)), end);
    }
};

struct RationalTest
{
    typedef vigra::Rational<int> R;

    void testGcdLcm()
    {
        shouldEqual(vigra::gcd(24, 18), 6);
        shouldEqual(vigra::lcm(6, 4), 12);
        shouldEqual(vigra::gcd(18, 24), 6);
        shouldEqual(vigra::lcm(4, 6), 12);
    }

    void testOStreamShifting()
    {
        std::ostringstream out;
        out << R(1,2);
        out << "Testing.." << R(42,23) << 3.141592653589793238 << std::endl;
    }

    void testOperators()
    {
        shouldEqual(R(3,4), R(3,4));
        shouldEqual(-R(3,4), R(-3,4));

        shouldEqual(R(3,4) + R(12,6), R(11,4));
        shouldEqual(R(3,4) - R(12,6), R(-5,4));
        shouldEqual(R(3,4) * R(12,6), R(3,2));
        shouldEqual(R(3,4) / R(12,6), R(3,8));

        shouldEqual(abs(R(-3,4)), R(3,4));
        shouldEqual(norm(R(-3,4)), R(3,4));
        shouldEqual(squaredNorm(R(-3,4)), R(9,16));

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

    void testConversion()
    {
        shouldEqual(vigra::rational_cast<R>(R(3,2)), R(3,2));
        shouldEqual(vigra::rational_cast<int>(R(3,2)), 1);
        shouldEqual(vigra::rational_cast<double>(R(3,2)), 1.5);
        shouldEqual(vigra::rational_cast<double>(1.5), 1.5);

        shouldEqual(R(vigra::Rational<short>((short)-2, (short)-4)), R(1,2));

        shouldEqual(R(3.5, 1e-4), R(7,2));
        shouldEqual(R(-3.5, 1e-4), R(-7,2));
        shouldEqual(R(0.123, 1e-4), R(123,1000));
        shouldEqual(R(-0.123, 1e-4), R(-123,1000));
        shouldEqual(R(0.123456, 1e-4), R(1235,10000));
        shouldEqual(R(0.123432, 1e-4), R(1234,10000));
        shouldEqual(R(-0.123456, 1e-4), R(-1235,10000));
        shouldEqual(R(-0.123432, 1e-4), R(-1234,10000));
    }

    void testFunctions()
    {
        shouldEqual(pow(R(1,2),2), R(1,4));
        shouldEqual(pow(R(2),-2), R(1,4));
        shouldEqual(pow(R(-1,2),2), R(1,4));
        shouldEqual(pow(R(-2),-2), R(1,4));
        shouldEqual(pow(R(-1,2),3), R(-1,8));
        shouldEqual(pow(R(-2),-3), R(-1,8));
        shouldEqual(pow(R(3),0), R(1));
        shouldEqual(pow(R(0),3), R(0));
        shouldEqual(pow(R(0),0), R(1));
        should(pow(R(0),-3).is_pinf());

        should(pow(R(1,0, false), 1).is_pinf());
        should(pow(R(-1,0, false), 1).is_ninf());
        shouldEqual(pow(R(1,0, false), -1), R(0));
        shouldEqual(pow(R(-1,0, false), -1), R(0));
        try { pow(R(1,0, false), 0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { pow(R(-1,0, false), 0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}

        shouldEqual(floor(R(2)), R(2));
        shouldEqual(floor(R(3,2)), R(1));
        shouldEqual(floor(R(1,2)), R(0));
        shouldEqual(floor(R(-1,2)), R(-1));
        shouldEqual(floor(R(1,-2)), R(-1));
        shouldEqual(floor(R(-3,2)), R(-2));
        shouldEqual(floor(R(-2)), R(-2));
        shouldEqual(floor(R(1,0,false)), R(1,0,false));
        shouldEqual(floor(R(-1,0,false)), R(-1,0,false));

        shouldEqual(ceil(R(2)), R(2));
        shouldEqual(ceil(R(3,2)), R(2));
        shouldEqual(ceil(R(1,2)), R(1));
        shouldEqual(ceil(R(-1,2)), R(0));
        shouldEqual(ceil(R(1,-2)), R(0));
        shouldEqual(ceil(R(-3,2)), R(-1));
        shouldEqual(ceil(R(-2)), R(-2));
        shouldEqual(ceil(R(1,0,false)), R(1,0,false));
        shouldEqual(ceil(R(-1,0,false)), R(-1,0,false));
    }

    void testInf()
    {
        R inf(2,0);
        R ninf(-2,0);

        should(inf.is_inf());
        should(inf.is_pinf());
        should(!inf.is_ninf());
        should(ninf.is_inf());
        should(ninf.is_ninf());
        should(!ninf.is_pinf());
        shouldEqual(inf.numerator(), 1);
        shouldEqual(ninf.numerator(), -1);

        should((inf + R(1)).is_pinf());
        should((inf + R(0)).is_pinf());
        should((inf + R(-1)).is_pinf());
        should((ninf + R(1)).is_ninf());
        should((ninf + R(0)).is_ninf());
        should((ninf + R(-1)).is_ninf());
        should((inf + 1).is_pinf());
        should((inf + 0).is_pinf());
        should((inf + (-1)).is_pinf());
        should((ninf + 1).is_ninf());
        should((ninf + 0).is_ninf());
        should((ninf + (-1)).is_ninf());
        should((inf + inf).is_pinf());
        should((ninf + ninf).is_ninf());
        shouldEqual((inf + R(3)).numerator(), 1);
        shouldEqual((ninf + R(3)).numerator(), -1);

        should((inf - R(1)).is_pinf());
        should((inf - R(0)).is_pinf());
        should((inf - R(-1)).is_pinf());
        should((ninf - R(1)).is_ninf());
        should((ninf - R(0)).is_ninf());
        should((ninf - R(-1)).is_ninf());
        should((inf - 1).is_pinf());
        should((inf - 0).is_pinf());
        should((inf - (-1)).is_pinf());
        should((ninf - 1).is_ninf());
        should((ninf - 0).is_ninf());
        should((ninf - (-1)).is_ninf());
        should((inf - ninf).is_pinf());
        should((ninf - inf).is_ninf());
        shouldEqual((inf - R(3)).numerator(), 1);
        shouldEqual((ninf - R(3)).numerator(), -1);

        should((inf * R(1)).is_pinf());
        should((inf * R(-1)).is_ninf());
        should((ninf * R(1)).is_ninf());
        should((ninf * R(-1)).is_pinf());
        should((inf * 1).is_pinf());
        should((inf * (-1)).is_ninf());
        should((ninf * 1).is_ninf());
        should((ninf * (-1)).is_pinf());
        should((inf * inf).is_pinf());
        should((inf * ninf).is_ninf());
        should((ninf * inf).is_ninf());
        should((ninf * ninf).is_pinf());
        shouldEqual((inf * R(3)).numerator(), 1);
        shouldEqual((ninf * R(3)).numerator(), -1);
        shouldEqual((inf * R(-3)).numerator(), -1);
        shouldEqual((ninf * R(-3)).numerator(), 1);

        should((inf / R(1)).is_pinf());
        should((inf / R(0)).is_pinf());
        should((inf / R(-1)).is_ninf());
        should((ninf / R(1)).is_ninf());
        should((ninf / R(0)).is_ninf());
        should((ninf / R(-1)).is_pinf());
        shouldEqual(R(1) / inf, R(0));
        shouldEqual(R(-1) / inf, R(0));
        shouldEqual(R(1) / ninf, R(0));
        shouldEqual(R(-1) / ninf, R(0));
        should((inf / 1).is_pinf());
        should((inf / 0).is_pinf());
        should((inf / (-1)).is_ninf());
        should((ninf / 1).is_ninf());
        should((ninf / 0).is_ninf());
        should((ninf / (-1)).is_pinf());

        shouldEqual(2 / inf, R(0));
        shouldEqual((-2) / inf, R(0));
        shouldEqual(2 / ninf, R(0));
        shouldEqual((-2) / ninf, R(0));
        shouldEqual((2 / inf).denominator(), 1);
        shouldEqual(((-2) / inf).denominator(), 1);
        shouldEqual((2 / ninf).denominator(), 1);
        shouldEqual(((-2) / ninf).denominator(), 1);

        shouldEqual((inf / R(3)).numerator(), 1);
        shouldEqual((ninf / R(3)).numerator(), -1);
        shouldEqual((inf / R(-3)).numerator(), -1);
        shouldEqual((ninf / R(-3)).numerator(), 1);

        should(inf == inf);
        should(!(inf != inf));
        should(!(inf < inf));
        should(inf <= inf);
        should(!(inf > inf));
        should(inf >= inf);
        should(ninf == ninf);
        should(!(ninf != ninf));
        should(!(ninf < ninf));
        should(ninf <= ninf);
        should(!(ninf > ninf));
        should(ninf >= ninf);
        should(inf != ninf);
        should(ninf != inf);
        should(inf > ninf);
        should(ninf < inf);
        should(!(inf < ninf));
        should(!(ninf > inf));

        should(inf != 0);
        should(ninf != 0);
        should(inf > 0);
        should(inf >= 0);
        should(ninf < 0);
        should(ninf <= 0);
        should(!(0 < ninf));
        should(!(0 > inf));
        should(!(0 <= ninf));
        should(!(0 >= inf));

        should(inf != R(1));
        should(ninf != R(1));
        should(inf > R(1));
        should(inf >= R(1));
        should(ninf < R(1));
        should(ninf <= R(1));
        should(!(R(1) < ninf));
        should(!(R(1) > inf));
        should(!(R(1) <= ninf));
        should(!(R(1) >= inf));

        try { inf + ninf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { ninf + inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { inf - inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { ninf - ninf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { inf * R(0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { ninf * R(0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { R(0) * inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { R(0) * ninf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { inf * 0; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { ninf * 0; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { 0 * inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { 0 * ninf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { inf / inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { inf / ninf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { ninf / inf; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { R(0) / R(0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { R(0) / 0; failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
        try { 0 / R(0); failTest("No exception thrown"); } catch(vigra::bad_rational &) {}
    }
};

struct FixedPointTest
{
    void testConstruction()
    {
        shouldEqual(vigra::fixedPoint(3).value, 3);
        shouldEqual(vigra::fixedPoint(-3).value, -3);
        shouldEqual(-vigra::fixedPoint(3).value, -3);

        shouldEqual((vigra::FixedPoint<3,4>(3).value), 3 << 4);
        shouldEqual((vigra::FixedPoint<3,4>(-3).value), -3 << 4);
        shouldEqual((-vigra::FixedPoint<3,4>(3).value), -3 << 4);

        shouldEqual((vigra::FixedPoint<3,4>(3.5).value), 56);
        shouldEqual((vigra::FixedPoint<3,4>(-3.5).value), -56);
        shouldEqual((-vigra::FixedPoint<3,4>(3.5).value), -56);

        try { vigra::FixedPoint<1, 8>(3.75); failTest("No exception thrown"); } catch(vigra::PreconditionViolation &) {}

        shouldEqual((vigra::NumericTraits<vigra::FixedPoint<1, 8> >::zero()).value, 0);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint<1, 8> >::one()).value, 1 << 8);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint<1, 8> >::max()).value, (1 << 9) - 1);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint<1, 8> >::min()).value, -((1 << 9) - 1));

        vigra::FixedPoint<2, 8> v(3.75);
        shouldEqual((vigra::FixedPoint<2, 8>(v).value), 15 << 6);
        shouldEqual((vigra::FixedPoint<3, 10>(v).value), 15 << 8);
        shouldEqual((vigra::FixedPoint<2, 2>(v).value), 15);
        shouldEqual((vigra::FixedPoint<2, 0>(v).value), 4);

        shouldEqual((vigra::FixedPoint<2, 8>(-v).value), -15 << 6);
        shouldEqual((vigra::FixedPoint<3, 10>(-v).value), -15 << 8);
        shouldEqual((vigra::FixedPoint<2, 2>(-v).value), -15);
        shouldEqual((vigra::FixedPoint<2, 0>(-v).value), -4);

        shouldEqual(vigra::fixed_point_cast<double>(v), 3.75);
        should((frac(v) == vigra::FixedPoint<0, 8>(0.75)));
        should((dual_frac(v) == vigra::FixedPoint<0, 8>(0.25)));
        should(vigra::floor(v) == 3);
        should(vigra::ceil(v) == 4);
        should(round(v) == 4);
        should(vigra::abs(v) == v);
        should((frac(-v) == vigra::FixedPoint<0, 8>(0.25)));
        should((dual_frac(-v) == vigra::FixedPoint<0, 8>(0.75)));
        should(vigra::floor(-v) == -4);
        should(vigra::ceil(-v) == -3);
        should(round(-v) == -4);
        should(vigra::abs(-v) == v);
        should(vigra::norm(-v) == v);
        should(vigra::squaredNorm(-v) == v*v);

        vigra::FixedPoint<3, 10> v1;
        shouldEqual((v1 = v).value, 15 << 8);
        shouldEqual((v1 = -v).value, -15 << 8);

        vigra::FixedPoint<2, 0> v2;
        shouldEqual((v2 = v).value, 4);
        shouldEqual((v2 = -v).value, -4);
    }

    void testComparison()
    {
        vigra::FixedPoint<3, 8> v1(3.75), v2(4);
        vigra::FixedPoint<2, 2> v3(3.75);
        should(v1 == v1);
        should(v1 == v3);
        should(!(v1 != v1));
        should(!(v1 != v3));
        should(v1 <= v1);
        should(v1 <= v3);
        should(!(v1 < v1));
        should(!(v1 < v3));
        should(v1 >= v1);
        should(v1 >= v3);
        should(!(v1 > v1));
        should(!(v1 > v3));

        should(v2 != v1);
        should(v2 != v3);
        should(!(v2 == v1));
        should(!(v2 == v3));
        should(!(v2 <= v1));
        should(!(v2 <= v3));
        should(!(v2 < v1));
        should(!(v2 < v3));
        should(v2 >= v1);
        should(v2 >= v3);
        should(v2 > v1);
        should(v2 > v3);
    }

    void testArithmetic()
    {
        vigra::FixedPoint<1, 16> t1(0.75), t2(0.25);
        signed char v1 = 1, v2 = 2, v4 = 4, v8 = 8;

        should((vigra::FixedPoint<1, 16>(t1) += t1) == (vigra::FixedPoint<1, 16>(1.5)));
        should((vigra::FixedPoint<1, 16>(t1) -= t1) == (vigra::FixedPoint<1, 16>(0.0)));
        should((vigra::FixedPoint<2, 16>(t1) *= t1) == (vigra::FixedPoint<1, 16>(9.0 / 16.0)));

        should(--t1 == (vigra::FixedPoint<1, 16>(-0.25)));
        should(t1 == (vigra::FixedPoint<1, 16>(-0.25)));
        should(++t1 == (vigra::FixedPoint<1, 16>(0.75)));
        should(t1 == (vigra::FixedPoint<1, 16>(0.75)));
        should(t1++ == (vigra::FixedPoint<1, 16>(0.75)));
        should(t1 == (vigra::FixedPoint<1, 16>(1.75)));
        should(t1-- == (vigra::FixedPoint<1, 16>(1.75)));
        should(t1 == (vigra::FixedPoint<1, 16>(0.75)));

        shouldEqual((t1 * vigra::fixedPoint(v1)).value, 3 << 14);
        shouldEqual((t2 * vigra::fixedPoint(v1)).value, 1 << 14);
        shouldEqual((-t1 * vigra::fixedPoint(v1)).value, -3 << 14);
        shouldEqual((-t2 * vigra::fixedPoint(v1)).value, -1 << 14);
        shouldEqual((t1 * -vigra::fixedPoint(v1)).value, -3 << 14);
        shouldEqual((t2 * -vigra::fixedPoint(v1)).value, -1 << 14);

        shouldEqual((vigra::FixedPoint<8, 2>(t1 * vigra::fixedPoint(v1))).value, 3);
        shouldEqual((vigra::FixedPoint<8, 2>(t2 * vigra::fixedPoint(v1))).value, 1);
        shouldEqual((vigra::FixedPoint<8, 2>(-t1 * vigra::fixedPoint(v1))).value, -3);
        shouldEqual((vigra::FixedPoint<8, 2>(-t2 * vigra::fixedPoint(v1))).value, -1);

        shouldEqual(vigra::floor(t1 * vigra::fixedPoint(v1) + t2 * vigra::fixedPoint(v2)), 1);
        shouldEqual(vigra::ceil(t1 * vigra::fixedPoint(v1) + t2 * vigra::fixedPoint(v2)), 2);
        shouldEqual(round(t1 * vigra::fixedPoint(v1) + t2 * vigra::fixedPoint(v2)), 1);
        shouldEqual(vigra::floor(t1 * vigra::fixedPoint(v4) + t2 * vigra::fixedPoint(v8)), 5);
        shouldEqual(vigra::ceil(t1 * vigra::fixedPoint(v4) + t2 * vigra::fixedPoint(v8)), 5);
        shouldEqual(round(t1 * vigra::fixedPoint(v4) + t2 * vigra::fixedPoint(v8)), 5);

        shouldEqual(vigra::floor(t1 * -vigra::fixedPoint(v1) - t2 * vigra::fixedPoint(v2)), -2);
        shouldEqual(vigra::ceil(t1 * -vigra::fixedPoint(v1) - t2 * vigra::fixedPoint(v2)), -1);
        shouldEqual(round(t1 * -vigra::fixedPoint(v1) - t2 * vigra::fixedPoint(v2)), -1);
        shouldEqual(vigra::floor(t1 * -vigra::fixedPoint(v4) - t2 * vigra::fixedPoint(v8)), -5);
        shouldEqual(vigra::ceil(t1 * -vigra::fixedPoint(v4) - t2 * vigra::fixedPoint(v8)), -5);
        shouldEqual(round(t1 * -vigra::fixedPoint(v4) - t2 * vigra::fixedPoint(v8)), -5);

        double d1 = 1.0 / 3.0, d2 = 1.0 / 7.0;
        vigra::FixedPoint<1, 24> r1(d1), r2(d2);
        vigra::FixedPoint<2, 24> r3;
        add(r1, r2, r3);
        shouldEqual(r3.value, (vigra::FixedPoint<2, 24>(d1 + d2)).value);
        sub(r1, r2, r3);
        shouldEqual(r3.value, (vigra::FixedPoint<2, 24>(d1 - d2)).value);
        mul(r1, r2, r3);
        shouldEqual(r3.value >> 2, (vigra::FixedPoint<2, 24>(d1 * d2)).value >> 2);
        
        for(int i = 0; i < 1024; ++i)
        {
            vigra::FixedPoint<4,5> fv1(i, vigra::FPNoShift);
            vigra::FixedPoint<5,4> fv2(i, vigra::FPNoShift);
            vigra::FixedPoint<5,5> fv3(i, vigra::FPNoShift);
            vigra::FixedPoint<6,6> fv4(i, vigra::FPNoShift);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv1)), vigra::floor(vigra::sqrt((double)fv1.value)) / 8.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv2)), vigra::floor(vigra::sqrt((double)fv2.value)) / 4.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv3)), vigra::floor(vigra::sqrt((double)fv3.value)) / 4.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv4)), vigra::floor(vigra::sqrt((double)fv4.value)) / 8.0);
        }
    }
};

struct LinalgTest
{
    typedef vigra::Matrix<double> Matrix;
    typedef Matrix::difference_type Shape;

    unsigned int size, iterations;

    LinalgTest()
    : size(50),
      iterations(5)
    {
        std::srand (0xdeadbeef);
    }

    void testOStreamShifting()
    {
        Matrix a = random_matrix (size, size);
        std::ostringstream out;
        out << a;
        out << "Testing.." << a << 42 << std::endl;
    }

    static double random_double ()
    {
        double ret = 2.0 * static_cast <double> (std::rand ()) / RAND_MAX - 1.0;
        return ret;
    }

    static Matrix random_matrix(unsigned int rows, unsigned int cols)
    {
        Matrix ret (rows, cols);
        for (unsigned int i = 0; i < rows; ++i)
            for (unsigned int j = 0; j < cols; ++j)
                ret (i, j) = random_double ();
        return ret;
    }

    static Matrix random_symmetric_matrix(unsigned int rows)
    {
        Matrix ret (rows, rows);
        for (unsigned int i = 0; i < rows; ++i)
            for (unsigned int j = i; j < rows; ++j)
                ret (j, i) = ret (i, j) = random_double ();
        return ret;
    }
    
    void testMatrix()
    {
        double data[] = {1.0, 5.0,
                         3.0, 2.0,
                         4.0, 7.0};
        double tref[] = {1.0, 3.0, 4.0,
                         5.0, 2.0, 7.0};
        double tref2[] = {1.0, 3.0,
                          5.0, 2.0};
        double idref[] = {1.0, 0.0, 0.0,
                          0.0, 1.0, 0.0,
                          0.0, 0.0, 1.0};
        std::string sref(" 1.0000  5.0000 \n 3.0000  2.0000 \n 4.0000  7.0000 \n");
        unsigned int r = 3, c = 2;

        Matrix a(r, c, data), zero(r, c);
        shouldEqual(a.rowCount(), r);
        shouldEqual(a.columnCount(), c);
        shouldEqual(a.elementCount(), r*c);
        shouldEqual(a.squaredNorm(), 104.0);
        shouldEqual(a.norm(), std::sqrt(104.0));
        shouldEqual(a.squaredNorm(), squaredNorm(a));
        shouldEqual(a.norm(), norm(a));
        shouldEqual(rowCount(a), r);
        shouldEqual(columnCount(a), c);
        
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(zero(i,j), 0.0);

        std::stringstream s;
        s << a;
        shouldEqual(s.str(), sref);

        for(unsigned int i=0, k=0; i<r; ++i)
        {
            Matrix::view_type ar = a.rowVector(i);
            shouldEqual(rowCount(ar), 1);
            shouldEqual(columnCount(ar), c);
            Matrix::view_type ar1 = rowVector(a, i);
            shouldEqual(rowCount(ar1), 1);
            shouldEqual(columnCount(ar1), c);
            for(unsigned int j=0; j<c; ++j, ++k)
            {
                shouldEqual(a(i,j), data[k]);
                shouldEqual(ar(0, j), data[k]);
                shouldEqual(ar1(0, j), data[k]);
            }
        }

        Matrix aa(r, c, tref, vigra::ColumnMajor);
        shouldEqual(aa.rowCount(), r);
        shouldEqual(aa.columnCount(), c);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(aa(i,j), a(i,j));

        Matrix b = a;
        shouldEqual(b.rowCount(), r);
        shouldEqual(b.columnCount(), c);
        shouldEqualSequence(a.begin(), a.end(), b.begin());

        b.init(0.0);
        should(b == zero);
        
        Matrix::iterator ib = b.begin();
        b = a;
        shouldEqual(ib, b.begin()); 
        shouldEqualSequence(a.begin(), a.end(), b.begin());

        b = 4.0 + a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0+data[k]);

        b = a + 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]+3.0);

        b = 4.0 - a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0-data[k]);

        b = a - 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]-3.0);

        b = 4.0 * a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0*data[k]);

        b = a * 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]*3.0);

        b = 4.0 / a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0/data[k]);

        b = a / 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k] / 3.0);

        b = a + a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 2.0 * data[k]);

        b = a - a;
        for(unsigned int i=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j)
                shouldEqual(b(i,j), 0.0);

        b = -a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), -data[k]);

        b = a * pointWise(a);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k] * data[k]);

        b = a / pointWise(a);
        for(unsigned int i=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j)
                shouldEqual(b(i,j), 1.0);

        b = pow(a, 2);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k] * data[k]);

        b = sqrt(a);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), sqrt(data[k]));

        b = sq(a);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), vigra::sq(data[k]));

        b = sign(a);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), vigra::sign(data[k]));

        Matrix at = transpose(a);
        shouldEqual(at.rowCount(), c);
        shouldEqual(at.columnCount(), r);
        for(unsigned int i=0, k=0; i<c; ++i)
        {
            Matrix::view_type ac = a.columnVector(i);
            shouldEqual(rowCount(ac), r);
            shouldEqual(columnCount(ac), 1);
            Matrix::view_type ac1 = columnVector(a, i);
            shouldEqual(rowCount(ac1), r);
            shouldEqual(columnCount(ac1), 1);
            for(unsigned int j=0; j<r; ++j, ++k)
            {
                shouldEqual(at(i,j), tref[k]);
                shouldEqual(ac(j,0), tref[k]);
                shouldEqual(ac1(j,0), tref[k]);
            }
        }
        
        double sn = squaredNorm(columnVector(a, 0));
        shouldEqual(sn, 26.0);
        shouldEqual(sn, dot(columnVector(a, 0), columnVector(a, 0)));
        shouldEqual(sn, dot(rowVector(at, 0), columnVector(a, 0)));
        shouldEqual(sn, dot(columnVector(a, 0), rowVector(at, 0)));
        shouldEqual(sn, dot(rowVector(at, 0), rowVector(at, 0)));
        shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(1,0)), a.subarray(Shape(0,0), Shape(0,1))));
        shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(0,1)), a.subarray(Shape(0,0), Shape(0,1))));
        shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(1,0)), a.subarray(Shape(0,0), Shape(1,0))));
        shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(0,1)), a.subarray(Shape(0,0), Shape(1,0))));

        Matrix a2(c, c, data);
        a2 = a2.transpose();
        for(unsigned int i=0, k=0; i<c; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(a2(i,j), tref2[k]);

        Matrix id = vigra::identityMatrix<double>(r);
        shouldEqual(id.rowCount(), r);
        shouldEqual(id.columnCount(), r);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<r; ++j, ++k)
                shouldEqual(id(i,j), idref[k]);

        Matrix d = diagonalMatrix(Matrix(r, 1, data));
        shouldEqual(d.rowCount(), r);
        shouldEqual(d.columnCount(), r);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<r; ++j, ++k)
                shouldEqual(d(i,j), idref[k]*data[i]);

        Matrix e(r*c, 1, data);
        shouldEqual(dot(transpose(e), e), e.squaredNorm());

        double dc1[] = {1.0, 1.0, 1.0},
               dc2[] = {1.2, 2.4, 3.6};
        Matrix c1(3,1, dc1), c2(3,1, dc2);
        Matrix cr = cross(c1, c2);
        shouldEqualTolerance(cr(0,0), 1.2, 1e-12);
        shouldEqualTolerance(cr(1,0), -2.4, 1e-12);
        shouldEqualTolerance(cr(2,0), 1.2, 1e-12);

        Matrix f(1, r*c - 1, tref);
        Matrix g = outer(e, f);
        shouldEqual(g.rowCount(), e.rowCount());
        shouldEqual(g.columnCount(), f.columnCount());
        for(unsigned int i=0; i<g.rowCount(); ++i)
            for(unsigned int j=0; j<g.columnCount(); ++j)
                shouldEqual(g(i,j), data[i]*tref[j]);

        Matrix h = transpose(a) * a;
        shouldEqual(h.rowCount(), c);
        shouldEqual(h.columnCount(), c);
        for(unsigned int i=0; i<c; ++i)
            for(unsigned int j=0; j<c; ++j)
                shouldEqual(h(i,j), dot(rowVector(at, i), columnVector(a, j)));

        should(isSymmetric(random_symmetric_matrix(10)));
        should(!isSymmetric(random_matrix(10, 10)));

        Matrix tm(2, 2, tref2);
        vigra::TinyVector<double, 2> tv(1.0, 2.0), tvrref(7.0, 9.0), tvlref(11.0, 7.0);
        shouldEqual(tm * tv, tvrref);
        shouldEqual(tv * tm, tvlref);
    }
 
    void testArgMinMax()
    {
        using namespace vigra::functor;
        
        double data[] = {1.0, 5.0,
                         3.0, 2.0,
                        -2.0, 4.0};
        unsigned int r = 3, c = 2;
        Matrix minmax(r, c, data);
        
        shouldEqual(argMin(minmax), 2);
        shouldEqual(argMax(minmax), 3);
        shouldEqual(argMinIf(minmax, Arg1() > Param(0.0)), 0);
        shouldEqual(argMinIf(minmax, Arg1() > Param(5.0)), -1);
        shouldEqual(argMaxIf(minmax, Arg1() < Param(5.0)), 5);
        shouldEqual(argMaxIf(minmax, Arg1() < Param(-2.0)), -1);
    }

    void testCholesky()
    {
        double epsilon = 1e-10;
        Matrix idref = vigra::identityMatrix<double>(size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            a = transpose(a) * a; // make a symmetric positive definite matrix
            Matrix l(size, size);
            choleskyDecomposition (a, l);
            Matrix ch = l * transpose(l);
            shouldEqualSequenceTolerance(ch.data(), ch.data()+size*size, a.data(), epsilon);
        }
    }

    void testQR()
    {
        double epsilon = 1e-10;
        Matrix idref = vigra::identityMatrix<double>(size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            Matrix r(size, size);
            Matrix q(size, size);
            qrDecomposition (a, q, r);
            Matrix id = transpose(q) * q;
            shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            Matrix qr = q * r;
            shouldEqualSequenceTolerance(qr.data(), qr.data()+size*size, a.data(), epsilon);
        }
    }

    void testLinearSolve()
    {
#if defined(__GNUC__) && __GNUC__ == 3
        double epsilon = 1e-8;
#else
        double epsilon = 1e-10;
#endif
        int size = 50;

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            Matrix b = random_matrix (size, 1);
            Matrix x(size, 1);

            should(linearSolve (a, b, x, "ExplicitQR"));
            Matrix ax = a * x;
            shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            should(linearSolve(a, b, x, "QR"));
            ax = a * x;
            shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            should(linearSolve(a, b, x, "SVD"));
            ax = a * x;
            shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);
            
            should(linearSolve(a, b, x, "NE"));
            ax = a * x;
            shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            Matrix c = transpose(a) * a; // make a symmetric positive definite matrix
            should(linearSolve (c, b, x, "Cholesky"));
            ax = c * x;
            shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);
        }
    }
    
    void testUnderdetermined()
    {   
        // test singular matrix
        Matrix a = vigra::identityMatrix<Matrix::value_type> (size);
        a(0,0) = 0; 
        Matrix b = random_matrix (size, 1);
        Matrix x(size, 1);
        should(!linearSolve (a, b, x, "Cholesky"));
        should(!linearSolve (a, b, x, "QR"));
        should(!linearSolve (a, b, x, "SVD"));
    }

    void testOverdetermined()
    {   
        double epsilon = 1e-12;
        
        unsigned int n = 5;
        unsigned int size = 1000;
        double noiseStdDev = 0.1;
        
        Matrix A(size, n), xs(n,1), xq(n,1), xn(n,1), r(size, 1);
        Matrix mean(1,1), stdDev(1, 1);

        for(unsigned int iter=0; iter<iterations; ++iter)
        {
            Matrix weights = random_matrix (n, 1);
            Matrix v = random_matrix (size, 1);
            
            // init rhs with Gaussian noise with zero mean and noiseStdDev
            Matrix rhs = 0.5*noiseStdDev*random_matrix (size, 1);
            for(unsigned int k=1; k<12; ++k)
                rhs += 0.5*noiseStdDev*random_matrix (size, 1);

            for(unsigned int k=0; k<size; ++k)
            {
                for(unsigned int l=0; l<n; ++l)
                {
                    A(k,l) = std::pow(v(k,0), double(l));
                    rhs(k,0) += weights(l,0)*A(k,l);
                }
            }       
            shouldEqual(linearSolve(A, rhs, xs, "SVD"), true);
            columnStatistics(A*xs-rhs, mean, stdDev);
            
            shouldEqualTolerance(mean(0,0), 0, 1e-12);
            shouldEqualTolerance(stdDev(0,0), noiseStdDev, 0.2);

            shouldEqual(linearSolveQR(A, rhs, xq), n);
            columnStatistics(A*xq-rhs, mean, stdDev);

            shouldEqualTolerance(mean(0,0), 0, 1e-12);
            shouldEqualTolerance(stdDev(0,0), noiseStdDev, 0.2);

            shouldEqual(linearSolve(A, rhs, xn, "ne"), true);
            columnStatistics(A*xn-rhs, mean, stdDev);
            
            shouldEqualTolerance(mean(0,0), 0, 1e-12);
            shouldEqualTolerance(stdDev(0,0), noiseStdDev, 0.2);
            
            shouldEqualSequenceTolerance(xs.data(), xs.data()+n, xq.data(), epsilon);
            shouldEqualSequenceTolerance(xs.data(), xs.data()+n, xn.data(), epsilon);
        }
    }

    void testIncrementalLinearSolve()
    {
#if defined(__GNUC__) && __GNUC__ == 3
        double epsilon = 1e-8;
#else
        double epsilon = 1e-10;
#endif
        int size = 50;

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            Matrix b = random_matrix (size, 1);
            Matrix x(size, 1);

            should(linearSolve(a, b, x, "QR"));
            
            {
                Matrix r(a), qtb(b), px(size,1), xx(size,1);
                vigra::ArrayVector<unsigned int> permutation(size);
                
                for(int k=0; k<size; ++k)
                {
                    // use Givens steps for a change (Householder steps like
                    //    should(vigra::linalg::detail::qrColumnHouseholderStep(k, r, qtb));
                    // work as well, but are already extensively tested within the QR algorithm)
                    should(vigra::linalg::detail::qrGivensStepImpl(k, r, qtb));
                    permutation[k] = k;
                }
                
                for(int k=0; k<size; ++k)
                {
                    int i = rand() % size, j = rand() % size;
                    if(i==j) continue;
                    
                    vigra::linalg::detail::upperTriangularCyclicShiftColumns(i, j, r, qtb, permutation);
                }
                should(vigra::linalg::linearSolveUpperTriangular(r, qtb, px));
                vigra::linalg::detail::inverseRowPermutation(px, xx, permutation);

                shouldEqualSequenceTolerance(x.data(), x.data()+size, xx.data(), epsilon);
            }

            {
                Matrix r(a), qtb(b), px(size,1), xx(size,1);
                vigra::ArrayVector<unsigned int> permutation(size);
                
                for(int k=0; k<size; ++k)
                {
                    // use Givens steps for a change (Householder steps like
                    //    should(vigra::linalg::detail::qrColumnHouseholderStep(k, r, qtb));
                    // work as well, but are already extensively tested within the QR algorithm)
                    should(vigra::linalg::detail::qrGivensStepImpl(k, r, qtb));
                    permutation[k] = k;
                }
                
                for(int k=0; k<size; ++k)
                {
                    int i = rand() % size, j = rand() % size;
                    vigra::linalg::detail::upperTriangularSwapColumns(i, j, r, qtb, permutation);
                }
                should(vigra::linalg::linearSolveUpperTriangular(r, qtb, px));
                vigra::linalg::detail::inverseRowPermutation(px, xx, permutation);
   
                shouldEqualSequenceTolerance(x.data(), x.data()+size, xx.data(), epsilon);
            }
        }
    }
    
    void testInverse()
    {
        double epsilon = 1e-10;
        Matrix idref = vigra::identityMatrix<double>(size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            Matrix id = a * inverse(a);
            shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            id = inverse(a) * a;
            shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            id = inverse(idref * a) * a; // test inverse(const TemporaryMatrix<T> &v)
            shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
        }

        double data[] = { 1.0, 0.0, 0.0, 0.0 };
        Matrix singular(2, 2, data);
        try {
            inverse(singular);
            failTest("inverse(singular) didn't throw an exception.");
        }
        catch(vigra::PreconditionViolation & c)
        {
            std::string expected("\nPrecondition violation!\ninverse(): matrix is not invertible.");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testSymmetricEigensystem()
    {
        double epsilon = 1e-8;
        Matrix idref = vigra::identityMatrix<double>(size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_symmetric_matrix (size);
            Matrix ew(size, 1);
            Matrix ev(size, size);
            should(symmetricEigensystem(a, ew, ev));
            Matrix id = ev * transpose(ev);
            shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            Matrix ae = ev * diagonalMatrix(ew) * transpose(ev);
            shouldEqualSequenceTolerance(ae.data(), ae.data()+size*size, a.data(), epsilon);
        }
    }

    void testNonsymmetricEigensystem()
    {
        double epsilon = 1e-8;
        Matrix idref = vigra::identityMatrix<double>(size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);
            vigra::Matrix<std::complex<double> > ew(size, 1);
            Matrix ev(size, size);
            should(nonsymmetricEigensystem(a, ew, ev));

            Matrix ewm(size, size);
            for(unsigned int k = 0; k < size; k++)
            {
                ewm(k, k) = ew(k, 0).real();
                if(ew(k, 0).imag() > 0.0)
                {
                    ewm(k, k+1) = ew(k, 0).imag();
                }
                else if(ew(k, 0).imag() < 0.0)
                {
                    ewm(k, k-1) = ew(k, 0).imag();
                }
            }
            Matrix ae = ev * ewm * inverse(ev);
            shouldEqualSequenceTolerance(ae.data(), ae.data()+size*size, a.data(), epsilon);
        }
    }
    
    void testDeterminant()
    {
        double ds2[] = {1, 2, 2, 1};
        double dns2[] = {1, 2, 3, 1};
        Matrix ms2(Shape(2,2), ds2);
        Matrix mns2(Shape(2,2), dns2);
        double eps = 1e-12;
        shouldEqualTolerance(determinant(ms2), -3.0, eps);
        shouldEqualTolerance(determinant(mns2), -5.0, eps);
        shouldEqualTolerance(logDeterminant(transpose(ms2)*ms2), std::log(9.0), eps);
        shouldEqualTolerance(logDeterminant(transpose(mns2)*mns2), std::log(25.0), eps);

        double ds3[] = {1, 2, 3, 2, 3, 1, 3, 1, 2};
        double dns3[] = {1, 2, 3, 5, 3, 1, 3, 1, 2};
        Matrix ms3(Shape(3,3), ds3);
        Matrix mns3(Shape(3,3), dns3);
        shouldEqualTolerance(determinant(ms3), -18.0, eps);
        shouldEqualTolerance(determinant(mns3), -21.0, eps);
        shouldEqualTolerance(determinant(transpose(ms3)*ms3, "Cholesky"), 324.0, eps);
        shouldEqualTolerance(determinant(transpose(mns3)*mns3, "Cholesky"), 441.0, eps);
        shouldEqualTolerance(logDeterminant(transpose(ms3)*ms3), std::log(324.0), eps);
        shouldEqualTolerance(logDeterminant(transpose(mns3)*mns3), std::log(441.0), eps);
    }
    
    void testSVD()
    {
        int m = 6, n = 4;
        Matrix a(m, n);
        for(int i1= 0; i1 < m; i1++)
            for(int i2= 0; i2 < n; i2++)
                a(i1, i2)= ((float)std::rand()-RAND_MAX/2.0)/10000.0;
	    Matrix u(m, n);
	    Matrix v(n, n);
	    Matrix S(n, 1);

	    unsigned int rank = singularValueDecomposition(a, u, S, v);
	    shouldEqual(rank, n);

        double eps = 1e-10;
	    shouldEqualTolerance(norm(a-u*diagonalMatrix(S)*transpose(v)), 0.0, eps);
	    shouldEqualTolerance(norm(vigra::identityMatrix<double>(4) - transpose(u)*u), 0.0, eps);
	    shouldEqualTolerance(norm(vigra::identityMatrix<double>(4) - transpose(v)*v), 0.0, eps);
	    shouldEqualTolerance(norm(vigra::identityMatrix<double>(4) - v*transpose(v)), 0.0, eps);
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

        // test only some polynomials, as the eigenvalue method is less robust
        add( testCase((&PolynomialTest<1, P1>::testPolynomialEigenvalueMethod)));
        add( testCase((&PolynomialTest<5, P1>::testPolynomialEigenvalueMethod)));
        add( testCase((&PolynomialTest<6, P1>::testPolynomialEigenvalueMethod)));
        add( testCase((&PolynomialTest<7, P1>::testPolynomialEigenvalueMethod)));

        add( testCase((&PolynomialTest<0, P2>::testPolynomial)));
        add( testCase((&PolynomialTest<1, P2>::testPolynomial)));
        add( testCase((&PolynomialTest<2, P2>::testPolynomial)));

        add( testCase(&HighOrderPolynomialTest::testPolynomial));
        add( testCase(&HighOrderPolynomialTest::testPolynomialEigenvalueMethod));

        add( testCase(&SplineTest<0>::testValues));
        add( testCase(&SplineTest<1>::testValues));
        add( testCase(&SplineTest<2>::testValues));
        add( testCase(&SplineTest<3>::testValues));
        add( testCase(&SplineTest<4>::testValues));
        add( testCase(&SplineTest<5>::testValues));
        add( testCase(&SplineTest<0>::testFixedPointValues));
        add( testCase(&SplineTest<1>::testFixedPointValues));
        add( testCase(&SplineTest<2>::testFixedPointValues));
        add( testCase(&SplineTest<3>::testFixedPointValues));
        add( testCase(&SplineTest<0>::testPrefilterCoefficients));
        add( testCase(&SplineTest<1>::testPrefilterCoefficients));
        add( testCase(&SplineTest<2>::testPrefilterCoefficients));
        add( testCase(&SplineTest<3>::testPrefilterCoefficients));
        add( testCase(&SplineTest<4>::testPrefilterCoefficients));
        add( testCase(&SplineTest<5>::testPrefilterCoefficients));
        add( testCase(&SplineTest<0>::testWeightMatrix));
        add( testCase(&SplineTest<1>::testWeightMatrix));
        add( testCase(&SplineTest<2>::testWeightMatrix));
        add( testCase(&SplineTest<3>::testWeightMatrix));
        add( testCase(&SplineTest<4>::testWeightMatrix));
        add( testCase(&SplineTest<5>::testWeightMatrix));

        add( testCase(&FunctionsTest::testGaussians));
        add( testCase(&FunctionsTest::intSquareRootTest));
        add( testCase(&FunctionsTest::testSpecialFunctions));
        add( testCase(&FunctionsTest::closeAtToleranceTest));
        add( testCase(&FunctionsTest::testArgMinMax));

        add( testCase(&RationalTest::testGcdLcm));
        add( testCase(&RationalTest::testOStreamShifting));
        add( testCase(&RationalTest::testOperators));
        add( testCase(&RationalTest::testConversion));
        add( testCase(&RationalTest::testFunctions));
        add( testCase(&RationalTest::testInf));

        add( testCase(&LinalgTest::testOStreamShifting));
        add( testCase(&LinalgTest::testMatrix));
        add( testCase(&LinalgTest::testArgMinMax));
        add( testCase(&LinalgTest::testCholesky));
        add( testCase(&LinalgTest::testQR));
        add( testCase(&LinalgTest::testLinearSolve));
        add( testCase(&LinalgTest::testUnderdetermined));
        add( testCase(&LinalgTest::testOverdetermined));
        add( testCase(&LinalgTest::testIncrementalLinearSolve));
        add( testCase(&LinalgTest::testInverse));
        add( testCase(&LinalgTest::testSymmetricEigensystem));
        add( testCase(&LinalgTest::testNonsymmetricEigensystem));
        add( testCase(&LinalgTest::testDeterminant));
        add( testCase(&LinalgTest::testSVD));

        add( testCase(&FixedPointTest::testConstruction));
        add( testCase(&FixedPointTest::testComparison));
        add( testCase(&FixedPointTest::testArithmetic));
    }
};

int main()
{
  try 
  {
    MathTestSuite test;

    int failed = test.run();

    std::cerr << test.report() << std::endl;

    return (failed != 0);
  }
  catch(std::exception & e)
  {
    std::cerr << "Unexpected exception: " << e.what() << "\n";
    return 1;
  }
}
