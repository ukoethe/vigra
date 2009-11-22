/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#if defined(__GNUC__) && __GNUC__ == 2 && __GNUC_MINOR__ == 95
// deactivate broken std::relops
#  define __SGI_STL_INTERNAL_RELOPS
#endif  // __GNUC__

#include <typeinfo>
#include <iostream>
#include <iomanip>
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
#include "vigra/random.hxx"

#define VIGRA_TOLERANCE_MESSAGE "If this test fails, please adjust the tolerance threshold and report\n" \
                       "your findings (including compiler information etc.) to the VIGRA mailing list:"

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

    void testSpecialIntegerFunctions()
    {
        for(vigra::Int32 i = 0; i < 1024; ++i)
        {
            shouldEqual(vigra::sqrti(i), (vigra::Int32)vigra::floor(vigra::sqrt((double)i)));
        }

        shouldEqual(vigra::roundi(0.0), 0);
        shouldEqual(vigra::roundi(1.0), 1);
        shouldEqual(vigra::roundi(1.1), 1);
        shouldEqual(vigra::roundi(1.6), 2);
        shouldEqual(vigra::roundi(-1.0), -1);
        shouldEqual(vigra::roundi(-1.1), -1);
        shouldEqual(vigra::roundi(-1.6), -2);

        vigra::UInt32 roundPower2[] = {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 0xffff, 0x7fffffff, 0x80000000, 0x80000001, 0xffffffff};
        vigra::UInt32 floorResult[] = {0, 1, 2, 2, 4, 4, 4, 8, 8, 8, 16, 0x8000, 0x40000000, 0x80000000, 0x80000000, 0x80000000};
        vigra::UInt32 ceilResult[] = {0, 1, 2, 4, 4, 8, 8, 8, 16, 16, 16, 0x10000, 0x80000000, 0x80000000, 0, 0};
        for(int i = 0; i < sizeof(roundPower2) / sizeof(vigra::UInt32); ++i)
        {
            shouldEqual(vigra::floorPower2(roundPower2[i]), floorResult[i]);
            shouldEqual(vigra::ceilPower2(roundPower2[i]), ceilResult[i]);
        }

        for(vigra::Int32 k=0; k<32; ++k)
        {
            shouldEqual(vigra::log2i(1 << k), k);
            shouldEqual(vigra::log2i((1 << k) + 1), k == 0 ? 1 : k);
            shouldEqual(vigra::log2i((1 << k) - 1), k-1);
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

struct FixedPoint16Test
{
    void testConstruction()
    {
        shouldEqual((vigra::FixedPoint16<3>(3).value), 3 << 12);
        shouldEqual((vigra::FixedPoint16<3>(-3).value), -3 << 12);
        shouldEqual((-vigra::FixedPoint16<3>(3).value), -3 << 12);

        shouldEqual((vigra::FixedPoint16<8>(3).value), 3 << 7);

        shouldEqual((vigra::FixedPoint16<3>(3.5).value), 7 << 11);
        shouldEqual((vigra::FixedPoint16<3>(-3.5).value), -(7 << 11));
        shouldEqual((-vigra::FixedPoint16<3>(3.5).value), -(7 << 11));

        shouldEqual((vigra::NumericTraits<vigra::FixedPoint16<4> >::zero()).value, 0);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint16<4> >::one()).value, 1 << 11);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint16<4> >::max()).value, (1 << 15) - 1);
        shouldEqual((vigra::NumericTraits<vigra::FixedPoint16<4> >::min()).value, -(1 << 15));

        shouldEqual((vigra::FixedPoint16<1, vigra::FPOverflowSaturate>(3.75).value), (1 << 15)-1);
        shouldEqual((vigra::FixedPoint16<1, vigra::FPOverflowSaturate>(-3.75).value), -(1 << 15));
        try { vigra::FixedPoint16<1, vigra::FPOverflowError>(3.75); failTest("No exception thrown"); } 
        catch(vigra::PreconditionViolation &) {}
        try { vigra::FixedPoint16<1, vigra::FPOverflowError>(-3.75); failTest("No exception thrown"); } 
        catch(vigra::PreconditionViolation &) {}

        vigra::FixedPoint16<4> v(3.75);
        shouldEqual((v.value), 15 << 9);
        shouldEqual(((-v).value), -15 << 9);
        shouldEqual((vigra::FixedPoint16<4>(v).value), 15 << 9);
        shouldEqual((vigra::FixedPoint16<6>(v).value), 15 << 7);
        shouldEqual((vigra::FixedPoint16<13>(v).value), 15);
        shouldEqual((vigra::FixedPoint16<15>(v).value), 4);

        shouldEqual((vigra::FixedPoint16<4>(-v).value), -15 << 9);
        shouldEqual((vigra::FixedPoint16<6>(-v).value), -15 << 7);
        shouldEqual((vigra::FixedPoint16<13>(-v).value), -15);
        shouldEqual((vigra::FixedPoint16<15>(-v).value), -4);

        shouldEqual(vigra::fixed_point_cast<double>(v), 3.75);
        shouldEqual(vigra::fixed_point_cast<double>(-v), -3.75);
        shouldEqual(frac(v), vigra::FixedPoint16<4>(0.75));
        shouldEqual(dual_frac(v), vigra::FixedPoint16<4>(0.25));
        shouldEqual(frac(-v), vigra::FixedPoint16<4>(0.25));
        shouldEqual(dual_frac(-v), vigra::FixedPoint16<4>(0.75));
        shouldEqual(vigra::floor(v), 3);
        shouldEqual(vigra::ceil(v), 4);
        shouldEqual(vigra::floor(-v), -4);
        shouldEqual(vigra::ceil(-v), -3);
        shouldEqual(round(v), 4);
        shouldEqual(round(-v), -4);
        shouldEqual(vigra::abs(v), v);
        shouldEqual(vigra::abs(-v), v);
        shouldEqual(vigra::norm(-v), v);
        shouldEqual(vigra::squaredNorm(-v), v*v);

        vigra::FixedPoint16<2> v1;
        shouldEqual((v1 = v).value, 15 << 11);
        shouldEqual((v1 = -v).value, -15 << 11);

        vigra::FixedPoint16<15> v2;
        shouldEqual((v2 = v).value, 4);
        shouldEqual((v2 = -v).value, -4);
    }

    void testComparison()
    {
        vigra::FixedPoint16<4> v1(3.75), v2(4);
        vigra::FixedPoint16<2> v3(3.75);
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
        typedef vigra::FixedPoint16<1> FP1;
        typedef vigra::FixedPoint16<2> FP2;
        typedef vigra::FixedPoint16<7> FP7;
        typedef vigra::FixedPoint16<8> FP8;
        typedef vigra::FixedPoint16<13> FP13;
        typedef vigra::FixedPoint16<15> FP15;
        
        FP1 t0(0), t1(0.75), t2(0.25);
        signed char v1 = 1, v2 = 2, v4 = 4, v8 = 8;

        shouldEqual(FP1(t1) += t1, FP1(1.5));
        shouldEqual(FP1(t1) -= t1, FP1(0.0));
        shouldEqual(FP1(t1) -= t2, FP1(0.5));
        shouldEqual(FP2(t1) *= t1, FP1(9.0 / 16.0));
        shouldEqual(FP2(t1) /= t2, FP2(3));
        shouldEqual(FP2(t1) /= t0, vigra::NumericTraits<FP2>::max());
        shouldEqual(FP2(-t1) /= t0, vigra::NumericTraits<FP2>::min());
        
        FP2 res;
        shouldEqual(add(t1, t1, res), FP2(1.5));
        shouldEqual(sub(t1, t1, res), FP2(0));
        shouldEqual(sub(t1, t2, res), FP2(0.5));
        shouldEqual(mul(t1, t1, res), FP2(9.0 / 16.0));
        shouldEqual(div(t1, t2, res), FP2(3));
        shouldEqual(div(t1, t0, res), vigra::NumericTraits<FP2>::max());
        shouldEqual(div(-t1, t0, res), vigra::NumericTraits<FP2>::min());

        shouldEqual(--t1, FP1(-0.25));
        shouldEqual(t1, FP1(-0.25));
        shouldEqual(++t1, FP1(0.75));
        shouldEqual(t1, FP1(0.75));
        shouldEqual(t1++, FP1(0.75));
        shouldEqual(t1, FP1(1.75));
        shouldEqual(t1--, FP1(1.75));
        shouldEqual(t1, FP1(0.75));

        shouldEqual((t1 * FP7(v1)).value, 3 << 6);
        shouldEqual((t2 * FP7(v1)).value, 1 << 6);
        shouldEqual((-t1 * FP7(v1)).value, -3 << 6);
        shouldEqual((-t2 * FP7(v1)).value, -1 << 6);
        shouldEqual((t1 * -FP7(v1)).value, -3 << 6);
        shouldEqual((t2 * -FP7(v1)).value, -1 << 6);

        shouldEqual((vigra::FixedPoint16<2, vigra::FPOverflowSaturate>(t1*FP7(v8)).value), (1 << 15)-1);
        shouldEqual((vigra::FixedPoint16<2, vigra::FPOverflowSaturate>(t1*FP7(-v8)).value), -(1 << 15));
        try { vigra::FixedPoint16<2, vigra::FPOverflowError>(t1*FP7(v8)); failTest("No exception thrown"); } 
        catch(vigra::PreconditionViolation &) {}
        try { vigra::FixedPoint16<2, vigra::FPOverflowError>(t1*FP7(-v8)); failTest("No exception thrown"); } 
        catch(vigra::PreconditionViolation &) {}

        shouldEqual((FP13(t1 * FP7(v1))).value, 3);
        shouldEqual((FP13(t2 * FP7(v1))).value, 1);
        shouldEqual((FP13(-t1 * FP7(v1))).value, -3);
        shouldEqual((FP13(-t2 * FP7(v1))).value, -1);

        shouldEqual((t1 * FP7(v4) + t2 * FP7(v8)).value, 5 << 8);
        shouldEqual((t1 * FP7(v1) + t2 * FP7(v2)).value, 5 << 6);

        shouldEqual(FP7(6) / FP7(3), FP7(2));
        shouldEqual(FP7(0.75) / FP7(0.25), FP7(3));
        shouldEqual(FP7(12) / FP7(48), FP7(0.25));
        shouldEqual(FP1(0.25) / FP7(2), FP7(0.125));
        shouldEqual(FP7(10) / FP1(0.25), FP7(40));
        shouldEqual(FP7(10) / t0, vigra::NumericTraits<FP7>::max());
        shouldEqual(FP7(-10) / t0, vigra::NumericTraits<FP7>::min());

        shouldEqual(vigra::floor(t1 * FP7(v1) + t2 * FP7(v2)), 1);
        shouldEqual(vigra::ceil(t1 * FP7(v1) + t2 * FP7(v2)), 2);
        shouldEqual(round(t1 * FP7(v1) + t2 * FP7(v2)), 1);
        shouldEqual(vigra::floor(t1 * FP7(v4) + t2 * FP7(v8)), 5);
        shouldEqual(vigra::ceil(t1 * FP7(v4) + t2 * FP7(v8)), 5);
        shouldEqual(round(t1 * FP7(v4) + t2 * FP7(v8)), 5);

        shouldEqual(vigra::floor(t1 * -FP7(v1) - t2 * FP7(v2)), -2);
        shouldEqual(vigra::ceil(t1 * -FP7(v1) - t2 * FP7(v2)), -1);
        shouldEqual(round(t1 * -FP7(v1) - t2 * FP7(v2)), -1);
        shouldEqual(vigra::floor(t1 * -FP7(v4) - t2 * FP7(v8)), -5);
        shouldEqual(vigra::ceil(t1 * -FP7(v4) - t2 * FP7(v8)), -5);
        shouldEqual(round(t1 * -FP7(v4) - t2 * FP7(v8)), -5);

        double d1 = 1.0 / 3.0, d2 = 1.0 / 7.0;
        FP1 r1(d1), r2(d2);
        FP2 r3;
        add(r1, r2, r3);
        shouldEqual(r3.value, FP2(d1 + d2).value);
        sub(r1, r2, r3);
        shouldEqual(r3.value, FP2(d1 - d2).value);
        mul(r1, r2, r3);
        shouldEqual(r3.value >> 2, FP2(d1 * d2).value >> 2);

        shouldEqual(vigra::sqrt(FP7(4)).value, 1 << 12);
        shouldEqual(vigra::sqrt(FP8(4)).value, 1 << 12);
        shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(FP7(4))), 2.0);
        shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(FP2(2.25))), 1.5);
        shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(FP8(6.25))), 2.5);

        for(int i = 0; i < 1024; ++i)
        {
            vigra::FixedPoint16<11> fv1(i, vigra::FPNoShift);
            vigra::FixedPoint16<10> fv2(i, vigra::FPNoShift);
            vigra::FixedPoint16<9>  fv3(i, vigra::FPNoShift);
            vigra::FixedPoint16<8>  fv4(i, vigra::FPNoShift);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv1)), vigra::floor(vigra::sqrt((double)(i << 14))) / 512.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv2)), vigra::floor(vigra::sqrt((double)(i << 15))) / 1024.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv3)), vigra::floor(vigra::sqrt((double)(i << 14))) / 1024.0);
            shouldEqual(vigra::fixed_point_cast<double>(vigra::sqrt(fv4)), vigra::floor(vigra::sqrt((double)(i << 15))) / 2048.0);
        }

        shouldEqual(vigra::atan2(FP1(0), FP1(1)), FP2(0));
        shouldEqual(vigra::atan2(FP1(0), FP1(-1)), FP2(M_PI));
        shouldEqual(vigra::atan2(FP1(1), FP1(0)), FP2(0.5*M_PI));
        shouldEqual(vigra::atan2(FP1(-1), FP1(0)), FP2(-0.5*M_PI));
        
        for(int i = -179; i < 180; ++i)
        {
            double angle = M_PI*i/180.0;
            double c = std::cos(angle), s = std::sin(angle);
            FP2 a = vigra::atan2(FP1(s), FP1(c));
            should(vigra::abs(i-vigra::fixed_point_cast<double>(a)/M_PI*180.0) < 0.5);
            a = vigra::atan2(FP15(30000.0*s), FP15(30000.0*c));
            should(vigra::abs(i-vigra::fixed_point_cast<double>(a)/M_PI*180.0) < 0.5);
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
        std::string sref("1.0000 5.0000 \n3.0000 2.0000 \n4.0000 7.0000 \n");
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

        Matrix one = zero + Matrix(r,c).init(1.0);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(one(i,j), 1.0);

        std::stringstream s;
        s << std::setprecision(4) << a;
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
        b += 4.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 7.0+data[k]);
        b += a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 7.0+2.0*data[k]);


        b = 4.0 - a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0-data[k]);
        b = a - 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]-3.0);
        b -= 4.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]-7.0);
        b -= a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), -7.0);

        b = 4.0 * a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0*data[k]);
        b = a * 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]*3.0);
        b *= 4.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]*12.0);
        b *= a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k]*data[k]*12.0);

        b = 4.0 / a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), 4.0/data[k]);
        b = a / 3.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k] / 3.0);
        b /= 4.0;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqual(b(i,j), data[k] / 12.0);
        b /= a;
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<c; ++j, ++k)
                shouldEqualTolerance(b(i,j), 1.0 / 12.0, 1e-12);

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
            shouldEqual(ac, subVector(ac, 0, r));
            shouldEqual(a.subarray(Shape(1, i), Shape(r-1, i+1)), subVector(ac, 1, r-1));
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
        
        shouldEqual(trace(a2), 3.0);

        Matrix id = vigra::identityMatrix<double>(r);
        shouldEqual(id.rowCount(), r);
        shouldEqual(id.columnCount(), r);
        for(unsigned int i=0, k=0; i<r; ++i)
            for(unsigned int j=0; j<r; ++j, ++k)
                shouldEqual(id(i,j), idref[k]);

        shouldEqual(trace(id), 3.0);

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
        for(int i=0; i<g.rowCount(); ++i)
            for(int j=0; j<g.columnCount(); ++j)
                shouldEqual(g(i,j), data[i]*tref[j]);

        Matrix h = transpose(a) * a;
        shouldEqual(h.rowCount(), c);
        shouldEqual(h.columnCount(), c);
        for(int i=0; i<(int)c; ++i)
            for(int j=0; j<(int)c; ++j)
                shouldEqual(h(i,j), dot(rowVector(at, i), columnVector(a, j)));

        should(isSymmetric(random_symmetric_matrix(10)));
        should(!isSymmetric(random_matrix(10, 10)));

        Matrix tm(2, 2, tref2);
        vigra::TinyVector<double, 2> tv(1.0, 2.0), tvrref(7.0, 9.0), tvlref(11.0, 7.0);
        shouldEqual(tm * tv, tvrref);
        shouldEqual(tv * tm, tvlref);

        Matrix rep = repeatMatrix(a, 2, 4);
        shouldEqual(rowCount(rep), 2*r);
        shouldEqual(columnCount(rep), 4*c);

        for(unsigned int l=0; l<4; ++l)
            for(unsigned int k=0; k<2; ++k)
                for(unsigned int j=0; j<c; ++j)
                    for(unsigned int i=0; i<r; ++i)
                        shouldEqual(rep(k*r+i, l*c+j), a(i,j));
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

    void testColumnAndRowStatistics()
    {
#if defined(__GNUC__)
        double epsilon = 1e-11, epsilon2 = 1e-8;
#else
        double epsilon = 1e-13, epsilon2 = 1e-10;
#endif

        Matrix rowMean(size, 1), columnMean(1, size);
        Matrix rowStdDev(size, 1), columnStdDev(1, size);
        Matrix rowNorm(size, 1), columnNorm(1, size);
        Matrix rowCovariance(size, size), columnCovariance(size, size);

        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_matrix (size, size);

            rowStatistics(a, rowMean, rowStdDev, rowNorm);
            columnStatistics(a, columnMean, columnStdDev, columnNorm);

            for(unsigned int k=0; k<size; ++k)
            {
                double rm = 0.0, cm = 0.0, rn = 0.0, cn = 0.0, rs = 0.0, cs = 0.0;
                for(unsigned int l=0; l<size; ++l)
                {
                    rm += a(k, l);
                    cm += a(l, k);
                    rn += vigra::sq(a(k, l));
                    cn += vigra::sq(a(l, k));
                }
                rm /= size;
                cm /= size;
                rn = std::sqrt(rn);
                cn = std::sqrt(cn);

                shouldEqualTolerance(rm, rowMean(k,0), epsilon);
                shouldEqualTolerance(cm, columnMean(0,k), epsilon);
                shouldEqualTolerance(rn, rowNorm(k,0), epsilon);
                shouldEqualTolerance(cn, columnNorm(0,k), epsilon);

                for(unsigned int l=0; l<size; ++l)
                {
                    rs += vigra::sq(a(k, l) - rm);
                    cs += vigra::sq(a(l, k) - cm);
                }
                rs = std::sqrt(rs / (size-1));
                cs = std::sqrt(cs / (size-1));

                shouldEqualTolerance(rs, rowStdDev(k,0), epsilon);
                shouldEqualTolerance(cs, columnStdDev(0,k), epsilon);
            }

            covarianceMatrixOfRows(a, rowCovariance);
            covarianceMatrixOfColumns(a, columnCovariance);
            Matrix rowCovarianceRef(size, size), columnCovarianceRef(size, size);
            for(unsigned int k=0; k<size; ++k)
            {
                for(unsigned int l=0; l<size; ++l)
                {
                    for(unsigned int m=0; m<size; ++m)
                    {
                        rowCovarianceRef(l, m) += (a(l, k) - rowMean(l, 0)) * (a(m, k) - rowMean(m, 0));
                        columnCovarianceRef(l, m) += (a(k, l) - columnMean(0, l)) * (a(k, m) - columnMean(0, m));
                    }
                }
            }
            rowCovarianceRef /= (size-1);
            columnCovarianceRef /= (size-1);

            shouldEqualSequenceTolerance(rowCovariance.data(), rowCovariance.data()+size*size, rowCovarianceRef.data(), epsilon2);
            shouldEqualSequenceTolerance(columnCovariance.data(), columnCovariance.data()+size*size, columnCovarianceRef.data(), epsilon2);
        }
    }

    void testColumnAndRowPreparation()
    {
        using vigra::ZeroMean;
        using vigra::UnitVariance;
        using vigra::UnitNorm;

#if defined(__GNUC__) && __GNUC__ == 3
        double epsilon = 1e-11, epsilon2 = 1e-8;
#else
        double epsilon = 1e-13, epsilon2 = 1e-10;
#endif

        Matrix rowMean(size, 1), columnMean(1, size);
        Matrix rowStdDev(size, 1), columnStdDev(1, size);
        Matrix rowNorm(size, 1), columnNorm(1, size);

        Matrix rowPrepared(size, size), columnPrepared(size, size);
        Matrix rowMeanPrepared(size, 1), columnMeanPrepared(1, size);
        Matrix rowStdDevPrepared(size, 1), columnStdDevPrepared(1, size);
        Matrix rowNormPrepared(size, 1), columnNormPrepared(1, size);
        Matrix rowOffset(size, 1), columnOffset(1, size);
        Matrix rowScaling(size, 1), columnScaling(1, size);

        Matrix zeroRowRef(size,1), zeroColRef(1, size);
        Matrix oneRowRef(size,1), oneColRef(1, size);
        oneRowRef.init(1.0);
        oneColRef.init(1.0);

        {
            Matrix a = random_matrix (size, size);

            columnStatistics(a, columnMean, columnStdDev, columnNorm);

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean);
            columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(columnStdDev.data(), columnStdDev.data()+size, columnStdDevPrepared.data(), epsilon);

            Matrix ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, UnitNorm);
            columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnNormPrepared.data(), epsilon);

            ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, UnitVariance);
            columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            columnMeanPrepared /= columnScaling;
            shouldEqualSequenceTolerance(columnMean.data(), columnMean.data()+size, columnMeanPrepared.data(), epsilon2);
            shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnStdDevPrepared.data(), epsilon);

            ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitVariance);
            columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnStdDevPrepared.data(), epsilon);

            ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitNorm);
            columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnNormPrepared.data(), epsilon);

            ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon2);

            rowStatistics(a, rowMean, rowStdDev, rowNorm);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean);
            rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(rowStdDev.data(), rowStdDev.data()+size, rowStdDevPrepared.data(), epsilon);

            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, UnitNorm);
            rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowNormPrepared.data(), epsilon);

            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, UnitVariance);
            rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            rowMeanPrepared /= rowScaling;
            shouldEqualSequenceTolerance(rowMean.data(), rowMean.data()+size, rowMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowStdDevPrepared.data(), epsilon);

            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitVariance);
            rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowStdDevPrepared.data(), epsilon);

            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitNorm);
            rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowNormPrepared.data(), epsilon);

            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon2);
        }

        {
            Matrix a(size, size, 2.0), aref(size, size, 1.0/std::sqrt((double)size));

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitVariance);
            shouldEqualSequence(a.data(), a.data()+size*size, columnPrepared.data());

            prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitNorm);
            shouldEqualSequenceTolerance(aref.data(), aref.data()+size*size, columnPrepared.data(), epsilon);
            Matrix ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitVariance);
            shouldEqualSequence(a.data(), a.data()+size*size, rowPrepared.data());

            prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitNorm);
            shouldEqualSequenceTolerance(aref.data(), aref.data()+size*size, rowPrepared.data(), epsilon);
            ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);
        }
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
#if defined(__GNUC__)
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

            should(linearSolve (a, b, x, "QR"));
            Matrix ax = a * x;
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

        {
            // square, rank-deficient system (compute minimum norm solution)
            double mdata[] = {1.0,  3.0,  7.0,
                             -1.0,  4.0,  4.0,
                              1.0, 10.0, 18.0};
            double rhsdata[] = { 5.0, 2.0, 12.0};
            double refdata[] = { 0.3850, -0.1103, 0.7066 };

            Matrix m(3,3,mdata), rhs(3,1,rhsdata), xx(3,1);

            shouldEqual(linearSolveQR(m, rhs, xx), 2);
            shouldEqualSequenceTolerance(refdata, refdata+3, xx.data(), 1e-3);
        }
        {
            // underdetermined, full-rank system (compute minimum norm solution)
            double mdata[] = {2.0, -3.0, 1.0, -6.0,
                              4.0,  1.0, 2.0,  9.0,
                              3.0,  1.0, 1.0,  8.0};
            double rhsdata[] = { -7.0, -7.0, -8.0};
            double refdata[] = { -3.26666666666667, 3.6, 5.13333333333333, -0.86666666666667 };

            Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            shouldEqual(linearSolveQR(m, rhs, xx), 3);
            shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-12);
        }
        {
            // underdetermined, rank-deficient, consistent system (compute minimum norm solution)
            double mdata[] = {1.0,  3.0, 3.0, 2.0,
                              2.0,  6.0, 9.0, 5.0,
                             -1.0, -3.0, 3.0, 0.0};
            double rhsdata[] = { 1.0, 5.0, 5.0};
            double refdata[] = { -0.211009, -0.633027, 0.963303, 0.110092 };

            Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            shouldEqual(linearSolveQR(m, rhs, xx), 2);
            shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-5);
        }
        {
            // underdetermined, rank-deficient, inconsistent system (compute minimum norm least squares solution)
            double mdata[] = {2.0, 1.0,  7.0, -7.0,
                             -3.0, 4.0, -5.0, -6.0,
                              1.0, 1.0,  4.0, -5.0};
            double rhsdata[] = { 2.0, 3.0, 2.0};
            double refdata[] = { -0.0627, 0.1561, -0.0321, -0.3427 };

            Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            shouldEqual(linearSolveQR(m, rhs, xx), 2);
            shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-3);
        }
    }

    void testOverdetermined()
    {
#if defined(__GNUC__)
        double epsilon = 1e-11;
#else
        double epsilon = 1e-12;
#endif

        unsigned int n = 5;
        unsigned int size = 1000;
        double noiseStdDev = 0.1;

        Matrix A(size, n), xs(n,1), xq(n,1), xn(n,1), r(size, 1);

        for(unsigned int iter=0; iter<iterations; ++iter)
        {
            // set up a linear regression problem for a polynomial of degree n
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

            // check that solution is indeed a minimum by
            // testing for zero derivative of the objective
            Matrix derivative = abs(transpose(A)*(A*xs - rhs));
            int absIndex = argMax(derivative);
            shouldEqualTolerance(derivative(absIndex,0), 0.0, epsilon);

            shouldEqual(linearSolveQR(A, rhs, xq), n);
            shouldEqualSequenceTolerance(xs.data(), xs.data()+n, xq.data(), epsilon);

            shouldEqual(linearSolve(A, rhs, xn, "ne"), true);
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

    void testSymmetricEigensystemAnalytic()
    {
        double epsilon = 1e-8;

        int size = 2;
        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_symmetric_matrix (size);
            Matrix ew(size, 1), ewref(size, 1);
            Matrix ev(size, size);
            symmetricEigensystem(a, ewref, ev);
            vigra::symmetric2x2Eigenvalues(
                a(0,0), a(0,1),
                a(1,1),
                &ew(0,0), &ew(1,0));
            shouldEqualSequenceTolerance(ew.data(), ew.data()+size, ewref.data(), epsilon);
        }

        size = 3;
        for(unsigned int i = 0; i < iterations; ++i)
        {
            Matrix a = random_symmetric_matrix (size);
            Matrix ew(size, 1), ewref(size, 1);
            Matrix ev(size, size);
            symmetricEigensystem(a, ewref, ev);
            vigra::symmetric3x3Eigenvalues<double>(
                a(0,0), a(0,1), a(0,2),
                a(1,1), a(1,2),
                a(2,2),
                &ew(0,0), &ew(1,0), &ew(2,0));
            shouldEqualSequenceTolerance(ew.data(), ew.data()+size, ewref.data(), epsilon);
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

#if defined(__GNUC__)
        double eps = 1e-9;
#else
        double eps = 1e-10;
#endif
   	    shouldEqualToleranceMessage(norm(a-u*diagonalMatrix(S)*transpose(v)), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
	    shouldEqualToleranceMessage(norm(vigra::identityMatrix<double>(4) - transpose(u)*u), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
	    shouldEqualToleranceMessage(norm(vigra::identityMatrix<double>(4) - transpose(v)*v), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
	    shouldEqualToleranceMessage(norm(vigra::identityMatrix<double>(4) - v*transpose(v)), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
    }
};

struct RandomTest
{
    void testTT800()
    {
        const unsigned int n = 50;
        unsigned int iref[n] = {
            3169973338U, 2724982910U,  347012937U, 1735893326U, 2282497071U,
            3975116866U,   62755666U,  500522132U,  129776071U, 1978109378U,
            4040131704U, 3800592193U, 3057303977U, 1468369496U,  370579849U,
            3630178833U,   51910867U,  819270944U,  476180518U,  190380673U,
            1370447020U, 1620916304U,  663482756U, 1354889312U, 4000276916U,
             868393086U, 1441698743U, 1086138563U, 1899869374U, 3717419747U,
            2455034041U, 2617437696U, 1595651084U, 4148285605U, 1860328467U,
             928897371U,  263340857U, 4091726170U, 2359987311U, 1669697327U,
            1882626857U, 1635656338U,  897501559U, 3233276032U,  373770970U,
            2950632840U, 2706386845U, 3294066568U, 3819538748U, 1902519841U };

        vigra::RandomTT800 random;
        for(int k=0; k<n; ++k)
            shouldEqual(random(), iref[k]);

        double fref[n] = {
              0.738067,   0.634460,   0.080795,   0.404169,   0.531435,
              0.925529,   0.014611,   0.116537,   0.030216,   0.460564,
              0.940666,   0.884894,   0.711834,   0.341881,   0.086282,
              0.845217,   0.012086,   0.190751,   0.110869,   0.044326,
              0.319082,   0.377399,   0.154479,   0.315460,   0.931387,
              0.202189,   0.335672,   0.252886,   0.442348,   0.865529,
              0.571607,   0.609420,   0.371516,   0.965848,   0.433141,
              0.216276,   0.061314,   0.952679,   0.549477,   0.388757,
              0.438333,   0.380831,   0.208966,   0.752806,   0.087025,
              0.686998,   0.630130,   0.766960,   0.889306,   0.442965 };
        vigra::RandomTT800 randomf;
        for(int k=0; k<n; ++k)
            should(vigra::abs(randomf.uniform() - fref[k]) < 2e-6);

        vigra::RandomTT800 randomr(vigra::RandomSeed);
    }

    void testMT19937()
    {
        const unsigned int n = 20, skip = 960, ilen = 4;
        unsigned int first[n] = {
             956529277U, 3842322136U, 3319553134U, 1843186657U, 2704993644U,
             595827513U,  938518626U, 1676224337U, 3221315650U, 1819026461U,
            2401778706U, 2494028885U,  767405145U, 1590064561U, 2766888951U,
            3951114980U, 2568046436U, 2550998890U, 2642089177U,  568249289U };
        unsigned int last[n] = {
            2396869032U, 1982500200U, 2649478910U,  839934727U, 3814542520U,
             918389387U,  995030736U, 2017568170U, 2621335422U, 1020082601U,
              24244213U, 2575242697U, 3941971804U,  922591409U, 2851763435U,
            2055641408U, 3695291669U, 2040276077U, 4118847636U, 3528766079U };

        vigra::RandomMT19937 random(0xDEADBEEF);
        for(int k=0; k<n; ++k)
            shouldEqual(random(), first[k]);
        for(int k=0; k<skip; ++k)
            random();
        for(int k=0; k<n; ++k)
            shouldEqual(random(), last[k]);

        for(int k=0; k<skip; ++k)
            should(random.uniformInt(31) < 31);

        random.seed(0xDEADBEEF);
        for(int k=0; k<n; ++k)
            shouldEqual(random(), first[k]);
        for(int k=0; k<skip; ++k)
            random();
        for(int k=0; k<n; ++k)
            shouldEqual(random(), last[k]);

        unsigned int firsta[n] = {
            1067595299U,  955945823U,  477289528U, 4107218783U, 4228976476U,
            3344332714U, 3355579695U,  227628506U,  810200273U, 2591290167U,
            2560260675U, 3242736208U,  646746669U, 1479517882U, 4245472273U,
            1143372638U, 3863670494U, 3221021970U, 1773610557U, 1138697238U };
        unsigned int lasta[n] = {
             123599888U,  472658308U, 1053598179U, 1012713758U, 3481064843U,
            3759461013U, 3981457956U, 3830587662U, 1877191791U, 3650996736U,
             988064871U, 3515461600U, 4089077232U, 2225147448U, 1249609188U,
            2643151863U, 3896204135U, 2416995901U, 1397735321U, 3460025646U };

        unsigned int init[ilen] = {0x123, 0x234, 0x345, 0x456};
        vigra::RandomMT19937 randoma(init, ilen);
        for(int k=0; k<n; ++k)
            shouldEqual(randoma(), firsta[k]);
        for(int k=0; k<skip; ++k)
            randoma();
        for(int k=0; k<n; ++k)
            shouldEqual(randoma(), lasta[k]);

        double ref53[n] = {
            0.76275444, 0.98670464, 0.27933125, 0.94218739, 0.78842173,
            0.92179002, 0.54534773, 0.38107717, 0.65286910, 0.22765212,
            0.74557914, 0.54708246, 0.42043117, 0.19189126, 0.70259889,
            0.77408120, 0.04605807, 0.69398269, 0.61711170, 0.10133577};
        for(int k=0; k<n; ++k)
            should(vigra::abs(randoma.uniform53()-ref53[k]) < 2e-8);

        randoma.seed(init, ilen);
        for(int k=0; k<n; ++k)
            shouldEqual(randoma(), firsta[k]);
        for(int k=0; k<skip; ++k)
            randoma();
        for(int k=0; k<n; ++k)
            shouldEqual(randoma(), lasta[k]);
    }

    void testRandomFunctors()
    {
        const unsigned int n = 50;
        unsigned int iref[n] = {
            3169973338U, 2724982910U,  347012937U, 1735893326U, 2282497071U,
            3975116866U,   62755666U,  500522132U,  129776071U, 1978109378U,
            4040131704U, 3800592193U, 3057303977U, 1468369496U,  370579849U,
            3630178833U,   51910867U,  819270944U,  476180518U,  190380673U,
            1370447020U, 1620916304U,  663482756U, 1354889312U, 4000276916U,
             868393086U, 1441698743U, 1086138563U, 1899869374U, 3717419747U,
            2455034041U, 2617437696U, 1595651084U, 4148285605U, 1860328467U,
             928897371U,  263340857U, 4091726170U, 2359987311U, 1669697327U,
            1882626857U, 1635656338U,  897501559U, 3233276032U,  373770970U,
            2950632840U, 2706386845U, 3294066568U, 3819538748U, 1902519841U };
        double fref[n] = {
              0.738067,   0.634460,   0.080795,   0.404169,   0.531435,
              0.925529,   0.014611,   0.116537,   0.030216,   0.460564,
              0.940666,   0.884894,   0.711834,   0.341881,   0.086282,
              0.845217,   0.012086,   0.190751,   0.110869,   0.044326,
              0.319082,   0.377399,   0.154479,   0.315460,   0.931387,
              0.202189,   0.335672,   0.252886,   0.442348,   0.865529,
              0.571607,   0.609420,   0.371516,   0.965848,   0.433141,
              0.216276,   0.061314,   0.952679,   0.549477,   0.388757,
              0.438333,   0.380831,   0.208966,   0.752806,   0.087025,
              0.686998,   0.630130,   0.766960,   0.889306,   0.442965 };
        double nref[n] = {
            1.35298, 0.764158, -0.757076, -0.173069, 0.0586711,
            0.794212, -0.483372, -0.0405762, 1.27956, -0.955101,
            -1.5062, -1.02069, -0.871562, -0.465495, -0.799888,
            -1.20286, -0.170944, 1.08383, 1.26832, 1.93807,
            -0.098183, 0.355986, -0.336965, -1.42996, 0.966012,
            -2.17195, -1.05422, -2.03724, -0.769992, 0.668851,
            -0.570259, 0.258217, 0.632492, 1.29755, 0.96869,
            -0.141918, -0.836236, -0.62337, 0.116509, -0.0314471,
            0.402451, -1.20504, -0.140861, -0.0765263, 1.06057,
            2.57671, 0.0299117, 0.471425, 1.59464, 1.37346};

        vigra::RandomTT800 random1;
        vigra::UniformRandomFunctor<> f1(random1);
        for(int k=0; k<n; ++k)
            should(vigra::abs(f1() - fref[k]) < 2e-6);

        vigra::RandomTT800 random2;
        vigra::UniformIntRandomFunctor<> f2(4, 34, random2, true);
        for(int k=0; k<n; ++k)
            shouldEqual(f2(), iref[k] % 31 + 4);

        vigra::RandomTT800 random3;
        vigra::UniformIntRandomFunctor<> f3(random3);
        for(int k=0; k<n; ++k)
            shouldEqual(f3(32), iref[k] % 32);

        vigra::RandomTT800 random4;
        vigra::NormalRandomFunctor<> f4(random4);
        for(int k=0; k<n; ++k)
            shouldEqualTolerance(f4(), nref[k], 1e-5);
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
        add( testCase(&FunctionsTest::testSpecialIntegerFunctions));
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
        add( testCase(&LinalgTest::testColumnAndRowStatistics));
        add( testCase(&LinalgTest::testColumnAndRowPreparation));
        add( testCase(&LinalgTest::testCholesky));
        add( testCase(&LinalgTest::testQR));
        add( testCase(&LinalgTest::testLinearSolve));
        add( testCase(&LinalgTest::testUnderdetermined));
        add( testCase(&LinalgTest::testOverdetermined));
        add( testCase(&LinalgTest::testIncrementalLinearSolve));
        add( testCase(&LinalgTest::testInverse));
        add( testCase(&LinalgTest::testSymmetricEigensystem));
        add( testCase(&LinalgTest::testNonsymmetricEigensystem));
        add( testCase(&LinalgTest::testSymmetricEigensystemAnalytic));
        add( testCase(&LinalgTest::testDeterminant));
        add( testCase(&LinalgTest::testSVD));

        add( testCase(&FixedPointTest::testConstruction));
        add( testCase(&FixedPointTest::testComparison));
        add( testCase(&FixedPointTest::testArithmetic));
        add( testCase(&FixedPoint16Test::testConstruction));
        add( testCase(&FixedPoint16Test::testComparison));
        add( testCase(&FixedPoint16Test::testArithmetic));

        add( testCase(&RandomTest::testTT800));
        add( testCase(&RandomTest::testMT19937));
        add( testCase(&RandomTest::testRandomFunctors));
    }
};

int main(int argc, char ** argv)
{
  try
  {
    MathTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cerr << test.report() << std::endl;

    return (failed != 0);
  }
  catch(std::exception & e)
  {
    std::cerr << "Unexpected exception: " << e.what() << "\n";
    return 1;
  }
}
