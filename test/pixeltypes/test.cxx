#include <iostream>
#include "unittest.hxx"

#include "vigra/accessor.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/rgbvalue.hxx"

using namespace vigra;

static float di[] = {1, 2, 4 };
static float df[] = {1.2, 2.4, 3.6 };

template <class BVector, class IVector, class FVector>
struct TinyVectorTest
{
    typedef BVector BV;
    typedef IVector IV;
    typedef FVector FV;

    BV bv0, bv1, bv3;
    IV iv0, iv1, iv3;
    FV fv0, fv1, fv3;

    template <class VECTOR>
    void printVector(VECTOR const & v)
    {
        std::cerr << "(";
        for(int i=0; i<v.size(); ++i)
            std::cerr << (float)v[i] << ", ";
        std::cerr << ")\n";
    }

    template <class VECTOR, class VALUE>
    bool equalValue(VECTOR const & v, VALUE const & vv)
    {
        for(int i=0; i<v.size(); ++i)
            if(v[i] != vv)
                return false;
        return true;
    }

    template <class VECTOR1, class VECTOR2>
    bool equalVector(VECTOR1 const & v1, VECTOR2 const & v2)
    {
        for(int i=0; i<v1.size(); ++i)
            if(v1[i] != v2[i])
                return false;
        return true;
    }

    template <class ITER1, class ITER2>
    bool equalIter(ITER1 i1, ITER1 i1end, ITER2 i2)
    {
        if(i1end - i1 != 3)
            return false;
        for(; i1<i1end; ++i1, ++i2)
            if(*i1 != *i2)
                return false;
        return true;
    }

    TinyVectorTest()
    : bv0(0), bv1(1), bv3(),
      iv0(0), iv1(1), iv3(),
      fv0(0.0), fv1(1.0), fv3()
    {
        bv3.init(df, df+3);
        iv3.init(df, df+3);
        fv3.init(df, df+3);
    }

    void testConstruction()
    {
        should(bv0.size() == 3);
        should(iv0.size() == 3);
        should(fv0.size() == 3);

        should(equalValue(bv0, 0));
        should(equalValue(iv0, 0));
        should(equalValue(fv0, 0.0f));

        should(equalValue(bv1, 1));
        should(equalValue(iv1, 1));
        should(equalValue(fv1, 1.0f));

        should(equalIter(bv3.begin(), bv3.end(), di));
        should(equalIter(iv3.begin(), iv3.end(), di));
        should(equalIter(fv3.begin(), fv3.end(), df));

        should(!equalVector(bv3, fv3));
        should(!equalVector(iv3, fv3));

        BV bv(fv3);
        should(equalIter(bv3.begin(), bv3.end(), bv.begin()));
        should(equalVector(bv3, bv));

        FV fv(iv3);
        should(equalIter(iv3.begin(), iv3.end(), fv.begin()));
        should(equalVector(iv3, fv));

        fv = fv3;
        should(equalIter(fv3.begin(), fv3.end(), fv.begin()));
        should(equalVector(fv3, fv));

        fv = bv3;
        should(equalIter(bv3.begin(), bv3.end(), fv.begin()));
        should(equalVector(bv3, fv));
    }

    void testComparison()
    {
        should(bv0 == bv0);
        should(iv0 == iv0);
        should(fv0 == fv0);
        should(iv0 == bv0);
        should(iv0 == fv0);
        should(fv0 == bv0);

        should(bv3 == bv3);
        should(iv3 == iv3);
        should(fv3 == fv3);
        should(iv3 == bv3);
        should(iv3 != fv3);
        should(fv3 != bv3);
    }

    void testArithmetic()
    {
        IV ivm3 = -iv3;
        FV fvm3 = -fv3;

        int mi[] = { -1, -2, -4};
        float mf[] = { -1.2, -2.4, -3.6 };

        should(equalIter(ivm3.begin(), ivm3.end(), mi));
        should(equalIter(fvm3.begin(), fvm3.end(), mf));

        IV iva3 = abs(ivm3);
        FV fva3 = abs(fvm3);
        should(equalVector(iv3, iva3));
        should(equalVector(fv3, fva3));

        int fmi[] = { -2, -3, -4};
        int fpi[] = { 1, 2, 3};
        IV ivi3 = floor(fvm3);
        should(equalIter(ivi3.begin(), ivi3.end(), fmi));
        ivi3 = -ceil(fv3);
        should(equalIter(ivi3.begin(), ivi3.end(), fmi));
        ivi3 = floor(fv3);
        should(equalIter(ivi3.begin(), ivi3.end(), fpi));
        ivi3 = -ceil(fvm3);
        should(equalIter(ivi3.begin(), ivi3.end(), fpi));

        should(bv1.squaredMagnitude() == 3);
        should(iv1.squaredMagnitude() == 3);
        should(fv1.squaredMagnitude() == 3.0);

        shouldEqualTolerance(fv3.squaredMagnitude(), (1.2f*1.2f + 2.4f*2.4f + 3.6f*3.6f), 1e-7f);

        should(dot(bv3, bv3) == bv3.squaredMagnitude());
        should(dot(iv3, bv3) == iv3.squaredMagnitude());
        should(dot(fv3, fv3) == fv3.squaredMagnitude());

        shouldEqualTolerance(VIGRA_CSTD::sqrt(
                (typename NumericTraits<typename BV::value_type>::RealPromote)
                dot(bv3, bv3)), bv3.magnitude(), 0.0);
        shouldEqualTolerance(VIGRA_CSTD::sqrt(
                (typename NumericTraits<typename IV::value_type>::RealPromote)
                dot(iv3, bv3)), iv3.magnitude(), 0.0);
        shouldEqualTolerance(VIGRA_CSTD::sqrt(
                (typename NumericTraits<typename FV::value_type>::RealPromote)
                dot(fv3, fv3)), fv3.magnitude(), 0.0f);

        BV bv = bv3;
        bv[2] = 200;
        should(dot(bv, bv) == 40005);
        should(bv.squaredMagnitude() == 40005);

        should(equalVector(bv3 - iv3, bv0));
        should(equalVector(fv3 - fv3, fv0));
        BV bvp = (bv3 + bv3)*0.5;
        FV fvp = (fv3 + fv3)*0.5;
        should(equalVector(bvp, bv3));
        should(equalVector(fvp, fv3));
        bvp = 2.0*bv3 - bv3;
        fvp = 2.0*fv3 - fv3;
        should(equalVector(bvp, bv3));
        should(equalVector(fvp, fv3));

        IV ivp = bv + bv;
        int ip1[] = {2, 4, 400};
        should(equalIter(ivp.begin(), ivp.end(), ip1));
        should(equalVector(bv0 - iv1, -iv1));

        bvp = bv3 / 2.0;
        fvp = bv3 / 2.0;
        int ip[] = {1, 1, 2};
        float fp[] = {0.5, 1.0, 2.0};
        should(equalIter(bvp.begin(), bvp.end(), ip));
        should(equalIter(fvp.begin(), fvp.end(), fp));
        fvp = fv3 / 2.0;
        float fp1[] = {0.6, 1.2, 1.8};
        should(equalIter(fvp.begin(), fvp.end(), fp1));
    }

    void testAccessor()
    {
        vigra::VectorAccessor<FV> v;
        FV  pfa[] = {fv3, fv1};
        FV * pf = pfa;

        should(v.size(pf) == 3);
        should(v.size(pf, 1) == 3);
        should(equalVector(v(pf), fv3));
        should(equalIter(v.begin(pf), v.end(pf), fv3.begin()));
        should(v.getComponent(pf, 2) == fv3[2]);
        v.setComponent(5.5, pf, 1);
        should(pf[0][1] == 5.5);
        should(equalVector(v(pf, 1), fv1));
        should(equalIter(v.begin(pf, 1), v.end(pf, 1), fv1.begin()));
        should(v.getComponent(pf, 1, 2) == fv1[2]);
        v.setComponent(5.5, pf, 1, 1);
        should(pf[1][1] == 5.5);
    }
};

struct RGBValueTest
: public TinyVectorTest<vigra::RGBValue<unsigned char>,
                        vigra::RGBValue<int>,
                        vigra::RGBValue<float> >
{
    typedef TinyVectorTest<vigra::RGBValue<unsigned char>,
                           vigra::RGBValue<int>,
                           vigra::RGBValue<float> > Base;

    RGBValueTest()
    : Base()
    {}


    void testRGBAccessors()
    {
        vigra::RGBAccessor<FV> rgb;
        vigra::RedAccessor<FV> red;
        vigra::GreenAccessor<FV> green;
        vigra::BlueAccessor<FV> blue;
        vigra::RGBToGrayAccessor<FV> gray;

        FV pfa[] = { FV(0.0), FV(1.0)};
        FV * pf = pfa;

        should(rgb(pf) == vigra::NumericTraits<FV>::zero());
        should(red(pf) == 0.0);
        should(green(pf) == 0.0);
        should(blue(pf) == 0.0);
        should(gray(pf) == 0.0);

        should(rgb(pf, 1) == vigra::NumericTraits<FV>::one());
        should(red(pf, 1) == 1.0);
        should(green(pf, 1) == 1.0);
        should(blue(pf, 1) == 1.0);
        should(gray(pf, 1) == 1.0);

        rgb.setRed(1.0, pf);
        rgb.setGreen(2.0, pf);
        rgb.setBlue(3.0, pf);
        should(red(pf) == 1.0);
        should(green(pf) == 2.0);
        should(blue(pf) == 3.0);

        red.set(4.0, pf);
        green.set(5.0, pf);
        blue.set(6.0, pf);
        should(rgb.red(pf) == 4.0);
        should(rgb.green(pf) == 5.0);
        should(rgb.blue(pf) == 6.0);
        should(vigra::abs(gray(pf) - 4.81) < 0.00001);

        rgb.setRed(7.0, pf, 1);
        rgb.setGreen(8.0, pf, 1);
        rgb.setBlue(9.0, pf, 1);
        should(red(pf, 1) == 7.0);
        should(green(pf, 1) == 8.0);
        should(blue(pf, 1) == 9.0);

        red.set(10.0, pf, 1);
        green.set(11.0, pf, 1);
        blue.set(12.0, pf, 1);
        should(rgb.red(pf, 1) == 10.0);
        should(rgb.green(pf, 1) == 11.0);
        should(rgb.blue(pf, 1) == 12.0);
        should(vigra::abs(gray(pf, 1) - 10.81) <0.00001);
    }

};

struct PixelTypesTestSuite
: public vigra::test_suite
{
    typedef TinyVectorTest<vigra::TinyVector<unsigned char, 3>,
                           vigra::TinyVector<int, 3>,
                           vigra::TinyVector<float, 3> > TinyVectorTest;

    PixelTypesTestSuite()
    : vigra::test_suite("PixelTypesTest")
    {
        add( testCase(&TinyVectorTest::testConstruction));
        add( testCase(&TinyVectorTest::testComparison));
        add( testCase(&TinyVectorTest::testArithmetic));
        add( testCase(&TinyVectorTest::testAccessor));

        add( testCase(&RGBValueTest::testConstruction));
        add( testCase(&RGBValueTest::testComparison));
        add( testCase(&RGBValueTest::testArithmetic));
        add( testCase(&RGBValueTest::testAccessor));
        add( testCase(&RGBValueTest::testRGBAccessors));
    }
};

int main()
{
    PixelTypesTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
