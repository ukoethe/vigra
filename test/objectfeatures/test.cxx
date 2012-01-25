#include <iostream>
#include <sstream>

#include "unittest.hxx"

#include "vigra/accessor.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/rgbvalue.hxx"

#include "vigra/coordinate_iterator.hxx"
#include "vigra/object_features.hxx"

#include "vigra/multi_pointoperators.hxx"

using namespace vigra;

using vigra::type_lists::use_template_list;
using vigra::Accumulators;

typedef float pixel_type;

typedef vigra::MultiArrayShape<2>::type  shape_2d;
typedef vigra::MultiArray<2, pixel_type> array_2d;

typedef vigra::RGBValue<double> rgb_type;

template <class T>
struct gen_accumulators
    : public use_template_list<Accumulators, T,
                               acc::Variance,
                               acc::Skewness,
                               acc::Kurtosis
                              >
{};

typedef gen_accumulators<pixel_type> pixel_accumulators;

typedef use_template_list<Accumulators, StridePairPointer<2, pixel_type>,
                          acc::Count,
                          acc::Sum,
                          acc::CoordSum
                         >
cv_test_accumulators;


typedef use_template_list<Accumulators, pixel_type,
                          acc::Count,
                          acc::Sum,
                          acc::Mean,
                          acc::Moment2
                         >
yy_test_accumulators;

typedef use_template_list<Accumulators, pixel_type,
                          acc::Mean,
                          acc::Moment2
                         >
yyy_test_accumulators;

struct count : public InitializerTag
{
    mutable unsigned long t;
    count() : t(0) {}
    unsigned long operator()() const
    {
        return t++;
    }
};

void test1()
{
    pixel_accumulators pix_acc;
    yyy_test_accumulators yy_acc;
    
	array_2d image(array_2d::size_type(10, 10));

    for (unsigned i = 0; i != 100; ++i)
        image(i / 10, i % 10) = i;

	inspectMultiArray(srcMultiArrayRange(image), pix_acc);

    shouldEqualTolerance(get<acc::Count>(pix_acc), 100, 1e-5);
    shouldEqualTolerance(get<acc::Mean>(pix_acc), 49.5, 1e-5);
    shouldEqualTolerance(get<acc::Moment2>(pix_acc), 83325, 1e-5);
    shouldEqualTolerance(get<acc::Moment2_2>(pix_acc), 83325, 1e-5);
    shouldEqualTolerance(get<acc::Variance>(pix_acc), 833.25, 1e-5);
    shouldEqualTolerance(get<acc::Skewness>(pix_acc), 0, 1e-5);
    shouldEqualTolerance(get<acc::Kurtosis>(pix_acc), 0.0179976, 1e-5);
}

void test2()
{
	array_2d stride_image2(array_2d::size_type(2, 3));
    initMultiArray(destMultiArrayRange(stride_image2), count());

    cv_test_accumulators cv_acc;

    inspectMultiArray(srcCoordinateMultiArrayRange(stride_image2), cv_acc);

    shouldEqualTolerance(get<acc::Count>(cv_acc), 6, 1e-5);
    shouldEqualTolerance(get<acc::Sum>(cv_acc), 15, 1e-5);
    { TinyVector<double, 2> x(3, 6); shouldEqualSequenceTolerance(x.begin(), x.end(), get<acc::CoordSum>(cv_acc).begin(), 1e-5); };

    rgb_type epstest = 1e-5;
    { TinyVector<double, 3> x(1e-05, 1e-05, 1e-05); shouldEqualSequenceTolerance(x.begin(), x.end(), epstest.begin(), 1e-5); };
}

struct ObjectFeaturesTest1
{
    void run1() { test1(); }
    void run2() { test2(); }
    
};

struct ObjectFeaturesTestSuite : public vigra::test_suite
{
    ObjectFeaturesTestSuite()
        : vigra::test_suite("ObjectFeaturesTestSuite")
    {
        add(testCase(&ObjectFeaturesTest1::run1));
        add(testCase(&ObjectFeaturesTest1::run2));
    }
};


int main(int argc, char** argv)
{
    ObjectFeaturesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
