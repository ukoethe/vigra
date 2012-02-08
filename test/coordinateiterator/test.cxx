#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <sstream>
#include "unittest.hxx"

#include "vigra/multi_pointoperators.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/coordinate_iterator.hxx"

#include "vigra/object_features.hxx"


using namespace vigra;

typedef float pixel_type;

typedef vigra::MultiArrayShape<2>::type  shape_2d;
typedef vigra::MultiArray<2, pixel_type> array_2d;
typedef vigra::MultiArrayShape<3>::type  shape_3d;
typedef vigra::MultiArray<3, pixel_type> array_3d;

typedef vigra::RGBValue<double> rgb_type;

typedef vigra::MultiArrayShape<2>::type  shape_rgb;
typedef vigra::MultiArray<2, rgb_type>   array_rgb;

using namespace vigra;

struct count : public InitializerTag
{
    mutable unsigned long t;
    count() : t(0) {}
    unsigned long operator()() const
    {
        return t++;
    }
};

void set_pr(std::ostream & os)
{
    os << std::fixed << std::setprecision(3);
}

struct echo
{
    std::ostringstream os;
    std::string empty;
    echo() { set_pr(os); }
    template<class X>
    void operator()(const X & x)
    {
        os.str(empty);
        os << '"' << x << '"' << ',';
        std::cout << os.str() << "\n";
    }
};

typedef const char* char_s;

struct echo_cmp : public echo
{
    typedef char_s* char_p;
    char_p p;
    echo_cmp() : p(0) {}
    echo_cmp & operator()(char_p r)
    {
        p = r;
        return *this;
    }
    template<class X>
    void operator()(const X & x)
    {
        os.str(empty);
        os << x;
        should (p && *p);
        if (p && *p)
        {
            shouldEqual(os.str(), *p);
            ++p;
        }
    }
};

void test1()
{
    StridePair<2>::coord_type c2_origin(10,  0);
    StridePair<2>::coord_type c2_stride( 1,  2);
    StridePair<2>::index_type c2_shape ( 2,  5);
    
    CoordinateMultiIterator<2>(c2_origin, c2_stride, c2_shape);

    coordinateMultiRange<2>(c2_shape, c2_stride, c2_origin);
    echo_cmp echo_f;
    char_s out_1[] =
    {
        "(10.000, 0.000)",
        "(11.000, 0.000)",
        "(10.000, 2.000)",
        "(11.000, 2.000)",
        "(10.000, 4.000)",
        "(11.000, 4.000)",
        "(10.000, 6.000)",
        "(11.000, 6.000)",
        "(10.000, 8.000)",
        "(11.000, 8.000)",
        0
    };
    (inspectMultiArray(coordinateMultiRange<2>(
                               c2_shape, c2_stride, c2_origin), echo_f(out_1)));

    StridePair<2>::coord_type c7_origin(10,  5.5);
    StridePair<2>::coord_type c7_stride( 1 / 7.0,  1 / 3.0);
    StridePair<2>::index_type c7_shape ( 7,  3);
    char_s out_2[] =
    {
        "(10.000, 5.500)",
        "(10.143, 5.500)",
        "(10.286, 5.500)",
        "(10.429, 5.500)",
        "(10.571, 5.500)",
        "(10.714, 5.500)",
        "(10.857, 5.500)",
        "(10.000, 5.833)",
        "(10.143, 5.833)",
        "(10.286, 5.833)",
        "(10.429, 5.833)",
        "(10.571, 5.833)",
        "(10.714, 5.833)",
        "(10.857, 5.833)",
        "(10.000, 6.167)",
        "(10.143, 6.167)",
        "(10.286, 6.167)",
        "(10.429, 6.167)",
        "(10.571, 6.167)",
        "(10.714, 6.167)",
        "(10.857, 6.167)",
        0
    };
    (inspectMultiArray(coordinateMultiRange<2>(
                               c7_shape, c7_stride, c7_origin), echo_f(out_2)));

    StridePair<3>::coord_type c3_origin(20,  -3,  1 / 3.0);
    StridePair<3>::coord_type c3_stride(0.5, 1/17.0,  1 / 6.0);
    StridePair<3>::index_type c3_shape (  2,  3,    3);
    char_s out_3[] =
    {
        "(20.000, -3.000, 0.333)",
        "(20.500, -3.000, 0.333)",
        "(20.000, -2.941, 0.333)",
        "(20.500, -2.941, 0.333)",
        "(20.000, -2.882, 0.333)",
        "(20.500, -2.882, 0.333)",
        "(20.000, -3.000, 0.500)",
        "(20.500, -3.000, 0.500)",
        "(20.000, -2.941, 0.500)",
        "(20.500, -2.941, 0.500)",
        "(20.000, -2.882, 0.500)",
        "(20.500, -2.882, 0.500)",
        "(20.000, -3.000, 0.667)",
        "(20.500, -3.000, 0.667)",
        "(20.000, -2.941, 0.667)",
        "(20.500, -2.941, 0.667)",
        "(20.000, -2.882, 0.667)",
        "(20.500, -2.882, 0.667)",
        0
    };
    (inspectMultiArray(coordinateMultiRange<3>(
                               c3_shape, c3_stride, c3_origin), echo_f(out_3)));

    StridePair<2>::coord_type c0_stride(0.5, 0.25);
    StridePair<2>::index_type c0_shape (  2,  3);
    char_s out_4[] =
    {
        "(0.000, 0.000)",
        "(0.500, 0.000)",
        "(0.000, 0.250)",
        "(0.500, 0.250)",
        "(0.000, 0.500)",
        "(0.500, 0.500)",
        0
    };
    (inspectMultiArray(coordinateMultiRange<2>(c0_shape, c0_stride),
                                echo_f(out_4)));

    StridePair<2>::index_type c1_shape (  2,  3);
    char_s out_5[] =
    {
        "(0.000, 0.000)",
        "(1.000, 0.000)",
        "(0.000, 1.000)",
        "(1.000, 1.000)",
        "(0.000, 2.000)",
        "(1.000, 2.000)",
        0
    };
    (inspectMultiArray(coordinateMultiRange<2>(c1_shape),
                                echo_f(out_5)));

    array_2d stride_image2(array_2d::size_type(2, 3));
    initMultiArray(destMultiArrayRange(stride_image2), count());
    char_s out_6[] =
    {
        "[0.000, (0.000, 0.000)]",
        "[1.000, (1.000, 0.000)]",
        "[2.000, (0.000, 1.000)]",
        "[3.000, (1.000, 1.000)]",
        "[4.000, (0.000, 2.000)]",
        "[5.000, (1.000, 2.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRange(stride_image2),
                                echo_f(out_6)));
    
    array_3d stride_image3(array_3d::size_type(2, 3, 4));
    initMultiArray(destMultiArrayRange(stride_image3), count());
    char_s out_7[] =
    {
        "[0.000, (0.000, 0.000, 0.000)]",
        "[1.000, (1.000, 0.000, 0.000)]",
        "[2.000, (0.000, 1.000, 0.000)]",
        "[3.000, (1.000, 1.000, 0.000)]",
        "[4.000, (0.000, 2.000, 0.000)]",
        "[5.000, (1.000, 2.000, 0.000)]",
        "[6.000, (0.000, 0.000, 1.000)]",
        "[7.000, (1.000, 0.000, 1.000)]",
        "[8.000, (0.000, 1.000, 1.000)]",
        "[9.000, (1.000, 1.000, 1.000)]",
        "[10.000, (0.000, 2.000, 1.000)]",
        "[11.000, (1.000, 2.000, 1.000)]",
        "[12.000, (0.000, 0.000, 2.000)]",
        "[13.000, (1.000, 0.000, 2.000)]",
        "[14.000, (0.000, 1.000, 2.000)]",
        "[15.000, (1.000, 1.000, 2.000)]",
        "[16.000, (0.000, 2.000, 2.000)]",
        "[17.000, (1.000, 2.000, 2.000)]",
        "[18.000, (0.000, 0.000, 3.000)]",
        "[19.000, (1.000, 0.000, 3.000)]",
        "[20.000, (0.000, 1.000, 3.000)]",
        "[21.000, (1.000, 1.000, 3.000)]",
        "[22.000, (0.000, 2.000, 3.000)]",
        "[23.000, (1.000, 2.000, 3.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRange(stride_image3),
                                echo_f(out_7)));

    array_rgb image_rgb(array_2d::size_type(2, 3));
    initMultiArray(destMultiArrayRange(image_rgb), count());
    char_s out_8[] =
    {
        "[(0.000, 0.000, 0.000), (0.000, 0.000)]",
        "[(1.000, 1.000, 1.000), (1.000, 0.000)]",
        "[(2.000, 2.000, 2.000), (0.000, 1.000)]",
        "[(3.000, 3.000, 3.000), (1.000, 1.000)]",
        "[(4.000, 4.000, 4.000), (0.000, 2.000)]",
        "[(5.000, 5.000, 5.000), (1.000, 2.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRange(image_rgb),
                                echo_f(out_8)));
    typedef AccessorTraits<rgb_type>::default_const_accessor
        default_rgb_access;
    char_s out_9[] =
    {
        "[(0.000, 0.000, 0.000), (0.000, 0.000)]",
        "[(1.000, 1.000, 1.000), (1.000, 0.000)]",
        "[(2.000, 2.000, 2.000), (0.000, 1.000)]",
        "[(3.000, 3.000, 3.000), (1.000, 1.000)]",
        "[(4.000, 4.000, 4.000), (0.000, 2.000)]",
        "[(5.000, 5.000, 5.000), (1.000, 2.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRangeAccessor(image_rgb,
                                                          default_rgb_access()),
                      echo_f(out_9)));
    char_s out_10[] =
    {
        "[0.000, (0.000, 0.000)]",
        "[1.000, (1.000, 0.000)]",
        "[2.000, (0.000, 1.000)]",
        "[3.000, (1.000, 1.000)]",
        "[4.000, (0.000, 2.000)]",
        "[5.000, (1.000, 2.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRangeAccessor(image_rgb,
                                                     GreenAccessor<rgb_type>()),
                      echo_f(out_10)));

    StridePair<2>::coord_type x2_origin(0.5,  2.2);
    const StridePair<2>::coord_type x2_stride(1.5, 0.4);
    char_s out_11[] =
    {
        "[0.000, (0.000, 0.000)]",
        "[1.000, (1.500, 0.000)]",
        "[2.000, (0.000, 0.400)]",
        "[3.000, (1.500, 0.400)]",
        "[4.000, (0.000, 0.800)]",
        "[5.000, (1.500, 0.800)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRange(stride_image2,
                                                             x2_stride),
                      echo_f(out_11)));
    char_s out_12[] =
    {
        "[0.000, (0.500, 2.200)]",
        "[1.000, (2.000, 2.200)]",
        "[2.000, (0.500, 2.600)]",
        "[3.000, (2.000, 2.600)]",
        "[4.000, (0.500, 3.000)]",
        "[5.000, (2.000, 3.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRange(stride_image2,
                                                   x2_stride,
                                                   x2_origin),
                      echo_f(out_12)));
    char_s out_13[] =
    {
        "[0.000, (0.000, 0.000)]",
        "[1.000, (1.500, 0.000)]",
        "[2.000, (0.000, 0.400)]",
        "[3.000, (1.500, 0.400)]",
        "[4.000, (0.000, 0.800)]",
        "[5.000, (1.500, 0.800)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRangeAccessor(image_rgb,
                                                    GreenAccessor<rgb_type>(),
                                                    x2_stride),
                      echo_f(out_13)));
    char_s out_14[] =
    {
        "[0.000, (0.500, 2.200)]",
        "[1.000, (2.000, 2.200)]",
        "[2.000, (0.500, 2.600)]",
        "[3.000, (2.000, 2.600)]",
        "[4.000, (0.500, 3.000)]",
        "[5.000, (2.000, 3.000)]",
        0
    };
    (inspectMultiArray(srcCoordinateMultiArrayRangeAccessor(image_rgb,
                                                    GreenAccessor<rgb_type>(),
                                                    x2_stride,
                                                    x2_origin),
                      echo_f(out_14)));
}

struct CoordinateIteratorTest1
{
    void run() { test1(); }
    
};

struct CoordinateIteratorTestSuite : public vigra::test_suite
{
    CoordinateIteratorTestSuite()
        : vigra::test_suite("CoordinateIteratorTestSuite")
    {
        add(testCase(&CoordinateIteratorTest1::run));
    }
};


int main(int argc, char** argv)
{
    CoordinateIteratorTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
