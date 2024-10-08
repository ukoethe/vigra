#ifndef RGB_BASIC_IMAGE_TEST_HXX
#define RGB_BASIC_IMAGE_TEST_HXX

#include "basic_image_test.hxx"

template<class Policy>
class RGBBasicImageTest
: public BasicImageTest<Policy>
{
public:
    typedef typename Policy::Image Image;
    typedef typename BasicImageTest<Policy>::ChildImage ChildImage;
    typedef typename Policy::value_type value_type;
    typename Policy::data_array_type data;

    RGBBasicImageTest()
    : data(Policy::getData())
    {}

    void testRed()
    {
        std::auto_ptr<Image> image(Policy::factory(2,3, data[0]));
        should((*image)(1,1)[0] == data[0].red());
        should((*image)[Diff2D(1,2)][0] == data[0].red());
    }
    void testGreen()
    {
        std::auto_ptr<Image> image(Policy::factory(2,3, data[3]));
        should((*image)(1,1)[1] == data[3].green());
        should((*image)[Diff2D(1,2)][1] == data[3].green());
    }
    void testBlue()
    {
        std::auto_ptr<Image> image(Policy::factory(2,3, data[5]));
        should((*image)(1,1)[2] == data[5].blue());
        should((*image)[Diff2D(1,2)][2] == data[5].blue());
    }
};


template<class Policy>
class RGBBasicImgTestSuite
: public BasicImageTestSuite<Policy>
{
public:
    // den Unterschied von RGBBasicImgTestSuite und RGBBasicImageTestSuite beachten
    RGBBasicImgTestSuite(const char * to_test_Image_name)
    : BasicImageTestSuite<Policy> (to_test_Image_name)
    {
        add ( testCase( &RGBBasicImageTest<Policy>::testRed));
        add ( testCase( &RGBBasicImageTest<Policy>::testGreen));
        add ( testCase( &RGBBasicImageTest<Policy>::testBlue));
    }
};

#endif // RGB_BASIC_IMAGE_TEST_HXX
