#ifndef RGB_IMAGEHIERARCHY_TEST_HXX
#define RGB_IMAGEHIERARCHY_TEST_HXX

#include "imagehierarchy_test.hxx"

template<class Policy>
class RGBImageHierarchyTest
: public ImageHierarchyTest<Policy>
{
public:
    typedef typename Policy::Image Image;
    typedef typename ImageHierarchyTest<Policy>::ChildImage ChildImage;
    typedef typename Policy::value_type value_type;
    typename Policy::data_array_type data;

    RGBImageHierarchyTest()
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
class RGBImgHierarchyTestSuite
: public ImageHierarchyTestSuite<Policy>
{
public:
    // den Unterschied von RGBImgHierarchyTestSuite und RGBImageHierarchyTestSuite beachten
    RGBImgHierarchyTestSuite(const char * to_test_Image_name)
    : ImageHierarchyTestSuite<Policy> (to_test_Image_name)
    {
        add ( testCase( &RGBImageHierarchyTest<Policy>::testRed));
        add ( testCase( &RGBImageHierarchyTest<Policy>::testGreen));
        add ( testCase( &RGBImageHierarchyTest<Policy>::testBlue));
    }
};
#endif // RGB_IMAGEHIERARCHY_TEST_HXX
