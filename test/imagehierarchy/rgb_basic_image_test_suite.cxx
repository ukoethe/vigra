#include "rgb_images_policy.hxx"
#include "rgb_basic_image_test.hxx"

RGBBasicImageTestSuite::RGBBasicImageTestSuite()
    : test_suite("RGBBasicImageTestSuite")
    {
        // den Unterschied zwischen RGBBasicImageTestSuite und RGBBasicImgTestSuite beachten !!!!
        add ( new RGBBasicImgTestSuite<RGBImagePolicy<vigra::IRGBImage> > ("vigra::IRGBImage"));
        add ( new RGBBasicImgTestSuite<RGBImagePolicy<vigra::FRGBImage> > ("vigra::FRGBImage"));
        add ( new RGBBasicImgTestSuite<RGBImagePolicy<vigra::DRGBImage> > ("vigra::DRGBImage"));
        add ( new RGBBasicImgTestSuite<RGBImagePolicy<vigra::BRGBImage> > ("vigra::BRGBImage"));
    }

// int main()
// {
//     RGBBasicImageTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
