#include "g++_relops_workaround.hxx"


#include "rgb_images_policy.hxx"
#include "rgb_imagehierarchy_test.hxx"

RGBImageHierarchyTestSuite::RGBImageHierarchyTestSuite()
    : test_suite("RGBImageHierarchyTestSuite")
    {
        // den Unterschied zwischen RGBImageHierarchyTestSuite und RGBImgHierarchyTestSuite beachten !!!!
        add ( new RGBImgHierarchyTestSuite<RGBImagePolicy<vigra::RGBImage> > ("vigra::RGBImage"));
    }

// int main()
// {
//     RGBImageHierarchyTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
