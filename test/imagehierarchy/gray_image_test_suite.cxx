#include "g++_relops_workaround.hxx"

#include <unittest.hxx>
#include "one_band_image_policy.hxx"
#include "imagehierarchy_test.hxx"

GrayImageTestSuite::GrayImageTestSuite()
    : vigra::test_suite(" GrayImageTestSuite")
    {
        add ( new ImageHierarchyTestSuite<OneBandImagePolicy<vigra::GrayImage> > ("vigra::GrayImage"));
    }

// int main()
// {
//     GrayImageTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
