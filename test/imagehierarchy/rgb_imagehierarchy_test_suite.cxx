#   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
#   define __STL_END_RELOPS_NAMESPACE }}
#   define __STD_RELOPS std::rel_ops
 
#include <stl_relops.h>

#   undef __STL_BEGIN_RELOPS_NAMESPACE 
#   undef __STL_END_RELOPS_NAMESPACE 
#   undef __STD_RELOPS

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
