#   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
#   define __STL_END_RELOPS_NAMESPACE }}
#   define __STD_RELOPS std::rel_ops
 
#include <stl_relops.h>

#   undef __STL_BEGIN_RELOPS_NAMESPACE 
#   undef __STL_END_RELOPS_NAMESPACE 
#   undef __STD_RELOPS

#include "vector_4_image_policy.hxx"
#include "imagehierarchy_test.hxx"

Vector4ImageHierarchyTestSuite::Vector4ImageHierarchyTestSuite()
    : vigra::test_suite(" Vector4ImageHierarchyTestSuite")
    {
        add ( new ImageHierarchyTestSuite<Vector4ImagePolicy<vigra::Vector4Image> > ("vigra::Vector4Image"));
    }

// int main()
// {
//     Vector4ImageHierarchyTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
