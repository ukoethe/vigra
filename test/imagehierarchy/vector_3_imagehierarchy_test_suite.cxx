#ifdef __GNUC__

#   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
#   define __STL_END_RELOPS_NAMESPACE }}
#   define __STD_RELOPS std::rel_ops
 
#include <stl_relops.h>

#   undef __STL_BEGIN_RELOPS_NAMESPACE 
#   undef __STL_END_RELOPS_NAMESPACE 
#   undef __STD_RELOPS

#endif

#include "vector_3_image_policy.hxx"
#include "imagehierarchy_test.hxx"

Vector3ImageHierarchyTestSuite::Vector3ImageHierarchyTestSuite()
    : vigra::test_suite(" Vector3ImageHierarchyTestSuite")
    {
        add ( new ImageHierarchyTestSuite<Vector3ImagePolicy<vigra::Vector3Image> > ("vigra::Vector3Image"));
    }

// int main()
// {
//     Vector3ImageHierarchyTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
