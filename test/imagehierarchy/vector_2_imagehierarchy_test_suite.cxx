#ifdef __GNUC__
#   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
#   define __STL_END_RELOPS_NAMESPACE }}
#   define __STD_RELOPS std::rel_ops
 
#include <stl_relops.h>

#   undef __STL_BEGIN_RELOPS_NAMESPACE 
#   undef __STL_END_RELOPS_NAMESPACE 
#   undef __STD_RELOPS

#endif

#include "vector_2_image_policy.hxx"
#include "imagehierarchy_test.hxx"

Vector2ImageHierarchyTestSuite::Vector2ImageHierarchyTestSuite()
    : vigra::test_suite(" Vector2ImageHierarchyTestSuite")
    {
        add ( new ImageHierarchyTestSuite<Vector2ImagePolicy<vigra::Vector2Image> > ("vigra::Vector2Image"));
    }

// int main()
// {
//     Vector2ImageHierarchyTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
