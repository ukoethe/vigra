#include "g++_relops_workaround.hxx"

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
