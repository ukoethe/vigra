#include "g++_relops_workaround.hxx"

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
