#include "g++_relops_workaround.hxx"

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
