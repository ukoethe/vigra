#include "g++_relops_workaround.hxx"

#include "vector_4_image_policy.hxx"
#include "basic_image_test.hxx"

Vector4BasicImageTestSuite::Vector4BasicImageTestSuite()
    : vigra::test_suite(" Vector4BasicImageTestSuite")
    {
        add ( new BasicImageTestSuite<Vector4ImagePolicy<vigra::FVector4Image> > ("vigra::FVector4Image"));
        add ( new BasicImageTestSuite<Vector4ImagePolicy<vigra::DVector4Image> > ("vigra::DVector4Image"));
    }

// int main()
// {
//     Vector4BasicImageTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
