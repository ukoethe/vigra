#include "vector_3_image_policy.hxx"
#include "basic_image_test.hxx"

Vector3BasicImageTestSuite::Vector3BasicImageTestSuite()
    : vigra::test_suite(" Vector3BasicImageTestSuite")
    {
        add ( new BasicImageTestSuite<Vector3ImagePolicy<vigra::FVector3Image> > ("vigra::FVector3Image"));
        add ( new BasicImageTestSuite<Vector3ImagePolicy<vigra::DVector3Image> > ("vigra::DVector3Image"));
    }

// int main()
// {
//     Vector3BasicImageTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
