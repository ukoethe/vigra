#include "vector_2_image_policy.hxx"
#include "basic_image_test.hxx"

Vector2BasicImageTestSuite ::Vector2BasicImageTestSuite()
    : vigra::test_suite(" Vector2BasicImageTestSuite")
    {
        add ( new BasicImageTestSuite<Vector2ImagePolicy<vigra::FVector2Image> > ("vigra::FVector2Image"));
        add ( new BasicImageTestSuite<Vector2ImagePolicy<vigra::DVector2Image> > ("vigra::DVector2Image"));
    }

// int main()
// {
//     Vector2BasicImageTestSuite suite;
//     int failed = suite.run();
//     std::cout << suite.report() << std::endl;
//     return (failed != 0);
// }
