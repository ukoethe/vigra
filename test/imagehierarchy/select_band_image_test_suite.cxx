#include "g++_relops_workaround.hxx"

#include "select_band_image_policy.hxx" 
#include "select_image_test.hxx"
#include "vector_2_image_policy.hxx" 
// #include "vector_3_image_policy.hxx" 
// #include "vector_4_image_policy.hxx" 
#include "rgb_images_policy.hxx" 
// #include "select_image_test.hxx" 
#include <unittest.hxx> 

/* Vorsicht !!! in select_image_test.hxx gibt's SelectBandImageTestSuite. Den Unterschied beachten!!!
*/
SelectImagesTestSuite::SelectImagesTestSuite()
    : vigra::test_suite("SelectBandImageTestSuite")
    {
        add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector2ImagePolicy<vigra::Vector2Image>, 0> > ("SelectBandImage an vigra::Vector2Image, selected band is 0"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector2ImagePolicy<vigra::Vector2Image>, 1> >("SelectBandImage an vigra::Vector2Image, selected band is 1"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 0> >("SelectBandImage an vigra::Vector3Image, selected band is 0"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 1> >("SelectBandImage an vigra::Vector3Image, selected band is 1"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 2> >("SelectBandImage an vigra::Vector3Image, selected band is 2"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 0> >("SelectBandImage an vigra::Vector4Image, selected band is 0"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 1> >("SelectBandImage an vigra::Vector4Image, selected band is 1"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 2> >("SelectBandImage an vigra::Vector4Image, selected band is 2"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 3> >("SelectBandImage an vigra::Vector4Image, selected band is 3"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 0> >("SelectBandImage an vigra::RGBImage, selected band is 0"));
//         add(new SelectBandImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 1> >("SelectBandImage an vigra::RGBImage, selected band is 1"));
        add(new SelectBandImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 2> >("SelectBandImage an vigra::RGBImage, selected band is 2"));
    }

// int main()
// {
//  
//     SelectImagesTestSuite itest;
// 
//     int failed = itest.run();
// 
//     std::cout << itest.report() << std::endl;
//     
//     return (failed != 0);
// }  
