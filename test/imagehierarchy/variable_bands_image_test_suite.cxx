#include "g++_relops_workaround.hxx"

#include "variable_bands_image_policy.hxx" 
#include "vector_2_image_policy.hxx" 
#include "vector_3_image_policy.hxx" 
#include "vector_4_image_policy.hxx" 
#include "rgb_images_policy.hxx" 
#include "one_band_image_policy.hxx" 
#include "variable_bands_image_test.hxx" 

/** Vorsicht !!! in select_image_test.hxx gibt's VariableBandImageTestSuite. Den Unterschied beachten!!!
*/
VariableImagesTestSuite::VariableImagesTestSuite()
    : vigra::test_suite("VariableBandImageTestSuite")
    {
        add( new VariableBandsImageTestSuite<VariableBandsImagePolicy<OneBandImagePolicy<vigra::GrayImage> > > ("VariableBandsImage an GrayImage"));
        add( new VariableBandsImageTestSuite<VariableBandsImagePolicy<RGBImagePolicy<vigra::RGBImage> > > ("VariableBandsImage an RGBImage"));
        add( new VariableBandsImageTestSuite<VariableBandsImagePolicy<Vector2ImagePolicy<vigra::Vector2Image> > > ("VariableBandsImage an Vector2Image"));
//        add( new VariableBandsImageTestSuite<VariableBandsImagePolicy<Vector3ImagePolicy<vigra::Vector3Image> > > ("VariableBandsImage an Vector3Image"));
        add( new VariableBandsImageTestSuite<VariableBandsImagePolicy<Vector4ImagePolicy<vigra::Vector4Image> > > ("VariableBandsImage an Vector4Image"));
    }

// int main()
// {
//  
//     VariableImagesTestSuite itest;
// 
//     int failed = itest.run();
// 
//     std::cout << itest.report() << std::endl;
//     
//     return (failed != 0);
// }  
