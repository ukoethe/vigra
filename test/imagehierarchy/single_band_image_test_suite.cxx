#include "g++_relops_workaround.hxx"

#include "single_band_image_policy.hxx" 
#include "one_band_image_policy.hxx" 
#include "single_band_image_test.hxx"
#include <unittest.hxx> 

/** Vorsicht !!! in select_image_test.hxx gibt's SingleBandImageTestSuite. Den Unterschied beachten!!!
*/
SingleImageTestSuite::SingleImageTestSuite()
    : vigra::test_suite("SingleBandImageTestSuite")
    {
        add( new SingleBandImageTestSuite<SingleBandImagePolicy<OneBandImagePolicy<vigra::GrayImage> > > ("SingleBandImage an vigra::GrayImage"));
    }

// int main()
// {
//  
//     SingleImageTestSuite itest;
// 
//     int failed = itest.run();
// 
//     std::cout << itest.report() << std::endl;
//     
//     return (failed != 0);
// }  
