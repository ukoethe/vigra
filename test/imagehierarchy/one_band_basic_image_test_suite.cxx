#include "g++_relops_workaround.hxx"

#include <unittest.h> 
#include "one_band_image_policy.hxx" 
#include "basic_image_test.hxx" 
// 
// struct OneBandBasicImageTestSuite
// : public vigra::test_suite
// {
OneBandBasicImageTestSuite::OneBandBasicImageTestSuite()
    : vigra::test_suite("OneBandBasicImageTestSuite")
    {
        add( new BasicImageTestSuite<OneBandImagePolicy<vigra::IImage> > ("vigra::IImage"));
        add( new BasicImageTestSuite<OneBandImagePolicy<vigra::FImage> > ("vigra::FImage")); 
        add( new BasicImageTestSuite<OneBandImagePolicy<vigra::BImage> > ("vigra::BImage"));
        add( new BasicImageTestSuite<OneBandImagePolicy<vigra::SImage> > ("vigra::SImage"));
        add( new BasicImageTestSuite<OneBandImagePolicy<vigra::DImage> > ("vigra::DImage"));
    }
// };
// 
// int main()
// {
//  
//     OneBandBasicImageTestSuite itest;
// 
//     int failed = itest.run();
// 
//     std::cout << itest.report() << std::endl;
//     
//     return (failed != 0);
// }  
