// #   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
// #   define __STL_END_RELOPS_NAMESPACE }}
// #   define __STD_RELOPS std::rel_ops
//  
// #include <stl_relops.h>
// 
// #   undef __STL_BEGIN_RELOPS_NAMESPACE 
// #   undef __STL_END_RELOPS_NAMESPACE 
// #   undef __STD_RELOPS

#include <iostream>
#include <unittest.h>
#include "policy.hxx"
#include "NewImHier.hxx"


using vigra::Diff2D;

#ifdef RGB_IMAGE
using vigra::FRGBImage;
using vigra::IRGBImage;
using vigra::BRGBImage;
using vigra::DRGBImage;
#endif //RGB_IMAGE

#ifdef nRGB_IMAGE
using vigra::FImage;
using vigra::IImage;
using vigra::BImage;
using vigra::SImage;
using vigra::DImage; 
#endif //nRGB_IMAGE

#ifdef VImage
using vigra::FVector2Image;
using vigra::FVector3Image;
using vigra::FVector4Image;
using vigra::DVector2Image;
using vigra::DVector3Image;
using vigra::DVector4Image;
#endif //VImage

#ifdef GRAY_IMAGE
using vigra::GrayImage;
#endif //GRAY_IMAGE

#ifdef FIXED_RGB_IMAGE
using vigra::RGBImage;
#endif //FIXED_RGB_IMAGE

#ifdef FIXED_V_IMAGE
using vigra::Vector2Image;
using vigra::Vector3Image;
using vigra::Vector4Image;
#endif //FIXED_V_IMAGE

#ifdef VARIABLE_BANDS_IMAGE
using vigra::VariableBandsImage;
#endif //VARIABLE_BANDS_IMAGE

#define createImageTest(Image) new TImageTestSuite<ImagePolicy<Image> >(#Image)
#define createRGBImageTest(Image) new RGBImageTestSuite<RGBImagePolicy<Image> >(#Image)
 
struct ImageTestSuite
: public vigra::test_suite
{
    ImageTestSuite()
    : vigra::test_suite("ImageTestSuite")
    {
#ifdef VImage  
        //add( new TImageTestSuite<Vector2ImagePolicy<FVector2Image> >("FVector2Image"));// fertig
        add( new TImageTestSuite<Vector3ImagePolicy<FVector3Image> >("FVector3Image"));// fertig
        //add( new TImageTestSuite<Vector4ImagePolicy<FVector4Image> >("FVector4Image"));// fertig
        //add( new TImageTestSuite<Vector2ImagePolicy<DVector2Image> >("DVector2Image"));// fertig
        add( new TImageTestSuite<Vector3ImagePolicy<DVector3Image> >("DVector3Image"));// fertig
        //add( new TImageTestSuite<Vector4ImagePolicy<DVector4Image> >("DVector4Image"));// fertig
#endif //VImage
   
#ifdef RGB_IMAGE
        add( createRGBImageTest(IRGBImage));
        add( createRGBImageTest(FRGBImage));
        add( createRGBImageTest(BRGBImage));
        add( createRGBImageTest(DRGBImage));
#endif //RGB_IMAGE
 
#if defined(nRGB_IMAGE)     

//         add( createImageTest(IImage));// fertig
//         add( createImageTest(FImage));// fertig 
//         add( createImageTest(BImage));// fertig
//         add( createImageTest(SImage));// fertig
        add( createImageTest(DImage));// fertig

#endif //nRGB_IMAGE  && ndefined(VARIABLE_BANDS_IMAGE)        

#ifdef GRAY_IMAGE

       add( createImageTest(GrayImage));// fertig       

#endif //GRAY_IMAGE

#ifdef FIXED_RGB_IMAGE

        add( new RGBImageTestSuite<RGBImagePolicy<RGBImage> >("RGBImage"));

#endif //FIXED_RGB_IMAGE 

#ifdef FIXED_V_IMAGE
        add( new TImageTestSuite<Vector2ImagePolicy<Vector2Image> >("Vector2Image"));// fertig
        add( new TImageTestSuite<Vector3ImagePolicy<Vector3Image> >("Vector3Image"));// fertig
        add( new TImageTestSuite<Vector4ImagePolicy<Vector4Image> >("Vector4Image"));// fertig
#endif //FIXED_V_IMAGE

#if defined(VARIABLE_BANDS_IMAGE)
  
//        add( new TImageTestSuite<VariableBandsImagePolicy<ImagePolicy<vigra::GrayImage> > > ("VariableBandsImage an GrayImage"));
//        add( new RGBImageTestSuite<VariableBandsImagePolicy<RGBImagePolicy<vigra::RGBImage> > > ("VariableBandsImage an RGBImage"));
//         add( new TImageTestSuite<VariableBandsImagePolicy<Vector2ImagePolicy<vigra::Vector2Image> > > ("VariableBandsImage an Vector2Image"));
//         add( new TImageTestSuite<VariableBandsImagePolicy<Vector3ImagePolicy<vigra::Vector3Image> > > ("VariableBandsImage an Vector3Image"));
//         add( new TImageTestSuite<VariableBandsImagePolicy<Vector4ImagePolicy<vigra::Vector4Image> > > ("VariableBandsImage an Vector4Image"));
        
#endif //VARIABLE_BANDS_IMAGE

#if defined(SINGLE_BAND_IMAGE)

        add( new TImageTestSuite<SingleBandImagePolicy<ImagePolicy<vigra::GrayImage> > > ("SingleBandImage"));// fertig

#endif //SINGLE_BAND_IMAGE
 
#ifdef SELECT_BAND_IMAGE
   
//        add(new TImageTestSuite<SelectBandImagePolicy<Vector2ImagePolicy<vigra::Vector2Image>, 0> > ("SelectBandImage an vigra::Vector2Image, selected band is 0"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector2ImagePolicy<vigra::Vector2Image>, 1> >("SelectBandImage an vigra::Vector2Image, selected band is 1"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 0> >("SelectBandImage an vigra::Vector3Image, selected band is 0"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 1> >("SelectBandImage an vigra::Vector3Image, selected band is 1"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector3ImagePolicy<vigra::Vector3Image>, 2> >("SelectBandImage an vigra::Vector3Image, selected band is 2"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 0> >("SelectBandImage an vigra::Vector4Image, selected band is 0"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 1> >("SelectBandImage an vigra::Vector4Image, selected band is 1"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 2> >("SelectBandImage an vigra::Vector4Image, selected band is 2"));
//         add(new TImageTestSuite<SelectBandImagePolicy<Vector4ImagePolicy<vigra::Vector4Image>, 3> >("SelectBandImage an vigra::Vector4Image, selected band is 3"));
//         add(new TImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 0> >("SelectBandImage an vigra::RGBImage, selected band is 0"));
//         add(new TImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 1> >("SelectBandImage an vigra::RGBImage, selected band is 1"));
//         add(new TImageTestSuite<SelectBandImagePolicy<RGBImagePolicy<vigra::RGBImage>, 2> >("SelectBandImage an vigra::RGBImage, selected band is 2"));

#endif SELECT_BAND_IMAGE
    }
};

int main()
{
 
    ImageTestSuite itest;

    int failed = itest.run();

    std::cout << itest.report() << std::endl;
    
    return (failed != 0);
}  
