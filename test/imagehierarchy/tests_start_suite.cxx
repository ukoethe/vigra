#   define __STL_BEGIN_RELOPS_NAMESPACE namespace std { namespace rel_ops {
#   define __STL_END_RELOPS_NAMESPACE }}
#   define __STD_RELOPS std::rel_ops
 
#include <stl_relops.h>

#   undef __STL_BEGIN_RELOPS_NAMESPACE 
#   undef __STL_END_RELOPS_NAMESPACE 
#   undef __STD_RELOPS

#include <unittest.h>
#include "one_band_image_policy.hxx"
#include "vector_2_image_policy.hxx"
#include "vector_3_image_policy.hxx"
#include "vector_4_image_policy.hxx"
#include "variable_bands_image_policy.hxx"
#include "rgb_images_policy.hxx"
#include "single_band_image_policy.hxx"
#include "select_band_image_policy.hxx"

struct TestsStartSuite
: public vigra::test_suite
{
    TestsStartSuite()
    : vigra::test_suite("ALL TESTS")
    {
         add( new VariableImagesTestSuite );
         add( new RGBImageHierarchyTestSuite );
         add( new RGBBasicImageTestSuite );
         add( new SelectImagesTestSuite );
         add( new SingleImageTestSuite );
         add( new GrayImageTestSuite );
         add( new OneBandBasicImageTestSuite );
         add( new Vector2ImageHierarchyTestSuite );
         add( new Vector3ImageHierarchyTestSuite );
         add( new Vector4ImageHierarchyTestSuite );
         add( new Vector2BasicImageTestSuite );
         add( new Vector3BasicImageTestSuite );
         add( new Vector4BasicImageTestSuite );
    }
};

int main()
{
 
    TestsStartSuite itest;

    int failed = itest.run();

    std::cout << itest.report() << std::endl;
    
    return (failed != 0);
}  
