#ifndef ONE_BAND_IMAGE_POLICY_HXX
#define ONE_BAND_IMAGE_POLICY_HXX

#include "test_policy_parent.hxx"

template<class ImageP>                                                  // bei dem template handelt es sich um ein EinBandImage, wie FImage, BImage, GrayImage usw.
class OneBandImagePolicy 
: public TestPolicy<ImageP> 
{

public:
    typedef ImageP                               Image;                      //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
    typedef ImageP                               ChildImage;                 //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
    typedef typename ImageP::PixelType           PixelType;
    typedef typename Image::value_type           value_type;
    typedef typename Image::value_type           child_value_type;
    typedef std::vector<value_type>              data_array_type;
    typedef std::vector<value_type>              child_data_array_type;
    
    static data_array_type getData()
    {
        
        static value_type variable = 1;
        
        variable = variable/100 + (variable*5)/1000 + variable;
        
        static value_type data[15];
        
        for(int i = 0; i<=14; i++)
        {
            data[i] = variable + i;
        }
        static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
        return data_vector;
    }
    
    static child_data_array_type getChildData()
    {   
        return getData();
    }
};

struct OneBandBasicImageTestSuite 
: public vigra::test_suite
{
    OneBandBasicImageTestSuite();
};

struct GrayImageTestSuite
: public vigra::test_suite
{
    GrayImageTestSuite();
};
#endif // ONE_BAND_IMAGE_POLICY_HXX
