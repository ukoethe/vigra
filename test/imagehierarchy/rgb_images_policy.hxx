#ifndef RGB__IMAGES_POLICY_HXX
#define RGB__IMAGES_POLICY_HXX

#include "test_policy_parent.hxx"
#include <unittest.hxx>

// bei dem template handelt es sich um ein RGBImage
template<class ImageP>                                                  
class RGBImagePolicy
: public TestPolicy<ImageP> 
{

public:
    typedef ImageP Image;
    typedef ImageP ChildImage;
    typedef typename ImageP::PixelType PixelType;
    typedef typename Image::value_type value_type;  
    typedef typename Image::value_type child_value_type;  
    typedef std::vector<value_type> data_array_type;
    typedef std::vector<value_type> child_data_array_type;
    typedef typename value_type::value_type type;
    
    static data_array_type getData()
    {
        static type variable = 1;
        variable = variable/100 + variable*5/1000 + variable;
        
        static value_type data[15];
        
        for(int i = 0; i <= 14 ; i ++)
        {
            data[i] = value_type((i+variable), (2*i), (2*i + variable));
        }
        
        static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
        return data_vector;
    }
    
    static child_data_array_type getChildData()
    {   
        return getData();
    }
};

struct RGBBasicImageTestSuite
: public vigra::test_suite
{
    RGBBasicImageTestSuite();
};


struct RGBImageHierarchyTestSuite
: public vigra::test_suite
{
    RGBImageHierarchyTestSuite();
};

#endif // RGB__IMAGES_POLICY_HXX
