#ifndef VECTOR_2_IMAGE_POLICY_HXX
#define VECTOR_2_IMAGE_POLICY_HXX

#include "test_policy_parent.hxx"
#include <unittest.h>

template<class Vector2ImageP>                                           // bei dem template handelt es sich um Vector2Image, der einer der folgenden Varianten einnehmen kann: FVector2Image, DVector2Image, Vector2Image
class Vector2ImagePolicy
: public TestPolicy<Vector2ImageP> 
{

public:
    typedef Vector2ImageP                       Image;                  // die zu testende Klasse
    typedef Vector2ImageP                       ChildImage;             // entspricht der zu testenden Klasse
    typedef typename Vector2ImageP::PixelType   PixelType;              // PixelType der zu testenden Klasse
    typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
    typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse
    typedef std::vector<value_type>             data_array_type;
    typedef std::vector<child_value_type>       child_data_array_type;
    typedef typename value_type::value_type     type;
 
    static data_array_type getData()
    {
        type frgb = 0.1;
        static value_type data[15];
    
        for(int i = 0; i <= 14 ; i ++)
        {
            data[i] = value_type((i+frgb), (2*i + frgb));
        }
        static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
        return data_vector;
    }
     
    static child_data_array_type getChildData()
    {
        return getData();
    }
};

struct Vector2ImageHierarchyTestSuite
: public vigra::test_suite
{
    Vector2ImageHierarchyTestSuite();
};

struct Vector2BasicImageTestSuite
: public vigra::test_suite
{
    Vector2BasicImageTestSuite();
};
#endif // VECTOR_2_IMAGE_POLICY_HXX
