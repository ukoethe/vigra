#ifndef VECTOR_4_IMAGE_POLICY_HXX
#define VECTOR_4_IMAGE_POLICY_HXX

#include "test_policy_parent.hxx"
#include <unittest.hxx>

template<class Vector4ImageP>                                           // bei dem template handelt es sich um Vector4Image, der einer der folgenden Varianten einnehmen kann: FVector4Image, DVector4Image, Vector4Image
class Vector4ImagePolicy
: public TestPolicy<Vector4ImageP>
{

public:
    typedef Vector4ImageP                       Image;                   // die zu testende Klasse
    typedef Vector4ImageP                       ChildImage;              // entspricht der zu testenden Klasse
    typedef typename Vector4ImageP::PixelType   PixelType;               // PixelType der zu testenden Klasse
    typedef typename Image::value_type          value_type;              // value_type der zu testenden Klasse
    typedef typename Image::value_type          child_value_type;        // entspricht dem value_type der zu testenden Klasse
    typedef std::vector<value_type>             data_array_type;
    typedef std::vector<value_type>             child_data_array_type;
    typedef typename value_type::value_type     type;
 
    static data_array_type getData()
    {
        type frgb = 0.1;
        static value_type data[15];
    
        for(int i = 0; i <= 14 ; i ++)
        {
            data[i] = value_type((i+frgb), (2*i), (2*i + frgb), (2*i + 2*frgb));
        }
        static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
        return data_vector;
    }
    
    static child_data_array_type getChildData()
    {        
        return getData();
    }
};

struct Vector4BasicImageTestSuite
: public vigra::test_suite
{
    Vector4BasicImageTestSuite();
};


struct Vector4ImageHierarchyTestSuite
: public vigra::test_suite
{
    Vector4ImageHierarchyTestSuite();
};


#endif // VECTOR_4_IMAGE_POLICY_HXX
