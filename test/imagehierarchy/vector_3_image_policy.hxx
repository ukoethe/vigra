#ifndef VECTOR_3_IMAGE_POLICY_HXX
#define VECTOR_3_IMAGE_POLICY_HXX

#include "test_policy_parent.hxx"
#include <unittest.h>

template<class Vector3ImageP>                                           // bei dem template handelt es sich um Vector3Image, der einer der folgenden Varianten einnehmen kann: FVector3Image, DVector3Image, Vector3Image
class Vector3ImagePolicy
: public TestPolicy<Vector3ImageP>
{

public:
    typedef Vector3ImageP                       Image;                  // die zu testende Klasse
    typedef Vector3ImageP                       ChildImage;             // entspricht der zu testenden Klasse
    typedef typename Vector3ImageP::PixelType   PixelType;              // PixelType der zu testenden Klasse
    typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
    typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse 
    typedef std::vector<value_type>             data_array_type;
    typedef std::vector<value_type>             child_data_array_type;
    typedef typename value_type::value_type     type;
    
    static data_array_type getData()
    {
        type frgb = 0.1;
        static value_type data[15];
    
        for(int i = 0; i <= 14 ; i ++)
        {
            data[i] = value_type((i+frgb), (2*i), (2*i + frgb));
        }
        static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
        return data_vector;
    }
    
    static child_data_array_type getChildData()
    {        
        return getData();
    }
};// end of class Vector3ImagePolicy

struct Vector3ImageHierarchyTestSuite
: public vigra::test_suite
{
    Vector3ImageHierarchyTestSuite();
};


struct Vector3BasicImageTestSuite
: public vigra::test_suite
{
    Vector3BasicImageTestSuite();
};


#endif // VECTOR_3_IMAGE_POLICY_HXX
