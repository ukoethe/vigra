#ifndef SELECT_BAND_IMAGE_POLICY_HXX
#define SELECT_BAND_IMAGE_POLICY_HXX

#include "NewImHier.hxx"
#include <unittest.h>

template <class MULTI_BAND_IMAGE_POLICY, int TO_SELECT_BAND>
class SelectBandImagePolicy
{
public:
    typedef vigra::SelectBandImage                          Image;                  // entspricht SelectBandImage
    typedef typename MULTI_BAND_IMAGE_POLICY::ChildImage    ChildImage;             // --||-- einer MultiBandKlasse, deren Obiekt an den Konstruktor der SelectBandImage Klasse uebergeben wurde 
    typedef vigra::SelectBandImage::PixelType               PixelType;              // --||-- GrayValue (float)
    typedef vigra::SelectBandImage::value_type              value_type;             // --||-- GrayValue (float)
    typedef typename ChildImage::value_type                 child_value_type;       // --||-- VectorProxy oder RGBValue, die beiden sind letzendlich TinyVektoren 
    typedef std::vector<value_type>                         data_array_type;        // Vektor von GrayValues (float)
    typedef std::vector<child_value_type>                   child_data_array_type;  // Vektor von VectorProxy oder RGBValue, die beiden sind letzendlich TinyVektoren
    
    static const int n = TO_SELECT_BAND;                                            // der Band der selektiert werden muss

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
        return MULTI_BAND_IMAGE_POLICY::getChildData();
    }
    
    static Image * factory()
    {
        return new Image(ChildImage(), TO_SELECT_BAND);
    }
    
    static Image * factory(int x, int y)
    {
        return new Image(ChildImage(x,y), TO_SELECT_BAND);
    }     
    
    static Image * factory(vigra::Diff2D size)
    {
        return new Image(ChildImage(size), TO_SELECT_BAND);
    }

    static Image * factory(int x, int y, child_value_type pixel)
    {
        return new Image(ChildImage(x, y, pixel), TO_SELECT_BAND);
    }
    
    static Image * factory(Image image)
    {
        return new Image(image);
    }
    
    static Image * factory(ChildImage image)
    {
        return new Image(ChildImage (image), TO_SELECT_BAND);
    }
    
    static Image * factory(typename ChildImage::InnerImage image)
    {
        return new Image(ChildImage(image), TO_SELECT_BAND);
    }
};

struct SelectImagesTestSuite
: public vigra::test_suite
{
    SelectImagesTestSuite();
};



#endif // SELECT_BAND_IMAGE_POLICY_HXX
