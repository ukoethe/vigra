#ifndef SINGLE_BAND_IMAGE_POLICY_HXX
#define SINGLE_BAND_IMAGE_POLICY_HXX

#include "NewImHier.hxx"
#include <unittest.hxx>

template<class IMAGEPOLICY>                                 // Bei der IMAGEPOLICY handelt es sich bis jetzt nur um GrayImage, spaeter soll SelectBandImage folgen
class SingleBandImagePolicy
{

/* Bei dem Image (in diesem Fall VariableBandsImage) handelt es sich um die zu testende Klasse
*  dagegen ChildImage ist die abgeleitete Klasse von Image, z.B. GrayImage oder SelectBandImage.
*  Die VariableBandsImage-Konstruktoren sind protected und somit nicht aufrufbar. Um dei Klasse zu testen,
*  muessen wir die Konstruktoren der abgeleiteten Klassen benutzen. Somit werden alle zur Verfuegung
*  stehende Funktionen getestet. 
*  In der Testklasse werden die ChildImage Konstruktoren aufgerufen, aber an einem Objekt
*  der Image Klasse (Polymorphie).
*/

public: 
    typedef typename IMAGEPOLICY::Image             ChildImage;             // abgeleitete Klasse der ImageKlasse. Es kann sich hierbei um GrayImage oder SelectBandImage handeln                  
    typedef vigra::SingleBandImage                  Image;                  // entspricht dem SingleBandImage
    typedef vigra::SingleBandImage::PixelType       PixelType;              // --||-- GrayValue (float)
    typedef typename ChildImage::value_type         value_type;             // --||-- GrayValue (float) 
    typedef typename ChildImage::value_type         child_value_type;       // --||-- GrayValue (float)  
    typedef std::vector<value_type>                 data_array_type;        // Vektor von GrayValues (float)
    typedef std::vector<value_type>                 child_data_array_type;  // Vektor von GrayValues (float)
    
    static  data_array_type getData()
    {
        return IMAGEPOLICY::getData();
    }
    
    static  child_data_array_type getChildData()
    {
        return IMAGEPOLICY::getChildData();
    }
    
    static ChildImage * factory()
    {
        return new ChildImage();
    }
    
    static ChildImage * factory(int x, int y)
    {
        return new ChildImage(x,y);
    }     
    
    static ChildImage * factory(vigra::Diff2D size)
    {
        return new ChildImage(size);
    }

    static ChildImage * factory(int x, int y, value_type pixel)
    {
        return new ChildImage(x, y, pixel);
    }
    
    static ChildImage * factory(ChildImage image)
    {
        return new ChildImage(image);
    }
    
    static ChildImage * factory(typename ChildImage::InnerImage image)
    {
        return new ChildImage(image);
    }
};

struct SingleImageTestSuite
: public vigra::test_suite
{
    SingleImageTestSuite();
};

#endif // SINGLE_BAND_IMAGE_POLICY_HXX
