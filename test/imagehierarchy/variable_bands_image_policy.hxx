#ifndef VARIABLE_BANDS_IMAGE_POLICY_HXX
#define VARIABLE_BANDS_IMAGE_POLICY_HXX

#include "NewImHier.hxx"
#include <unittest.hxx>

/**An das template wird eine Policy-Klasse der Bilder der Imagehierarchie uebergeben
*  die moeglichen sind: ImagePolicy<GrayImage>, RGBImagePolicy<RGBImage>,
*  Vector2ImagePolicy<Vector2Image>, Vector3ImagePolicy<Vector3Image>, Vector4ImagePolicy<Vector4Image>.
*/
template<class IMAGEPOLICY>
class VariableBandsImagePolicy
{

/** Bei dem Image (in diesem Fall VariableBandsImage) handelt es sich um die zu testende Klasse
*  dagegen ChildImage ist die abgeleitete Klasse von Image, z.B. GrayImage oder FixedRGBImage.
*  Die VariableBandsImage-Konstruktoren sind protected und somit nicht aufrufbar. Um die Klasse zu testen,
*  muessen wir die Konstruktoren der abgeleiteten Klassen benutzen. Somit werden alle zur Verfuegung
*  stehende Funktionen getestet. 
*  In der Testklasse werden die ChildImage Konstruktoren aufgerufen, aber an einem Objekt
*  der VariableBandsImage Klasse (Polymorphie).
*/

public: 
    typedef typename IMAGEPOLICY::Image                 ChildImage;             // abgeleitete Klasse der ImageKlasse, es kann sich hierbei um GrayImage, SelectBandImage, SingleBandImage, FixedRGBImage oder beliebiges FixedBandImage (Vector2Image, Vector3Image usw.) handeln                
    typedef vigra::VariableBandsImage                   Image;                  // VariableBandsImage
    typedef vigra::VariableBandsImage::PixelType        PixelType;              // VektorProxy
    typedef typename ChildImage::value_type             value_type;             // 
    typedef typename ChildImage::value_type             child_value_type;       
    typedef std::vector<value_type>                     data_array_type;
    typedef std::vector<value_type>                     child_data_array_type;
    
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

struct VariableImagesTestSuite
: public vigra::test_suite
{
    VariableImagesTestSuite();
};

#endif // VARIABLE_BANDS_IMAGE_POLICY_HXX
