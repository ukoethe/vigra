#ifndef VIGRA_TEST_POLICY_PARENT_HXX
#define VIGRA_TEST_POLICY_PARENT_HXX

#include "NewImHier.hxx"

/** Ist eine Oberklasse aller ImagePolicys. ImagePolicys stellen die Besonderheiten der einzelner
* Bilder dar, z.B. Vector2ImagePolicy stellt die Besonderheiten der Bilder Vector2Image, DVector2Image FVector2Image,
* bei denen die Pixel TinyVectoren sind. Zusaetzlich wird der Vektor von entsprechenden Pixeln zu Verfuegung
* gestellt, so dass die entsprechende Pixeln brauchen nicht in der Testklasse ImageTest erzeugt zu 
* werden. Da aber die Methoden factory() fast in allen Policy - Klassen gleich sind wurde deswegen die 
* abstrakte Oberklasse geschaffen, damit die Methoden brauchen nicht mehrmals auf gleiche Weise in 
* unterschiedlichen Unterklassen wiederhollt zu werden.
*/
template<class ToTestImage>
class TestPolicy
{
public:
    typedef ToTestImage                          Image;                      //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
    typedef ToTestImage                          ChildImage;                 //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
    typedef typename ToTestImage::PixelType      PixelType;
    typedef typename Image::value_type           value_type;
    typedef typename Image::value_type           child_value_type;
    typedef std::vector<value_type>              data_array_type;
    typedef std::vector<value_type>              child_data_array_type;
    
    /** factory() - Methoden sind dazu, da um den Aufruf des richtigen Konstruktors
    * zu kapseln. Insbesondere geht es darum, dass fuer Tests von SingleBandImage
    * und VariableBandsImage die Konstruktoren der Unterklassen aufgerufen werden
    * muessen, da sie selbst abstrakte Klassen sind, also man kann nicht die Instanzen
    * von Ihnen erzeugen. Bei allen anderen Bildern muessen aber die Konstruktoren 
    * der Testklasse selber aufgerufen werden. Aus diesem Grunde auch die Unterscheidung
    * zwischen Image und ChildImage. 
    */
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

#ifdef IMAGEHIERARCHY    
    static ChildImage * factory(typename ChildImage::InnerImage image)
    {
        return new ChildImage(image);
    }
#endif IMAGEHIERARCHY     

};// end of class TestPolicy

#endif // VIGRA_TEST_POLICY_PARENT_HXX

