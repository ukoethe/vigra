#include <iostream>
#include <unittest.h>
#include <functional>
#include <algorithm>
#include "vigra/stdimage.hxx"
#include "vigra/rgbvalue.hxx"
#include "NewImHier.hxx"

#if defined(GRAY_IMAGE) || defined(nRGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE)
    #include "one_band_image_policy.hxx"
#endif

#if defined(FIXED_RGB_IMAGE) || defined(RGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE) || defined(SELECT_BAND_IMAGE)
    #include "rgb_images_policy.hxx"
#endif

#if defined(FIXED_V_IMAGE)  ||defined(VImage)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE) || defined (SELECT_BAND_IMAGE)
    #include "vector_2_image_policy.hxx"
    #include "vector_3_image_policy.hxx"
    #include "vector_4_image_policy.hxx"
#endif

#ifdef VARIABLE_BANDS_IMAGE
    #include "variable_bands_image_policy.hxx"
#endif

#ifdef SINGLE_BAND_IMAGE
    #include "single_band_image_policy.hxx"
#endif

#ifdef SELECT_BAND_IMAGE
    #include "select_band_image_policy.hxx"
#endif 
using vigra::Diff2D;

std::ostream & operator<<(std::ostream & o, vigra::ConstVectorProxy const & p)
{
    if(p.size() == 1)
    {
        o << p[0];
    }
    else
    {
        o << "(";
        int i;
        for (i=0; i<p.size()-1; ++i)
            o << p[i] << ", ";
        o << p[i] << ")";
    }
    return o;
}

/** Vergleicht zwei Bilder auf aequivalenz der Pixel. Die Pixel muessen die GLEICHE 
* ANZAHL der Baender haben
*/
template <class _InputIter1, class _InputIter2>
inline bool equalPixels(_InputIter1 __first1, _InputIter1 __last1, _InputIter2 __first2) 
{
  for ( ; __first1 != __last1; ++__first1, ++__first2)
    if (!(*__first1 == *__first2))
      return false;
  return true;
}

/** Vergleicht den gleichen Durchlauf von Iterator und ScanOrderIterator. Die Pixel muessen die GLEICHE 
* ANZAHL der Baender haben
*/
template <class Iterator, class ScanOrderIterator>
inline bool equalIteratorRun(Iterator ul, Iterator lr, ScanOrderIterator f) 
{
  for(; ul.y < lr.y; ++ul.y)
  {
    Iterator c = ul;
    for(; c.x < lr.x; ++c.x, ++f)
    {
        if(*c != *f)
            return false;
    }
  }
  return true;
}

/** Untersucht zwei vergleichbare Bider auf Groesse- und PixelIdentitaet
*  aeusserste Vorsicht: die Bilder muessen vergleichbar sein!!!! Es sind Bilder
*  der gleichen Klasse vergleichbar und Bilder mit gleichen Anzahl der Baender !!!
*/

template<class Image, class ComparableImage>
bool equal(Image const & image1, ComparableImage const & image2)
{
    if (image1.size() == image2.size())
        return image1.height() == image2.height() &&
               image1.width() == image2.width() &&
               equalPixels(image1.begin(), image1.end(), image2.begin());
    else return false;
};

/** Entspricht der "normalen" binary_function aus std, wurde aber aus Kompatibilitaetsgruenden
* umbenannt und hier mitaufgenommen.
*/
template <class _Arg1, class _Arg2, class _Result>
struct Pixels_binary_function {
  typedef _Arg1 first_argument_type;
  typedef _Arg2 second_argument_type;
  typedef _Result result_type;
}; 

/** Vergleicht zwei Pixel auf Gleichheit.Entspricht der "normalen" not_equal_to aus std, 
* wurde aber aus Kompatibilitaetsgruenden umbenannt und hier mitaufgenommen.
*/
template <class _Tp>
struct Pixels_not_equal_to : public Pixels_binary_function<_Tp,_Tp,bool> 
{
  bool operator()(const _Tp& __x, const _Tp& __y) const { return !(__x == __y); }
};

// /** Gibt den Inhalt eines Vektors auf der Konsole aus
// */
// template <class Vector>
// void printVector(Vector to_print_vector)
// {
//     for (int i = 0; i < to_print_vector.size() ; i++ )
//         cout <<"[" << i << "]: " << to_print_vector[i] << ", ";
//     cout << endl;
// }

// /** Ist eine Oberklasse aller ImagePolicys. ImagePolicys stellen die Besonderheiten der einzelner
// * Bilder dar, z.B. Vector2ImagePolicy stellt die Besonderheiten der Bilder Vector2Image, DVector2Image FVector2Image,
// * bei denen die Pixel TinyVectoren sind. Zusaetzlich wird der Vektor von entsprechenden Pixeln zu Verfuegung
// * gestellt, so dass die entsprechende Pixeln brauchen nicht in der Testklasse ImageTest erzeugt zu 
// * werden. Da aber die Methoden factory() fast in allen Policy - Klassen gleich sind wurde deswegen die 
// * abstrakte Oberklasse geschaffen, damit die Methoden brauchen nicht mehrmals auf gleiche Weise in 
// * unterschiedlichen Unterklassen wiederhollt zu werden.
// */
// template<class ToTestImage>
// class TestPolicy
// {
// public:
//     typedef ToTestImage                          Image;                      //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef ToTestImage                          ChildImage;                 //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef typename ToTestImage::PixelType      PixelType;
//     typedef typename Image::value_type           value_type;
//     typedef typename Image::value_type           child_value_type;
//     typedef std::vector<value_type>              data_array_type;
//     typedef std::vector<value_type>              child_data_array_type;
//     
//     /** factory() - Methoden sind dazu, da um den Aufruf des richtigen Konstruktors
//     * zu kapseln. Insbesondere geht es darum, dass fuer Tests von SingleBandImage
//     * und VariableBandsImage die Konstruktoren der Unterklassen aufgerufen werden
//     * muessen, da sie selbst abstrakte Klassen sind, also man kann nicht die Instanzen
//     * von Ihnen erzeugen. Bei allen anderen Bildern muessen aber die Konstruktoren 
//     * der Testklasse selber aufgerufen werden. Aus diesem Grunde auch die Unterscheidung
//     * zwischen Image und ChildImage. 
//     */
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
// 
// #ifdef IMAGEHIERARCHY    
//     static ChildImage * factory(typename ChildImage::InnerImage image)
//     {
//         return new ChildImage(image);
//     }
// #endif IMAGEHIERARCHY     
//  
// };

//#if defined(GRAY_IMAGE) || defined(nRGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE)
// template<class ImageP>                                                  // bei dem template handelt es sich um ein EinBandImage, wie FImage, BImage, GrayImage usw.
// class ImagePolicy
// {
// 
// public:
//     typedef ImageP                          Image;                      //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef ImageP                          ChildImage;                 //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef typename ImageP::PixelType      PixelType;
//     typedef typename Image::value_type      value_type;  
//     typedef typename Image::value_type      child_value_type;  
//     typedef std::vector<value_type>         data_array_type;
//     typedef std::vector<value_type>         child_data_array_type;
//     
//     
//     static data_array_type getData()
//     {
//         
//         static value_type variable = 1;
//         
//         variable = variable/100 + (variable*5)/1000 + variable;
//         
//         static value_type data[15];
//         
//         for(int i = 0; i<=14; i++)
//         {
//             data[i] = variable + i;
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {   
//         return getData();
//     }
//     
//     /** factory() - Methoden sind dazu da um den Aufruf des richtigen Konstruktors
//     * zu kapseln. Insbesondere geht es darum, dass fuer Tests von SingleBandImage
//     * und VariableBandsImage die Konstruktoren der Unterklassen aufgerufen werden
//     * muessen, da sie selbst abstrakte Klassen sind, also man kann nicht die Instanzen
//     * von Ihnen erzeugen. 
//     */
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
// 
// #ifdef IMAGEHIERARCHY    
//     static ChildImage * factory(typename ChildImage::InnerImage image)
//     {
//         return new ChildImage(image);
//     }
// #endif IMAGEHIERARCHY     
// };
// template<class ImageP>                                                  // bei dem template handelt es sich um ein EinBandImage, wie FImage, BImage, GrayImage usw.
// class ImagePolicy 
// : public TestPolicy<ImageP> 
// {
// 
// public:
//     typedef ImageP                               Image;                      //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef ImageP                               ChildImage;                 //es handelt sich um ein EinBandImage, wie FImage, BImage, GrayImage usw. ChildImage und Image sind aequivalent.
//     typedef typename ImageP::PixelType           PixelType;
//     typedef typename Image::value_type           value_type;
//     typedef typename Image::value_type           child_value_type;
//     typedef std::vector<value_type>              data_array_type;
//     typedef std::vector<value_type>              child_data_array_type;
//     
//     static data_array_type getData()
//     {
//         
//         static value_type variable = 1;
//         
//         variable = variable/100 + (variable*5)/1000 + variable;
//         
//         static value_type data[15];
//         
//         for(int i = 0; i<=14; i++)
//         {
//             data[i] = variable + i;
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {   
//         return getData();
//     }
// };
//#endif // GRAY_IMAGE ||nRGB_IMAGE

//#if defined(FIXED_RGB_IMAGE) || defined(RGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE) || defined(SELECT_BAND_IMAGE)
// template<class ImageP>                                                  // bei dem template handelt es sich um ein RGBImage
// class RGBImagePolicy
// {
// 
// public:
//     typedef ImageP Image;
//     typedef ImageP ChildImage;
//     typedef typename ImageP::PixelType PixelType;
//     typedef typename Image::value_type value_type;  
//     typedef typename Image::value_type child_value_type;  
//     typedef std::vector<value_type> data_array_type;
//     typedef std::vector<value_type> child_data_array_type;
//     typedef typename value_type::value_type type;
//     
//     static data_array_type getData()
//     {
//         static type variable = 1;
//         variable = variable/100 + variable*5/1000 + variable;
//         
//         static value_type data[15];
//         
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+variable), (2*i), (2*i + variable));
//         }
//         
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {   
//         return getData();
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
// 
//     #ifdef IMAGEHIERARCHY    
//         static ChildImage * factory(typename ChildImage::InnerImage image)
//         {
//             return new ChildImage(image);
//         }
//     #endif IMAGEHIERARCHY     
// };
// template<class ImageP>                                                  // bei dem template handelt es sich um ein RGBImage
// class RGBImagePolicy
// : public TestPolicy<ImageP> 
// {
// 
// public:
//     typedef ImageP Image;
//     typedef ImageP ChildImage;
//     typedef typename ImageP::PixelType PixelType;
//     typedef typename Image::value_type value_type;  
//     typedef typename Image::value_type child_value_type;  
//     typedef std::vector<value_type> data_array_type;
//     typedef std::vector<value_type> child_data_array_type;
//     typedef typename value_type::value_type type;
//     
//     static data_array_type getData()
//     {
//         static type variable = 1;
//         variable = variable/100 + variable*5/1000 + variable;
//         
//         static value_type data[15];
//         
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+variable), (2*i), (2*i + variable));
//         }
//         
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {   
//         return getData();
//     }
// };
//#endif //FIXED_RGB_IMAGE || RGB_IMAGE  || defined (VARIABLE_BANDS_IMAGE)

//#if defined(FIXED_V_IMAGE)  ||defined(VImage)  || defined (VARIABLE_BANDS_IMAGE) || defined(SINGLE_BAND_IMAGE) || defined (SELECT_BAND_IMAGE)
// 
// template<class Vector2ImageP>                                           // bei dem template handelt es sich um Vector2Image, der einer der folgenden Varianten einnehmen kann: FVector2Image, DVector2Image, Vector2Image
// class Vector2ImagePolicy
// {
// 
// public:
//     typedef Vector2ImageP                       Image;                  // die zu testende Klasse
//     typedef Vector2ImageP                       ChildImage;             // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector2Image::PixelType      PixelType;
// #else                                                                   // PixelType der zu testenden Klasse
//     typedef typename Vector2ImageP::PixelType   PixelType;
// #endif
//     typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<child_value_type>       child_data_array_type;
//     typedef typename value_type::value_type     type;
//  
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i + frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//      
//     static child_data_array_type getChildData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i + frgb));
//         }
//         static child_data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
//     
//     #ifdef IMAGEHIERARCHY                                               // Da diese Policy wie fuer Tests von Vector2Image der imagehierarchy, als auch von FVector2Image, DVector2Image  der BasicImage muss diese Methode auskommentiert werden, da BasicImage kein InnerImage besitzt 
//         /** factory() - Methoden sind dazu da um den Aufruf des richtigen Konstruktors
//         * zu kapseln. Insbesondere geht es darum, dass fuer Tests von SingleBandImage
//         * und VariableBandsImage die Konstruktoren der Unterklassen aufgerufen werden
//         * muessen, da sie selbst abstrakte Klassen sind, also man kann nicht die Instanzen
//         * von Ihnen erzeugen. 
//         */
//         static ChildImage * factory(typename ChildImage::InnerImage image)
//         {
//             return new ChildImage(image);
//         }
//     #endif IMAGEHIERARCHY       
// };

// template<class Vector2ImageP>                                           // bei dem template handelt es sich um Vector2Image, der einer der folgenden Varianten einnehmen kann: FVector2Image, DVector2Image, Vector2Image
// class Vector2ImagePolicy
// : public TestPolicy<Vector2ImageP> 
// {
// 
// public:
//     typedef Vector2ImageP                       Image;                  // die zu testende Klasse
//     typedef Vector2ImageP                       ChildImage;             // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector2Image::PixelType      PixelType;
// #else                                                                   // PixelType der zu testenden Klasse
//     typedef typename Vector2ImageP::PixelType   PixelType;
// #endif
//     typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<child_value_type>       child_data_array_type;
//     typedef typename value_type::value_type     type;
//  
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i + frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//      
//     static child_data_array_type getChildData()
//     {
//         return getData();
//     }
// };
// 
// template<class Vector3ImageP>                                           // bei dem template handelt es sich um Vector3Image, der einer der folgenden Varianten einnehmen kann: FVector3Image, DVector3Image, Vector3Image
// class Vector3ImagePolicy
// {
// 
// public:
//     typedef Vector3ImageP                       Image;                  // die zu testende Klasse
//     typedef Vector3ImageP                       ChildImage;             // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector3Image::PixelType      PixelType;
// #else                                                                   // PixelType der zu testenden Klasse
//     typedef typename Vector3ImageP::PixelType   PixelType;
// #endif    
//     typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse 
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<value_type>             child_data_array_type;
//     typedef typename value_type::value_type     type;
//     
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i), (2*i + frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {        
//         return getData();
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory( ChildImage image)
//     {
//         return new ChildImage(image);
//     }
//     
// #ifdef IMAGEHIERARCHY                                                   // Da diese Policy wie fuer Tests  der  Klasse imagehierarchy, als auch der BasicImage muss diese Methode auskommentiert werden, da BasicImage kein InnerImage besitzt    
//     static ChildImage * factory(typename ChildImage::InnerImage image)
//     {
//         return new ChildImage(image);
//     }
// #endif IMAGEHIERARCHY      
// };

// template<class Vector3ImageP>                                           // bei dem template handelt es sich um Vector3Image, der einer der folgenden Varianten einnehmen kann: FVector3Image, DVector3Image, Vector3Image
// class Vector3ImagePolicy
// : public TestPolicy<Vector3ImageP>
// {
// 
// public:
//     typedef Vector3ImageP                       Image;                  // die zu testende Klasse
//     typedef Vector3ImageP                       ChildImage;             // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector3Image::PixelType      PixelType;
// #else                                                                   // PixelType der zu testenden Klasse
//     typedef typename Vector3ImageP::PixelType   PixelType;
// #endif    
//     typedef typename Image::value_type          value_type;             // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;       // entspricht dem value_type der zu testenden Klasse 
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<value_type>             child_data_array_type;
//     typedef typename value_type::value_type     type;
//     
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i), (2*i + frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {        
//         return getData();
//     }
// };// end of class Vector3ImagePolicy
// 
// template<class Vector4ImageP>                                           // bei dem template handelt es sich um Vector4Image, der einer der folgenden Varianten einnehmen kann: FVector4Image, DVector4Image, Vector4Image
// class Vector4ImagePolicy
// {
// 
// public:
//     typedef Vector4ImageP                       Image;                   // die zu testende Klasse
//     typedef Vector4ImageP                       ChildImage;              // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector4Image::PixelType      PixelType;
// #else                                                                    // PixelType der zu testenden Klasse
//     typedef typename Vector4ImageP::PixelType   PixelType;
// #endif
//     typedef typename Image::value_type          value_type;              // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;        // entspricht dem value_type der zu testenden Klasse
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<value_type>             child_data_array_type;
//     typedef typename value_type::value_type     type;
//  
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i), (2*i + frgb), (2*i + 2*frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {        
//         return getData();
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
//     
//     #ifdef IMAGEHIERARCHY                                                    // Da diese Policy wie fuer Tests  der  Klasse imagehierarchy, als auch der BasicImage muss diese Methode auskommentiert werden, da BasicImage kein InnerImage besitzt    
//         static ChildImage * factory(typename ChildImage::InnerImage image)
//         {
//             return new ChildImage(image);
//         }
//     #endif IMAGEHIERARCHY       
// };
// 
// template<class Vector4ImageP>                                           // bei dem template handelt es sich um Vector4Image, der einer der folgenden Varianten einnehmen kann: FVector4Image, DVector4Image, Vector4Image
// class Vector4ImagePolicy
// : public TestPolicy<Vector4ImageP>
// {
// 
// public:
//     typedef Vector4ImageP                       Image;                   // die zu testende Klasse
//     typedef Vector4ImageP                       ChildImage;              // entspricht der zu testenden Klasse
// #ifndef VImage
//     typedef vigra::Vector4Image::PixelType      PixelType;
// #else                                                                    // PixelType der zu testenden Klasse
//     typedef typename Vector4ImageP::PixelType   PixelType;
// #endif
//     typedef typename Image::value_type          value_type;              // value_type der zu testenden Klasse
//     typedef typename Image::value_type          child_value_type;        // entspricht dem value_type der zu testenden Klasse
//     typedef std::vector<value_type>             data_array_type;
//     typedef std::vector<value_type>             child_data_array_type;
//     typedef typename value_type::value_type     type;
//  
//     static data_array_type getData()
//     {
//         type frgb = 0.1;
//         static value_type data[15];
//     
//         for(int i = 0; i <= 14 ; i ++)
//         {
//             data[i] = value_type((i+frgb), (2*i), (2*i + frgb), (2*i + 2*frgb));
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {        
//         return getData();
//     }
// };
//#endif // FIXED_V_IMAGE ||VImage  || defined (VARIABLE_BANDS_IMAGE)

//#ifdef VARIABLE_BANDS_IMAGE

// /**An das template wird eine Policy-Klasse der Bilder der Imagehierarchie uebergeben
// *  die moeglichen sind: ImagePolicy<GrayImage>, RGBImagePolicy<RGBImage>,
// *  Vector2ImagePolicy<Vector2Image>, Vector3ImagePolicy<Vector3Image>, Vector4ImagePolicy<Vector4Image>.
// */
// template<class IMAGEPOLICY>
// class VariableBandsImagePolicy
// {
// 
// /** Bei dem Image (in diesem Fall VariableBandsImage) handelt es sich um die zu testende Klasse
// *  dagegen ChildImage ist die abgeleitete Klasse von Image, z.B. GrayImage oder FixedRGBImage.
// *  Die VariableBandsImage-Konstruktoren sind protected und somit nicht aufrufbar. Um die Klasse zu testen,
// *  muessen wir die Konstruktoren der abgeleiteten Klassen benutzen. Somit werden alle zur Verfuegung
// *  stehende Funktionen getestet. 
// *  In der Testklasse werden die ChildImage Konstruktoren aufgerufen, aber an einem Objekt
// *  der VariableBandsImage Klasse (Polymorphie).
// */
// 
// public: 
//     typedef typename IMAGEPOLICY::Image                 ChildImage;             // abgeleitete Klasse der ImageKlasse, es kann sich hierbei um GrayImage, SelectBandImage, SingleBandImage, FixedRGBImage oder beliebiges FixedBandImage (Vector2Image, Vector3Image usw.) handeln                
//     typedef vigra::VariableBandsImage                   Image;                  // VariableBandsImage
//     typedef vigra::VariableBandsImage::PixelType        PixelType;              // VektorProxy
//     typedef typename ChildImage::value_type             value_type;             // 
//     typedef typename ChildImage::value_type             child_value_type;       
//     typedef std::vector<value_type>                     data_array_type;
//     typedef std::vector<value_type>                     child_data_array_type;
//     
//     static  data_array_type getData()
//     {
//         return IMAGEPOLICY::getData();
//     }
//     
//     static  child_data_array_type getChildData()
//     {
//         return IMAGEPOLICY::getChildData();
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
//     
//     static ChildImage * factory(typename ChildImage::InnerImage image)
//     {
//         return new ChildImage(image);
//     }
// 
// };
//#endif //VARIABLE_BANDS_IMAGE

//#ifdef SINGLE_BAND_IMAGE

// template<class IMAGEPOLICY>                                 // Bei der IMAGEPOLICY handelt es sich bis jetzt nur um GrayImage, spaeter soll SelectBandImage folgen
// class SingleBandImagePolicy
// {
// 
// /* Bei dem Image (in diesem Fall VariableBandsImage) handelt es sich um die zu testende Klasse
// *  dagegen ChildImage ist die abgeleitete Klasse von Image, z.B. GrayImage oder SelectBandImage.
// *  Die VariableBandsImage-Konstruktoren sind protected und somit nicht aufrufbar. Um dei Klasse zu testen,
// *  muessen wir die Konstruktoren der abgeleiteten Klassen benutzen. Somit werden alle zur Verfuegung
// *  stehende Funktionen getestet. 
// *  In der Testklasse werden die ChildImage Konstruktoren aufgerufen, aber an einem Objekt
// *  der Image Klasse (Polymorphie).
// */
// 
// public: 
//     typedef typename IMAGEPOLICY::Image             ChildImage;             // abgeleitete Klasse der ImageKlasse. Es kann sich hierbei um GrayImage oder SelectBandImage handeln                  
//     typedef vigra::SingleBandImage                  Image;                  // entspricht dem SingleBandImage
//     typedef vigra::SingleBandImage::PixelType       PixelType;              // --||-- GrayValue (float)
//     typedef typename ChildImage::value_type         value_type;             // --||-- GrayValue (float) 
//     typedef typename ChildImage::value_type         child_value_type;       // --||-- GrayValue (float)  
//     typedef std::vector<value_type>                 data_array_type;        // Vektor von GrayValues (float)
//     typedef std::vector<value_type>                 child_data_array_type;  // Vektor von GrayValues (float)
//     
//     static  data_array_type getData()
//     {
//         return IMAGEPOLICY::getData();
//     }
//     
//     static  child_data_array_type getChildData()
//     {
//         return IMAGEPOLICY::getChildData();
//     }
//     
//     static ChildImage * factory()
//     {
//         return new ChildImage();
//     }
//     
//     static ChildImage * factory(int x, int y)
//     {
//         return new ChildImage(x,y);
//     }     
//     
//     static ChildImage * factory(Diff2D size)
//     {
//         return new ChildImage(size);
//     }
// 
//     static ChildImage * factory(int x, int y, value_type pixel)
//     {
//         return new ChildImage(x, y, pixel);
//     }
//     
//     static ChildImage * factory(ChildImage image)
//     {
//         return new ChildImage(image);
//     }
//     
//     static ChildImage * factory(typename ChildImage::InnerImage image)
//     {
//         return new ChildImage(image);
//     }
// };
//#endif //SINGLE_BAND_IMAGE

//#ifdef SELECT_BAND_IMAGE
// template <class MULTI_BAND_IMAGE_POLICY, int TO_SELECT_BAND>
// class SelectBandImagePolicy
// {
// public:
//     typedef vigra::SelectBandImage                          Image;                  // entspricht SelectBandImage
//     typedef typename MULTI_BAND_IMAGE_POLICY::ChildImage    ChildImage;             // --||-- einer MultiBandKlasse, deren Obiekt an den Konstruktor der SelectBandImage Klasse uebergeben wurde 
//     typedef vigra::SelectBandImage::PixelType               PixelType;              // --||-- GrayValue (float)
//     typedef vigra::SelectBandImage::value_type              value_type;             // --||-- GrayValue (float)
//     typedef typename ChildImage::value_type                 child_value_type;       // --||-- VectorProxy oder RGBValue, die beiden sind letzendlich TinyVektoren 
//     typedef std::vector<value_type>                         data_array_type;        // Vektor von GrayValues (float)
//     typedef std::vector<child_value_type>                   child_data_array_type;  // Vektor von VectorProxy oder RGBValue, die beiden sind letzendlich TinyVektoren
//     
//     static const int n = TO_SELECT_BAND;                                            // der Band der selektiert werden muss
// 
//     static data_array_type getData()
//     {
//         
//         static value_type variable = 1;
//         variable = variable/100 + (variable*5)/1000 + variable;
//         static value_type data[15];
//         for(int i = 0; i<=14; i++)
//         {
//             data[i] = variable + i;
//         }
//         static data_array_type data_vector(data, data+sizeof(data)/sizeof(value_type));
//         return data_vector;
//     }
//     
//     static child_data_array_type getChildData()
//     {        
//         return MULTI_BAND_IMAGE_POLICY::getChildData();
//     }
//     
//     static Image * factory()
//     {
//         return new Image(ChildImage(), TO_SELECT_BAND);
//     }
//     
//     static Image * factory(int x, int y)
//     {
//         return new Image(ChildImage(x,y), TO_SELECT_BAND);
//     }     
//     
//     static Image * factory(Diff2D size)
//     {
//         return new Image(ChildImage(size), TO_SELECT_BAND);
//     }
// 
//     static Image * factory(int x, int y, child_value_type pixel)
//     {
//         return new Image(ChildImage(x, y, pixel), TO_SELECT_BAND);
//     }
//     
//     static Image * factory(Image image)
//     {
//         return new Image(image);
//     }
//     
//     static Image * factory(ChildImage image)
//     {
//         return new Image(ChildImage (image), TO_SELECT_BAND);
//     }
//     
//     static Image * factory(typename ChildImage::InnerImage image)
//     {
//         return new Image(ChildImage(image), TO_SELECT_BAND);
//     }
// };
//#endif //SELECT_BAND_IMAGE

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////                                    TESTKLASSE                                 ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class Policy>
class ImageTest
{
public:
    typedef typename Policy::Image              Image;              // zu testende Klasse z.B. GrayImage, VariableBandsImage usw.
    typedef typename Policy::ChildImage         ChildImage;         // unterscheidet sich von zu testender Klasse Image, nur wenn einer der abstrakten Klassen VariableBandImage oder SingleBandImage oder die Klasse SelectBandImage getestet wird, sonst entspricht es der Klasse Image.  
    typedef typename Policy::value_type         value_type;         // value_type der zu testenden Klasse
    typedef typename Policy::child_value_type   child_value_type;   // value_type der Klasse an der die zu testende Klasse getestet wird, also z.B. VariableBandsImage wird am Beispiel von Vector2Image getestet dann ist child_value_type die value_type von Vector2Image
    
    typename Policy::data_array_type data;
    typename Policy::child_data_array_type child_data;
 
    std::auto_ptr<Image> image0_;
    std::auto_ptr<Image> image1_;
    
    typename Image::Accessor acc0_;
    typename Image::Accessor acc1_;
    
    typename Image::ScanOrderIterator scanOrderIter0_;
    typename Image::ConstScanOrderIterator constScanOrderIter0_;
    typename Image::ScanOrderIterator scanOrderIter1_;
    typename Image::ConstScanOrderIterator constScanOrderIter1_;
    
    typename Image::traverser traverserIter0_;
    typename Image::const_traverser constTraverserIter0_;
    typename Image::traverser traverserIter1_;
    typename Image::const_traverser constTraverserIter1_;

    typename Image::traverser::row_iterator rowIter1_;
    typename Image::const_traverser::row_iterator constRowIter1_;

    typename Image::traverser::column_iterator columnIter1_;
    typename Image::const_traverser::column_iterator constColumnIter1_;
    
    ImageTest()
    : image0_(Policy::factory()),
      image1_(Policy::factory(3,5)), 
      data(Policy::getData()),
      child_data(Policy::getChildData())
    {
        acc0_ = image0_->accessor();
        acc1_ = image1_->accessor();

        scanOrderIter0_ = image0_->begin();
        scanOrderIter1_ = image1_->begin();

        acc1_.set(data[0], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[1], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[2], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[3], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[4], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[5], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[6], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[7], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[8], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[9], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[10], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[11], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[12], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[13], scanOrderIter1_);
        ++scanOrderIter1_;
        acc1_.set(data[14], scanOrderIter1_);
        ++scanOrderIter1_;

        should(scanOrderIter1_ == image1_->end());
    }
    
    /**  Testet den Default-Konstruktor der zu testenden Imageklasse
    */
    void testImageDefaultConstuctor()
    {
        /*
        *  image0_ wurde mit dem Default-Konstruktor erzeugt, es hat die Groesse Diff2D(0, 0)
        */
        should(image0_->height() == 0);
        should(image0_->width() == 0);
        should(image0_->size() == vigra::Diff2D(0,0));
        should(!image0_->isInside(vigra::Diff2D(-1,0)));
        should(!image0_->isInside(vigra::Diff2D(1,0)));
        should(!image0_->isInside(vigra::Diff2D(0,1)));
        should(!image0_->isInside(vigra::Diff2D(0,-1)));
    }
    
    /**  Testet den "(WIDTH, HEIGHT)"-Konstruktor der zu testenden Imageklasse
    */
    void testImageIntConstuctor()
    {
        should(image1_->height() == 5);
        should(image1_->width() == 3);
        should(image1_->size() == vigra::Diff2D(3,5));
    }
    
    /**  Testet den "(Diff2D)"-Konstruktor der zu testenden Imageklasse
    */
    void testImage2DConstuctor()
    {
        std::auto_ptr<Image> image1(Policy::factory(Diff2D(2,3)));
        should(image1->height() == 3);
        should(image1->width() == 2);
        should(image1->size() == vigra::Diff2D(2,3));
        
        std::auto_ptr<Image> image2(Policy::factory(Diff2D(0, 0)));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
    }
    
    /**  Testet den "(WIDTH, HEIGHT, PIXELTYPE)"-Konstruktor der zu testenden Imageklasse
    */
    void testImageIntPixelConstuctor()
    {
        std::auto_ptr<Image> image1(Policy::factory(2,3, child_data[0]));
        should(image1->height() == 3);
        should(image1->width() == 2);
        should(image1->size() == vigra::Diff2D(2,3));
    
    #ifdef SELECT_BAND_IMAGE                                                // Bei SelectBandImage wird nur der selektierte Band mit dem Pixel child_data[0] initialisiert   
        should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[0][Policy::n])));
    #else
        should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[0])));    
    #endif SELECT_BAND_IMAGE
        
        std::auto_ptr<Image> image2(Policy::factory(0, 0, child_data[1]));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
    
    #ifdef SELECT_BAND_IMAGE                                                // Bei SelectBandImage wird nur der selektierte Band mit dem Pixel child_data[1] initialisiert   
        should(image2->end() == find_if(image2->begin(), image2->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[1][Policy::n])));
    #else
        should(image2->end() == find_if(image2->begin(), image2->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[1])));    
    #endif SELECT_BAND_IMAGE
    }
    
    /** Testet den Copy Konstruktor ( Image(Image img) ).
    */
    void testCopyConstructor()
    {
        
        #ifndef SELECT_BAND_IMAGE 
            ChildImage * image1 = new ChildImage(3, 5, child_data[2]);
            std::auto_ptr<Image> image1_copy(Policy::factory(*image1));
            
            ChildImage image0(0, 0, child_data[3]);
            std::auto_ptr<Image> image0_copy(Policy::factory(image0));
        #else
            std::auto_ptr<Image> image1_copy(Policy::factory(*image1_));
            std::auto_ptr<Image> image1 = image1_;
            
            std::auto_ptr<Image> image0_copy(Policy::factory(*image0_));
        #endif  SELECT_BAND_IMAGE 

            should(image1_copy->height() == 5);
            should(image1_copy->width() == 3);
            should(image1_copy->size() == vigra::Diff2D(3,5));
            should(equalPixels(image1->begin(),image1->end(), image1_copy->begin()));
            should(image0_copy->height() == 0);
            should(image0_copy->width() == 0);
            should(image0_copy->size() == vigra::Diff2D(0,0));
            should(equalPixels(image0_->begin(), image0_->end(), image0_copy->begin()));
     }
    
//     /** testet die Konstruktorden der imagehierarchy Klasse, an die die InnerImage (BasicImage<float>)
//     * uebergeben wird, also Image (InnerImage i)
//     */
//     void testInnerImageConstructor()
//     {
//         std::auto_ptr<Image> image1(Policy::factory(new typename ChildImage::InnerImage(3, 4, child_data[0]))); 
//         should(image1->height() == 4);
//         should(image1->width() == 3);
//         should(image1->size() == vigra::Diff2D(3,4));
//     #ifdef SELECT_BAND_IMAGE    
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[0][Policy::n])));
//     #else
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[0])));    
//     #endif SELECT_BAND_IMAGE    
//         
//         std::auto_ptr<Image> image2(Policy::factory(new typename ChildImage::InnerImage(0, 0, child_data[1])));
//         should(image2->height() == 0);
//         should(image2->width() == 0);
//         should(image2->size() == vigra::Diff2D(0,0));
//     #ifdef SELECT_BAND_IMAGE    
//         should(image2->end() == find_if(image2->begin(), image2->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[1][Policy::n])));
//     #else
//         should(image2->end() == find_if(image2->begin(), image2->end(), bind2nd(Pixels_not_equal_to<value_type>(), child_data[1])));    
//     #endif SELECT_BAND_IMAGE
//     }
//     
//     /** testet die clone() Methode der Klasse aus imagehierarchy
//     */
//     void testClone()
//     {
//         /*
//         *  Im Falle der Aenderungsvornehmungen an einem der Images, werden die Aenderungen auch nur an einem sichtbar
//         */
//         std::auto_ptr<typename Image::CloneType> image1(image1_->clone()); 
//         should(equal(*image1, *image1_));                                               
//         should(equalPixels(image1_->begin(),image1_->end(), image1->begin()));             
// #ifndef SELECT_BAND_IMAGE
//         should(&(*image1) != &(*image1_));                                              
// #else        
//         should(static_cast<vigra::SingleBandImage *> (&(*image1)) != static_cast<vigra::SingleBandImage *> (&(*image1_)));
// #endif SELECT_BAND_IMAGE
//         
//         /* Aenderung mit der init-Funktion
//         */
//         image1->init(data[5]); 
//         should((*image1_->begin()) != static_cast<typename Image::PixelType> (data[5]));
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
//         
//         image1_->init(data[6]);
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
//         
//         std::auto_ptr<typename Image::CloneType> image0(image0_->clone());
//         should(equal(*image0, *image0_));
//         should(equalPixels(image0_->begin(),image0_->end(), image0->begin()));image1_->init(data[6]);
// #ifndef SELECT_BAND_IMAGE
//         should(&(*image0) != &(*image0_));
// #else  
//         should(static_cast<vigra::SingleBandImage *> (&(*image0)) != static_cast<vigra::SingleBandImage *> (&(*image0_)));
// #endif SELECT_BAND_IMAGE
// 
//         /* Aenderung mit dem Zuweisungsoperator -> siehe in testShallowCopy()
//         */
// #if !defined(VARIABLE_BANDS_IMAGE) && !defined(SINGLE_BAND_IMAGE)               //Bei VariableBandsImage kann man den Zuweisungsoperator nicht testen, da sie eine abstrakte Klasse ist
//         (*image0_) = (*image1_);                              
//         should(image0->size() == vigra::Diff2D(0,0));
//         should(image0_->size() == vigra::Diff2D(3,5));
//         should(!equal(*image0, *image0_));
//         should(equal(*image1_, *image0_));
// #endif //!VARIABLE_BANDS_IMAGE && ! SIMGLE_BAND_IMAGE   
//     }
//     
//     /** testet die shallowCopy() Methode der Klassen aus imagehierarchy.
//     */
//     void testShallowCopy()
//     {
//         /*
//         *  Im Falle der Aenderungsvornehmungen an einem der Images, werden die Aenderungen an beiden sichtbar
//         *  Es gibt nur zwei Aenderungsmoeglichkeiten init() und Zuweisung eines anderen Images
//         *  bei der Zuweisung werden die Zeiger von den eigenen Daten auf die anderen umgeleitet,
//         *  und das entspricht nicht dem Sinn der shallowCopy
//         */
//         std::auto_ptr<Image> image1(image1_->shallowCopy());
//         should(equal(*image1, *image1_));     
//         should(&(*image1) != &(*image1_));
//         
//         /* Aenderung mit der init-Funktion
//         */
//         image1->init(data[7]);
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
//         
//         image1->init(data[8]);
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
//         
//         /* Eine shallowCopy zeigt auf die selben Daten des kopierten Objektes
//         */
//         std::auto_ptr<Image> image1Copy(image1->shallowCopy());
//         should(equal(*image1Copy, *image1_));
//         should(&(*image1Copy) != &(*image1_));
//         
//         image1Copy->init(data[9]);
//         should(image1Copy->end() == find_if(image1Copy->begin(), image1Copy->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
//         should(image1->end() == find_if(image1->begin(), image1->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
//         
//         std::auto_ptr<Image> image0(image0_->shallowCopy());
//         should(equal(*image0, *image0_));
//         should(&(*image0) != &(*image0_));
//         
// #if !defined(VARIABLE_BANDS_IMAGE) && !defined(SINGLE_BAND_IMAGE)
//         (*image0_) = (*image1_);
//         should(equal(*image1_, *image0_));
//         should(!equal(*image1_, *image0));
//         should(!equal(*image0_, *image0));
// #endif //!VARIABLE_BANDS_IMAGE && !SINGLE_BAND_IMAGE
//     }
//     
//     /** Testet die Methode actualSize(), die die Groesse des kompletten Bildes liefert und nicht der ROI
//     */
//     void testActualSize()
//     {
//         should(image1_->size() == image1_->actualSize());
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         should(image1_->size() != image1_->actualSize());
//         should(image1_->size() == Diff2D(1,3));
//     }
//     
//     /** Testet die Methode actualWidth(), die die Breite des kompletten Bildes liefert und nicht der ROI
//     */
//     void testActualWidth()
//     {
//         should(image1_->width() == image1_->actualWidth());
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         should(image1_->width() == 1);
//         should(image1_->width() != image1_->actualWidth());
//     }
//     
//     /** Testet die Methode actualHeight(), die die Hoehe des kompletten Bildes liefert und nicht der ROI
//     */
//     void testActualHeight()
//     {
//         should(image1_->height() == image1_->actualHeight());
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         should(image1_->height() == 3);
//         should(image1_->height() != image1_->actualHeight());
//     }
// 
//     /** Testet die Methode actualBands(), die die Anzahl der Baender eines Bildes liefert 
//     */
//     void testActualBands()
//     {
//         #ifndef SELECT_BAND_IMAGE
//             should(image1_->bands() == image1_->actualBands());
//         #else
//             should(image1_->bands() != image1_->actualBands());
//             int act_band = image1_->actualBands();
//             Policy::n == 0 ? image1_->setROI(1) : image1_->setROI(Policy::n - 1);
//             should(image1_->actualBands() != 1);
//             should(image1_->actualBands() == act_band);
//         #endif SELECT_BAND_IMAGE    
//     }
//     
//     /** Testet die Methode bands(), die die Anzahl der selektirerten Baender eines Bildes liefert.
//     * Die Anzahl der selektirerten Baender ist bei allen Bildern gleich dem Anzahl der Baender des Bildes,
//     * die Ausnahme bildet SelectBandImage dort ist immer nur ein Band selektiert.
//     */
//     void testBands()
//     {
//         #ifndef SELECT_BAND_IMAGE
//             should(image1_->bands() == image1_->actualBands());
//         #else
//             should(image1_->bands() != image1_->actualBands());
//             should(image1_->bands() == 1);
//             Policy::n == 0 ? image1_->setROI(1) : image1_->setROI(Policy::n - 1);
//             should(image1_->bands() == 1);            
//         #endif SELECT_BAND_IMAGE    
//     }
//     
//     /** testet die Methode roiUpperLeft(), die die Koordinaten der linken
//     * oberen Ecke von ROI liefert.
//     */
//     void testRoiUpperLeft()
//     {
//         should(image1_->roiUpperLeft() == Diff2D(0, 0));
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         should(image1_->roiUpperLeft() == Diff2D(1, 1));
//         should((*image1_)[image1_->roiUpperLeft()] != (*image1_->upperLeft()));
//     }
// 
//     /** testet die Methode setROI(), die eine ROI auf einem Bild setzt.
//     */
//     void testSetROI()
//     {
//         std::auto_ptr<Image> image1_copy(image1_->shallowCopy());               //dient der Vergleichsmoeglichkeit
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         
//         should(image1_->width() == (image1_->actualWidth() - 2));
//         should(image1_->height() == (image1_->actualHeight() - 2));
//         should(image1_->size() == (image1_->actualSize() - Diff2D(2, 2)));
//         should((*image1_copy)(1,1) == (*image1_)(0,0));
//         
//         /*  Differenz zweier ScanOrderIteratoren ist eine int_Zahl
//         */
//         should((image1_->end() - image1_->begin()) == 3);                       // Anzahl der Pixel in der ROI
//         should((image1_copy->end() - image1_copy->begin()) == 15);              // Anzahl der Pixel im Hauptbild
//                 
//         /*  Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
//         *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
//         *   dabei ist an der shallowCopy genau so die ROI gesetzt.
//         */
//         std::auto_ptr<Image> image1_ROI_copy(image1_->shallowCopy());
//          
//         should(!(image1_ROI_copy->size() == image1_ROI_copy->actualSize()));
//         should(image1_ROI_copy->size() == image1_->size());
//         
//         /*  Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
//         *   mit der Groesse und Pixelinitialisierung der ROI
//         */
//         std::auto_ptr<typename Image::CloneType> image1_ROI_clone(image1_->clone());
//         
//         should(image1_ROI_clone->size() == image1_ROI_clone->actualSize());
//         should(image1_ROI_clone->size() == image1_->size());
//                    
//         for(int x = 0; x < image1_ROI_clone->width(); x++)
//             for(int y = 0; y < image1_ROI_clone->height(); y++)
//                 should((*image1_ROI_clone)(x,y) == (*image1_copy)(x + 1,y + 1));
//                 
//           
//         /*  Wenn man die Aenderungen an der ROI vornimmt, aendert sich auch die shallowCopy des Hauptbildes
//         *   aber verschoben um roiUpperLeft
//         */
//         image1_->init(data[7]);
//         
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
//                     
//         for(int x = 0; x < image1_->width(); x++)
//             for(int y = 0; y < image1_->height(); y++)
//                 should((*image1_)(x,y) == static_cast<typename Policy::PixelType>(data[7]));
// 
//                  
//         for(int x = 1; x < (image1_->width() + 1); x++)                                 //um 1 verschoben
//             for(int y = 1; y < (image1_->height() + 1); y++)                            //um 1 verschoben
//                 should((*image1_copy)(x,y) == static_cast<typename Policy::PixelType>(data[7]));
//                 
//         
//         for(int x = 0; x < image1_->width(); x++)                                       //Verhalten einer ROI
//             for(int y = 0; y < image1_->height(); y++)
//                 should((*image1_ROI_copy)(x,y) == static_cast<typename Policy::PixelType>(data[7]));
//                 
//         
//         /* Es wird ein neues Bild erzeugt mit unterschiedlichen Pixeldaten
//         */
//         std::auto_ptr<Image> image1_new(Policy::factory(3, 5));
//         typename Image::Accessor acc;
//         acc = image1_new->accessor();
//         scanOrderIter1_ = image1_new->begin();
//         for(int i = 0; i < data.size(); i++)
//             {   
//                 acc.set(data[i], scanOrderIter1_);
//                 ++scanOrderIter1_; 
//             }
//         
//         std::auto_ptr<Image> image1_new_copy(image1_new->shallowCopy());
//         
//         /*  Setzt man eine ROI auf der shallowCopy so ist sie auch nur in der Copy gesetzt.
//         *   Die Veraenderungen in der ROI sind aber in der Copie und im image1_new sichtbar,
//         *   im image1_new nur verschoben
//         */
//         image1_new_copy->setROI(Diff2D(1,1), Diff2D(1,3));
//         
//         should(image1_new->size() == image1_new->actualSize());
//         should(image1_new_copy->size() != image1_new_copy->actualSize());
//         
//         (*image1_new_copy)(0,0) = data[14];
//         
//         should((*image1_new_copy)(0,0) == static_cast<typename Policy::PixelType>(data[14]));
//         should((*image1_new)(1,1) == (*image1_new_copy)(0,0));
//         
//         /*  Von der ROI kann man in den negativen Bereich zurueckgreifen, da dort das Bild immernoch definiert ist
//         *   Vorsicht! der negative Zugrif darf nicht ausserhalb des eigentlichen Bildes erfolgen
//         */
//         should((image1_new_copy->upperLeft())[Diff2D(-1,-1)] == (*image1_new)(0,0));
//         
//         /*  Eine ROI ist auch auf die Groesse des Bildes ruecksetzbar
//         *   Vorsicht! Die Grenzueberschreitungen beachten!
//         */
//         image1_new_copy->setROI(Diff2D(-1,-1), image1_new_copy->actualSize());                      //Ruecksetzen
//         
//         should(image1_new_copy->size() == image1_new_copy->actualSize());
//         should((image1_new->upperLeft())[Diff2D(0,0)] == (*image1_new_copy)[image1_new_copy->roiUpperLeft()]);
//         
//         
//         /*  Setzt man die ROI mit neagtiven Werten in size so wird der negative Wert vom lowerRight des Bildes abgezogen
//         *   in unserem Falle sind Folgende ROIsetzungen gleich
//         */
//         image1_new_copy->setROI(Diff2D(1,1), Diff2D(-1,-1));
//         should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
//         should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
//         should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
//         
//         image1_new_copy->setROI(Diff2D(-1,-1), image1_new_copy->actualSize());                      //Ruecksetzen
//         image1_new_copy->setROI(Diff2D(1,1), Diff2D(1,-1));                                         //neu setzen
//         
//         should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
//         should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
//         should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
//         
//         image1_new_copy->resetROI();                                                                //Ruecksetzen
//         image1_new_copy->setROI(Diff2D(1,1), Diff2D(-1,3));                                         //neu setzen
//         
//         should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
//         should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
//         should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
//         
//      }//end of testSetRoi()
//      
//      /** testet die Methode resetROI(), die die gesetzte ROI wieder aufhebt
//      */
//      void testResetROI()
//      {
//         std::auto_ptr<typename Image::CloneType> image1(image1_->clone()); 
//         image1_->setROI(Diff2D(1,1), Diff2D(-1,3));
//         image1_->resetROI();  
//         should(equal(*image1, *image1_));
//      }
//      
//      /** Testet die setROI(int band) Methode der Klasse SelectBandImage. Die Methode aendert den 
//      * selektierten Band. 
//      */
//      void testSetROIBand()
//      {
//         should(image1_->getSelectedBand() == Policy::n);
//         image1_->setROI(Policy::n == 0 ? 1 : (Policy::n - 1));
//         should(image1_->getSelectedBand() != Policy::n);
//         should(image1_->getSelectedBand() == Policy::n == 0 ? 1 : (Policy::n - 1));
//         testSetROI();
//      }
//      
//      /** Testet die getSelectedBand() Methode der Klasse SelectBandImage. Die Methode liefert den 
//      * selektierten Band. 
//      */
//      void testGetSelectedBand()
//      {
//         should(image1_->getSelectedBand() == Policy::n);
//         image1_->setROI(Policy::n == 0 ? 1 : (Policy::n - 1));
//         should(image1_->getSelectedBand() != Policy::n);
//         should(image1_->getSelectedBand() == Policy::n == 0 ? 1 : (Policy::n - 1));
//      }
//          
//      void testIsInsideROI()
//      {
//         should(image1_->isInsideROI((image1_->actualSize())- Diff2D(1,1)));
//         
//         image1_->setROI(Diff2D(1,1), Diff2D(1,3));
//         
//         bool ergebnis = true;
//         for ( int i = 1 ; i < 2 ; i++)
//             for (int j = 1; j < 4; j ++)
//                 ergebnis &= image1_->isInsideROI(vigra::Diff2D(i,j));
// 
//         ergebnis &= !image1_->isInsideROI(vigra::Diff2D(0,0)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(2,0)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(0,4)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(2,4)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(2,2)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(0,2)) &&
//                      !image1_->isInsideROI(vigra::Diff2D(2,0)); 
//                      
//         should(ergebnis);
//      }
//      
    /** Testet den Zuweisungsoperator bei dem auf den beiden Seiten das zu testende Bild steht
    */
    void testOperatorAssignmentImage()
    {
        std::auto_ptr<Image> image(Policy::factory());
        (*image) = (*image1_);
        should(equal(*image, *image1_));
        should(&image != &image1_);
        
        (*image) = (*image0_);
        should(!equal(*image, *image1_));
        should(equal(*image, *image0_));
        should( &(*image) !=  &(*image0_));
    }
//     
//     /** Testet den Zuweisungsoperator bei dem auf der linken Seiten ein Pixel steht
//     */
//     void testOperatorAssignmentPixel()
//     {
//         (*image1_) = data[4];
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[4])));
//         
//         (*image0_) = data[5];
//         should(image0_->end() == find_if(image0_->begin(), image0_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
//     }
//     
    /** Testet die init() Methode
    */
    void testInit()
    {
        image1_->init(data[6]);
        should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
        image1_->init(data[7]);
        should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
        should(image1_->end() != find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
        
        image0_->init(data[8]);
        should(image0_->end() == find_if(image0_->begin(), image0_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
    }
    
//     void testResizeInt()
//     {
//         image1_->resize(6, 10);
//         should(image1_->height() == 10);
//         should(image1_->width() == 6);
//         should(image1_->size() == vigra::Diff2D(6,10));
//         
//         image0_->resize(2, 3);
//         should(image0_->height() == 3);
//         should(image0_->width() == 2);
//         should(image0_->size() == vigra::Diff2D(2,3));          
//     }
//     
//     void testResize2D()
//     {
//         image1_->resize(vigra::Diff2D(7,8));
//         should(image1_->height() == 8);
//         should(image1_->width() == 7);
//         should(image1_->size() == vigra::Diff2D(7,8));
//         
//         image0_->resize(vigra::Diff2D(1, 4));
//         should(image0_->height() == 4);
//         should(image0_->width() == 1);
//         should(image0_->size() == vigra::Diff2D(1,4));         
//     }
//     
//     void testResizeIntInit()
//     {
//         image1_->resize(4, 6, data[9]);
//         should(image1_->height() == 6);
//         should(image1_->width() == 4);
//         should(image1_->size() == vigra::Diff2D(4,6));
//         should(image1_->end() == find_if(image1_->begin(), image1_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
//         
//         image0_->resize(1, 4, data[10]);
//         should(image0_->height() == 4);
//         should(image0_->width() == 1);
//         should(image0_->size() == vigra::Diff2D(1,4));
//         should(image0_->end() == find_if(image0_->begin(), image0_->end(), bind2nd(Pixels_not_equal_to<value_type>(), data[10])));      
//     }
//     
//     /** testet die Methode resizeCopy(BasicImage img) an Instanzen der Klasse BasicImage.
//     * Imagehierarchie enthaelt die Methode nicht !!! Also, wird es auch nicht getestet
//     * fuer Imagehierarchie.
//     */
//     void testResizeCopy()
//     {
//         image0_->resizeCopy(*image1_);
//         should(equal(*image0_, *image1_));
//         should(&(*image0_) != &(*image1_));
//         
//         std::auto_ptr<Image> image1(new Image(*image1_));
//         image1->resize(6, 10, data[11]);
//         
//         image1_->resizeCopy(*image1);
//         should(equal(*image1, *image1_));
//         should(!equal(*image0_, *image1_));
//         should(& (*image0_ )!= & (*image1_));
//         should(& (*image1) != & (*image1_));
//     }

    void testWidth()
    {
        should(3 == image1_->width());
        should(0 == image0_->width());
        
        #ifndef IMAGEHIERARCHY                                              // da resize() in Imagehierarchy nicht implementiert ist, so kann auch die Veraenderung nicht getestet werden
            image1_->resize(4,6);
            should(4 == image1_->width());
            image0_->resize(5,7);
            should(5 == image0_->width());
        #endif IMAGEHIERARCHY
    } 
    
    void testHeight()
    {
        should(5 == image1_->height());
        should(0 == image0_->height());
        
        #ifndef IMAGEHIERARCHY                                              // da resize() in Imagehierarchy nicht implementiert ist, so kann auch die Veraenderung nicht getestet werden
            image1_->resize(4,6);
            should(6 == image1_->height());
            image0_->resize(7,8);
            should(8 == image0_->height());
        #endif IMAGEHIERARCHY
    } 
    
    void testSize()
    {
        should(image1_->size() == vigra::Diff2D(3,5));
        should(image0_->size() == vigra::Diff2D(0,0));
        
        #ifndef IMAGEHIERARCHY                                              // da resize() in Imagehierarchy nicht implementiert ist, so kann auch die Veraenderung nicht getestet werden
            image1_->resize(4,6);
            should(vigra::Diff2D(4,6) == image1_->size());
            image1_->resize(7,8);
            should(vigra::Diff2D(7,8) == image1_->size());
        #endif IMAGEHIERARCHY
    }
    
    void testIsInside()
    {
        bool ergebnis = true;
        for ( int i = 0 ; i < 3 ; i++)
            for (int j = 0; j < 5; j ++)
                ergebnis &= image1_->isInside(vigra::Diff2D(i,j));

        ergebnis &= !image1_->isInside(vigra::Diff2D(-1,-1)) &&
                     !image1_->isInside(vigra::Diff2D(-1,2)) &&
                     !image1_->isInside(vigra::Diff2D(2,-1)) &&
                     !image1_->isInside(vigra::Diff2D(-1,5)) &&
                     !image1_->isInside(vigra::Diff2D(3,-1)) &&
                     !image1_->isInside(vigra::Diff2D(3,2)) &&
                     !image1_->isInside(vigra::Diff2D(3,5)) &&
                     !image1_->isInside(vigra::Diff2D(2,5));
                     
        should(ergebnis);
        
        ergebnis = !image0_->isInside(vigra::Diff2D(-1,-1)) &&
                   !image0_->isInside(vigra::Diff2D(-1,0)) &&
                   !image0_->isInside(vigra::Diff2D(0,-1)) &&
                   !image0_->isInside(vigra::Diff2D(1,0)) &&
                   !image0_->isInside(vigra::Diff2D(0,1)) &&
                   !image0_->isInside(vigra::Diff2D(1,1)) &&
                   !image0_->isInside(vigra::Diff2D(1,-1)) &&
                   !image0_->isInside(vigra::Diff2D(-1,1)) &&
                   !image0_->isInside(vigra::Diff2D(0,0));
        
        should(ergebnis);
    }
    
    /** testet den operator[](Diff2D)
    */
    void testOperatorIndex2D()
    {
        typename Image::ScanOrderIterator k = image1_->begin();
        for ( int y = 0 ; y < 5 ; y++)
            for (int x = 0; x < 3; x++, k++)
                should((*image1_)[Diff2D(x,y)] == *k);
    }
    
    /** test der Funktion "PixelType const& operator[] (Diff2D const &d) const"
    */
    void testOperatorIndex2DConst()
    {
        std::auto_ptr<Image> const & cimage = image1_;
        typename Image::ScanOrderIterator k = cimage->begin();
        for ( int y = 0 ; y < 5 ; y++)
            for (int x = 0; x < 3; x++, k++)
                should((*cimage)[Diff2D(x,y)] == *k);
    }
    
    /** testet den operator(int x, int y).
    */
    void testOperatorFunctionCallInt()
    {
        for ( int y = 0 ; y < 5 ; y++)
            for (int x = 0; x < 3; x++)
                should((*image1_)(x,y) == (*image1_)[Diff2D(x,y)]);
    }
    
    /** testet den const operator()(int x, int y) const.
    */   
    void testOperatorFunctionCallIntConst()
    {
        std::auto_ptr<Image> const & cimage = image1_;
        for ( int y = 0 ; y < 5 ; y++)
            for (int x = 0; x < 3; x++)
                should((*cimage)(x,y) == (*cimage)[Diff2D(x,y)]);
    }
    
    /** testet den operator[](int dy)
    */
    void testOperatorDoubleIndex()
    {
        for ( int x = 0 ; x < 3 ; x++)
            for (int y = 0; y < 5; y++)
                should((*image1_)[y][x] == (*image1_)(x,y));
    }
    
    /** testet den const operator[](int dy) 
    */
    void testOperatorDoubleIndexConst()
    {
        std::auto_ptr<Image> const & cimage = image1_;
        for ( int x = 0 ; x < 3 ; x++)
            for (int y = 0; y < 5; y++)
                should((*cimage)[y][x] == (*image1_)(x,y));
    }
    
    /** testet die upperLeft() - Funktion der zu testenden Imageklasse
    */
    void testUpperLeft()
    {
        /*
        * upperLeft liefert einen Iterator zurueck, um den Iterator 
        * auf Richtigkeit zu ueberpruefen entreferenziere ich ihn, er soll
        * dann eine !!!!SomePixelType!!!! zuruekgeben und das sollte derselbe Typ und value
        * sein wie bei Index2d(0,0) vom anderen Image. 
        * Eigentlich werden PixelTypes vergliechen, Da der Iterator keine
        * Vergleichsmoeglichkeiten bietet
        */
        should((image1_->upperLeft())[2][1] == (*image1_)(1, 2));
        should(*(image1_->upperLeft() + image1_->size() - Diff2D(1, 1)) == (*image1_)(2, 4));
        should(image1_->upperLeft()[Diff2D(2,3)] == (*image1_)[Diff2D(2,3)]);

        (*image1_->upperLeft()) = data[3];
        should((*image1_->upperLeft()) == static_cast<typename Policy::PixelType>(data[3]));
    }

    void testLowerRight()
    {
        should(*(image1_->lowerRight() - Diff2D(1,1)) == (*image1_)(image1_->width() - 1 , image1_->height() - 1));   
        should(*(image1_->lowerRight() - Diff2D(3, 5)) == *(image1_->upperLeft()));
    }
    
    void testConstUpperLeft()
    {
        Image const * cimage = image1_.get();
        should(*(cimage->upperLeft()) == (*cimage)[Diff2D(0,0)]);
        should(cimage->upperLeft()[2][1] == (*cimage)(1, 2));
        should(*(cimage->upperLeft() + Diff2D(2, 4)) == (*cimage)(2, 4));
        should(cimage->upperLeft()[Diff2D(2,3)] == (*cimage)[Diff2D(2,3)]);       
    }
    
    void testConstLowerRight()
    {
        Image const * cimage = image1_.get();
        should(*(cimage->lowerRight() - Diff2D(1,1)) == (*cimage)(image1_->width() - 1 , image1_->height() - 1));         
        should(*(cimage->lowerRight() - Diff2D(3, 5)) == *(cimage->upperLeft()));
    }    
        
    void  testBegin()
    { 
        should(*image1_->begin() == *image1_->upperLeft());

        //An dieser Stelle sollte begin()++ angesetzt werden funktioniert leider nicht (nur in BasicImage)                 //TO DO
        should(*(image1_->begin() + 1) == (*image1_)[Diff2D(1,0)]);
        should(*(image1_->begin() + image1_->width()) == (*image1_)[Diff2D(0,1)]);
        should(*(image1_->begin() + ((image1_->width()*(image1_->height()) - 1))) == (*image1_)[Diff2D(2,4)]);
    }

    void testEnd()
    {       
        should(*(image1_->end() - 1) == (*image1_)(image1_->width() - 1 , image1_->height() - 1));
        should(*(image1_->end() -(image1_->width()*image1_->height())) == (*image1_)(0, 0));    
    }

    void  testConstBegin()
    {           
        Image const * cimage = image1_.get();
        should(*cimage->begin() == *cimage->upperLeft());
        should(*(cimage->begin() + 1) == (*cimage)[Diff2D(1,0)]);
        should(*(cimage->begin() + image1_->width()) == (*cimage)[Diff2D(3,0)]);
        should(*(cimage->begin() + ((image1_->width()*image1_->height()) - 1)) == (*cimage)[Diff2D(2,4)]);
    }

    void testConstEnd()
    {           
        Image const * cimage = image1_.get();
        should(*(cimage->end() - 1) == (*cimage)(image1_->width() - 1 , image1_->height() - 1));   
        should(*(cimage->end() -(image1_->width()*image1_->height())) == (*cimage)(0, 0));    
    }
    
    /** Testet den Accessor
    */   
    void testAccessor()
    {
        typename Image::Accessor acc = image1_->accessor();
        
        traverserIter1_ = image1_->upperLeft();
        should(acc(traverserIter1_) == acc1_(traverserIter1_));
        acc.set(data[6],traverserIter1_);
        should(acc(traverserIter1_) == static_cast<typename Policy::PixelType>(data[6]));
    } 
    
    /** Testet den konstanten Accessor
    */
    void testAccessorConst()
    {
        image1_->init(data[5]);
        Image const * cimage = image1_.get();                                           //get ist eine Funktion des std::auto_ptr<...>
        typename Image::ConstAccessor acc = cimage->accessor();
        typename Image::ConstIterator i = cimage->upperLeft();
        should(acc(i) == static_cast<typename Policy::PixelType>(data[5]));
    }
    
    /** Testet richtige Funktionsweise aller Iteratoren und Zugriffsmoeglichkeiten (Accessor, Diff2D usw.) auf 
    * ein Pixel des Bildes
    */
    void testAllAccessAndSetMethods()
    {
        scanOrderIter1_ = image1_->begin();
        traverserIter1_ = image1_->upperLeft();
        rowIter1_ = traverserIter1_.rowIterator();
        
        for(int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 3; j++)
            {   
                should ((*image1_)(j,i) == (*image1_)[vigra::Diff2D(j,i)]);
                should ((*image1_)(j,i) == (*image1_)[i][j]);
                should ((*image1_)(j,i) == *scanOrderIter1_ );                            
                should ((*image1_)(j,i) == *traverserIter1_ ); 
                should ((*image1_)(j,i) == *rowIter1_ ); 
                
                (*image1_)(j,i) = data[7];
                should ((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[7]));
             
                should ((*image1_)(j,i) == (*image1_)[vigra::Diff2D(j,i)]);
                should ((*image1_)(j,i) == (*image1_)[i][j]);
                should ((*image1_)(j,i) == acc1_(scanOrderIter1_ ));
                should ((*image1_)(j,i) == acc1_(traverserIter1_ ));
                should ((*image1_)(j,i) == acc1_(rowIter1_));
                
                (*image1_)[Diff2D(j,i)] = data[8];
                should ((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[8]));
                  
                should ((*image1_)(j,i) == acc1_(scanOrderIter1_));
                should ((*image1_)(j,i) == acc1_(rowIter1_));
                should ((*image1_)(j,i) == acc1_(traverserIter1_));
                
                (*image1_)[i][j] = data[9];
                should ((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[9]));
                
                should ((*image1_)(j,i) == (*image1_)[vigra::Diff2D(j,i)]);
                should ((*image1_)(j,i) == (*image1_)[i][j]);
                should((*image1_)(j,i) == acc1_(scanOrderIter1_));                        
                should((*image1_)(j,i) == acc1_(traverserIter1_));
                should((*image1_)(j,i) == acc1_(rowIter1_));
                
                *traverserIter1_ = data[10];
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[10]));
                
                traverserIter1_[Diff2D(0,0)] = data[11];
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[11]));
                
                acc1_.set(data[12], traverserIter1_);
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[12]));
                
                acc1_.set(data[13], traverserIter1_, Diff2D(0,0));
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[13]));
                 
                *rowIter1_ = data[10];
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[10]));
                
                rowIter1_[0] = data[0];
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[0]));
                
                acc1_.set(data[1], rowIter1_);
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[1]));
                
                acc1_.set(data[2], rowIter1_, 0);
                should((*image1_)(j,i) == static_cast<typename Policy::PixelType>(data[2]));

                scanOrderIter1_++;
                traverserIter1_.x++;
                rowIter1_++;
            }
            traverserIter1_.x = image1_->upperLeft().x;
            traverserIter1_.y++;
            rowIter1_ = traverserIter1_.rowIterator();
        }
        
        /** tests columnIterator */
        traverserIter1_ = image1_->upperLeft();
        columnIter1_ = traverserIter1_.columnIterator();
        for(int i = 0; i < 3; i++)  
        {
            for (int j = 0; j < 5; j++)
            {
                should((*image1_)(i,j) == acc1_(columnIter1_));
                acc1_.set(data[3], columnIter1_);
                should((*image1_)(i,j) == acc1_(columnIter1_));
                columnIter1_++;
                traverserIter1_.y++;
            
            }
            traverserIter1_.y = image1_->upperLeft().y;
            traverserIter1_.x++;
            columnIter1_ = traverserIter1_.columnIterator();
        }         
    } 
       
};

#if defined(FIXED_RGB_IMAGE) ||defined(RGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE)
    /** Klasse testet die zusaetzlichen Methoden der RGBBilder, die restlichen Methoden, die 
    * auch in anderen Bildern enthalten sind, werden mittels der Oberklasse ImageTest getestet.
    */
    template <class Policy>
    class RGBImageTest : public ImageTest <Policy>
    {
    public:
        typedef typename Policy::Image Image;
        typedef typename ImageTest<Policy>::ChildImage ChildImage;
        typedef typename Policy::value_type value_type;
        typename Policy::data_array_type data;

        RGBImageTest()
        : data(Policy::getData())
        {
        }

        void testRed()
        {
            std::auto_ptr<Image> image(Policy::factory(2,3, data[0]));
            should((*image)(1,1)[0] == data[0].red());
            should((*image)[Diff2D(1,2)][0] == data[0].red());
        }    
        void testGreen()
        {
            std::auto_ptr<Image> image(Policy::factory(2,3, data[3]));
            should((*image)(1,1)[1] == data[3].green());
            should((*image)[Diff2D(1,2)][1] == data[3].green());
        }    
        void testBlue()
        {
            std::auto_ptr<Image> image(Policy::factory(2,3, data[5]));
            should((*image)(1,1)[2] == data[5].blue());
            should((*image)[Diff2D(1,2)][2] == data[5].blue());
        } 
    };
#endif // FIXED_RGB_IMAGE ||RGB_IMAGE  || defined (VARIABLE_BANDS_IMAGE)

template <class POLICY>
struct TImageTestSuite
: public vigra::test_suite
{
    TImageTestSuite(const char * to_test_Image_name)
    :vigra::test_suite(to_test_Image_name)
    {
        add( testCase( &ImageTest<POLICY>::testImageDefaultConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntConstuctor));
        add( testCase( &ImageTest<POLICY>::testImage2DConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntPixelConstuctor));
        add( testCase( &ImageTest<POLICY>::testCopyConstructor));

#ifdef IMAGEHIERARCHY
        add( testCase( &ImageTest<POLICY>::testInnerImageConstructor));
        add( testCase( &ImageTest<POLICY>::testClone));
        add( testCase( &ImageTest<POLICY>::testShallowCopy));
        add( testCase( &ImageTest<POLICY>::testSetROI));
        add( testCase( &ImageTest<POLICY>::testResetROI));
        add( testCase( &ImageTest<POLICY>::testActualSize));
        add( testCase( &ImageTest<POLICY>::testActualWidth));
        add( testCase( &ImageTest<POLICY>::testActualHeight));
        add( testCase( &ImageTest<POLICY>::testActualBands));
        add( testCase( &ImageTest<POLICY>::testBands));
        add( testCase( &ImageTest<POLICY>::testRoiUpperLeft));
        add( testCase( &ImageTest<POLICY>::testIsInsideROI));
#else      
        add( testCase( &ImageTest<POLICY>::testOperatorAssignmentPixel)); 
        add( testCase( &ImageTest<POLICY>::testResizeInt));
        add( testCase( &ImageTest<POLICY>::testResize2D));
        add( testCase( &ImageTest<POLICY>::testResizeIntInit));
        add( testCase( &ImageTest<POLICY>::testResizeCopy));
#endif  IMAGEHIERARCHY 
       
#if !defined(VARIABLE_BANDS_IMAGE) && !defined(SINGLE_BAND_IMAGE)                       // da VariableBandsImage eine abstrakte Klasse ist, kann man auch nicht den Zuweisungsoperator testen, da der default Konstruktor aufgerufen wird.
        add( testCase( &ImageTest<POLICY>::testOperatorAssignmentImage));
#endif //!VARIABLE_BANDS_IMAGE && !SINGLE_BAND_IMAGE

#ifdef SELECT_BAND_IMAGE                                                                // die Klasse SelectBandImage enthaelt zusaetzliche Methoden setROI(int band) und getSelectedBand(). Die Methode setROI(int band) aendert den selektierten Band und getSelectedBand() liefert den selektierten Band. Und hier werden deren Tests aufgerufen
        add( testCase( &ImageTest<POLICY>::testSetROIBand));
        add( testCase( &ImageTest<POLICY>::testGetSelectedBand));
#endif SELECT_BAND_IMAGE

        add( testCase( &ImageTest<POLICY>::testInit));
        add( testCase( &ImageTest<POLICY>::testWidth));
        add( testCase( &ImageTest<POLICY>::testHeight));
        add( testCase( &ImageTest<POLICY>::testSize)); 
        add( testCase( &ImageTest<POLICY>::testIsInside));
        add( testCase( &ImageTest<POLICY>::testOperatorIndex2D));
        add( testCase( &ImageTest<POLICY>::testOperatorIndex2DConst));
        add( testCase( &ImageTest<POLICY>::testOperatorFunctionCallInt));
        add( testCase( &ImageTest<POLICY>::testOperatorFunctionCallIntConst));
        add( testCase( &ImageTest<POLICY>::testUpperLeft));
        add( testCase( &ImageTest<POLICY>::testLowerRight));
        add( testCase( &ImageTest<POLICY>::testConstUpperLeft));
        add( testCase( &ImageTest<POLICY>::testConstLowerRight));
        add( testCase( &ImageTest<POLICY>::testBegin));
        add( testCase( &ImageTest<POLICY>::testEnd));
        add( testCase( &ImageTest<POLICY>::testConstBegin));
        add( testCase( &ImageTest<POLICY>::testConstEnd));
        add( testCase( &ImageTest<POLICY>::testOperatorDoubleIndex));
        add( testCase( &ImageTest<POLICY>::testOperatorDoubleIndexConst));
        add( testCase( &ImageTest<POLICY>::testAccessor));
        add( testCase( &ImageTest<POLICY>::testAccessorConst));
        add( testCase( &ImageTest<POLICY>::testAllAccessAndSetMethods));
    }
};//end of struct TImageTestSuite

#if defined(FIXED_RGB_IMAGE) ||defined(RGB_IMAGE)  || defined (VARIABLE_BANDS_IMAGE)
    /** testet die drei extra Methoden von RGBImages
    */
    template <class POLICY>
    struct RGBImageTestSuite
    :public TImageTestSuite<POLICY>
    {
        RGBImageTestSuite(const char * to_test_Image_name)
        : TImageTestSuite<POLICY>(to_test_Image_name)
        {
            add( testCase( &RGBImageTest<POLICY>::testRed));
            add( testCase( &RGBImageTest<POLICY>::testGreen));
            add( testCase( &RGBImageTest<POLICY>::testBlue));
        }
    };// end of struct RGBImageTestSuite
#endif // FIXED_RGB_IMAGE ||RGB_IMAGE  || defined (VARIABLE_BANDS_IMAGE)
