#ifndef SELECT_IMAGE_TEST_HXX
#define SELECT_IMAGE_TEST_HXX

#include "imagehierarchy_test.hxx"

template<class Policy>
class SelectBandImageTest
: public ImageHierarchyTest<Policy>
{
public:
    typedef typename Policy::Image              Image;              // zu testende Klasse z.B. GrayImage, VariableBandsImage usw.
    typedef typename Policy::ChildImage         ChildImage;         // unterscheidet sich von zu testender Klasse Image, nur wenn einer der abstrakten Klassen VariableBandImage oder SingleBandImage oder die Klasse SelectBandImage getestet wird, sonst entspricht es der Klasse Image.  
    typedef typename Policy::value_type         value_type;         // value_type der zu testenden Klasse
    typedef typename Policy::child_value_type   child_value_type;   // value_type der Klasse an der die zu testende Klasse getestet wird, also z.B. VariableBandsImage wird am Beispiel von Vector2Image getestet dann ist child_value_type die value_type von Vector2Image
    
    /**  Testet den "(WIDTH, HEIGHT, PIXELTYPE)"-Konstruktor der zu testenden Imageklasse
    */
    void testImageIntPixelConstuctor()
    {
        std::auto_ptr<Image> image1(Policy::factory(2,3, child_data[0]));
        should(image1->height() == 3);
        should(image1->width() == 2);
        should(image1->size() == vigra::Diff2D(2,3));
    
        // Bei SelectBandImage wird nur der selektierte Band mit dem Pixel child_data[0] initialisiert   
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[0][Policy::n])));
        
        std::auto_ptr<Image> image2(Policy::factory(0, 0, child_data[1]));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
    
        // Bei SelectBandImage wird nur der selektierte Band mit dem Pixel child_data[1] initialisiert   
        should(image2->end() == std::find_if(image2->begin(), image2->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[1][Policy::n])));
    }
    
    /** Testet den Copy Konstruktor ( Image(Image img) ).
    */
    void testCopyConstructor()
    {
        std::auto_ptr<Image> image1_copy(Policy::factory(*image1_));
        std::auto_ptr<Image> image1 = image1_;

        std::auto_ptr<Image> image0_copy(Policy::factory(*image0_));

        should(image1_copy->height() == 5);
        should(image1_copy->width() == 3);
        should(image1_copy->size() == vigra::Diff2D(3,5));
        should(equalPixels(image1->begin(),image1->end(), image1_copy->begin()));
        should(image0_copy->height() == 0);
        should(image0_copy->width() == 0);
        should(image0_copy->size() == vigra::Diff2D(0,0));
        should(equalPixels(image0_->begin(), image0_->end(), image0_copy->begin()));
     }
    
    /** testet die Konstruktorden der imagehierarchy Klasse, an die die InnerImage (BasicImage<float>)
    * uebergeben wird, also Image (InnerImage i)
    */
    void testInnerImageConstructor()
    {
        std::auto_ptr<Image> image1(Policy::factory(new typename ChildImage::InnerImage(3, 4, child_data[0]))); 
        should(image1->height() == 4);
        should(image1->width() == 3);
        should(image1->size() == vigra::Diff2D(3,4));
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[0][Policy::n])));

        std::auto_ptr<Image> image2(Policy::factory(new typename ChildImage::InnerImage(0, 0, child_data[1])));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
        should(image2->end() == std::find_if(image2->begin(), image2->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[1][Policy::n])));
    }
    
    /** testet die clone() Methode der Klasse aus imagehierarchy
    */
    void testClone()
    {
        /*
        *  Im Falle der Aenderungsvornehmungen an einem der Images, werden die Aenderungen auch nur an einem sichtbar
        */
        std::auto_ptr<typename Image::CloneType> image1(image1_->clone()); 
        should(equal(*image1, *image1_));                                               
        should(equalPixels(image1_->begin(),image1_->end(), image1->begin()));             
        should(static_cast<vigra::SingleBandImage *> (&(*image1)) != static_cast<vigra::SingleBandImage *> (&(*image1_)));
        
        /* Aenderung mit der init-Funktion
        */
        image1->init(data[5]); 
        should((*image1_->begin()) != static_cast<typename Image::PixelType> (data[5]));
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
        
        image1_->init(data[6]);
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
        
        std::auto_ptr<typename Image::CloneType> image0(image0_->clone());
        should(equal(*image0, *image0_));
        should(equalPixels(image0_->begin(),image0_->end(), image0->begin()));image1_->init(data[6]);
        should(static_cast<vigra::SingleBandImage *> (&(*image0)) != static_cast<vigra::SingleBandImage *> (&(*image0_)));

        /* Aenderung mit dem Zuweisungsoperator -> siehe in testShallowCopy()
        */
        (*image0_) = (*image1_);                              
        should(image0->size() == vigra::Diff2D(0,0));
        should(image0_->size() == vigra::Diff2D(3,5));
        should(!equal(*image0, *image0_));
        should(equal(*image1_, *image0_));
    }

    /** Testet die Methode actualBands(), die die Anzahl der Baender eines Bildes liefert 
    */
    void testActualBands()
    {
        should(image1_->bands() != image1_->actualBands());
        int act_band = image1_->actualBands();
        Policy::n == 0 ? image1_->setROI(1) : image1_->setROI(Policy::n - 1);
        should(image1_->actualBands() != 1);
        should(image1_->actualBands() == act_band);
    }
    
    /** Testet die Methode bands(), die die Anzahl der selektirerten Baender eines Bildes liefert.
    * Die Anzahl der selektirerten Baender ist bei allen Bildern gleich dem Anzahl der Baender des Bildes,
    * die Ausnahme bildet SelectBandImage dort ist immer nur ein Band selektiert.
    */
    void testBands()
    {
        should(image1_->bands() != image1_->actualBands());
        should(image1_->bands() == 1);
        Policy::n == 0 ? image1_->setROI(1) : image1_->setROI(Policy::n - 1);
        should(image1_->bands() == 1);
    }
    
    /** Testet die setROI(int band) Methode der Klasse SelectBandImage. Die Methode aendert den 
    * selektierten Band. 
    */
    void testSetROIBand()
    {
    should(image1_->getSelectedBand() == Policy::n);
    image1_->setROI(Policy::n == 0 ? 1 : (Policy::n - 1));
    should(image1_->getSelectedBand() != Policy::n);
    should(image1_->getSelectedBand() == ((Policy::n == 0) ? 1 : (Policy::n - 1)));
    testSetROI();
    }

    /** Testet die getSelectedBand() Methode der Klasse SelectBandImage. Die Methode liefert den 
    * selektierten Band. 
    */
    void testGetSelectedBand()
    {
    should(image1_->getSelectedBand() == Policy::n);
    image1_->setROI(Policy::n == 0 ? 1 : (Policy::n - 1));
    should(image1_->getSelectedBand() != Policy::n);
    should(image1_->getSelectedBand() == ((Policy::n == 0) ? 1 : (Policy::n - 1)));
    }

};

template <class POLICY>
struct SelectBandImageTestSuite
: vigra::test_suite
{ 
    SelectBandImageTestSuite(const char * to_test_Image_name)
    : vigra::test_suite(to_test_Image_name)
    {
        add( testCase( &ImageTest<POLICY>::testImageDefaultConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntConstuctor));
        add( testCase( &ImageTest<POLICY>::testImage2DConstuctor));
        add( testCase( &SelectBandImageTest<POLICY>::testImageIntPixelConstuctor));
        add( testCase( &SelectBandImageTest<POLICY>::testCopyConstructor));
        add( testCase( &ImageTest<POLICY>::testOperatorAssignmentImage));
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
        add( testCase( &SelectBandImageTest<POLICY>::testInnerImageConstructor));
        add( testCase( &SelectBandImageTest<POLICY>::testClone));
        add( testCase( &ImageHierarchyTest<POLICY>::testShallowCopy));
        add( testCase( &ImageHierarchyTest<POLICY>::testSetROI));
        add( testCase( &ImageHierarchyTest<POLICY>::testResetROI));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualSize));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualWidth));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualHeight));
        add( testCase( &SelectBandImageTest<POLICY>::testActualBands));
        add( testCase( &SelectBandImageTest<POLICY>::testBands));
        add( testCase( &ImageHierarchyTest<POLICY>::testRoiUpperLeft));
        add( testCase( &ImageHierarchyTest<POLICY>::testIsInsideROI));
        add( testCase( &SelectBandImageTest<POLICY>::testSetROIBand));
        add( testCase( &SelectBandImageTest<POLICY>::testGetSelectedBand));
    }
};//end of struct SelectBandImageTestSuite
#endif // SELECT_IMAGE_TEST_HXX
