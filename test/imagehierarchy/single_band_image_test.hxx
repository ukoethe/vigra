#ifndef SINGLE_BAND_IMAGE_TEST_HXX
#define SINGLE_BAND_IMAGE_TEST_HXX

#include "imagehierarchy_test.hxx"

template<class Policy>
class SingleBandImageTest
: public ImageHierarchyTest<Policy>
{    
public:
    typedef typename Policy::Image              Image;              // zu testende Klasse z.B. GrayImage, SingleBandImage usw.
    typedef typename Policy::ChildImage         ChildImage;         // unterscheidet sich von zu testender Klasse Image, nur wenn einer der abstrakten Klassen VariableBandImage oder SingleBandImage oder die Klasse SelectBandImage getestet wird, sonst entspricht es der Klasse Image.  
    typedef typename Policy::value_type         value_type;         // value_type der zu testenden Klasse
    typedef typename Policy::child_value_type   child_value_type;   // value_type der Klasse an der die zu testende Klasse getestet wird, also z.B. SingleBandImage wird am Beispiel von Vector2Image getestet dann ist child_value_type die value_type von Vector2Image
    
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
        should(&(*image1) != &(*image1_));                                              
        
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
        should(&(*image0) != &(*image0_));
    }
    
    /** testet die shallowCopy() Methode der Klassen aus imagehierarchy.
    */
    void testShallowCopy()
    {
        /*
        *  Im Falle der Aenderungsvornehmungen an einem der Images, werden die Aenderungen an beiden sichtbar
        *  Es gibt nur zwei Aenderungsmoeglichkeiten init() und Zuweisung eines anderen Images
        *  bei der Zuweisung werden die Zeiger von den eigenen Daten auf die anderen umgeleitet,
        *  und das entspricht nicht dem Sinn der shallowCopy
        */
        std::auto_ptr<Image> image1(image1_->shallowCopy());
        should(equal(*image1, *image1_));     
        should(&(*image1) != &(*image1_));
        
        /* Aenderung mit der init-Funktion
        */
        image1->init(data[7]);
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
        
        image1->init(data[8]);
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
        
        /* Eine shallowCopy zeigt auf die selben Daten des kopierten Objektes
        */
        std::auto_ptr<Image> image1Copy(image1->shallowCopy());
        should(equal(*image1Copy, *image1_));
        should(&(*image1Copy) != &(*image1_));
        
        image1Copy->init(data[9]);
        should(image1Copy->end() == std::find_if(image1Copy->begin(), image1Copy->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
        
        std::auto_ptr<Image> image0(image0_->shallowCopy());
        should(equal(*image0, *image0_));
        should(&(*image0) != &(*image0_));
    }
};

template <class POLICY>
struct SingleBandImageTestSuite
: vigra::test_suite
{
    SingleBandImageTestSuite(const char * to_test_Image_name)
    : vigra::test_suite(to_test_Image_name)
    {
        add( testCase( &ImageTest<POLICY>::testImageDefaultConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntConstuctor));
        add( testCase( &ImageTest<POLICY>::testImage2DConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntPixelConstuctor));
        add( testCase( &ImageTest<POLICY>::testCopyConstructor));
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
        add( testCase( &ImageHierarchyTest<POLICY>::testInnerImageConstructor));
        add( testCase( &SingleBandImageTest<POLICY>::testClone));
        add( testCase( &SingleBandImageTest<POLICY>::testShallowCopy));
        add( testCase( &ImageHierarchyTest<POLICY>::testSetROI));
        add( testCase( &ImageHierarchyTest<POLICY>::testResetROI));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualSize));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualWidth));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualHeight));
        add( testCase( &ImageHierarchyTest<POLICY>::testActualBands));
        add( testCase( &ImageHierarchyTest<POLICY>::testBands));
        add( testCase( &ImageHierarchyTest<POLICY>::testRoiUpperLeft));
        add( testCase( &ImageHierarchyTest<POLICY>::testIsInsideROI));
    }
};//end of struct SingleBandImageTestSuite
#endif // SINGLE_BAND_IMAGE_TEST_HXX
