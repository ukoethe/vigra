#ifndef BASIC_IMAGE_TEST_HXX
#define BASIC_IMAGE_TEST_HXX

#include "parent_test_class.hxx"

template<class Policy>
class BasicImageTest
: public ImageTest<Policy>
{
public:
    typedef typename Policy::Image              Image;              // zu testende Klasse z.B. GrayImage, VariableBandsImage usw.
    typedef typename Policy::ChildImage         ChildImage;         // unterscheidet sich von zu testender Klasse Image, nur wenn einer der abstrakten Klassen VariableBandImage oder SingleBandImage oder die Klasse SelectBandImage getestet wird, sonst entspricht es der Klasse Image.  
    typedef typename Policy::value_type         value_type;         // value_type der zu testenden Klasse
    typedef typename Policy::child_value_type   child_value_type;   // value_type der Klasse an der die zu testende Klasse getestet wird, also z.B. VariableBandsImage wird am Beispiel von Vector2Image getestet dann ist child_value_type die value_type von Vector2Image
    
    
    /** Testet den Zuweisungsoperator bei dem auf der linken Seiten ein Pixel steht
    */
    void testOperatorAssignmentPixel()
    {
        (*image1_) = data[4];
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[4])));
        
        (*image0_) = data[5];
        should(image0_->end() == std::find_if(image0_->begin(), image0_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[5])));
    }
    
    void testResizeInt()
    {
        image1_->resize(6, 10);
        should(image1_->height() == 10);
        should(image1_->width() == 6);
        should(image1_->size() == vigra::Diff2D(6,10));
        
        image0_->resize(2, 3);
        should(image0_->height() == 3);
        should(image0_->width() == 2);
        should(image0_->size() == vigra::Diff2D(2,3));          
    }
    
    void testResize2D()
    {
        image1_->resize(vigra::Diff2D(7,8));
        should(image1_->height() == 8);
        should(image1_->width() == 7);
        should(image1_->size() == vigra::Diff2D(7,8));
        
        image0_->resize(vigra::Diff2D(1, 4));
        should(image0_->height() == 4);
        should(image0_->width() == 1);
        should(image0_->size() == vigra::Diff2D(1,4));         
    }
    
    void testResizeIntInit()
    {
        image1_->resize(4, 6, data[9]);
        should(image1_->height() == 6);
        should(image1_->width() == 4);
        should(image1_->size() == vigra::Diff2D(4,6));
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[9])));
        
        image0_->resize(1, 4, data[10]);
        should(image0_->height() == 4);
        should(image0_->width() == 1);
        should(image0_->size() == vigra::Diff2D(1,4));
        should(image0_->end() == std::find_if(image0_->begin(), image0_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[10])));      
    }
    
    /** testet die Methode resizeCopy(BasicImage img) an Instanzen der Klasse BasicImage.
    * Imagehierarchie enthaelt die Methode nicht !!! Also, wird es auch nicht getestet
    * fuer Imagehierarchie.
    */
    void testResizeCopy()
    {
        image0_->resizeCopy(*image1_);
        should(equal(*image0_, *image1_));
        should(&(*image0_) != &(*image1_));
        
        std::auto_ptr<Image> image1(new Image(*image1_));
        image1->resize(6, 10, data[11]);
        
        image1_->resizeCopy(*image1);
        should(equal(*image1, *image1_));
        should(!equal(*image0_, *image1_));
        should(& (*image0_ )!= & (*image1_));
        should(& (*image1) != & (*image1_));
    }
};// end of class BasicImageTest

template <class POLICY>
struct BasicImageTestSuite
: public ImageTestSuite<POLICY>
{
    BasicImageTestSuite(const char * to_test_Image_name)
    : ImageTestSuite<POLICY>(to_test_Image_name)
    {
        add( testCase( &BasicImageTest<POLICY>::testOperatorAssignmentPixel)); 
        add( testCase( &BasicImageTest<POLICY>::testResizeInt));
        add( testCase( &BasicImageTest<POLICY>::testResize2D));
        add( testCase( &BasicImageTest<POLICY>::testResizeIntInit));
        add( testCase( &BasicImageTest<POLICY>::testResizeCopy));
    }
};//end of struct BasicImageTestSuite

#endif // BASIC_IMAGE_TEST_HXX
