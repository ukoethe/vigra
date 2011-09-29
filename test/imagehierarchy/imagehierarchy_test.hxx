#ifndef IMAGEHIERARCHY_TEST_HXX
#define IMAGEHIERARCHY_TEST_HXX

#include "parent_test_class.hxx"

template<class Policy>
class ImageHierarchyTest
: public ImageTest<Policy>
{    
public:
    typedef typename Policy::Image              Image;              // zu testende Klasse z.B. GrayImage, VariableBandsImage usw.
    typedef typename Policy::ChildImage         ChildImage;         // unterscheidet sich von zu testender Klasse Image, nur wenn einer der abstrakten Klassen VariableBandImage oder SingleBandImage oder die Klasse SelectBandImage getestet wird, sonst entspricht es der Klasse Image.  
    typedef typename Policy::value_type         value_type;         // value_type der zu testenden Klasse
    typedef typename Policy::child_value_type   child_value_type;   // value_type der Klasse an der die zu testende Klasse getestet wird, also z.B. VariableBandsImage wird am Beispiel von Vector2Image getestet dann ist child_value_type die value_type von Vector2Image
    
    /** testet die Konstruktorden der imagehierarchy Klasse, an die die InnerImage (BasicImage<float>)
    * uebergeben wird, also Image (InnerImage i)
    */
    void testInnerImageConstructor()
    {
        std::auto_ptr<Image> image1(Policy::factory(new typename ChildImage::InnerImage(3, 4, child_data[0]))); 
        should(image1->height() == 4);
        should(image1->width() == 3);
        should(image1->size() == vigra::Diff2D(3,4));
        should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[0])));    
        
        std::auto_ptr<Image> image2(Policy::factory(new typename ChildImage::InnerImage(0, 0, child_data[1])));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
        should(image2->end() == std::find_if(image2->begin(), image2->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[1])));    
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

        /* Aenderung mit dem Zuweisungsoperator -> siehe in testShallowCopy()
        */
        (*image0_) = (*image1_);                              
        should(image0->size() == vigra::Diff2D(0,0));
        should(image0_->size() == vigra::Diff2D(3,5));
        should(!equal(*image0, *image0_));
        should(equal(*image1_, *image0_));
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
        
        (*image0_) = (*image1_);
        should(equal(*image1_, *image0_));
        should(!equal(*image1_, *image0));
        should(!equal(*image0_, *image0));
    }
    
    /** Testet die Methode actualSize(), die die Groesse des kompletten Bildes liefert und nicht der ROI
    */
    void testActualSize()
    {
        should(image1_->size() == image1_->actualSize());
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        should(image1_->size() != image1_->actualSize());
        should(image1_->size() == Diff2D(1,3));
    }
    
    /** Testet die Methode actualWidth(), die die Breite des kompletten Bildes liefert und nicht der ROI
    */
    void testActualWidth()
    {
        should(image1_->width() == image1_->actualWidth());
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        should(image1_->width() == 1);
        should(image1_->width() != image1_->actualWidth());
    }
    
    /** Testet die Methode actualHeight(), die die Hoehe des kompletten Bildes liefert und nicht der ROI
    */
    void testActualHeight()
    {
        should(image1_->height() == image1_->actualHeight());
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        should(image1_->height() == 3);
        should(image1_->height() != image1_->actualHeight());
    }

    /** Testet die Methode actualBands(), die die Anzahl der Baender eines Bildes liefert 
    */
    void testActualBands()
    {
        should(image1_->bands() == image1_->actualBands());
    }
    
    /** Testet die Methode bands(), die die Anzahl der selektirerten Baender eines Bildes liefert.
    * Die Anzahl der selektirerten Baender ist bei allen Bildern gleich dem Anzahl der Baender des Bildes,
    * die Ausnahme bildet SelectBandImage dort ist immer nur ein Band selektiert.
    */
    void testBands()
    {
        should(image1_->bands() == image1_->actualBands());
    }
    
    /** testet die Methode roiUpperLeft(), die die Koordinaten der linken
    * oberen Ecke von ROI liefert.
    */
    void testRoiUpperLeft()
    {
        should(image1_->roiUpperLeft() == Diff2D(0, 0));
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        should(image1_->roiUpperLeft() == Diff2D(1, 1));
        should((*image1_)[image1_->roiUpperLeft()] != (*image1_->upperLeft()));
    }

    /** testet die Methode setROI(), die eine ROI auf einem Bild setzt.
    */
    void testSetROI()
    {
        std::auto_ptr<Image> image1_copy(image1_->shallowCopy());               //dient der Vergleichsmoeglichkeit
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        
        should(image1_->width() == (image1_->actualWidth() - 2));
        should(image1_->height() == (image1_->actualHeight() - 2));
        should(image1_->size() == (image1_->actualSize() - Diff2D(2, 2)));
        should((*image1_copy)(1,1) == (*image1_)(0,0));
        
        /*  Differenz zweier ScanOrderIteratoren ist eine int_Zahl
        */
        should((image1_->end() - image1_->begin()) == 3);                       // Anzahl der Pixel in der ROI
        should((image1_copy->end() - image1_copy->begin()) == 15);              // Anzahl der Pixel im Hauptbild
                
        /*  Erzeugt man eine shallowCopy eines Bildes an dem ROI gesetz ist, dann ist
        *   die shallowCopy eine Copy des gesamten Bildes und nicht nur der ROI
        *   dabei ist an der shallowCopy genau so die ROI gesetzt.
        */
        std::auto_ptr<Image> image1_ROI_copy(image1_->shallowCopy());
         
        should(!(image1_ROI_copy->size() == image1_ROI_copy->actualSize()));
        should(image1_ROI_copy->size() == image1_->size());
        
        /*  Erzeugt man einen Clone eines Bildes an dem ROI gesetz ist, dann kriegt man ein Bild
        *   mit der Groesse und Pixelinitialisierung der ROI
        */
        std::auto_ptr<typename Image::CloneType> image1_ROI_clone(image1_->clone());
        
        should(image1_ROI_clone->size() == image1_ROI_clone->actualSize());
        should(image1_ROI_clone->size() == image1_->size());
                   
        for(int x = 0; x < image1_ROI_clone->width(); x++)
            for(int y = 0; y < image1_ROI_clone->height(); y++)
                should((*image1_ROI_clone)(x,y) == (*image1_copy)(x + 1,y + 1));
                
          
        /*  Wenn man die Aenderungen an der ROI vornimmt, aendert sich auch die shallowCopy des Hauptbildes
        *   aber verschoben um roiUpperLeft
        */
        image1_->init(data[7]);
        
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
                    
        for(int x = 0; x < image1_->width(); x++)
            for(int y = 0; y < image1_->height(); y++)
                should((*image1_)(x,y) == static_cast<typename Policy::PixelType>(data[7]));

                 
        for(int x = 1; x < (image1_->width() + 1); x++)                                 //um 1 verschoben
            for(int y = 1; y < (image1_->height() + 1); y++)                            //um 1 verschoben
                should((*image1_copy)(x,y) == static_cast<typename Policy::PixelType>(data[7]));
                
        
        for(int x = 0; x < image1_->width(); x++)                                       //Verhalten einer ROI
            for(int y = 0; y < image1_->height(); y++)
                should((*image1_ROI_copy)(x,y) == static_cast<typename Policy::PixelType>(data[7]));
                
        
        /* Es wird ein neues Bild erzeugt mit unterschiedlichen Pixeldaten
        */
        std::auto_ptr<Image> image1_new(Policy::factory(3, 5));
        typename Image::Accessor acc;
        acc = image1_new->accessor();
        scanOrderIter1_ = image1_new->begin();
        for(unsigned int i = 0; i < data.size(); i++)
            {   
                acc.set(data[i], scanOrderIter1_);
                ++scanOrderIter1_; 
            }
        
        std::auto_ptr<Image> image1_new_copy(image1_new->shallowCopy());
        
        /*  Setzt man eine ROI auf der shallowCopy so ist sie auch nur in der Copy gesetzt.
        *   Die Veraenderungen in der ROI sind aber in der Copie und im image1_new sichtbar,
        *   im image1_new nur verschoben
        */
        image1_new_copy->setROI(Diff2D(1,1), Diff2D(1,3));
        
        should(image1_new->size() == image1_new->actualSize());
        should(image1_new_copy->size() != image1_new_copy->actualSize());
        
        (*image1_new_copy)(0,0) = data[14];
        
        should((*image1_new_copy)(0,0) == static_cast<typename Policy::PixelType>(data[14]));
        should((*image1_new)(1,1) == (*image1_new_copy)(0,0));
        
        /*  Von der ROI kann man in den negativen Bereich zurueckgreifen, da dort das Bild immernoch definiert ist
        *   Vorsicht! der negative Zugrif darf nicht ausserhalb des eigentlichen Bildes erfolgen
        */
        should((image1_new_copy->upperLeft())[Diff2D(-1,-1)] == (*image1_new)(0,0));
        
        /*  Eine ROI ist auch auf die Groesse des Bildes ruecksetzbar
        *   Vorsicht! Die Grenzueberschreitungen beachten!
        */
        image1_new_copy->setROI(Diff2D(-1,-1), image1_new_copy->actualSize());                      //Ruecksetzen
        
        should(image1_new_copy->size() == image1_new_copy->actualSize());
        should((image1_new->upperLeft())[Diff2D(0,0)] == (*image1_new_copy)[image1_new_copy->roiUpperLeft()]);
        
        
        /*  Setzt man die ROI mit neagtiven Werten in size so wird der negative Wert vom lowerRight des Bildes abgezogen
        *   in unserem Falle sind Folgende ROIsetzungen gleich
        */
        image1_new_copy->setROI(Diff2D(1,1), Diff2D(-1,-1));
        should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
        should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
        should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
        
        image1_new_copy->setROI(Diff2D(-1,-1), image1_new_copy->actualSize());                      //Ruecksetzen
        image1_new_copy->setROI(Diff2D(1,1), Diff2D(1,-1));                                         //neu setzen
        
        should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
        should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
        should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
        
        image1_new_copy->resetROI();                                                                //Ruecksetzen
        image1_new_copy->setROI(Diff2D(1,1), Diff2D(-1,3));                                         //neu setzen
        
        should((*image1_new_copy)(image1_new_copy->width() - 1, image1_new_copy->height() - 1) == (*image1_new)(1,3));
        should(image1_new_copy->lowerRight()[Diff2D(0,0)] == (*image1_new)(image1_new->width() - 1,image1_new->height() - 1));
        should(equalIteratorRun(image1_new_copy->upperLeft(), image1_new_copy->lowerRight(), image1_new_copy->begin()));
     }//end of testSetRoi()
     
     /** testet die Methode resetROI(), die die gesetzte ROI wieder aufhebt
     */
     void testResetROI()
     {
        std::auto_ptr<typename Image::CloneType> image1(image1_->clone()); 
        image1_->setROI(Diff2D(1,1), Diff2D(-1,3));
        image1_->resetROI();  
        should(equal(*image1, *image1_));
     }
     
     /** Testet die setROI(int band) Methode der Klasse SelectBandImage. Die Methode aendert den 
     * selektierten Band. 
     */
     void testSetROIBand()
     {
        should(image1_->getSelectedBand() == Policy::n);
        image1_->setROI(Policy::n == 0 ? 1 : (Policy::n - 1));
        should(image1_->getSelectedBand() != Policy::n);
        should(image1_->getSelectedBand() == Policy::n == 0 ? 1 : (Policy::n - 1));
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
        should(image1_->getSelectedBand() == Policy::n == 0 ? 1 : (Policy::n - 1));
     }
         
     void testIsInsideROI()
     {
        should(image1_->isInsideROI((image1_->actualSize())- Diff2D(1,1)));
        
        image1_->setROI(Diff2D(1,1), Diff2D(1,3));
        
        bool ergebnis = true;
        for ( int i = 1 ; i < 2 ; i++)
            for (int j = 1; j < 4; j ++)
                ergebnis &= image1_->isInsideROI(vigra::Diff2D(i,j));

        ergebnis &= !image1_->isInsideROI(vigra::Diff2D(0,0)) &&
                     !image1_->isInsideROI(vigra::Diff2D(2,0)) &&
                     !image1_->isInsideROI(vigra::Diff2D(0,4)) &&
                     !image1_->isInsideROI(vigra::Diff2D(2,4)) &&
                     !image1_->isInsideROI(vigra::Diff2D(2,2)) &&
                     !image1_->isInsideROI(vigra::Diff2D(0,2)) &&
                     !image1_->isInsideROI(vigra::Diff2D(2,0)); 
                     
        should(ergebnis);
     }
};// end of class ImageHierarchyTest

template <class POLICY>
struct ImageHierarchyTestSuite
: public ImageTestSuite<POLICY>
{
    ImageHierarchyTestSuite(const char * to_test_Image_name)
    : ImageTestSuite<POLICY>(to_test_Image_name)
    {
        add( testCase( &ImageHierarchyTest<POLICY>::testInnerImageConstructor));
        add( testCase( &ImageHierarchyTest<POLICY>::testClone));
        add( testCase( &ImageHierarchyTest<POLICY>::testShallowCopy));
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
};//end of struct ImageHierarchyTestSuite
#endif // IMAGEHIERARCHY_TEST_HXX
