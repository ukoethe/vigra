#ifndef PARENT_TEST_CLASS_HXX
#define PARENT_TEST_CLASS_HXX

#include <unittest.h>

using vigra::Diff2D;
// 
// std::ostream & operator<<(std::ostream & o, vigra::ConstVectorProxy const & p)
// {
//     if(p.size() == 1)
//     {
//         o << p[0];
//     }
//     else
//     {
//         o << "(";
//         int i;
//         for (i=0; i<p.size()-1; ++i)
//             o << p[i] << ", ";
//         o << p[i] << ")";
//     }
//     return o;
// }

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
		should(image1->end() == std::find_if(image1->begin(), image1->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[0])));    
        
        std::auto_ptr<Image> image2(Policy::factory(0, 0, child_data[1]));
        should(image2->height() == 0);
        should(image2->width() == 0);
        should(image2->size() == vigra::Diff2D(0,0));
        should(image2->end() == std::find_if(image2->begin(), image2->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), child_data[1])));    
    }
    
    /** Testet den Copy Konstruktor ( Image(Image img) ).
    */
    void testCopyConstructor()
    {
        ChildImage * image1 = new ChildImage(3, 5, child_data[2]);
        std::auto_ptr<Image> image1_copy(Policy::factory(*image1));

        ChildImage image0(0, 0, child_data[3]);
        std::auto_ptr<Image> image0_copy(Policy::factory(image0));

        should(image1_copy->height() == 5);
        should(image1_copy->width() == 3);
        should(image1_copy->size() == vigra::Diff2D(3,5));
        should(equalPixels(image1->begin(),image1->end(), image1_copy->begin()));
        should(image0_copy->height() == 0);
        should(image0_copy->width() == 0);
        should(image0_copy->size() == vigra::Diff2D(0,0));
        should(equalPixels(image0_->begin(), image0_->end(), image0_copy->begin()));
     }
     
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
    
    /** Testet die init() Methode
    */
    void testInit()
    {
        image1_->init(data[6]);
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
        image1_->init(data[7]);
        should(image1_->end() == std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[7])));
        should(image1_->end() != std::find_if(image1_->begin(), image1_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[6])));
        
        image0_->init(data[8]);
        should(image0_->end() == std::find_if(image0_->begin(), image0_->end(), std::bind2nd(Pixels_not_equal_to<value_type>(), data[8])));
    }
    
    void testWidth()
    {
        should(3 == image1_->width());
        should(0 == image0_->width());
    } 
    
    void testHeight()
    {
        should(5 == image1_->height());
        should(0 == image0_->height());
    } 
    
    void testSize()
    {
        should(image1_->size() == vigra::Diff2D(3,5));
        should(image0_->size() == vigra::Diff2D(0,0));
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
       
};// end of class ImageTest

template <class POLICY>
struct ImageTestSuite
: public vigra::test_suite
{
    ImageTestSuite(const char * to_test_Image_name)
    : vigra::test_suite(to_test_Image_name)
    {
        add( testCase( &ImageTest<POLICY>::testImageDefaultConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntConstuctor));
        add( testCase( &ImageTest<POLICY>::testImage2DConstuctor));
        add( testCase( &ImageTest<POLICY>::testImageIntPixelConstuctor));
        add( testCase( &ImageTest<POLICY>::testCopyConstructor));
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
    }
};//end of struct TImageTestSuite
#endif // PARENT_TEST_CLASS_HXX
