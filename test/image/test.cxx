#include <iostream>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/basicimageview.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/impex.hxx"
#include "vigra/imagecontainer.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/pixelneighborhood.hxx"
#include "vigra/contourcirculator.hxx"

using vigra::Diff2D;
using namespace vigra;

unsigned char * testData(unsigned char)
{
    static unsigned char data[] = {1,2,3,4,5,6,7,8,9};
    return data;
}

double * testData(double)
{
    static double data[] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9};
    return data;
}

RGBValue<unsigned char> * testData(RGBValue<unsigned char>)
{
    typedef RGBValue<unsigned char> BRGB;
    static BRGB data[] = {
        BRGB(1,1,1),
        BRGB(2,2,2),
        BRGB(3,3,3),
        BRGB(4,4,4),
        BRGB(5,5,5),
        BRGB(6,6,6),
        BRGB(7,7,7),
        BRGB(8,8,8),
        BRGB(9,9,9)
    };
    return data;
}

RGBValue<float> * testData(RGBValue<float>)
{
    typedef vigra::RGBValue<float> FRGB;
    static FRGB data[] = {
        FRGB(1.1, 1.1, 1.1),
        FRGB(2.2, 2.2, 2.2),
        FRGB(3.3, 3.3, 3.3),
        FRGB(4.4, 4.4, 4.4),
        FRGB(5.5, 5.5, 5.5),
        FRGB(6.6, 6.6, 6.6),
        FRGB(7.7, 7.7, 7.7),
        FRGB(8.8, 8.8, 8.8),
        FRGB(9.9, 9.9, 9.9)
    };
    return data;
}

template <class IMAGE>
struct ImageTest
{
    typedef IMAGE Image;
    typedef typename Image::value_type value_type;
    value_type internalMemory[9];
    value_type * data;

    ImageTest(IMAGE const & image)
    : img(image), data(testData(value_type()))
    {
        typename Image::Accessor acc = img.accessor();
        typename Image::iterator i = img.begin();
        
        acc.set(data[0], i);
        ++i;
        acc.set(data[1], i);
        ++i;
        acc.set(data[2], i);
        ++i;
        acc.set(data[3], i);
        ++i;
        acc.set(data[4], i);
        ++i;
        acc.set(data[5], i);
        ++i;
        acc.set(data[6], i);
        ++i;
        acc.set(data[7], i);
        ++i;
        acc.set(data[8], i);
        ++i;
        should(i == img.end());
    }

    template <class Iterator>
    void scanImage(Iterator ul, Iterator lr)
    {
        Iterator y = ul;
        Iterator x = ul;
        typename Image::Accessor acc = img.accessor();

        should(acc(x) == data[0]);
        ++x.x;
        should(acc(x) == data[1]);
        ++x.x;
        should(acc(x) == data[2]);
        ++x.x;
        should(x.x == lr.x);

        ++y.y;
        x = y;
        should(acc(x) == data[3]);
        ++x.x;
        should(acc(x) == data[4]);
        ++x.x;
        should(acc(x) == data[5]);
        ++y.y;
        x = y;
        should(acc(x) == data[6]);
        ++x.x;
        should(acc(x) == data[7]);
        ++x.x;
        should(acc(x) == data[8]);
        ++y.y;
        should(y.y == lr.y);

        y = ul;
        should(acc(y, vigra::Diff2D(1,1)) == data[4]);
    }

    template <class Iterator>
    void scanRows(Iterator r1, Iterator r2, Iterator r3, int w)
    {
        Iterator end = r1 + w;
        typename Image::Accessor acc = img.accessor();

        should(acc(r1) == data[0]);
        ++r1;
        should(acc(r1) == data[1]);
        ++r1;
        should(acc(r1) == data[2]);
        ++r1;
        should(r1 == end);

        end = r2 + w;
        should(acc(r2) == data[3]);
        ++r2;
        should(acc(r2) == data[4]);
        ++r2;
        should(acc(r2) == data[5]);
        ++r2;
        should(r2 == end);

        end = r3 + w;
        should(acc(r3) == data[6]);
        ++r3;
        should(acc(r3) == data[7]);
        ++r3;
        should(acc(r3) == data[8]);
        ++r3;
        should(r3 == end);
    }

    template <class Iterator>
    void scanColumns(Iterator c1, Iterator c2, Iterator c3, int h)
    {
        Iterator end = c1 + h;
        typename Image::Accessor acc = img.accessor();

        should(acc(c1) == data[0]);
        ++c1;
        should(acc(c1) == data[3]);
        ++c1;
        should(acc(c1) == data[6]);
        ++c1;
        should(c1 == end);

        end = c2 + h;
        should(acc(c2) == data[1]);
        ++c2;
        should(acc(c2) == data[4]);
        ++c2;
        should(acc(c2) == data[7]);
        ++c2;
        should(c2 == end);

        end = c3 + h;
        should(acc(c3) == data[2]);
        ++c3;
        should(acc(c3) == data[5]);
        ++c3;
        should(acc(c3) == data[8]);
        ++c3;
        should(c3 == end);
    }

    void testIterator()
    {
        typename Image::Iterator ul = img.upperLeft();
        typename Image::Iterator lr = img.lowerRight();

        scanImage(ul, lr);
        scanRows(ul.rowIterator(), (ul+Diff2D(0,1)).rowIterator(),
                 (ul+Diff2D(0,2)).rowIterator(), img.width());
        scanColumns(ul.columnIterator(), (ul+Diff2D(1,0)).columnIterator(),
                 (ul+Diff2D(2,0)).columnIterator(), img.height());

        typename Image::Accessor acc = img.accessor();
        typename Image::iterator i = img.begin();
        should(acc(i, 4) == data[4]);
    }

    void testIndex()
    {
        for (int y=0; y<img.height(); ++y)
        {
            for(int x=0; x<img.width(); ++x)
            {
                shouldEqual(data[y*img.width()+x], img[Diff2D(x,y)]);
                shouldEqual(data[y*img.width()+x], img(x,y));
            }
        }
    }

    void copyImage()
    {
        typedef typename Image::value_type Value;

        Image img1(img);
        typename Image::iterator i = img.begin();
        typename Image::iterator i1 = img1.begin();
        typename Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }

        img.init(0);
        for(; i != img.end(); ++i)
        {
			should(acc(i) == Value(0));
        }
        img(1,1) = Value(200);
        img1 = img;
        i = img.begin();
        i1 = img1.begin();

        for(; i != img.end(); ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }
    }

    Image img;
};

template <class IMAGE>
struct BasicImageTest
: public ImageTest<IMAGE>
{
    BasicImageTest()
    : ImageTest<IMAGE>(IMAGE(Diff2D(3,3)))
    {}

    // next lines needed due to gcc 2.95 bug
    void testIterator()
    {
        ImageTest<IMAGE>::testIterator();
    }

    void testIndex()
    {
        ImageTest<IMAGE>::testIndex();
    }

    void copyImage()
    {
        ImageTest<IMAGE>::copyImage();
    }

    void swapImage()
    {
        IMAGE other(1,1);
        other(0,0) = this->data[2];
        
        this->img.swap(other);
        
        shouldEqual(this->img.width(), 1);
        shouldEqual(this->img.height(), 1);
        shouldEqual(this->img(0,0), this->data[2]);
        shouldEqual(other.width(), 3);
        shouldEqual(other.height(), 3);
        shouldEqual(other(2,2), this->data[8]);
    }
};

template <class IMAGE>
struct BasicImageViewTest
: public ImageTest<IMAGE>
{
    
    BasicImageViewTest()
    : ImageTest<IMAGE>(IMAGE(ImageTest<IMAGE>::internalMemory, Diff2D(3,3)))
    {}

    // next lines needed due to gcc 2.95 bug
    void testIterator()
    {
        ImageTest<IMAGE>::testIterator();
    }

    void testIndex()
    {
        ImageTest<IMAGE>::testIndex();
    }

    void copyImage()
    {
        ImageTest<IMAGE>::copyImage();
    }
};

template <class T>
struct StridedImageIteratorTest
{
    T * data_;
    
    StridedImageIteratorTest()
    : data_(testData(T()))
    {}

    void testIterator()
    {
        int w = 3, h = 3;
        int xskip = 2, yskip = 2;
        int ws = w / xskip + 1, hs = h / yskip + 1;
        
        StridedImageIterator<T> ul(data_, w, xskip, yskip);
        StridedImageIterator<T> lr = ul + Diff2D(ws, hs);
        
        shouldEqual(ws, lr.x - ul.x);
        shouldEqual(hs, lr.y - ul.y);
        shouldEqual(Diff2D(ws, hs), lr - ul);
        
        StridedImageIterator<T> x = ul;
        typename StridedImageIterator<T>::row_iterator r = ul.rowIterator();
        typename StridedImageIterator<T>::row_iterator rend = r + ws;
        shouldEqual(data_[0], *x);
        shouldEqual(data_[0], *r);
        should(x.x < lr.x);
        should(r < rend);
        ++x.x;
        ++r;
        shouldEqual(data_[2], *x);
        shouldEqual(data_[2], *r);
        should(x.x < lr.x);
        should(r < rend);
        ++x.x;
        ++r;
        should(x.x == lr.x);
        should(r == rend);
        
        ++ul.y;
        x = ul;
        r = ul.rowIterator();
        rend = r + ws;
        shouldEqual(data_[6], *x);
        shouldEqual(data_[6], *r);
        should(x.x < lr.x);
        should(r < rend);
        ++x.x;
        ++r;
        shouldEqual(data_[8], *x);
        shouldEqual(data_[8], *r);
        should(x.x < lr.x);
        should(r < rend);
        ++x.x;
        ++r;
        should(x.x == lr.x);
        should(r == rend);

        ++ul.y;
        should(ul.y == lr.y);
    }
};


struct ImageTestSuite
: public vigra::test_suite
{
    ImageTestSuite()
    : vigra::test_suite("ImageTestSuite")
    {
        add( testCase( &BasicImageTest<BasicImage<unsigned char> >::testIterator));
        add( testCase( &BasicImageTest<BasicImage<unsigned char> >::testIndex));
        add( testCase( &BasicImageTest<BasicImage<unsigned char> >::copyImage));
        add( testCase( &BasicImageTest<BasicImage<unsigned char> >::swapImage));
        add( testCase( &BasicImageTest<BasicImage<double> >::testIterator));
        add( testCase( &BasicImageTest<BasicImage<double> >::testIndex));
        add( testCase( &BasicImageTest<BasicImage<double> >::copyImage));
        add( testCase( &BasicImageTest<BasicImage<double> >::swapImage));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<unsigned char> > >::testIterator));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<unsigned char> > >::testIndex));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<unsigned char> > >::copyImage));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<unsigned char> > >::swapImage));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<float> > >::testIterator));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<float> > >::testIndex));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<float> > >::copyImage));
        add( testCase( &BasicImageTest<BasicImage<RGBValue<float> > >::swapImage));
        add( testCase( &BasicImageViewTest<BasicImageView<unsigned char> >::testIterator));
        add( testCase( &BasicImageViewTest<BasicImageView<unsigned char> >::testIndex));
        add( testCase( &BasicImageViewTest<BasicImageView<unsigned char> >::copyImage));
        add( testCase( &BasicImageViewTest<BasicImageView<double> >::testIterator));
        add( testCase( &BasicImageViewTest<BasicImageView<double> >::testIndex));
        add( testCase( &BasicImageViewTest<BasicImageView<double> >::copyImage));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<unsigned char> > >::testIterator));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<unsigned char> > >::testIndex));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<unsigned char> > >::copyImage));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<float> > >::testIterator));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<float> > >::testIndex));
        add( testCase( &BasicImageViewTest<BasicImageView<RGBValue<float> > >::copyImage));
        add( testCase( &StridedImageIteratorTest<unsigned char>::testIterator));
        add( testCase( &StridedImageIteratorTest<RGBValue<float> >::testIterator));
    }
};

struct CompareFunctor
{
	double sumDifference_;

	CompareFunctor(): sumDifference_(0) {}

	void operator()(const float &a, const float &b)
    { sumDifference_+=  VIGRA_CSTD::abs(a-b); }

    double operator()()
		{ return sumDifference_; }
};

struct ImageContainerTests
{
	ImageImportInfo info;
	int w, h;
	FImage lennaImage;

	ImageContainerTests()
		: info("../../images/lenna.xv"),
		  w(info.width()), h(info.height()),
		  lennaImage(w, h)
	{
		importImage(info, destImage(lennaImage));
	}

	void initArrayWithImageTest()
	{
		ImageArray<FImage> threeLennas(3, lennaImage);
		CompareFunctor cmp;
		inspectTwoImages(srcImageRange(threeLennas[0]), srcImage(threeLennas[2]), cmp);
		shouldEqual(cmp(), 0.0);

		Diff2D newsize(50, 50);
		threeLennas.resizeImages(newsize);
		for (ImageArray<FImage>::iterator it= threeLennas.begin();
			 it!= threeLennas.end(); it++)
			shouldEqual((*it).size(), newsize);
	}

	void initArrayWithSizeTest()
	{
		Diff2D testsize(50, 50);
		ImageArray<FImage> ia(6, testsize);

		for (unsigned int i=0; i<ia.size(); i++)
			shouldEqual(ia[i].size(), testsize);

		ImageArray<FImage> ia2(ia.begin(), ia.end());
		shouldEqual(ia2.imageSize(), testsize);

		ia2.erase(ia2.begin()+1);
		shouldEqual(ia2.size(), ia.size()-1);

		ia.clear();
		shouldEqual(ia.size(), 0);
	}
};

struct ImageContainerTestSuite
: public vigra::test_suite
{
    ImageContainerTestSuite()
    : vigra::test_suite("ImageContainerTestSuite")
    {
        add( testCase( &ImageContainerTests::initArrayWithImageTest ) );
        add( testCase( &ImageContainerTests::initArrayWithSizeTest ) );
    }
};

struct NeighborhoodCirculatorTest
{
    typedef vigra::NeighborhoodCirculator<vigra::BImage::Iterator, vigra::EightNeighborCode>
    EightCirculator;
    typedef vigra::NeighborhoodCirculator<vigra::BImage::Iterator, vigra::FourNeighborCode>
    FourCirculator;

    vigra::BImage img;
    EightCirculator eightCirc;
    FourCirculator fourCirc;

    NeighborhoodCirculatorTest()
        : img(4, 4),
          eightCirc(img.upperLeft() + vigra::Diff2D(1, 1)),
          fourCirc(img.upperLeft() + vigra::Diff2D(1, 1))
    {
        for(int y= 0; y<img.height(); y++)
            for(int x= 0; x<img.width(); x++)
                img(x, y)= y*img.width() + x;
    }

    void testInit()
    {
        shouldEqual(*eightCirc, 6);
        should(!eightCirc.isDiagonal());
        shouldEqual(*fourCirc, 6);
        should(!fourCirc.isDiagonal());
    }

    void testEightCirculateForward()
    {
        eightCirc++;
        shouldEqual(*eightCirc, 2);
        eightCirc++;
        shouldEqual(*eightCirc, 1);
        eightCirc++;
        shouldEqual(*eightCirc, 0);
        eightCirc++;
        shouldEqual(*eightCirc, 4);
        eightCirc++;
        shouldEqual(*eightCirc, 8);
        eightCirc++;
        shouldEqual(*eightCirc, 9);
        eightCirc++;
        shouldEqual(*eightCirc, 10);
        eightCirc++;
        shouldEqual(*eightCirc, 6);
    }

    void testEightCirculateReverse()
    {
        eightCirc--;
        shouldEqual(*eightCirc, 10);
        eightCirc--;
        shouldEqual(*eightCirc, 9);
        eightCirc--;
        shouldEqual(*eightCirc, 8);
        eightCirc--;
        shouldEqual(*eightCirc, 4);
        eightCirc--;
        shouldEqual(*eightCirc, 0);
        eightCirc--;
        shouldEqual(*eightCirc, 1);
        eightCirc--;
        shouldEqual(*eightCirc, 2);
        eightCirc--;
        shouldEqual(*eightCirc, 6);
    }

    void testFourCirculateForward()
    {
        fourCirc++;
        shouldEqual(*fourCirc, 1);
        fourCirc++;
        shouldEqual(*fourCirc, 4);
        fourCirc++;
        shouldEqual(*fourCirc, 9);
        fourCirc++;
        shouldEqual(*fourCirc, 6);
    }

    void testFourCirculateReverse()
    {
        fourCirc--;
        shouldEqual(*fourCirc, 9);
        fourCirc--;
        shouldEqual(*fourCirc, 4);
        fourCirc--;
        shouldEqual(*fourCirc, 1);
        fourCirc--;
        shouldEqual(*fourCirc, 6);
    }

    void testIsDiagonal()
    {
        for(int i=0; i<10; i++, eightCirc++, fourCirc++)
        {
            if(i%2)
                should(eightCirc.isDiagonal());
            else
                should(!eightCirc.isDiagonal());
            should(!fourCirc.isDiagonal());
        }
    }

    void testEquality()
    {
        EightCirculator eightCirc2 = eightCirc;
        should(eightCirc == eightCirc2);
        eightCirc2++;
        should(eightCirc != eightCirc2);
        eightCirc2++;
        should(eightCirc != eightCirc2);
        eightCirc++;
        should(eightCirc != eightCirc2);
        eightCirc2--;
        should(eightCirc == eightCirc2);

        FourCirculator fourCirc2(img.upperLeft() + vigra::Diff2D(1, 1));
        should(fourCirc == fourCirc2);
        fourCirc--;
        should(fourCirc != fourCirc2);

        eightCirc2 = eightCirc + 3;
        eightCirc += 3;
        should(eightCirc == eightCirc2);

        fourCirc2 = fourCirc + 3;
        fourCirc += 3;
        should(fourCirc == fourCirc2);
    }

    void testTurning()
    {
        for(int i=0; i<4; i++)
        {
            shouldEqual(*fourCirc, *eightCirc);
            fourCirc.turnRight();
            eightCirc.turnRight();
        }

        for(int i=0; i<4; i++)
        {
            shouldEqual(*fourCirc, *eightCirc);
            fourCirc.turnLeft();
            eightCirc.turnLeft();
        }

        fourCirc.turnLeft();
        fourCirc.turnLeft();
        eightCirc.turnRound();
        shouldEqual(*fourCirc, *eightCirc);

        eightCirc.turnRight();
        eightCirc.turnRight();
        fourCirc.turnRound();
        shouldEqual(*fourCirc, *eightCirc);
    }

    void testMoving()
    {
        eightCirc.swapCenterNeighbor();
        shouldEqual(*eightCirc, 5); // looking west from 6 now
        eightCirc++;
        shouldEqual(*eightCirc, 9);
        eightCirc++;
        shouldEqual(*eightCirc, 10);
        eightCirc++;
        shouldEqual(*eightCirc, 11);
        eightCirc++;
        shouldEqual(*eightCirc, 7); // looking east again

        eightCirc+= 4; // looking west again
        eightCirc.moveCenterToNeighbor();
        shouldEqual(*eightCirc, 4);
    }

    void testMiscellaneous()
    {
        // test direction()
        should(fourCirc.direction() == vigra::FourNeighborCode::East);

        // test operator -
        EightCirculator eightCirc2 = eightCirc;
        for(int i=0; i<7; i++, eightCirc2++)
            shouldEqual(eightCirc2 - eightCirc, i);

        // test operator[]
        should(eightCirc2[1] == *eightCirc);

        // test base()
        eightCirc += vigra::EightNeighborCode::SouthEast - eightCirc.direction();
        eightCirc.moveCenterToNeighbor();
        should(eightCirc.base() - img.upperLeft() == vigra::Diff2D(3, 3));
        eightCirc.turnRound();
        eightCirc.swapCenterNeighbor();
        should(eightCirc.base() - img.upperLeft() == vigra::Diff2D(2, 2));
        eightCirc.turnRound();
        should(eightCirc.base() - img.upperLeft() == vigra::Diff2D(0, 0));
    }
};

struct CrackContourCirculatorTest
{
    typedef vigra::CrackContourCirculator<vigra::BImage::traverser> CrackCirc;

    vigra::BImage img;

    CrackContourCirculatorTest()
            : img(8,8)
    {
        static int imdata[] = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 1, 0, 0, 1, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 0,
            0, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 0, 1, 1, 0, 0, 0,
            0, 0, 1, 1, 1, 1, 1, 0,
            0, 0, 1, 1, 1, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        };
    
        for(int i = 0; i<img.width()*img.height(); ++i)
            img.begin()[i] = imdata[i];
    }

    void testInit()
    {
        CrackCirc crackCirc(img.upperLeft() + vigra::Diff2D(1, 1));
        CrackCirc end = crackCirc;
        
        should(crackCirc.pos() == vigra::Diff2D(0, 0));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(0, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(0, 2));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(0, 3));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(1, 3));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(2, 3));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(2, 4));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(1, 4));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(1, 5));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(1, 6));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(2, 6));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(3, 6));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 6));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 6));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 5));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(6, 5));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(6, 4));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 4));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 4));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 3));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 2));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 2));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(6, 2));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(6, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(5, 0));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 0));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(4, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(3, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(2, 1));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(2, 0));
        ++crackCirc;
        should(crackCirc.pos() == vigra::Diff2D(1, 0));
        ++crackCirc;
        should(crackCirc == end);
    }

};

struct NeighborhoodCirculatorTestSuite
: public vigra::test_suite
{
    NeighborhoodCirculatorTestSuite()
    : vigra::test_suite("NeighborhoodCirculatorTestSuite")
    {
        add( testCase( &NeighborhoodCirculatorTest::testInit));
        add( testCase( &NeighborhoodCirculatorTest::testEightCirculateForward));
        add( testCase( &NeighborhoodCirculatorTest::testEightCirculateReverse));
        add( testCase( &NeighborhoodCirculatorTest::testFourCirculateForward));
        add( testCase( &NeighborhoodCirculatorTest::testFourCirculateReverse));
        add( testCase( &NeighborhoodCirculatorTest::testIsDiagonal));
        add( testCase( &NeighborhoodCirculatorTest::testEquality));
        add( testCase( &NeighborhoodCirculatorTest::testTurning));
        add( testCase( &NeighborhoodCirculatorTest::testMoving));
        add( testCase( &NeighborhoodCirculatorTest::testMiscellaneous));
   }
};

struct CrackContourCirculatorTestSuite
: public vigra::test_suite
{
    CrackContourCirculatorTestSuite()
    : vigra::test_suite("CrackContourCirculatorTestSuite")
    {
        add( testCase( &CrackContourCirculatorTest::testInit));
   }
};

struct ImageTestCollection
: public vigra::test_suite
{
    ImageTestCollection()
    : vigra::test_suite("ImageTestCollection")
    {
        add( new ImageTestSuite);
        add( new ImageContainerTestSuite);
        add( new NeighborhoodCirculatorTestSuite);
        add( new CrackContourCirculatorTestSuite);
   }
};

int main()
{
    ImageTestCollection test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

