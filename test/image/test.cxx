#include <iostream>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/imageiterator.hxx"
#include "vigra/impex.hxx"
#include "vigra/imagecontainer.hxx"
#include "vigra/inspectimage.hxx"

using vigra::Diff2D;
using namespace vigra;

template <class IMAGE>
struct ImageTest
{
    typedef IMAGE Image;
    static typename Image::value_type data[];

    ImageTest()
    : img(3,3)
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

    void testBasicImageIterator()
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

    void testImageIterator()
    {
        vigra::ImageIterator<typename Image::value_type> ul(img.begin(), img.width());
        vigra::ImageIterator<typename Image::value_type> lr = ul + img.size();
        scanImage(ul, lr);
        scanRows(ul.rowIterator(), (ul+Diff2D(0,1)).rowIterator(),
                 (ul+Diff2D(0,2)).rowIterator(), img.width());
        scanColumns(ul.columnIterator(), (ul+Diff2D(1,0)).columnIterator(),
                 (ul+Diff2D(2,0)).columnIterator(), img.height());
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

typedef ImageTest<vigra::BImage> BImageTest;

unsigned char BImageTest::data[] = {1,2,3,4,5,6,7,8,9};

typedef ImageTest<vigra::DImage> DImageTest;

double DImageTest::data[] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9};

typedef ImageTest<vigra::BRGBImage> BRGBImageTest;
typedef vigra::RGBValue<unsigned char> BRGB;
BRGB BRGBImageTest::data[] = {
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

typedef ImageTest<vigra::FRGBImage> FRGBImageTest;
typedef vigra::RGBValue<float> FRGB;
FRGB FRGBImageTest::data[] = {
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


struct ImageTestSuite
: public vigra::test_suite
{
    ImageTestSuite()
    : vigra::test_suite("ImageTestSuite")
    {
        add( testCase( &BImageTest::testBasicImageIterator));
        add( testCase( &BImageTest::testImageIterator));
        add( testCase( &BImageTest::copyImage));
        add( testCase( &DImageTest::testBasicImageIterator));
        add( testCase( &DImageTest::testImageIterator));
        add( testCase( &DImageTest::copyImage));
        add( testCase( &BRGBImageTest::testBasicImageIterator));
        add( testCase( &BRGBImageTest::testImageIterator));
        add( testCase( &BRGBImageTest::copyImage));
        add( testCase( &FRGBImageTest::testBasicImageIterator));
        add( testCase( &FRGBImageTest::testImageIterator));
        add( testCase( &FRGBImageTest::copyImage));
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

// for shouldEqual:
std::ostream & operator <<(std::ostream &s, Diff2D size)
{
	s << "Diff2D(" << size.x << ", " << size.y << ")";
    return s;
}

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

		for (int i=0; i<ia.size(); i++)
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

int main()
{
    ImageTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

	ImageContainerTestSuite containerTest;
	
	failed = failed || containerTest.run();

	std::cout << containerTest.report() << std::endl;

    return (failed != 0);
}

