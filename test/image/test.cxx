#include <iostream>
#include "unittest.h"
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/impex.hxx"

struct BImageTest
{
	typedef BImage Image;

	BImageTest()
	: img(3,3)
	{
		Image::Accessor acc = img.accessor();
		Image::ScanOrderIterator i = img.begin();

		acc.set(1, i);
		++i;
		acc.set(2, i);
		++i;
		acc.set(3, i);
		++i;
		acc.set(4, i);
		++i;
		acc.set(5, i);
		++i;
		acc.set(6, i);
		++i;
		acc.set(7, i);
		++i;
		acc.set(8, i);
		++i;
		acc.set(9, i);
		++i;
		should(i == img.end());
	}

	void scanImage()
	{
		Image::Iterator y = img.upperLeft();
		Image::Iterator lr = img.lowerRight();
		Image::Iterator x = y;
		Image::Accessor acc = img.accessor();

		should(acc(x) == 1);
		++x.x;
		should(acc(x) == 2);
		++x.x;
		should(acc(x) == 3);
		++x.x;
		should(x.x == lr.x);

		++y.y;
		x = y;
		should(acc(x) == 4);
		++x.x;
		should(acc(x) == 5);
		++x.x;
		should(acc(x) == 6);
		++y.y;
		x = y;
		should(acc(x) == 7);
		++x.x;
		should(acc(x) == 8);
		++x.x;
		should(acc(x) == 9);
		++y.y;
		should(y.y == lr.y);

		y = img.upperLeft();
		should(acc(y, Diff2D(1,1)) == 5);

		Image::ScanOrderIterator i = img.begin();
		should(acc(i, 4) == 5);
	}

	void copyImage()
	{
		Image img1 = img;
		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	void storeLoadViffImage()
	{
                ViffImage * viff = createViffImage(srcImageRange(img));
		writeViffImage("test.xv", viff);

		freeViffImage(viff);
		viff = 0;

		viff = readViffImage("test.xv");

		should(viff != 0);
		should(viff->row_size == 3);
		should(viff->col_size == 3);
		should(viff->num_data_bands == 1);
		should(viff->data_storage_type == VFF_TYP_1_BYTE);

		Image img1(3,3);

		importViffImage(viff, destImage(img1));
		freeViffImage(viff);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}


#if 0
	void storeLoadGifImage()
	{
                LugImage * lug = createLugImage(srcImageRange(img));
		writeLugImage("test.gif", lug);

		freeLugImage(lug);
		lug = 0;

		lug = readLugImage("test.gif");

		should(lug != 0);
		should(lug->xsize == 3);
		should(lug->ysize == 3);

		Image img1(3,3);

		importLugImage(lug, destImage(img1));
		freeLugImage(lug);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}
#endif /* #if 0 */

	void storeLoadGifImage()
	{
                exportImage(srcImageRange(img), ImageExportInfo("test.gif"));
                
                ImageImportInfo info("test.gif");

		should(info.isGrayscale());
                should(info.width() == 3);
		should(info.height() == 3);

		Image img1(3,3);

		importImage(info, destImage(img1));

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	Image img;
};

struct DImageTest
{
	typedef DImage Image;

	DImageTest()
	: img(3,3)
	{
		Image::Accessor acc = img.accessor();
		Image::ScanOrderIterator i = img.begin();

		acc.set(1.1, i);
		++i;
		acc.set(2.2, i);
		++i;
		acc.set(3.3, i);
		++i;
		acc.set(4.4, i);
		++i;
		acc.set(5.5, i);
		++i;
		acc.set(6.6, i);
		++i;
		acc.set(7.7, i);
		++i;
		acc.set(8.8, i);
		++i;
		acc.set(9.9, i);
		++i;
		should(i == img.end());
	}

	void scanImage()
	{
		Image::Iterator y = img.upperLeft();
		Image::Iterator lr = img.lowerRight();
		Image::Iterator x = y;
		Image::Accessor acc = img.accessor();

		should(acc(x) == 1.1);
		++x.x;
		should(acc(x) == 2.2);
		++x.x;
		should(acc(x) == 3.3);
		++x.x;
		should(x.x == lr.x);

		++y.y;
		x = y;
		should(acc(x) == 4.4);
		++x.x;
		should(acc(x) == 5.5);
		++x.x;
		should(acc(x) == 6.6);
		++y.y;
		x = y;
		should(acc(x) == 7.7);
		++x.x;
		should(acc(x) == 8.8);
		++x.x;
		should(acc(x) == 9.9);
		++y.y;
		should(y.y == lr.y);

		y = img.upperLeft();
		should(acc(y, Diff2D(1,1)) == 5.5);

		Image::ScanOrderIterator i = img.begin();
		should(acc(i, 4) == 5.5);
	}

	void copyImage()
	{
		Image img1 = img;
		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	void storeLoadImage()
	{
                ViffImage * viff = createViffImage(srcImageRange(img));
		writeViffImage("test.xv", viff);

		freeViffImage(viff);
		viff = 0;

		viff = readViffImage("test.xv");

		should(viff != 0);
		should(viff->row_size == 3);
		should(viff->col_size == 3);
		should(viff->num_data_bands == 1);
		should(viff->data_storage_type == VFF_TYP_DOUBLE);

		Image img1(3,3);

		importViffImage(viff, destImage(img1));
		freeViffImage(viff);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	Image img;
};

struct BRGBImageTest
{
	typedef BRGBImage Image;

	BRGBImageTest()
	: img(3,3), c1(1, 2, 3), c2(4, 5, 6), c3(7, 8, 9)
	{
		Image::Accessor acc = img.accessor();
		Image::ScanOrderIterator i = img.begin();

		acc.set(c1, i);
		++i;
		acc.set(c2, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c2, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c1, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c1, i);
		++i;
		acc.set(c2, i);
		++i;
		should(i == img.end());
	}

	void scanImage()
	{
		Image::Iterator y = img.upperLeft();
		Image::Iterator lr = img.lowerRight();
		Image::Iterator x = y;
		Image::Accessor acc = img.accessor();

		should(acc(x) == c1);
		++x.x;
		should(acc(x) == c2);
		++x.x;
		should(acc(x) == c3);
		++x.x;
		should(x.x == lr.x);

		++y.y;
		x = y;
		should(acc(x) == c2);
		++x.x;
		should(acc(x) == c3);
		++x.x;
		should(acc(x) == c1);
		++y.y;
		x = y;
		should(acc(x) == c3);
		++x.x;
		should(acc(x) == c1);
		++x.x;
		should(acc(x) == c2);
		++y.y;
		should(y.y == lr.y);

		y = img.upperLeft();
		should(acc(y, Diff2D(1,1)) == c3);

		Image::ScanOrderIterator i = img.begin();
		should(acc(i, 4) == c3);
	}

	void copyImage()
	{
		Image img1 = img;
		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	void storeLoadViffImage()
	{
                ViffImage * viff = createViffImage(srcImageRange(img));
		writeViffImage("test.xv", viff);

		freeViffImage(viff);
		viff = 0;

		viff = readViffImage("test.xv");

		should(viff != 0);
		should(viff->row_size == 3);
		should(viff->col_size == 3);
		should(viff->num_data_bands == 3);
		should(viff->data_storage_type == VFF_TYP_1_BYTE);

		Image img1(3,3);

		importViffImage(viff, destImage(img1));
		freeViffImage(viff);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}


#if 0
	void storeLoadGifImage()
	{
		LugImage * lug = readLugImage("lenna.gif");

		should(lug != 0);

		Image img1(lug->xsize, lug->ysize);
                
		importLugImage(lug, destImage(img1));
		freeLugImage(lug);
                
                lug = createLugImage(srcImageRange(img1));
                write_rgb_file("test.rgb", lug);

                read_rgb_file("test.rgb", lug, img1.width(), img1.height());
                ViffImage * viff = readViffImage("lenna.xv");
                
		should(lug != 0);
		should(viff != 0);
		should(lug->xsize == viff->row_size);
		should(lug->ysize == viff->col_size);

		Image img2(lug->xsize, lug->ysize);

		importLugImage(lug, destImage(img1));
		freeLugImage(lug);

		importViffImage(viff, destImage(img2));
		freeViffImage(viff);

		Image::ScanOrderIterator i1 = img1.begin();
		Image::ScanOrderIterator i2 = img2.begin();
		Image::Accessor acc = img1.accessor();

		for(; i1 != img1.end(); ++i1, ++i2)
		{
			should(acc(i1) == acc(i2));
		}
	}

#endif /* #if 0 */

	void storeLoadGifImage()
	{
		ImageImportInfo info("lenna.gif");
                
                should(info.isColor());

		Image img1(info.width(), info.height());
                
		importImage(info, destImage(img1));
               
                exportImage(srcImageRange(img1), ImageExportInfo("test.ras"));

                ViffImage * viff = readViffImage("lenna.xv");
		should(viff != 0);
                
                ImageImportInfo info1("test.ras");
                
		should(info1.isColor());
		should(info1.width() == viff->row_size);
		should(info1.height() == viff->col_size);
		should(info1.width() == info.width());
		should(info1.height() == info.height());

		Image img2(info1.width(), info1.height());

		importImage(info1, destImage(img1));

		importViffImage(viff, destImage(img2));
		freeViffImage(viff);

		Image::ScanOrderIterator i1 = img1.begin();
		Image::ScanOrderIterator i2 = img2.begin();
		Image::Accessor acc = img1.accessor();

		for(; i1 != img1.end(); ++i1, ++i2)
		{
			should(acc(i1) == acc(i2));
		}
	}

	Image img;
	Image::value_type c1, c2, c3;
};

struct FRGBImageTest
{
	typedef FRGBImage Image;

	FRGBImageTest()
	: img(3,3), c1(1.1, 2.2, 3.3), c2(4.4, 5.5, 6.6), c3(7.7, 8.8, 9.9)
	{
		Image::Accessor acc = img.accessor();
		Image::ScanOrderIterator i = img.begin();

		acc.set(c1, i);
		++i;
		acc.set(c2, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c2, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c1, i);
		++i;
		acc.set(c3, i);
		++i;
		acc.set(c1, i);
		++i;
		acc.set(c2, i);
		++i;
		should(i == img.end());
	}

	void scanImage()
	{
		Image::Iterator y = img.upperLeft();
		Image::Iterator lr = img.lowerRight();
		Image::Iterator x = y;
		Image::Accessor acc = img.accessor();

		should(acc(x) == c1);
		++x.x;
		should(acc(x) == c2);
		++x.x;
		should(acc(x) == c3);
		++x.x;
		should(x.x == lr.x);

		++y.y;
		x = y;
		should(acc(x) == c2);
		++x.x;
		should(acc(x) == c3);
		++x.x;
		should(acc(x) == c1);
		++y.y;
		x = y;
		should(acc(x) == c3);
		++x.x;
		should(acc(x) == c1);
		++x.x;
		should(acc(x) == c2);
		++y.y;
		should(y.y == lr.y);

		y = img.upperLeft();
		should(acc(y, Diff2D(1,1)) == c3);

		Image::ScanOrderIterator i = img.begin();
		should(acc(i, 4) == c3);
	}

	void copyImage()
	{
		Image img1 = img;
		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	void storeLoadImage()
	{
        ViffImage * viff = createViffImage(srcImageRange(img));
		writeViffImage("test.xv", viff);

		freeViffImage(viff);
		viff = 0;

		viff = readViffImage("test.xv");

		should(viff != 0);
		should(viff->row_size == 3);
		should(viff->col_size == 3);
		should(viff->num_data_bands == 3);
		should(viff->data_storage_type == VFF_TYP_FLOAT);

		Image img1(3,3);

		importViffImage(viff, destImage(img1));
		freeViffImage(viff);

		Image::ScanOrderIterator i = img.begin();
		Image::ScanOrderIterator i1 = img1.begin();
		Image::Accessor acc = img.accessor();

		for(; i != img.end(); ++i, ++i1)
		{
			should(acc(i) == acc(i1));
		}
	}

	Image img;
	Image::value_type c1, c2, c3;
};

struct ImageTestSuite
: public TestSuite
{
	ImageTestSuite()
	: TestSuite("ImageTestSuite")
	{
		add( testCase( &BImageTest::scanImage));
		add( testCase( &BImageTest::copyImage));
		add( testCase( &BImageTest::storeLoadViffImage));
		add( testCase( &BImageTest::storeLoadGifImage));
		add( testCase( &DImageTest::scanImage));
		add( testCase( &DImageTest::copyImage));
		add( testCase( &DImageTest::storeLoadImage));
		add( testCase( &BRGBImageTest::scanImage));
		add( testCase( &BRGBImageTest::copyImage));
		add( testCase( &BRGBImageTest::storeLoadViffImage));
		add( testCase( &BRGBImageTest::storeLoadGifImage));
		add( testCase( &FRGBImageTest::scanImage));
		add( testCase( &FRGBImageTest::copyImage));
		add( testCase( &FRGBImageTest::storeLoadImage));
	}
};

int main()
{
    ImageTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

