#include <iostream>
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/impex.hxx"
#include "unittest.h"


class ByteImageExportImportTest
{
    typedef BImage Image;
  public:
    ByteImageExportImportTest()
    {
        ViffImage * viff = readViffImage("lenna.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        int w = viff->row_size;
        int h = viff->col_size;

        img.resize(w, h);

        importViffImage(viff, destImage(img));

        freeViffImage(viff);
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.gif"));
        
        ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }

    void testEmptyGIF()
    {
        Image img(100,100);
        
        img = 0;
        
        exportImage(srcImageRange(img), ImageExportInfo("res.gif"));
        
        ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = res.begin();
	Image::Accessor acc = res.accessor();

	for(; i != res.end(); ++i)
	{
            should(acc(i) == 0);
	}
    }

    void testJPEG()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.jpg").setCompression("100"));
        
        ImageImportInfo info("res.jpg");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	float sum = 0;
        for(; i != img.end(); ++i, ++i1)
	{
            sum += abs(acc(i) - acc(i1));
	}
        should(sum / (info.width() * info.height()) < 0.1);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.tif").setCompression("RunLength"));
        
        ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testBMP()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.bmp"));
        
        ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testSUN()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.ras"));
        
        ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF1()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.xv"));
        
        ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF2()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.foo").setFileType("VIFF"));
        
        ViffImage * viff = readViffImage("res.foo");
        
        should(viff != 0);
        should(viff->row_size == img.width());
        should(viff->col_size == img.height());
        
        Image res(viff->row_size, viff->col_size);
        
        importViffImage(viff, destImage(res));
        freeViffImage(viff);
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    Image img;
};

class ByteRGBImageExportImportTest
{
    typedef BRGBImage Image;
  public:
    ByteRGBImageExportImportTest()
    {
        ViffImage * viff = readViffImage("lennargb.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        int w = viff->row_size;
        int h = viff->col_size;

        img.resize(w, h);

        importViffImage(viff, destImage(img));

        freeViffImage(viff);
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.gif"));
        
        ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testJPEG()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.jpg").setCompression("100"));
        
        ImageImportInfo info("res.jpg");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	float sum = 0;
        for(; i != img.end(); ++i, ++i1)
	{
            sum += (acc(i) - acc(i1)).magnitude();
	}
        should(sum / (info.width() * info.height()) < 2.0);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.tif").setCompression("RunLength"));
        
        ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testBMP()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.bmp"));
        
        ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testSUN()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.ras"));
        
        ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF1()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.xv"));
        
        ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF2()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.foo").setFileType("VIFF"));
        
        ViffImage * viff = readViffImage("res.foo");
        
        should(viff != 0);
        should(viff->row_size == img.width());
        should(viff->col_size == img.height());
        
        Image res(viff->row_size, viff->col_size);
        
        importViffImage(viff, destImage(res));
        freeViffImage(viff);
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    Image img;
};

class FloatImageExportImportTest
{
    typedef FImage Image;
  public:
    FloatImageExportImportTest()
    {
        ViffImage * viff = readViffImage("lenna.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        int w = viff->row_size;
        int h = viff->col_size;

        img.resize(w, h);

        importViffImage(viff, destImage(img));

        freeViffImage(viff);
        
        viff = readViffImage("lennafloat.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        reread.resize(w, h);

        importViffImage(viff, destImage(reread));

        freeViffImage(viff);
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.gif"));
        
        ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testJPEG()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.jpg").setCompression("100"));
        
        ImageImportInfo info("res.jpg");
        
        should(info.width() == reread.width());
        should(info.height() == reread.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	float sum = 0;
        for(; i != reread.end(); ++i, ++i1)
	{
            sum += abs(acc(i) - acc(i1));
	}
        should(sum / (info.width() * info.height()) < 0.1);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.tif").setCompression("RunLength"));
        
        ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testBMP()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.bmp"));
        
        ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testSUN()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.ras"));
        
        ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.xv"));
        
        ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    Image img, reread;
};

class FloatRGBImageExportImportTest
{
    typedef FRGBImage Image;
  public:
    FloatRGBImageExportImportTest()
    {
        ViffImage * viff = readViffImage("lennargb.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        int w = viff->row_size;
        int h = viff->col_size;

        img.resize(w, h);

        importViffImage(viff, destImage(img));

        freeViffImage(viff);
        
        viff = readViffImage("lennafloatrgb.xv");
        
        shouldMsg(viff != 0, "Unable to read input image");

        reread.resize(w, h);

        importViffImage(viff, destImage(reread));

        freeViffImage(viff);
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.gif"));
        
        ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testJPEG()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.jpg").setCompression("100"));
        
        ImageImportInfo info("res.jpg");
        
        should(info.width() == reread.width());
        should(info.height() == reread.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	float sum = 0;
        for(; i != reread.end(); ++i, ++i1)
	{
            sum += (acc(i) - acc(i1)).magnitude();
	}
        should(sum / (info.width() * info.height()) < 2.0);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    ImageExportInfo("res.tif").setCompression("RunLength"));
        
        ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testBMP()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.bmp"));
        
        ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testSUN()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.ras"));
        
        ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	for(; i != reread.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    void testVIFF()
    {
        exportImage(srcImageRange(img), ImageExportInfo("res.xv"));
        
        ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	for(; i != img.end(); ++i, ++i1)
	{
            should(acc(i) == acc(i1));
	}
    }
    
    Image img, reread;
};

class ImageExportImportFailureTest
{
  public:
    ImageExportImportFailureTest()
    : img(3,3)
    {}
    
    void testGIFExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.gif"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testGIFImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.gif");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    void testJPEGExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.jpg"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testJPEGImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.jpg");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    void testTIFFExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.tif"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testTIFFImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.tif");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    void testBMPExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.bmp"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testBMPImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.bmp");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    void testSUNExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.ras"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testSUNImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.ras");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    void testVIFFExport()
    {
        bool caught = false;
        
        try
        {
            exportImage(srcImageRange(img), ImageExportInfo("intentionalFailure/res.xv"));
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
    
    void testVIFFImport()
    {
        bool caught = false;
        
        try
        {
            ImageImportInfo info("foo.xv");
        }
        catch(PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    BImage img;
};

class ImageImportExportTestSuite
: public TestSuite
{
  public:
    ImageImportExportTestSuite()
    : TestSuite("ImageImportExportTestSuite")
    {
        add( testCase(&ByteImageExportImportTest::testGIF));
        add( testCase(&ByteImageExportImportTest::testEmptyGIF));
        add( testCase(&ByteImageExportImportTest::testJPEG));
        add( testCase(&ByteImageExportImportTest::testTIFF));
        add( testCase(&ByteImageExportImportTest::testBMP));
        add( testCase(&ByteImageExportImportTest::testSUN));
        add( testCase(&ByteImageExportImportTest::testVIFF1));
        add( testCase(&ByteImageExportImportTest::testVIFF2));

        add( testCase(&ByteRGBImageExportImportTest::testGIF));
        add( testCase(&ByteRGBImageExportImportTest::testJPEG));
        add( testCase(&ByteRGBImageExportImportTest::testTIFF));
        add( testCase(&ByteRGBImageExportImportTest::testBMP));
        add( testCase(&ByteRGBImageExportImportTest::testSUN));
        add( testCase(&ByteRGBImageExportImportTest::testVIFF1));
        add( testCase(&ByteRGBImageExportImportTest::testVIFF2));

        add( testCase(&FloatImageExportImportTest::testGIF));
        add( testCase(&FloatImageExportImportTest::testJPEG));
        add( testCase(&FloatImageExportImportTest::testTIFF));
        add( testCase(&FloatImageExportImportTest::testBMP));
        add( testCase(&FloatImageExportImportTest::testSUN));
        add( testCase(&FloatImageExportImportTest::testVIFF));

        add( testCase(&FloatRGBImageExportImportTest::testGIF));
        add( testCase(&FloatRGBImageExportImportTest::testJPEG));
        add( testCase(&FloatRGBImageExportImportTest::testTIFF));
        add( testCase(&FloatRGBImageExportImportTest::testBMP));
        add( testCase(&FloatRGBImageExportImportTest::testSUN));
        add( testCase(&FloatRGBImageExportImportTest::testVIFF));

        add( testCase(&ImageExportImportFailureTest::testGIFExport));
        add( testCase(&ImageExportImportFailureTest::testGIFImport));
        add( testCase(&ImageExportImportFailureTest::testJPEGExport));
        add( testCase(&ImageExportImportFailureTest::testJPEGImport));
        add( testCase(&ImageExportImportFailureTest::testTIFFExport));
        add( testCase(&ImageExportImportFailureTest::testTIFFImport));
        add( testCase(&ImageExportImportFailureTest::testBMPExport));
        add( testCase(&ImageExportImportFailureTest::testBMPImport));
        add( testCase(&ImageExportImportFailureTest::testSUNExport));
        add( testCase(&ImageExportImportFailureTest::testSUNImport));
        add( testCase(&ImageExportImportFailureTest::testVIFFExport));
        add( testCase(&ImageExportImportFailureTest::testVIFFImport));
    }
};

int main()
{
    ImageImportExportTestSuite test;
    
    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
