#include <iostream>
#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/impex.hxx"
#include "unittest.h"

using namespace vigra;

class ByteImageExportImportTest
{
    typedef vigra::BImage Image;
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
    
    void testFile(const char * filename);
    
    void testGIF()
    {
        testFile("res.gif");
    }

    void testEmptyGIF()
    {
        Image img(100,100);
        
        img = 0;
        
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.gif"));
        
        vigra::ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
                    vigra::ImageExportInfo("res.jpg").setCompression("100"));
        
        vigra::ImageImportInfo info("res.jpg");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = img.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = img.accessor();

	float sum = 0;
        for(; i != img.end(); ++i, ++i1)
	{
            sum += std::abs(acc(i) - acc(i1));
	}
        should(sum / (info.width() * info.height()) < 0.1);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    vigra::ImageExportInfo("res.tif").setCompression("RunLength"));
        
        vigra::ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        testFile("res.bmp");
    }
    
    void testPGM()
    {
        testFile("res.pgm");
    }
    
    void testPNM()
    {
        testFile("res.pnm");
    }
    
    void testSUN()
    {
        testFile("res.ras");
    }
    
    void testVIFF1()
    {
        testFile("res.xv");
    }
    
    void testVIFF2()
    {
        exportImage(srcImageRange(img), 
                    vigra::ImageExportInfo("res.foo").setFileType("VIFF"));
        
        ViffImage * viff = readViffImage("res.foo");
        
        should(viff != 0);
        should((int)viff->row_size == img.width());
        should((int)viff->col_size == img.height());
        
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

void ByteImageExportImportTest::testFile(const char * filename)
{
    exportImage(srcImageRange(img), vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == img.width());
    should(info.height() == img.height());
    should(info.isGrayscale());
    should(info.pixelType() == vigra::ImageImportInfo::UINT8);

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


class ByteRGBImageExportImportTest
{
    typedef vigra::BRGBImage Image;
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
    
    void testFile(const char * fileName);
    
    void testGIF()
    {
        testFile("res.gif");
    }
    
    void testJPEG()
    {
        exportImage(srcImageRange(img), 
                    vigra::ImageExportInfo("res.jpg").setCompression("100"));
        
        vigra::ImageImportInfo info("res.jpg");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
                    vigra::ImageExportInfo("res.tif").setCompression("RunLength"));
        
        vigra::ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        testFile("res.bmp");
    }
    
    void testPPM()
    {
        testFile("res.ppm");
    }
    
    void testPNM()
    {
        testFile("res.pnm");
    }
    
    void testSUN()
    {
        testFile("res.ras");
    }
    
    void testVIFF1()
    {
        testFile("res.xv");
    }
    
    void testVIFF2()
    {
        exportImage(srcImageRange(img), 
                    vigra::ImageExportInfo("res.foo").setFileType("VIFF"));
        
        ViffImage * viff = readViffImage("res.foo");
        
        should(viff != 0);
        should((int)viff->row_size == img.width());
        should((int)viff->col_size == img.height());
        
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

void ByteRGBImageExportImportTest::testFile(const char * fileName)
{
    exportImage(srcImageRange(img), vigra::ImageExportInfo(fileName));

    vigra::ImageImportInfo info(fileName);

    should(info.width() == img.width());
    should(info.height() == img.height());
    should(info.isColor());
    should(info.pixelType() == vigra::ImageImportInfo::UINT8);

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
    
class FloatImageExportImportTest
{
    typedef vigra::FImage Image;
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.gif"));
        
        vigra::ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
                    vigra::ImageExportInfo("res.jpg").setCompression("100"));
        
        vigra::ImageImportInfo info("res.jpg");
        
        should(info.width() == reread.width());
        should(info.height() == reread.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
        Image res(info.width(), info.height());
        
        importImage(info, destImage(res));
        
	Image::ScanOrderIterator i = reread.begin();
	Image::ScanOrderIterator i1 = res.begin();
	Image::Accessor acc = reread.accessor();

	float sum = 0;
        for(; i != reread.end(); ++i, ++i1)
	{
            sum += std::abs(acc(i) - acc(i1));
	}
        should(sum / (info.width() * info.height()) < 0.1);
    }

    void testTIFF()
    {
        exportImage(srcImageRange(img), 
                    vigra::ImageExportInfo("res.tif").setCompression("LZW"));
        
        vigra::ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.bmp"));
        
        vigra::ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.ras"));
        
        vigra::ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.xv"));
        
        vigra::ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
        
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
    typedef vigra::FRGBImage Image;
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.gif"));
        
        vigra::ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
                    vigra::ImageExportInfo("res.jpg").setCompression("100"));
        
        vigra::ImageImportInfo info("res.jpg");
        
        should(info.width() == reread.width());
        should(info.height() == reread.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
                    vigra::ImageExportInfo("res.tif").setCompression("LZW"));
        
        vigra::ImageImportInfo info("res.tif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.bmp"));
        
        vigra::ImageImportInfo info("res.bmp");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.ras"));
        
        vigra::ImageImportInfo info("res.ras");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
        
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
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.xv"));
        
        vigra::ImageImportInfo info("res.xv");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
        
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.gif"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.gif");
        }
        catch(vigra::PostconditionViolation &)
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.jpg"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.jpg");
        }
        catch(vigra::PostconditionViolation &)
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.tif"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.tif");
        }
        catch(vigra::PostconditionViolation &)
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.bmp"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.bmp");
        }
        catch(vigra::PostconditionViolation &)
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.ras"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.ras");
        }
        catch(vigra::PostconditionViolation &)
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
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.xv"));
        }
        catch(vigra::PostconditionViolation &)
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
            vigra::ImageImportInfo info("foo.xv");
        }
        catch(vigra::PostconditionViolation &)
        {
            caught = true;
        }
        should(caught == true);
    }
        
    vigra::BImage img;
};

class ImageImportExportTestSuite
: public vigra::test_suite
{
  public:
    ImageImportExportTestSuite()
    : vigra::test_suite("ImageImportExportTestSuite")
    {
        add( testCase(&ByteImageExportImportTest::testGIF));
        add( testCase(&ByteImageExportImportTest::testEmptyGIF));
        add( testCase(&ByteImageExportImportTest::testJPEG));
        add( testCase(&ByteImageExportImportTest::testTIFF));
        add( testCase(&ByteImageExportImportTest::testBMP));
        add( testCase(&ByteImageExportImportTest::testPGM));
        add( testCase(&ByteImageExportImportTest::testPNM));
        add( testCase(&ByteImageExportImportTest::testSUN));
        add( testCase(&ByteImageExportImportTest::testVIFF1));
        add( testCase(&ByteImageExportImportTest::testVIFF2));

        add( testCase(&ByteRGBImageExportImportTest::testGIF));
        add( testCase(&ByteRGBImageExportImportTest::testJPEG));
        add( testCase(&ByteRGBImageExportImportTest::testTIFF));
        add( testCase(&ByteRGBImageExportImportTest::testBMP));
        add( testCase(&ByteRGBImageExportImportTest::testPPM));
        add( testCase(&ByteRGBImageExportImportTest::testPNM));
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
