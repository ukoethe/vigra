#include <iostream>
#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/impex.hxx"
#include "unittest.hxx"

using namespace vigra;

class ByteImageExportImportTest
{
    typedef vigra::BImage Image;
  public:
    ByteImageExportImportTest()
    {
        vigra::ImageImportInfo info("lenna.xv");

        int w = info.width();
        int h = info.height();

        img.resize(w, h);

        importImage(info, destImage(img));
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
        
        vigra::ImageImportInfo info("res.foo");

        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        should(info.getFileType() == std::string("VIFF"));
        
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
    
    Image img;
};

void ByteImageExportImportTest::testFile(const char * filename)
{
    exportImage(srcImageRange(img), vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == img.width());
    should(info.height() == img.height());
    should(info.isGrayscale());
#ifndef NEWIMPEX
    should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
    should(info.pixelType() == std::string("UINT8"));
#endif

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
        vigra::ImageImportInfo info("lennargb.xv");

        int w = info.width();
        int h = info.height();

        img.resize(w, h);

        importImage(info, destImage(img));
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
        
        vigra::ImageImportInfo info("res.foo");

        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        should(info.getFileType() == std::string("VIFF"));
        
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
    
    Image img;
};

void ByteRGBImageExportImportTest::testFile(const char * fileName)
{
    exportImage(srcImageRange(img), vigra::ImageExportInfo(fileName));

    vigra::ImageImportInfo info(fileName);

    should(info.width() == img.width());
    should(info.height() == img.height());
    should(info.isColor());
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif

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
        vigra::ImageImportInfo info("lenna.xv");

        int w = info.width();
        int h = info.height();

        img.resize(w, h);

        importImage(info, destImage(img));
        
        vigra::ImageImportInfo rinfo("lennafloat.xv");

        reread.resize(w, h);

        importImage(rinfo, destImage(reread));
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.gif"));
        
        vigra::ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isGrayscale());
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
#else
        should(info.pixelType() == std::string("FLOAT"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
#else
        should(info.pixelType() == std::string("FLOAT"));
#endif
        
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
        vigra::ImageImportInfo info("lennargb.xv");

        int w = info.width();
        int h = info.height();

        img.resize(w, h);

        importImage(info, destImage(img));
        
        vigra::ImageImportInfo rinfo("lennafloatrgb.xv");

        reread.resize(w, h);

        importImage(rinfo, destImage(reread));
    }
    
    void testGIF()
    {
        exportImage(srcImageRange(img), vigra::ImageExportInfo("res.gif"));
        
        vigra::ImageImportInfo info("res.gif");
        
        should(info.width() == img.width());
        should(info.height() == img.height());
        should(info.isColor());
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
#else
        should(info.pixelType() == std::string("FLOAT"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::UINT8);
#else
        should(info.pixelType() == std::string("UINT8"));
#endif
        
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
#ifndef NEWIMPEX
        should(info.pixelType() == vigra::ImageImportInfo::FLOAT);
#else
        should(info.pixelType() == std::string("FLOAT"));
#endif
        
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
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.gif"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.gif'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testGIFImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.gif");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.gif'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testJPEGExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.jpg"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.jpg'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testJPEGImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.jpg");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.jpg'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testTIFFExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.tif"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.tif'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testTIFFImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.tif");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.tif'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testBMPExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.bmp"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.bmp'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testBMPImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.bmp");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.bmp'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testPNMExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.pnm"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.pnm'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testPNMImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.pnm");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.pnm'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testSUNExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.ras"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.ras'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testSUNImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.ras");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.ras'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
        
    void testVIFFExport()
    {
        try
        {
            exportImage(srcImageRange(img), vigra::ImageExportInfo("intentionalFailure/res.xv"));
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'intentionalFailure/res.xv'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
    }
    
    void testVIFFImport()
    {
        try
        {
            vigra::ImageImportInfo info("foo.xv");
            failTest("Failed to throw exception.");
        }
        catch(vigra::PreconditionViolation & e)
        {
            char const * expected = "\nPrecondition violation!\n"
                                    "Unable to open file 'foo.xv'.";
            should(strncmp(expected, e.what(), strlen(expected)) == 0);
        }
        catch(vigra::PostconditionViolation & e)  // for old impex
        {}
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
#ifndef NEWIMPEX
        add( testCase(&ByteImageExportImportTest::testGIF));
        add( testCase(&ByteImageExportImportTest::testEmptyGIF));
#endif
        add( testCase(&ByteImageExportImportTest::testJPEG));
        add( testCase(&ByteImageExportImportTest::testTIFF));
        add( testCase(&ByteImageExportImportTest::testBMP));
        add( testCase(&ByteImageExportImportTest::testPGM));
        add( testCase(&ByteImageExportImportTest::testPNM));
        add( testCase(&ByteImageExportImportTest::testSUN));
        add( testCase(&ByteImageExportImportTest::testVIFF1));
        add( testCase(&ByteImageExportImportTest::testVIFF2));

#ifndef NEWIMPEX
        add( testCase(&ByteRGBImageExportImportTest::testGIF));
#endif
        add( testCase(&ByteRGBImageExportImportTest::testJPEG));
        add( testCase(&ByteRGBImageExportImportTest::testTIFF));
        add( testCase(&ByteRGBImageExportImportTest::testBMP));
        add( testCase(&ByteRGBImageExportImportTest::testPPM));
        add( testCase(&ByteRGBImageExportImportTest::testPNM));
        add( testCase(&ByteRGBImageExportImportTest::testSUN));
        add( testCase(&ByteRGBImageExportImportTest::testVIFF1));
        add( testCase(&ByteRGBImageExportImportTest::testVIFF2));

#ifndef NEWIMPEX
        add( testCase(&FloatImageExportImportTest::testGIF));
#endif
        add( testCase(&FloatImageExportImportTest::testJPEG));
        add( testCase(&FloatImageExportImportTest::testTIFF));
        add( testCase(&FloatImageExportImportTest::testBMP));
        add( testCase(&FloatImageExportImportTest::testSUN));
        add( testCase(&FloatImageExportImportTest::testVIFF));

#ifndef NEWIMPEX
        add( testCase(&FloatRGBImageExportImportTest::testGIF));
#endif
        add( testCase(&FloatRGBImageExportImportTest::testJPEG));
        add( testCase(&FloatRGBImageExportImportTest::testTIFF));
        add( testCase(&FloatRGBImageExportImportTest::testBMP));
        add( testCase(&FloatRGBImageExportImportTest::testSUN));
        add( testCase(&FloatRGBImageExportImportTest::testVIFF));

#ifndef NEWIMPEX
        add( testCase(&ImageExportImportFailureTest::testGIFExport));
        add( testCase(&ImageExportImportFailureTest::testGIFImport));
#endif
        add( testCase(&ImageExportImportFailureTest::testJPEGExport));
        add( testCase(&ImageExportImportFailureTest::testJPEGImport));
        add( testCase(&ImageExportImportFailureTest::testTIFFExport));
        add( testCase(&ImageExportImportFailureTest::testTIFFImport));
        add( testCase(&ImageExportImportFailureTest::testBMPExport));
        add( testCase(&ImageExportImportFailureTest::testBMPImport));
        add( testCase(&ImageExportImportFailureTest::testPNMExport));
        add( testCase(&ImageExportImportFailureTest::testPNMImport));
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
