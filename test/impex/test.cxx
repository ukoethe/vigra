#include <iostream>
#include <cmath>
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "unittest.hxx"

using namespace vigra;

class ByteImageExportImportTest
{
    typedef vigra::BImage Image;

public:

    ByteImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lenna.xv");

        const int w = info.width ();
        const int h = info.height ();

        img.resize (w, h);

        importImage (info, destImage (img));
    }

    void testListFormatsExtensions()
    {
        const std::string formats = impexListFormats();
        const std::string extensions = impexListExtensions();

        shouldEqual(formats, "BMP GIF JPEG PNG PNM SUN TIFF VIFF");
        shouldEqual(extensions, "bmp gif jpeg jpg pbm pgm png pnm ppm ras tif tiff xv");
    }

    void testIsImage()
    {
        should(isImage("lenna.xv"));
        should(!isImage("Makefile"));
    }

    void testFile (const char *filename);

    void testGIF ()
    {
        vigra::ImageExportInfo exportinfo ("res.gif");
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.gif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        float sum = 0;
        for (; i != img.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));

        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testJPEG ()
    {
        vigra::ImageExportInfo exportinfo ("res.jpg");
        exportinfo.setCompression ("100");
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        float sum = 0;
        for (; i != img.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));

        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testTIFF ()
    {
        vigra::ImageExportInfo exportinfo ("res.tif");
        exportinfo.setCompression ("RunLength");
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    void testBMP ()
    {
        testFile ("res.bmp");
    }

    void testPGM ()
    {
        testFile ("res.pgm");
    }

    void testPNM ()
    {
        testFile ("res.pnm");
    }

    void testPNM2 ()
    {
        vigra::ImageExportInfo exportinfo ("res.pgm");
        exportinfo.setCompression ("ASCII");
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.pgm");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("PNM"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    void testPNG ()
    {
        testFile ("res.png");
    }

    void testSUN ()
    {
        testFile ("res.ras");
    }

    void testVIFF1 ()
    {
        testFile ("res.xv");
    }

    void testVIFF2 ()
    {
        vigra::ImageExportInfo exportinfo ("res.foo");
        exportinfo.setFileType ("VIFF");
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.foo");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("VIFF"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    Image img;
};

void
ByteImageExportImportTest::testFile (const char *filename)
{
    exportImage (srcImageRange (img), vigra::ImageExportInfo (filename));

    vigra::ImageImportInfo info (filename);

    should (info.width () == img.width ());
    should (info.height () == img.height ());
    should (info.isGrayscale ());
    should (info.getPixelType () == std::string ("UINT8"));

    Image res (info.width (), info.height ());

    importImage (info, destImage (res));

    Image::ScanOrderIterator i = img.begin ();
    Image::ScanOrderIterator i1 = res.begin ();
    Image::Accessor acc = img.accessor ();

    for (; i != img.end (); ++i, ++i1)
        {
            if (acc (i) != acc (i1))
                {
                    std::cerr << acc (i) - acc (i1) << std::endl;
                    exit (1);
                }
            should (acc (i) == acc (i1));
        }
}


class ByteRGBImageExportImportTest
{
    typedef vigra::BRGBImage Image;
public:
    ByteRGBImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lennargb.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize (w, h);

        importImage (info, destImage (img));
    }

    void testFile (const char *fileName);

    void testGIF ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("resrgb.gif"));

        vigra::ImageImportInfo info ("resrgb.gif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.width (), info.height ());
        Image ref (info.width (), info.height ());

        importImage (info, destImage (res));
        importImage (vigra::ImageImportInfo("lenna_gifref.xv"), destImage (ref));

        Image::ScanOrderIterator i = ref.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = ref.accessor ();

        float sum = 0;
        for (; i != ref.end (); ++i, ++i1)
                sum += (acc (i) - acc (i1)).magnitude ();

        should (sum / (info.width () * info.height ()) < 3.0);  // use rather large tolerance to make the 
                                                                // test portable
    }

    void testJPEG ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        float sum = 0;
        for (; i != img.end (); ++i, ++i1)
            {
                sum += (acc (i) - acc (i1)).magnitude ();
            }
        should (sum / (info.width () * info.height ()) < 2.0);
    }

    void testTIFF ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").
                     setCompression ("RunLength"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            {
                should (acc (i) == acc (i1));
            }
    }

    void testBMP ()
    {
        testFile ("res.bmp");
    }

    void testPPM ()
    {
        testFile ("res.ppm");
    }

    void testPNM ()
    {
        testFile ("res.pnm");
    }


    void testPNM2 ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.ppm").setCompression ("ASCII"));

        vigra::ImageImportInfo info ("res.ppm");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("PNM"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    void testPNG ()
    {
        testFile ("res.png");
    }

    void testSUN ()
    {
        testFile ("res.ras");
    }

    void testVIFF1 ()
    {
        testFile ("res.xv");
    }

    void testVIFF2 ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.foo").setFileType ("VIFF"));

        vigra::ImageImportInfo info ("res.foo");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("VIFF"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    Image img;
};

void
ByteRGBImageExportImportTest::testFile (const char *fileName)
{
    exportImage (srcImageRange (img), vigra::ImageExportInfo (fileName));

    vigra::ImageImportInfo info (fileName);

    should (info.width () == img.width ());
    should (info.height () == img.height ());
    should (info.isColor ());
    should (info.getPixelType () == std::string ("UINT8"));

    Image res (info.width (), info.height ());

    importImage (info, destImage (res));

    Image::ScanOrderIterator i = img.begin ();
    Image::ScanOrderIterator i1 = res.begin ();
    Image::Accessor acc = img.accessor ();

    for (; i != img.end (); ++i, ++i1)
        should (acc (i) == acc (i1));
}

class PNGInt16Test
{
  public:
    void testByteOrder()
    {
        SImage i(1,1);
        i(0,0) = 1;
        exportImage(srcImageRange(i), ImageExportInfo("res.png"));
        ImageImportInfo info("res.png");
        shouldEqual(info.width(), 1);
        shouldEqual(info.height(), 1);
        shouldEqual(info.numBands(), 1);
        shouldEqual(info.isGrayscale(), true);
        shouldEqual(std::string(info.getPixelType()), std::string("INT16"));
        i(0,0) = 0;
        importImage(info, destImage(i));
        shouldEqual(i(0,0), 1);

        BasicImage<RGBValue<short> > rgb(1,1);
        rgb(0,0) = RGBValue<short>(1,2,3);
        exportImage(srcImageRange(rgb), ImageExportInfo("res.png"));
        ImageImportInfo rgbinfo("res.png");
        shouldEqual(rgbinfo.width(), 1);
        shouldEqual(rgbinfo.height(), 1);
        shouldEqual(rgbinfo.numBands(), 3);
        shouldEqual(std::string(rgbinfo.getPixelType()), std::string("INT16"));
        rgb(0,0) = RGBValue<short>(0,0,0);
        importImage(rgbinfo, destImage(rgb));
        shouldEqual(rgb(0,0), RGBValue<short>(1,2,3));
    }
};

class FloatImageExportImportTest
{
    typedef vigra::FImage Image;

public:

    FloatImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lenna.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize (w, h);

        importImage (info, destImage (img));

        vigra::ImageImportInfo rinfo ("lennafloat.xv");

        reread.resize (w, h);

        importImage (rinfo, destImage (reread));
    }

    void testGIF()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.gif"));

        vigra::ImageImportInfo info ("res.gif");

        should (info.width () == reread.width ());
        should (info.height () == reread.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        float sum = 0;
        for (; i != reread.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));
        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testJPEG ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == reread.width ());
        should (info.height () == reread.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        float sum = 0;
        for (; i != reread.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));
        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testTIFF ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();
        
        shouldEqualSequence(i, img.end(), i1);

 /*       for (; i != img.end (); ++i, ++i1)
            {
                shouldEqual(acc (i), acc (i1));
            }
 */
    }

    void testBMP ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.bmp"));

        vigra::ImageImportInfo info ("res.bmp");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        for (; i != reread.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }

    void testSUN ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.ras"));

        vigra::ImageImportInfo info ("res.ras");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        for (; i != reread.end (); ++i, ++i1)
            {
                should (acc (i) == acc (i1));
            }
    }

    void testVIFF ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            {
                should (acc (i) == acc (i1));
            }
    }

    Image img, reread;
};

class FloatRGBImageExportImportTest
{
    typedef vigra::FRGBImage Image;
    Image img, reread;

public:

    FloatRGBImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lennargb.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize (w, h);

        importImage (info, destImage (img));

        vigra::ImageImportInfo rinfo ("lennafloatrgb.xv");

        reread.resize (w, h);

        importImage (rinfo, destImage (reread));
    }

    void testJPEG ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == reread.width ());
        should (info.height () == reread.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += (acc (i) - acc (i1)).magnitude ();
        should (sum / (info.width () * info.height ()) < 2.0f);
    }

    void testTIFF ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (acc (i), acc (i1));
    }

    void testBMP ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.bmp"));

        vigra::ImageImportInfo info ("res.bmp");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += (acc (i) - acc (i1)).magnitude ();
        should (sum / (info.width () * info.height ()) < 2.0f);
    }

    void testSUN ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.ras"));

        vigra::ImageImportInfo info ("res.ras");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += (acc (i) - acc (i1)).magnitude ();
        should (sum / (info.width () * info.height ()) < 2.0f);
    }

    void testVIFF ()
    {
        exportImage (srcImageRange (img), vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            should (acc (i) == acc (i1));
    }
};

class Vector4ExportImportTest
{
  public:
  
    typedef vigra::FVector4Image Image;
    typedef vigra::BasicImage<TinyVector<unsigned char, 4> > BImage;
    Image img, reread;
    BImage breread, breference;

public:

    Vector4ExportImportTest ()
    : img(2,3),
      reread(2,3),
      breread(2,3),
      breference(2,3)
    {
        double scale = 255.0 / 11.0;
        double offset = 5.5;
        for(int y = 0; y<3; ++y)
        {
            for(int x=0; x<2; ++x)
            {
                img(x,y)[0] = 2*y+x + 0.5;
                img(x,y)[1] = -img(x,y)[0];
                img(x,y)[2] = 0.0;
                img(x,y)[3] = 0.5;
                for(int b=0; b<4; ++b)
                {
                    breference(x,y)[b] = 
                        NumericTraits<unsigned char>::fromRealPromote(scale*(img(x,y)[b]+offset));
                }
            }
        }
    }

    void failingTest (char const * filename)
    {
        try 
        {
            exportImage( srcImageRange(img), vigra::ImageExportInfo( filename ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::PreconditionViolation & e ) 
        {
            std::string expected = "\nPrecondition violation!\n";
            expected += "exportImage(): file format does not support requested number of bands (color channels)";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testJPEG ()
    {
        failingTest("res.jpg");
    }

    void testGIF ()
    {
        failingTest("res.gif");
    }

    void testBMP ()
    {
        failingTest("res.bmp");
    }

    void testPNM ()
    {
        failingTest("res.pnm");
    }

    void testSUN ()
    {
        failingTest("res.ras");
    }

    void testVIFF ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("FLOAT"));

        importImage (info, destImage (reread));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = reread.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (acc (i), acc (i1));
    }

    void testTIFF ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("FLOAT"));

        importImage (info, destImage (reread));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = reread.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (acc (i), acc (i1));
    }

    void testPNG ()
    {
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.png"));

        vigra::ImageImportInfo info ("res.png");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("UINT8"));

        importImage (info, destImage (breread));

        BImage::ScanOrderIterator i = breference.begin ();
        BImage::ScanOrderIterator i1 = breread.begin ();
        BImage::Accessor acc = breference.accessor ();

        for (; i != breference.end (); ++i, ++i1)
        {
            should ((acc (i)- acc (i1)).magnitude() <= 1.0);
        }
    }
};


class ImageExportImportFailureTest
{
    vigra::BImage img;

public:

    ImageExportImportFailureTest()
        : img( 3, 3 )
    {}

    // gif

    void testGIFExport()
    {
        testExport("gif");
    }

    void testGIFImport()
    {
        testImport("gif");
    }

    // jpeg

    void testJPEGExport()
    {
        testExport("jpg");
    }

    void testJPEGImport()
    {
        testImport("jpg");
    }

    // tiff

    void testTIFFExport()
    {
        testExport("tiff");
    }

    void testTIFFImport()
    {
        testImport("tiff");
    }

    // viff

    void testVIFFExport()
    {
        testExport("xv");
    }

    void testVIFFImport()
    {
        testImport("xv");
    }

    // sun

    void testSUNExport()
    {
        testExport("ras");
    }

    void testSUNImport()
    {
        testImport("ras");
    }

    // pnm

    void testPNMExport()
    {
        testExport("pnm");
    }

    void testPNMImport()
    {
        testImport("pnm");
    }

    // png

    void testPNGExport()
    {
        testExport("png");
    }

    void testPNGImport()
    {
        testImport("png");
    }

    // bmp

    void testBMPExport()
    {
        testExport("bmp");
    }

    void testBMPImport()
    {
        testImport("bmp");
    }

    // test implementation

    void testImport( const char * fext )
    {
        std::string fname = "foo.";
        fname += fext;
        try {
            vigra::ImageImportInfo info( fname.c_str() );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::PreconditionViolation & e ) {
            std::string expected = "\nPrecondition violation!\n";
            expected += "Unable to open file '";
            expected += fname;
            expected += "'.";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testExport( const char * fext )
    {
        std::string fname = "intentionalFailure/foo.";
        fname += fext;
        try {
            exportImage( srcImageRange(img), vigra::ImageExportInfo( fname.c_str() ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::PreconditionViolation & e ) {
            std::string expected = "\nPrecondition violation!\n";
            expected += "Unable to open file '";
            expected += fname;
            expected += "'.";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }
};

struct ImageImportExportTestSuite : public vigra::test_suite
{
    ImageImportExportTestSuite()
        : vigra::test_suite("ImageImportExportTestSuite")
    {
        // general tests
        add(testCase(&ByteImageExportImportTest::testListFormatsExtensions));
        add(testCase(&ByteImageExportImportTest::testIsImage));
        
        // grayscale byte images
        add(testCase(&ByteImageExportImportTest::testGIF));
        add(testCase(&ByteImageExportImportTest::testJPEG));
        add(testCase(&ByteImageExportImportTest::testTIFF));
        add(testCase(&ByteImageExportImportTest::testBMP));
        add(testCase(&ByteImageExportImportTest::testPGM));
        add(testCase(&ByteImageExportImportTest::testPNM));
        add(testCase(&ByteImageExportImportTest::testPNM2));
        add(testCase(&ByteImageExportImportTest::testPNG));
        add(testCase(&ByteImageExportImportTest::testSUN));
        add(testCase(&ByteImageExportImportTest::testVIFF1));
        add(testCase(&ByteImageExportImportTest::testVIFF2));

        // rgb byte images
        add(testCase(&ByteRGBImageExportImportTest::testGIF));
        add(testCase(&ByteRGBImageExportImportTest::testJPEG));
        add(testCase(&ByteRGBImageExportImportTest::testTIFF));
        add(testCase(&ByteRGBImageExportImportTest::testBMP));
        add(testCase(&ByteRGBImageExportImportTest::testPPM));
        add(testCase(&ByteRGBImageExportImportTest::testPNM));
        add(testCase(&ByteRGBImageExportImportTest::testPNM2));
        add(testCase(&ByteRGBImageExportImportTest::testPNG));
        add(testCase(&ByteRGBImageExportImportTest::testSUN));
        add(testCase(&ByteRGBImageExportImportTest::testVIFF1));
        add(testCase(&ByteRGBImageExportImportTest::testVIFF2));

        // 16-bit PNG
        add(testCase(&PNGInt16Test::testByteOrder));

        // grayscale float images
        add(testCase(&FloatImageExportImportTest::testGIF));
        add(testCase(&FloatImageExportImportTest::testJPEG));
        add(testCase(&FloatImageExportImportTest::testTIFF));
        add(testCase(&FloatImageExportImportTest::testBMP));
        add(testCase(&FloatImageExportImportTest::testSUN));
        add(testCase(&FloatImageExportImportTest::testVIFF));

        // rgb float images
        add(testCase(&FloatRGBImageExportImportTest::testJPEG));
        add(testCase(&FloatRGBImageExportImportTest::testTIFF));
        add(testCase(&FloatRGBImageExportImportTest::testBMP));
        add(testCase(&FloatRGBImageExportImportTest::testSUN));
        add(testCase(&FloatRGBImageExportImportTest::testVIFF));

        // 4-band images
        add(testCase(&Vector4ExportImportTest::testJPEG));
        add(testCase(&Vector4ExportImportTest::testGIF));
        add(testCase(&Vector4ExportImportTest::testBMP));
        add(testCase(&Vector4ExportImportTest::testPNM));
        add(testCase(&Vector4ExportImportTest::testSUN));
        add(testCase(&Vector4ExportImportTest::testVIFF));
        add(testCase(&Vector4ExportImportTest::testTIFF));
        add(testCase(&Vector4ExportImportTest::testPNG));

        // failure tests
        add(testCase(&ImageExportImportFailureTest::testGIFExport));
        add(testCase(&ImageExportImportFailureTest::testGIFImport));
        add(testCase(&ImageExportImportFailureTest::testJPEGExport));
        add(testCase(&ImageExportImportFailureTest::testJPEGImport));
        add(testCase(&ImageExportImportFailureTest::testTIFFExport));
        add(testCase(&ImageExportImportFailureTest::testTIFFImport));
        add(testCase(&ImageExportImportFailureTest::testBMPExport));
        add(testCase(&ImageExportImportFailureTest::testBMPImport));
        add(testCase(&ImageExportImportFailureTest::testPNMExport));
        add(testCase(&ImageExportImportFailureTest::testPNMImport));
        add(testCase(&ImageExportImportFailureTest::testPNGExport));
        add(testCase(&ImageExportImportFailureTest::testPNGImport));
        add(testCase(&ImageExportImportFailureTest::testSUNExport));
        add(testCase(&ImageExportImportFailureTest::testSUNImport));
        add(testCase(&ImageExportImportFailureTest::testVIFFExport));
        add(testCase(&ImageExportImportFailureTest::testVIFFImport));
    }
};

int main ()
{
    ImageImportExportTestSuite test;
    const int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed != 0;
}
