/************************************************************************/
/*                                                                      */
/*           Copyright 2004-2012 by Ullrich Koethe                      */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/impexalpha.hxx"
#include "unittest.hxx"
#include "vigra/multi_array.hxx"

using namespace vigra;

template <class Image>
void failCodec(Image const & img, ImageExportInfo const & info)
{
    try {
        exportImage (srcImageRange (img), info);
        failTest( "Failed to throw exception." );
    }
    catch( vigra::PreconditionViolation & e )
    {
        std::string expected = "\nPrecondition violation!\n";
        expected += "did not find a matching codec for the given file extension";
        const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
        should(rc);
    }
}

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

        const char * refFormats = "BMP "
#if defined(HasEXR)
        "EXR "
#endif
        "GIF HDR "
#if defined(HasJPEG)
        "JPEG "
#endif
#if defined(HasPNG)
        "PNG "
#endif
        "PNM SUN "
#if defined(HasTIFF)
        "TIFF "
#endif
        "VIFF";
        shouldEqual(formats, refFormats);

        const char * refExtensions = "bmp "
#if defined(HasEXR)
        "exr "
#endif
        "gif hdr "
#if defined(HasJPEG)
        "jpeg jpg "
#endif
        "pbm pgm "
#if defined(HasPNG)
        "png "
#endif
        "pnm ppm ras "
#if defined(HasTIFF)
        "tif tiff "
#endif
        "xv";
        shouldEqual(extensions, refExtensions);
    }

    void testIsImage()
    {
        should(isImage("lenna.xv"));
        should(!isImage("no-image.txt"));
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
#if !defined(HasJPEG)
        failCodec(img, exportinfo);
#else
        exportinfo.setCompression ("JPEG QUALITY=100");
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
#endif
    }

    void testTIFF ()
    {
        vigra::ImageExportInfo exportinfo ("res.tif");
#if !defined(HasTIFF)
        failCodec(img, exportinfo);
#else
        exportinfo.setCompression ("LZW");
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
#endif
    }

    void testTIFFSequence()
    {
#if defined(HasTIFF)
        for (int i=0; i < 3; ++i)
        {
            std::string fileName = std::string("lenna_") + vigra::asString(i) + ".tif";
            vigra::ImageImportInfo ininfo (fileName.c_str());
            Image inimg(ininfo.width(), ininfo.height());
            importImage(ininfo, destImage(inimg));
            vigra::ImageExportInfo outinfo ("resseq.tif", i==0?"w":"a");
            exportImage(srcImageRange(inimg), outinfo);
        }

        for (int j=0; j < 3; ++j)
        {
            vigra::ImageImportInfo ininfo ("resseq.tif", j);
            Image inimg(ininfo.width(), ininfo.height());
            std::string fileName = std::string("lenna_") + vigra::asString(j) + ".tif";
            importImage(ininfo, destImage(inimg));
            vigra::ImageImportInfo originfo (fileName.c_str());
            Image origimg(originfo.width(), originfo.height());
            importImage(originfo, destImage(origimg));

            Image::ScanOrderIterator it = inimg.begin ();
            Image::ScanOrderIterator it1 = origimg.begin ();
            Image::Accessor acc = inimg.accessor ();
            for (; it != inimg.end (); ++it, ++it1)
                should (acc (it) == acc (it1));
        }
#endif
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
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        testFile ("res.png");
#endif
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
        should (acc (i) == acc (i1));
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

        double sum = 0.0;
        for (; i != ref.end (); ++i, ++i1)
                sum += (acc (i) - acc (i1)).magnitude ();

        should (sum / (info.width () * info.height ()) < 4.0);  // use rather large tolerance to make the
                                                                // test portable
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

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

        double sum = 0.0;
        for (; i != img.end (); ++i, ++i1)
            {
                sum += (acc (i) - acc (i1)).magnitude ();
            }
        should (sum / (info.width () * info.height ()) < 2.0);
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").
                     setCompression ("LZW"));

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
#endif
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
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        testFile ("res.png");
#endif
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

class CanvasSizeTest
{
  public:
    void testTIFFCanvasSize ()
    {
        vigra::ImageExportInfo exportinfo ("res.tif");
        FRGBImage img(1, 1);
#if !defined(HasTIFF)
        failCodec(img, exportinfo);
#else
        img(0,0) = 1;
        exportinfo.setCompression ("LZW");
        Size2D canvasSize(3, 8);
        exportinfo.setCanvasSize (canvasSize);
        exportImage (srcImageRange (img), exportinfo);

        vigra::ImageImportInfo info ("res.tif");

        should (info.getCanvasSize () == canvasSize);
#endif
    }
};

class PNGInt16Test
{
  public:
    void testByteOrder()
    {
        UInt16Image i(1,1);
        i(0,0) = 1;
        exportImage(srcImageRange(i), ImageExportInfo("res.png"));
        ImageImportInfo info("res.png");
        shouldEqual(info.width(), 1);
        shouldEqual(info.height(), 1);
        shouldEqual(info.numBands(), 1);
        shouldEqual(info.isGrayscale(), true);
        shouldEqual(std::string(info.getPixelType()), std::string("UINT16"));
        i(0,0) = 0;
        importImage(info, destImage(i));
        shouldEqual(i(0,0), 1);

        // DGSW: Note that this produces a PNG standard conformant image
        //       but both Imagemagick 'identify' and photoshop CS2 see
        //       the data incorrectly
        BasicImage<RGBValue<unsigned short> > rgb(1,1);
        // Using unsigned values 0xff01, 0xfff1, 0xfffd
        rgb(0,0) = RGBValue<unsigned short>(65281,65521,65533);
        // Using unsigned values 0x7f01, 0x7ff1, 0x7ffd
        // rgb(0,0) = RGBValue<unsigned short>(32513,32753,32765);
        exportImage(srcImageRange(rgb), ImageExportInfo("res.png"));
        ImageImportInfo rgbinfo("res.png");
        shouldEqual(rgbinfo.width(), 1);
        shouldEqual(rgbinfo.height(), 1);
        shouldEqual(rgbinfo.numBands(), 3);
        shouldEqual(std::string(rgbinfo.getPixelType()), std::string("UINT16"));
        rgb(0,0) = RGBValue<short>(0,0,0);
        importImage(rgbinfo, destImage(rgb));
        shouldEqual(rgb(0,0), RGBValue<unsigned short>(65281,65521,65533));
//        shouldEqual(rgb(0,0), RGBValue<unsigned short>(32513,32753,32765));
    }
};

class FloatImageExportImportTest
{
    typedef vigra::DImage Image;
    std::string rereadType;

public:

    FloatImageExportImportTest ()
    : rereadType("DOUBLE")
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

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));
        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

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

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));
        should (sum / (info.width () * info.height ()) < 0.1);
#endif
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo ("res.png"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res_.png"));

        vigra::ImageImportInfo info ("res_.png");

        should (info.width () == reread.width ());
        should (info.height () == reread.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = reread.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = reread.accessor ();

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += std::abs (acc (i) - acc (i1));
        should (sum / (info.width () * info.height ()) < 0.1);
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == rereadType);

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        //Image::Accessor acc = img.accessor ();

        shouldEqualSequence(i, img.end(), i1);
#endif
    }

    void testTIFFForcedRange ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.width ());
        should (info.height () == img.height ());
        should (info.isGrayscale ());
        should (info.getPixelType () == rereadType);

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqualTolerance(acc (i) / 255.0, acc (i1) - 1.0, 1e-12);
#endif
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
        should (info.getPixelType () == rereadType);

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
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

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
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif"));
#else
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
#endif
    }

    void testTIFFForcedRange ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0,255,1,2));

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
        {
            shouldEqualTolerance(acc.red(i)/255.0f, acc.red(i1)-1.0f, 1e-4);
            shouldEqualTolerance(acc.green(i)/255.0f, acc.green(i1)-1.0f, 1e-4);
            shouldEqualTolerance(acc.blue(i)/255.0f, acc.blue(i1)-1.0f, 1e-4);
        }
#endif
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

    void testHDR ()
    {
        vigra::ImageExportInfo exi("res.hdr");

        exportImage (srcImageRange (img), exi );

        vigra::ImageImportInfo info ("res.hdr");

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

};

class Vector4ExportImportTest
{
  public:

    typedef vigra::FVector4Image Image;
    typedef vigra::BasicImage<TinyVector<UInt8, 4> > BImage;
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
                img(x,y)[0] = 2*y+x + 0.5f;
                img(x,y)[1] = -img(x,y)[0];
                img(x,y)[2] = 0.0;
                img(x,y)[3] = 0.5;
                for(int b=0; b<4; ++b)
                {
                    breference(x,y)[b] =
                        NumericTraits<UInt8>::fromRealPromote(scale*(img(x,y)[b]+offset));
                }
            }
        }
    }

    void failingTest (char const * filename,
                      char const * message = "exportImage(): file format does not support requested number of bands (color channels)")
    {
        try
        {
            exportImage( srcImageRange(img), vigra::ImageExportInfo( filename ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::PreconditionViolation & e )
        {
            std::string expected = "\nPrecondition violation!\n";
            expected += message;
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failingTest("res.jpg", "did not find a matching codec for the given file extension");
#else
        failingTest("res.jpg");
#endif
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
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif"));
#else
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
#endif
    }

    void testEXR ()
    {
#if !defined(HasEXR)
        failCodec(img, vigra::ImageExportInfo ("res.exr"));
#else
        exportImage (srcImageRange (img),
                     vigra::ImageExportInfo ("res.exr"));

        vigra::ImageImportInfo info ("res.exr");

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
#endif
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
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
#endif
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
#if !defined(HasJPEG)
        testExport("jpg", "did not find a matching codec for the given file extension");
#else
        testExport("jpg");
#endif
    }

    void testJPEGImport()
    {
        testImport("jpg");
    }

    // tiff

    void testTIFFExport()
    {
#if !defined(HasTIFF)
        testExport("tiff", "did not find a matching codec for the given file extension");
#else
        testExport("tiff");
#endif
    }

    void testTIFFImport()
    {
        testImport("tiff");
    }

    // exr

    void testEXRExport()
    {
#if !defined(HasEXR)
        testExport("exr", "did not find a matching codec for the given file extension");
#else
        testExport("exr");
#endif
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
#if !defined(HasPNG)
        testExport("png", "did not find a matching codec for the given file extension");
#else
        testExport("png");
#endif
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

    void testExport( const char * fext ,
                     const char * message = 0)
    {
        std::string fname = "intentionalFailure/foo.";
        fname += fext;
        try {
            exportImage( srcImageRange(img), vigra::ImageExportInfo( fname.c_str() ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::PreconditionViolation & e ) {
            std::string expected = "\nPrecondition violation!\n";
            if(message)
            {
                expected += message;
            }
            else
            {
                expected += "Unable to open file '";
                expected += fname;
                expected += "'.";
            }
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }
};

class GrayscaleImportExportAlphaTest
{
public:
    GrayscaleImportExportAlphaTest()
    {
#if defined(HasTIFF)
        vigra::ImageImportInfo info("lenna_masked_gray.tif");
        
        image_.resize(info.size());
        alpha_.resize(info.size());
        importImageAlpha(info, destImage(image_), destImage(alpha_));
#else
        image_.resize(Size2D(20,10));
        alpha_.resize(image_.size());
        
        image_.init(10);
        alpha_.init(255);
#endif
    }

    void testFile(const char* filename);

    void testTIFF()
    {
        const char filename[] = "res.tif";

#if defined(HasTIFF)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

    void testPNG()
    {
        const char filename[] = "res.png";

#if defined(HasPNG)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

private:
    BImage image_;
    BImage alpha_;
};

void
GrayscaleImportExportAlphaTest::testFile(const char* filename)
{
    exportImageAlpha(srcImageRange(image_), srcImage(alpha_), vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == image_.width());
    should(info.height() == image_.height());
    should(!info.isColor());
    should(!strcmp(info.getPixelType(), "UINT8"));
    should(info.numBands() == 2);

    should(info.width() == alpha_.width());
    should(info.height() == alpha_.height());
    should(info.numExtraBands() == 1);

    BImage image(info.size());
    BImage alpha(info.size());

    importImageAlpha(info, destImage(image), destImage(alpha));

    for (BImage::const_iterator x = alpha_.begin(), xx = alpha_.begin(); x != alpha_.end(); ++x, ++xx)
    {
        should(*x == 255);
        should(*x == *xx);
    }

    for (BImage::const_iterator x = image_.begin(), xx = image.begin(); x != image_.end(); ++x, ++xx)
    {
        should(*x == *xx);
    }
}

class RGBImportExportAlphaTest
{
public:
    RGBImportExportAlphaTest()
    {
#if defined(HasTIFF)
        vigra::ImageImportInfo info("lenna_masked_color.tif");

        image_.resize(info.size());
        alpha_.resize(info.size());
        importImageAlpha(info, destImage(image_), destImage(alpha_));
#else
        image_.resize(Size2D(20,10));
        alpha_.resize(image_.size());
        
        image_.init(RGBValue<unsigned char>(10,20,30));
        alpha_.init(255);
#endif
    }

    void testFile(const char* filename);

    void testTIFF()
    {
        const char filename[] = "res.tif";

#if defined(HasTIFF)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

    void testPNG()
    {
        const char filename[] = "res.png";

#if defined(HasPNG)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

private:
    BRGBImage image_;
    BImage alpha_;
};

void
RGBImportExportAlphaTest::testFile(const char* filename)
{
    exportImageAlpha(srcImageRange(image_), srcImage(alpha_), vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == image_.width());
    should(info.height() == image_.height());
    should(info.isColor());
    should(!strcmp(info.getPixelType(), "UINT8"));
    should(info.numBands() == 4);

    should(info.width() == alpha_.width());
    should(info.height() == alpha_.height());
    should(info.numExtraBands() == 1);

    BRGBImage image(info.size());
    BImage alpha(info.size());

    importImageAlpha(info, destImage(image), destImage(alpha));

    for (BImage::const_iterator x = alpha_.begin(), xx = alpha_.begin(); x != alpha_.end(); ++x, ++xx)
    {
        should(*x == 255);
        should(*x == *xx);
    }

    for (BRGBImage::const_iterator x = image_.begin(), xx = image.begin(); x != image_.end(); ++x, ++xx)
    {
        should(*x == *xx);
    }
}

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
        add(testCase(&ByteImageExportImportTest::testTIFFSequence));
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

#if defined(HasPNG)
        // 16-bit PNG
        add(testCase(&PNGInt16Test::testByteOrder));
#endif

        add(testCase(&CanvasSizeTest::testTIFFCanvasSize));

        // grayscale float images
        add(testCase(&FloatImageExportImportTest::testGIF));
        add(testCase(&FloatImageExportImportTest::testJPEG));
        add(testCase(&FloatImageExportImportTest::testPNG));
        add(testCase(&FloatImageExportImportTest::testTIFF));
        add(testCase(&FloatImageExportImportTest::testTIFFForcedRange));
        add(testCase(&FloatImageExportImportTest::testBMP));
        add(testCase(&FloatImageExportImportTest::testSUN));
        add(testCase(&FloatImageExportImportTest::testVIFF));

        // 4-band images
        add(testCase(&Vector4ExportImportTest::testJPEG));
        add(testCase(&Vector4ExportImportTest::testGIF));
        add(testCase(&Vector4ExportImportTest::testBMP));
        add(testCase(&Vector4ExportImportTest::testPNM));
        add(testCase(&Vector4ExportImportTest::testSUN));
        add(testCase(&Vector4ExportImportTest::testVIFF));
        add(testCase(&Vector4ExportImportTest::testTIFF));
        add(testCase(&Vector4ExportImportTest::testEXR));
        add(testCase(&Vector4ExportImportTest::testPNG));

        // rgb float images
        add(testCase(&FloatRGBImageExportImportTest::testJPEG));
        add(testCase(&FloatRGBImageExportImportTest::testTIFF));
        add(testCase(&FloatRGBImageExportImportTest::testTIFFForcedRange));
        add(testCase(&FloatRGBImageExportImportTest::testBMP));
        add(testCase(&FloatRGBImageExportImportTest::testSUN));
        add(testCase(&FloatRGBImageExportImportTest::testVIFF));
        add(testCase(&FloatRGBImageExportImportTest::testHDR));

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

        // alpha-channel tests
        add(testCase(&GrayscaleImportExportAlphaTest::testTIFF));
        add(testCase(&GrayscaleImportExportAlphaTest::testPNG));
        add(testCase(&RGBImportExportAlphaTest::testTIFF));
        add(testCase(&RGBImportExportAlphaTest::testPNG));
    }
};


int main (int argc, char ** argv)
{
    ImageImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
