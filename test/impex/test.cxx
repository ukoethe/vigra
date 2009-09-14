/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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
#include "unittest.hxx"
#include "vigra/hdf5impex.hxx"
#include "vigra/multi_array.hxx"
#include "time.h"

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

        shouldEqual(formats, "BMP GIF HDR JPEG PNG PNM SUN TIFF VIFF");
        shouldEqual(extensions, "bmp gif hdr jpeg jpg pbm pgm png pnm ppm ras tif tiff xv");
    }

    void testIsImage()
    {
        should(isImage("lenna.xv"));
#ifdef _MSC_VER
        should(!isImage("impex.vcproj"));
#else
        should(!isImage("Makefile"));
#endif
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

        double sum = 0.0;
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

        double sum = 0.0;
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
        should (info.getPixelType () == rereadType);

        Image res (info.width (), info.height ());

        importImage (info, destImage (res));

        Image::ScanOrderIterator i = img.begin ();
        Image::ScanOrderIterator i1 = res.begin ();
        Image::Accessor acc = img.accessor ();

        shouldEqualSequence(i, img.end(), i1);
    }

    void testTIFFForcedRange ()
    {
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

    void testTIFFForcedRange ()
    {
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

class HDF5ExportImportTest
{

public:

    HDF5ExportImportTest()
    {}

	void testUnstridedHDF5ExportImport()
	{
		// export and import data from and to unstrided array

		char hdf5File[] = "testfile1.hdf5";

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		hid_t out_file_1 = H5Fcreate(hdf5File, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		//C++: H5File out_file_1("testfile.hdf5", H5F_ACC_TRUNC);
		// ...data 1
		char hdf5group_1[] = "group/subgroup/subsubgroup/data1";
		writeToHDF5File(out_file_1, hdf5group_1, out_data_1);
		writeToHDF5File2(out_file_1, "group/subgroup/subsubgroup/data1_B", out_data_1);
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data2";
		writeToHDF5File(out_file_1, hdf5group_2, out_data_2);
		H5Fclose(out_file_1);
		//C++: out_file_1.close();

		// import
		// ...data 1
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,int> in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_1, in_data_1);
		// ...data 2
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2);

		// compare content
		// ...data 1
		should (in_data_1 == out_data_1);
		// ...data 2
		should (in_data_2 == out_data_2);
	}


	void testStridedHDF5ExportImport1()
	{
		// export data from strided arrays and import to unstrided arrays

		char hdf5File[] = "testfile2.hdf5";

		// int data in 2 dimensions (partly negative)
		MultiArray<3,int> out_data_3(MultiArrayShape<3>::type(2, 3, 4));
        // initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_3.data () [i] = i;
		//std::cout << "Test (0,0,0), (1,0,0), (0,1,0), (0,0,1): " << out_data_3(0,0,0) << " " << out_data_3(1,0,0) << " " << out_data_3(0,1,0) << " " << out_data_3(0,0,1) << std::endl;
		// bind inner dimension to test if strided data impex works
		MultiArrayView<2,int,StridedArrayTag> out_data_4(out_data_3.bindInner(1));
		// bind dimension in the middle to test if strided data impex works
		MultiArrayView<2,int,StridedArrayTag> out_data_5(out_data_3.bindAt(1, 1));

		// export the two sets
		hid_t out_file_2 = H5Fcreate(hdf5File, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		//C++: H5File out_file_2("testfile2.hdf5", H5F_ACC_TRUNC);
		char hdf5group_4[] = "group/subgroup/data4";
		//std::cout << "data (0,0),(0,1),(1,0): " << out_data_4(0,0) << " " << out_data_4(0,1) << " " << out_data_4(1,0) << std::endl;

		clock_t t1,t2;
		t1 = clock();
		writeToHDF5File(out_file_2, hdf5group_4, out_data_4);
		t2 = clock();
		std::cout << "Time 1: " << ((float)(t2-t1)) << std::endl;
		t1 = clock();
		writeToHDF5File2(out_file_2, "group/subgroup/data4_B", out_data_4);
		t2 = clock();
		std::cout << "Time 2: " << ((float)(t2-t1)) << std::endl;

		char hdf5group_5[] = "group/subgroup/data5";
		//std::cout << "data (0,0),(0,1),(1,0): " << out_data_4(0,0) << " " << out_data_4(0,1) << " " << out_data_4(1,0) << std::endl;
		writeToHDF5File(out_file_2, hdf5group_5, out_data_5);
		writeToHDF5File(out_file_2, "group/subgroup/data5_B", out_data_5);
		H5Fclose(out_file_2);
		//C++: out_file_2.close();

		// import test: copy data to unstrided array
		HDF5ImportInfo infoHDF5_4(hdf5File, hdf5group_4);
        MultiArray<2,int> in_data_4(MultiArrayShape<2>::type(infoHDF5_4.shapeOfDimension(0), infoHDF5_4.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_4, in_data_4);
		HDF5ImportInfo infoHDF5_5(hdf5File, hdf5group_5);
        MultiArray<2,int> in_data_5(MultiArrayShape<2>::type(infoHDF5_5.shapeOfDimension(0), infoHDF5_5.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_5, in_data_5);
		// compare content
		should (in_data_4 == out_data_4);
		should (in_data_5 == out_data_5);
	}

	void testStridedHDF5ExportImport2()
	{
		// export data from unstrided arrays and import to strided arrays

		char hdf5File[] = "testfile3.hdf5";

		// int data in 2 dimensions (partly negative)
		MultiArray<3,int> out_data_6(MultiArrayShape<3>::type(2, 3, 4));
        // initialize the array to the test data
        for (int i = 0; i < 24; ++i)
            out_data_6.data () [i] = i;

		// export the two sets
		hid_t out_file_3 = H5Fcreate(hdf5File, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		char hdf5group_6[] = "group/subgroup/data6";
		writeToHDF5File(out_file_3, hdf5group_6, out_data_6);
		H5Fclose(out_file_3);

		// import test: copy data to strided array
		HDF5ImportInfo infoHDF5_6(hdf5File, hdf5group_6);
        MultiArray<4,int> in_data_6a(MultiArrayShape<4>::type(3, infoHDF5_6.shapeOfDimension(0), infoHDF5_6.shapeOfDimension(1), infoHDF5_6.shapeOfDimension(2)));
		in_data_6a.init(42);
		MultiArrayView<3,int,StridedArrayTag> in_data_6b(in_data_6a.bindInner(0));
        loadFromHDF5File(infoHDF5_6, in_data_6b);
		// compare content
		should (in_data_6b == out_data_6);
	}

	void testAppendHDF5ExportImport()
	{
		// write data to file, close file, open file, write more data, compare

		char hdf5File[] = "testfile4.hdf5";

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		hid_t out_file_1 = H5Fcreate(hdf5File, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		// ...data 1
		char hdf5group_1[] = "group/subgroup/subsubgroup/data1";
		writeToHDF5File(out_file_1, hdf5group_1, out_data_1);
		H5Fclose(out_file_1);

		// append to existing file
		hid_t out_file_2 = H5Fopen(hdf5File, H5F_ACC_RDWR, H5P_DEFAULT);
		// ...data 2
		char hdf5group_2[] = "group/subgroup/data2";
		writeToHDF5File(out_file_2, hdf5group_2, out_data_2);
		H5Fclose(out_file_2);

		// import
		// ...data 1
		HDF5ImportInfo infoHDF5_1(hdf5File, hdf5group_1);
        MultiArray<2,int> in_data_1(MultiArrayShape<2>::type(infoHDF5_1.shapeOfDimension(0), infoHDF5_1.shapeOfDimension(1)));
        loadFromHDF5File(infoHDF5_1, in_data_1);
		// ...data 2
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2);

		// compare content
		// ...data 1
		should (in_data_1 == out_data_1);
		// ...data 2
		should (in_data_2 == out_data_2);
	}

	void testOverwriteHDF5ExportImport()
	{
		// write data to file, close file, open file, overwrite data, compare

		char hdf5File[] = "testfile5.hdf5";

		// data 1: int data in 2 dimensions (partly negative)
		MultiArray<2,int> out_data_1(MultiArrayShape<2>::type(10, 11));
        // ...initialize the array to the test data
        for (int i = 0; i < 110; ++i)
            out_data_1.data () [i] = i - 55;
		//std::cout << "Test (0,0), (0,1), (1,0): " << out_data_1(0,0) << " " << out_data_1(0,1) << " " << out_data_1(1,0) << " " << std::endl;

		// data 2: double data in 4 dimensions (partly negative)
		MultiArray<4,double> out_data_2(MultiArrayShape<4>::type(10, 2, 3, 4));
        // ...initialize the array to the test data
        for (int i = 0; i < 240; ++i)
            out_data_2.data () [i] = i + (std::rand() / (double)RAND_MAX) - 120;

		// export
		hid_t out_file_1 = H5Fcreate(hdf5File, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		// ...data 1
		char hdf5group_1[] = "group/subgroup/data1";
		writeToHDF5File(out_file_1, hdf5group_1, out_data_1);
		H5Fclose(out_file_1);

		// overwrite in existing file
		// ...data 2
		hid_t out_file_2 = H5Fopen(hdf5File, H5F_ACC_RDWR, H5P_DEFAULT);
		char hdf5group_2[] = "group/subgroup/data1";
		writeToHDF5File(out_file_2, hdf5group_1, out_data_2);
		H5Fclose(out_file_2);

		// import
		HDF5ImportInfo infoHDF5_2(hdf5File, hdf5group_2);
        MultiArray<4,double> in_data_2(MultiArrayShape<4>::type(infoHDF5_2.shapeOfDimension(0), infoHDF5_2.shapeOfDimension(1), infoHDF5_2.shapeOfDimension(2), infoHDF5_2.shapeOfDimension(3)));
		loadFromHDF5File(infoHDF5_2, in_data_2);

		// compare content
		should (in_data_2 == out_data_2);
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
    }
};

struct HDF5ImportExportTestSuite : public vigra::test_suite
{
    HDF5ImportExportTestSuite()
        : vigra::test_suite("HDF5ImportExportTestSuite")
    {
		// general tests
        add(testCase(&HDF5ExportImportTest::testUnstridedHDF5ExportImport));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport1));
        add(testCase(&HDF5ExportImportTest::testStridedHDF5ExportImport2));
        add(testCase(&HDF5ExportImportTest::testAppendHDF5ExportImport));
        add(testCase(&HDF5ExportImportTest::testOverwriteHDF5ExportImport));
	}
};


int main (int argc, char ** argv)
{
    ImageImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    HDF5ImportExportTestSuite test2;
    const int failed2 = test2.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test2.report() << std::endl;

	return failed*failed2 != 0;
}
