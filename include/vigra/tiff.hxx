/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
#ifndef VIGRA_TIFF_HXX
#define VIGRA_TIFF_HXX

#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/tiff.h"

namespace vigra {

/** \defgroup TIFFImpex Import/export of the TIFF format

    TIFF conversion and file export/import.
    
    Normally, you need not call the TIFF functions directly. They are
    available much more conveniently via \ref importImage() and \ref exportImage() 
    
    TIFF (Tagged Image File Format) is a very versatile image format - 
    one can store different pixel types (byte, integer, float, double) and
    color models (black and white, RGB, mapped RGB, other color systems). 
    For more details and information on how to create a TIFF image,
    refer to the TIFF documentation at 
    <a href="http://www.libtiff.org/">http://www.libtiff.org/</a> for details.
*/
//@{

/********************************************************/
/*                                                      */
/*                     importTiffImage                  */
/*                                                      */
/********************************************************/

/** \brief Convert given TiffImage into image specified by iterator range.

    Accessors are used to write the data.    
    This function calls \ref tiffToScalarImage() or \ref tiffToRGBImage(), depending on 
    the accessor's value_type.

    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importTiffImage(TiffImage * tiff, ImageIterator iter, Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importTiffImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    uint32 w, h;
    TiffImage * tiff = TIFFOpen("tiffimage.tiff", "r");
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    
    vigra::BImage img(w,h);
    
    vigra::importTiffImage(tiff, destImage(img));
    
    TIFFClose(tiff);
    \endcode
    
    <b> Required Interface:</b>
    
    see \ref tiffToScalarImage() and \ref tiffToRGBImage()
    
    <b> Preconditions:</b>
    
    see \ref tiffToScalarImage() and \ref tiffToRGBImage()
    
*/
template <class ImageIterator, class Accessor>
inline void
importTiffImage(TiffImage * tiff, ImageIterator iter, Accessor a)
{
    typedef typename 
        NumericTraits<typename Accessor::value_type>::isScalar
        isScalar;
    importTiffImage(tiff, iter, a, isScalar());
}

template <class ImageIterator, class Accessor>
inline void
importTiffImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
{
    importTiffImage(tiff, dest.first, dest.second);
}

template <class ImageIterator, class Accessor>
inline void
importTiffImage(TiffImage * tiff, ImageIterator iter, Accessor a, VigraTrueType)
{
    tiffToScalarImage(tiff, iter, a);
}

template <class ImageIterator, class Accessor>
inline void
importTiffImage(TiffImage * tiff, ImageIterator iter, Accessor a, VigraFalseType)
{
    tiffToRGBImage(tiff, iter, a);
}

/********************************************************/
/*                                                      */
/*                    tiffToScalarImage                 */
/*                                                      */
/********************************************************/

/** \brief Convert single-band TiffImage to scalar image.

    This function uses accessors to write the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        tiffToScalarImage(TiffImage * tiff, ImageIterator iter, Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        tiffToScalarImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    
    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    uint32 w, h;
    uint16 photometric
    TiffImage * tiff = TIFFOpen("tiffimage.tiff", "r");
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tiff_, TIFFTAG_PHOTOMETRIC, &photometric);
        
    if(photometric != PHOTOMETRIC_MINISWHITE &&
       photometric != PHOTOMETRIC_MINISBLACK)
    {
        // not a scalar image - handle error
    }
    
    vigra::BImage img(w,h);
    
    vigra::tiffToScalarImage(tiff, destImage(img));
    
    TIFFClose(tiff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    <unsigned char, short, long, float, double> value;
    
    Accessor accessor;
               
    accessor.set(value, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    ImageIterator must refer to a large enough image.
    
    \code
    uint16 sampleFormat, samplesPerPixel, bitsPerSample, photometric;
           
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);

    sampleFormat != SAMPLEFORMAT_VOID
    samplesPerPixel == 1
    photometric == PHOTOMETRIC_MINISWHITE ||
       photometric == PHOTOMETRIC_MINISBLACK
    bitsPerSample == 1 || 
       bitsPerSample == 8 || 
       bitsPerSample == 16 || 
       bitsPerSample == 32 || 
       bitsPerSample == 64
    
    \endcode
    
*/
template <class ImageIterator, class Accessor>
void
tiffToScalarImage(TiffImage * tiff, ImageIterator iter, Accessor a)
{
    vigra_precondition(tiff != 0, 
             "tiffToScalarImage(TiffImage *, ScalarImageIterator): " 
             "NULL pointer to input data.");
    
    uint16 sampleFormat = 1, bitsPerSample, 
           fillorder, samplesPerPixel, photometric;
    uint32 w,h;
    
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_FILLORDER, &fillorder);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    
    vigra_precondition(photometric == PHOTOMETRIC_MINISWHITE ||
                 photometric == PHOTOMETRIC_MINISBLACK, 
             "tiffToScalarImage(TiffImage *, ScalarImageIterator): " 
             "Image isn't grayscale.");
    
    vigra_precondition(samplesPerPixel == 1, 
             "tiffToScalarImage(TiffImage *, ScalarImageIterator): " 
             "Image is multiband, not scalar.");
    
    vigra_precondition(sampleFormat != SAMPLEFORMAT_VOID,
             "tiffToScalarImage(TiffImage *, ScalarImageIterator): " 
             "undefined pixeltype (SAMPLEFORMAT_VOID).");

    ImageIterator yd(iter);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    int offset, scale, max, min;
    if(photometric == PHOTOMETRIC_MINISWHITE)
    {
        min = 255;
        max = 0;
        scale = -1;
        offset = 255;
    }
    else
    {
        scale = 1;
        offset = 0;
        min = 0;
        max = 255;
    }
    
    try{
        switch(sampleFormat)
        {
          case SAMPLEFORMAT_UINT:
          {
            switch (bitsPerSample)
            {
              case 1:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        if(fillorder == FILLORDER_MSB2LSB)
                        {
                            a.set(((((uint8 *)buf)[x/8] >> (7 - x%8)) & 1) ? max : min, xd);
                        }
                        else
                        {
                            a.set(((((uint8 *)buf)[x/8] >> (x%8)) & 1) ? max : min, xd);
                        }
                    }
                }
                break;
              }
              case 8:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(offset + scale*((uint8 *)buf)[x], xd);
                    }
                }
                break;
              }
              case 16:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((uint16 *)buf)[x], xd);
                    }
                }
                break;
              }
              case 32:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((uint32 *)buf)[x], xd);
                    }
                }
                break;
              }
              default:
                vigra_fail("tiffToScalarImage(TiffImage *, ScalarImageIterator): "
                     "unsupported number of bits per pixel");
            }
            break;
          }
          case SAMPLEFORMAT_INT:
          {
            switch (bitsPerSample)
            {
              case 1:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        if(fillorder == FILLORDER_MSB2LSB)
                        {
                            a.set(((((int8 *)buf)[x/8] >> (7 - x%8)) & 1) ? max : min, xd);
                        }
                        else
                        {
                            a.set(((((int8 *)buf)[x/8] >> (x%8)) & 1) ? max : min, xd);
                        }
                    }
                }
                break;
              }
              case 8:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(offset + scale*((uint8 *)buf)[x], xd);
                    }
                }
                break;
              }
              case 16:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((int16 *)buf)[x], xd);
                    }
                }
                break;
              }
              case 32:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((int32 *)buf)[x], xd);
                    }
                }
                break;
              }
              default:
                vigra_fail("tiffToScalarImage(TiffImage *, ScalarImageIterator): "
                     "unsupported number of bits per pixel");
            }
            break;
          }
          case SAMPLEFORMAT_IEEEFP:
          {
            switch (bitsPerSample)
            {
              case sizeof(float)*8:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((float *)buf)[x], xd);
                    }
                }
                break;
              }
              case sizeof(double)*8:
              {
                for(unsigned int y=0; y<h; ++y, ++yd.y)
                {
                    TIFFReadScanline(tiff, buf, y);
                    ImageIterator xd(yd);

                    for(unsigned int x=0; x<w; ++x, ++xd.x)
                    {
                        a.set(((double *)buf)[x], xd);
                    }
                }
                break;
              }
              default:
                vigra_fail("tiffToScalarImage(TiffImage *, ScalarImageIterator): "
                     "unsupported number of bits per pixel");
            }
            break;
          }
          default:
          {
            // should never happen
            vigra_fail("tiffToScalarImage(TiffImage *, ScalarImageIterator): " 
                 "internal error.");
          }
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class ImageIterator, class Accessor>
void
tiffToScalarImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
{
    tiffToScalarImage(tiff, dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                      tiffToRGBImage                  */
/*                                                      */
/********************************************************/

/** \brief Convert RGB (3-band or color-mapped) TiffImage 
    to RGB image.
    
    This function uses \ref RGBAccessor to write the data.
    A RGBImageIterator is an iterator which is associated with a
    RGBAccessor.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        tiffToRGBImage(TiffImage * tiff, RGBImageIterator iter, RGBAccessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        tiffToRGBImage(TiffImage * tiff, pair<RGBImageIterator, RGBAccessor> dest)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    uint32 w, h;
    uint16 photometric
    TiffImage * tiff = TIFFOpen("tiffimage.tiff", "r");
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tiff_, TIFFTAG_PHOTOMETRIC, &photometric);
        
    if(photometric != PHOTOMETRIC_RGB &&
       photometric != PHOTOMETRIC_PALETTE)
    {
        // not an RGB image - handle error
    }
    
    vigra::BRGBImage img(w, h);
    
    vigra::tiffToRGBImage(tiff, destImage(img));
    
    TIFFClose(tiff);
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    <unsigned char, short, long, float, double> rvalue, gvalue, bvalue;
    
    RGBAccessor accessor;
                           
    accessor.setRed(rvalue, upperleft);
    accessor.setGreen(gvalue, upperleft);
    accessor.setBlue(bvalue, upperleft);
    \endcode
    
    <b> Preconditions:</b>
    
    ImageIterator must refer to a large enough image.
    
    \code
    uint16 sampleFormat, samplesPerPixel, bitsPerSample, photometric;
           
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);

    sampleFormat != SAMPLEFORMAT_VOID
    samplesPerPixel == 3 // unlass photometric == PHOTOMETRIC_PALETTE
    photometric == PHOTOMETRIC_RGB ||
       photometric == PHOTOMETRIC_PALETTE
    bitsPerSample == 1 || 
       bitsPerSample == 8 || 
       bitsPerSample == 16 || 
       bitsPerSample == 32 || 
       bitsPerSample == 64
    \endcode
    
*/
template <class RGBImageIterator, class RGBAccessor>
void
tiffToRGBImage(TiffImage * tiff, RGBImageIterator iter, RGBAccessor a)
{
    vigra_precondition(tiff != 0,
              "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
          "NULL pointer to input data.");
    
    uint16 sampleFormat = 1, bitsPerSample, 
           samplesPerPixel, planarConfig, photometric;
    uint32 w,h;
    
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);
    TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planarConfig);
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    
    vigra_precondition(photometric == PHOTOMETRIC_RGB ||
                 photometric == PHOTOMETRIC_PALETTE, 
             "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
             "Image isn't RGB.");
    
    vigra_precondition(sampleFormat != SAMPLEFORMAT_VOID,
             "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
             "undefined pixeltype (SAMPLEFORMAT_VOID).");
        
    RGBImageIterator yd(iter);
    
    switch (photometric)
    {
      case PHOTOMETRIC_PALETTE:
      {
        uint32 * raster = new uint32[w*h];
        try
        {
            if (!TIFFReadRGBAImage(tiff, w, h, raster, 0)) 
            {
                vigra_fail(
                  "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
                  "unable to read image data.");
            }
          
            for(unsigned int y=0; y<h; ++y, ++yd.y)
            {
                RGBImageIterator xd(yd);

                for(unsigned int x=0; x<w; ++x, ++xd.x)
                {
                    a.setRed(TIFFGetR(raster[x+y*w]), xd);
                    a.setGreen(TIFFGetG(raster[x+y*w]), xd);
                    a.setBlue(TIFFGetB(raster[x+y*w]), xd);
                }
            }
        }
        catch(...)
        {
            delete[] raster;
            throw;
        }
        delete[] raster;
        break;
      }
      case PHOTOMETRIC_RGB:
      {
        vigra_precondition(samplesPerPixel == 3,
                 "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
                 "number of samples per pixel must be 3.");
        
        int bufsize = TIFFScanlineSize(tiff);
        tdata_t * bufr = new tdata_t[bufsize];
        tdata_t * bufg = new tdata_t[bufsize];
        tdata_t * bufb = new tdata_t[bufsize];
        
        int offset = (planarConfig == PLANARCONFIG_CONTIG) ? 3 : 1;
        
        try
        {
            switch(sampleFormat)
            {
              case SAMPLEFORMAT_UINT:
              {
                switch (bitsPerSample)
                {
                  case 8:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        uint8 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (uint8 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (uint8 *)bufr;
                            pg = (uint8 *)bufg;
                            pb = (uint8 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  case 16:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        uint16 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (uint16 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (uint16 *)bufr;
                            pg = (uint16 *)bufg;
                            pb = (uint16 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  case 32:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        uint32 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (uint32 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (uint32 *)bufr;
                            pg = (uint32 *)bufg;
                            pb = (uint32 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  default:
                  {
                    vigra_fail("tiffToRGBImage(TiffImage *, RGBImageIterator): "
                         "unsupported number of bits per pixel");
                  }
                }
                break;
              }
              case SAMPLEFORMAT_INT:
              {
                switch (bitsPerSample)
                {
                  case 8:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        int8 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (int8 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (int8 *)bufr;
                            pg = (int8 *)bufg;
                            pb = (int8 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  case 16:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        int16 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (int16 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (int16 *)bufr;
                            pg = (int16 *)bufg;
                            pb = (int16 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  case 32:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        int32 *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (int32 *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (int32 *)bufr;
                            pg = (int32 *)bufg;
                            pb = (int32 *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  default:
                    vigra_fail("tiffToRGBImage(TiffImage *, RGBImageIterator): "
                         "unsupported number of bits per pixel");
                }
                break;
              }
              case SAMPLEFORMAT_IEEEFP:
              {
                switch (bitsPerSample)
                {
                  case sizeof(float)*8:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        float *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (float *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (float *)bufr;
                            pg = (float *)bufg;
                            pb = (float *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  case sizeof(double)*8:
                  {
                    for(unsigned int y=0; y<h; ++y, ++yd.y)
                    {
                        double *pr, *pg, *pb;
                        
                        if(planarConfig == PLANARCONFIG_CONTIG)
                        {
                            TIFFReadScanline(tiff, bufr, y);
                            pr = (double *)bufr;
                            pg = pr+1;
                            pb = pg+1;
                        }
                        else
                        {
                            TIFFReadScanline(tiff, bufr, y, 0);
                            TIFFReadScanline(tiff, bufg, y, 1);
                            TIFFReadScanline(tiff, bufb, y, 2);
                            pr = (double *)bufr;
                            pg = (double *)bufg;
                            pb = (double *)bufb;
                        }
                        
                        RGBImageIterator xd(yd);

                        for(unsigned int x=0; x<w; ++x, ++xd.x, pr+=offset, pg+=offset, pb+=offset)
                        {
                            a.setRed(*pr, xd);
                            a.setGreen(*pg, xd);
                            a.setBlue(*pb, xd);
                        }
                    }
                    break;
                  }
                  default:
                    vigra_fail("tiffToRGBImage(TiffImage *, RGBImageIterator): "
                         "unsupported number of bits per pixel");
                }
                break;
              }
              default:
              {
                // should never happen
                vigra_fail("tiffToRGBImage(TiffImage *, RGBImageIterator): " 
                     "internal error.");
              }
          }
        }
        catch(...)
        {
            delete[] bufr;
            delete[] bufg;
            delete[] bufb;
            throw;
        }
        delete[] bufr;
        delete[] bufg;
        delete[] bufb;
        
        break;
      }
      default:
      {
        // should never happen
        vigra_fail(
          "tiffToRGBImage(TiffImage *, RGBImageIterator): " 
          "internal error.");
      }
    }
}

template <class ImageIterator, class VectorComponentAccessor>
void
tiffToRGBImage(TiffImage * tiff, pair<ImageIterator, VectorComponentAccessor> dest)
{
    tiffToRGBImage(tiff, dest.first, dest.second);
}

template <class T>
struct CreateTiffImage;

/********************************************************/
/*                                                      */
/*                     createTiffImage                  */
/*                                                      */
/********************************************************/

/** \brief Create a TiffImage from the given iterator range.

    Type and size of the TiffImage are determined by the input image. 
    Currently, the function can create scalar images and RGB images of type 
    unsigned char, short, int, float, and double.
    This function uses accessors to read the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline TiffImage *
        createTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                        Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline TiffImage *
        createTiffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    vigra::BImage img(width, height);
    
    ...
    
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    vigra::createTiffImage(srcImageRange(img), tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    Accessor accessor;
                           
    accessor(upperleft);   // result written into TiffImage
    \endcode
    
*/
template <class ImageIterator, class Accessor>
inline void
createTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
{
    CreateTiffImage<typename Accessor::value_type>::
        exec(upperleft, lowerright, a, tiff);
}

template <class ImageIterator, class Accessor>
inline void
createTiffImage(triple<ImageIterator, ImageIterator, Accessor> src, TiffImage * tiff)
{
    createTiffImage(src.first, src.second, src.third, tiff);
}

/********************************************************/
/*                                                      */
/*                createScalarTiffImage                 */
/*                                                      */
/********************************************************/

/** \brief Create a single-band TiffImage from the given scalar image.

    Type and size of the TiffImage are determined by the input image 
    (may be one of unsigned char, short, int, float, or double).
    This function uses accessors to read the data.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline TiffImage *
        createScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                  Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline TiffImage *
        createScalarTiffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    vigra::BImage img(width, height);
    
    ...
    
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    vigra::createScalarTiffImage(srcImageRange(img), tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    Accessor accessor;
                           
    accessor(upperleft);   // result written into TiffImage
    \endcode
    
*/
template <class ImageIterator, class Accessor>
inline void
createScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
{
    CreateTiffImage<typename Accessor::value_type>::
        exec(upperleft, lowerright, a, tiff);
}

template <class ImageIterator, class Accessor>
inline void
createScalarTiffImage(triple<ImageIterator, ImageIterator, Accessor> src, TiffImage * tiff)
{
    createScalarTiffImage(src.first, src.second, src.third, tiff);
}

template <class ImageIterator, class Accessor>
void
createBScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    ImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            uint8 * p = (uint8 *)buf;
            ImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x)
            {
                p[x] = a(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class ImageIterator, class Accessor>
void
createShortScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    ImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            int16 * p = (int16 *)buf;
            ImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x)
            {
                p[x] = a(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class ImageIterator, class Accessor>
void
createIScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    ImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            int32 * p = (int32 *)buf;
            ImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x)
            {
                p[x] = a(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class ImageIterator, class Accessor>
void
createFScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(float)*8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    ImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            float * p = (float *)buf;
            ImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x)
            {
                p[x] = a(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class ImageIterator, class Accessor>
void
createDScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(double)*8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    ImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            double * p = (double *)buf;
            ImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x)
            {
                p[x] = a(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <>
struct CreateTiffImage<unsigned char>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createBScalarTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<short>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createShortScalarTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<int>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createIScalarTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<float>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createFScalarTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<double>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createDScalarTiffImage(upperleft, lowerright, a, tiff);
    }
};

/********************************************************/
/*                                                      */
/*                  createRGBTiffImage                  */
/*                                                      */
/********************************************************/

/** \brief Create a 3-band TiffImage from the given RGB image.

    Type and size of the TiffImage are determined by the input image 
    (may be one of unsigned char, int, float, or double).
    This function uses \ref RGBAccessor to read the data. A
    RGBImageIterator is an iterator that is associated with a
    RGBAccessor.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        TiffImage *
        createRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a)
                }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        inline TiffImage *
        createRGBTiffImage(triple<RGBImageIterator, RGBImageIterator, RGBAccessor> src)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="tiff_8hxx-source.html">vigra/tiff.hxx</a>"
    
    \code
    vigra::BRGBImage img(width, height);
    
    ...
    
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    vigra::createRGBTiffImage(srcImageRange(img), tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode
    
    <b> Required Interface:</b>
    
    \code
    ImageIterator upperleft;
    RGBAccessor accessor;
                           
    accessor.red(upperleft);     // result written into TiffImage
    accessor.green(upperleft);   // result written into TiffImage
    accessor.blue(upperleft);    // result written into TiffImage
    \endcode
    
*/
template <class RGBImageIterator, class RGBAccessor>
inline void
createRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a, TiffImage * tiff)
{
    CreateTiffImage<typename RGBAccessor::value_type>::
        exec(upperleft, lowerright, a, tiff);
}

template <class RGBImageIterator, class RGBAccessor>
inline void
createRGBTiffImage(triple<RGBImageIterator, RGBImageIterator, RGBAccessor> src, TiffImage * tiff)
{
    createRGBTiffImage(src.first, src.second, src.third, tiff);
}

template <class RGBImageIterator, class RGBAccessor>
void
createBRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    RGBImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            uint8 * pr = (uint8 *)buf;
            uint8 * pg = pr+1;
            uint8 * pb = pg+1;
            
            RGBImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x, pr+=3, pg+=3, pb+=3)
            {
                *pr = a.red(xs);
                *pg = a.green(xs);
                *pb = a.blue(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class RGBImageIterator, class RGBAccessor>
void
createShortRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    RGBImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            uint16 * pr = (uint16 *)buf;
            uint16 * pg = pr+1;
            uint16 * pb = pg+1;
            
            RGBImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x, pr+=3, pg+=3, pb+=3)
            {
                *pr = a.red(xs);
                *pg = a.green(xs);
                *pb = a.blue(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class RGBImageIterator, class RGBAccessor>
void
createIRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    RGBImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            uint32 * pr = (uint32 *)buf;
            uint32 * pg = pr+1;
            uint32 * pb = pg+1;
            
            RGBImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x, pr+=3, pg+=3, pb+=3)
            {
                *pr = a.red(xs);
                *pg = a.green(xs);
                *pb = a.blue(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class RGBImageIterator, class RGBAccessor>
void
createFRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(float)*8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    RGBImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            float * pr = (float *)buf;
            float * pg = pr+1;
            float * pb = pg+1;
            
            RGBImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x, pr+=3, pg+=3, pb+=3)
            {
                *pr = a.red(xs);
                *pg = a.green(xs);
                *pb = a.blue(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <class RGBImageIterator, class RGBAccessor>
void
createDRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright, 
                                   RGBAccessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, sizeof(double)*8);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    
    int bufsize = TIFFScanlineSize(tiff);
    tdata_t * buf = new tdata_t[bufsize];
    
    RGBImageIterator ys(upperleft);
    
    try
    {
        for(int y=0; y<h; ++y, ++ys.y)
        {
            double * pr = (double *)buf;
            double * pg = pr+1;
            double * pb = pg+1;
            
            RGBImageIterator xs(ys);
            
            for(int x=0; x<w; ++x, ++xs.x, pr+=3, pg+=3, pb+=3)
            {
                *pr = a.red(xs);
                *pg = a.green(xs);
                *pb = a.blue(xs);
            }
            TIFFWriteScanline(tiff, buf, y);
        }
    }
    catch(...)
    {
        delete[] buf;
        throw;
    }
    delete[] buf;
}

template <>
struct CreateTiffImage<RGBValue<unsigned char> >
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createBRGBTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<RGBValue<short> >
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createShortRGBTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<RGBValue<int> >
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createIRGBTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<RGBValue<float> >
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createFRGBTiffImage(upperleft, lowerright, a, tiff);
    }
};

template <>
struct CreateTiffImage<RGBValue<double> >
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createDRGBTiffImage(upperleft, lowerright, a, tiff);
    }
};


//@}

} // namespace vigra


#endif /* VIGRA_TIFF_HXX */
