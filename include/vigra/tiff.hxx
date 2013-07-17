/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_TIFF_HXX
#define VIGRA_TIFF_HXX

#include "utilities.hxx"
#include "numerictraits.hxx"
#include "rgbvalue.hxx"
#include "multi_shape.hxx"

extern "C"
{
#include <tiff.h>
#include <tiffio.h>
}

namespace vigra {

typedef TIFF TiffImage;

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

/** \brief Read a given TIFF image.

    This function calls \ref tiffToScalarImage() or \ref tiffToRGBImage(), depending on 
    the destinations's value_type. Usually, it is better to use \ref importImage().
    importTiffImage() should only be used if explicit access to the TIFF object
    <tt>TiffImage</tt> is required.

    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S>
        void
        importTiffImage(TiffImage * tiff, MultiArrayView<2, T, S> dest);
    }
    \endcode
    
    \deprecatedAPI{importTiffImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importTiffImage(TiffImage * tiff, ImageIterator iter, Accessor a)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importTiffImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra
    
    \code
    uint32 w, h;
    TiffImage * tiff = TIFFOpen("tiffimage.tiff", "r");
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    
    MultiArray<2, unsigned char> img(w,h);
    
    importTiffImage(tiff, img);
    
    TIFFClose(tiff);
    \endcode
    
    <b> Required Interface:</b>
    
    see \ref tiffToScalarImage() and \ref tiffToRGBImage()
    
    <b> Preconditions:</b>
    
    see \ref tiffToScalarImage() and \ref tiffToRGBImage()
    
*/
doxygen_overloaded_function(template <...> void importTiffImage)

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

template <class T, class S>
inline void
importTiffImage(TiffImage * tiff, MultiArrayView<2, T, S> dest)
{
    importTiffImage(tiff, destImage(dest));
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

    Note that unexpected results can occur when the destination pixel type is weaker than the pixel type
    in the file (e.g. when a <tt>float</tt> file is imported into a <tt>unsigned char</tt> image).

    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        tiffToScalarImage(TiffImage * tiff, ImageIterator iter, Accessor a)
    }
    \endcode
    
    \deprecatedAPI{tiffToScalarImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        tiffToScalarImage(TiffImage * tiff, ImageIterator iter, Accessor a)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        tiffToScalarImage(TiffImage * tiff, pair<ImageIterator, Accessor> dest)
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra

    \code
    uint32 w, h;
    uint16 photometric;
    TiffImage * tiff = TIFFOpen("tiffimage.tiff", "r");
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tiff_, TIFFTAG_PHOTOMETRIC, &photometric);
        
    if(photometric != PHOTOMETRIC_MINISWHITE &&
       photometric != PHOTOMETRIC_MINISBLACK)
    {
        // not a scalar image - handle error
    }
    
    MultiArray<2, unsigned char> img(w,h);
    
    tiffToScalarImage(tiff, img);
    
    TIFFClose(tiff);
    \endcode

    \deprecatedUsage{tiffToScalarImage}
    \code
    uint32 w, h;
    uint16 photometric;
    
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
    \deprecatedEnd
    
    <b> Preconditions:</b>
    
    The output array must have the correct shape.
    
    \code
    uint16 sampleFormat, samplesPerPixel, bitsPerSample, photometric;
           
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);

    sampleFormat != SAMPLEFORMAT_VOID
    samplesPerPixel == 1
    photometric == PHOTOMETRIC_MINISWHITE || photometric == PHOTOMETRIC_MINISBLACK
    bitsPerSample == 1 || 
       bitsPerSample == 8 || 
       bitsPerSample == 16 || 
       bitsPerSample == 32 || 
       bitsPerSample == 64
    \endcode
*/
doxygen_overloaded_function(template <...> void tiffToScalarImage)

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

/** \brief Import a RGB (3-band or color-mapped) TiffImage 
    into a RGB image.

    Note that unexpected results can occur when the destination pixel type is weaker than the pixel type
    in the file (e.g. when a <tt>float</tt> file is imported into a <tt>unsigned char</tt> image).
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        tiffToRGBImage(TiffImage * tiff, RGBImageIterator iter, RGBAccessor a)
    }
    \endcode
    
    \deprecatedAPI{tiffToRGBImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        tiffToRGBImage(TiffImage * tiff, RGBImageIterator iter, RGBAccessor a)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        void
        tiffToRGBImage(TiffImage * tiff, pair<RGBImageIterator, RGBAccessor> dest)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra

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
    
    MultiArray<2, RGBValue<unsigned char> > img(w, h);
    
    tiffToRGBImage(tiff, img);
    
    TIFFClose(tiff);
    \endcode

    \deprecatedUsage{tiffToRGBImage}
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
    \deprecatedEnd
    
    <b> Preconditions:</b>
    
    The destination image must have the appropriate size.
    
    \code
    uint16 sampleFormat, samplesPerPixel, bitsPerSample, photometric;
           
    TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samplesPerPixel);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &photometric);

    sampleFormat != SAMPLEFORMAT_VOID
    samplesPerPixel == 3 // unless photometric == PHOTOMETRIC_PALETTE
    photometric == PHOTOMETRIC_RGB ||
       photometric == PHOTOMETRIC_PALETTE
    bitsPerSample == 1 || 
       bitsPerSample == 8 || 
       bitsPerSample == 16 || 
       bitsPerSample == 32 || 
       bitsPerSample == 64
    \endcode
*/
doxygen_overloaded_function(template <...> void tiffToRGBImage)

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
                typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                typename RGBImageIterator::row_iterator rowend = rowit + w;
                for(int x=0; rowit<rowend; ++rowit,++x )
                {
                    uint32 rast = raster[x+y*w];
                    a.setRGB(TIFFGetR(rast),TIFFGetG(rast),TIFFGetB(rast),rowit);
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
                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                                                                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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

                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
                        
                        typename RGBImageIterator::row_iterator rowit = yd.rowIterator();
                        typename RGBImageIterator::row_iterator rowend = rowit + w;
                        for(; rowit<rowend; ++rowit, pr+=offset, pg+=offset, pb+=offset)
                            a.setRGB(*pr,*pg, *pb, rowit);
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
   
    Usually, it is better to use \ref exportImage(). createTiffImage() should only be used if explicit access to the TIFF object
    <tt>TiffImage</tt> is required.    
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T, class S>
        void
        createTiffImage(MultiArrayView<2, T, S> const & src, TiffImage * tiff);
    }
    \endcode
    
    \deprecatedAPI{createTiffImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        TiffImage *
        createTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                        Accessor a)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        TiffImage *
        createTiffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> img(width, height);
    
    ...
    
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    createTiffImage(img, tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode

    \deprecatedUsage{createTiffImage}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void createTiffImage)

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

template <class T, class S>
inline void
createTiffImage(MultiArrayView<2, T, S> const & src, TiffImage * tiff)
{
    createTiffImage(srcImageRange(src), tiff);
}

/********************************************************/
/*                                                      */
/*                createScalarTiffImage                 */
/*                                                      */
/********************************************************/

/** \brief Create a single-band TiffImage from the given scalar image.

    Type and size of the TiffImage are determined by the input image 
    (may be one of unsigned char, short, int, float, or double).
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        TiffImage *
        createScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                  Accessor a)
    }
    \endcode
    
    \deprecatedAPI{createScalarTiffImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        TiffImage *
        createScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                  Accessor a)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        TiffImage *
        createScalarTiffImage(triple<ImageIterator, ImageIterator, Accessor> src)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> img(width, height);
    ...
   
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    createScalarTiffImage(img, tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode

    \deprecatedUsage{createScalarTiffImage}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void createScalarTiffImage)

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
createUShortScalarTiffImage(ImageIterator upperleft, ImageIterator lowerright, 
                                 Accessor a, TiffImage * tiff)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, 16);
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
            uint16 * p = (uint16 *)buf;
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
struct CreateTiffImage<unsigned short>
{
    template <class ImageIterator, class Accessor>
    static void
    exec(ImageIterator upperleft, ImageIterator lowerright, 
                      Accessor a, TiffImage * tiff)
    {
        createUShortScalarTiffImage(upperleft, lowerright, a, tiff);
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
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        TiffImage *
        createRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a)
                }
    \endcode
    
    \deprecatedAPI{createRGBTiffImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        TiffImage *
        createRGBTiffImage(RGBImageIterator upperleft, RGBImageIterator lowerright,
                   RGBAccessor a)
                }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class RGBImageIterator, class RGBAccessor>
        TiffImage *
        createRGBTiffImage(triple<RGBImageIterator, RGBImageIterator, RGBAccessor> src)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tiff.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, RGBValue<unsigned char> > img(width, height);
    ...
    
    TiffImage * tiff = TIFFOpen(("tiffimage.tiff", "w");

    createRGBTiffImage(img, tiff);

    TIFFClose(tiff);   // implicitly writes the image to the disk
    \endcode

    \deprecatedUsage{createRGBTiffImage}
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
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void createRGBTiffImage)

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
