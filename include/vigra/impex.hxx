/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.1.0, Dec 06 2000 )                                    */
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


#ifndef VIGRA_IMPEX_HXX
#define VIGRA_IMPEX_HXX

#include <string>
#include "vigra/impex.h"
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/tiff.hxx"

namespace vigra {

typedef unsigned int
   (*ImageExportFunctionPointer)(VigraImpexImageInfo *image_info,VigraImpexImage *image);

typedef VigraImpexImage * (*ImageImportFunctionPointer)(VigraImpexImageInfo *image_info);

struct ImageFileTypeInfo
{
    char const * typeTag;
    char const * fileExtension;
    char const * fileMagicString;
    int lengthOfMagicString;
    ImageExportFunctionPointer exportFunction;
    ImageImportFunctionPointer importFunction;

};

/** \addtogroup VigraImpex VIGRA's Image Import/Export Facilities
    
    supports GIF, TIFF, JPEG, BMP, PNM (PBM, PGM, PPM), SunRaster, KHOROS-VIFF formats
**/
//@{

/** \brief List the image formats VIGRA can read and write.

    This is useful for creating error messages if VIGRA encounters an
    image format it doesn't recognize.

    <b> Usage:</b>

    <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
    Namespace: vigra

    \code
    std::cout << "supported formats: " << vigra::impexListFormats() << std::endl;
    \endcode

**/
std::string impexListFormats();

/** \brief Test whether a file is an image format known to VIGRA.

    This checks the first few bytes of the file and compares them with the 
    "magic strings" of each recognized image format.

    <b> Usage:</b>

    <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
    Namespace: vigra

    \code
    std::cout << "is image: " << vigra::isImage("foo.bmp") << std::endl;
    \endcode

**/
bool isImage(char const * filename);

/********************************************************/
/*                                                      */
/*                   ImageExportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function exportImage().

    See \ref exportImage() for usage example. This object must
    be used to define the properties of an image to be written to disk.
**/
class ImageExportInfo
{
  public:
        /** Construct ImageExportInfo object.
            The image will be stored under the given filename.
            The file type will be guessed from the extension unless overridden
            by \ref setFileType(). Recognized extensions: '.bmp', '.gif', '.jpeg',
            '.jpg', '.p7', '.pbm', '.pgm', '.pnm', '.ppm',
            '.ras', '.tif', '.tiff', '.xv'. JPEG and TIFF are only available when
            libjpeg and libtiff are installed.

        **/
    ImageExportInfo(char const * filename);

        /** Store image as given file type. This will override any type guessed
            from the file name's extension. Recognized file types:

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>

            With the exception of TIFF and VIFF, all file types store 1 byte (gray scale
            and mapped RGB) or 3 bytes (RGB) per pixel.

            TIFF and VIFF are aditionally able to store
            short and long integers (2 or 4 bytes) and real values (32 bit float and 64 bit
            double) without conversion. So you will need to use TIFF or VIFF if you need to store
            images with high accuracy (the appropriate type to write is automatically derived
            from the image type to be exported). However, many other programs using TIFF
            (e.g. ImageMagick) have not implemented support for those pixel types.
            So don't be surprised if the generated TIFF is not readable in some cases.
            If this happens, convert the image to 'unsigned char' or 'RGBValue<unsigned char>'
            prior to exporting.

        **/
    ImageExportInfo & setFileType(char const * filetype);

        /** Set compression type. This will be ignored if the given compression
            is not supported by the specified file type. Recognized strings:
            "LZW", "RunLength", "1" ... "100". A number is interpreted
            as the compression quality for JPEG compression. JPEG compression is supported
            by the JPEG and TIFF formats.

        **/
    ImageExportInfo & setCompression(char const * compression);

    bool isViff() const
    {
        return (filetype_ != 0) ? strcmp(filetype_->typeTag, "VIFF") == 0 : false;
    }
    bool isTiff() const
    {
        return (filetype_ != 0) ? strcmp(filetype_->typeTag, "TIFF") == 0 : false;
    }
    char const * fileName() const { return filename_.c_str(); }
    void guessFiletypeFromExtension(char const * filename);
    void initImageInfo(VigraImpexImageInfo &) const;

    ImageExportFunctionPointer exportFunction() const;

    std::string filename_;
    ImageFileTypeInfo * filetype_;
    VigraImpexCompressionType compression_;
    int quality_;
};

/********************************************************/
/*                                                      */
/*                   ImageImportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function importImage().

    See \ref importImage() for a usage
    example. This object must be used to read an image from disk and
    enquire about its properties.
**/
class ImageImportInfo
{
  public:
    enum ColorSpace { GRAY = VigraImpexGRAYColorspace,
                      RGB = VigraImpexRGBColorspace,
                      UNDEF = VigraImpexUndefinedColorspace };
    enum PixelType { UINT8, INT16, INT32, FLOAT, DOUBLE };

        /** Construct ImageImportInfo object.
            The image with the given filename is read into memory.
            The file type will be determined by the first few bytes of the
            file (magic number). Recognized file types:

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>

        **/
    ImageImportInfo(char const * filename);

    ~ImageImportInfo();

    void loadImage(char const * filename);
    void deleteInternalImages();

        /** Returns true if the image is a KHOROS VIFF image.
        **/
    bool isViff() const
    {
        return viff_ != 0;
    }
        /** Returns true if the image is a TIFF image.
        **/
    bool isTiff() const
    {
        return tiff_ != 0;
    }

    char const * fileName() const { return filename_.c_str(); }

    bool fileIsLoaded() const { return viff_|| tiff_ || impex_ ; }

        /** Get width of the image.
        **/
    int width() const
    {
        vigra_precondition(fileIsLoaded(), "ImageImportInfo::width(): No image loaded");
        return width_;
    }

        /** Get height of the image.
        **/
    int height() const
    {
        vigra_precondition(fileIsLoaded(), "ImageImportInfo::height(): No image loaded");
        return height_;
    }

        /** Get size of the image.
        **/
    Diff2D size() const
    {
        return Diff2D(width(), height());
    }

        /** Returns true if the image is gray scale.
        **/
    bool isGrayscale() const
    {
        return colorSpace() == GRAY;
    }

        /** Returns true if the image is colored (RGB).
        **/
    bool isColor() const
    {
        return colorSpace() == RGB;
    }

        /** Query the color space of the image. Possible values are:
            <DL>
            <DT>"ImageImportInfo::GRAY"<DD> grayscale image
            <DT>"ImageImportInfo::RGB"<DD> RGB image
            <DT>"ImageImportInfo::UNDEF"<DD> unsupported color space
            </DL>
        **/
    ColorSpace colorSpace() const
    {
        vigra_precondition(fileIsLoaded(), "ImageImportInfo::colorSpace(): No image loaded");
        return colorspace_;
    }

        /** Query the pixel type of the image. Possible values are:
            <DL>
            <DT>"ImageImportInfo::UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"ImageImportInfo::UINT16"<DD> 16-bit signed integer (short)
            <DT>"ImageImportInfo::UINT32"<DD> 32-bit signed integer (long)
            <DT>"ImageImportInfo::FLOAT"<DD> 32-bit floating point (float)
            <DT>"ImageImportInfo::DOUBLE"<DD> 64-bit floating point (double)
            </DL>
        **/
    PixelType pixelType() const
    {
        vigra_precondition(fileIsLoaded(), "ImageImportInfo::pixelType(): No image loaded");
        return pixelType_;
    }

        /** Returns true if the image has 1 byte per pixel (gray) or
            3 bytes per pixel (RGB).
        **/
    bool isByte() const
    {
        return pixelType() == UINT8;
    }

    static ImageFileTypeInfo * findFileTypeFromMagicString(char const * filename);

    const char * getFileType() const { return filetype_->typeTag; }

    std::string filename_;
    ImageFileTypeInfo * filetype_;
    ColorSpace colorspace_;
    int width_, height_;
    PixelType pixelType_;
    ViffImage * viff_;
    TiffImage * tiff_;
    VigraImpexImage * impex_;

  private:
    // do not copy this class
    ImageImportInfo(ImageImportInfo const &) {}
    ImageImportInfo & operator=(ImageImportInfo const &) { return *this; }
};

VigraImpexImage * readVigraImpexImage(char const * filename);

template <class ImageIterator, class Accessor>
void internalImportRGBImage(VigraImpexImage const * image,
                      ImageIterator dul, Accessor ad)
{
    vigra_precondition(image != 0,
                 "importImage(): Null pointer to image data");

    int w = image->columns;
    int h = image->rows;

    ImageIterator dlr = dul + Diff2D(w,h);

    switch (image->c_class)
    {
      case VigraImpexDirectClass:
      {
        VigraImpexRunlengthPacket const * p = image->pixels;
        int run_length = p->length;

        for (; dul.y != dlr.y; ++dul.y)
        {
            ImageIterator dx = dul;
            for(; dx.x != dlr.x; ++dx.x)
            {
                ad.setRed(p->red, dx);
                ad.setGreen(p->green, dx);
                ad.setBlue(p->blue, dx);

                --run_length;
                if(run_length < 0)
                {
                    ++p;
                    run_length = p->length;
                }
            }
        }

        break;
      }

      case VigraImpexPseudoClass:
      {
        VigraImpexRunlengthPacket const * p=image->pixels;
        VigraImpexColorPacket const * colormap = image->colormap;
        int run_length = p->length;

        for (; dul.y != dlr.y; ++dul.y)
        {
            ImageIterator dx = dul;
            for(; dx.x != dlr.x; ++dx.x)
            {
                ad.setRed(colormap[p->index].red, dx);
                ad.setGreen(colormap[p->index].green, dx);
                ad.setBlue(colormap[p->index].blue, dx);

                --run_length;
                if(run_length < 0)
                {
                    ++p;
                    run_length = p->length;
                }
            }
        }
        break;
      }

      default:
        vigra_precondition(0, "importImage(): Unknown image class.");
    }
}

template <class ImageIterator, class Accessor>
void internalImportScalarImage(VigraImpexImage const * image,
                               ImageIterator dul, Accessor ad)
{
    vigra_precondition(image != 0,
                 "importImage(): Null pointer to image data");

    int w = image->columns;
    int h = image->rows;

    ImageIterator dlr = dul + Diff2D(w,h);

    switch (image->c_class)
    {
      case VigraImpexDirectClass:
      {
        VigraImpexRunlengthPacket const * p = image->pixels;
        int run_length = p->length;

        for (; dul.y != dlr.y; ++dul.y)
        {
            ImageIterator dx = dul;
            for(; dx.x != dlr.x; ++dx.x)
            {
                ad.set(p->red, dx);

                --run_length;
                if(run_length < 0)
                {
                    ++p;
                    run_length = p->length;
                }
            }
        }

        break;
      }

      case VigraImpexPseudoClass:
      {
        VigraImpexRunlengthPacket const * p=image->pixels;
        VigraImpexColorPacket const * colormap = image->colormap;
        int run_length = p->length;

        for (; dul.y != dlr.y; ++dul.y)
        {
            ImageIterator dx = dul;
            for(; dx.x != dlr.x; ++dx.x)
            {
                ad.set(colormap[p->index].red, dx);

                --run_length;
                if(run_length < 0)
                {
                    ++p;
                    run_length = p->length;
                }
            }
        }
        break;
      }

      default:
        vigra_precondition(0, "importImage(): Unknown image class.");
    }
}

/********************************************************/
/*                                                      */
/*                     importImage                      */
/*                                                      */
/********************************************************/

/** \brief Read an image, given an \ref vigra::ImageImportInfo object.

    This function uses code adapted from
    <a href="http://www.imagemagick.org/">ImageMagick 4.0</a>
    &copy; 1998 E. I. du Pont de Nemours and Company

    This function uses code adapted from
    <a href="http://www.khoral.com/">KHOROS 1.0</a>
    &copy; 1990, University of New Mexico

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importImage(ImageImportInfo const & image, ImageIterator iter, Accessor a)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        inline void
        importImage(ImageImportInfo const & image, pair<ImageIterator, Accessor> dest)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
    Namespace: vigra

    \code

    vigra::ImageImportInfo info("myimage.gif");

    if(info.isGrayscale())
    {
        // create byte image of appropriate size
        vigra::BImage in(info.width(), info.height());

        vigra::importImage(info, destImage(in)); // read the image
        ...
    }
    else
    {
        // create byte RGB image of appropriate size
        vigra::BRGBImage in(info.width(), info.height());

        vigra::importImage(info, destImage(in)); // read the image
        ...
    }

    \endcode

    <b> Preconditions:</b>

    <UL>

    <LI> the image file must be readable
    <LI> the file type must be one of

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>
    </UL>
**/
template <class ImageIterator, class Accessor>
inline void
importImage(ImageImportInfo const & image, ImageIterator iter, Accessor a)
{
    if(image.isViff())
    {
        importViffImage(image.viff_, iter, a);
    }
    else if(image.isTiff())
    {
        importTiffImage(image.tiff_, iter, a);
    }
    else
    {
        typedef typename
            NumericTraits<typename Accessor::value_type>::isScalar
            isScalar;
        internalImportImage(image.impex_, iter, a, isScalar());
    }
}

template <class ImageIterator, class Accessor>
inline void
importImage(ImageImportInfo const & image, pair<ImageIterator, Accessor> dest)
{
    importImage(image, dest.first, dest.second);
}

template <class ImageIterator, class Accessor>
inline void
internalImportImage(VigraImpexImage const * image,
                    ImageIterator iter, Accessor a, VigraTrueType)
{
    internalImportScalarImage(image, iter, a);
}

template <class ImageIterator, class Accessor>
inline void
internalImportImage(VigraImpexImage const * image,
                   ImageIterator iter, Accessor a, VigraFalseType)
{
    internalImportRGBImage(image, iter, a);
}

template <class SrcIterator, class SrcAccessor>
void internalScalarExportImage(
                 SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info, VigraTrueType)
{
    VigraImpexImageInfo image_info;

    info.initImageInfo(image_info);

    VigraImpexImage * image;

    image = vigraImpexAllocateImage(&image_info);

    vigra_postcondition(image != 0, "exportImage(): Unable to allocate memory");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    image->c_class = VigraImpexDirectClass;
    image->columns = w;
    image->rows = h;
    image->packets = image->columns*image->rows;
    image->pixels = (VigraImpexRunlengthPacket *)
                     malloc(image->packets*sizeof(VigraImpexRunlengthPacket));

    if (image->pixels == 0)
    {
        vigraImpexDestroyImage(image);
        vigraImpexDestroyImageInfo(&image_info);
        vigra_postcondition(0, "exportImage(): Unable to allocate memory");
    }

    VigraImpexRunlengthPacket * p = image->pixels;
    SrcIterator sy = sul;

    for (; sy.y != slr.y; ++sy.y)
    {
        SrcIterator sx = sy;
        for(; sx.x != slr.x; ++sx.x, ++p)
        {
            p->red = sget(sx);
            p->green = sget(sx);
            p->blue = sget(sx);
            p->length = 0;
        }
    }

    ImageExportFunctionPointer exportFunction = info.exportFunction();
    vigra_precondition(exportFunction != 0,
                 "exportImage(): Unknown file type");

    int status = (*exportFunction)(&image_info, image);

    vigraImpexDestroyImage(image);
    vigraImpexDestroyImageInfo(&image_info);

    if(status == 0)
        vigra_postcondition(0, "exportImage(): write vigra_failed");
}

template <class SrcIterator, class SrcAccessor>
void internalScalarExportImage(
                 SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info,
                 VigraFalseType)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    SrcIterator sy = sul;

    typedef typename SrcAccessor::value_type ValueType;

    ValueType min = sget(sul);
    ValueType max = sget(sul);

    for (; sy.y != slr.y; ++sy.y)
    {
        SrcIterator sx = sy;
        for(; sx.x != slr.x; ++sx.x)
        {
            if(sget(sx) < min) min = sget(sx);
            if(max < sget(sx)) max = sget(sx);
        }
    }

    BImage tmp(w,h);

    ValueType offset = -min;

    typename NumericTraits<ValueType>::RealPromote
        scale = 255.0 / (max - min);

    sy = sul;
    BImage::Iterator dy = tmp.upperLeft();

    for (; sy.y != slr.y; ++sy.y, ++dy.y)
    {
        SrcIterator sx = sy;
        BImage::Iterator dx = dy;

        for(; sx.x != slr.x; ++sx.x, ++dx.x)
        {
            *dx = BImage::value_type(scale * (sget(sx) + offset));
        }
    }

    internalScalarExportImage(tmp.upperLeft(), tmp.lowerRight(), tmp.accessor(),
                           info, VigraTrueType());
}

template <class SrcIterator, class SrcAccessor>
void internalRGBExportImage(
                 SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info, VigraTrueType)
{
    VigraImpexImageInfo image_info;

    info.initImageInfo(image_info);

    VigraImpexImage * image;

    image = vigraImpexAllocateImage(&image_info);

    vigra_postcondition(image != 0, "exportImage(): Unable to allocate memory");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    image->c_class = VigraImpexDirectClass;
    image->columns = w;
    image->rows = h;
    image->packets = image->columns*image->rows;
    image->pixels = (VigraImpexRunlengthPacket *)
                     malloc(image->packets*sizeof(VigraImpexRunlengthPacket));

    if (image->pixels == 0)
    {
        vigraImpexDestroyImage(image);
        vigraImpexDestroyImageInfo(&image_info);
        vigra_postcondition(0, "exportImage(): Unable to allocate memory");
    }

    VigraImpexRunlengthPacket * p = image->pixels;
    SrcIterator sy = sul;

    for (; sy.y != slr.y; ++sy.y)
    {
        SrcIterator sx = sy;
        for(; sx.x != slr.x; ++sx.x, ++p)
        {
            p->red = sget.red(sx);
            p->green = sget.green(sx);
            p->blue = sget.blue(sx);
            p->length = 0;
        }
    }

    ImageExportFunctionPointer exportFunction = info.exportFunction();
    vigra_precondition(exportFunction != 0,
                 "exportImage(): Unknown file type");

    int status = (*exportFunction)(&image_info, image);

    vigraImpexDestroyImage(image);
    vigraImpexDestroyImageInfo(&image_info);

    if(status == 0)
        vigra_postcondition(0, "exportImage(): write vigra_failed");
}

template <class SrcIterator, class SrcAccessor>
void internalRGBExportImage(
                 SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info,
                 VigraFalseType)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    SrcIterator sy = sul;

    typedef typename SrcAccessor::value_type RGBValueType;
    typedef typename RGBValueType::value_type RGBComponent;

    RGBComponent min = sget.red(sul);
    RGBComponent max = sget.red(sul);

    for (; sy.y != slr.y; ++sy.y)
    {
        SrcIterator sx = sy;
        for(; sx.x != slr.x; ++sx.x)
        {
            if(sget.red(sx) < min) min = sget.red(sx);
            if(sget.green(sx) < min) min = sget.green(sx);
            if(sget.blue(sx) < min) min = sget.blue(sx);

            if(max < sget.red(sx)) max = sget.red(sx);
            if(max < sget.green(sx)) max = sget.green(sx);
            if(max < sget.blue(sx)) max = sget.blue(sx);
        }
    }

    BRGBImage tmp(w,h);

    RGBValueType offset(-min, -min, -min);

    typename NumericTraits<RGBComponent>::RealPromote
        scale = 255.0 / (max - min);

    sy = sul;
    BRGBImage::Iterator dy = tmp.upperLeft();

    for (; sy.y != slr.y; ++sy.y, ++dy.y)
    {
        SrcIterator sx = sy;
        BRGBImage::Iterator dx = dy;

        for(; sx.x != slr.x; ++sx.x, ++dx.x)
        {
            *dx = BRGBImage::value_type(scale * (sget(sx) + offset));
        }
    }

    internalRGBExportImage(tmp.upperLeft(), tmp.lowerRight(), tmp.accessor(),
                           info, VigraTrueType());
}

template <int N>
struct InternalIsByte
{
    typedef VigraFalseType Ret;
};

template <>
struct InternalIsByte<1>
{
    typedef VigraTrueType Ret;
};


template <class SrcIterator, class SrcAccessor>
inline
void internalExportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info, VigraTrueType)
{
    typedef typename SrcAccessor::value_type Value;
    typedef typename InternalIsByte<sizeof(Value)>::Ret IsByte;
    internalScalarExportImage(sul, slr, sget, info, IsByte());
}

template <class SrcIterator, class SrcAccessor>
inline
void internalExportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info, VigraFalseType)
{
    typedef typename SrcAccessor::value_type Value;
    typedef typename Value::value_type Component;
    typedef typename InternalIsByte<sizeof(Component)>::Ret IsByte;
    internalRGBExportImage(sul, slr, sget, info, IsByte());
}

/********************************************************/
/*                                                      */
/*                     exportImage                      */
/*                                                      */
/********************************************************/

/** \brief Write an image, given an \ref vigra::ImageExportInfo object.

    This function uses code adapted from
    <a href="http://www.wizards.dupont.com/cristy/ImageMagick.html">ImageMagick 4.0</a>
    &copy; 1998 E. I. du Pont de Nemours and Company

    This function uses code adapted from
    <a href="http://www.khoral.com/">KHOROS 1.0</a>
    &copy; 1990, University of New Mexico

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         ImageExportInfo const & info)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         ImageExportInfo const & info)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
    Namespace: vigra

    \code


    vigra::BRGBImage out(w, h);
    ...

    // write as JPEG image, using compression quality 80
    vigra::exportImage(srcImageRange(out),
                      vigra::ImageExportInfo("myimage.jpg").setCompression("80"));

    \endcode

    <b> Preconditions:</b>

    <UL>

    <LI> the image file must be writable
    <LI> the file type must be one of

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>

    </UL>
**/
template <class SrcIterator, class SrcAccessor>
inline
void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                 ImageExportInfo const & info)
{
    if(info.isViff())
    {
        ViffImage * viff = createViffImage(sul, slr, sget);
        int status = writeViffImage((char *)info.fileName(), viff);
        freeViffImage(viff);

        if(status == 0)
            vigra_postcondition(0, "exportImage(): write failed");
    }
    else if(info.isTiff())
    {
        TiffImage * tiff = TIFFOpen((char *)info.fileName(), "w");
        vigra_postcondition(tiff != 0,
               "exportImage(): Unable to open image");

        switch(info.compression_)
        {
            case VigraImpexLZWCompression:
                TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
                break;
            case VigraImpexRunlengthEncodedCompression:
                TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
                break;
            case VigraImpexJPEGCompression:
                TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_JPEG);
                TIFFSetField(tiff, TIFFTAG_JPEGQUALITY, info.quality_);
                break;
            default:
                TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        }
        createTiffImage(sul, slr, sget, tiff);
        TIFFClose(tiff);
    }
    else
    {
        typedef typename
            NumericTraits<typename SrcAccessor::value_type>::isScalar isScalar;
        internalExportImage(sul, slr, sget, info, isScalar());
    }
}

template <class SrcIterator, class SrcAccessor>
inline
void exportImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                 ImageExportInfo const & info)
{
    exportImage(src.first, src.second, src.third, info);
}

//@}

} // namespace vigra

#endif /* VIGRA_IMPEX_HXX */
