/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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
 

#ifndef VIGRA_IMPEX_HXX
#define VIGRA_IMPEX_HXX

#include <string>
#include "vigra/impex.h"
#include "vigra/stdimage.hxx"
#include "vigra/viff.hxx"
#include "vigra/tiff.hxx"

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

std::string vigraImpexListFormats();

/** @name VIGRA's Image Import/Export Facilities
    @memo supports GIF, TIFF, JPEG, BMP, SunRaster, KHOROS-VIFF formats
**/
//@{

/********************************************************/
/*                                                      */
/*                   ImageExportInfo                    */
/*                                                      */
/********************************************************/

/** Argument object for the function \Ref{exportImage} (see there for usage
    example). This object must
    be used to define the properties of an image to be written to disk. 
    
    @memo Argument object for the function \Ref{exportImage}
**/
class ImageExportInfo
{
  public:
        /** Construct ImageExportInfo object.
            The image will be stored under the given filename.
            The file type will be guessed from the extension unless overridden
            by \Ref{setFileType}. Recognized extensions: .gif .jpg .jpeg
            .bmp .ras .tif .tiff .xv. JPEG and TIFF are only available when
            libjpeg and libtiff are installed. 

            @memo 
        **/
    ImageExportInfo(char const * filename);
    
        /** Store image as given file type. This will override any type guessed
            from the file name's extension. Recognized types: "GIF", "JPEG",
            "BMP" (Windows bitmap), "SUN" (SUN raster format), "TIFF", 
            "VIFF" (KHOROS format). JPEG and TIFF are only available when
            libjpeg and libtiff are installed.
            
            With the exception of VIFF, all file types store 1 byte (gray scale
            and mapped RGB) or 3 bytes (RGB) per pixel. VIFF is able to store
            long integers (4 bytes) and real values (float and double) without
            conversion. So you will need to use VIFF if you need to store
            images with high accuracy.

            @memo 
        **/
    ImageExportInfo & setFileType(char const * filetype);
    
        /** Set compression type. This will be ignored if the given compression
            is not supported by the specified file type. Recognized strings:
            "LZW", "RunLength", "1" ... "100". A number is interpreted 
            as the compression quality for JPEG compression.

            @memo 
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

/** Argument object for the function \Ref{importImage} (see there for usage
    example). This object must be used to read an image from disk and
    enquire about its properties.
    
    @memo Argument object for the function \Ref{importImage}
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
            file (magic number). Recognized file types: "GIF", "JPEG",
            "BMP" (Windows bitmap), "SUN" (SUN raster format), "TIFF", 
            "VIFF" (KHOROS format).  JPEG and TIFF are only available when
            libjpeg and libtiff are installed. 

            @memo 
        **/
    ImageImportInfo(char const * filename);

    ~ImageImportInfo();
    
    void loadImage(char const * filename);
    void deleteInternalImages();
    
        /** Returns true if the image has type KHOROS VIFF image.

            @memo 
        **/
    bool isViff() const 
    { 
        return viff_ != 0; 
    }
    bool isTiff() const 
    { 
        return tiff_ != 0; 
    }
    char const * fileName() const { return filename_.c_str(); }
    
    bool fileIsLoaded() const { return viff_|| tiff_ || impex_ ; }

        /** Get width of the image.

            @memo 
        **/
    int width() const
    {
        precondition(fileIsLoaded(), "ImageImportInfo::width(): No image loaded");
        return width_;
    }
        /** Get height of the image.

            @memo 
        **/
    int height() const
    {
        precondition(fileIsLoaded(), "ImageImportInfo::height(): No image loaded");
        return height_;
    }
        /** Get size of the image.

            @memo 
        **/
    Diff2D size() const
    {
        return Diff2D(width(), height());
    }
        /** Returns true if the image is gray scale.

            @memo 
        **/
    bool isGrayscale() const
    {
        return colorSpace() == GRAY;
    }
        /** Returns true if the image is colored (RGB).

            @memo 
        **/
    bool isColor() const
    {
        return colorSpace() == RGB;
    }
    
    ColorSpace colorSpace() const 
    {
        precondition(fileIsLoaded(), "ImageImportInfo::colorSpace(): No image loaded");
        return colorspace_; 
    }
    
    PixelType pixelType() const 
    { 
        precondition(fileIsLoaded(), "ImageImportInfo::pixelType(): No image loaded");
        return pixelType_; 
    }
    
        /** Returns true if the image 1 byte per pixel (gray) or
            3 bytes per pixel (RGB).

            @memo 
        **/
    bool isByte() const
    {
        return pixelType() == UINT8;
    }
    
    void findFileTypeFromMagicString(char const * filename);
    
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
    precondition(image != 0,
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
        precondition(0, "importImage(): Unknown image class.");
    }
}

template <class ImageIterator, class Accessor>
void internalImportScalarImage(VigraImpexImage const * image, 
                               ImageIterator dul, Accessor ad)
{
    precondition(image != 0,
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
        precondition(0, "importImage(): Unknown image class.");
    }
}

/********************************************************/
/*                                                      */
/*                     importImage                      */
/*                                                      */
/********************************************************/

/** Read an image, given an \Ref{ImageImportInfo} object.

    This function uses code adapted from 
    \URL[ImageMagick 4.0]{http://www.wizards.dupont.com/cristy/ImageMagick.html}
    &copy; 1998 E. I. du Pont de Nemours and Company

    This function uses code adapted from 
    \URL[KHOROS 1.0]{http://www.khoral.com/}
    &copy; 1990, University of New Mexico

    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class ImageIterator, class Accessor>
    void
    importImage(ImageImportInfo const & image, ImageIterator iter, Accessor a)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class ImageIterator, class Accessor>
    inline void
    importImage(ImageImportInfo const & image, pair<ImageIterator, Accessor> dest)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/impex.hxx]{../include/vigra/impex.hxx}
    
    \begin{verbatim}
    
    ImageImportInfo info("myimage.gif");
    
    if(info.isGrayscale())
    {
        // create byte image of appropriate size
        BImage in(info.width(), info.height()); 
        
        importImage(info, destImage(in)); // read the image
        ...
    }
    else
    {
        // create byte RGB image of appropriate size
        BRGBImage in(info.width(), info.height()); 
        
        importImage(info, destImage(in)); // read the image
        ...
    }
    
    \end{verbatim}

    {\bf Preconditions:}
    
    \begin{itemize}
    
    \item the image file must be readable 
    \item the file type must be one of "GIF", "JPEG",
            "BMP" (Windows bitmap), "SUN" (SUN raster format), "TIFF", 
            "VIFF" (KHOROS format).  JPEG and TIFF are only available when
            libjpeg and libtiff are installed. 
    
    \end{itemize}
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
    
    postcondition(image != 0, "exportImage(): Unable to allocate memory");
    
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
        postcondition(0, "exportImage(): Unable to allocate memory");
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
    precondition(exportFunction != 0, 
                 "exportImage(): Unknown file type");
    
    int status = (*exportFunction)(&image_info, image);
    
    vigraImpexDestroyImage(image);
    vigraImpexDestroyImageInfo(&image_info);
    
    if(status == 0)
        postcondition(0, "exportImage(): write failed");
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
    
    postcondition(image != 0, "exportImage(): Unable to allocate memory");
    
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
        postcondition(0, "exportImage(): Unable to allocate memory");
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
    precondition(exportFunction != 0, 
                 "exportImage(): Unknown file type");
    
    int status = (*exportFunction)(&image_info, image);
    
    vigraImpexDestroyImage(image);
    vigraImpexDestroyImageInfo(&image_info);
    
    if(status == 0)
        postcondition(0, "exportImage(): write failed");
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

/** Write an image, given an \Ref{ImageExportInfo} object.

    This function uses code adapted from 
    \URL[ImageMagick 4.0]{http://www.wizards.dupont.com/cristy/ImageMagick.html}
    &copy; 1998 E. I. du Pont de Nemours and Company

    This function uses code adapted from 
    \URL[KHOROS 1.0]{http://www.khoral.com/}
    &copy; 1990, University of New Mexico

    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor>
    void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                     ImageExportInfo const & info)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor>
    void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                     ImageExportInfo const & info)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/impex.hxx]{../include/vigra/impex.hxx}
    
    \begin{verbatim}
    
    
    BRGBImage out(w, h); 
    ...
        
    // write as JPEG image, using compression quality 80
    exportImage(srcImageRange(out), 
                ImageExportInfo("myimage.jpg").setCompression("80")); 
    
    \end{verbatim}

    {\bf Preconditions:}
    
    \begin{itemize}
    
    \item the image file must be writable 
    \item the file type must be one of "GIF", "JPEG",
            "BMP" (Windows bitmap), "SUN" (SUN raster format), "TIFF", 
            "VIFF" (KHOROS format).  JPEG and TIFF are only available when
            libjpeg and libtiff are installed. 
    
    \end{itemize}
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
            postcondition(0, "exportImage(): write failed");
    }
    else if(info.isTiff())
    {
        TiffImage * tiff = TIFFOpen((char *)info.fileName(), "w");
        postcondition(tiff != 0, 
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


#endif /* VIGRA_IMPEX_HXX */
