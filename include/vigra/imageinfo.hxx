/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*               Copyright 2001-2002 by Gunnar Kedenburg                */
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

#ifndef VIGRA_IMAGEINFO_HXX
#define VIGRA_IMAGEINFO_HXX

#include <memory>
#include <string>
#include "vigra/utilities.hxx"
#include "vigra/codec.hxx"

namespace vigra
{
/** \addtogroup VigraImpex Image Import/Export Facilities
    
    supports GIF, TIFF, JPEG, BMP, PNM (PBM, PGM, PPM), PNG, SunRaster, KHOROS-VIFF formats
**/
//@{

    /** \brief List the image formats VIGRA can read and write.

        This is useful for creating error messages if VIGRA encounters an
        image format it doesn't recognize.

        <b> Usage:</b>

        <b>\#include</b> "<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>"<br>
        Namespace: vigra

        \code
        std::cout << "supported formats: " << vigra::impexListFormats() << std::endl;
        \endcode

    **/
std::string impexListFormats();

    /** \brief List the file extension VIGRA understands.

        This is useful for creating file dialogs that only list image files
        VIGRA can actually import.

        <b> Usage:</b>

        <b>\#include</b> "<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>"<br>
        Namespace: vigra

        \code
        std::cout << "supported extensions: " << vigra::impexListExtensions() << std::endl;
        \endcode

    **/
std::string impexListExtensions();

/** \brief Test whether a file is an image format known to VIGRA.

    This checks the first few bytes of the file and compares them with the
    "magic strings" of each recognized image format.

    <b> Usage:</b>

    <b>\#include</b> "<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>"<br>
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
    See \ref exportImage() for usage example. This object must be used
    to define the properties of an image to be written to disk.

    <b>\#include</b> "<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>"<br>
    Namespace: vigra
**/
class ImageExportInfo
{
  public:
        /** Construct ImageExportInfo object.
        
            The image will be stored under the given filename.
            The file type will be guessed from the extension unless overridden
            by \ref setFileType(). Recognized extensions: '.bmp', '.gif',
            '.jpeg', '.jpg', '.p7', '.png', '.pbm', '.pgm', '.pnm', '.ppm', '.ras',
            '.tif', '.tiff', '.xv'.
            JPEG support requires libjpeg, PNG support requires libpng, and
            TIFF support requires libtiff.
        **/
    ImageExportInfo( const char * );

    const char * getFileName() const;

    /** Store image as given file type.
    
        This will override any type guessed
        from the file name's extension. Recognized file types:

        <DL>
        <DT>"BMP"<DD> Microsoft Windows bitmap image file.
        <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
        <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format;
        compressed 24-bit color (only available if libjpeg is installed).
        <DT>"PNG"<DD> Portable Network Graphic
        (only available if libpng is installed).
        <DT>"PBM"<DD> Portable bitmap format (black and white).
        <DT>"PGM"<DD> Portable graymap format (gray scale).
        <DT>"PNM"<DD> Portable anymap.
        <DT>"PPM"<DD> Portable pixmap format (color).
        <DT>"SUN"<DD> SUN Rasterfile.
        <DT>"TIFF"<DD> Tagged Image File Format.
        (only available if libtiff is installed.)
        <DT>"VIFF"<DD> Khoros Visualization image file.
        </DL>

        With the exception of TIFF and VIFF, all file types store
        1 byte (gray scale and mapped RGB) or 3 bytes (RGB) per
        pixel.

        TIFF and VIFF are aditionally able to store short and long
        integers (2 or 4 bytes) and real values (32 bit float and
        64 bit double) without conversion. So you will need to use
        TIFF or VIFF if you need to store images with high
        accuracy (the appropriate type to write is automatically
        derived from the image type to be exported). However, many
        other programs using TIFF (e.g. ImageMagick) have not
        implemented support for those pixel types.  So don't be
        surprised if the generated TIFF is not readable in some
        cases.  If this happens, convert the image to 'unsigned
        char' or 'RGBValue\<unsigned char\>' prior to exporting.
    **/
    ImageExportInfo & setFileType( const char * );
    const char * getFileType() const;

    /** Set compression type.
    
        Recognized strings: "LZW",
        "RunLength", "1" ... "100". A number is interpreted as the
        compression quality for JPEG compression. JPEG compression is
        supported by the JPEG and TIFF formats.
    **/
    ImageExportInfo & setCompression( const char * );
    const char * getCompression() const;

    /** Set the pixel type of the image. Possible values are:
        <DL>
        <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
        <DT>"INT16"<DD> 16-bit signed integer (short)
        <DT>"INT32"<DD> 32-bit signed integer (long)
        <DT>"FLOAT"<DD> 32-bit floating point (float)
        <DT>"DOUBLE"<DD> 64-bit floating point (double)
        </DL>
    **/
    ImageExportInfo & setPixelType( const char * );

    /** Get the pixel type of the image. Possible values are:
        <DL>
        <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
        <DT>"INT16"<DD> 16-bit signed integer (short)
        <DT>"INT32"<DD> 32-bit signed integer (long)
        <DT>"FLOAT"<DD> 32-bit floating point (float)
        <DT>"DOUBLE"<DD> 64-bit floating point (double)
        </DL>
    **/
    const char * getPixelType() const;

        /** Set the image resolution in horizontal direction
         **/
    ImageExportInfo & setXResolution( float );

        /** Set the image resolution in vertical direction
         **/
    ImageExportInfo & setYResolution( float );
    
    float getXResolution() const;
    float getYResolution() const;

  private:
    std::string m_filename, m_filetype, m_pixeltype, m_comp;
    float m_x_res, m_y_res;
};

// return an encoder for a given ImageExportInfo object
std::auto_ptr<Encoder> encoder( const ImageExportInfo & info );


/********************************************************/
/*                                                      */
/*                   ImageImportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function importImage().
See \ref importImage() for a usage example. This object must be
used to read an image from disk and enquire about its properties.

<b>\#include</b> "<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>"<br>
Namespace: vigra
**/
class ImageImportInfo
{
  public:
    enum PixelType { UINT8, INT16, INT32, FLOAT, DOUBLE };

        /** Construct ImageImportInfo object.
        
            The image with the given filename is read into memory.
            The file type will be determined by the first few bytes of the
            file (magic number). Recognized file types:

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format
            (only available if libjpeg is installed).
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"PNG"<DD> Portable Network Graphics
            (only available if libpng is installed).
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format.
            (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>
        **/
    ImageImportInfo( const char *  );
    
    const char * getFileName() const;
    
        /** Get the file type of the image associated with this
            info object.
            
            See ImageImportInfo::ImageImportInfo for a list of the
            available file types.
        **/
    const char * getFileType() const;

        /** Get width of the image.
         **/
    int width() const;

        /** Get height of the image.
         **/
    int height() const;

        /** Get the number bands in the image.
         **/
    int numBands() const;

        /** Get size of the image.
         **/
    Size2D size() const;

        /** Returns true if the image is gray scale.
         **/
    bool isGrayscale() const;

        /** Returns true if the image is colored (RGB).
         **/
    bool isColor() const;

        /** Query the pixel type of the image.
        
            Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
        **/
    const char * getPixelType() const;

        /// deprecated: use getPixelType()
    PixelType pixelType() const;

        /** Returns true if the image has 1 byte per pixel (gray) or
            3 bytes per pixel (RGB).
        **/
    bool isByte() const;

        /** Returns the image resolution in horizontal direction
         **/
    float getXResolution() const;

        /** Returns the image resolution in vertical direction
         **/
    float getYResolution() const;

  private:
    std::string m_filename, m_filetype, m_pixeltype;
    int m_width, m_height, m_num_bands;
    float m_x_res, m_y_res;
};

// return a decoder for a given ImageImportInfo object
std::auto_ptr<vigra::Decoder> decoder( const ImageImportInfo & info );

} // namespace vigra

#endif // VIGRA_IMAGEINFO_HXX
