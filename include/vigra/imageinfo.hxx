/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*               Copyright 2001-2002 by Gunnar Kedenburg                */
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

/* Modifications by Pablo d'Angelo
 * updated to vigra 1.4 by Douglas Wilkins
 * as of 18 Febuary 2006:
 *  - Added UINT16 and UINT32 pixel types.
 *  - Added support for obtaining extra bands beyond RGB.
 *  - Added support for a position field that indicates the start of this
 *    image relative to some global origin.
 *  - Added support for x and y resolution fields.
 *  - Added support for ICC profiles
 */

#ifndef VIGRA_IMAGEINFO_HXX
#define VIGRA_IMAGEINFO_HXX

#include <memory>
#include <string>
#include "config.hxx"
#include "error.hxx"
#include "diff2d.hxx"
#include "codec.hxx"
#include "array_vector.hxx"
#include "multi_iterator.hxx"

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

        <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
        Namespace: vigra

        \code
        std::cout << "supported formats: " << vigra::impexListFormats() << std::endl;
        \endcode

    **/
VIGRA_EXPORT std::string impexListFormats();

    /** \brief List the file extension VIGRA understands.

        This is useful for creating file dialogs that only list image files
        VIGRA can actually import.

        <b> Usage:</b>

        <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
        Namespace: vigra

        \code
        std::cout << "supported extensions: " << vigra::impexListExtensions() << std::endl;
        \endcode

    **/
VIGRA_EXPORT std::string impexListExtensions();

/** \brief Test whether a file is an image format known to VIGRA.

    This checks the first few bytes of the file and compares them with the
    "magic strings" of each recognized image format.

    <b> Usage:</b>

    <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
    Namespace: vigra

    \code
    std::cout << "is image: " << vigra::isImage("foo.bmp") << std::endl;
    \endcode

**/
VIGRA_EXPORT bool isImage(char const * filename);

/********************************************************/
/*                                                      */
/*                   ImageExportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function exportImage().

    See \ref exportImage() for usage example. This object must be used
    to define the properties of an image to be written to disk.

    <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
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
            '.tif', '.tiff', '.xv', '.hdr'.
            JPEG support requires libjpeg, PNG support requires libpng, and
            TIFF support requires libtiff.
         **/
    VIGRA_EXPORT ImageExportInfo( const char * );
    VIGRA_EXPORT ~ImageExportInfo();

        /** Set image file name.
        
            The file type will be guessed from the extension unless overridden
            by \ref setFileType(). Recognized extensions: '.bmp', '.gif',
            '.jpeg', '.jpg', '.p7', '.png', '.pbm', '.pgm', '.pnm', '.ppm', '.ras',
            '.tif', '.tiff', '.xv', '.hdr'.
            JPEG support requires libjpeg, PNG support requires libpng, and
            TIFF support requires libtiff.
         **/
    VIGRA_EXPORT ImageExportInfo & setFileName(const char * filename);
    VIGRA_EXPORT const char * getFileName() const;

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

            With the exception of TIFF, VIFF, PNG, and PNM all file types store
            1 byte (gray scale and mapped RGB) or 3 bytes (RGB) per
            pixel.

            PNG can store UInt8 and UInt16 values, and supports 1 and 3 channel
            images. One additional alpha channel is also supported.

            PNM can store 1 and 3 channel images with UInt8, UInt16 and UInt32
            values in each channel.

            TIFF and VIFF are aditionally able to store short and long
            integers (2 or 4 bytes) and real values (32 bit float and
            64 bit double) without conversion. So you will need to use
            TIFF or VIFF if you need to store images with high
            accuracy (the appropriate type to write is automatically
            derived from the image type to be exported). However, many
            other programs using TIFF (e.g. ImageMagick) have not
            implemented support for those pixel types.  So don't be
            surprised if the generated TIFF is not readable in some
            cases.  If this happens, export the image as 'unsigned
            char' or 'RGBValue\<unsigned char\>' by calling
            \ref ImageExportInfo::setPixelType().

            Support to reading and writing ICC color profiles is
            provided for TIFF, JPEG, and PNG images.
         **/
    VIGRA_EXPORT ImageExportInfo & setFileType( const char * );
    VIGRA_EXPORT const char * getFileType() const;

        /** Set compression type.

            Recognized strings: "" (no compression), "LZW",
            "RunLength", "1" ... "100". A number is interpreted as the
            compression quality for JPEG compression. JPEG compression is
            supported by the JPEG and TIFF formats. "LZW" is only available
            if libtiff was installed with LZW enabled. By default, libtiff came
            with LZW disabled due to Unisys patent enforcement. In this case,
            VIGRA stores the image uncompressed.

                Valid Compression for TIFF files:
                  JPEG    jpeg compression, call setQuality as well!
                  RLE     runlength compression
                  LZW     lzw compression
                  DEFLATE deflate compression
         **/
    VIGRA_EXPORT ImageExportInfo & setCompression( const char * );
    VIGRA_EXPORT const char * getCompression() const;

        /** Set the pixel type of the image file. Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>

            <b>Usage:</b>
            
            \code
            FImage img(w,h);

            // by default, float images are exported with pixeltype float
            // when the target format support this type, i.e. is TIFF or VIFF.
            exportImage(srcImageRange(img), ImageExportInfo("asFloat.tif"));

            // if this is not desired, force a different pixeltype
            exportImage(srcImageRange(img), ImageExportInfo("asByte.tif").setPixelType("UINT8"));
            \endcode
         **/
    VIGRA_EXPORT ImageExportInfo & setPixelType( const char * );

        /** Get the pixel type of the image. Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT const char * getPixelType() const;
    
    VIGRA_EXPORT ImageExportInfo & setForcedRangeMapping(double fromMin, double fromMax,
                                                     double toMin, double toMax);    
    VIGRA_EXPORT bool hasForcedRangeMapping() const;
    VIGRA_EXPORT double getFromMin() const;
    VIGRA_EXPORT double getFromMax() const;
    VIGRA_EXPORT double getToMin() const;
    VIGRA_EXPORT double getToMax() const;
    
        /** Set the image resolution in horizontal direction
         **/
    VIGRA_EXPORT ImageExportInfo & setXResolution( float );
    VIGRA_EXPORT float getXResolution() const;

        /** Set the image resolution in vertical direction
         **/
    VIGRA_EXPORT ImageExportInfo & setYResolution( float );
    VIGRA_EXPORT float getYResolution() const;

        /** Set the position of the upper Left corner on a global
            canvas.

            Currently only supported by TIFF and PNG files.

            The offset is encoded in the XPosition and YPosition TIFF tags.

            @param pos     position of the upper left corner in pixels
                           (must be >= 0)
         **/
    VIGRA_EXPORT ImageExportInfo & setPosition(const Diff2D & pos);

        /** Get the position of the upper left corner on
            a global canvas.
         **/
    VIGRA_EXPORT Diff2D getPosition() const;

        /**
          ICC profiles (handled as raw data so far).
          see getICCProfile()/setICCProfile()
         **/
    typedef ArrayVector<unsigned char> ICCProfile;

        /** Returns a reference to the ICC profile.
         */
    VIGRA_EXPORT const ICCProfile & getICCProfile() const;

        /** Sets the ICC profile.
            ICC profiles are currently supported by TIFF, PNG and JPEG images.
            (Otherwise, the profile data is silently ignored.)
         **/
    VIGRA_EXPORT ImageExportInfo & setICCProfile(const ICCProfile & profile);

  private:
    std::string m_filename, m_filetype, m_pixeltype, m_comp;
    float m_x_res, m_y_res;
    Diff2D m_pos;
    ICCProfile m_icc_profile;
    double fromMin_, fromMax_, toMin_, toMax_;
};

// return an encoder for a given ImageExportInfo object
VIGRA_EXPORT std::auto_ptr<Encoder> encoder( const ImageExportInfo & info );

/********************************************************/
/*                                                      */
/*                   ImageImportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function importImage().

See \ref importImage() for a usage example. This object must be
used to read an image from disk and enquire about its properties.

<b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
Namespace: vigra
**/
class ImageImportInfo
{
  public:
    enum PixelType { UINT8, INT16, UINT16, INT32, UINT32, FLOAT, DOUBLE };

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
    VIGRA_EXPORT ImageImportInfo( const char *  );
    VIGRA_EXPORT ~ImageImportInfo();

    VIGRA_EXPORT const char * getFileName() const;

        /** Get the file type of the image associated with this
            info object.

            See ImageImportInfo::ImageImportInfo for a list of the
            available file types.
         **/
    VIGRA_EXPORT const char * getFileType() const;

        /** Get width of the image.
         **/
    VIGRA_EXPORT int width() const;

        /** Get height of the image.
         **/
    VIGRA_EXPORT int height() const;

        /** Get the total number of bands in the image.
         **/
    VIGRA_EXPORT int numBands() const;

        /** Get the number of extra (non color) bands in the image.
         ** Usually these are the alpha channels.
         **/
    VIGRA_EXPORT int numExtraBands() const;

        /** Get size of the image.
         **/
    VIGRA_EXPORT Size2D size() const;

        /** Get size of the image in a form compatible to MultiArray.
         **/
	VIGRA_EXPORT MultiArrayShape<2>::type shape() const;

        /** Returns true if the image is gray scale.
         **/
    VIGRA_EXPORT bool isGrayscale() const;

        /** Returns true if the image is colored (RGB).
         **/
    VIGRA_EXPORT bool isColor() const;

        /** Query the pixel type of the image.

            Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT const char * getPixelType() const;

        /** Query the pixel type of the image.

            Same as getPixelType(), but the result is returned as a 
            ImageImportInfo::PixelType enum. This is useful to implement
            a switch() on the pixel type.

            Possible values are:
            <DL>
            <DT>UINT8<DD> 8-bit unsigned integer (unsigned char)
            <DT>INT16<DD> 16-bit signed integer (short)
            <DT>UINT16<DD> 16-bit unsigned integer (unsigned short)
            <DT>INT32<DD> 32-bit signed integer (long)
            <DT>UINT32<DD> 32-bit unsigned integer (unsigned long)
            <DT>FLOAT<DD> 32-bit floating point (float)
            <DT>DOUBLE<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT PixelType pixelType() const;

        /** Returns true if the image has 1 byte per pixel (gray) or
            3 bytes per pixel (RGB).
         **/
    VIGRA_EXPORT bool isByte() const;

        /** Returns the layer offset of the current image, if there is one
         **/
    VIGRA_EXPORT Diff2D getPosition() const;

        /** Returns the image resolution in horizontal direction
         **/
    VIGRA_EXPORT float getXResolution() const;

        /** Returns the image resolution in vertical direction
         **/
    VIGRA_EXPORT float getYResolution() const;

        /**
          ICC profiles (handled as raw data so far).
          see getICCProfile()/setICCProfile()
         **/
    typedef ArrayVector<unsigned char> ICCProfile;

        /** Returns a reference to the ICC profile.

           Note: The reference will become invalid when the
           ImageImportInfo object has been destroyed.
         **/
    VIGRA_EXPORT const ICCProfile & getICCProfile() const;

  private:
    std::string m_filename, m_filetype, m_pixeltype;
    int m_width, m_height, m_num_bands, m_num_extra_bands;
    float m_x_res, m_y_res;
    Diff2D m_pos;
    ICCProfile m_icc_profile;
};

// return a decoder for a given ImageImportInfo object
VIGRA_EXPORT std::auto_ptr<Decoder> decoder( const ImageImportInfo & info );

//@}

} // namespace vigra

#endif // VIGRA_IMAGEINFO_HXX
