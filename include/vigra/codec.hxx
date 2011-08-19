/************************************************************************/
/*                                                                      */
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
 *  - Added support for ICC Profiles
 */

#ifndef VIGRA_CODEC_HXX
#define VIGRA_CODEC_HXX

#include <memory>
#include <string>
#include <vector>

#include "array_vector.hxx"
#include "config.hxx"
#include "diff2d.hxx"
#include "sized_int.hxx"

// possible pixel types:
// "undefined", "UINT8", "UINT16", "INT16", "UINT32", "INT32", "FLOAT", "DOUBLE"

// possible compression types:
// "undefined", "RLE", "LZW", "LOSSLESS", "JPEG", "DEFLATE"

// possible file types:
// "undefined", "TIFF", "VIFF", "JPEG", "PNG", "PNM", "BMP", "SUN", "XPM", "EXR"

// possible name extensions:
// "undefined", "tif", "tiff", "jpg", "jpeg", "png", "pnm", "bmp", "sun",
// "xpm", "exr" (also capital forms)

namespace vigra
{
    template <class T>
    struct TypeAsString
    {
        static std::string result() { return "undefined"; }
    };

    template <>
    struct TypeAsString<Int8>
    {
        static std::string result() { return "INT8"; }
    };

    template <>
    struct TypeAsString<UInt8>
    {
        static std::string result() { return "UINT8"; }
    };

    template <>
    struct TypeAsString<Int16>
    {
        static std::string result() { return "INT16"; }
    };

    template <>
    struct TypeAsString<UInt16>
    {
        static std::string result() { return "UINT16"; }
    };

    template <>
    struct TypeAsString<Int32>
    {
        static std::string result() { return "INT32"; }
    };

    template <>
    struct TypeAsString<UInt32>
    {
        static std::string result() { return "UINT32"; }
    };

    template <>
    struct TypeAsString<float>
    {
        static std::string result() { return "FLOAT"; }
    };

    template <>
    struct TypeAsString<double>
    {
        static std::string result() { return "DOUBLE"; }
    };


    // codec description
    struct CodecDesc
    {
        std::string fileType;
        std::vector<std::string> pixelTypes;
        std::vector<std::string> compressionTypes;
        std::vector<std::vector<char> > magicStrings;
        std::vector<std::string> fileExtensions;
        std::vector<int> bandNumbers;
    };

    // Decoder and Encoder are virtual types that define a common
    // interface for all image file formats impex supports.

    struct Decoder
    {
        virtual ~Decoder() {};
        virtual void init( const std::string & ) = 0;

        // initialize with an image index. For codecs that do not support this feature, the standard init is called.
        virtual void init( const std::string & fileName, unsigned int imageIndex)
        {
          init(fileName);
        }

        virtual void close() = 0;
        virtual void abort() = 0;

        virtual std::string getFileType() const = 0;
        virtual std::string getPixelType() const = 0;

        virtual unsigned int getNumImages() const
        {
          return 1;
        }

        virtual void setImageIndex(unsigned int)
        {
        }

        virtual unsigned int getImageIndex() const
        {
          return 0;
        }

        virtual unsigned int getWidth() const = 0;
        virtual unsigned int getHeight() const = 0;
        virtual unsigned int getNumBands() const = 0;
        virtual unsigned int getNumExtraBands() const
        {
            return 0;
        }

        virtual vigra::Diff2D getPosition() const
        {
            return vigra::Diff2D();
        }

        virtual float getXResolution() const
        {
            return 0.0f;
        }
        virtual float getYResolution() const
        {
            return 0.0f;
        }

        virtual vigra::Size2D getCanvasSize() const
        {
            return vigra::Size2D(this->getWidth(), this->getHeight());
        }

        virtual unsigned int getOffset() const = 0;

        virtual const void * currentScanlineOfBand( unsigned int ) const = 0;
        virtual void nextScanline() = 0;

        typedef ArrayVector<unsigned char> ICCProfile;

        const ICCProfile & getICCProfile() const
        {
            return iccProfile_;
        }

        ICCProfile iccProfile_;
    };

    struct Encoder
    {
        virtual ~Encoder() {};
        virtual void init( const std::string & ) = 0;

        // initialize with file access mode. For codecs that do not support this feature, the standard init is called.
        virtual void init( const std::string & fileName, const std::string & mode )
        {
          init(fileName);
        }

        virtual void close() = 0;
        virtual void abort() = 0;

        virtual std::string getFileType() const = 0;
        virtual unsigned int getOffset() const = 0;

        virtual void setWidth( unsigned int ) = 0;
        virtual void setHeight( unsigned int ) = 0;
        virtual void setNumBands( unsigned int ) = 0;
        virtual void setCompressionType( const std::string &, int = -1 ) = 0;
        virtual void setPixelType( const std::string & ) = 0;
        virtual void finalizeSettings() = 0;

        virtual void setPosition( const vigra::Diff2D & /*pos*/ )
        {
        }
        virtual void setCanvasSize( const vigra::Size2D & size)
        {
        }
        virtual void setXResolution( float /*xres*/ )
        {
        }
        virtual void setYResolution( float /*yres*/ )
        {
        }

        typedef ArrayVector<unsigned char> ICCProfile;

        virtual void setICCProfile(const ICCProfile & /* data */)
        {
        }

        virtual void * currentScanlineOfBand( unsigned int ) = 0;
        virtual void nextScanline() = 0;

        struct TIFFCompressionException {};
    };

    // codec factory for registration at the codec manager

    struct CodecFactory
    {
        virtual CodecDesc getCodecDesc() const = 0;
        virtual std::auto_ptr<Decoder> getDecoder() const = 0;
        virtual std::auto_ptr<Encoder> getEncoder() const = 0;
        virtual ~CodecFactory() {};
    };

    // factory functions to encapsulate the codec managers
    //
    // codecs are selected according to the following order:
    // - (if provided) the FileType
    // - (in case of decoders) the file's magic string
    // - the filename extension

    VIGRA_EXPORT std::auto_ptr<Decoder>
    getDecoder( const std::string &, const std::string & = "undefined", unsigned int = 0 );

    VIGRA_EXPORT std::auto_ptr<Encoder>
    getEncoder( const std::string &, const std::string & = "undefined", const std::string & = "w" );

    VIGRA_EXPORT std::string
    getEncoderType( const std::string &, const std::string & = "undefined" );

    // functions to query the capabilities of certain codecs

    VIGRA_EXPORT std::vector<std::string> queryCodecPixelTypes( const std::string & );

    VIGRA_EXPORT bool negotiatePixelType( std::string const & codecname,
                 std::string const & srcPixeltype, std::string & destPixeltype);

    VIGRA_EXPORT bool isPixelTypeSupported( const std::string &, const std::string & );

    VIGRA_EXPORT bool isBandNumberSupported( const std::string &, int bands );
}

#endif // VIGRA_CODEC_HXX
