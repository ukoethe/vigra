/************************************************************************/
/*                                                                      */
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

#ifndef VIGRA_CODEC_HXX
#define VIGRA_CODEC_HXX

#include <memory>
#include <string>
#include <vector>

// possible pixel types:
// "undefined", "UINT8", "INT16", "INT32", "FLOAT", "DOUBLE"

// possible compression types:
// "undefined", "RLE", "LZW", "LOSSLESS", "JPEG"

// possible file types:
// "undefined", "TIFF", "VIFF", "JPEG", "PNG", "PNM", "BMP", "SUN", "XPM"

// possible name extensions:
// "undefined", "tif", "tiff", "jpg", "jpeg", "png", "pnm", "bmp", "sun",
// "xpm" (also capital forms)

namespace vigra
{
    // codec description

    struct CodecDesc
    {
        std::string fileType;
        std::vector<std::string> pixelTypes;
        std::vector<std::string> compressionTypes;
        std::vector<std::vector<char> > magicStrings;
        std::vector<std::string> fileExtensions;
    };

    // Decoder and Encoder are pure virtual types that define a common
    // interface for all image file formats impex supports.

    struct Decoder
    {
        virtual ~Decoder() {};
        virtual void init( const std::string & ) = 0;
        virtual void close() = 0;
        virtual void abort() = 0;

        virtual std::string getFileType() const = 0;
        virtual std::string getPixelType() const = 0;

        virtual unsigned int getWidth() const = 0;
        virtual unsigned int getHeight() const = 0;
        virtual unsigned int getNumBands() const = 0;
        virtual unsigned int getOffset() const = 0;

        virtual const void * currentScanlineOfBand( unsigned int ) const = 0;
        virtual void nextScanline() = 0;
    };

    struct Encoder
    {
        virtual ~Encoder() {};
        virtual void init( const std::string & ) = 0;
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

        virtual void * currentScanlineOfBand( unsigned int ) = 0;
        virtual void nextScanline() = 0;
    };

    // codec factory for registration at the codec manager

    struct CodecFactory
    {
        virtual CodecDesc getCodecDesc() const = 0;
        virtual std::auto_ptr<Decoder> getDecoder() const = 0;
        virtual std::auto_ptr<Encoder> getEncoder() const = 0;
    };

    // factory functions to encapsulate the codec managers
    //
    // codecs are selected according to the following order:
    // - (if provided) the FileType
    // - (in case of decoders) the file's magic string
    // - the filename extension

    std::auto_ptr<Decoder>
    getDecoder( const std::string &, const std::string & = "undefined" );

    std::auto_ptr<Encoder>
    getEncoder( const std::string &, const std::string & = "undefined" );

    // functions to query the capabilities of certain codecs

    std::vector<std::string> queryCodecPixelTypes( const std::string & );

    bool isPixelTypeSupported( const std::string &, const std::string & );
}

#endif // VIGRA_CODEC_HXX
