/************************************************************************/
/*                                                                      */
/*               Copyright 2001-2002 by Pablo d'Angelo                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.4.0, Dec 21 2005 )                                    */
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

#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "hdr.hxx"
#include "auto_file.hxx"
#include "void_vector.hxx"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "rgbe.h"

extern "C"
{
#include "rgbe.h"
}

namespace vigra {

    CodecDesc HDRCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "HDR";

        // init pixel types
        desc.pixelTypes.resize(1);
        desc.pixelTypes[0] = "FLOAT";

        // init compression types
        desc.compressionTypes.resize(1);
        desc.compressionTypes[0] = "NONE";

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = '#';
        desc.magicStrings[0][1] = '?';
        desc.magicStrings[0][2] = 'R';
        desc.magicStrings[0][3] = 'A';

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "hdr";

        desc.bandNumbers.resize(1);
        desc.bandNumbers[0] = 3;

        return desc;
    }

    std::auto_ptr<Decoder> HDRCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new HDRDecoder() );
    }

    std::auto_ptr<Encoder> HDRCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new HDREncoder() );
    }

    class HDRCodecImpl
    {

    protected:

        std::string pixeltype;

        vigra_rgbe_header_info rgbe_h;

        int width, height;
        int samples_per_pixel;

    public:

        HDRCodecImpl();
        ~HDRCodecImpl();
    };

    HDRCodecImpl::HDRCodecImpl()
        : pixeltype("FLOAT")
    {
        samples_per_pixel = 3;
    }

    HDRCodecImpl::~HDRCodecImpl()
    {
    }

    class HDRDecoderImpl : public HDRCodecImpl
    {
        friend class HDRDecoder;

        // data sink
        auto_file infile;
#ifdef DEBUG_HDR
        auto_file dbgFile;
#endif

        // image container
        void_vector<float> scanline;
        int scanline_idx;

    public:

        HDRDecoderImpl( const std::string & filename );
        ~HDRDecoderImpl();

        const void * currentScanlineOfBand( unsigned int band ) const;
        void nextScanline();


    };

    HDRDecoderImpl::HDRDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        // Returns the layer
    : infile( filename.c_str(), "rb" )
#else
    : infile( filename.c_str(), "r" )
#endif
{
        // read width and height
        VIGRA_RGBE_ReadHeader(infile.get() ,&width,&height,&rgbe_h);

        scanline.resize(samples_per_pixel*width);
        scanline_idx = 0;
    }


    HDRDecoderImpl::~HDRDecoderImpl()
    {
    }


    const void *
    HDRDecoderImpl::currentScanlineOfBand( unsigned int band ) const
    {
        return &(scanline[band]);
    }

    void HDRDecoderImpl::nextScanline()
    {
        auto_file & f = const_cast<auto_file &>(infile);
        VIGRA_RGBE_ReadPixels_RLE(f.get(), scanline.data(), width, 1);
#ifdef DEBUG_HDR
        for (int i=0; i < width; i++) {
            fprintf(dbgFile.get(), "%f ", scanline[i]);
        }
        fprintf(dbgFile.get(), "\n");
#endif
    }

    void HDRDecoder::init( const std::string & filename )
    {
        pimpl = new HDRDecoderImpl(filename);
    }

    HDRDecoder::~HDRDecoder()
    {
        delete pimpl;
    }

    std::string HDRDecoder::getFileType() const
    {
        return "HDR";
    }

    unsigned int HDRDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int HDRDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int HDRDecoder::getNumBands() const
    {
        return pimpl->samples_per_pixel;
    }

    std::string HDRDecoder::getPixelType() const
    {
        return pimpl->pixeltype;
    }

    unsigned int HDRDecoder::getOffset() const
    {
        return pimpl->samples_per_pixel;
    }

    const void * HDRDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->currentScanlineOfBand(band);
    }

    void HDRDecoder::nextScanline()
    {
        pimpl->nextScanline();
    }

    void HDRDecoder::close() {}
    void HDRDecoder::abort() {}

    // this encoder always writes interleaved tiff files
    class HDREncoderImpl : public HDRCodecImpl
    {
        friend class HDREncoder;

        // data sink
        auto_file file;

        // image container
        void_vector<float> scanline;
        bool finalized;

    public:
        // ctor, dtor

        HDREncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        // Returns the layer
            : file( filename.c_str(), "wb" ),
#else
            : file( filename.c_str(), "w" ),
#endif
              finalized(false)
        {
        }

        // methods

        void setCompressionType( const std::string &, int );
        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int band )
        {
            return scanline.data() + band;
        }

        void nextScanline()
        {
            // save one scanline
            if (VIGRA_RGBE_WritePixels_RLE(file.get(), scanline.begin(), width, 1) != VIGRA_RGBE_RETURN_SUCCESS)
            {
                vigra_fail("HDREncoder: Could not write scanline");
            }
        }
    };

    void HDREncoderImpl::setCompressionType( const std::string & comp,
                                              int quality = -1 )
    {
    }

    void HDREncoderImpl::finalizeSettings()
    {
        rgbe_h.valid=-1;
        strcpy(rgbe_h.programtype,"RADIANCE");
        rgbe_h.gamma=1.0;
        rgbe_h.exposure=1.0;

        scanline.resize(samples_per_pixel*width);

        if (VIGRA_RGBE_WriteHeader(file.get(), width, height, &rgbe_h) != VIGRA_RGBE_RETURN_SUCCESS ) {
            vigra_fail("HDREncoder: Could not write header");
        }
        finalized = true;
    }

    void HDREncoder::init( const std::string & filename )
    {
        pimpl = new HDREncoderImpl(filename);
    }

    HDREncoder::~HDREncoder()
    {
        delete pimpl;
    }

    std::string HDREncoder::getFileType() const
    {
        return "HDR";
    }

    unsigned int HDREncoder::getOffset() const
    {
        return pimpl->samples_per_pixel;
    }

    void HDREncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void HDREncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void HDREncoder::setNumBands( unsigned int bands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if (bands != 3) {
            vigra_fail("HDREncoder: can only save 3 channel images");
        }
        pimpl->samples_per_pixel = 3;
    }

    void HDREncoder::setCompressionType( const std::string & comp,
                                          int quality )
    {
    }

    void HDREncoder::setPixelType( const std::string & pixeltype )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( pixeltype != "FLOAT" )
            vigra_fail( "internal error: pixeltype not supported." );
        pimpl->pixeltype = "FLOAT";
    }

    void HDREncoder::setPosition( const vigra::Diff2D & pos )
    {
    }

    void HDREncoder::setXResolution( float xres )
    {
    }

    void HDREncoder::setYResolution( float yres )
    {
    }

    void HDREncoder::finalizeSettings()
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->finalizeSettings();
    }

    void * HDREncoder::currentScanlineOfBand( unsigned int band )
    {
        return pimpl->currentScanlineOfBand(band);
    }

    void HDREncoder::nextScanline()
    {
        pimpl->nextScanline();
    }

    void HDREncoder::close() {}
    void HDREncoder::abort() {}
}

