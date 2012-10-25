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
 * as of 4 March 2006:
 *  - Added support for x and y resolution fields.
 *  - Added ICC support.
 */

#ifdef HasJPEG

#include <stdexcept>
#include <csetjmp>
#include "vigra/config.hxx"
#include "void_vector.hxx"
#include "error.hxx"
#include "auto_file.hxx"
#include "jpeg.hxx"

extern "C" {

#include <jpeglib.h>
#include "iccjpeg.h"

} // extern "C"

namespace {

struct JPEGCodecErrorManager
{
    jpeg_error_mgr pub;
    std::jmp_buf buf;
};

} // namespace

extern "C"
{

static void JPEGCodecLongjumper( j_common_ptr info )
{
    (*info->err->output_message)(info);
    JPEGCodecErrorManager * error = reinterpret_cast< JPEGCodecErrorManager * >(info->err);
    std::longjmp( error->buf, 1 );
}

} // extern "C"

namespace vigra
{
    CodecDesc JPEGCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "JPEG";

        // init pixel types
        desc.pixelTypes.resize(1);
        desc.pixelTypes[0] = "UINT8";

        // init compression types
        desc.compressionTypes.resize(1);
        desc.compressionTypes[0] = "JPEG";

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(3);
        desc.magicStrings[0][0] = '\377';
        desc.magicStrings[0][1] = '\330';
        desc.magicStrings[0][2] = '\377';

        // init file extensions
        desc.fileExtensions.resize(2);
        desc.fileExtensions[0] = "jpg";
        desc.fileExtensions[1] = "jpeg";

        desc.bandNumbers.resize(2);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 3;

        return desc;
    }

    VIGRA_UNIQUE_PTR<Decoder> JPEGCodecFactory::getDecoder() const
    {
        return VIGRA_UNIQUE_PTR<Decoder>( new JPEGDecoder() );
    }

    VIGRA_UNIQUE_PTR<Encoder> JPEGCodecFactory::getEncoder() const
    {
        return VIGRA_UNIQUE_PTR<Encoder>( new JPEGEncoder() );
    }

    class JPEGCodecImpl
    {
        // extend the jpeg_error_mgr by a jump buffer

    public:

        // attributes

        JPEGCodecErrorManager err;

    };

    class JPEGDecoderImplBase : public JPEGCodecImpl
    {
    public:

        // attributes
        jpeg_decompress_struct info;

        JPEGDecoderImplBase()
        {
            // create the decompression struct
            jpeg_create_decompress(&info);
        }

        virtual ~JPEGDecoderImplBase()
        {
            // delete the decompression struct
            jpeg_destroy_decompress(&info);
        }
    };

    class JPEGEncoderImplBase : public JPEGCodecImpl
    {
    public:

        // attributes
        jpeg_compress_struct info;

        JPEGEncoderImplBase()
        {
            // create the decompression struct
            jpeg_create_compress(&info);
        }

        virtual ~JPEGEncoderImplBase()
        {
            // delete the decompression struct
            jpeg_destroy_compress(&info);
        }
    };

    struct JPEGDecoderImpl : public JPEGDecoderImplBase
    {
        // attributes

        auto_file file;
        void_vector<JSAMPLE> bands;
        unsigned int width, height, components, scanline;

        // icc profile, if available
        UInt32 iccProfileLength;
        const unsigned char *iccProfilePtr;

        // ctor, dtor
        JPEGDecoderImpl( const std::string & filename );
        ~JPEGDecoderImpl();

        // methods

        void init();
    };

    JPEGDecoderImpl::JPEGDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : file( filename.c_str(), "rb" ),
#else
        : file( filename.c_str(), "r" ),
#endif
          bands(0), scanline(0), iccProfileLength(0), iccProfilePtr(NULL)
    {
        // setup setjmp() error handling
        info.err = jpeg_std_error( ( jpeg_error_mgr * ) &err );
        err.pub.error_exit = &JPEGCodecLongjumper;

        // setup the data source
        if (setjmp(err.buf)) {
            vigra_fail( "error in jpeg_stdio_src()" );
        }
        jpeg_stdio_src( &info, file.get() );
        // prepare for icc profile
        setup_read_icc_profile(&info);
    }

    void JPEGDecoderImpl::init()
    {
        // read the header
        if (setjmp(err.buf))
            vigra_fail( "error in jpeg_read_header()" );
        jpeg_read_header( &info, TRUE );

        // extract ICC profile
        JOCTET *iccBuf;
        unsigned int iccLen;
        if (read_icc_profile(&info, &iccBuf, &iccLen)) {
            iccProfileLength = iccLen;
            iccProfilePtr = iccBuf;
        }

        // start the decompression
        if (setjmp(err.buf))
            vigra_fail( "error in jpeg_start_decompress()" );
        jpeg_start_decompress(&info);

        // transfer interesting header information
        width = info.output_width;
        height = info.output_height;
        components = info.output_components;

        // alloc memory for a single scanline
        bands.resize( width * components );

        // set colorspace
        info.jpeg_color_space = components == 1 ? JCS_GRAYSCALE : JCS_RGB;
    }

    JPEGDecoderImpl::~JPEGDecoderImpl()
    {
        if (iccProfilePtr && iccProfileLength)
            free((void *)iccProfilePtr);
    }

    void JPEGDecoder::init( const std::string & filename )
    {
        pimpl = new JPEGDecoderImpl(filename);
        pimpl->init();
        if(pimpl->iccProfileLength)
        {
            Decoder::ICCProfile iccData(
                pimpl->iccProfilePtr,
                pimpl->iccProfilePtr + pimpl->iccProfileLength);
            iccProfile_.swap(iccData);
        }
    }

    JPEGDecoder::~JPEGDecoder()
    {
        delete pimpl;
    }

    std::string JPEGDecoder::getFileType() const
    {
        return "JPEG";
    }

    unsigned int JPEGDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int JPEGDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int JPEGDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string JPEGDecoder::getPixelType() const
    {
        return "UINT8";
    }

    unsigned int JPEGDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * JPEGDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->bands.data() + band;
    }

    void JPEGDecoder::nextScanline()
    {
        // check if there are scanlines left at all, eventually read one
        JSAMPLE * band = pimpl->bands.data();
        if ( pimpl->info.output_scanline < pimpl->info.output_height ) {
            if (setjmp(pimpl->err.buf))
                vigra_fail( "error in jpeg_read_scanlines()" );
            jpeg_read_scanlines( &pimpl->info, &band, 1 );
        }
    }

    void JPEGDecoder::close()
    {
        // finish any pending decompression
        if (setjmp(pimpl->err.buf))
            vigra_fail( "error in jpeg_finish_decompress()" );
        jpeg_finish_decompress(&pimpl->info);
    }

    void JPEGDecoder::abort() {}

    struct JPEGEncoderImpl : public JPEGEncoderImplBase
    {
        // attributes

        auto_file file;
        void_vector<JSAMPLE> bands;
        unsigned int width, height, components, scanline;
        int quality;

        // icc profile, if available
        Encoder::ICCProfile iccProfile;

        // state
        bool finalized;

        // ctor, dtor

        JPEGEncoderImpl( const std::string & filename );
        ~JPEGEncoderImpl();

        // methods

        void finalize();
    };

    JPEGEncoderImpl::JPEGEncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : file( filename.c_str(), "wb" ),
#else
        : file( filename.c_str(), "w" ),
#endif
          scanline(0), quality(-1), finalized(false)
    {
        // setup setjmp() error handling
        info.err = jpeg_std_error( ( jpeg_error_mgr * ) &err );
        err.pub.error_exit = &JPEGCodecLongjumper;

        // setup the data dest
        if (setjmp(err.buf)) {
            vigra_fail( "error in jpeg_stdio_dest()" );
        }
        jpeg_stdio_dest( &info, file.get() );
    }

    JPEGEncoderImpl::~JPEGEncoderImpl()
    {
    }

    void JPEGEncoderImpl::finalize()
    {
        VIGRA_IMPEX_FINALIZED(finalized);

        // alloc memory for a single scanline
        bands.resize( width * components );
        finalized = true;

        // init the compression
        info.image_width = width;
        info.image_height = height;
        info.input_components = components;

        // rgb or gray can be assumed here
        info.in_color_space = components == 1 ? JCS_GRAYSCALE : JCS_RGB;
        info.X_density = 100;
        info.Y_density = 100;

        // set defaults based upon the set values
        if (setjmp(err.buf))
            vigra_fail( "error in jpeg_set_defaults()" );
        jpeg_set_defaults(&info);

        // set the quality level
        if ( quality != -1 ) {
            if (setjmp(err.buf))
                vigra_fail( "error in jpeg_set_quality()" );
            jpeg_set_quality( &info, quality, TRUE );
        }

        // enhance the quality a little bit
        for ( unsigned int i = 0; i < MAX_COMPONENTS; ++i ) {
            info.comp_info[i].h_samp_factor = 1;
            info.comp_info[i].v_samp_factor = 1;
        }
#ifdef ENTROPY_OPT_SUPPORTED
        info.optimize_coding = TRUE;
#endif
        info.dct_method = JDCT_FLOAT;

        // start the compression
        if (setjmp(err.buf))
            vigra_fail( "error in jpeg_start_compress()" );
        jpeg_start_compress( &info, TRUE );

        if (iccProfile.size()) {
            write_icc_profile(&info, iccProfile.begin(), (unsigned int)iccProfile.size());
        }
    }

    void JPEGEncoder::init( const std::string & filename )
    {
        pimpl = new JPEGEncoderImpl(filename);
    }

    JPEGEncoder::~JPEGEncoder()
    {
        delete pimpl;
    }

    std::string JPEGEncoder::getFileType() const
    {
        return "JPEG";
    }

    void JPEGEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void JPEGEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void JPEGEncoder::setNumBands( unsigned int bands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->components = bands;
    }

    void JPEGEncoder::setCompressionType( const std::string & comp,
                                          int quality )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( comp == "LOSSLESS" )
            vigra_fail( "lossless encoding is not supported by your jpeg library." );
        if ( comp == "JPEG_ARITH" )
#ifdef C_ARITH_CODING_SUPPORTED
            pimpl->info.arith_code = TRUE;
#else
            vigra_fail( "arithmetic encoding is not supported by your jpeg library." );
#endif
        pimpl->quality = quality;
    }

    void JPEGEncoder::setPixelType( const std::string & pixelType )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( pixelType != "UINT8" )
            vigra_precondition( false, "only UINT8 pixels are supported." );
    }

    unsigned int JPEGEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void JPEGEncoder::setICCProfile(const ICCProfile & data)
    {
        pimpl->iccProfile = data;
    }

    void JPEGEncoder::finalizeSettings()
    {
        pimpl->finalize();
    }

    void * JPEGEncoder::currentScanlineOfBand( unsigned int band )
    {
        return pimpl->bands.data() + band;
    }

    void JPEGEncoder::nextScanline()
    {
        // check if there are scanlines left at all, eventually write one
        JSAMPLE * band = pimpl->bands.data();
        if ( pimpl->info.next_scanline < pimpl->info.image_height ) {
            if (setjmp(pimpl->err.buf))
                vigra_fail( "error in jpeg_write_scanlines()" );
            jpeg_write_scanlines( &pimpl->info, &band, 1 );
        }
    }

    void JPEGEncoder::close()
    {
        // finish any pending compression
        if (setjmp(pimpl->err.buf))
            vigra_fail( "error in jpeg_finish_compress()" );
        jpeg_finish_compress( &pimpl->info );
    }

    void JPEGEncoder::abort() {}
}

#endif // HasJPEG
