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

#ifdef HasPNG

#include <stdexcept>
#include <iostream>
#include "void_vector.hxx"
#include "auto_file.hxx"
#include "png.hxx"
#include "error.hxx"

extern "C"
{
#include <png.h>
}

#if PNG_LIBPNG_VER < 10201
#error "please update your libpng to at least 1.2.1"
#endif

// TODO: per-scanline reading/writing

namespace vigra {

    CodecDesc PngCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "PNG";

        // init pixel types
        desc.pixelTypes.resize(2);
        desc.pixelTypes[0] = "UINT8";
        desc.pixelTypes[1] = "INT16";

        // init compression types
        desc.compressionTypes.resize(1);
        desc.compressionTypes[0] = "LOSSLESS";

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = '\x89';
        desc.magicStrings[0][1] = 'P';
        desc.magicStrings[0][2] = 'N';
        desc.magicStrings[0][3] = 'G';

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "png";

        desc.bandNumbers.resize(4);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 2;
        desc.bandNumbers[2] = 3;
        desc.bandNumbers[3] = 4;
        
        return desc;
    }

    std::auto_ptr<Decoder> PngCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new PngDecoder() );
    }

    std::auto_ptr<Encoder> PngCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new PngEncoder() );
    }
    
    namespace {
        std::string png_error_message;
    }

    // called on fatal errors
    static void PngError( png_structp png_ptr, png_const_charp error_msg )
    {
        png_error_message = std::string(error_msg);
        longjmp( png_ptr->jmpbuf, 1 );
    }

    // called on non-fatal errors
    static void PngWarning( png_structp png_ptr, png_const_charp warning_msg )
    {
        std::cerr << warning_msg << std::endl;
    }

    struct PngDecoderImpl
    {
        // data source
        auto_file file;

        // data container
        void_vector_base bands;

        // this is where libpng stores its state
        png_structp png;
        png_infop info;

        // image header fields
        png_uint_32 width, height, components;
        int bit_depth, color_type;

        // scanline counter
        int scanline;

        // ctor, dtor
        PngDecoderImpl( const std::string & filename );
        ~PngDecoderImpl();

        // methods
        void init();
    };

    PngDecoderImpl::PngDecoderImpl( const std::string & filename )
#ifdef _MSC_VER
        : file( filename.c_str(), "rb" ), 
#else
        : file( filename.c_str(), "r" ), 
#endif
          bands(0), scanline(-1)
    {
        png_error_message = "";
        // check if the file is a png file
        const unsigned int sig_size = 8;
        png_byte sig[sig_size];
        std::fread( sig, sig_size, 1, file.get() );
        const int no_png = png_sig_cmp( sig, 0, sig_size );
        vigra_precondition( !no_png, "given file is not a png file.");

        // create png read struct with user defined handlers
        png = png_create_read_struct( PNG_LIBPNG_VER_STRING, NULL,
                                      &PngError, &PngWarning );
        vigra_postcondition( png != 0, "could not create the read struct." );

        // create info struct
        if (setjmp(png->jmpbuf)) {
            png_destroy_read_struct( &png, &info, NULL );
            vigra_postcondition( false, png_error_message.insert(0, "error in png_create_info_struct(): ").c_str() );
        }
        info = png_create_info_struct(png);
        vigra_postcondition( info != 0, "could not create the info struct." );

        // init png i/o
        if (setjmp(png->jmpbuf)) {
            png_destroy_read_struct( &png, &info, NULL );
            vigra_postcondition( false, png_error_message.insert(0, "error in png_init_io(): ").c_str() );
        }
        png_init_io( png, file.get() );

        // specify that the signature was already read
        if (setjmp(png->jmpbuf)) {
            png_destroy_read_struct( &png, &info, NULL );
            vigra_postcondition( false, png_error_message.insert(0, "error in png_set_sig_bytes(): ").c_str() );
        }
        png_set_sig_bytes( png, sig_size );

    }

    PngDecoderImpl::~PngDecoderImpl()
    {
        png_destroy_read_struct( &png, &info, NULL );
    }

    void PngDecoderImpl::init()
    {
        // read all chunks up to the image data
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_read_info(): ").c_str() );
        png_read_info( png, info );

        // pull over the header fields
        int interlace_method, compression_method, filter_method;
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_get_IHDR(): ").c_str() );
        png_get_IHDR( png, info, &width, &height, &bit_depth, &color_type,
                      &interlace_method, &compression_method, &filter_method );
        
        // transform palette to rgb
        if ( color_type == PNG_COLOR_TYPE_PALETTE) {
            if (setjmp(png->jmpbuf))
                vigra_postcondition( false, png_error_message.insert(0, "error in png_palette_to_rgb(): ").c_str() );
            png_set_palette_to_rgb(png);
            color_type = PNG_COLOR_TYPE_RGB;
            bit_depth = 8;
        }

        // expand gray values to at least one byte size
        if ( color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8 ) {
            if (setjmp(png->jmpbuf))
                vigra_postcondition( false,png_error_message.insert(0, "error in png_set_gray_1_2_4_to_8(): ").c_str());
            png_set_gray_1_2_4_to_8(png);
            bit_depth = 8;
        }


#if 0
        // strip alpha channel
        if ( color_type & PNG_COLOR_MASK_ALPHA ) {
            if (setjmp(png->jmpbuf))
                vigra_postcondition( false, png_error_message.insert(0, "error in png_set_strip_alpha(): ").c_str() );
            png_set_strip_alpha(png);
            color_type ^= PNG_COLOR_MASK_ALPHA;
        }
#endif /* #if 0 */


        // find out the number of components
        switch (color_type) {
        case PNG_COLOR_TYPE_GRAY:
            components = 1;
            break;
        case PNG_COLOR_TYPE_GRAY_ALPHA:
            components = 2;
            break;
        case PNG_COLOR_TYPE_RGB:
            components = 3;
            break;
        case PNG_COLOR_TYPE_RGB_ALPHA:
            components = 4;
            break;
        default:
            vigra_fail( "internal error: illegal color type." );
        }

#if 0
        // gamma correction changes the pixels, this is unwanted.

        // image gamma
        double image_gamma = 0.45455;
        if ( png_get_valid( png, info, PNG_INFO_gAMA ) ) {
            if (setjmp(png->jmpbuf))
                vigra_postcondition( false, png_error_message.insert(0, "error in png_get_gAMA(): ").c_str() );
            png_get_gAMA( png, info, &image_gamma );
        }

        // screen gamma
        double screen_gamma = 2.2;

        // set gamma correction
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_set_gamma(): ").c_str() );
        png_set_gamma( png, screen_gamma, image_gamma );
#endif

        // update png library state to reflect any changes that were made
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_read_update_info(): ").c_str() );
        png_read_update_info( png, info );

        const unsigned int size = width * height * components
            * ( bit_depth >> 3 );
        const unsigned int row_stride = size / height;

        // prepare the bands vector
        typedef void_vector< unsigned char > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);
        cbands.resize(size);

        // prepare the row pointers
        void_vector<png_bytep> row_pointers(height);
        for ( unsigned int i = 0; i < height; ++i )
            row_pointers[i] = cbands.data() + row_stride * i;

        // read the whole image
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_read_image(): ").c_str() );
        png_read_image( png, row_pointers.begin() );
    }

    void PngDecoder::init( const std::string & filename )
    {
        pimpl = new PngDecoderImpl(filename);
        pimpl->init();
    }

    PngDecoder::~PngDecoder()
    {
        delete pimpl;
    }

    std::string PngDecoder::getFileType() const
    {
        return "PNG";
    }

    unsigned int PngDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int PngDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int PngDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string PngDecoder::getPixelType() const
    {
        switch (pimpl->bit_depth) {
        case 8:
            return "UINT8";
        case 16:
            return "INT16";
        default:
            vigra_fail( "internal error: illegal pixel type." );
        }
        return "";
    }

    unsigned int PngDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * PngDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        const unsigned int index = pimpl->width * pimpl->components
            * pimpl->scanline + band;
        switch (pimpl->bit_depth) {
        case 8:
            {
                typedef void_vector< unsigned char > bands_type;
                const bands_type & bands
                    = static_cast< const bands_type & >(pimpl->bands);
                return bands.data() + index;
            }
        case 16:
            {
                typedef void_vector<short> bands_type;
                const bands_type & bands
                    = static_cast< const bands_type & >(pimpl->bands);
                return bands.data() + index;
            }
        default:
            vigra_fail( "internal error: illegal bit depth." );
        }
        return 0;
    }

    void PngDecoder::nextScanline()
    {
        ++(pimpl->scanline);
    }

    void PngDecoder::close() {}

    void PngDecoder::abort() {}

    struct PngEncoderImpl
    {
        // data sink
        auto_file file;

        // data container
        void_vector_base bands;

        // this is where libpng stores its state
        png_structp png;
        png_infop info;

        // image header fields
        png_uint_32 width, height, components;
        int bit_depth, color_type;

        // scanline counter
        int scanline;

        // state
        bool finalized;

        // ctor, dtor
        PngEncoderImpl( const std::string & filename );
        ~PngEncoderImpl();

        // methods
        void finalize();
        void write();
    };

    PngEncoderImpl::PngEncoderImpl( const std::string & filename )
#ifdef _MSC_VER
        : file( filename.c_str(), "wb" ), 
#else
        : file( filename.c_str(), "w" ), 
#endif
          bands(0),
          scanline(0), finalized(false)
    {
        png_error_message = "";
        // create png struct with user defined handlers
        png = png_create_write_struct( PNG_LIBPNG_VER_STRING, NULL, 
                                       &PngError, &PngWarning );
        vigra_postcondition( png != 0, "could not create the write struct." );

        // create info struct
        if (setjmp(png->jmpbuf)) {
            png_destroy_write_struct( &png, &info );
            vigra_postcondition( false, png_error_message.insert(0, "error in png_info_struct(): ").c_str() );
        }
        info = png_create_info_struct(png);
        if ( !info ) {
            png_destroy_write_struct( &png, &info );
            vigra_postcondition( false, png_error_message.insert(0, "could not create the info struct.: ").c_str() );
        }

        // init png i/o
        if (setjmp(png->jmpbuf)) {
            png_destroy_write_struct( &png, &info );
            vigra_postcondition( false, png_error_message.insert(0, "error in png_init_io(): ").c_str() );
        }
        png_init_io( png, file.get() );
    }

    PngEncoderImpl::~PngEncoderImpl()
    {
        png_destroy_write_struct( &png, &info );
    }

    void PngEncoderImpl::finalize()
    {
        // write the IHDR
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_set_IHDR(): ").c_str() );
        png_set_IHDR( png, info, width, height, bit_depth, color_type,
                      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                      PNG_FILTER_TYPE_DEFAULT );

        // write the info struct
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_write_info(): ").c_str() );
        png_write_info( png, info );

        // prepare the bands
        bands.resize( ( bit_depth >> 3 ) * width * components * height );

        // enter finalized state
        finalized = true;
    }

    void PngEncoderImpl::write()
    {
        // prepare row pointers
        png_uint_32 row_stride = ( bit_depth >> 3 ) * width * components;
        void_vector<png_byte *>  row_pointers(height);
        typedef void_vector<png_byte> vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);
        png_byte * mover = cbands.data();
        for( png_uint_32 i = 0; i < height; ++i ) {
            row_pointers[i] = mover;
            mover += row_stride;
        }
        // write the whole image
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_write_image(): ").c_str() );
        png_write_image( png, row_pointers.begin() );
        if (setjmp(png->jmpbuf))
            vigra_postcondition( false, png_error_message.insert(0, "error in png_write_end(): ").c_str() );
        png_write_end(png, info);
    }

    void PngEncoder::init( const std::string & filename )
    {
        pimpl = new PngEncoderImpl(filename);
    }

    PngEncoder::~PngEncoder()
    {
        delete pimpl;
    }

    std::string PngEncoder::getFileType() const
    {
        return "PNG";
    }

    void PngEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void PngEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void PngEncoder::setNumBands( unsigned int bands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( bands == 1 )
            pimpl->color_type = PNG_COLOR_TYPE_GRAY;
        else if ( bands == 2 )
            pimpl->color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        else if ( bands == 3 )
            pimpl->color_type = PNG_COLOR_TYPE_RGB;
        else if ( bands == 4 )
            pimpl->color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        else
            vigra_fail( "internal error: number of components not supported." );
        pimpl->components = bands;
    }

    void PngEncoder::setCompressionType( const std::string & comp,
                                         int quality )
    {
        // nothing is settable => do nothing
    }

    void PngEncoder::setPixelType( const std::string & pixelType )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( pixelType == "UINT8" )
            pimpl->bit_depth = 8;
        else if ( pixelType == "INT16" )
            pimpl->bit_depth = 16;
        else
            vigra_fail( "internal error: pixeltype not supported." );
    }

    unsigned int PngEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void PngEncoder::finalizeSettings()
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->finalize();
    }

    void * PngEncoder::currentScanlineOfBand( unsigned int band )
    {
        const unsigned int index = pimpl->width * pimpl->components
            * pimpl->scanline + band;
        switch (pimpl->bit_depth) {
        case 8:
            {
                typedef void_vector< unsigned char > bands_type;
                bands_type & bands
                    = static_cast< bands_type & >(pimpl->bands);
                return bands.data() + index;
            }
        case 16:
            {
                typedef void_vector<short> bands_type;
                bands_type & bands
                    = static_cast< bands_type & >(pimpl->bands);
                return bands.data() + index;
            }
        default:
            vigra_fail( "internal error: illegal bit depth." );
        }
        return 0;
    }

    void PngEncoder::nextScanline()
    {
        ++(pimpl->scanline);
    }

    void PngEncoder::close()
    {
        pimpl->write();
    }

    void PngEncoder::abort() {}
}

#endif // HasPNG
