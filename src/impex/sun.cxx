/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Gunnar Kedenburg                     */
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

#include <fstream>
#include <stdexcept>
#include "error.hxx"
#include "byteorder.hxx"
#include "void_vector.hxx"
#include "sun.hxx"

// magic number
#define RAS_MAGIC         0x59A66A95
#define RAS_MAGIC_REVERSE 0x956AA659

// raster map type
#define RMT_NONE       0
#define RMT_EQUAL_RGB  1
#define RMT_RAW        2

// raster type
#define RT_OLD         0
#define RT_STANDARD    1
#define RT_ENCODED     2
#define RT_FORMAT_RGB  3

namespace vigra {

    CodecDesc SunCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "SUN";

        // init pixel types
        desc.pixelTypes.resize(1);
        desc.pixelTypes[0] = "UINT8";

        // init compression types
        desc.compressionTypes.resize(0);

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = '\x59';
        desc.magicStrings[0][1] = '\xA6';
        desc.magicStrings[0][2] = '\x6A';
        desc.magicStrings[0][3] = '\x95';

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "ras";

        return desc;
    }

    std::auto_ptr<Decoder> SunCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new SunDecoder() );
    }

    std::auto_ptr<Encoder> SunCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new SunEncoder() );
    }

    struct SunHeader
    {
        typedef unsigned long field_type;

        // attributes

        field_type width, height, depth, length, type, maptype, maplength;

        // methods

        void from_stream( std::ifstream & stream, const byteorder & bo );
        void to_stream( std::ofstream & stream, const byteorder & bo );
    };

    void SunHeader::from_stream( std::ifstream & stream, const byteorder & bo )
    {
        read_field( stream, bo, width );
        read_field( stream, bo, height );
        read_field( stream, bo, depth );
        read_field( stream, bo, length );
        read_field( stream, bo, type );
        read_field( stream, bo, maptype );
        read_field( stream, bo, maplength );
    }

    void SunHeader::to_stream( std::ofstream & stream, const byteorder & bo )
    {
        write_field( stream, bo, width );
        write_field( stream, bo, height );
        write_field( stream, bo, depth );
        write_field( stream, bo, length );
        write_field( stream, bo, type );
        write_field( stream, bo, maptype );
        write_field( stream, bo, maplength );
    }

    struct SunDecoderImpl
    {
        // attributes

        SunHeader header;
        std::ifstream stream;
        byteorder bo;
        void_vector< unsigned char > maps, bands;
        unsigned int components, row_stride;
        bool recode;

        // methods

        unsigned char get_val( unsigned int i )
        {
            return ( header.depth == 1 ) ?
                bands[ i / 8 ] & ( 1 << ( 7 - ( i % 8 ) ) ) : bands[i];
        }

        void read_scanline();

        // ctor

        SunDecoderImpl( const std::string & filename );
    };

    SunDecoderImpl::SunDecoderImpl( const std::string & filename )
#ifdef _MSC_VER
        : stream( filename.c_str(), std::ios::binary ), 
#else
        : stream( filename.c_str() ), 
#endif
          bo("big endian"), 
          maps(0), 
          bands(0),
          recode(false)
    {
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
        
        // read the magic number, adjust byte order if necessary
        SunHeader::field_type magic;
        read_field( stream, bo, magic );
        if ( magic == RAS_MAGIC_REVERSE ) {
            bo.set("little endian");
        } else {
            vigra_precondition( magic == RAS_MAGIC,
                                "the stored magic number is invalid" );
        }

        // read the header
        header.from_stream( stream, bo );

        // byte encoded files are not supported
        vigra_precondition( header.type != RT_ENCODED,
                            "ras byte encoding is not supported" );

        // calculate the row stride and adjust the bands vector
        // row stride is rounded up to the next multiple of 16 bits
        row_stride = ( 2 * header.width * ( header.depth / 8 ) + 1 ) / 2;
        bands.resize(row_stride);

        // read the color map, if there is one
        if ( header.maptype != RMT_NONE && header.maplength != 0 ) {

            // read the maps
            maps.resize(header.maplength);
            read_array( stream, bo, maps.data(), header.maplength );

            // map to gray or rgb
            recode = true;
        }

        // regenerate header.length, if needed
        if ( header.length == 0 )
            header.length = header.height * row_stride;

        // figure out the number of components
        switch (header.depth) {
        case 1:
            components = 1;
            // expand bi-level images to byte images
            recode = true;
            break;
        case 8:
            components = 1;
            break;
        case 24:
            components = 3;
            break;
        default:
            vigra_precondition( false, "number of bands is unsupported" );
        }
    }

    void SunDecoderImpl::read_scanline()
    {
        // read the scanline
        read_array( stream, bo, bands.data(), row_stride );

        // recode if necessary
        if (recode) {
            void_vector< unsigned char > recode_bands;
            if ( header.depth == 1 ) {
                // expand to unsigned char
                recode_bands.resize(header.width);
                for ( unsigned int i = 0; i < header.width; ++i )
                    recode_bands[i] = get_val(i);
            } else if ( header.maptype == RMT_EQUAL_RGB ) {
                recode_bands.resize( header.maptype == RMT_EQUAL_RGB ?
                                     3 * header.width : header.width );
                const unsigned int mapstride = header.maplength / 3;
                for ( unsigned int i = 0; i < header.width; ++i ) {
                    recode_bands[ 3 * i ] = maps[get_val(i)];
                    recode_bands[ 3 * i + 1 ]
                        = maps[ mapstride + get_val(i) ];
                    recode_bands[ 3 * i + 2 ]
                        = maps[ 2 * mapstride + get_val(i) ];
                }
            } else if ( header.maptype == RMT_RAW ) {
                recode_bands.resize(bands.size());
                for ( unsigned int i = 0; i < header.width; ++i )
                    recode_bands[i] = maps[get_val(i)];
            }
            swap_void_vector( recode_bands, bands );
        }

        // bgr -> rgb
        if ( header.type == RT_STANDARD && components == 3 ) {
            // rgb -> bgr
            void_vector< unsigned char > recode_bands(bands.size());
            for ( unsigned int i = 0; i < header.width; ++i ) {
                recode_bands[ 3 * i ] = bands[ 3 * i + 2 ];
                recode_bands[ 3 * i + 1 ] = bands[ 3 * i + 1 ];
                recode_bands[ 3 * i + 2 ] = bands[ 3 * i ];
            }
            swap_void_vector( recode_bands, bands );
        }
    }

    void SunDecoder::init( const std::string & filename )
    {
        pimpl = new SunDecoderImpl( filename );
    }

    SunDecoder::~SunDecoder()
    {
        delete pimpl;
    }

    std::string SunDecoder::getFileType() const
    {
        return "SUN";
    }

    unsigned int SunDecoder::getWidth() const
    {
        return pimpl->header.width;
    }

    unsigned int SunDecoder::getHeight() const
    {
        return pimpl->header.height;
    }

    unsigned int SunDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string SunDecoder::getPixelType() const
    {
        return "UINT8";
    }

    unsigned int SunDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * SunDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->bands.data() + band;
    }

    void SunDecoder::nextScanline()
    {
        pimpl->read_scanline();
    }

    void SunDecoder::close() {}
    void SunDecoder::abort() {}

    struct SunEncoderImpl
    {
        // attributes

        SunHeader header;
        std::ofstream stream;
        byteorder bo;
        void_vector< unsigned char > bands;
        unsigned int components, row_stride;
        bool finalized;

        // methods

        void finalize();
        void write_scanline();

        // ctor

        SunEncoderImpl( const std::string & filename );
    };

    SunEncoderImpl::SunEncoderImpl( const std::string & filename )
#ifdef _MSC_VER
        : stream( filename.c_str(), std::ios::binary ), 
#else
        : stream( filename.c_str() ), 
#endif
          bo("big endian"),
          bands(0), finalized(false)
    {
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
        // write the magic number
        SunHeader::field_type magic = RAS_MAGIC;
        write_field( stream, bo, magic );
    }

    void SunEncoderImpl::finalize()
    {
        // color depth
        vigra_precondition( components == 1 || components == 3,
                            "number of bands is not supported" );

        header.depth = components << 3;

        // calculate the row stride and adjust the bands vector
        // row stride is rounded up to the next multiple of 16 bits
        row_stride = ( 2 * header.width * ( header.depth / 8 ) + 1 ) / 2;
        bands.resize(row_stride);

        // set bands memory to zero
        for ( unsigned int i = 0; i < row_stride; ++i )
            bands[i] = 0;

        // set the band length
        header.length = header.height * row_stride;

        // standard format
        header.type = RT_STANDARD;

        // no colormap
        header.maptype = RMT_NONE;
        header.maplength = 0;

        // write the header
        header.to_stream( stream, bo );
    }

    void SunEncoderImpl::write_scanline()
    {
        if ( components == 3 ) {
            // rgb -> bgr
            void_vector< unsigned char > recode_bands(bands.size());
            for ( unsigned int i = 0; i < header.width; ++i ) {
                recode_bands[ 3 * i ] = bands[ 3 * i + 2 ];
                recode_bands[ 3 * i + 1 ] = bands[ 3 * i + 1 ];
                recode_bands[ 3 * i + 2 ] = bands[ 3 * i ];
            }
            swap_void_vector( recode_bands, bands );
        }
        write_array( stream, bo, bands.data(), row_stride );
    }

    void SunEncoder::init( const std::string & filename )
    {
        pimpl = new SunEncoderImpl(filename);
    }

    SunEncoder::~SunEncoder()
    {
        delete pimpl;
    }

    std::string SunEncoder::getFileType() const
    {
        return "SUN";
    }

    void SunEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.width = width;
    }

    void SunEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.height = height;
    }

    void SunEncoder::setNumBands( unsigned int numBands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->components = numBands;
    }

    void SunEncoder::setCompressionType( const std::string & comp,
                                         int quality )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    }

    void SunEncoder::setPixelType( const std::string & pixeltype )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        vigra_precondition( pixeltype == "UINT8",
                            "SunEncoder::setPixelType(): "
                            "SUN raster supports only the UINT8 pixeltype" );
    }

    unsigned int SunEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void SunEncoder::finalizeSettings()
    {
        pimpl->finalize();
        pimpl->finalized = true;
    }

    void * SunEncoder::currentScanlineOfBand( unsigned int band )
    {
        return pimpl->bands.data() + band;
    }

    void SunEncoder::nextScanline()
    {
        pimpl->write_scanline();
    }

    void SunEncoder::close()
    {
        nextScanline();
    }

    void SunEncoder::abort() {}
}
