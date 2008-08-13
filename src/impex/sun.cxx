/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
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

#include <fstream>
#include <stdexcept>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
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

        desc.bandNumbers.resize(2);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 3;
        
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
        typedef UInt32 field_type;

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
        void_vector< UInt8 > maps, bands;
        UInt32 components, row_stride;
        bool recode;

        // methods

        void read_scanline();

        // ctor

        SunDecoderImpl( const std::string & filename );
    };

    SunDecoderImpl::SunDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : stream( filename.c_str(), std::ios::binary ), 
#else
        : stream( filename.c_str() ), 
#endif
          bo ("big endian"), 
          maps (0), 
          bands (0),
          recode (false)
    {
        if (!stream.good ())
        {
            std::string msg ("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition (0, msg.c_str ());
        }
        
        // read the magic number, adjust byte order if necessary
        SunHeader::field_type magic;
        read_field (stream, bo, magic);
        if (magic == RAS_MAGIC_REVERSE)
            bo.set("little endian");
        else
            vigra_precondition
                (magic == RAS_MAGIC, "the stored magic number is invalid");

        // read the header
        header.from_stream (stream, bo);

        // byte encoded files are not supported
        vigra_precondition (header.type != RT_ENCODED,
                            "ras byte encoding is not supported");

        // calculate the row stride and adjust the bands vector
        // row stride is rounded up to the next multiple of 16 bits
        row_stride = (2*header.width*(header.depth/8)+1)/2;
        bands.resize (row_stride);

        // read the color map, if there is one
        if (header.maptype != RMT_NONE) {
            vigra_precondition
                (header.maplength != 0,
                 "mapping requested, but color maps have zero length");
            maps.resize (header.maplength);
            read_array (stream, bo, maps.data (), header.maplength);
        }

        // compute the header length, if it is not set.
        if (header.length == 0)
            header.length = header.height * row_stride;

        // find out if recoding is necessary.
        if (header.maptype != RMT_NONE || header.depth == 1)
            recode = true;
        else
            recode = false;

        // find out the number of components.
        if (header.depth == 24 || header.maptype == RMT_EQUAL_RGB)
            components = 3;
        else
            components = 1;

        // sanity check on the depth
        vigra_precondition
            (header.depth == 1 || header.depth == 8 || header.depth == 24,
             "unsupported color depth");
    }

    void SunDecoderImpl::read_scanline()
    {
        // read the scanline
        read_array (stream, bo, bands.data (), row_stride);

        // recode if necessary
        if (recode) {

            void_vector <UInt8> recode_bands;

            if (header.depth == 1) {

                // expand to UInt8.
                recode_bands.resize (header.width);
                for (unsigned int i = 0; i < header.width; ++i) {

                    // there are eight pixels in each byte.
                    const UInt8 b = bands [i/8];
                    recode_bands [i] = b >> i%8 & 0x01;
                }

                // commit.
                swap_void_vector (recode_bands, bands);
            }

            // color map the scanline.
            if (header.maptype == RMT_EQUAL_RGB) {

                // map from UInt8 to rgb
                recode_bands.resize (3*header.width);
                const unsigned int mapstride = header.maplength/3;
                UInt8 *recode_mover = recode_bands.data ();
                for (unsigned int i = 0; i < header.width; ++i) {
                    // find out the pointer to the red color
                    UInt8 *map_mover = maps.data () + bands [i];
                    // red
                    *recode_mover++ = *map_mover;
                    map_mover += mapstride;
                    // green
                    *recode_mover++ = *map_mover;
                    map_mover += mapstride;
                    // blue
                    *recode_mover++ = *map_mover;
                }
                
            } else if (header.maptype == RMT_RAW) {

                // map from UInt8 to UInt8
                recode_bands.resize (header.width);
                for (unsigned int i = 0; i < header.width; ++i)
                    recode_bands [i] = maps [bands [i]];
            }

            // commit.
            swap_void_vector (recode_bands, bands);
        }

        // swap the color components of a BGR image to RGB.
        // i really don't know the exact condition for this.
        if (header.type == RT_STANDARD && header.maptype != RMT_EQUAL_RGB
            && components == 3) {

            void_vector <UInt8> recode_bands (3*header.width);
            for (unsigned int i = 0; i < header.width; ++i) {
                recode_bands [3*i]   = bands [3*i+2];
                recode_bands [3*i+1] = bands [3*i+1];
                recode_bands [3*i+2] = bands [3*i];
            }

            // commit.
            swap_void_vector (recode_bands, bands);
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
        void_vector< UInt8 > bands;
        UInt32 components, row_stride;
        bool finalized;

        // methods

        void finalize();
        void write_scanline();

        // ctor

        SunEncoderImpl( const std::string & filename );
    };

    SunEncoderImpl::SunEncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
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
            void_vector< UInt8 > recode_bands(bands.size());
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
