/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2004 by Gunnar Kedenburg                */
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

#include <iostream>
#include <fstream>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "void_vector.hxx"
#include "byteorder.hxx"
#include "bmp.hxx"

// Windows Bitmap 3.0

#define BMP_RLE8 1
#define BMP_RLE4 2

namespace vigra {

CodecDesc BmpCodecFactory::getCodecDesc() const
{
    CodecDesc desc;

    // init file type
    desc.fileType = "BMP";

    // init pixel types
    desc.pixelTypes.resize(1);
    desc.pixelTypes[0] = "UINT8";

    // init compression types
    desc.compressionTypes.resize(1);
    desc.compressionTypes[0] = "RLE";

    // init magic strings
    desc.magicStrings.resize(1);
    desc.magicStrings[0].resize(2);
    desc.magicStrings[0][0] = 'B';
    desc.magicStrings[0][1] = 'M';

    // init file extensions
    desc.fileExtensions.resize(1);
    desc.fileExtensions[0] = "bmp";

    desc.bandNumbers.resize(2);
    desc.bandNumbers[0] = 1;
    desc.bandNumbers[1] = 3;

    return desc;
}

std::auto_ptr<Decoder> BmpCodecFactory::getDecoder() const
{
    return std::auto_ptr<Decoder>( new BmpDecoder() );
}

std::auto_ptr<Encoder> BmpCodecFactory::getEncoder() const
{
    return std::auto_ptr<Encoder>( new BmpEncoder() );
}

struct BmpFileHeader
{
    // attributes

    // the magic number
    UInt16 magic;
    // size of the whole file
    UInt32 size;
    // offset (from this field) to the data
    UInt32 offset;

    // ctor

    BmpFileHeader();

    // methods

    void from_stream( std::ifstream & stream, byteorder & bo );
    void to_stream( std::ofstream & stream, byteorder & bo );
};

BmpFileHeader::BmpFileHeader()
{
    // "BM"
    magic = 0x4D42;
}

void BmpFileHeader::from_stream( std::ifstream & stream, byteorder & bo )
{
    UInt16 filemagic;
    read_field( stream, bo, filemagic );
    vigra_precondition( filemagic == magic, "magic value is incorrect." );
    read_field( stream, bo, size );
    stream.seekg( 4, std::ios::cur );
    read_field( stream, bo, offset );
}

void BmpFileHeader::to_stream( std::ofstream & stream, byteorder & bo )
{
    write_field( stream, bo, magic );
    write_field( stream, bo, size );
    for( unsigned int i = 0; i < 4; ++i )
        stream.put(0); // clear reserved fields
    write_field( stream, bo, offset );
}

struct BmpInfoHeader
{
    // attributes

    // size of this header in the file
    UInt32 info_size;
    // image dimensions
    Int32 width, height;
    // number of planes, always set to one
    UInt16 planes;
    // bits per pixel
    UInt16 bit_count;
    // compression type
    UInt32 compression;
    // image size in bytes, may be zero for 24 bpp images
    UInt32 image_size;
    // image resolution
    Int32 x_pixels_per_meter, y_pixels_per_meter;
    // number of used colors, may be zero
    UInt32 clr_used;
    // number of important colors, may be zero
    UInt32 clr_important;

    // methods

    void from_stream( std::ifstream & stream, byteorder & bo );
    void to_stream( std::ofstream & stream, byteorder & bo );
};

void BmpInfoHeader::from_stream( std::ifstream & stream, byteorder & bo )
{
    const UInt32 info_impl_size = 40;
    read_field( stream, bo, info_size );
    vigra_precondition( info_size >= info_impl_size,
                        "info header has a wrong size" );
    read_field( stream, bo, width );
    vigra_precondition( width > 0, "width must be > 0" );
    read_field( stream, bo, height );
    vigra_precondition( height > 0, "height must be > 0" );
    read_field( stream, bo, planes );
    vigra_precondition( planes == 1, "planes must be 1" );
    read_field( stream, bo, bit_count );
    vigra_precondition( bit_count == 1 || bit_count == 4 || bit_count == 8
                        || bit_count == 24, "illegal bit count" );
    read_field( stream, bo, compression );
    read_field( stream, bo, image_size );
    vigra_precondition( image_size != 0 || bit_count == 24,
                        "illegal image size" );
    read_field( stream, bo, x_pixels_per_meter );
    read_field( stream, bo, y_pixels_per_meter );
    read_field( stream, bo, clr_used );
    const unsigned int max_colors = 1 << bit_count;
    vigra_precondition( clr_used <= max_colors,
                        "used colors field invalid" );
    read_field( stream, bo, clr_important );
    vigra_precondition( clr_important <= max_colors,
                        "important colors field invalid" );
    // skip any padding
    stream.seekg( info_size - info_impl_size, std::ios::cur );
}

void BmpInfoHeader::to_stream( std::ofstream & stream, byteorder & bo )
{
    write_field( stream, bo, info_size );
    write_field( stream, bo, width );
    write_field( stream, bo, height );
    write_field( stream, bo, planes = 1 );
    write_field( stream, bo, bit_count );
    write_field( stream, bo, compression );
    write_field( stream, bo, image_size );
    write_field( stream, bo, x_pixels_per_meter );
    write_field( stream, bo, y_pixels_per_meter );
    write_field( stream, bo, clr_used );
    write_field( stream, bo, clr_important );
}

struct BmpDecoderImpl
{
    // attributes

    // data source
    std::ifstream stream;

    // bmp headers
    BmpFileHeader file_header;
    BmpInfoHeader info_header;

    // image containers
    void_vector< UInt8 > pixels;
    void_vector< UInt8 > map;

    int scanline;

    // this flag is figured out from the info header's bit_count
    // field and the colormap, if there is one.
    bool grayscale;

    // the data is read when the first scanline is requested.
    bool data_read;

    // methods

    void read_data ();
    void read_colormap ();
    void read_1bit_data ();
    void read_rle4_data ();
    void read_4bit_data ();
    void read_rle8_data ();
    void read_8bit_data ();
    void read_rgb_data ();

    // ctor

    BmpDecoderImpl( const std::string & filename );
};


// reads the header.
BmpDecoderImpl::BmpDecoderImpl( const std::string & filename )
    :
#ifdef VIGRA_NEED_BIN_STREAMS
      stream (filename.c_str (), std::ios::binary),
#else
      stream (filename.c_str ()),
#endif
      scanline(-1)
{
    if( !stream.good() )
    {
        std::string msg("Unable to open file '");
        msg += filename;
        msg += "'.";
        vigra_precondition(0, msg.c_str());
    }
    byteorder bo( "little endian" );

    // read the header
    file_header.from_stream( stream, bo );
    info_header.from_stream( stream, bo );

    // read the map, if there is one, to determine if this is a grayscale
    // or rgb image.
    grayscale = false; // default is rgb
    if (info_header.bit_count != 24)
        read_colormap ();

    data_read = false;
}

void BmpDecoderImpl::read_data ()
{
    // read (and eventually map) the bands
    switch (info_header.bit_count) {
    case 1:
        read_1bit_data ();
        break;
    case 4:
        if (info_header.compression) {
            read_rle4_data ();
        } else {
            read_4bit_data ();
        } break;
    case 8:
        if (info_header.compression) {
            read_rle8_data ();
        } else {
            read_8bit_data ();
        } break;
    case 24:
        read_rgb_data ();
        break;
    }
    data_read = true;
}

void BmpDecoderImpl::read_colormap ()
{
    const unsigned int num_colors = 1 << info_header.bit_count;
    map.resize( 3 * num_colors );
    grayscale = true;
    for ( unsigned int i = 0; i < num_colors; ++i ) {
        map[ 3 * i + 2 ] = stream.get();
        map[ 3 * i + 1 ] = stream.get();
        map[ 3 * i     ] = stream.get();
        stream.get(); // skip unused byte
        grayscale = grayscale && ( map[ 3 * i ] == map[ 3 * i + 1 ] );
        grayscale = grayscale && ( map[ 3 * i + 1 ] == map[ 3 * i + 2 ] );
    }
}

void BmpDecoderImpl::read_1bit_data ()
{
    // these are sizes for the final image
    const unsigned int ncomp = grayscale ? 1 : 3;
    const unsigned int line_size = info_header.width * ncomp;
    const unsigned int image_size = info_header.height * line_size;
    int c = 0;

    // seek to the data
    stream.seekg( file_header.offset, std::ios::beg );

    // make room for the pixels
    pixels.resize(image_size);

    // padding after each line (must end on a 32-bit boundary)
    unsigned int pad_size = ( ( info_header.width + 7 ) / 8 ) % 4;
    if ( pad_size > 0 )
        pad_size = 4 - pad_size;

    // setup the base pointer at one line after the end
    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover = base;

    // read scanlines from bottom to top
    for ( int y = info_header.height - 1; y >= 0; --y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // read the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {

            // eventually get a new byte from the stream
            if ( x % 8 == 0 )
                c = stream.get();

            // get the color bit
            const UInt8 index = ( c >> ( 7 - ( x % 8 ) ) ) & 0x01;

            // map and assign the pixel
            const UInt8 * map_base = map.data() + 3 * index;
            for ( unsigned int i = 0; i < ncomp; ++i )
                mover[i] = map_base[i];
            mover += ncomp;
        }

        // advance to the next 32 bit boundary
        stream.seekg( pad_size, std::ios::cur );
    }
}

void BmpDecoderImpl::read_4bit_data ()
{
    // these are sizes for the final image
    const unsigned int ncomp = grayscale ? 1 : 3;
    const unsigned int line_size = info_header.width * ncomp;
    const unsigned int image_size = info_header.height * line_size;
    int c = 0;

    // seek to the data
    stream.seekg( file_header.offset, std::ios::beg );

    // make room for the pixels
    pixels.resize(image_size);

    // padding after each line (must end on a 32-bit boundary)
    unsigned int pad_size = ( ( info_header.width + 1 ) / 2 ) % 4;
    if ( pad_size > 0 )
        pad_size = 4 - pad_size;

    // setup the base pointer at one line after the end
    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover = base;

    // read scanlines from bottom to top
    for ( int y = info_header.height - 1; y >= 0; --y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // read the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {

            // eventually get a new byte from the stream
            if ( x % 2 == 0 )
                c = stream.get();

            // get the color index
            const UInt8 index = ( c >> ( 1 - ( x % 2 ) ) ) & 0x0f;

            // map and assign the pixel
            const UInt8 * map_base = map.data() + 3 * index;
            for ( unsigned int i = 0; i < ncomp; ++i )
                mover[i] = map_base[i];
            mover += ncomp;
        }

        // advance to the next 32 bit boundary
        stream.seekg( pad_size, std::ios::cur );
    }
}

void BmpDecoderImpl::read_rle4_data ()
{
    // these are sizes for the final image
    const unsigned int ncomp = grayscale ? 1 : 3;
    const unsigned int line_size = info_header.width * ncomp;
    const unsigned int image_size = info_header.height * line_size;

    // seek to the data
    stream.seekg (file_header.offset, std::ios::beg);

    // make room for the pixels
    pixels.resize (image_size);

    // setup the base pointer at the beginning of the last line.
    UInt8 * base = pixels.data() + image_size - line_size;
    UInt8 * mover = base;

    // set the image's background color to black.
    VIGRA_CSTD::memset (pixels.data (), 0, image_size);

    // read the rle codes one by one.
    int c1, c2;
    bool painting = true;

    int x = 0;
    int y = 0;

    while (painting) {

        c1 = stream.get ();
        c2 = stream.get ();

        if (c1 == 0x00) {

            // escape sequences
            switch (c2) {
            case 0x00:
                {
                    // end of line

                    // finished painting this line, move the pointer to
                    // the beginning of the previous line, which is the
                    // next bmp line.
                    // we need to go back to the beginning of this line,
                    // and another additional line back.
                    mover -= line_size + ncomp*x;

                    // so we are at the beginning of the next bmp line.
                    x = 0;
                    ++y;

                    break;
                }
            case 0x01:
                {
                    // end of bitmap
                    painting = false;
                    break;
                }
            case 0x02:
                {
                    // "delta" movement

                    // handle the end-of-line case.
                    if (x == info_header.width) {
                        mover -= line_size + ncomp*x;
                        x = 0;
                        ++y;
                    }

                    unsigned int dx = stream.get ();
                    unsigned int dy = stream.get ();
                    int nx = x + dx;

                    // x movement.
                    if (nx > info_header.width) {
                        dy += nx / info_header.width + 1;
                        nx %= info_header.width;
                    }
                    mover += ncomp*(nx-x);
                    x = nx;

                    // y movement.
                    if (dy != 0) {
                        mover -= line_size*dy;
                        y += dy;
                    }

                    break;
                }
            default:

                // absolute mode. paint the following c2 nibbles.
                // then, eventually skip one byte to respect 16 bit boundaries.
                int k = 0;

                while (true) {
                    const int c = stream.get ();

                    // the high-order and lower-order nibbles
                    const int nh = (c & 0xf0) >> 4;
                    const int nl = (c & 0x0f);

                    // paint the the higher-order nibble.
                    const UInt8 *map_base_h = map.data () + 3*nh;
                    unsigned int j;
                    for (j = 0; j < ncomp; ++j)
                        mover [j] = map_base_h [j];
                    mover += ncomp;
                    if (++k >= c2)
                        break;

                    // paint the lower-order nibble.
                    const UInt8 *map_base_l = map.data () + 3*nl;
                    for (j = 0; j < ncomp; ++j)
                        mover [j] = map_base_l [j];
                    mover += ncomp;
                    if (++k >= c2)
                        break;
                }

                if (c2 % 2)
                    stream.get ();
            }

        } else {

            // plain rle: repeat the next byte mapped, c1 times.
            // a line break may not happen here.
            for (int i = 0; i < c1; ++i) {
                // the high-order and lower-order nibbles
                const int nh = (c2 & 0xf0) >> 4;
                const int nl = (c2 & 0x0f);
                // paint the higher-order nibble.
                const UInt8 *map_base_h = map.data () + 3*nh;
                unsigned int j;
                for (j = 0; j < ncomp; ++j)
                    mover [j] = map_base_h [j];
                mover += ncomp;
                // paint the lower-order nibble.
                const UInt8 *map_base_l = map.data () + 3*nl;
                for (j = 0; j < ncomp; ++j)
                    mover [j] = map_base_l [j];
                mover += ncomp;
            }
            x += c1;
        }
    }
}

void BmpDecoderImpl::read_8bit_data ()
{
    // these are sizes for the final image
    const unsigned int ncomp = grayscale ? 1 : 3;
    const unsigned int line_size = info_header.width * ncomp;
    const unsigned int image_size = info_header.height * line_size;

    // seek to the data
    stream.seekg( file_header.offset, std::ios::beg );

    // make room for the pixels
    pixels.resize(image_size);

    // padding after each line (must end on a 32-bit boundary)
    unsigned int pad_size = info_header.width % 4;
    if ( pad_size > 0 )
        pad_size = 4 - pad_size;

    // setup the base pointer at one line after the end
    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover = base;

    // read scanlines from bottom to top
    for ( int y = info_header.height - 1; y >= 0; --y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // read the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {

            // get the color byte
            const int index = stream.get();

            // map and assign the pixel
            const UInt8 * map_base = map.data() + 3 * index;
            for ( unsigned int i = 0; i < ncomp; ++i )
                mover[i] = map_base[i];
            mover += ncomp;
        }

        // advance to the next 32 bit boundary
        stream.seekg( pad_size, std::ios::cur );
    }
}

void BmpDecoderImpl::read_rle8_data ()
{
    // these are sizes for the final image
    const unsigned int ncomp = grayscale ? 1 : 3;
    const unsigned int line_size = info_header.width * ncomp;
    const unsigned int image_size = info_header.height * line_size;

    // seek to the data
    stream.seekg (file_header.offset, std::ios::beg);

    // make room for the pixels
    pixels.resize (image_size);

    // setup the base pointer at the beginning of the last line.
    UInt8 * base = pixels.data() + image_size - line_size;
    UInt8 * mover = base;

    // set the image's background color to black.
    VIGRA_CSTD::memset (pixels.data (), 0, image_size);

    // read the rle codes one by one.
    int c1, c2;
    bool painting = true;

    int x = 0;
    int y = 0;

    while (painting) {

        c1 = stream.get ();
        c2 = stream.get ();

        if (c1 == 0x00) {

            // escape sequences
            switch (c2) {
            case 0x00:
                {
                    // end of line

                    // finished painting this line, move the pointer to
                    // the beginning of the previous line, which is the
                    // next bmp line.
                    // we need to go back to the beginning of this line,
                    // and another additional line back.
                    mover -= line_size + ncomp*x;

                    // so we are at the beginning of the next bmp line.
                    x = 0;
                    ++y;

                    break;
                }
            case 0x01:
                {
                    // end of bitmap
                    painting = false;
                    break;
                }
            case 0x02:
                {
                    // "delta" movement

                    // handle the end-of-line case.
                    if (x == info_header.width) {
                        mover -= line_size + ncomp*x;
                        x = 0;
                        ++y;
                    }

                    unsigned int dx = stream.get ();
                    unsigned int dy = stream.get ();
                    int nx = x + dx;

                    // x movement.
                    if (nx > info_header.width) {
                        dy += nx / info_header.width + 1;
                        nx %= info_header.width;
                    }
                    mover += ncomp*(nx-x);
                    x = nx;

                    // y movement.
                    if (dy != 0) {
                        mover -= line_size*dy;
                        y += dy;
                    }
                    break;
                }
            default:

                // absolute mode. paint the following c2 bytes.
                // then, eventually skip one byte to respect 16 bit boundaries.
                for (int k = 0; k < c2; ++k) {
                    const int c = stream.get ();
                    const UInt8 *map_base = map.data () + 3*c;
                    for (unsigned int j = 0; j < ncomp; ++j)
                        mover [j] = map_base [j];
                    mover += ncomp;
                }
                if (c2 % 2)
                    stream.get ();
            }

        } else {

            // plain rle: repeat c2 mapped, c1 times.
            // a line break may not happen here.
            for (int i = 0; i < c1; ++i) {
                const UInt8 *map_base = map.data () + 3*c2;
                for (unsigned int j = 0; j < ncomp; ++j)
                    mover [j] = map_base [j];
                mover += ncomp;
            }
            x += c1;
        }
    }
}

void BmpDecoderImpl::read_rgb_data ()
{
    const unsigned int line_size = 3 * info_header.width;
    const unsigned int image_size = info_header.height * line_size;

    // seek to the data
    stream.seekg( file_header.offset, std::ios::beg );

    // make room for the pixels
    pixels.resize(image_size);

    // padding after each scanline.
    // citing Kevin D. Quitt's mail to wotsit.org:
    // In RGB encoding (no compression), when using 8 bits per pixel, lines
    // must start on a long-word boundary (i.e., low two bits zero).
    const unsigned int pad_size = (line_size % 4) ? 4 - (line_size % 4) : 0;

    // setup the base pointer at one line after the end
    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover;

    // read scanlines from bottom to top
    for ( int y = info_header.height - 1; y >= 0; --y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // read the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {
            mover[2] = stream.get(); // B
            mover[1] = stream.get(); // G
            mover[0] = stream.get(); // R
            mover += 3;
        }

        // advance to the next 32 bit boundary
        stream.seekg( pad_size, std::ios::cur );
    }
}

void BmpDecoder::init( const std::string & filename )
{
    pimpl = new BmpDecoderImpl( filename.c_str() );
}

BmpDecoder::~BmpDecoder()
{
    delete pimpl;
}

std::string BmpDecoder::getFileType() const
{
    return "BMP";
}

unsigned int BmpDecoder::getWidth() const
{
    return pimpl->info_header.width;
}

unsigned int BmpDecoder::getHeight() const
{
    return pimpl->info_header.height;
}

unsigned int BmpDecoder::getNumBands() const
{
    return pimpl->grayscale ? 1 : 3;
}

std::string BmpDecoder::getPixelType() const
{
    return "UINT8";
}

unsigned int BmpDecoder::getOffset() const
{
    return pimpl->grayscale ? 1 : 3;
}

const void * BmpDecoder::currentScanlineOfBand( unsigned int band ) const
{
    if (!pimpl->data_read)
        pimpl->read_data ();
    return pimpl->pixels.data() + band
        + pimpl->scanline * ( pimpl->grayscale ? 1 : 3 )
        * pimpl->info_header.width;
}

void BmpDecoder::nextScanline()
{
    ++(pimpl->scanline);
}

void BmpDecoder::close() {}

void BmpDecoder::abort() {}

struct BmpEncoderImpl
{
    // attributes

    // bmp headers
    BmpFileHeader file_header;
    BmpInfoHeader info_header;

    // output stream
    byteorder bo;
    std::ofstream stream;

    // image container
    void_vector< UInt8 > pixels;

    int scanline;

    // this flag is set when the number of bands is set to one
    bool grayscale;

    // finalized settings
    bool finalized;

    // ctor

    BmpEncoderImpl( const std::string & );

    // methods

    void finalize();
    void write();
    void write_colormap();
    void write_8bit_data();
    void write_rgb_data();
};

BmpEncoderImpl::BmpEncoderImpl( const std::string & filename )
    : bo( "little endian" ),
#ifdef VIGRA_NEED_BIN_STREAMS
      stream( filename.c_str(), std::ios::binary ),
#else
      stream( filename.c_str() ),
#endif
      scanline(0), finalized(false)
{
    if( !stream.good() )
    {
        std::string msg("Unable to open file '");
        msg += filename;
        msg += "'.";
        vigra_precondition(0, msg.c_str());
    }
}

void BmpEncoderImpl::finalize()
{
    if ( grayscale ) {

        unsigned int line_size = 3 * info_header.width;
        unsigned int pad_size = info_header.width % 4;
        if ( pad_size > 0 )
            pad_size = 4 - pad_size;
        line_size += pad_size;
        const unsigned int cmap_size = 4 * 256;
        file_header.size = 10 + 40 + cmap_size
            + line_size * info_header.height;
        file_header.offset = 10 + 40 + 4 + cmap_size;
        info_header.info_size = 40;
        info_header.planes = 1;
        info_header.bit_count = 8;
        info_header.compression = 0;
        info_header.image_size = line_size * info_header.height;
        info_header.x_pixels_per_meter = 0;
        info_header.y_pixels_per_meter = 0;
        info_header.clr_used = 256;
        info_header.clr_important = 256;

    } else {

        // write the header fields
        file_header.size = 10 + 40
            + 3 * info_header.width * info_header.height;
        file_header.offset = 10 + 40 + 4;
        info_header.info_size = 40;
        info_header.planes = 1;
        info_header.bit_count = 24;
        info_header.compression = 0;
        info_header.image_size = 0;
        info_header.x_pixels_per_meter = 0;
        info_header.y_pixels_per_meter = 0;
        info_header.clr_used = 0;
        info_header.clr_important = 0;
    }

    // resize the vector, now that width and height are clear
    pixels.resize( ( grayscale ? 1 : 3 ) * info_header.width * info_header.height );

    // ok, finalize.
    finalized = true;
}

void BmpEncoderImpl::write()
{
    // write the headers
    file_header.to_stream( stream, bo );
    info_header.to_stream( stream, bo );

    if (grayscale) {
        write_colormap();
        write_8bit_data();
    } else {
        write_rgb_data();
    }
}

void BmpEncoderImpl::write_colormap()
{
    const unsigned int num_colors = 256;
    for ( unsigned int i = 0; i < num_colors; ++i ) {
        const int c = i;
        stream.put(c);
        stream.put(c);
        stream.put(c);
        stream.put(0); // skip unused byte
    }
}

void BmpEncoderImpl::write_8bit_data()
{
    const unsigned int line_size = info_header.width;
    const unsigned int image_size = info_header.height * line_size;
    unsigned int pad_size = info_header.width % 4;
    if ( pad_size > 0 )
        pad_size = 4 - pad_size;

    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover = base;

    // write scanlines, the scanlines are already in bottom-to-top
    // order.
    for ( int y = 0; y < info_header.height; ++y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // write the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {
            stream.put(*mover);
            ++mover;
        }

        // pad
        for ( unsigned int p = 0; p < pad_size; ++p )
            stream.put(0);
    }
}

void BmpEncoderImpl::write_rgb_data()
{
    const unsigned int line_size = 3 * info_header.width;
    const unsigned int image_size = info_header.height * line_size;
    const unsigned int pad_size = (line_size % 4) ? 4 - (line_size % 4) : 0;

    UInt8 * base = pixels.data() + image_size;
    UInt8 * mover = base;

    // write scanlines, the scanlines are already in bottom-to-top
    // order.
    for ( int y = 0; y < info_header.height; ++y ) {

        // set the base pointer to the beginning of the line to be read
        base -= line_size;
        mover = base;

        // write the line from left to right
        for ( int x = 0; x < info_header.width; ++x ) {
            stream.put(mover[2]); // B
            stream.put(mover[1]); // G
            stream.put(mover[0]); // R
            mover += 3;
        }

        // pad
        for ( unsigned int p = 0; p < pad_size; ++p )
	    stream.put(0);
    }
}

void BmpEncoder::init( const std::string & filename )
{
    pimpl = new BmpEncoderImpl(filename);
}

BmpEncoder::~BmpEncoder()
{
    delete pimpl;
}

std::string BmpEncoder::getFileType() const
{
    return "BMP";
}

void BmpEncoder::setWidth( unsigned int width )
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    pimpl->info_header.width = width;
}

void BmpEncoder::setHeight( unsigned int height )
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    pimpl->info_header.height = height;
}

void BmpEncoder::setNumBands( unsigned int bands )
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    vigra_precondition( ( bands == 1 ) || ( bands == 3 ),
                        "bmp supports only rgb and grayscale images" );
    pimpl->grayscale = ( bands == 1 );
}

void BmpEncoder::setCompressionType( const std::string & /* comp */, int /* quality */)
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
}

void BmpEncoder::setPixelType( const std::string & pixelType )
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    vigra_precondition( pixelType == "UINT8",
                        "bmp supports only the UINT8 pixeltype" );
}

unsigned int BmpEncoder::getOffset() const
{
    return pimpl->grayscale ? 1 : 3;
}

void BmpEncoder::finalizeSettings()
{
    VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    pimpl->finalize();
}

void * BmpEncoder::currentScanlineOfBand( unsigned int band )
{
    if (pimpl->grayscale) {
        const unsigned int row_stride = pimpl->info_header.width;
        return pimpl->pixels.data() + ( row_stride * pimpl->scanline );
    } else {
        const unsigned int row_stride = 3 * pimpl->info_header.width;
        return pimpl->pixels.data() + ( row_stride * pimpl->scanline ) + band;
    }
}

void BmpEncoder::nextScanline()
{
    ++(pimpl->scanline);
}

void BmpEncoder::close()
{
    pimpl->write();
}

void BmpEncoder::abort() {}
}
