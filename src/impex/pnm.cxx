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

#include <cmath>
#include <cctype>
// under Solaris, this is a macro leading to syntax errors:
#undef isspace
#include <iostream>
#include <fstream>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "void_vector.hxx"
#include "pnm.hxx"
#include "byteorder.hxx"

namespace vigra {

    CodecDesc PnmCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "PNM";

        // init pixel types
        desc.pixelTypes.resize(3);
        desc.pixelTypes[0] = "UINT8";
        desc.pixelTypes[1] = "UINT16";
        desc.pixelTypes[2] = "UINT32";

        // init compression types
        desc.compressionTypes.resize(3);
        desc.compressionTypes[0] = "ASCII";
        desc.compressionTypes[1] = "RAW";
        desc.compressionTypes[2] = "BILEVEL";

        // init magic strings
        desc.magicStrings.resize(6);
        desc.magicStrings[0].resize(2);
        desc.magicStrings[0][0] = 'P';
        desc.magicStrings[0][1] = '1';
        desc.magicStrings[1].resize(2);
        desc.magicStrings[1][0] = 'P';
        desc.magicStrings[1][1] = '2';
        desc.magicStrings[2].resize(2);
        desc.magicStrings[2][0] = 'P';
        desc.magicStrings[2][1] = '3';
        desc.magicStrings[3].resize(2);
        desc.magicStrings[3][0] = 'P';
        desc.magicStrings[3][1] = '4';
        desc.magicStrings[4].resize(2);
        desc.magicStrings[4][0] = 'P';
        desc.magicStrings[4][1] = '5';
        desc.magicStrings[5].resize(2);
        desc.magicStrings[5][0] = 'P';
        desc.magicStrings[5][1] = '6';

        // init file extensions
        desc.fileExtensions.resize(4);
        desc.fileExtensions[0] = "pnm";
        desc.fileExtensions[1] = "pbm";
        desc.fileExtensions[2] = "pgm";
        desc.fileExtensions[3] = "ppm";

        desc.bandNumbers.resize(2);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 3;

        return desc;
    }

    std::auto_ptr<Decoder> PnmCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new PnmDecoder() );
    }

    std::auto_ptr<Encoder> PnmCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new PnmEncoder() );
    }

    struct PnmDecoderImpl
    {
        // data source
        std::ifstream stream;

        // image container
        void_vector_base bands;

        // data storage
        bool raw, bilevel;

        // image dimensions
        unsigned int width, height, components;

        // pixel type
        std::string pixeltype;

        // skip whitespace
        void skip_whitespace();

        // scanline reading
        void read_bilevel_ascii_scanline();
        void read_ascii_scanline();
        void read_bilevel_raw_scanline();
        void read_raw_scanline();

        void read_raw_scanline_uchar();
        void read_raw_scanline_ushort();
        void read_raw_scanline_uint();

        // skip whitespace and comment blocks
        void skip();

        // ctor
        PnmDecoderImpl( const std::string & );
    };

    void PnmDecoderImpl::skip_whitespace()
    {
        while ( std::isspace( stream.peek() ) )
            stream.get();
    }

    void PnmDecoderImpl::skip()
    {
        // skip whitespace
        skip_whitespace();

        // skip comments
        while ( stream.peek() == '#' ) {
            // skip line
            while( stream.peek() != '\n' )
                stream.get();
            // skip whitespace
            skip_whitespace();
        }
    }

    void PnmDecoderImpl::read_bilevel_ascii_scanline()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // read and store
        for( unsigned int i = 0; i < width * components; ++i ) {
            skip_whitespace();
            cbands[i] = ( stream.get() - '0' ) * 255;
        }
    }

    void PnmDecoderImpl::read_ascii_scanline()
    {
        // XXX implement this also for other pixel types than unsigned char

        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // read and store
        int x;
        for( unsigned int i = 0; i < width * components; ++i ) {
            skip_whitespace();
            stream >> x;
            cbands[i] = x;
        }
    }

    void PnmDecoderImpl::read_bilevel_raw_scanline()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // read and store
        UInt8 buf = 0;
        const unsigned int size = width / 8;  // XXX wrong assumption
        for( unsigned int i = 0; i < size; ++i ) {
            stream.read( reinterpret_cast< char * >(&buf), 1 );
            for( unsigned int j = 0; j < 8; ++j )
                cbands[ 8 * i + j ] = ( ( buf >> j ) & 1 ) * 255;
        }
    }

    void PnmDecoderImpl::read_raw_scanline()
    {
        if ( pixeltype == std::string("UINT8") ) {
            read_raw_scanline_uchar();
        }
        if ( pixeltype == std::string("UINT16") ) {
            read_raw_scanline_ushort();
        }
        if ( pixeltype == std::string("UINT32") ) {
            read_raw_scanline_uint();
        }
    }

    void PnmDecoderImpl::read_raw_scanline_uchar()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // read and store
        stream.read( reinterpret_cast< char * >(cbands.data()),
                     width * components );
    }

    void PnmDecoderImpl::read_raw_scanline_ushort()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt16 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // read and store, need to swap bytes.
        byteorder bo( "big endian" );
        read_array( stream, bo, reinterpret_cast< UInt16 * >(cbands.data()),
                    width * components );
    }

    void PnmDecoderImpl::read_raw_scanline_uint()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt32 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        byteorder bo( "big endian" );
        read_array( stream, bo, reinterpret_cast< UInt32 * >(cbands.data()),
                    width * components );
    }

    // reads the header.
    PnmDecoderImpl::PnmDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : stream( filename.c_str(), std::ios::binary )
#else
        : stream( filename.c_str() )
#endif
    {
        long maxval = 1;
        char type;

        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }

        // read the pnm header
        vigra_postcondition( stream.get() == 'P', "bad magic number" );

        // read the type, find out if the file is raw or ascii
        type = stream.get();

        switch (type) {
        case '1': // plain bitmap
            raw = false;
            bilevel = true;
            components = 1;
            maxval = 1;
            pixeltype = "UINT8";
            break;
        case '2': // plain graymap
            raw = false;
            bilevel = false;
            components = 1;
            break;
        case '3': // plain pixmap
            raw = false;
            bilevel = false;
            components = 3;
            break;
        case '4': // raw bitmap
            raw = true;
            bilevel = true;
            components = 1;
            maxval = 1;
            pixeltype = "UINT8";
            break;
        case '5': // raw graymap
            raw = true;
            bilevel = false;
            components = 1;
            maxval = 255;
            pixeltype = "UINT8";
            break;
        case '6': // raw pixmap
            raw = true;
            bilevel = false;
            components = 3;
            maxval = 255;
            pixeltype = "UINT8";
            break;
        default:
            vigra_precondition( false, "unknown magic number in file" );
        }

        // read width, height and maxval
        skip();
        stream >> width;
        skip();
        stream >> height;

        // bitmaps implicitly have maxval 1
        if ( type != '1' && type != '4' ) {
            skip();
            stream >> maxval;
        }

        // select a pixeltype depending on maxval
        int bits = 0;
        do
        {
            ++bits;
            maxval >>= 1;
        }
        while(maxval > 0);

        vigra_precondition( bits >= 0, "the file's maxval field is corrupt" );
        if ( bits <= 8 )
            pixeltype = "UINT8";
        else if ( bits <= 16 )
            pixeltype = "UINT16";
        else if ( bits <= 32 )
            pixeltype = "UINT32";
        else
            vigra_precondition( false,
                                "the file's maxval field is too large" );

        // adjust the buffer size
        if ( pixeltype == "UINT8" )
            bands.resize( width * components );
        else if ( pixeltype == "UINT16" )
            bands.resize( width * components * 2 );
        else if ( pixeltype == "UINT32" )
            bands.resize( width * components * 4 );

#ifdef DEBUG
        // print stats
        std::cerr << width << "x" << height << "x" << components
                  << "x " << pixeltype << std::endl;
#endif

        // advance to the beginning of the "data section"
        if (raw == false)
          skip();
        else
        {
#if defined(__GNUC__) && __GNUC__ == 2
          typedef streamoff streamOffset;
#else
          typedef std::ifstream::off_type streamOffset;
#endif
          {
              UInt32 seekOffset = width * height * components;
              if ( pixeltype == "UINT8" )
                  seekOffset *= 1;
              else if ( pixeltype == "UINT16" )
                  seekOffset *= 2;
              else if ( pixeltype == "UINT32" )
                  seekOffset *= 4;

              stream.seekg( -static_cast<streamOffset>(seekOffset), std::ios::end );
          }
        }
    }

    void PnmDecoder::init( const std::string & filename )
    {
        pimpl = new PnmDecoderImpl( filename.c_str() );
    }

    PnmDecoder::~PnmDecoder()
    {
        delete pimpl;
    }

    std::string PnmDecoder::getFileType() const
    {
        return "PNM";
    }

    unsigned int PnmDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int PnmDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int PnmDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string PnmDecoder::getPixelType() const
    {
        return pimpl->pixeltype;
    }

    unsigned int PnmDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * PnmDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        if ( pimpl->pixeltype == "UINT8" ) {
            typedef void_vector< UInt8 > bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + band;
        } else if ( pimpl->pixeltype == "UINT16" ) {
            typedef void_vector<UInt16> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + band;
        } else if ( pimpl->pixeltype == "UINT32" ) {
            typedef void_vector<UInt32> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + band;
        }
        vigra_precondition( false, "internal error: unknown pixeltype" );
        return 0; // this is not reached
    }

    void PnmDecoder::nextScanline()
    {
        if ( pimpl->raw ) {
            if ( pimpl->bilevel ) pimpl->read_bilevel_raw_scanline();
            else pimpl->read_raw_scanline();
        } else {
            if ( pimpl->bilevel ) pimpl->read_bilevel_ascii_scanline();
            else pimpl->read_ascii_scanline();
        }
    }

    void PnmDecoder::close()
    {}

    void PnmDecoder::abort()
    {}

    struct PnmEncoderImpl
    {
        // data source
        std::ofstream stream;

        // image container
        void_vector_base bands;

        // raw data storage
        bool raw, bilevel;

        // finalized settings
        bool finalized;

        // image dimensions
        unsigned int width, height, components;
        unsigned int maxval;

        // current scanline number
        unsigned int scanline;

        // pixel type
        std::string pixeltype;

        // writing
        void write_bilevel_ascii();
        void write_ascii();
        void write_bilevel_raw();
        void write_raw();

        // ctor
        PnmEncoderImpl( const std::string & );
    };

    PnmEncoderImpl::PnmEncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : stream( filename.c_str(), std::ios::binary ),
#else
        : stream( filename.c_str() ),
#endif
          raw(true), bilevel(false), finalized(false), scanline(0)
    {
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
    }

    void PnmEncoder::init( const std::string & filename )
    {
        pimpl = new PnmEncoderImpl(filename);
    }

    PnmEncoder::~PnmEncoder()
    {
        delete pimpl;
    }

    std::string PnmEncoder::getFileType() const
    {
        return "PNM";
    }

    void PnmEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void PnmEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void PnmEncoder::setNumBands( unsigned int bands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->components = bands;
    }

    void PnmEncoder::setCompressionType( const std::string & comp, int /* quality */)
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( comp == "ASCII" )
            pimpl->raw = false;
        else if ( comp == "RAW" )
            pimpl->raw = true;
        else if ( comp == "BILEVEL" )
            pimpl->bilevel = true;
    }

    void PnmEncoder::setPixelType( const std::string & pixelType )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->pixeltype = pixelType;
    }

    unsigned int PnmEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void PnmEncoder::finalizeSettings()
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->finalized = true;

        if ( pimpl->pixeltype == "INT32" )
          pimpl->raw = false;

        // write magic number
        pimpl->stream << "P";
        if ( pimpl->components == 1 ) {
            if ( pimpl->bilevel ) {
                // bitmap
                pimpl->stream << (pimpl->raw ? "4" : "1");
            } else {
                // graymap
                pimpl->stream << (pimpl->raw ? "5" : "2");
            }
        } else if ( pimpl->components == 3 ) {
            // pixmap
            pimpl->stream << (pimpl->raw ? "6" : "3");
        } else
            vigra_precondition( false, "number of bands is not supported" );
        pimpl->stream << "\n";

        // write advertisement
        pimpl->stream << "# generated by the VIGRA library\n";

        // write width and height
        pimpl->stream << pimpl->width << " " << pimpl->height << std::endl;

        // allocate image memory
        if ( pimpl->pixeltype == "UINT8" )
            pimpl->bands.resize( pimpl->height * pimpl->width
                                 * pimpl->components );
        else if ( pimpl->pixeltype == "UINT16" )
            pimpl->bands.resize( 2 * pimpl->height * pimpl->width
                                 * pimpl->components );
        else if ( pimpl->pixeltype == "UINT32" )
            pimpl->bands.resize( 4 * pimpl->height * pimpl->width
                                 * pimpl->components );
    }

    void * PnmEncoder::currentScanlineOfBand( unsigned int band )
    {
        const unsigned int row_stride = pimpl->width * pimpl->components;
        if ( pimpl->pixeltype == "UINT8" ) {
            typedef void_vector< UInt8 > bands_type;
            bands_type & bands
                = static_cast< bands_type & >(pimpl->bands);
            return bands.data() + pimpl->scanline * row_stride + band;
        } else if ( pimpl->pixeltype == "UINT16" ) {
            typedef void_vector<UInt16> bands_type;
            bands_type & bands
                = static_cast< bands_type & >(pimpl->bands);
            return bands.data() + pimpl->scanline * row_stride + band;
        } else if ( pimpl->pixeltype == "UINT32" ) {
            typedef void_vector<UInt32> bands_type;
            bands_type & bands
                = static_cast< bands_type & >(pimpl->bands);
            return bands.data() + pimpl->scanline * row_stride + band;
        }
        vigra_postcondition( false, "internal error" );
        return 0;
    }

    void PnmEncoderImpl::write_bilevel_ascii()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // write and store
        UInt8 * iter = cbands.data();
        for( unsigned int i = 0; i < height; ++i ) {
            for ( unsigned int j = 0; j < width; ++j ) {
                for ( unsigned int k = 0; k < components; ++k ) {
                    const int value = *iter / 255;
                    ++iter;
                    stream << value + '0' << " ";
                }
                stream << " "; // separate pixels with an extra space
            }
            stream << std::endl; // separate lines with a newline
        }
    }

    void PnmEncoderImpl::write_ascii()
    {
        if ( pixeltype == "UINT8" ) {

          typedef void_vector< UInt8 > vector_type;
          vector_type & cbands = static_cast< vector_type & >(bands);
          UInt8 * iter = cbands.data();
          for( unsigned int i = 0; i < height; ++i ) {
            for ( unsigned int j = 0; j < width; ++j ) {
              for ( unsigned int k = 0; k < components; ++k ) {
                const int value = *iter;
                ++iter;
                stream << value << " ";
              }
              stream << " "; // separate pixels with an extra space
            }
            stream << std::endl; // separate lines with a newline
          }

        } else if ( pixeltype == "UINT16" ) {

          typedef void_vector<UInt16> vector_type;
          vector_type & cbands = static_cast< vector_type & >(bands);
          UInt16 * iter = cbands.data();
          for( unsigned int i = 0; i < height; ++i ) {
            for ( unsigned int j = 0; j < width; ++j ) {
              for ( unsigned int k = 0; k < components; ++k ) {
                const int value = *iter;
                ++iter;
                stream << value << " ";
              }
              stream << " "; // separate pixels with an extra space
            }
            stream << std::endl; // separate lines with a newline
          }

        } else if ( pixeltype == "UINT32" ) {

          typedef void_vector<UInt32> vector_type;
          vector_type & cbands = static_cast< vector_type & >(bands);
          UInt32 * iter = cbands.data();
          for( unsigned int i = 0; i < height; ++i ) {
            for ( unsigned int j = 0; j < width; ++j ) {
              for ( unsigned int k = 0; k < components; ++k ) {
                const int value = *iter;
                ++iter;
                stream << value << " ";
              }
              stream << " "; // separate pixels with an extra space
            }
            stream << std::endl; // separate lines with a newline
          }

        }

        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        vector_type & cbands = static_cast< vector_type & >(bands);

        // write and store
        int x;
        for( unsigned int i = 0; i < width * components; ++i ) {
            x = cbands[i];
            stream << x << " ";
        }
    }

    void PnmEncoderImpl::write_bilevel_raw()
    {
        // cast the bands to the correct type
        typedef void_vector< UInt8 > vector_type;
        //vector_type & cbands = static_cast< vector_type & >(bands);

        // XXX
    }

    void PnmEncoderImpl::write_raw()
    {
        if ( pixeltype == "UINT8" ) {

            // cast the bands to the correct type
            typedef void_vector< UInt8 > vector_type;
            vector_type & cbands = static_cast< vector_type & >(bands);

            // write and store
            stream.write( reinterpret_cast< char * >(cbands.data()),
                          height * width * components );

        } else if ( pixeltype == "UINT16" ) {

            // cast the bands to the correct type
            typedef void_vector<UInt16> vector_type;
            vector_type & cbands = static_cast< vector_type & >(bands);

            // write and store
            byteorder bo( "big endian" );

            write_array( stream, bo, cbands.data(),
                         height * width * components );
        } else {
          vigra_postcondition( false, "internal error" );
        }
    }

    void PnmEncoder::nextScanline()
    {
        ++(pimpl->scanline);
    }

    void PnmEncoder::close()
    {
        if (!pimpl->bilevel) {

            // find out maxval and print it into the stream
            UInt32 maxval = 0;
            if ( pimpl->pixeltype == "UINT8" ) {
                void_vector< UInt8 > & cbands
                    = static_cast< void_vector< UInt8 > & >
                    (pimpl->bands);
                for( UInt8 * iter = cbands.begin();
                     iter < cbands.end(); ++iter )
                    if ( *iter > maxval ) maxval = *iter;
            } else if ( pimpl->pixeltype == "UINT16" ) {
                void_vector<UInt16> & cbands
                    = static_cast< void_vector<UInt16> & >(pimpl->bands);
                for( UInt16 * iter = cbands.begin();
                     iter < cbands.end(); ++iter )
                    if ( *iter > maxval ) maxval = *iter;
            } else if ( pimpl->pixeltype == "UINT32" ) {
                void_vector<UInt32> & cbands
                    = static_cast< void_vector<UInt32> & >(pimpl->bands);
                for( UInt32 * iter = cbands.begin();
                     iter < cbands.end(); ++iter )
                    if ( *iter > maxval ) maxval = *iter;
            }
            pimpl->stream << maxval << std::endl;

            // print the data
            if ( pimpl->raw )
                pimpl->write_raw();
            else
                pimpl->write_ascii();

        } else {

            // print the data
            if ( pimpl->raw )
                pimpl->write_bilevel_raw();
            else
                pimpl->write_bilevel_ascii();
        }
    }

    void PnmEncoder::abort() {}

} // namespace vigra
