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
 *
 * Modification by Andrew Mihal, 27 October 2004:
 *  - Modified encoder to better estimate the number of rows per strip.
 *  - Modified decoder to use the scanline interface - the strip-based
 *    interface hogs memory when the rows/strip value is large.
 *  - Added support for ICC profiles
 * Andrew Mihal's modifications are covered by the VIGRA license.
 */

#ifdef HasTIFF

#ifdef _MSC_VER
// NB (jbeda): tiffio.h is going to include this anyway.  Let's include
// it now so that we can control how it comes in.  Namely, we want
// to get our version that doesn't set the evil min/max macros.
#include "vigra/windows.h"
#endif

#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "tiff.hxx"
#include <iostream>
#include <iomanip>
#include <sstream>

extern "C"
{

#include <tiff.h>
#include <tiffio.h>
#include <tiffvers.h>

static void vigraWarningHandler(char* module, char* fmt, va_list ap)
{
    static const std::string ignore("Unknown field with tag");
    if(ignore.compare(0, ignore.size(), fmt, ignore.size()) == 0)
        return;
    
    if (module != NULL)
    {
        static const std::string ignore("TIFFFetchNormalTag");
        if(ignore.compare(module) == 0)
            return;
            
        fprintf(stderr, "%s: ", module);
    }
    
    fprintf(stderr, "Warning, ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, ".\n");
}

} // extern "C"

namespace vigra {

    CodecDesc TIFFCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "TIFF";

        // init pixel types
        desc.pixelTypes.resize(9);
        desc.pixelTypes[0] = "BILEVEL";
        desc.pixelTypes[1] = "UINT8";
        desc.pixelTypes[2] = "INT8";
        desc.pixelTypes[3] = "UINT16";
        desc.pixelTypes[4] = "INT16";
        desc.pixelTypes[5] = "UINT32";
        desc.pixelTypes[6] = "INT32";
        desc.pixelTypes[7] = "FLOAT";
        desc.pixelTypes[8] = "DOUBLE";

        // init compression types
        desc.compressionTypes.resize(6);
        desc.compressionTypes[0] = "NONE";
        desc.compressionTypes[1] = "RLE";
        desc.compressionTypes[2] = "PACKBITS";
        desc.compressionTypes[3] = "JPEG";
        desc.compressionTypes[4] = "LZW";
        desc.compressionTypes[5] = "DEFLATE";

        // init magic strings
#if TIFFLIB_VERSION > 20070712
        desc.magicStrings.resize(3);
#else
        desc.magicStrings.resize(2);
#endif
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = '\115';
        desc.magicStrings[0][1] = '\115';
        desc.magicStrings[0][2] = '\000';
        desc.magicStrings[0][3] = '\052';
        desc.magicStrings[1].resize(4);
        desc.magicStrings[1][0] = '\111';
        desc.magicStrings[1][1] = '\111';
        desc.magicStrings[1][2] = '\052';
        desc.magicStrings[1][3] = '\000';
#if TIFFLIB_VERSION > 20070712
        // magic for bigtiff
        desc.magicStrings[2].resize(4);
        desc.magicStrings[2][0] = '\111';
        desc.magicStrings[2][1] = '\111';
        desc.magicStrings[2][2] = '\053';
        desc.magicStrings[2][3] = '\000';
#endif

        // init file extensions
        desc.fileExtensions.resize(2);
        desc.fileExtensions[0] = "tif";
        desc.fileExtensions[1] = "tiff";

        desc.bandNumbers.resize(4);
        desc.bandNumbers[0] = 1;
        desc.bandNumbers[1] = 2;
        desc.bandNumbers[2] = 3;
        desc.bandNumbers[3] = 4;

        return desc;
    }

    VIGRA_UNIQUE_PTR<Decoder> TIFFCodecFactory::getDecoder() const
    {
        // use 'NULL' to silence all warnings
        TIFFSetWarningHandler((TIFFErrorHandler)&vigraWarningHandler);

        return VIGRA_UNIQUE_PTR<Decoder>( new TIFFDecoder() );
    }

    VIGRA_UNIQUE_PTR<Encoder> TIFFCodecFactory::getEncoder() const
    {
        return VIGRA_UNIQUE_PTR<Encoder>( new TIFFEncoder() );
    }

    class TIFFCodecImpl
    {

    protected:

        std::string pixeltype;

        TIFF * tiff;
        tdata_t * stripbuffer;
        tstrip_t strip;

        uint32 stripindex, stripheight;
        uint32 width, height;
        uint16 samples_per_pixel, bits_per_sample,
            photometric, planarconfig, fillorder, extra_samples_per_pixel;
        float x_resolution, y_resolution;
        Diff2D position;
        Size2D canvasSize;

        Decoder::ICCProfile iccProfile;

    public:

        TIFFCodecImpl();
        ~TIFFCodecImpl();
    };

    TIFFCodecImpl::TIFFCodecImpl()
        : pixeltype("undefined")
    {
        tiff = 0;
        stripbuffer = 0;
        strip = 0;
        stripindex = 0;
        planarconfig = PLANARCONFIG_CONTIG;
        x_resolution = 0;
        y_resolution = 0;
        extra_samples_per_pixel = 0;
   }

    TIFFCodecImpl::~TIFFCodecImpl()
    {
        if ( planarconfig == PLANARCONFIG_SEPARATE ) {
            if ( stripbuffer != 0 ) {
                for( unsigned int i = 0; i < samples_per_pixel; ++i )
                    if ( stripbuffer[i] != 0 )
                        _TIFFfree(stripbuffer[i]);
                delete[] stripbuffer;
            }
        } else {
            if ( stripbuffer != 0 ) {
                if ( stripbuffer[0] != 0 )
                    _TIFFfree(stripbuffer[0]);
                delete[] stripbuffer;
            }
        }

        if ( tiff != 0 )
            TIFFClose(tiff);
    }

    class TIFFDecoderImpl : public TIFFCodecImpl
    {
        friend class TIFFDecoder;

        unsigned int scanline;

        std::string get_pixeltype_by_sampleformat() const;
        std::string get_pixeltype_by_datatype() const;

    public:

        TIFFDecoderImpl( const std::string & filename );

        void init( unsigned int imageIndex );

        unsigned int getNumImages();
        void setImageIndex( unsigned int index );
        unsigned int getImageIndex();

        const void * currentScanlineOfBand( unsigned int band ) const;
        void nextScanline();
    };

    TIFFDecoderImpl::TIFFDecoderImpl( const std::string & filename )
    {
        tiff = TIFFOpen( filename.c_str(), "r" );

        if ( !tiff ) {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }

        scanline = 0;
    }

    std::string TIFFDecoderImpl::get_pixeltype_by_sampleformat() const
    {
        uint16 sampleformat;

        if ( TIFFGetField( tiff, TIFFTAG_SAMPLEFORMAT, &sampleformat ) ) {

            switch (sampleformat) {

            case SAMPLEFORMAT_UINT:
                // added by dangelo, support for UINT 16 & 32 bit
                switch (bits_per_sample) {
                case 8:
                    return "UINT8";
                case 16:
                    return "UINT16";
                case 32:
                    return "UINT32";
                }
                break;

            case SAMPLEFORMAT_INT:
                switch (bits_per_sample) {
                case 8:
                    return "INT8";
                case 16:
                    return "INT16";
                case 32:
                    return "INT32";
                }
                break;

            case SAMPLEFORMAT_IEEEFP:
                switch (bits_per_sample) {
                case 32:
                    return "FLOAT";
                case 64:
                    return "DOUBLE";
                }
                break;
            }
        }
        return "undefined";
    }

    std::string TIFFDecoderImpl::get_pixeltype_by_datatype() const
    {
        uint16 datatype;

        if ( TIFFGetField( tiff, TIFFTAG_DATATYPE, &datatype ) ) {
            // dangelo: correct parsing of INT/UINT (given in tiff.h)
            switch (datatype) {
            case TIFF_BYTE:
                return "UINT8";
            case TIFF_SBYTE:
                return "INT8";
            case TIFF_SHORT:
                return "UINT16";
            case TIFF_SSHORT:
                return "INT16";
            case TIFF_LONG:
                return "UINT32";
            case TIFF_SLONG:
                return "INT32";
            case TIFF_FLOAT:
                return "FLOAT";
            case TIFF_DOUBLE:
                return "DOUBLE";
            }
        }

        return "undefined";
    }

    void TIFFDecoderImpl::init(unsigned int imageIndex)
    {
        // set image directory, if necessary:
        if (imageIndex != TIFFCurrentDirectory(tiff))
        {
            if (!TIFFSetDirectory( tiff, (tdir_t)(imageIndex)) )
                vigra_fail( "Invalid TIFF image index" );
        }

        // read width and height
        TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &width );
        TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &height );

        // check for tiled TIFFs
        uint32 tileWidth, tileHeight;
        if( TIFFGetField( tiff, TIFFTAG_TILEWIDTH, &tileWidth ) &&
            TIFFGetField( tiff, TIFFTAG_TILELENGTH, &tileHeight ) )
            vigra_precondition( (tileWidth == width) && (tileHeight == height),
                                "TIFFDecoderImpl::init(): "
                                "Cannot read tiled TIFFs (not implemented)." );

        // find out strip heights
        stripheight = 1; // now using scanline interface instead of strip interface

        // get samples_per_pixel
        samples_per_pixel = 0;
        extra_samples_per_pixel = 0;
        if ( !TIFFGetFieldDefaulted( tiff, TIFFTAG_SAMPLESPERPIXEL,
                                     &samples_per_pixel ) )
            vigra_fail( "TIFFDecoderImpl::init(): Samples per pixel not set."
                        " A suitable default was not found." );

        // read extra samples (# of alpha channels)
        uint16 *extra_sample_types=0;
        if (TIFFGetField( tiff, TIFFTAG_EXTRASAMPLES,
                          &extra_samples_per_pixel, &extra_sample_types )!=1)
        {
            extra_samples_per_pixel=0;
        }

        if (extra_samples_per_pixel > 0) {
            for (int i=0; i< extra_samples_per_pixel; i++) {
                if (extra_sample_types[i] ==  EXTRASAMPLE_ASSOCALPHA)
                {
                    std::cerr << "WARNING: TIFFDecoderImpl::init(): associated alpha treated"
                                 " as unassociated alpha!" << std::endl;
                }
            }
        }

        // get photometric
        if ( !TIFFGetFieldDefaulted( tiff, TIFFTAG_PHOTOMETRIC,
                                     &photometric ) )
            vigra_fail( "TIFFDecoderImpl::init(): Photometric tag is not set."
                        " A suitable default was not found." );

        // check photometric preconditions
        switch ( photometric )
        {
            case PHOTOMETRIC_MINISWHITE:
            case PHOTOMETRIC_MINISBLACK:
            case PHOTOMETRIC_PALETTE:
            {
                if ( samples_per_pixel - extra_samples_per_pixel != 1 )
                vigra_fail("TIFFDecoderImpl::init():"
                                " Photometric tag does not fit the number of"
                                " samples per pixel." );
                break;
            }
            case PHOTOMETRIC_RGB:
            {
                if ( samples_per_pixel > 3 && extra_samples_per_pixel == 0 ) {
                    // file probably lacks the extra_samples tag
                    extra_samples_per_pixel = samples_per_pixel - 3;
                }
                if ( samples_per_pixel - extra_samples_per_pixel != 3 )
                    vigra_fail("TIFFDecoderImpl::init():"
                                    " Photometric tag does not fit the number of"
                                    " samples per pixel." );
                break;
            }
            case PHOTOMETRIC_LOGL:
            case PHOTOMETRIC_LOGLUV:
            {
                uint16 tiffcomp;
                TIFFGetFieldDefaulted( tiff, TIFFTAG_COMPRESSION, &tiffcomp );
                if (tiffcomp != COMPRESSION_SGILOG && tiffcomp != COMPRESSION_SGILOG24)
                    vigra_fail("TIFFDecoderImpl::init():"
                                    " Only SGILOG compression is supported for"
                                    " LogLuv TIFF."
                    );
                TIFFSetField(tiff, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
                break;
            }
        }

        // get planarconfig
        if ( samples_per_pixel > 1 ) {
            if ( !TIFFGetFieldDefaulted( tiff, TIFFTAG_PLANARCONFIG,
                                         &planarconfig ) )
                vigra_fail( "TIFFDecoderImpl::init(): Planarconfig is not"
                            " set. A suitable default was not found." );
        }

        // get bits per pixel
        if ( !TIFFGetField( tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample ) )
        {
            std::cerr << "Warning: no TIFFTAG_BITSPERSAMPLE, using 8 bits per sample.\n";
            bits_per_sample = 8;
        }
        // get pixeltype
        if ( bits_per_sample != 1 ) {

            // try the sampleformat tag
            pixeltype = get_pixeltype_by_sampleformat();

            if ( pixeltype == "undefined" ) {

                // try the (obsolete) datatype tag
                pixeltype = get_pixeltype_by_datatype();

                if ( pixeltype == "undefined" ) {
                    // ERROR: no useable pixeltype found..
                    // imagemagick can write files without it..
                    // try to guess a suitable one here.
                    switch(bits_per_sample)
                    {
                        case 8:
                            pixeltype = "UINT8";
                            break;
                        case 16:
                            pixeltype = "UINT16";
                            break;
                        case 32:
                            pixeltype = "UINT32"; // prefer int over float
                            break;
                        case 64:
                            pixeltype =  "DOUBLE";
                            break;
                        default:
                            vigra_fail( "TIFFDecoderImpl::init(): Sampleformat or Datatype tag undefined and guessing sampletype from Bits per Sample failed." );
                            break;
                    }
                    if(bits_per_sample != 8) // issue the warning only for non-trivial cases
                        std::cerr << "Warning: no TIFFTAG_SAMPLEFORMAT or TIFFTAG_DATATYPE, "
                                     "guessing pixeltype '" << pixeltype << "'.\n";
                }
            }

        } else {

            // if each sample is 1 bit long
            pixeltype = "BILEVEL";

            // get fillorder
            if ( !TIFFGetField( tiff, TIFFTAG_FILLORDER, &fillorder ) )
                fillorder = FILLORDER_MSB2LSB;
        }

        // make sure the LogLuv has correct pixeltype because only float is supported
        if (photometric == PHOTOMETRIC_LOGLUV) {
            pixeltype = "FLOAT";
            samples_per_pixel = 3;
        }

        // other fields
        uint16 u16value;
        uint32 u32value;
        float unitLength = 1.0f;
        if (TIFFGetField( tiff, TIFFTAG_RESOLUTIONUNIT, &u16value )) {
            switch (u16value) {
            case RESUNIT_NONE:
                unitLength = 1.0f;
                break;
            case RESUNIT_INCH:
                unitLength = 1.0f;
                break;
            case RESUNIT_CENTIMETER:
                unitLength = 1.0f/2.54f;
                break;
            default:
                vigra_fail("Unknown resolution unit");
            }
        }

        float fvalue;
        if (TIFFGetField( tiff, TIFFTAG_XRESOLUTION, &fvalue )) {
            x_resolution = fvalue / unitLength;
        }

        if (TIFFGetField( tiff, TIFFTAG_YRESOLUTION, &fvalue )) {
            y_resolution = fvalue / unitLength;
        }

        // XPosition
        if (TIFFGetField( tiff, TIFFTAG_XPOSITION, &fvalue )) {
            fvalue = fvalue * x_resolution;
            position.x = (int)floor(fvalue + 0.5);
        }
        // YPosition
        if (TIFFGetField( tiff, TIFFTAG_YPOSITION, &fvalue )) {
            fvalue = fvalue * y_resolution;
            position.y = (int)floor(fvalue + 0.5);
        }

        // canvas size
        if (TIFFGetField( tiff, TIFFTAG_PIXAR_IMAGEFULLWIDTH, &u32value )) {
            canvasSize.x = u32value;
        }
        if (TIFFGetField( tiff, TIFFTAG_PIXAR_IMAGEFULLLENGTH, &u32value )) {
            canvasSize.y = u32value;
        }

        if ((uint32)canvasSize.x < position.x + width || (uint32)canvasSize.y < position.y + height)
        {
            canvasSize.x = canvasSize.y = 0;
        }

        // ICC Profile
        UInt32 iccProfileLength = 0;
        const unsigned char *iccProfilePtr = NULL;
        if(TIFFGetField(tiff, TIFFTAG_ICCPROFILE,
                        &iccProfileLength, &iccProfilePtr)
           && iccProfileLength)
        {
            Decoder::ICCProfile iccData(
                iccProfilePtr, iccProfilePtr + iccProfileLength);
            iccProfile.swap(iccData);
        }

        // allocate data buffers
        const unsigned int stripsize = TIFFScanlineSize(tiff);
        if ( planarconfig == PLANARCONFIG_SEPARATE ) {
            stripbuffer = new tdata_t[samples_per_pixel];
            for( unsigned int i = 0; i < samples_per_pixel; ++i ) {
                stripbuffer[i] = 0;
            }
            for( unsigned int i = 0; i < samples_per_pixel; ++i ) {
                stripbuffer[i] = _TIFFmalloc(stripsize);
                if(stripbuffer[i] == 0)
                    throw std::bad_alloc();
            }
        } else {
            stripbuffer = new tdata_t[1];
            stripbuffer[0] = 0;
            stripbuffer[0] = _TIFFmalloc(stripsize < width ? width : stripsize);
            if(stripbuffer[0] == 0)
                throw std::bad_alloc();
        }

        // let the codec read a new strip
        stripindex = stripheight;
    }

    const void *
    TIFFDecoderImpl::currentScanlineOfBand( unsigned int band ) const
    {
        if ( bits_per_sample == 1 ) {
            const unsigned int n = TIFFScanlineSize(tiff);
            UInt8 * const startpointer
                = static_cast< UInt8 * >(stripbuffer[0]);
            UInt8 * bytepointer = startpointer;

            bytepointer += n-1;
            for (int byte = n-1 ; byte >= 0; --byte)
            {
                UInt8 currentByte = *bytepointer;
                --bytepointer;

                UInt8 * bitpointer = startpointer;
                bitpointer += byte * 8;

                for (unsigned char bit = 7; bit < 8; --bit)
                {
                    *bitpointer = ((currentByte & (1 << bit)) ? photometric : 1 - photometric);
                    ++bitpointer;
                    if (byte * 8 + 7 - bit == width - 1) break;
                }
            }
            // XXX probably right
            return startpointer + ( stripindex * width ) / 8;
        } else {
            if ( planarconfig == PLANARCONFIG_SEPARATE ) {
                UInt8 * const buf
                    = static_cast< UInt8 * >(stripbuffer[band]);
                return buf + ( stripindex * width ) * ( bits_per_sample / 8 );
            } else {
                UInt8 * const buf
                    = static_cast< UInt8 * >(stripbuffer[0]);
                return buf + ( band + stripindex * width * samples_per_pixel )
                    * ( bits_per_sample / 8 );
            }
        }
    }

    void TIFFDecoderImpl::nextScanline()
    {
        // eventually read a new strip
        if ( ++stripindex >= stripheight ) {
            stripindex = 0;

            if ( planarconfig == PLANARCONFIG_SEPARATE ) {
                const tsize_t size = TIFFScanlineSize(tiff);
                for( unsigned int i = 0; i < samples_per_pixel; ++i )
                    TIFFReadScanline(tiff, stripbuffer[i], scanline++, size);
            } else {
                TIFFReadScanline( tiff, stripbuffer[0], scanline++, 0);
            }

            // XXX handle bilevel images

            // invert grayscale images that interpret 0 as white
            if ( photometric == PHOTOMETRIC_MINISWHITE &&
                 samples_per_pixel == 1 && pixeltype == "UINT8" ) {

                UInt8 * buf = static_cast< UInt8 * >(stripbuffer[0]);
                const unsigned int n = TIFFScanlineSize(tiff);

                // invert every pixel
                for ( unsigned int i = 0; i < n; ++i, ++buf )
                    *buf = 0xff - *buf;
            }
        }
    }

    void TIFFDecoder::init( const std::string & filename, unsigned int imageIndex=0 )
    {
        pimpl = new TIFFDecoderImpl(filename);
        pimpl->init(imageIndex);
        iccProfile_ = pimpl->iccProfile;
    }

    TIFFDecoder::~TIFFDecoder()
    {
        delete pimpl;
    }

    std::string TIFFDecoder::getFileType() const
    {
        return "TIFF";
    }

    unsigned int TIFFDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int TIFFDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int TIFFDecoder::getNumBands() const
    {
        return pimpl->samples_per_pixel;
    }

    unsigned int TIFFDecoder::getNumExtraBands() const
    {
        return pimpl->extra_samples_per_pixel;
    }

    unsigned int TIFFDecoder::getNumImages() const
    {
        return pimpl->getNumImages();
    }

    void TIFFDecoder::setImageIndex(unsigned int imageIndex)
    {
        pimpl->setImageIndex(imageIndex);
    }

    unsigned int TIFFDecoder::getImageIndex() const
    {
        return pimpl->getImageIndex();
    }

    vigra::Diff2D TIFFDecoder::getPosition() const
    {
        return pimpl->position;
    }

    vigra::Size2D TIFFDecoder::getCanvasSize() const
    {
        return pimpl->canvasSize;
    }

    float TIFFDecoder::getXResolution() const
    {
        return pimpl->x_resolution;
    }

    float TIFFDecoder::getYResolution() const
    {
        return pimpl->y_resolution;
    }

    std::string TIFFDecoder::getPixelType() const
    {
        return pimpl->pixeltype;
    }

    unsigned int TIFFDecoder::getOffset() const
    {
        return pimpl->planarconfig == PLANARCONFIG_SEPARATE ?
            1 : pimpl->samples_per_pixel;
    }

    unsigned int
    TIFFDecoderImpl::getNumImages()
    {
        unsigned int currIndex = getImageIndex();
        TIFFSetDirectory(tiff, 0);
        int numPages = 1;
        while (TIFFReadDirectory(tiff)) numPages++;
        TIFFSetDirectory(tiff, currIndex);
        return numPages;
    }

    void
    TIFFDecoderImpl::setImageIndex( unsigned int imageIndex )
    {
        init(imageIndex);
    }

    unsigned int
    TIFFDecoderImpl::getImageIndex()
    {
        return TIFFCurrentDirectory(tiff);
    }

    const void * TIFFDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->currentScanlineOfBand(band);
    }

    void TIFFDecoder::nextScanline()
    {
        pimpl->nextScanline();
    }

    void TIFFDecoder::close() {}
    void TIFFDecoder::abort() {}

    // this encoder always writes interleaved tiff files
    class TIFFEncoderImpl : public TIFFCodecImpl
    {
        friend class TIFFEncoder;

        // attributes

        unsigned short tiffcomp;
        bool finalized;

    public:

        // ctor, dtor

        TIFFEncoderImpl( const std::string & filename, const std::string & mode )
            : tiffcomp(COMPRESSION_LZW), finalized(false)
        {
            tiff = TIFFOpen( filename.c_str(), mode.c_str() );
            if (!tiff)
            {
                std::string msg("Unable to open file '");
                msg += filename;
                msg += "'.";
                vigra_precondition( false, msg.c_str() );
            }

            planarconfig = PLANARCONFIG_CONTIG;
        }

        // methods

        void setCompressionType( const std::string &, int );
        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int band ) const
        {
            const unsigned int atomicbytes = bits_per_sample >> 3;
            UInt8 * buf = ( UInt8 * ) stripbuffer[0];
            return buf + atomicbytes *
                ( width * samples_per_pixel * stripindex + band );
        }

        void nextScanline()
        {
            // compute the number of rows in the current strip
            unsigned int rows = ( strip + 1 ) * stripheight > height ?
                height - strip * stripheight : stripheight;

            if ( ++stripindex >= rows ) {

                // write next strip
                stripindex = 0;

                int success = TIFFWriteEncodedStrip( tiff, strip++, stripbuffer[0],
                                       TIFFVStripSize( tiff, rows ) );
                if(success == -1 && tiffcomp != COMPRESSION_NONE)
                {
                    throw Encoder::TIFFCompressionException(); // retry without compression
                }

                vigra_postcondition(success != -1,
                        "exportImage(): Unable to write TIFF data.");
            }
        }
    };

    void TIFFEncoderImpl::setCompressionType( const std::string & comp,
                                              int quality = -1 )
    {
        // if any compression type is set that we do not support,
        // the expected behavior is to do nothing
        if ( comp == "NONE" )
            tiffcomp = COMPRESSION_NONE;
        else if ( ( comp == "JPEG" ) && ( quality != -1 ) )
            tiffcomp = COMPRESSION_OJPEG;
        else if ( comp == "RLE" || comp == "RunLength")
            tiffcomp = COMPRESSION_CCITTRLE;
        else if ( comp == "PACKBITS")
            tiffcomp = COMPRESSION_PACKBITS;
        else if ( comp == "LZW" )
            tiffcomp = COMPRESSION_LZW;
        else if ( comp == "DEFLATE" )
            tiffcomp = COMPRESSION_DEFLATE;
    }

    void TIFFEncoderImpl::finalizeSettings()
    {
        // decide if we should write Grey, or RGB files
        // all additional channels are treated as extra samples
        if (samples_per_pixel < 3) {
            extra_samples_per_pixel = samples_per_pixel - 1;
        } else {
            extra_samples_per_pixel = samples_per_pixel - 3;
        }

        // set some fields
        TIFFSetField( tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
        TIFFSetField( tiff, TIFFTAG_IMAGEWIDTH, width );
        TIFFSetField( tiff, TIFFTAG_IMAGELENGTH, height );
        // TIFFDefaultStripSize tries for 8kb strips! Laughable!
        // This will do a 1MB strip for 8-bit images,
        // 2MB strip for 16-bit, and so forth.
        unsigned int estimate =
            (unsigned int)std::max(static_cast<UIntBiggest>(1),
                                  (static_cast<UIntBiggest>(1)<<20) / (width * samples_per_pixel));
        TIFFSetField( tiff, TIFFTAG_ROWSPERSTRIP,
                      stripheight = TIFFDefaultStripSize( tiff, estimate ) );
        TIFFSetField( tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel );
        TIFFSetField( tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
        TIFFSetField( tiff, TIFFTAG_COMPRESSION, tiffcomp );

        // subfile descriptor
        TIFFSetField( tiff, TIFFTAG_SUBFILETYPE, 0);

        // set pixel type
        if ( pixeltype == "BILEVEL" ) {
            // no sampleformat applies for bilevel
            bits_per_sample = 1;
        } else if ( pixeltype == "UINT8" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT );
            bits_per_sample = 8;
        } else if ( pixeltype == "INT16" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT );
            bits_per_sample = 16;
        } else if ( pixeltype == "UINT16" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT );
            bits_per_sample = 16;
        } else if ( pixeltype == "INT32" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT );
            bits_per_sample = 32;
        } else if ( pixeltype == "UINT32" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT );
            bits_per_sample = 32;
        } else if ( pixeltype == "FLOAT" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
            bits_per_sample = 32;
        } else if ( pixeltype == "DOUBLE" ) {
            TIFFSetField( tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
            bits_per_sample = 64;
        }
        TIFFSetField( tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample );

       if (extra_samples_per_pixel > 0) {
              uint16 * types = new  uint16[extra_samples_per_pixel];
           for ( int i=0; i < extra_samples_per_pixel; i++ ) {
              types[i] = EXTRASAMPLE_UNASSALPHA;
            }
            TIFFSetField( tiff, TIFFTAG_EXTRASAMPLES, extra_samples_per_pixel,
                          types );
                delete[] types;
        }

        // set photometric
        if ( samples_per_pixel - extra_samples_per_pixel == 1 )
            TIFFSetField( tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
        else if ( samples_per_pixel - extra_samples_per_pixel == 3 )
            TIFFSetField( tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB );

        // set resolution
        if ( x_resolution > 0) {
            TIFFSetField( tiff, TIFFTAG_XRESOLUTION, x_resolution );
        }
        if ( y_resolution > 0 ) {
            TIFFSetField( tiff, TIFFTAG_YRESOLUTION, y_resolution );
        }
        if (x_resolution > 0 || y_resolution > 0) {
            TIFFSetField( tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH );
        }

        // save position, if available
        if (position.x >= 0 && position.y >= 0 &&
            x_resolution > 0 && y_resolution > 0)
        {
            TIFFSetField( tiff, TIFFTAG_XPOSITION, position.x / x_resolution);
            TIFFSetField( tiff, TIFFTAG_YPOSITION, position.y / y_resolution);
        }

        if ((uint32)canvasSize.x >= position.x + width
            && (uint32)canvasSize.y >= position.y + height)
        {
            TIFFSetField( tiff, TIFFTAG_PIXAR_IMAGEFULLWIDTH, canvasSize.x);
            TIFFSetField( tiff, TIFFTAG_PIXAR_IMAGEFULLLENGTH, canvasSize.y);
        }

        // Set ICC profile, if available.
        if (iccProfile.size()) {
            TIFFSetField(tiff, TIFFTAG_ICCPROFILE,
                         iccProfile.size(), iccProfile.begin());
        }

        // alloc memory
        stripbuffer = new tdata_t[1];
        stripbuffer[0] = 0;
        stripbuffer[0] = _TIFFmalloc( TIFFStripSize(tiff) );
        if(stripbuffer[0] == 0)
            throw std::bad_alloc();

        finalized = true;
    }

    void TIFFEncoder::init( const std::string & filename, const std::string & mode )
    {
        pimpl = new TIFFEncoderImpl(filename, mode);
    }

    TIFFEncoder::~TIFFEncoder()
    {
        delete pimpl;
    }

    std::string TIFFEncoder::getFileType() const
    {
        return "TIFF";
    }

    void TIFFEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void TIFFEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void TIFFEncoder::setNumBands( unsigned int bands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->samples_per_pixel = bands;
    }

    void TIFFEncoder::setCompressionType( const std::string & comp,
                                          int quality )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->setCompressionType( comp, quality );
    }

    void TIFFEncoder::setPixelType( const std::string & pixeltype )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->pixeltype = pixeltype;
    }

    void TIFFEncoder::setPosition( const vigra::Diff2D & pos )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->position = pos;
    }

    void TIFFEncoder::setCanvasSize( const vigra::Size2D & size )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->canvasSize = size;
    }

    void TIFFEncoder::setXResolution( float xres )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->x_resolution = xres;
    }

    void TIFFEncoder::setYResolution( float yres )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->y_resolution = yres;
    }

    unsigned int TIFFEncoder::getOffset() const
    {
        return pimpl->samples_per_pixel;
    }

    void TIFFEncoder::finalizeSettings()
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->finalizeSettings();
    }

    void * TIFFEncoder::currentScanlineOfBand( unsigned int band )
    {
        return pimpl->currentScanlineOfBand(band);
    }

    void TIFFEncoder::nextScanline()
    {
        pimpl->nextScanline();
    }

    void TIFFEncoder::setICCProfile(const ICCProfile & data)
    {
        pimpl->iccProfile = data;
    }

    void TIFFEncoder::close() {}
    void TIFFEncoder::abort() {}
}

#endif // HasTIFF
