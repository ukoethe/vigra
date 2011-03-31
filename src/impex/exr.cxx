/************************************************************************/
/*                                                                      */
/*               Copyright 2007 by Pablo d'Angelo                       */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifdef HasEXR

#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "void_vector.hxx"
#include "auto_file.hxx"
#include "exr.hxx"
#include "byteorder.hxx"
#include "error.hxx"
#include <stdexcept>
#include <iostream>

#include <ImfRgbaFile.h>
#include <ImfCRgbaFile.h>
#include <ImfStandardAttributes.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>

using namespace Imf;
using namespace Imath;

namespace vigra {

    CodecDesc ExrCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "EXR";

        // init pixel types
        desc.pixelTypes.resize(1);
        desc.pixelTypes[0] = "FLOAT";

        // init compression types
#if defined(IMF_B44_COMPRESSION) && defined(IMF_B44A_COMPRESSION)
        desc.compressionTypes.resize(7);
#else
        desc.compressionTypes.resize(5);
#endif
        desc.compressionTypes[0] = "NONE";
        desc.compressionTypes[1] = "ZIP";
        desc.compressionTypes[2] = "RLE";
        desc.compressionTypes[3] = "PIZ";
        desc.compressionTypes[4] = "PXR24";
#if defined(IMF_B44_COMPRESSION) && defined(IMF_B44A_COMPRESSION)
        desc.compressionTypes[5] = "B44";
        desc.compressionTypes[6] = "B44A";
#endif

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(4);
        desc.magicStrings[0][0] = '\x76';
        desc.magicStrings[0][1] = '\x2f';
        desc.magicStrings[0][2] = '\x31';
        desc.magicStrings[0][3] = '\x01';

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "exr";

        desc.bandNumbers.resize(1);
        desc.bandNumbers[0] = 4;

        return desc;
    }

    std::auto_ptr<Decoder> ExrCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new ExrDecoder() );
    }

    std::auto_ptr<Encoder> ExrCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new ExrEncoder() );
    }

    struct ExrDecoderImpl
    {
        std::string filename;
        // data source
        // auto_file file;

        // this is where libopenexr stores its state
        RgbaInputFile file;

        // data container
        ArrayVector<Rgba> pixels;
        ArrayVector<float> bands;

        // scanline counter
        int scanline;

        int width;
        int height;

        int components;
        int extra_components;
        Diff2D position;
        Size2D canvasSize;

        float x_resolution, y_resolution;

        // ctor, dtor
        ExrDecoderImpl( const std::string & filename );
        ~ExrDecoderImpl();

        // methods
        void init();
        void nextScanline();
    };

    ExrDecoderImpl::ExrDecoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
        : file( filename.c_str() ),
#else
        : file( filename.c_str() ),
#endif
          bands(0),
          scanline(-1), width(0), height(0),
          components(4), extra_components(1),
          x_resolution(0), y_resolution(0)
    {
    }

    ExrDecoderImpl::~ExrDecoderImpl()
    {
    }

    void ExrDecoderImpl::init()
    {
        // setup framebuffer
        Box2i dw = file.header().dataWindow();
        width  = dw.max.x - dw.min.x + 1;
        height = dw.max.y - dw.min.y + 1;

        position.x = dw.min.x;
        scanline = dw.min.y;
        position.y = dw.min.y;

        dw = file.header().displayWindow();
        canvasSize.x = dw.max.x+1;
        canvasSize.y = dw.max.y+1;

        // allocate data buffers
        pixels.resize(width);
        bands.resize(4*width);
    }

    void ExrDecoderImpl::nextScanline()
    {
        file.setFrameBuffer (pixels.data() - position.x - scanline * width, 1, width);
        file.readPixels (scanline, scanline);
        scanline++;
        // convert scanline to float
        float * dest = bands.begin();
        for (int i=0; i < width; i++) {
            *dest++ = pixels[i].r;
            *dest++ = pixels[i].g;
            *dest++ = pixels[i].b;
            *dest++ = pixels[i].a;
        }
    }

    void ExrDecoder::init( const std::string & filename )
    {
        pimpl = new ExrDecoderImpl(filename);
        pimpl->init();
    }

    ExrDecoder::~ExrDecoder()
    {
        delete pimpl;
    }

    std::string ExrDecoder::getFileType() const
    {
        return "EXR";
    }

    unsigned int ExrDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int ExrDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int ExrDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    unsigned int ExrDecoder::getNumExtraBands() const
    {
        return pimpl->extra_components;
    }

    float ExrDecoder::getXResolution() const
    {
        return pimpl->x_resolution;
    }

    float ExrDecoder::getYResolution() const
    {
        return pimpl->y_resolution;
    }

    Diff2D ExrDecoder::getPosition() const
    {
        return pimpl->position;
    }

    Size2D ExrDecoder::getCanvasSize() const
    {
        return pimpl->canvasSize;
    }

    std::string ExrDecoder::getPixelType() const
    {
        return "FLOAT";
    }

    unsigned int ExrDecoder::getOffset() const
    {
        return pimpl->components;
    }

    const void * ExrDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        return pimpl->bands.begin() + band;
    }

    void ExrDecoder::nextScanline()
    {
        pimpl->nextScanline();
    }

    void ExrDecoder::close() {}

    void ExrDecoder::abort() {}

    struct ExrEncoderImpl
    {
        std::string filename;
        // data sink
        RgbaOutputFile *file;

        // data container
        ArrayVector<float> bands;
        ArrayVector<Rgba> pixels;

        // image header fields
        int width, height, components;
        int extra_components;
        int bit_depth, color_type;
        Compression exrcomp;

        // scanline counter
        int scanline;

        // state
        bool finalized;

        // image layer position
        Diff2D position;
        Size2D canvasSize;

        // resolution
        float x_resolution, y_resolution;

        // ctor, dtor
        ExrEncoderImpl( const std::string & filename );
        ~ExrEncoderImpl();

        // methods
        void nextScanline();
        void finalize();
        void setCompressionType( const std::string &, int );
        void close();
    };

    ExrEncoderImpl::ExrEncoderImpl( const std::string & filename )
        : filename(filename), file(0), bands(0),
          exrcomp(PIZ_COMPRESSION), scanline(0), finalized(false),
          x_resolution(0), y_resolution(0)
    {
    }

    ExrEncoderImpl::~ExrEncoderImpl()
    {
        if (file)
            delete file;
    }

    void ExrEncoderImpl::finalize()
    {
        // prepare the bands
        bands.resize( 4 * width );
        pixels.resize(width);

        // set proper position
        Imath::Box2i displayWindow;
        if (canvasSize.x < width + position.x ||
            canvasSize.y < height + position.y)
        {
            displayWindow.min.x = 0;
            displayWindow.min.y = 0;
            displayWindow.max.x = width+position.x -1;
            displayWindow.max.y = height+position.y-1;
        } else {
            displayWindow.min.x = 0;
            displayWindow.min.y = 0;
            displayWindow.max.x = canvasSize.x -1;
            displayWindow.max.y = canvasSize.y -1;
        }
        Imath::Box2i dataWindow (Imath::V2i (position.x , position.y),
                                 Imath::V2i (width+position.x -1, height+position.y-1));
        Header header(displayWindow, dataWindow, 1, Imath::V2f(0, 0), 1, INCREASING_Y, exrcomp);
        file = new RgbaOutputFile(filename.c_str(), header, WRITE_RGBA);
        // enter finalized state
        finalized = true;
    }

    void ExrEncoderImpl::nextScanline()
    {
        // check if there are scanlines left at all, eventually write one
        if ( scanline < height ) {
            float * src = bands.data();
            for (int i=0; i < width; i++) {
                // convert to half
                pixels[i].r = *src++;
                pixels[i].g = *src++;
                pixels[i].b = *src++;
                pixels[i].a = *src++;
            }
            file->setFrameBuffer (pixels.begin() - position.x -(scanline+position.y)*width, 1, width);
            file->writePixels (1);
        }
        scanline++;
    }

    void ExrEncoderImpl::setCompressionType( const std::string & comp, int quality){
       if (comp == "NONE")
           exrcomp = NO_COMPRESSION;
       else if (comp == "ZIP")
           exrcomp = ZIP_COMPRESSION;
       else if (comp == "RLE")
           exrcomp = RLE_COMPRESSION;
       else if (comp == "PIZ")
           exrcomp = PIZ_COMPRESSION;
       else if (comp == "PXR24")
           exrcomp = PXR24_COMPRESSION;
#if defined(IMF_B44_COMPRESSION) && defined(IMF_B44A_COMPRESSION)
       else if (comp == "B44")
           exrcomp = B44_COMPRESSION;
       else if (comp == "B44A")
    	   exrcomp = B44A_COMPRESSION;
#endif
    }

    void ExrEncoderImpl::close()
    {
        delete file;
        file = 0;
    }

    void ExrEncoder::init( const std::string & filename )
    {
        pimpl = new ExrEncoderImpl(filename);
    }

    ExrEncoder::~ExrEncoder()
    {
        delete pimpl;
    }

    std::string ExrEncoder::getFileType() const
    {
        return "EXR";
    }

    void ExrEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->width = width;
    }

    void ExrEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->height = height;
    }

    void ExrEncoder::setNumBands( unsigned int bands )
    {
        if ( bands != 4 )
            vigra_fail( "internal error: number of components not supported." );
        pimpl->components = bands;
    }

    void ExrEncoder::setCompressionType( const std::string & comp,
                                         int quality )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->setCompressionType(comp, quality);
    }

    void ExrEncoder::setPosition( const Diff2D & pos )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->position = pos;
    }

    void ExrEncoder::setCanvasSize( const Size2D & size )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->canvasSize = size;
    }

    void ExrEncoder::setXResolution( float xres )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->x_resolution = xres;
    }

    void ExrEncoder::setYResolution( float yres )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->y_resolution = yres;
    }

    void ExrEncoder::setPixelType( const std::string & pixelType )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        if ( pixelType != "FLOAT" )
            vigra_fail( "internal error: pixeltype not supported." );
    }

    unsigned int ExrEncoder::getOffset() const
    {
        return pimpl->components;
    }

    void ExrEncoder::finalizeSettings()
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->finalize();
    }

    void * ExrEncoder::currentScanlineOfBand( unsigned int band )
    {
        return pimpl->bands.begin() + band;
    }

    void ExrEncoder::nextScanline()
    {
        // write scanline
        pimpl->nextScanline();
    }

    void ExrEncoder::close()
    {
        pimpl->close();
    }

    void ExrEncoder::abort() {}
}

#endif // HasEXR
