/************************************************************************/
/*                                                                      */
/*               Copyright 2007 by Pablo d'Angelo                       */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.5.0, Dec 07 2006 )                                    */
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

#ifndef VIGRA_IMPEX_EXR_HXX
#define VIGRA_IMPEX_EXR_HXX

#include "vigra/codec.hxx"

// EXR - OpenEXR

namespace vigra {

    struct ExrDecoderImpl;
    struct ExrEncoderImpl;

    struct ExrCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        VIGRA_UNIQUE_PTR<Decoder> getDecoder() const;
        VIGRA_UNIQUE_PTR<Encoder> getEncoder() const;
    };

    class ExrDecoder : public Decoder
    {
        ExrDecoderImpl * pimpl;

    public:

        ExrDecoder() : pimpl(0) {}

        ~ExrDecoder();

        void init( const std::string & );
        void close();
        void abort();

        std::string getFileType() const;
        std::string getPixelType() const;

        unsigned int getWidth() const;
        unsigned int getHeight() const;
        unsigned int getNumBands() const;
        unsigned int getNumExtraBands() const;
        float getXResolution() const;
        float getYResolution() const;
        Diff2D getPosition() const;
        Size2D getCanvasSize() const;

        unsigned int getOffset() const;

        const void * currentScanlineOfBand( unsigned int ) const;
        void nextScanline();
    };

    class ExrEncoder : public Encoder
    {
        ExrEncoderImpl * pimpl;

    public:

        ExrEncoder() : pimpl(0) {}

        ~ExrEncoder();

        void init( const std::string & );
        void close();
        void abort();

        std::string getFileType() const;
        unsigned int getOffset() const;

        void setWidth( unsigned int );
        void setHeight( unsigned int );
        void setNumBands( unsigned int );
        void setCompressionType( const std::string &, int = -1 );
        void setPixelType( const std::string & );

        void setPosition( const Diff2D & pos );
        void setCanvasSize( const Size2D & pos );
        void setXResolution( float xres );
        void setYResolution( float yres );

        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int );
        void nextScanline();
    };
}

#endif // VIGRA_IMPEX_EXR_HXX
