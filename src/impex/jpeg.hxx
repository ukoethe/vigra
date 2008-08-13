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

#ifndef VIGRA_IMPEX_JPEG_HXX
#define VIGRA_IMPEX_JPEG_HXX

#include <vector>
#include "vigra/codec.hxx"

namespace vigra {

    struct JPEGCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        std::auto_ptr<Decoder> getDecoder() const;
        std::auto_ptr<Encoder> getEncoder() const;
    };

    struct JPEGDecoderImpl;
    struct JPEGEncoderImpl;

    class JPEGDecoder : public Decoder
    {
        JPEGDecoderImpl * pimpl;

    public:

        JPEGDecoder() : pimpl(0) {}

        ~JPEGDecoder();

        std::string getFileType() const;
        unsigned int getWidth() const;
        unsigned int getHeight() const;
        unsigned int getNumBands() const;

        const void * currentScanlineOfBand( unsigned int ) const;
        void nextScanline();

        std::string getPixelType() const;
        unsigned int getOffset() const;

        void init( const std::string & );
        void close();
        void abort();

    };

    class JPEGEncoder : public Encoder
    {
        JPEGEncoderImpl * pimpl;

    public:

        JPEGEncoder() : pimpl(0) {}

        ~JPEGEncoder();

        std::string getFileType() const;
        void setWidth( unsigned int );
        void setHeight( unsigned int );
        void setNumBands( unsigned int );

        void setICCProfile(const ICCProfile & data);

        void setCompressionType( const std::string &, int = -1 );
        void setPixelType( const std::string & );
        unsigned int getOffset() const;

        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int );
        void nextScanline();

        void init( const std::string & );
        void close();
        void abort();
    };
}

#endif // VIGRA_IMPEX_JPEG_HXX
