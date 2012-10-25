/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Pablo d'Angelo                  */
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


#ifndef VIGRA_IMPEX_HDR_HXX
#define VIGRA_IMPEX_HDR_HXX

#include <vector>
#include "vigra/diff2d.hxx"
#include "vigra/codec.hxx"

namespace vigra {

    struct HDRCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        VIGRA_UNIQUE_PTR<Decoder> getDecoder() const;
        VIGRA_UNIQUE_PTR<Encoder> getEncoder() const;
    };

    class HDRDecoderImpl;
    class HDREncoderImpl;

    class HDRDecoder : public Decoder
    {
        HDRDecoderImpl * pimpl;

    public:

        HDRDecoder() : pimpl(0) {}

        ~HDRDecoder();

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

    class HDREncoder : public Encoder
    {
        HDREncoderImpl * pimpl;

    public:

        HDREncoder() : pimpl(0) {}

        ~HDREncoder();

        std::string getFileType() const;
        void setWidth( unsigned int );
        void setHeight( unsigned int );
        void setNumBands( unsigned int );

        void setCompressionType( const std::string &, int = -1 );
        void setPixelType( const std::string & );

        void setPosition( const vigra::Diff2D & pos );
        void setXResolution( float xres );
        void setYResolution( float yres );

        unsigned int getOffset() const;

        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int );
        void nextScanline();

        void init( const std::string & );
        void close();
        void abort();
    };
}

#endif // VIGRA_IMPEX_HDR_HXX
