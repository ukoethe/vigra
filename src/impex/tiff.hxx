/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Gunnar Kedenburg                */
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

#ifndef VIGRA_IMPEX_TIFF_HXX
#define VIGRA_IMPEX_TIFF_HXX

#include <vector>
#include "vigra/codec.hxx"

namespace vigra {

    struct TIFFCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        std::auto_ptr<Decoder> getDecoder() const;
        std::auto_ptr<Encoder> getEncoder() const;
    };

    class TIFFDecoderImpl;
    class TIFFEncoderImpl;

    class TIFFDecoder : public Decoder
    {
        TIFFDecoderImpl * pimpl;

    public:

        TIFFDecoder() : pimpl(0) {}

        ~TIFFDecoder();

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

    class TIFFEncoder : public Encoder
    {
        TIFFEncoderImpl * pimpl;

    public:

        TIFFEncoder() : pimpl(0) {}

        ~TIFFEncoder();

        std::string getFileType() const;
        void setWidth( unsigned int );
        void setHeight( unsigned int );
        void setNumBands( unsigned int );

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

#endif // VIGRA_IMPEX_TIFF_HXX
