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

#ifndef VIGRA_IMPEX_SUN_HXX
#define VIGRA_IMPEX_SUN_HXX

#include "vigra/codec.hxx"

// SUN Rasterfile format

namespace vigra {

    struct SunDecoderImpl;
    struct SunEncoderImpl;

    struct SunCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        std::auto_ptr<Decoder> getDecoder() const;
        std::auto_ptr<Encoder> getEncoder() const;
    };

    class SunDecoder : public Decoder
    {
        SunDecoderImpl * pimpl;

    public:

        SunDecoder() : pimpl(0) {}

        ~SunDecoder();
        void init( const std::string & );
        void close();
        void abort();

        std::string getFileType() const;
        std::string getPixelType() const;

        unsigned int getWidth() const;
        unsigned int getHeight() const;
        unsigned int getNumBands() const;
        unsigned int getOffset() const;

        const void * currentScanlineOfBand( unsigned int ) const;
        void nextScanline();
    };

    class SunEncoder : public Encoder
    {
        SunEncoderImpl * pimpl;

    public:

        SunEncoder() : pimpl(0) {}

        ~SunEncoder();
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
        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int );
        void nextScanline();
    };
}

#endif // VIGRA_IMPEX_SUN_HXX
