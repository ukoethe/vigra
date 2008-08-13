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

#ifndef VIGRA_IMPEX_VIFF_HXX
#define VIGRA_IMPEX_VIFF_HXX

#include "vigra/codec.hxx"

// VIFF - Khoros Visualization/Image File Format

namespace vigra {

    struct ViffCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        std::auto_ptr<Decoder> getDecoder() const;
        std::auto_ptr<Encoder> getEncoder() const;
    };

    struct ViffDecoderImpl;
    struct ViffEncoderImpl;

    class ViffDecoder : public Decoder
    {
        ViffDecoderImpl * pimpl;

    public:

        ViffDecoder() : pimpl(0) {}

        ~ViffDecoder();
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

    class ViffEncoder : public Encoder
    {
        ViffEncoderImpl * pimpl;

    public:

        ViffEncoder() : pimpl(0) {}

        ~ViffEncoder();
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

#endif // VIGRA_IMPEX_VIFF_HXX
