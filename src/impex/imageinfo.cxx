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

#include <algorithm>
#include <iterator>
#include <strstream>
#include <vector>
#include <iterator>
#include "vigra/imageinfo.hxx"
#include "codecmanager.hxx"

namespace vigra
{
    // build a string from a sequence.
    template< class iterator >
    std::string stringify( const iterator & start, const iterator & end )
    {
        typedef typename std::iterator_traits<iterator>::value_type value_type;
        std::ostrstream out;
        iterator last = end - 1;
        std::copy( start, last,
                   std::ostream_iterator<value_type>( out, " " ) );
        out << *last << (char)0;
        return out.str();
    }

    void validate_filetype( std::string filetype )
    {
        vigra_precondition( codecManager().fileTypeSupported(filetype),
                            "given file type is not supported" );
    }

    std::string impexListFormats()
    {
        std::vector<std::string> ft = codecManager().supportedFileTypes();
        return stringify( ft.begin(), ft.end() );
    }

    std::string impexListExtensions()
    {
        std::vector<std::string> ft = codecManager().supportedFileExtensions();
        return stringify( ft.begin(), ft.end() );
    }

    // class ImageExportInfo

    ImageExportInfo::ImageExportInfo( const char * filename )
        : m_filename(filename)
    {}

    ImageExportInfo & ImageExportInfo::setFileType( const char * filetype )
    {
        m_filetype = filetype;
        return *this;
    }

    ImageExportInfo & ImageExportInfo::setCompression( const char * comp )
    {
        m_comp = comp;
        return *this;
    }

    const char * ImageExportInfo::getFileName() const
    {
        return m_filename.c_str();
    }

    const char * ImageExportInfo::getFileType() const
    {
        return m_filetype.c_str();
    }

    ImageExportInfo & ImageExportInfo::setPixelType( const char * s )
    {
        m_pixeltype = s;
        return *this;
    }

    const char * ImageExportInfo::getPixelType() const
    {
        return m_pixeltype.c_str();
    }

    const char * ImageExportInfo::getCompression() const
    {
        return m_comp.c_str();
    }

    float ImageExportInfo::getXResolution() const
    {
        return m_x_res;
    }

    float ImageExportInfo::getYResolution() const
    {
        return m_y_res;
    }

    ImageExportInfo & ImageExportInfo::setXResolution( float val )
    {
        m_x_res = val;
        return *this;
    }

    ImageExportInfo & ImageExportInfo::setYResolution( float val )
    {
        m_y_res = val;
        return *this;
    }

    // return an encoder for a given ImageExportInfo object
    std::auto_ptr<Encoder> encoder( const ImageExportInfo & info )
    {
        std::auto_ptr<Encoder> enc;

        std::string filetype = info.getFileType();
        if ( filetype != "" ) {
            validate_filetype(filetype);
            std::auto_ptr<Encoder> enc2
                = getEncoder( std::string( info.getFileName() ), filetype );
            enc = enc2;
        } else {
            std::auto_ptr<Encoder> enc2
                = getEncoder( std::string( info.getFileName() ) );
            enc = enc2;
        }

        std::string comp = info.getCompression();
        if ( comp != "" ) {

            // check for JPEG compression
            int quality = -1;
            std::istrstream compstream(comp.c_str());
            compstream >> quality;
            if ( quality != -1 ) {
                enc->setCompressionType( "JPEG", quality );
                return enc;
            }

            // leave any other compression type to the codec
            enc->setCompressionType(comp);
        }

        std::string pixel_type = info.getPixelType();
        if ( pixel_type != "" ) {
            enc->setPixelType(pixel_type);
        }

        return enc;
    }

    // class ImageImportInfo

    ImageImportInfo::ImageImportInfo( const char * filename )
        : m_filename(filename)
    {
        std::auto_ptr<Decoder> decoder = getDecoder(m_filename);

        m_filetype = decoder->getFileType();
        m_pixeltype = decoder->getPixelType();
        m_width = decoder->getWidth();
        m_height = decoder->getHeight();
        m_num_bands = decoder->getNumBands();

        decoder->abort(); // there probably is no better way than this
    }

    const char * ImageImportInfo::getFileName() const
    {
        return m_filename.c_str();
    }

    const char * ImageImportInfo::getFileType() const
    {
        return m_filetype.c_str();
    }

    const char * ImageImportInfo::getPixelType() const
    {
        return m_pixeltype.c_str();
    }

    ImageImportInfo::PixelType ImageImportInfo::pixelType() const
    {
       const std::string pixeltype=getPixelType();
       if (pixeltype == "UINT8")
         return UINT8;
       if (pixeltype == "INT16")
         return INT16;
       if (pixeltype == "INT32")
         return INT32;
       if (pixeltype == "FLOAT")
         return FLOAT;
       if (pixeltype == "DOUBLE")
         return DOUBLE;
       vigra_fail( "internal error: unknown pixel type" );
       return ImageImportInfo::PixelType();
    }

    int ImageImportInfo::width() const
    {
        return m_width;
    }

    int ImageImportInfo::height() const
    {
        return m_height;
    }

    int ImageImportInfo::numBands() const
    {
        return m_num_bands;
    }

    Diff2D ImageImportInfo::size() const
    {
        return Diff2D( m_width, m_height );
    }

    bool ImageImportInfo::isGrayscale() const
    {
        return m_num_bands == 1;
    }

    bool ImageImportInfo::isColor() const
    {
        return m_num_bands == 3;
    }

    bool ImageImportInfo::isByte() const
    {
        return m_pixeltype == "UINT8";
    }

    float ImageImportInfo::getXResolution() const
    {
        return m_x_res;
    }

    float ImageImportInfo::getYResolution() const
    {
        return m_y_res;
    }

    // return a decoder for a given ImageImportInfo object
    std::auto_ptr<Decoder> decoder( const ImageImportInfo & info )
    {
        std::string filetype = info.getFileType();
        validate_filetype(filetype);
        return getDecoder( std::string( info.getFileName() ), filetype );
    }

} // namespace vigra
