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

#include <fstream>
#include <algorithm>
#include <cctype> // std::tolower
#ifdef DEBUG
#include <iostream>
#endif

#include "vigra/config.hxx"
#include "vigra/error.hxx"
#include "codecmanager.hxx"

// the codecs
#include "jpeg.hxx"
#include "tiff.hxx"
#include "viff.hxx"
#include "sun.hxx"
#include "png.hxx"
#include "pnm.hxx"
#include "bmp.hxx"
#include "gif.hxx"
#include "hdr.hxx"
#include "exr.hxx"

namespace vigra
{
    // singleton pattern
    CodecManager & CodecManager::manager()
    {
        static CodecManager manager;
        return manager;
    }

    CodecManager::CodecManager()
    {
#ifdef HasPNG
        import( new PngCodecFactory() );
#endif
#ifdef HasJPEG
        import( new JPEGCodecFactory() );
#endif
#ifdef HasTIFF
        import( new TIFFCodecFactory() );
#endif
#ifdef HasEXR
        import( new ExrCodecFactory() );
#endif
        import( new SunCodecFactory() );
        import( new PnmCodecFactory() );
        import( new ViffCodecFactory() );
        import( new BmpCodecFactory() );
        import( new GIFCodecFactory() );
        import( new HDRCodecFactory() );
    }

    CodecManager::~CodecManager() {
       // release previously allocated codecs
       // (use erase ideom similar to
       // S. Meyers' "Effective STL", Item 9)
       for (std::map< std::string, CodecFactory * >::iterator i
             = factoryMap.begin();
            i != factoryMap.end();
            /* nothing */ ) {
         delete (*i).second;
         factoryMap.erase(i++);
       }
    }

    // add an encoder to the stores
    void CodecManager::import( CodecFactory * cf )
    {
        CodecDesc desc = cf->getCodecDesc();

        // fill extension map
        const std::vector<std::string> & ext = desc.fileExtensions;
        typedef std::vector<std::string>::const_iterator iter_type;
        for( iter_type iter = ext.begin(); iter < ext.end(); ++iter )
            extensionMap[*iter] = desc.fileType;

        // fill magic vector
        for( VIGRA_CSTD::size_t i = 0; i < desc.magicStrings.size(); ++i )
            magicStrings.push_back( std::pair<std::vector<char>, std::string>
                                    ( desc.magicStrings[i],desc.fileType ) );
        // fill factory map
        factoryMap[desc.fileType] = cf;
    }

    // find out which pixel types a given codec supports
    std::vector<std::string>
    CodecManager::queryCodecPixelTypes( const std::string & filetype ) const
    {
        std::map< std::string, CodecFactory * >::const_iterator result
            = factoryMap.find( filetype );
        std::string message("queryCodecPixelTypes(): codec '");
        message += filetype + "' does not exist";
        vigra_precondition( result != factoryMap.end(), message.c_str());

        return result->second->getCodecDesc().pixelTypes;
    }

    // find out which pixel types a given codec supports
    std::vector<int>
    CodecManager::queryCodecBandNumbers( const std::string & filetype ) const
    {
        std::map< std::string, CodecFactory * >::const_iterator result
            = factoryMap.find( filetype );
        vigra_precondition( result != factoryMap.end(),
        "the codec that was queried for its pixeltype does not exist" );

        return result->second->getCodecDesc().bandNumbers;
    }

    // find out if a given file type is supported
    bool CodecManager::fileTypeSupported( const std::string & fileType )
    {
        std::map< std::string, CodecFactory * >::const_iterator search
            = factoryMap.find( fileType );
        return ( search != factoryMap.end() );
    }

    // find out which file types are supported
    std::vector<std::string> CodecManager::supportedFileTypes()
    {
        std::vector<std::string> fileTypes;
        std::map< std::string, CodecFactory * >::const_iterator iter
            = factoryMap.begin();
        while ( iter != factoryMap.end() ) {
            fileTypes.push_back( iter->first );
            ++iter;
        }
        std::sort(fileTypes.begin(), fileTypes.end());
        return fileTypes;
    }

    // find out which file extensions are supported
    std::vector<std::string> CodecManager::supportedFileExtensions()
    {
        std::vector<std::string> fileExtensions;
        std::map< std::string, std::string >::const_iterator iter
            = extensionMap.begin();
        while ( iter != extensionMap.end() ) {
            fileExtensions.push_back( iter->first );
            ++iter;
        }
        std::sort(fileExtensions.begin(), fileExtensions.end());
        return fileExtensions;
    }

    std::string
    CodecManager::getFileTypeByMagicString( const std::string & filename )
        const
    {
        // support for reading the magic string from stdin has been dropped
        // it was not guaranteed to work by the Standard

        // get the magic string
        const unsigned int magiclen = 4;
        char fmagic[magiclen];
#ifdef VIGRA_NEED_BIN_STREAMS
        std::ifstream stream(filename.c_str(), std::ios::binary);
#else
        std::ifstream stream(filename.c_str());
#endif
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
        stream.read( fmagic, magiclen );
        stream.close();

        // compare with the known magic strings
        typedef std::vector< std::pair< std::vector<char>, std::string > >
            magic_type;
        for( magic_type::const_iterator iter = magicStrings.begin();
             iter < magicStrings.end(); ++iter ) {
            const std::vector<char> & magic = iter->first;
            if ( std::equal( magic.begin(), magic.end(), fmagic ) )
                return iter->second;
        }

        // did not find a matching string
        return std::string();
    }

    // look up decoder from the list, then return it
    std::auto_ptr<Decoder>
    CodecManager::getDecoder( const std::string & filename,
                              const std::string & filetype,
                              unsigned int imageindex ) const
    {
        std::string fileType = filetype;

        if ( fileType == "undefined" ) {

            fileType = getFileTypeByMagicString(filename);
            vigra_precondition( !fileType.empty(),
                                "did not find a matching file type." );

#ifdef DEBUG
            std::cerr << "detected " << fileType
                      << " file format by magicstring of " << filename
                      << std::endl;
#endif
        }

        // return a codec factory by the file type
        std::map< std::string, CodecFactory * >::const_iterator search
            = factoryMap.find(fileType);
        vigra_precondition( search != factoryMap.end(),
        "did not find a matching codec for the given filetype" );

        // okay, we can return a decoder
        std::auto_ptr<Decoder> dec = search->second->getDecoder();
        dec->init(filename, imageindex);
        return dec;
    }

    // look up encoder from the list, then return it
    std::string
    CodecManager::getEncoderType( const std::string & filename,
                                  const std::string & fType,
                                  const std::string & mode ) const
    {
        std::string fileType = fType;

        if ( fileType == "" || fileType == "undefined" ) {
            // look up the file type by the file extension
            std::string ext
                = filename.substr( filename.find_last_of(".") + 1 );
            std::transform( ext.begin(), ext.end(), ext.begin(), (int (*)(int))&std::tolower );
            std::map< std::string, std::string >::const_iterator search
                = extensionMap.find(ext);
            vigra_precondition( search != extensionMap.end(),
            "did not find a matching codec for the given file extension" );
            // at this point, we have found a valid fileType
            fileType = search->second;
        }

        return fileType;
    }

    // look up encoder from the list, then return it
    std::auto_ptr<Encoder>
    CodecManager::getEncoder( const std::string & filename,
                              const std::string & fType,
                              const std::string & mode ) const
    {
        std::string fileType = getEncoderType(filename, fType);

        // return a codec factory by the file type
        std::map< std::string, CodecFactory * >::const_iterator search
            = factoryMap.find( fileType );
        vigra_precondition( search != factoryMap.end(),
        "did not find a matching codec for the given filetype" );

        // okay, we can return an encoder
        std::auto_ptr<Encoder> enc = search->second->getEncoder();
        enc->init(filename, mode);
        return enc;
    }

    // get a decoder
    std::auto_ptr<Decoder>
    getDecoder( const std::string & filename, const std::string & filetype, unsigned int imageindex )
    {
        return codecManager().getDecoder( filename, filetype, imageindex );
    }

    // get an encoder type
    std::string
    getEncoderType( const std::string & filename, const std::string & filetype )
    {
        return codecManager().getEncoderType( filename, filetype );
    }

    // get an encoder
    std::auto_ptr<Encoder>
    getEncoder( const std::string & filename, const std::string & filetype, const std::string & mode )
    {
        return codecManager().getEncoder( filename, filetype, mode );
    }

    std::vector<std::string>
    queryCodecPixelTypes( const std::string & codecname )
    {
        return codecManager().queryCodecPixelTypes(codecname);
    }


    // return true if downcasting is required, false otherwise
    bool negotiatePixelType( std::string const & codecname,
                 std::string const & srcPixeltype, std::string & destPixeltype )
    {
        std::vector<std::string> ptypes
            = codecManager().queryCodecPixelTypes(codecname);

        std::vector<std::string>::iterator pend;
        if(destPixeltype != "")
        {
            pend = std::find(ptypes.begin(), ptypes.end(), destPixeltype);
            if(pend == ptypes.end())
            {
                std::string msg("exportImage(): file type ");
                msg += codecname + " does not support requested pixel type "
                                          + destPixeltype + ".";
                vigra_precondition(false, msg.c_str());
            }
            ++pend;
        }
        else
        {
            pend = ptypes.end();
        }

        std::vector<std::string>::const_iterator result
                                   = std::find( ptypes.begin(), pend, srcPixeltype );

        if( result == pend)
        {
            if(destPixeltype == "")
                destPixeltype = "UINT8";
            // must always downcast
            return true;
        }
        else
        {
            if(destPixeltype == "")
                destPixeltype = srcPixeltype;
            // don't downcast
            return false;
        }
    }

    bool isPixelTypeSupported( const std::string & codecname,
                               const std::string & pixeltype )
    {
        std::vector<std::string> ptypes
            = codecManager().queryCodecPixelTypes(codecname);
        std::vector<std::string>::const_iterator result
            = std::find( ptypes.begin(), ptypes.end(), pixeltype );
        return ( result != ptypes.end() );
    }

    bool isBandNumberSupported( const std::string & codecname,
                                int bands )
    {
        std::vector<int> bandNumbers
            = codecManager().queryCodecBandNumbers(codecname);
        if(bandNumbers[0] == 0)
            return true; // any band number supported
        std::vector<int>::const_iterator result
            = std::find( bandNumbers.begin(), bandNumbers.end(), bands );
        return ( result != bandNumbers.end() );
    }

} // namespace vigra
