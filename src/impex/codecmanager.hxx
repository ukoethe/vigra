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

#ifndef VIGRA_IMPEX_CODECMANAGER_HXX
#define VIGRA_IMPEX_CODECMANAGER_HXX

#include <map>
#include <memory>
#include "vigra/codec.hxx"

namespace vigra
{
    // the codec manager object is the global codec store
    // from which codecs can be allocated.

    class CodecManager
    {
        // the codec store
        std::vector<std::pair<std::vector<char>, std::string> > magicStrings;
        std::map< std::string, std::string > extensionMap;
        std::map< std::string, CodecFactory * > factoryMap;

    public:

        // singleton pattern
        static CodecManager & manager();

        // add an encoder to the stores
        void import( CodecFactory * cf );

        // find out if a given file type is supported
        bool fileTypeSupported( const std::string & fileType );

        std::vector<std::string>
        queryCodecPixelTypes( const std::string & codecName ) const;

        std::vector<int>
        queryCodecBandNumbers( const std::string & codecName ) const;

        // find out which file types are supported
        std::vector<std::string> supportedFileTypes();

        // find out which file extensions are supported
        std::vector<std::string> supportedFileExtensions();

        // look up decoder from the list, then return it
        std::auto_ptr<Decoder>
        getDecoder( const std::string & fileName,
                    const std::string & fileType = "undefined" ) const;

        // look up encoder from the list, then return it
        std::auto_ptr<Encoder>
        getEncoder( const std::string & fileName,
                    const std::string & fileType = "undefined" ) const;

        // try to figure out the correct file type
        std::string getFileTypeByMagicString( const std::string & filename ) const;

    private:

        // this will only be called by the singleton pattern
        CodecManager();

    }; // class CodecManager

    // singleton pattern
    inline CodecManager & codecManager()
    {
        return CodecManager::manager();
    }

} // namespace vigra

#endif // VIGRA_IMPEX_CODECMANAGER_HXX
