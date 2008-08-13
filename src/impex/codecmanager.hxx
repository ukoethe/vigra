/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Gunnar Kedenburg                */
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

        // look up encoder type from the list
        std::string
        getEncoderType( const std::string & fileName,
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
        
        ~CodecManager();

    }; // class CodecManager

    // singleton pattern
    inline CodecManager & codecManager()
    {
        return CodecManager::manager();
    }

} // namespace vigra

#endif // VIGRA_IMPEX_CODECMANAGER_HXX
