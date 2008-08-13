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

#ifndef VIGRA_IMPEX_AUTO_FILE_HXX
#define VIGRA_IMPEX_AUTO_FILE_HXX

#include <string>
#include <cstdio>
#include "vigra/error.hxx"

namespace vigra
{
    // resource aquisition is initialisation

    struct auto_file
    {
        FILE * m_file;

        auto_file( const char * name, const char * mode )
        : m_file(0)
        {
            m_file = VIGRA_CSTD::fopen( name, mode );
            if(!m_file)
            {
                std::string msg("Unable to open file '");
                msg += name;
                msg += "'.";
                vigra_precondition(0, msg.c_str());
            }
        }

        ~auto_file()
        {
            if(m_file)
                VIGRA_CSTD::fclose(m_file);
        }

        FILE * get()
        {
            return m_file;
        }

    private:

        // this can not be default constructed, copy constructed or assigned
        auto_file();
        auto_file( const auto_file & );
        auto_file & operator=( const auto_file & );
    };
}

#endif // VIGRA_IMPEX_AUTO_FILE_HXX
