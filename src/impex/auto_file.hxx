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
            m_file = std::fopen( name, mode );
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
                std::fclose(m_file);
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
