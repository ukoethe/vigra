/************************************************************************/
/*                                                                      */
/*    Copyright 2002-2013 by Gunnar Kedenburg and Ullrich Koethe        */
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
/* Modifications by Pablo d'Angelo
 * updated to vigra 1.4 by Douglas Wilkins
 * as of 18 Febuary 2006:
 *  - Added UINT16 and UINT32 pixel types.
 *  - Added support for obtaining extra bands beyond RGB.
 *  - Added support for a position field that indicates the start of this
 *    image relative to some global origin.
 *  - Added support for x and y resolution fields.
 *  - Added support for ICC profiles
 */

#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "vigra/array_vector.hxx"
#include "vigra/imageinfo.hxx"
#include "codecmanager.hxx"

#if defined(_MSC_VER)
#  include "vigra/windows.h"
#  include <io.h>
#else
#  include <dirent.h>
#  include <unistd.h>
#endif

#ifdef HasHDF5
# include "vigra/hdf5impex.hxx"
#endif

namespace vigra
{

// build a string from a sequence.
#if defined(_MSC_VER) && (_MSC_VER < 1300)
template <class iterator>
std::string stringify (const iterator &start, const iterator &end)
{
    return stringifyImpl(start, end, *start);
}

template <class iterator, class Value>
std::string stringifyImpl (const iterator &start, const iterator &end, Value const &)
{
    std::ostringstream out;
    // do not place a space character after the last sequence element.
    std::copy (start, end - 1,
               std::ostream_iterator <Value> (out, " "));
    out << *(end-1);
    return out.str ();
}

#else

template <class iterator>
std::string stringify (const iterator &start, const iterator &end)
{
    typedef typename std::iterator_traits<iterator>::value_type value_type;
    std::ostringstream out;
    // do not place a space character after the last sequence element.
    std::copy (start, end - 1,
               std::ostream_iterator <value_type> (out, " "));
    out << *(end-1);
    return out.str ();
}

#endif // _MSC_VER < 1300

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

bool isImage(char const * filename)
{
#ifdef _MSC_VER
    return _access(filename, 0) != -1 && CodecManager::manager().getFileTypeByMagicString(filename) != "";
#else
    return access(filename, F_OK) == 0 && CodecManager::manager().getFileTypeByMagicString(filename) != "";
#endif
}

} // namespace vigra
