/************************************************************************/
/*                                                                      */
/*               Copyright 2013-2014 by Ullrich Koethe                  */
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


#ifndef VIGRA_COMPRESSION_HXX
#define VIGRA_COMPRESSION_HXX

#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"
#include <vector>

namespace vigra {

enum CompressionMethod {  ZLIB=-1,     // default compression using zlib
                          ZLIB_NONE=0, // no compression using zlib
                          ZLIB_FAST=1, // fastest compression using zlib
                          ZLIB_BEST=9, // highest compression using zlib
                          LZ4          // very fast LZ4 algorithm
                       };

/** Compress the source buffer.

    The destination array will be resized as required.
*/
VIGRA_EXPORT void compress(char const * source, std::size_t size, ArrayVector<char> & dest, CompressionMethod method);
VIGRA_EXPORT void compress(char const * source, std::size_t size, std::vector<char> & dest, CompressionMethod method);

/** Uncompress the source buffer when the uncompressed size is known.

    The destination buffer must be allocated to the correct size.
*/
VIGRA_EXPORT void uncompress(char const * source, std::size_t srcSize, 
                             char * dest, std::size_t destSize, CompressionMethod method);


} // namespace vigra

#endif // VIGRA_COMPRESSION_HXX
