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

#include <algorithm>
#include "vigra/compression.hxx"
#include "lz4.h"

#ifdef HasZLIB
#include <zlib.h>
#endif

namespace vigra {

std::size_t compressImpl(char const * source, std::size_t srcSize, 
                         ArrayVector<char> & buffer,
                         CompressionMethod method)
{
    switch(method)
    {
      case NO_COMPRESSION:
      {
        ArrayVector<char>(source, source+srcSize).swap(buffer);
        return srcSize;
      }
      case ZLIB:
      case ZLIB_NONE:
      case ZLIB_FAST:
      case ZLIB_BEST:
      {
    #ifdef HasZLIB
        uLong destSize = ::compressBound(srcSize);
        buffer.resize(destSize);
        int res = ::compress2((Bytef *)buffer.data(), &destSize, (Bytef *)source, srcSize, method);
        vigra_postcondition(res == Z_OK, "compress(): zlib compression failed.");                    
        return destSize;
    #else
        vigra_precondition(false, "compress(): VIGRA was compiled without ZLIB compression.");
        return 0;
    #endif
      }
      case DEFAULT_COMPRESSION:
      case LZ4:
      {
        std::size_t destSize = ::LZ4_compressBound(srcSize);
        buffer.resize(destSize);
        destSize = ::LZ4_compress(source, buffer.data(), srcSize);
        vigra_postcondition(destSize > 0, "compress(): lz4 compression failed.");
        return destSize;
      }

#if 0  // currently unsupported
      case SNAPPY:
      {
    #ifdef HasSNAPPY
        std::size_t destSize = snappy::MaxCompressedLength(srcSize);
        buffer.resize(destSize);
        snappy::RawCompress(source, srcSize, buffer.data(), &destSize);
        return destSize;
    #else
        vigra_precondition(false, "compress(): VIGRA was compiled without SNAPPY compression.");
        return 0;
    #endif
      }
      case LZO:
      {
    #ifdef HasLZO
        static const int workmemSize = 
                  (LZO1X_1_MEM_COMPRESS + sizeof(lzo_align_t) - 1) / sizeof(lzo_align_t);
        ArrayVector<lzo_align_t> wrkmem(workmemSize);
        
        lzo_uint destSize = srcSize + srcSize / 16 + 64 + 3;
        buffer.resize(destSize);
        int res = ::lzo1x_1_compress((const lzo_bytep)source, srcSize, 
                                     (lzo_bytep)buffer.data(), &destSize, 
                                     wrkmem.data());
        return destSize;
    #else
        vigra_precondition(false, "compress(): VIGRA was compiled without LZO compression.");
        return 0;
    #endif
      }
#endif
      default:
        vigra_precondition(false, "compress(): Unknown compression method.");
    }
    return 0;
}

void compress(char const * source, std::size_t size, ArrayVector<char> & dest, CompressionMethod method)
{
    ArrayVector<char> buffer;
    std::size_t destSize = compressImpl(source, size, buffer, method);
    dest.insert(dest.begin(), buffer.data(), buffer.data() + destSize);
}

void compress(char const * source, std::size_t size, std::vector<char> & dest, CompressionMethod method)
{
    ArrayVector<char> buffer;
    std::size_t destSize = compressImpl(source, size, buffer, method);
    dest.insert(dest.begin(), buffer.data(), buffer.data() + destSize);
}

void uncompress(char const * source, std::size_t srcSize, 
                char * dest, std::size_t destSize, CompressionMethod method)
{
    switch(method)
    {
      case NO_COMPRESSION:
      {
        std::copy(source, source+srcSize, dest);
        break;
      }
      case ZLIB:
      case ZLIB_NONE:
      case ZLIB_FAST:
      case ZLIB_BEST:
      {
    #ifdef HasZLIB
        uLong destLen = destSize;
        int res = ::uncompress((Bytef *)dest, &destLen, (Bytef *)source, srcSize);
        vigra_postcondition(res == Z_OK, "uncompress(): zlib decompression failed.");
    #else
        vigra_precondition(false, "uncompress(): VIGRA was compiled without ZLIB compression.");
    #endif
        break;
      }
      case DEFAULT_COMPRESSION:
      case LZ4:
      {
        int sourceLen = ::LZ4_decompress_fast(source, dest, destSize);
        vigra_postcondition(sourceLen == srcSize, "uncompress(): lz4 decompression failed.");
        break;
      }
      
#if 0 // currently unsupported
      case SNAPPY:
      {
    #ifdef HasSNAPPY
        snappy::RawUncompress(source, srcSize, dest);
    #else
        vigra_precondition(false, "compress(): VIGRA was compiled without SNAPPY compression.");
    #endif
        break;
      }
      case LZO:
      {
    #ifdef HasLZO
        lzo_uint destLen = destSize;
        int res = ::lzo1x_decompress((const lzo_bytep)source, srcSize,
                                     (lzo_bytep)dest, &destLen, NULL);
        vigra_postcondition(res == LZO_E_OK, "uncompress(): lzo decompression failed.");
    #else
        vigra_precondition(false, "compress(): VIGRA was compiled without LZO compression.");
    #endif
        break;
      }
#endif
      default:
        vigra_precondition(false, "uncompress(): Unknown compression method.");
    }
}

/** Uncompress a data buffer when the uncompressed size is unknown.

    The destination array will be resized as required.
*/
void uncompress(char const * source, std::size_t size, ArrayVector<char> & dest, CompressionMethod method);
void uncompress(char const * source, std::size_t size, std::vector<char> & dest, CompressionMethod method);


} // namespace vigra
