/************************************************************************/
/*                                                                      */
/*     Copyright 2012-2014 by Ullrich Koethe and Thorben Kroeger        */
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

/* benchmark results for a simple loop 'if(iter.get<1>() != count++)'

    ********************
    image size: 200^3, chunk size: 64^3, i.e. chunk count: 4^3
    times in msec, excluding time to store file on disk
    
    win64/vs2012 (koethe-laptop):            uint8     float    double
    plain array                                 18        18        18
    chunked array (all in cache)                25        26        26
    thread-safe chunked (all in cache)          27        28        29
    thread-safe chunked (1 slice in cache)      29        33        39
    thread-safe chunked (1 row in cache)        45        48        52
    chunked (initial creation, all in cache)    33        43        57
    
    linux/gcc 4.7.3 (birdofprey):            uint8     float    double
    plain array                                 16        20        21
    chunked array (all in cache)                17        23        24
    thread-safe chunked (all in cache)          19        24        25
    thread-safe chunked (1 slice in cache)      20        29        34
    thread-safe chunked (1 row in cache)        24        33        39
    chunked (initial creation, all in cache)    22        34        48
    
    OS X 10.7:                               uint8     float    double
    plain array                                 11        22        24
    chunked array (all in cache)                --        --        --
    thread-safe chunked (all in cache)          20        25        26
    thread-safe chunked (1 slice in cache)      23        37        46
    thread-safe chunked (1 row in cache)        32        50        56
    chunked (initial creation, all in cache)    34        56        77

    **********************
    image size: 400^3, chunk size: 127^3, i.e. chunk count: 4^3
    times in msec, excluding time to store file on disk
    
    win64/vs2012 (koethe-laptop):            uint8     float    double
    plain array                                130       130       130
    chunked array (all in cache)               190       190       200
    thread-safe chunked (all in cache)         190       200       210
    thread-safe chunked (1 slice in cache)     210       235       280
    thread-safe chunked (1 row in cache)       240       270       300
    chunked (initial creation, all in cache)   230       300       400
    
    linux/gcc 4.7.3 (birdofprey):            uint8     float    double
    plain array                                130       162       165
    chunked array (all in cache)               131       180       184
    thread-safe chunked (all in cache)         135       183       188
    thread-safe chunked (1 slice in cache)     146       218       258
    thread-safe chunked (1 row in cache)       154       229       270
    chunked (initial creation, all in cache)   173       269       372

    ***********************
    Compression:
    * I tried ZLIB, LZO, SNAPPY, LZ4, LZFX and FASTLZ (with compression levels 1 -- faster
      and level 2 -- higher compression). There are also QuickLZ and LZMAT which claim 
      to be fast, but they are under a GPL license.
    * ZLIB compresses best, but is quite slow even at compression level 1 
      (times are in ms and include compression and decompression).
                 byte   float   double
        ZLIB      121    3100     5800
        ZLIB1      79    1110     1190
        LZO        43     104      280
        SNAPPY     46      71      305
        LZ4        42      70      283
        LZFX       76     278      330
        FASTLZ1    52     280      309
        FASTLZ1    53     286      339
    * The fast compression algorithms are unable to compress the float array
      and achieve ~50% for the double array, whereas ZLIB achieves 32% and 16%
      respectively (at the fastest compression level 1, it is still 33% and 17%
      respectively). LZFX cannot even compress the byte data (probably a bug?).
      Average compression ratios for the byte array are
        ZLIB:    2.3%
        ZLIB1:   4.6%
        LZO:     5.6%
        SNAPPY:  9.8%
        LZ4:     9.7%
        FASTLZ1: 7.6%
        FASTLZ2: 7.9%
    * LZO is under GPL (but there is a Java implementation under Apache license at
      http://svn.apache.org/repos/asf/hadoop/common/tags/release-0.19.2/src/core/org/apache/hadoop/io/compress/lzo/)
      The others are BSD and MIT (FASTLZ).
    * Snappy doesn't support Windows natively, but porting is simple (see my github repo)
    * The source code for LZO, LZ4, LZFX, and FASTLZ can simply be copied to VIGRA, 
      but LZO's GPL license is unsuitable.
    * HDF5 compression is already sufficient at level 1 (4-15%,
      higher levels don't lead to big gains) and only a factor 3-10 slower 
      than without compression.
*/

#ifndef VIGRA_MULTI_ARRAY_CHUNKED_HXX
#define VIGRA_MULTI_ARRAY_CHUNKED_HXX

#include "metaprogramming.hxx"
#include "multi_array.hxx"
#include "threading.hxx"
#include "hdf5impex.hxx"
#include "compression.hxx"

#ifdef _WIN32
# include "windows.h"
#else
# include <fcntl.h>
# include <stdlib.h>
# include <unistd.h>
# include <sys/stat.h>
# include <sys/mman.h>
#endif

double timeit = 0;

namespace vigra {

#ifdef __APPLE__
    #define VIGRA_NO_SPARSE_FILE
#endif

#ifdef _WIN32

inline 
void winErrorToException(std::string message = "") 
{ 
    LPVOID lpMsgBuf;
    DWORD dw = GetLastError(); 

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL );

    message += (char*)lpMsgBuf;
    LocalFree(lpMsgBuf);
    
    throw std::runtime_error(message);
}

inline 
std::string winTempFileName(std::string path = "") 
{ 
    if(path == "")
    {
        TCHAR default_path[MAX_PATH];
        if(!GetTempPath(MAX_PATH, default_path))
            winErrorToException("winTempFileName(): ");
        path = default_path;
    }
    
    TCHAR name[MAX_PATH];  
    if(!GetTempFileName(path.c_str(), TEXT("vigra"), 0, name))
        winErrorToException("winTempFileName(): ");

    return std::string(name);
}

inline 
std::size_t winClusterSize()
{
    SYSTEM_INFO info;
    ::GetSystemInfo(&info); 
    return info.dwAllocationGranularity;
}

#endif

namespace {

#ifdef _WIN32
std::size_t mmap_alignment = winClusterSize();
#else
std::size_t mmap_alignment = sysconf(_SC_PAGE_SIZE);
#endif

} // anonymous namespace

template <class T>
struct ChunkedMemory;

template <unsigned int N, class T>
class ChunkedArray;

template <unsigned int N, class T>
class ChunkedHandle;


template <unsigned int N, class T>
struct ChunkShape;

template <class T>
struct ChunkShape<1, T>
{
    static const unsigned int bits0 = 18;
    static const unsigned int mask0 = (1 << bits0) - 1;
    
    static void chunkIndex(Shape1 const & p, Shape1 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits0;
    }
    
    static std::size_t chunkOffset(Shape1 const & p, Shape1 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits0) * s[0];
    }
    
    static std::size_t offsetInChunk(Shape1 const & p, Shape1 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & mask0) * s[0];
    }

    static Shape1 chunkArrayShape(Shape1 shape)
    {
        shape[0] = (shape[0] + mask0) >> bits0;
        return shape;
    }

    static Shape1 standardChunkShape()
    {
        return Shape1(1 << bits0);
    }
};

template <class T>
struct ChunkShape<2, T>
{
    static const unsigned int bits0 = 9;
    static const unsigned int bits1 = 9;
    static const unsigned int mask0 = (1 << bits0) - 1;
    static const unsigned int mask1 = (1 << bits1) - 1;
    
    static void chunkIndex(Shape2 const & p, Shape2 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits0;
        b[1] = (UI)p[1] >> bits1;
    }
    
    static std::size_t chunkOffset(Shape2 const & p, Shape2 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits0) * s[0] +
               ((UI)p[1] >> bits1) * s[1];
    }
    
    static std::size_t offsetInChunk(Shape2 const & p, Shape2 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & mask0) * s[0] +
               ((UI)p[1] & mask1) * s[1];
    }

    static Shape2 chunkArrayShape(Shape2 shape)
    {
        shape[0] = (shape[0] + mask0) >> bits0;
        shape[1] = (shape[1] + mask1) >> bits1;
        return shape;
    }

    static Shape2 standardChunkShape()
    {
        return Shape2(1 << bits0, 1 << bits1);
    }
};

template <class T>
struct ChunkShape<3, T>
{
    static const unsigned int bits0 = 6;
    static const unsigned int bits1 = 6;
    static const unsigned int bits2 = 6;
    static const unsigned int mask0 = (1 << bits0) - 1;
    static const unsigned int mask1 = (1 << bits1) - 1;
    static const unsigned int mask2 = (1 << bits2) - 1;
    
    static void chunkIndex(Shape3 const & p, Shape3 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits0;
        b[1] = (UI)p[1] >> bits1;
        b[2] = (UI)p[2] >> bits2;
    }
    
    static std::size_t chunkOffset(Shape3 const & p, Shape3 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits0) * s[0] +
               ((UI)p[1] >> bits1) * s[1] +
               ((UI)p[2] >> bits2) * s[2];
    }
    
    static std::size_t offsetInChunk(Shape3 const & p, Shape3 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & mask0) * s[0] +
               ((UI)p[1] & mask1) * s[1] +
               ((UI)p[2] & mask2) * s[2];
    }

    static Shape3 chunkArrayShape(Shape3 shape)
    {
        shape[0] = (shape[0] + mask0) >> bits0;
        shape[1] = (shape[1] + mask1) >> bits1;
        shape[2] = (shape[2] + mask2) >> bits2;
        return shape;
    }

    static Shape3 standardChunkShape()
    {
        return Shape3(1 << bits0, 1 << bits1, 1 << bits2);
    }
};

template <class T>
struct ChunkShape<4, T>
{
    static const unsigned int bits0 = 6;
    static const unsigned int bits1 = 6;
    static const unsigned int bits2 = 4;
    static const unsigned int bits3 = 2;
    static const unsigned int mask0 = (1 << bits0) - 1;
    static const unsigned int mask1 = (1 << bits1) - 1;
    static const unsigned int mask2 = (1 << bits2) - 1;
    static const unsigned int mask3 = (1 << bits3) - 1;
    
    static void chunkIndex(Shape4 const & p, Shape4 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits0;
        b[1] = (UI)p[1] >> bits1;
        b[2] = (UI)p[2] >> bits2;
        b[3] = (UI)p[3] >> bits3;
    }
    
    static std::size_t chunkOffset(Shape4 const & p, Shape4 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits0) * s[0] +
               ((UI)p[1] >> bits1) * s[1] +
               ((UI)p[2] >> bits2) * s[2] +
               ((UI)p[3] >> bits3) * s[3];
    }
    
    static std::size_t offsetInChunk(Shape4 const & p, Shape4 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & mask0) * s[0] +
               ((UI)p[1] & mask1) * s[1] +
               ((UI)p[2] & mask2) * s[2] +
               ((UI)p[3] & mask3) * s[3];
    }

    static Shape4 chunkArrayShape(Shape4 shape)
    {
        shape[0] = (shape[0] + mask0) >> bits0;
        shape[1] = (shape[1] + mask1) >> bits1;
        shape[2] = (shape[2] + mask2) >> bits2;
        shape[3] = (shape[3] + mask3) >> bits3;
        return shape;
    }

    static Shape4 standardChunkShape()
    {
        return Shape4(1 << bits0, 1 << bits1, 1 << bits2, 1 << bits3);
    }
};

template <unsigned int N, class T>
class ChunkBase
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    
    ChunkBase()
    : pointer_(0),
      shape_(),
      strides_(),
      cache_next_(0),
      refcount_(0)
    {}
    
    ChunkBase(shape_type const & shape)
    : pointer_(0),
      shape_(shape),
      strides_(detail::defaultStride(shape)),
      cache_next_(0),
      refcount_(0)
    {}
    
    ChunkBase(ChunkBase const & rhs)
    : pointer_(0),
      shape_(rhs.shape_),
      strides_(rhs.strides_),
      cache_next_(0),
      refcount_(0)
    {}
    
    ChunkBase & operator=(ChunkBase const & rhs)
    {
        if(this != &rhs)
        {
            pointer p = pointer_.load(threading::memory_order_release);
            pointer rp = rhs.pointer_.load(threading::memory_order_release);
            if(p == 0 && rp == 0)
            {
                shape_ = rhs.shape_;
                strides_ = rhs.strides_;
            }
            else if(p != 0 && rp != 0)
            {
                vigra_precondition(shape_ == rhs.shape_, // implies stride equality
                    "ChunkBase::operator=(): shape mismatch.");
                std::copy(rp, rp + rhs.size(), p);
            }
            else
            {
                vigra_precondition(0,
                    "ChunkBase::operator=(): invalid assignment.");
            }
        }
        return *this;
    }
    
    MultiArrayIndex size() const
    {
        return prod(shape_);
    }
    
    threading::atomic<pointer> pointer_;
    shape_type shape_, strides_;
    ChunkBase * cache_next_;
    mutable threading::atomic<int> refcount_;
};

template <unsigned int N, class T>
class ChunkedSubarrayCopy
: public MultiArray<N, T>
{
  public:
    typedef ChunkedArray<N, T>  CArray;
    typedef MultiArray<N, T>        SArray;
    typedef typename  MultiArrayView<N, T>::difference_type shape_type;
  
    ChunkedSubarrayCopy(shape_type const & s)
    : MultiArray<N, T>(s),
      array_(0),
      offset_(),
      write_on_destruction_(false)
    {}
    
    ~ChunkedSubarrayCopy()
    {
        writeBack();
    }
    
    void writeBack()
    {
        if(!write_on_destruction_ || !array_)
            return;
        typename CArray::iterator i   = array_->begin().restrictToSubarray(offset_, offset_ + this->shape()),
                                  end = i.getEndIterator();
        typename SArray::iterator j   = this->begin();
        
        for(; i != end; ++i, ++j)
            *i = *j;
    }
  
    CArray * array_;
    shape_type offset_;
    bool write_on_destruction_;
};

template <unsigned int N, class T>
class ChunkedArrayBase
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    typedef ChunkedSubarrayCopy<N, T>  SubarrayType;
    typedef ChunkBase<N, T> Chunk;

    virtual ~ChunkedArrayBase()
    {}
    
    virtual pointer loadChunk(Chunk * chunk) = 0;
    
    virtual void unloadChunk(Chunk * chunk) = 0;
    
    virtual Chunk * lookupChunk(shape_type const & index) = 0;
    
    virtual void unrefChunk(ChunkedHandle<N, T> * h) = 0;
    
    virtual void unrefChunks(shape_type const & start, shape_type const & stop) = 0;
    
    virtual pointer getChunk(shape_type const & point, 
                             shape_type & strides, shape_type & upper_bound, 
                             ChunkedHandle<N, T> * h) = 0;

    virtual void copySubarray(shape_type const & start, 
                              ChunkedSubarrayCopy<N, T> & subarray, bool write_on_destruction = false) = 0;
};

template <unsigned int N, class T>
class MultiArrayView<N, T, ChunkedArrayTag>
{
  public:
    typedef typename MultiArrayShape<N>::type shape_type;
    typedef T* pointer;
    typedef T& reference;
    
    class Chunk
    {
      public:
        
        pointer pointer_;
        shape_type strides_;
    };
    
    struct UnrefProxy
    {
        UnrefProxy(shape_type const & start, shape_type const & stop,
                   ChunkedArrayBase<N, T> * array)
        : start_(start),
          stop_(stop),
          array_(array)
        {}
        
        ~UnrefProxy()
        {
            if(array_)
                array_->unrefChunks(start_, stop_);
        }
        
        shape_type start_, stop_;
        ChunkedArrayBase<N, T> * array_;
    };
    
    MultiArrayView(shape_type const & shape)
    : shape_(shape)
    {}
    
    reference operator[](shape_type p)
    {
        p += offset_;
        Chunk * chunk = chunks_.data() + ChunkShape<N, T>::chunkOffset(p, chunks_.stride());
        return *(chunk->pointer_ + ChunkShape<N, T>::offsetInChunk(p, chunk->strides_));
    }

    MultiArray<N, Chunk> chunks_;
    shape_type shape_, offset_;
    std::shared_ptr<UnrefProxy> unref_;
};

/*
The present implementation uses a memory-mapped sparse file to store the chunks.
A sparse file is created on Linux using the O_TRUNC flag (this seems to be 
the default file behavior on Linux anyway), and on Windows by
calling DeviceIoControl(file_handle, FSCTL_SET_SPARSE,...) after file creation.

We can automatically delete the file upon closing. On Windows, this happens
if the file was opened with FILE_FLAG_DELETE_ON_CLOSE. (It may be useful to
combine this with the flag FILE_ATTRIBUTE_TEMPORARY, which tells the OS to
avoid writing the file to disk if possible. However, judging from the timings,
something is still written, or cleanup takes considerable time.)
On Linux, you can call fileno(tmpfile()) for the same purpose.

Alternatives are:
* Keep the array in memory, but compress unused chunks.
* Don't create a file explicitly, but use the swap file instead. This is 
  achieved on Linux by mmap(..., MAP_PRIVATE | MAP_ANONYMOUS, -1, ...), 
  on Windows by calling CreateFileMapping(INVALID_HANDLE_VALUE, ...).
   * On Linux, the memory must not be unmapped because this
     looses the data. In fact, anonymous mmap() is very similar to 
     malloc(), and there is probably no good reason to use anonymous mmap().
   * On Windows, this is much faster, because the OS will never try to 
     actually write something to disk (unless swapping is necessary).
* Back chunks by HDF5 chunks, possibly using on-the-fly compression. This
  is in particular useful for existing HDF5 files.
* Back chunks by HDF5 datasets. This can be combined with compression 
  (both explicit and on-the-fly) or with memory mapping (using the 
  function H5Dget_offset() to get the offset from the beginning of the file).
    
FIXME:
* backends:
   * allocators are not used
   * HDF5 only works for scalar types so far
   * HDF5 must support read-only and read/write mode
   * temp file arrays in swap (just an API addition to the constructor)
   * support TIFF chunked reading
* the array implementations should go into cxx files in src/impex
  * this requires implementation of the low-level functions independently of dtype
    (use 'char *' and multiply shape and stride with sizeof(T))
  * don't forget to increment the soversion after the change
  * alternative: provide 'config_local.hxx' with flags for available packages
* decide on chunk locking policies for array views (in particular, for index access)
  * array view has functions fetch()/release() (better names?) to lock/unlock 
    _all_ chunks in the view
  * release() is automatically called in the destructor
  * it should be possible to call fetch in the constructor via a flag,
    but should the constructor fetch by default?
  * how should fetch() handle the case when the cache is too small
    * throw an exception?
    * silently enlarge the cache?
    * temporarily enlarge the cache?
    * provide an option to control the behavior?
  * also provide copySubarray() with ReadOnly and ReadWrite flags, where
    ReadWrite copies the subarray back in the destructor or on demand
    * locking is only required while each slice is copied
    * the copy functions can use normal array views and iterators
    * the ReadWrite version can store a checksum for each chunk (or part
      of a chunk) to detect collisions on write
    * use shared pointers to support memory management of the subarrays?
* find efficient ways to support slicing and transposition in the indexing
  functions of a view. 
  1. possibility: each view contains
      * an index object 'bound_index_' with original dimension whose values denote 
        coordinates of bound axes and offsets for unbound coordinates
      * a permutation object 'permutation_' with dimension of the view that maps 
        view coordinates to original coordinates
      * that is:
        operator[](index)
        {
            shape_type full_index(bound_index_);
            for(int k=0; k<N_view; ++k)
                full_index[permutation_[k]] += index[k];
            split full_index into chunk part and local part
            look up chunk
            return pixel
        }
      * maybe this is faster if it is combined with the stride computation?
      * an optimization for unsliced arrays is desirable
  2. possibility:
      * add the point offset to the low-dimensional index
      * split low-dimensional index into chunk part and local part
      * look up chunk
      * determine scalar pointer offset from local part and strides plus a 
        chunk-specific correction that can be stored in a 3^N array 
        - but can we efficiently determine where to look in that offset array?
  3. possibility:
      * don't care about speed - require copySubarray() if indexing should
        be fast
*/
template <unsigned int N, class T>
class ChunkedArray
: public ChunkedArrayBase<N, T>
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    // typedef typename CoupledHandleType<N, ChunkedMemory<T> >::type  P1;
    // typedef CoupledScanOrderIterator<N, P1>                         iterator;
    // typedef typename CoupledIteratorType<N, ChunkedMemory<T>>::type iterator;
    // typedef typename iterator::value_type                           handle;
    typedef StridedScanOrderIterator<N, ChunkedMemory<T>, T&, T*>   iterator;
    typedef ChunkedSubarrayCopy<N, T>                                  SubarrayType;
    typedef ChunkBase<N, T> Chunk;
    typedef MultiArrayView<N, T, ChunkedArrayTag>                   ViewType;
        
    ChunkedArray(int cache_max = 0)
    : shape_(0),
      chunk_shape_(ChunkShape<N, T>::standardChunkShape()),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      cache_max_size_(cache_max)
    {}
        
    ChunkedArray(shape_type const & shape, int cache_max = 0)
    : shape_(shape),
      chunk_shape_(ChunkShape<N, T>::standardChunkShape()),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      cache_max_size_(cache_max)
    {}
    
    virtual ~ChunkedArray()
    {}
    
    virtual pointer loadChunk(Chunk * chunk) = 0;
    
    virtual void unloadChunk(Chunk * chunk)
    {}
    
    virtual Chunk * lookupChunk(shape_type const & index) = 0;
    
    virtual void unrefChunk(ChunkedHandle<N, T> * h)
    {
        if(h->chunk_)
        {
            static_cast<Chunk*>(h->chunk_)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
            h->chunk_ = 0;
        }
    }
    
    virtual pointer getChunk(shape_type const & point, 
                             shape_type & strides, shape_type & upper_bound, 
                             ChunkedHandle<N, T> * h)
    {
        if(h->chunk_)
        {
            static_cast<Chunk*>(h->chunk_)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
            h->chunk_ = 0;
        }
        
        shape_type global_point = point + h->offset_;

        if(!this->isInside(global_point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N, T>::chunkIndex(global_point, chunkIndex);
        ChunkBase<N, T> * chunk = lookupChunk(chunkIndex);
        
        // Obtain a reference to the current chunk.
        // We use a simple spin-lock here because it is very fast in case of success
        // and failures (i.e. a collisions with another thread) are presumably 
        // very rare.
        while(true)
        {
            int rc = chunk->refcount_.load(threading::memory_order_acquire);
            if(rc < 0)
            {
                // cache management in progress => try again later
                threading::this_thread::yield();
            }
            else if(chunk->refcount_.compare_exchange_weak(rc, rc+1, threading::memory_order_seq_cst))
            {
                // success
                break;
            }
        }
        
        T * p = chunk->pointer_.load(threading::memory_order_acquire);
        if(p == 0)
        {
            // apply the double-checked locking pattern
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            p = chunk->pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                p = loadChunk(chunk);
                
                if(cache_max_size_ > 0)
                {
                    // insert in queue of mapped chunks
                    if(!cache_first_)
                        cache_first_ = chunk;
                    if(cache_last_)
                        cache_last_->cache_next_ = chunk;
                    cache_last_ = chunk;
                    ++cache_size_;
                    
                    // do cache management if cache is full
                    // (note that we still hold the cache_lock_)
                    if(cache_size_ > cache_max_size_)
                    {
                        // remove first element from queue
                        Chunk * c = cache_first_;
                        cache_first_ = c->cache_next_;
                        if(cache_first_ == 0)
                            cache_last_ = 0;
                        
                        // check if the refcount is zero and temporarily set it 
                        // to -1 in order to lock the chunk while cache management 
                        // is done
                        int rc = 0;
                        if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_acquire))
                        {
                            // refcount was zero => can unload
                            unloadChunk(c);
                            c->cache_next_ = 0;
                            --cache_size_;
                            c->refcount_.store(0, std::memory_order_release);
                        }
                        else
                        {
                            // refcount was non-zero => reinsert chunk into queue
                            if(!cache_first_)
                                cache_first_ = c;
                            if(cache_last_)
                                cache_last_->cache_next_ = c;
                            cache_last_ = c;
                        }
                    }
                }
            }
        }

        strides = chunk->strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_ - h->offset_;
        std::size_t offset = ChunkShape<N, T>::offsetInChunk(global_point, strides);
        h->chunk_ = chunk;
        return p + offset;
    }
    
    virtual void refChunks(shape_type const & start, shape_type const & stop)
    {
        vigra_precondition((shape_type() <= start).all() && (start < stop).all() && (stop < shape()).all(),
                           "ChunkedArray::copySubarray(): subarray out of bounds.");
        
        shape_type chunkStart(SkipInitialization), chunkStop(SkipInitialization);
        ChunkShape<N, T>::chunkIndex(start, chunkStart);
        ChunkShape<N, T>::chunkIndex(stop-shape_type(1), chunkStop);
        chunkStop += shape_type(1);
        
        threading::lock_guard<threading::mutex> guard(cache_lock_);
        
        MultiCoordinateIterator<N> i(chunkStart, chunkStop),
                                   end(i.getEndIterator());
        for(; i != end; ++i)
        {
            Chunk * chunk = lookupChunk(*i);
            chunk->refcount_.fetch_add(1, threading::memory_order_relaxed);
            if(chunk->pointer_.load(threading::memory_order_relaxed) == 0)
            {
                loadChunk(chunk);
                if(cache_max_size_ > 0)
                {
                    // insert in queue of mapped chunks
                    if(!cache_first_)
                        cache_first_ = chunk;
                    if(cache_last_)
                        cache_last_->cache_next_ = chunk;
                    cache_last_ = chunk;
                    ++cache_size_;
                }
            }
        }
    }
    
    virtual void unrefChunks(shape_type const & start, shape_type const & stop)
    {
        vigra_precondition((shape_type() <= start).all() && (start < stop).all() && (stop < shape()).all(),
                           "ChunkedArray::copySubarray(): subarray out of bounds.");
        
        shape_type chunkStart(SkipInitialization), chunkStop(SkipInitialization);
        ChunkShape<N, T>::chunkIndex(start, chunkStart);
        ChunkShape<N, T>::chunkIndex(stop-shape_type(1), chunkStop);
        chunkStop += shape_type(1);
        
        MultiCoordinateIterator<N> i(chunkStart, chunkStop),
                                   end(i.getEndIterator());
        for(; i != end; ++i)
        {
            Chunk * chunk = lookupChunk(*i);
            lookupChunk(*i)->refcount_.fetch_sub(1, threading::memory_order_acquire);
        }
        
        if(cache_max_size_ > 0)
        {
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            Chunk * c = cache_first_;
            while(cache_size_ > cache_max_size_ && c != 0)
            {
                Chunk * n = c->cache_next_;
                int rc = 0;
                if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_acquire))
                {
                    // refcount was zero => can unload
                    unloadChunk(c);
                    c->cache_next_ = 0;
                    --cache_size_;
                    c->refcount_.store(0, std::memory_order_release);
                }
                c = n;
            }
        }
    }
    
    void subarray(shape_type const & start, ViewType & view)
    {
        shape_type stop = start + view.shape_;
        
        vigra_precondition((shape_type() <= start).all() && (start < stop).all() && (stop < shape()).all(),
                           "ChunkedArray::copySubarray(): subarray out of bounds.");
        
        shape_type chunkStart(SkipInitialization), chunkStop(SkipInitialization);
        ChunkShape<N, T>::chunkIndex(start, chunkStart);
        ChunkShape<N, T>::chunkIndex(stop-shape_type(1), chunkStop);
        chunkStop += shape_type(1);
        
        view.chunks_.reshape(chunkStop-chunkStart);
        view.offset_ = start - chunkStart * this->chunk_shape_;
        
        typedef typename ViewType::UnrefProxy UP;
        view.unref_ = std::shared_ptr<UP>(new UP(start, stop, this));
        
        threading::lock_guard<threading::mutex> guard(cache_lock_);
        
        MultiCoordinateIterator<N> i(chunkStart, chunkStop),
                                   end(i.getEndIterator());
        for(; i != end; ++i)
        {
            Chunk * chunk = lookupChunk(*i);
            chunk->refcount_.fetch_add(1, threading::memory_order_relaxed);
            T * p = chunk->pointer_.load(threading::memory_order_relaxed);
            if(p == 0)
            {
                p = loadChunk(chunk);
                if(cache_max_size_ > 0)
                {
                    // insert in queue of mapped chunks
                    if(!cache_first_)
                        cache_first_ = chunk;
                    if(cache_last_)
                        cache_last_->cache_next_ = chunk;
                    cache_last_ = chunk;
                    ++cache_size_;
                }
            }
            typename ViewType::Chunk * vc = &view.chunks_[*i - chunkStart];
            vc->pointer_ = p;
            vc->strides_ = chunk->strides_;
        }
    }
    
    virtual void copySubarray(shape_type const & start, 
                              ChunkedSubarrayCopy<N, T> & subarray, bool write_on_destruction = false)
    {
        shape_type stop   = start + subarray.shape();
        
        vigra_precondition((shape_type() <= start).all() && (start < stop).all() && (stop < shape()).all(),
                           "ChunkedArray::copySubarray(): subarray out of bounds.");
                           
        iterator i(begin().restrictToSubarray(start, stop)),
                 end(i.getEndIterator());
        typename SubarrayType::iterator j = subarray.begin();
        
        for(; i != end; ++i, ++j)
        {
           *j = *i;
        }
        subarray.array_ = this;
        subarray.write_on_destruction_ = write_on_destruction;
        subarray.offset_ = start;
    }
            
    shape_type const & shape() const
    {
        return shape_;
    }
    
    bool isInside (shape_type const & p) const
    {
        for(int d=0; d<N; ++d)
            if(p[d] < 0 || p[d] >= shape_[d])
                return false;
        return true;
    }
    
    iterator begin()
    {
        return createCoupledIterator(*this);
    }
    
    iterator end()
    {
        return begin().getEndIterator();
    }
    
    // FIXME: this must be re-implemented in derived classes
    virtual void reshape(shape_type const & shape)
    {
        vigra_precondition(shape_ == shape_type(),
            "ChunkedArray::reshape(): can only reshape an empty array.");
        shape_ = shape;
    }
    
    shape_type shape_, chunk_shape_;
    Chunk * cache_first_, * cache_last_;
    std::size_t cache_size_, cache_max_size_;
    threading::mutex cache_lock_;
};

/** Returns a CoupledScanOrderIterator to simultaneously iterate over image m1 and its coordinates. 
 */
template <unsigned int N, class T>
typename ChunkedArray<N, T>::iterator
createCoupledIterator(ChunkedArray<N, T> & m)
{
    typedef typename ChunkedArray<N, T>::iterator    IteratorType;
    typedef typename IteratorType::handle_type           P1;
    typedef typename P1::base_type                       P0;
    
    return IteratorType(P1(m, 
                        P0(m.shape())));
}

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayFull
: public ChunkedArray<N, T>
{
  public:
        
    typedef MultiArray<N, T, Alloc> Storage;
    typedef typename Storage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayFull(shape_type const & shape, Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape),
      upper_bound_(shape + shape_type(1)),
      array_(shape, alloc)
    {
        this->chunk_shape_ = shape + shape_type(2);
    }
    
    ~ChunkedArrayFull()
    {}
    
    virtual pointer loadChunk(ChunkBase<N, T> *)
    {
        vigra_fail("ChunkedArrayFull::loadChunk() must not be called.");
        return 0;
    }
    
    virtual ChunkBase<N,T> * lookupChunk(shape_type const &)
    {
        vigra_fail("ChunkedArrayFull::lookupChunk() must not be called.");
        return 0;
    }
    
    virtual pointer ptr(shape_type const & point, 
                        shape_type & strides, shape_type & upper_bound, 
                        ChunkBase<N, T> **)
    {
        if(!this->isInside(point))
            return 0;

        strides = array_.stride();
        upper_bound = upper_bound_;
        return &array_[point];
    }
    
    Storage array_;  // a contiguous array
    shape_type upper_bound_;
};

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayLazy
: public ChunkedArray<N, T>
{
  public:
    
    class Chunk
    : public ChunkBase<N, T>
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Chunk(Alloc const & alloc = Alloc())
        : ChunkBase<N, T>(),
          alloc_(alloc)
        {}
        
        ~Chunk()
        {
            deallocate();
        }
    
        void reshape(shape_type const & shape)
        {
            vigra_precondition(this->pointer_ == (void*)0,
                "ChunkedArrayLazy::Chunk::reshape(): chunk was already allocated.");
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(shape);
        }
        
        pointer allocate()
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(!p)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                p = alloc_.allocate (s);
                try {
                    for (; i < s; ++i)
                        alloc_.construct (p + i, T());
                }
                catch (...) {
                    for (Size j = 0; j < i; ++j)
                        alloc_.destroy (p + j);
                    alloc_.deallocate (p, s);
                    throw;
                }
                this->pointer_.store(p, threading::memory_order_release);
            }
            return p;
        }
        
        void deallocate()
        {
            pointer p = this->pointer_.exchange(0, threading::memory_order_release);
            if(p)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                for (; i < s; ++i)
                    alloc_.destroy (p + i);
                alloc_.deallocate (p, s);
            }
        }
        
        Alloc alloc_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayLazy(shape_type const & shape, Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape),
      outer_array_(ChunkShape<N, T>::chunkArrayShape(shape), Chunk(alloc))
    {
        // set shape of the chunks
        typename ChunkStorage::iterator i   = outer_array_.begin(), 
                                        end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->reshape(min(this->chunk_shape_, 
                           shape - i.point()*this->chunk_shape_));
        }
    }
    
    ~ChunkedArrayLazy()
    {}
    
    virtual pointer loadChunk(ChunkBase<N, T> * chunk)
    {
        return static_cast<Chunk *>(chunk)->allocate();
    }
    
    virtual Chunk * lookupChunk(shape_type const & index)
    {
        return &outer_array_[index];
    }
    
    ChunkStorage outer_array_;  // the array of chunks
};

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayCompressed
: public ChunkedArray<N, T>
{
  public:
    
    class Chunk
    : public ChunkBase<N, T>
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Chunk()
        : ChunkBase<N, T>(),
          compressed_()
        {}
        
        ~Chunk()
        {
            deallocate();
        }
    
        void reshape(shape_type const & shape)
        {
            vigra_precondition(this->pointer_ == (void*)0,
                "ChunkedArrayCompressed::Chunk::reshape(): chunk was already allocated.");
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(shape);
        }
        
        pointer allocate()
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                p = alloc_.allocate (s);
                try {
                    for (; i < s; ++i)
                        alloc_.construct (p + i, T());
                }
                catch (...) {
                    for (Size j = 0; j < i; ++j)
                        alloc_.destroy (p + j);
                    alloc_.deallocate (p, s);
                    throw;
                }
                this->pointer_.store(p, threading::memory_order_release);
            }
            return p;
        }
        
        void deallocate()
        {
            pointer p = this->pointer_.exchange(0, threading::memory_order_release);
            if(p != 0)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                for (; i < s; ++i)
                    alloc_.destroy (p + i);
                alloc_.deallocate (p, s);
            }
            compressed_.clear();
        }
                
        void compress(CompressionMethod method)
        {
            pointer p = this->pointer_.exchange(0, threading::memory_order_release);
            if(p != 0)
            {
                vigra_invariant(compressed_.size() == 0,
                    "ChunkedArrayCompressed::Chunk::compress(): compressed and uncompressed pointer are both non-zero.");

                ::vigra::compress((char const *)p, this->size()*sizeof(T), compressed_, method);

                // std::cerr << "compression ratio: " << double(compressed_.size())/(this->size()*sizeof(T)) << "\n";
                std::size_t i = 0,
                            s = this->size();
                for (; i < s; ++i)
                    alloc_.destroy (p + i);
                alloc_.deallocate(p, (typename Alloc::size_type)this->size());
            }
        }
        
        pointer uncompress(CompressionMethod method)
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                if(compressed_.size())
                {
                    p = alloc_.allocate((typename Alloc::size_type)this->size());

                    ::vigra::uncompress(compressed_.data(), compressed_.size(), 
                                        (char*)p, this->size()*sizeof(T), method);
                    compressed_.clear();
                    this->pointer_.store(p, threading::memory_order_release);
                }
                else
                {
                    p = allocate();
                }
            }
            else
            {
                vigra_invariant(compressed_.size() == 0,
                    "ChunkedArrayCompressed::Chunk::uncompress(): compressed and uncompressed pointer are both non-zero.");
            }
            return p;
        }
        
        ArrayVector<char> compressed_;
        Alloc alloc_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayCompressed(shape_type const & shape, int cache_max = 0, CompressionMethod method=LZ4)
    : ChunkedArray<N, T>(shape, cache_max ? cache_max : max(outer_array_.shape())*2),
      outer_array_(ChunkShape<N, T>::chunkArrayShape(shape)),
      compression_method_(method)
    {
        // set shape of the chunks
        typename ChunkStorage::iterator i   = outer_array_.begin(), 
                                        end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->reshape(min(this->chunk_shape_, 
                           shape - i.point()*this->chunk_shape_));
        }
    }
    
    ~ChunkedArrayCompressed()
    {
        std::cerr << "    final cache size: " << this->cache_size_ << "\n";
    }
    
    virtual pointer loadChunk(ChunkBase<N, T> * chunk)
    {
        return static_cast<Chunk *>(chunk)->uncompress(compression_method_);
    }
    
    virtual void unloadChunk(ChunkBase<N, T> * chunk)
    {
        static_cast<Chunk *>(chunk)->compress(compression_method_);
    }
    
    virtual Chunk * lookupChunk(shape_type const & index)
    {
        return &outer_array_[index];
    }
        
    ChunkStorage outer_array_;  // the array of chunks

    CompressionMethod compression_method_;
};

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayHDF5
: public ChunkedArray<N, T>
{
  public:
    
    class Chunk
    : public ChunkBase<N, T>
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Chunk()
        : ChunkBase<N, T>(),
          array_(0)
        {}
        
        ~Chunk()
        {
            write();
        }
    
        void reshape(shape_type const & shape, shape_type const & start, 
                     ChunkedArrayHDF5 * array)
        {
            vigra_precondition(this->pointer_ == (void*)0,
                "ChunkedArrayCompressed::Chunk::reshape(): chunk was already allocated.");
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(shape);
            start_ = start;
            array_ = array;
        }
        
        void write()
        {
            pointer p = this->pointer_.exchange(0, threading::memory_order_release);
            if(p != 0)
            {
                herr_t status = array_->file_.writeBlock(array_->dataset_, start_, storage_);
                vigra_postcondition(status >= 0,
                    "ChunkedArrayHDF5: write to dataset failed.");
                storage_.swap(MultiArray<N, T>());
            }
        }
        
        pointer read()
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                storage_.reshape(this->shape_);
                herr_t status = array_->file_.readBlock(array_->dataset_, start_, this->shape_, storage_);
                vigra_postcondition(status >= 0,
                    "ChunkedArrayHDF5: read from dataset failed.");
                p = storage_.data();
                this->pointer_.store(p, threading::memory_order_release);
            }
            return p;
        }
        
        MultiArray<N, T> storage_;
        shape_type start_;
        ChunkedArrayHDF5 * array_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     shape_type const & shape, int cache_max = 0, int compression = 0)
    : ChunkedArray<N, T>(cache_max ? cache_max : max(outer_array_.shape())*2),
      file_(file),
      dataset_(),
      outer_array_()
    {
        if(file_.existsDataset(dataset))
        {
            ArrayVector<hsize_t> fileShape(file.getDatasetShape(dataset));
            
            // FIXME: the following checks don't work when T == TinyVector
            vigra_precondition(fileShape.size() == N,
                "ChunkedArrayHDF5(file, dataset, shape): dataset has wrong dimension.");
            for(unsigned int k=0; k<N; ++k)
                vigra_precondition(fileShape[k] == shape[k],
                    "ChunkedArrayHDF5(file, dataset, shape): shape mismatch between dataset and shape argument.");
            dataset_ = file.getDatasetHandleShared(dataset);
        }
        else
        {
            dataset_ = file_.createDataset(dataset, shape, T(), this->chunk_shape_, compression);
        }
        init(shape);
    }
    
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     int cache_max = 0, int compression = 0)
    : ChunkedArray<N, T>(cache_max ? cache_max : max(outer_array_.shape())*2),
      file_(file),
      dataset_(),
      outer_array_()
    {
        vigra_precondition(file.existsDataset(dataset),
            "ChunkedArrayHDF5(file, dataset): dataset does not exist.");
        dataset_ = file.getDatasetHandleShared(dataset);
            
        ArrayVector<hsize_t> fileShape(file.getDatasetShape(dataset));
        // FIXME: the following checks don't work when T == TinyVector
        vigra_precondition(fileShape.size() == N,
            "ChunkedArrayHDF5(file, dataset): dataset has wrong dimension.");
        shape_type shape(fileShape.begin());
        init(shape);
    }
    
    void init(shape_type const & shape)
    {
        reshape(shape);
        outer_array_.reshape(ChunkShape<N, T>::chunkArrayShape(shape));
        if(this->cache_max_size_ == 0)
            this->cache_max_size_ = max(outer_array_.shape())*2;
        
        // set shape of the chunks
        typename ChunkStorage::iterator i   = outer_array_.begin(), 
                                        end = outer_array_.end();
        for(; i != end; ++i)
        {
            shape_type start = i.point()*this->chunk_shape_;
            i->reshape(min(this->chunk_shape_, shape - start),
                       start,
                       this);
        }
    }
    
    ~ChunkedArrayHDF5()
    {
        std::cerr << "    final cache size: " << this->cache_size_ << "\n";
        
        // make sure that chunks are written to disk before the destructor of 
        // file_ is called
        outer_array_.swap(ChunkStorage());
    }
    
    void flushToDisk()
    {
        file_.flushToDisk();
    }
    
    virtual pointer loadChunk(ChunkBase<N, T> * chunk)
    {
        return static_cast<Chunk *>(chunk)->read();
    }
    
    virtual void unloadChunk(ChunkBase<N, T> * chunk)
    {
        static_cast<Chunk *>(chunk)->write();
    }
    
    virtual Chunk * lookupChunk(shape_type const & index)
    {
        return &outer_array_[index];
    }
    
    HDF5File file_;
    HDF5HandleShared dataset_;
    
    ChunkStorage outer_array_;  // the array of chunks
};

template <unsigned int N, class T>
class ChunkedArrayTmpFile
: public ChunkedArray<N, T>
{
  public:
#ifdef _WIN32
    typedef HANDLE FileHandle;
#else
    typedef int FileHandle;
#endif    
    
    class Chunk
    : public ChunkBase<N, T>
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Chunk()
        : ChunkBase<N, T>()
        {}
        
        ~Chunk()
        {
            unmap();
        }
        
        void reshape(shape_type const & shape, std::size_t alignment)
        {
            file_ = 0;
            offset_ = 0;
            std::size_t size = prod(shape)*sizeof(T);
            std::size_t mask = alignment - 1;
            alloc_size_ = (size + mask) & ~mask;
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(this->shape_);
        }
        
        void setFile(FileHandle file, std::ptrdiff_t offset)
        {
            file_ = file;
            offset_ = offset;
        }
        
        std::size_t size() const
        {
            return prod(this->shape_);
        }
        
        std::size_t alloc_size() const
        {
            return alloc_size_;
        }
        
        pointer map(std::ptrdiff_t offset)
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(!p)
            {
                if(offset_ < 0) 
                    offset_ = offset;
            #ifdef _WIN32
                static const std::size_t bits = sizeof(DWORD)*8,
                                         mask = (std::size_t(1) << bits) - 1;
                p = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                           std::size_t(offset_) >> bits, std::size_t(offset_) & mask, alloc_size_);
                if(!p)
                    winErrorToException("ChunkedArrayChunk::map(): ");
            #else
                p = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                  file_, std::size_t(offset_));
                if(!p)
                    throw std::runtime_error("ChunkedArrayChunk::map(): mmap() failed.");
            #endif
                this->pointer_.store(p, threading::memory_order_release);
            }
            return p;
        }
                
        void unmap()
        {
            void * p = this->pointer_.exchange(0, threading::memory_order_release);
            if(p)
            {
        #ifdef _WIN32
                ::UnmapViewOfFile(p);
        #else
                munmap(p, alloc_size_);
        #endif
            }
        }
        
        std::ptrdiff_t offset_;
        std::size_t alloc_size_;
        FileHandle file_;
    };

    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayTmpFile(shape_type const & shape, int cache_max=0, std::string const & path = "")
    : ChunkedArray<N, T>(shape, cache_max ? cache_max : max(outer_array_.shape())*2),
      outer_array_(ChunkShape<N, T>::chunkArrayShape(shape)),
      file_size_(),
      file_capacity_()
    {
        // set shape of the chunks
        typename ChunkStorage::iterator i = outer_array_.begin(), 
             end = outer_array_.end();
        std::size_t size = 0;
        for(; i != end; ++i)
        {
            i->reshape(min(this->chunk_shape_, shape - i.point()*this->chunk_shape_),
                       mmap_alignment);
            size += i->alloc_size();
        }
        
        std::cerr << "    file size: " << size << "\n";

    #ifdef VIGRA_NO_SPARSE_FILE
        file_capacity_ = 4*prod(this->chunk_shape_)*sizeof(T);
    #else
        file_capacity_ = size;
    #endif
    
    #ifdef _WIN32
        // create a temp file
        file_ = ::CreateFile(winTempFileName(path).c_str(), GENERIC_READ | GENERIC_WRITE,
                             0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE, NULL);
        if (file_ == INVALID_HANDLE_VALUE) 
            winErrorToException("ChunkedArrayTmpFile(): ");
            
        // make it a sparse file
        DWORD dwTemp;
        if(!::DeviceIoControl(file_, FSCTL_SET_SPARSE, NULL, 0, NULL, 0, &dwTemp, NULL))
            winErrorToException("ChunkedArrayTmpFile(): ");

        // place the data in the swap file
        // file_ = INVALID_HANDLE_VALUE;
        
        // resize and memory-map the file
        static const std::size_t bits = sizeof(LONG)*8, mask = (std::size_t(1) << bits) - 1;
        mappedFile_ = CreateFileMapping(file_, NULL, PAGE_READWRITE, 
                                        file_capacity_ >> bits, file_capacity_ & mask, NULL);
        if(!mappedFile_)
            winErrorToException("ChunkedArrayTmpFile(): ");
    #else
        mappedFile_ = file_ = fileno(tmpfile());
        if(file_ == -1)
            throw std::runtime_error("ChunkedArrayTmpFile(): unable to open file.");
        lseek(file_, file_capacity_-1, SEEK_SET);
        if(write(file_, "0", 1) == -1)
            throw std::runtime_error("ChunkedArrayTmpFile(): unable to resize file.");
    #endif
        
        // tell the chunks about the memory mapping
        i = outer_array_.begin();
        int offset = 0;
        for(; i != end; ++i)
        {
    #ifdef VIGRA_NO_SPARSE_FILE
            i->setFile(mappedFile_, -1);
    #else
            i->setFile(mappedFile_, offset);
            offset += i->alloc_size();
    #endif
        }
    }
    
    ~ChunkedArrayTmpFile()
    {
        std::cerr << "    final cache size: " << this->cache_size_ << "\n";
        unmap();
    #ifdef _WIN32
        ::CloseHandle(mappedFile_);
        ::CloseHandle(file_);
    #else
        ::close(file_);
    #endif
    }
    
    virtual pointer loadChunk(ChunkBase<N, T> * chunk_base)
    {
        Chunk * chunk = static_cast<Chunk *>(chunk_base);
        std::ptrdiff_t offset = file_size_;
    #ifdef VIGRA_NO_SPARSE_FILE
        std::size_t chunk_size = chunk->alloc_size();
        if(chunk->offset_ < 0)
        {
            if(offset + chunk_size > file_capacity_)
            {
                file_capacity_ *= 2;
                if(lseek(file_, file_capacity_-1, SEEK_SET) == -1)
                    throw std::runtime_error("ChunkedArrayTmpFile(): unable to reset file size.");
                if(write(file_, "0", 1) == -1)
                    throw std::runtime_error("ChunkedArrayTmpFile(): unable to resize file.");
            }
            file_size_ += chunk_size;
        }
    #endif
        return chunk->map(offset);
    }
    
    virtual void unloadChunk(ChunkBase<N, T> * chunk)
    {
        static_cast<Chunk *>(chunk)->unmap();
    }
    
    virtual Chunk * lookupChunk(shape_type const & index)
    {
        return &outer_array_[index];
    }
    
    void unmap()
    {
        typename ChunkStorage::iterator  i = outer_array_.begin(), 
             end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->unmap();
        }
    }
    
    ChunkStorage outer_array_;  // the array of chunks

    FileHandle file_, mappedFile_;  // the file back-end
    std::size_t file_size_, file_capacity_;
};

template <unsigned int N, class T>
class ChunkedHandle
{
  public:
    typedef ChunkedArray<N, T>             array_type;
    typedef ChunkBase<N, T>                    chunk_type;
    typedef typename MultiArrayShape<N>::type  shape_type;
    
    ChunkedHandle()
    : offset_(),
      chunk_(0)
    {}
    
    ChunkedHandle(ChunkedHandle const & other)
    : offset_(other.offset_),
      chunk_(0)
    {}
    
    ChunkedHandle & operator=(ChunkedHandle const & other)
    {
        offset_ = other.offset_;
        chunk_ = 0;
        return *this;
    }
    
    shape_type offset_;
    void * chunk_;
};

    /*
        The handle must store a pointer to a chunk because the chunk knows 
        about memory menagement, and to an array view because it knows about
        subarrays and slices.
        
        Perhaps we can reduce this to a single pointer or otherwise reduce 
        the handle memory to make it faster?
    */
template <class T, class NEXT>
class CoupledHandle<ChunkedMemory<T>, NEXT>
: public NEXT,
  public ChunkedHandle<NEXT::dimensions, T>
{
public:
    typedef NEXT                                  base_type;
    typedef ChunkedHandle<NEXT::dimensions, T>    base_type2;
    typedef CoupledHandle<ChunkedMemory<T>, NEXT> self_type;
    
    static const int index =                      NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef ChunkedArrayBase<dimensions, T>       array_type;
    typedef ChunkBase<dimensions, T>              chunk_type;
    typedef ChunkShape<dimensions, T>             chunk_shape;
    typedef T                                     value_type;
    typedef value_type *                          pointer;
    typedef value_type const *                    const_pointer;
    typedef value_type &                          reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;
    
    CoupledHandle()
    : base_type(),
      base_type2(),
      pointer_(), 
      array_()
    {}

    CoupledHandle(CoupledHandle const & other)
    : base_type(other),
      base_type2(other),
      pointer_(),
      array_(other.array_)
    {
        pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
    }

    CoupledHandle(array_type & array, NEXT const & next)
    : base_type(next),
      base_type2(),
      pointer_(), 
      array_(&array)
    {
        pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
    }

    ~CoupledHandle()
    {
        // deref the present chunk
        array_->unrefChunk(this);
    }

    CoupledHandle & operator=(CoupledHandle const & other)
    {
        if(this != &other)
        {
            // deref the present chunk
            array_->unrefChunk(this);
            base_type::operator=(other);
            base_type2::operator=(other);
            array_ = other.array_;
            pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
        }
        return *this;
    }
    
    using base_type::point;
    using base_type::shape;

    inline void incDim(int dim) 
    {
        base_type::incDim(dim);
        pointer_ += strides_[dim];
        if(point()[dim] == upper_bound_[dim])
        {
            // if(point()[dim] < shape()[dim])
                pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
        }
    }

    inline void decDim(int dim) 
    {
        base_type::decDim(dim);
        pointer_ -= strides_[dim];
        if(point()[dim] < upper_bound_[dim] - array_->chunk_shape_[dim])
        {
            // if(point()[dim] >= 0)
                pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
        }
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        base_type::addDim(dim, d);
        if(point()[dim] < shape()[dim] && point()[dim] >= 0)
            pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
    }

    inline void add(shape_type const & d) 
    {
        base_type::add(d);
        pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
    }
    
    template<int DIMENSION>
    inline void increment() 
    {
        // incDim(DIMENSION);
        base_type::template increment<DIMENSION>();
        pointer_ += strides_[DIMENSION];
        if(point()[DIMENSION] == upper_bound_[DIMENSION])
        {
            if(point()[DIMENSION] > shape()[DIMENSION])
                // this invariant check prevents the compiler from optimizing stupidly
                // (it makes a difference of a factor of 2!)                
                vigra_invariant(false, "CoupledHandle<ChunkedMemory<T>>: internal error.");
            else 
                pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
        }
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        // decDim(DIMENSION);
        base_type::template decrement<DIMENSION>();
        pointer_ -= strides_[DIMENSION];
        if(point()[DIMENSION] < upper_bound_[DIMENSION] - array_->chunk_shape_[DIMENSION])
        {
            if(point()[DIMENSION] < -1)
                // this invariant check prevents the compiler from optimizing stupidly
                // (it makes a difference of a factor of 2!)                
                vigra_invariant(false, "CoupledHandle<ChunkedMemory<T>>: internal error.");
            else
                pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
        }
    }
    
    template<int DIMENSION>
    inline void increment(MultiArrayIndex d) 
    {
        addDim(DIMENSION, d);
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex d) 
    {
        addDim(DIMENSION, -d);
    }
    
    void restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        base_type::restrictToSubarray(start, end);
        this->offset_ += start;
        pointer_ = array_->getChunk(point(), strides_, upper_bound_, this);
    }

    // ptr access
    reference operator*()
    {
        return *pointer_;
    }

    const_reference operator*() const
    {
        return *pointer_;
    }

    pointer operator->()
    {
        return pointer_;
    }

    const_pointer operator->() const
    {
        return pointer_;
    }

    pointer ptr()
    {
        return pointer_;
    }

    const_pointer ptr() const
    {
        return pointer_;
    }

    pointer pointer_;
    shape_type strides_, upper_bound_;
    array_type * array_;
};

} // namespace vigra

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HXX */
