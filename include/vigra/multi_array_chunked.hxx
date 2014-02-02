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
    compression:
    * SNAPPY is very fast, but cannot compress the float data
      (it achieves 5-10% for uint8 and 50% for double)
    * ZLIB achieves at least a factor of 5 better compression,
      but is more than an order of magnitude slower
    * HDF5 compression is already sufficient at level 1 (4-15%,
      higher levels don't lead to big gains) and a factor 3-10 slower 
      than without compression.
*/

#ifndef VIGRA_MULTI_ARRAY_CHUNKED_HXX
#define VIGRA_MULTI_ARRAY_CHUNKED_HXX

#include "metaprogramming.hxx"
#include "multi_array.hxx"
#include "threading.hxx"
#include "hdf5impex.hxx"

#ifdef _WIN32
# include "windows.h"
#else
# include <fcntl.h>
# include <stdlib.h>
# include <unistd.h>
# include <sys/stat.h>
# include <sys/mman.h>
#endif

#include <zlib.h>
#include <snappy.h>

#define VIGRA_CHUNKED_ARRAY_THREAD_SAFE

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

template <unsigned int N>
struct ChunkShape;

template <>
struct ChunkShape<1>
{
    static const unsigned int bits = 18;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct ChunkShape<2>
{
    static const unsigned int bits = 9;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct ChunkShape<3>
{
    static const unsigned int bits = 6;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
    
    static std::size_t unmangle(Shape3 const & s, Shape3 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)s[0] >> bits;
        UI o = (UI)s[0] & mask;
        b[1] = (UI)s[1] >> bits;
        o += ((UI)s[1] & mask) << bits;
        b[2] = (UI)s[2] >> bits;
        o += ((UI)s[2] & mask) << (2*bits);
        return o;
    }
    
    static void chunkIndex(Shape3 const & p, Shape3 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits;
        b[1] = (UI)p[1] >> bits;
        b[2] = (UI)p[2] >> bits;
    }
    
    static std::size_t chunkOffset(Shape3 const & p, Shape3 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits) * s[0] +
               ((UI)p[1] >> bits) * s[1] +
               ((UI)p[2] >> bits) * s[2];
    }
    
    static std::size_t offset(Shape3 const & p, Shape3 const & s)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & mask) * s[0] +
               ((UI)p[1] & mask) * s[1] +
               ((UI)p[2] & mask) * s[2];
    }
};

static int globalCount = 0;

template <class Shape>
Shape outerShape(Shape shape)
{
    static const int N = Shape::static_size;
    for(int k=0; k<N; ++k)
        shape[k] = (shape[k] + ChunkShape<N>::mask) >> ChunkShape<N>::bits;
    return shape;
}

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
      refcount_(0)
    {}
    
    ChunkBase(shape_type const & shape)
    : pointer_(0),
      shape_(shape),
      strides_(detail::defaultStride(shape)),
      refcount_(0)
    {}
    
    ChunkBase(ChunkBase const & rhs)
    : pointer_(0),
      shape_(rhs.shape_),
      strides_(rhs.strides_),
      refcount_(0)
    {}
    
    ChunkBase & operator=(ChunkBase const & rhs)
    {
        if(this != &rhs)
        {
            pointer p = pointer_.load(threading::memory_order_acquire);
            pointer rp = rhs.pointer_.load(threading::memory_order_acquire);
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
                pointer_.store(p, threading::memory_order_release);
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
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
    mutable threading::atomic<int> refcount_;
#else
    mutable int refcount_;
#endif
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
* HDF5 only works for scalar types so far
* public API for temp file arrays in swap missing
* support ZLIB compression levels
* allocators are not used
* the array implementations should go into cxx files in src/impex
  * this requires implementation of the low-level functions independently of dtype
    (use 'char *' and multiply shape and stride with sizeof(T))
* decide chunk locking policies for array views (in particular, for index access)
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
class ChunkedArrayBase
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
        
    ChunkedArrayBase()
    : shape_(0),
      chunk_shape_(ChunkShape<N>::value)
    {}
        
    ChunkedArrayBase(shape_type const & shape)
    : shape_(shape),
      chunk_shape_(ChunkShape<N>::value)
    {}
    
    virtual ~ChunkedArrayBase()
    {}
    
    virtual pointer ptr(shape_type const & point, 
                        shape_type & strides, shape_type & upper_bound, 
                        ChunkBase<N, T> ** chunk) = 0;
    
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
    
    void reshape(shape_type const & shape)
    {
        shape_ = shape;
    }
    
    shape_type shape_, chunk_shape_;
};

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayFull
: public ChunkedArrayBase<N, T>
{
  public:
        
    typedef MultiArray<N, T, Alloc> Storage;
    typedef typename Storage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayFull(shape_type const & shape, Alloc const & alloc = Alloc())
    : ChunkedArrayBase<N, T>(shape),
      upper_bound_(shape + shape_type(1)),
      array_(shape, alloc)
    {
        this->chunk_shape_ = shape + shape_type(2);
    }
    
    ~ChunkedArrayFull()
    {}
    
    // void map()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->map();
        // }
    // }
    
    // void unmap()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->unmap();
        // }
    // }

    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
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
class ChunkedArray
: public ChunkedArrayBase<N, T>
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
            vigra_precondition(pointer_ == (void*)0,
                "ChunkedArray::Chunk::reshape(): chunk was already allocated.");
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
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                for (; i < s; ++i)
                    alloc_.destroy (p + i);
                alloc_.deallocate (p, s);
                this->pointer_.store(0, threading::memory_order_release);
            }
        }
        
        Alloc alloc_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArray(shape_type const & shape, Alloc const & alloc = Alloc())
    : ChunkedArrayBase<N, T>(shape),
      outer_array_(outerShape(shape), Chunk(alloc))
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
    
    ~ChunkedArray()
    {}
    
    // void map()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->map();
        // }
    // }
    
    // void unmap()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->unmap();
        // }
    // }

    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & point, shape_type & strides, shape_type & upper_bound, ChunkBase<N, T> **)
    {
        if(!this->isInside(point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(point, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
                
        T * p = chunk.pointer_.load(threading::memory_order_acquire);
        if(p == 0)
        {
            // apply the double-checked locking pattern
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            p = chunk.allocate();
        }

        strides = chunk.strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(point, strides);
        return p + offset;
    }
    
    ChunkStorage outer_array_;  // the array of chunks
    threading::mutex cache_lock_;
};

enum ChunkCompression { ZLIB, SNAPPY };

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayCompressed
: public ChunkedArrayBase<N, T>
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
          compressed_(),
          cache_next_(0)
        {}
        
        ~Chunk()
        {
            deallocate();
        }
    
        void reshape(shape_type const & shape)
        {
            vigra_precondition(pointer_ == (void*)0,
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
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p != 0)
            {
                typedef typename Alloc::size_type Size;
                Size i = 0,
                     s = this->size();
                for (; i < s; ++i)
                    alloc_.destroy (p + i);
                alloc_.deallocate (p, s);
                this->pointer_.store(0, threading::memory_order_release);
            }
            compressed_.clear();
        }
        
        void compress(ChunkCompression method)
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p != 0)
            {
                vigra_invariant(compressed_.size() == 0,
                    "ChunkedArrayCompressed::Chunk::compress(): compressed and uncompressed pointer are both non-zero.");

                switch(method)
                {
                  case ZLIB:
                  {
                    uLong sourceLen = this->size()*sizeof(T),
                          destLen   = compressBound(sourceLen);
                    ArrayVector<Bytef> buffer(destLen);
                    int res = ::compress(buffer.data(), &destLen, (Bytef *)p, sourceLen);
                    vigra_postcondition(res == Z_OK,
                        "ChunkedArrayCompressed::Chunk::compress(): zlib compression failed.");                    
                    compressed_.insert(compressed_.begin(), 
                                       (char*)buffer.data(), (char*)buffer.data() + destLen);
                    break;
                  }
                  case SNAPPY:
                  {
                    size_t sourceLen = this->size()*sizeof(T),
                           destLen;
                    ArrayVector<char> buffer(snappy::MaxCompressedLength(sourceLen));
                    snappy::RawCompress((const char*)p, sourceLen, &buffer[0], &destLen);
                    compressed_.insert(compressed_.begin(), buffer.data(), buffer.data() + destLen);
                    break;
                  }
                }

                // std::cerr << "compression ratio: " << double(destLen) / sourceLen << "\n";
                alloc_.deallocate(p, (typename Alloc::size_type)this->size());
                this->pointer_.store(0, threading::memory_order_release);
            }
        }
        
        pointer uncompress(ChunkCompression method)
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                if(compressed_.size())
                {
                    p = alloc_.allocate((typename Alloc::size_type)this->size());

                    switch(method)
                    {
                      case ZLIB:
                      {
                        uLong sourceLen = compressed_.size(),
                              destLen   = this->size()*sizeof(T);
                        int res = ::uncompress((Bytef *)p, &destLen, (Bytef *)compressed_.data(), sourceLen);
                        vigra_postcondition(res == Z_OK,
                            "ChunkedArrayCompressed::Chunk::uncompress(): zlib decompression failed.");
                        break;
                      }
                      case SNAPPY:
                      {
                        snappy::RawUncompress(compressed_.data(), compressed_.size(), (char *)p);
                        break;
                      }
                    }
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
        Chunk * cache_next_;
        Alloc alloc_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayCompressed(shape_type const & shape, int cache_max = 0, ChunkCompression method=SNAPPY)
    : ChunkedArrayBase<N, T>(shape),
      outer_array_(outerShape(shape)),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      cache_max_size_(cache_max ? cache_max : max(outer_array_.shape())*2),
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
        std::cerr << "    final cache size: " << cache_size_ << "\n";
    }
    
    // void map()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->map();
        // }
    // }
    
    // void unmap()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->unmap();
        // }
    // }

    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & point, shape_type & strides, shape_type & upper_bound, ChunkBase<N, T> ** chunk_ptr)
    {
        if(*chunk_ptr)
        {
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
            (*chunk_ptr)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
#else
            --(*chunk_ptr)->refcount_;
#endif
            *chunk_ptr = 0;
        }
        
        if(!this->isInside(point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(point, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
        
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
        // Obtain a reference to the current chunk.
        // We use a simple spin-lock here because it is very fast in case of success
        // and failures (i.e. a collisions with another thread) are presumably 
        // very rare.
        while(true)
        {
            int rc = chunk.refcount_.load(threading::memory_order_acquire);
            if(rc < 0)
            {
                // cache management in progress => try again later
                threading::this_thread::yield();
            }
            else if(chunk.refcount_.compare_exchange_weak(rc, rc+1, threading::memory_order_seq_cst))
            {
                // success
                break;
            }
        }
#else
        ++chunk.refcount_;
#endif
        
        T * p = chunk.pointer_.load(threading::memory_order_acquire);
        if(p == 0)
        {
            // apply the double-checked locking pattern
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            p = chunk.pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                p = chunk.uncompress(compression_method_);
                
                // insert in queue of mapped chunks
                if(!cache_first_)
                    cache_first_ = &chunk;
                if(cache_last_)
                    cache_last_->cache_next_ = &chunk;
                cache_last_ = &chunk;
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
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                    if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_seq_cst))
#else
                    if(c->refcount_ == 0)
#endif
                    {
                        // refcount was zero => can compress
                        c->compress(compression_method_);
                        c->cache_next_ = 0;
                        --cache_size_;
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                        c->refcount_.store(0, std::memory_order_release);
#endif
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

        strides = chunk.strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(point, strides);
        *chunk_ptr = &chunk;
        return p + offset;
    }
    
    ChunkStorage outer_array_;  // the array of chunks

    Chunk * cache_first_, * cache_last_;  // cache of loaded chunks
    std::size_t cache_size_, cache_max_size_;
    ChunkCompression compression_method_;
    threading::mutex cache_lock_;
};

template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayHDF5
: public ChunkedArrayBase<N, T>
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
          array_(0),
          cache_next_(0)
        {}
        
        ~Chunk()
        {
            write();
        }
    
        void reshape(shape_type const & shape, shape_type const & start, 
                     ChunkedArrayHDF5 * array)
        {
            vigra_precondition(pointer_ == (void*)0,
                "ChunkedArrayCompressed::Chunk::reshape(): chunk was already allocated.");
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(shape);
            start_ = start;
            array_ = array;
        }
        
        void write()
        {
            pointer p = this->pointer_.load(threading::memory_order_acquire);
            if(p != 0)
            {
                herr_t status = array_->file_.writeBlock(array_->dataset_, start_, storage_);
                vigra_postcondition(status >= 0,
                    "ChunkedArrayHDF5: write to dataset failed.");
                storage_.swap(MultiArray<N, T>());
                this->pointer_.store(0, threading::memory_order_release);
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
        Chunk * cache_next_;
    };
    
    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     shape_type const & shape, int cache_max = 0, int compression = 0)
    : ChunkedArrayBase<N, T>(),
      file_(file),
      dataset_(),
      outer_array_(),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      cache_max_size_(cache_max ? cache_max : max(outer_array_.shape())*2)
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
    : ChunkedArrayBase<N, T>(),
      file_(file),
      dataset_(),
      outer_array_(),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      cache_max_size_(cache_max ? cache_max : max(outer_array_.shape())*2)
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
        outer_array_.reshape(outerShape(shape));
        if(cache_max_size_ == 0)
            cache_max_size_ = max(outer_array_.shape())*2;
        
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
        std::cerr << "    final cache size: " << cache_size_ << "\n";
        
        // make sure that chunks are written to disk before the destructor of 
        // file_ is called
        outer_array_.swap(ChunkStorage());
    }
    
    void flushToDisk()
    {
        file_.flushToDisk();
    }
    
    // void map()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->map();
        // }
    // }
    
    // void unmap()
    // {
        // typename ChunkStorage::iterator  i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
        // {
            // i->unmap();
        // }
    // }

    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & point, shape_type & strides, shape_type & upper_bound, ChunkBase<N, T> ** chunk_ptr)
    {
        if(*chunk_ptr)
        {
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
            (*chunk_ptr)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
#else
            --(*chunk_ptr)->refcount_;
#endif
            *chunk_ptr = 0;
        }
        
        if(!this->isInside(point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(point, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
        
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
        // Obtain a reference to the current chunk.
        // We use a simple spin-lock here because it is very fast in case of success
        // and failures (i.e. a collisions with another thread) are presumably 
        // very rare.
        while(true)
        {
            int rc = chunk.refcount_.load(threading::memory_order_acquire);
            if(rc < 0)
            {
                // cache management in progress => try again later
                threading::this_thread::yield();
            }
            else if(chunk.refcount_.compare_exchange_weak(rc, rc+1, threading::memory_order_seq_cst))
            {
                // success
                break;
            }
        }
#else
        ++chunk.refcount_;
#endif
        
        T * p = chunk.pointer_.load(threading::memory_order_acquire);
        if(p == 0)
        {
            // apply the double-checked locking pattern
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            p = chunk.pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                p = chunk.read();
                
                // insert in queue of mapped chunks
                if(!cache_first_)
                    cache_first_ = &chunk;
                if(cache_last_)
                    cache_last_->cache_next_ = &chunk;
                cache_last_ = &chunk;
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
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                    if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_seq_cst))
#else
                    if(c->refcount_ == 0)
#endif
                    {
                        // refcount was zero => can compress
                        c->write();
                        c->cache_next_ = 0;
                        --cache_size_;
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                        c->refcount_.store(0, std::memory_order_release);
#endif
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

        strides = chunk.strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(point, strides);
        *chunk_ptr = &chunk;
        return p + offset;
    }
    
    HDF5File file_;
    HDF5HandleShared dataset_;
    
    ChunkStorage outer_array_;  // the array of chunks

    Chunk * cache_first_, * cache_last_;  // cache of loaded chunks
    std::size_t cache_size_, cache_max_size_;

    threading::mutex cache_lock_;
};

template <unsigned int N, class T>
class ChunkedArrayTmpFile
: public ChunkedArrayBase<N, T>
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
        : ChunkBase<N, T>(),
          cache_next_(0)
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
            pointer p = this->pointer_.load();
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
        Chunk * cache_next_;
    };

    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayTmpFile(shape_type const & shape, int cache_max=0, std::string const & path = "")
    : ChunkedArrayBase<N, T>(shape),
      outer_array_(outerShape(shape)),
      file_size_(),
      file_capacity_(),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      // cache_max_size_(max(outer_array_.shape())*2)
      cache_max_size_(cache_max ? cache_max : max(outer_array_.shape())*2)
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
        std::cerr << "    final cache size: " << cache_size_ << "\n";
        unmap();
    #ifdef _WIN32
        ::CloseHandle(mappedFile_);
        ::CloseHandle(file_);
    #else
        close(file_);
    #endif
    }
    
    void map()
    {
        typename ChunkStorage::iterator  i = outer_array_.begin(), 
             end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->map();
        }
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

    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & point, shape_type & strides, shape_type & upper_bound, ChunkBase<N, T> ** chunk_ptr)
    {
        // USETICTOC;
        // TIC;
        if(*chunk_ptr)
        {
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
            int old_rc = (*chunk_ptr)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
            if(old_rc == 0)
                vigra_invariant(0, "refcount got negative!");
#else
            --(*chunk_ptr)->refcount_;
#endif
            *chunk_ptr = 0;
        }
        
        if(!this->isInside(point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(point, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
        
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
        // Obtain a reference to the current chunk.
        // We use a simple spin-lock here because it is very fast in case of success
        // and failures (i.e. a collisions with another thread) are presumably 
        // very rare.
        while(true)
        {
            int rc = chunk.refcount_.load(threading::memory_order_acquire);
            if(rc < 0)
            {
                // cache management in progress => try again later
                threading::this_thread::yield();
            }
            else if(chunk.refcount_.compare_exchange_weak(rc, rc+1, threading::memory_order_seq_cst))
            {
                // success
                break;
            }
        }
#else
        ++chunk.refcount_;
#endif
        
        T * p = chunk.pointer_.load(threading::memory_order_acquire);
        if(p == 0)
        {
            // apply the double-checked locking pattern
            threading::lock_guard<threading::mutex> guard(cache_lock_);
            p = chunk.pointer_.load(threading::memory_order_acquire);
            if(p == 0)
            {
                std::ptrdiff_t offset = file_size_;
            #ifdef VIGRA_NO_SPARSE_FILE
                std::size_t chunk_size = chunk.alloc_size();
                if(chunk.offset_ < 0)
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
                p = chunk.map(offset);
                
                // insert in queue of mapped chunks
                if(!cache_first_)
                    cache_first_ = &chunk;
                if(cache_last_)
                    cache_last_->cache_next_ = &chunk;
                cache_last_ = &chunk;
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
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                    if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_seq_cst))
#else
                    if(c->refcount_ == 0)
#endif
                    {
                        // refcount was zero => can unmap
                        c->unmap();
                        c->cache_next_ = 0;
                        --cache_size_;
#ifdef VIGRA_CHUNKED_ARRAY_THREAD_SAFE
                        c->refcount_.store(0, std::memory_order_release);
#endif
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

        strides = chunk.strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(point, strides);
        *chunk_ptr = &chunk;
        // timeit += TOCN;
        return p + offset;
    }
    
    ChunkStorage outer_array_;  // the array of chunks

    FileHandle file_, mappedFile_;  // the file back-end
    std::size_t file_size_, file_capacity_;
    
    Chunk * cache_first_, * cache_last_;  // cache of loaded chunks
    std::size_t cache_size_, cache_max_size_;
    threading::mutex cache_lock_;
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
: public NEXT
{
public:
    typedef NEXT                                  base_type;
    typedef CoupledHandle<ChunkedMemory<T>, NEXT> self_type;
    
    static const int index =                      NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef ChunkedArrayBase<dimensions, T>           array_type;
    typedef ChunkBase<dimensions, T>              chunk_type;
    typedef ChunkShape<dimensions>                chunk_shape;
    typedef T                                     value_type;
    typedef value_type *                          pointer;
    typedef value_type const *                    const_pointer;
    typedef value_type &                          reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;
    
    typedef pointer (*SetPointerFct)(array_type * array, shape_type const & p, 
                                     shape_type & strides, shape_type & border, chunk_type ** chunk);
    
    static SetPointerFct setPointer;
    
    static pointer setPointerFct(array_type * array, shape_type const & p, 
                                 shape_type & strides, shape_type & border, chunk_type ** chunk)
    {
        return array->ptr(p, strides, border, chunk);
    }

    CoupledHandle()
    : base_type(),
      pointer_(), 
      array_()
    {}

    CoupledHandle(CoupledHandle const & other)
    : base_type(other),
      pointer_(), 
      array_(other.array_),
      chunk_()
    {
        pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
    }

    CoupledHandle(array_type & array, NEXT const & next)
    : base_type(next),
      pointer_(), 
      array_(&array),
      chunk_()
    {
        pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
    }

    ~CoupledHandle()
    {
        // deref the present chunk
        array_->ptr(shape_type(-1), strides_, upper_bound_, &chunk_);
    }

    CoupledHandle & operator=(CoupledHandle const & other)
    {
        if(this != &other)
        {
            // deref the present chunk
             array_->ptr(shape_type(-1), strides_, upper_bound_, &chunk_);
             base_type::operator=(other);
             array_ = other.array_;
             pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
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
                pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
        }
    }

    inline void decDim(int dim) 
    {
        base_type::decDim(dim);
        pointer_ -= strides_[dim];
        if(point()[dim] < upper_bound_[dim] - array_->chunk_shape_[dim])
        {
            // if(point()[dim] >= 0)
                pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
        }
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        base_type::addDim(dim, d);
        if(point()[dim] < shape()[dim] && point()[dim] >= 0)
            pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
    }

    inline void add(shape_type const & d) 
    {
        base_type::add(d);
        pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
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
                pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
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
                pointer_ = array_->ptr(point(), strides_, upper_bound_, &chunk_);
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
    
    // void restrictToSubarray(shape_type const & start, shape_type const & end)
    // {
        // view_.unsafePtr() += dot(start, strides_);
        // base_type::restrictToSubarray(start, end);
    // }

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
    chunk_type * chunk_;
};

template <class T, class NEXT>
typename CoupledHandle<ChunkedMemory<T>, NEXT>::SetPointerFct 
CoupledHandle<ChunkedMemory<T>, NEXT>::setPointer = &CoupledHandle<ChunkedMemory<T>, NEXT>::setPointerFct;

// template <class T, class NEXT>
// class CoupledHandle<ChunkedMemory<T>, NEXT>
// : public NEXT
// {
// public:
    // typedef NEXT                                  base_type;
    // typedef CoupledHandle<ChunkedMemory<T>, NEXT> self_type;
    
    // static const int index =                      NEXT::index + 1;    // index of this member of the chain
    // static const unsigned int dimensions =        NEXT::dimensions;

    // typedef ChunkedArrayBase<dimensions, T>           array_type;
    // typedef ChunkShape<dimensions>                chunk_shape;
    // typedef T                                     value_type;
    // typedef value_type *                          pointer;
    // typedef value_type const *                    const_pointer;
    // typedef value_type &                          reference;
    // typedef value_type const &                    const_reference;
    // typedef typename base_type::shape_type        shape_type;
    
    // typedef T* (*SetPointerFct)(shape_type const & p, array_type * array, shape_type & strides, shape_type & border);
    
    // static SetPointerFct setPointer;
    
    // static T* setPointerFct(shape_type const & p, array_type * array, shape_type & strides, shape_type & border)
    // {
        // return array->ptr(p, strides, border);
    // }

    // CoupledHandle()
    // : base_type(),
      // pointer_(), 
      // array_()
    // {}

    // CoupledHandle(array_type & array, NEXT const & next)
    // : base_type(next),
      // pointer_((*setPointer)(shape_type(), &array, strides_, border_)), 
      // array_(&array)
    // {}
    
    // using base_type::point;
    // using base_type::shape;

    // inline void incDim(int dim) 
    // {
        // base_type::incDim(dim);
        // // if((point()[dim] & chunk_shape::mask) == 0)
        // // {
            // // if(point()[dim] < shape()[dim])
                // // pointer_ = (*setPointer)(point(), array_, strides_, border_);
        // // }
        // // else
        // // {
            // // pointer_ += strides_[dim];
        // // }
        // pointer_ += strides_[dim];
        // if(point()[dim] == border_[dim])
        // {
            // if(point()[dim] < shape()[dim])
                // // pointer_ = (*setPointer)(point(), array_, strides_, border_);
                // pointer_ = array_->ptr(point(), strides_, border_);
        // }
    // }

    // // inline void decDim(int dim) 
    // // {
        // // view_.unsafePtr() -= strides_[dim];
        // // base_type::decDim(dim);
    // // }

    // inline void addDim(int dim, MultiArrayIndex d) 
    // {
        // base_type::addDim(dim, d);
        // if(point()[dim] < shape()[dim])
            // pointer_ = (*setPointer)(point(), array_, strides_, border_);
        // // if(dim == 0)
            // // std::cerr << point() << " " << pointer_ << " " << strides_ << " " << border_ << "\n";
    // }

    // // inline void add(shape_type const & d) 
    // // {
        // // view_.unsafePtr() += dot(d, strides_);
        // // base_type::add(d);
    // // }
    
    // // template<int DIMENSION>
    // // inline void increment() 
    // // {
        // // view_.unsafePtr() += strides_[DIMENSION];
        // // base_type::template increment<DIMENSION>();
    // // }
    
    // // template<int DIMENSION>
    // // inline void decrement() 
    // // {
        // // view_.unsafePtr() -= strides_[DIMENSION];
        // // base_type::template decrement<DIMENSION>();
    // // }
    
    // // // TODO: test if making the above a default case of the this hurts performance
    // // template<int DIMENSION>
    // // inline void increment(MultiArrayIndex offset) 
    // // {
        // // view_.unsafePtr() += offset*strides_[DIMENSION];
        // // base_type::template increment<DIMENSION>(offset);
    // // }
    
    // // template<int DIMENSION>
    // // inline void decrement(MultiArrayIndex offset) 
    // // {
        // // view_.unsafePtr() -= offset*strides_[DIMENSION];
        // // base_type::template decrement<DIMENSION>(offset);
    // // }
    
    // // void restrictToSubarray(shape_type const & start, shape_type const & end)
    // // {
        // // view_.unsafePtr() += dot(start, strides_);
        // // base_type::restrictToSubarray(start, end);
    // // }

    // // ptr access
    // reference operator*()
    // {
        // return *pointer_;
    // }

    // const_reference operator*() const
    // {
        // return *pointer_;
    // }

    // pointer operator->()
    // {
        // return pointer_;
    // }

    // const_pointer operator->() const
    // {
        // return pointer_;
    // }

    // pointer ptr()
    // {
        // return pointer_;
    // }

    // const_pointer ptr() const
    // {
        // return pointer_;
    // }

    // pointer pointer_;
    // shape_type strides_, border_;
    // array_type * array_;
// };

// template <class T, class NEXT>
// typename CoupledHandle<ChunkedMemory<T>, NEXT>::SetPointerFct 
// CoupledHandle<ChunkedMemory<T>, NEXT>::setPointer = &CoupledHandle<ChunkedMemory<T>, NEXT>::setPointerFct;


} // namespace vigra

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HXX */
