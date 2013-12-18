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

#ifndef VIGRA_MULTI_ARRAY_CHUNKED_HXX
#define VIGRA_MULTI_ARRAY_CHUNKED_HXX

#include "metaprogramming.hxx"
#include "multi_array.hxx"
#include "threading.hxx"

#ifdef _WIN32
# include "windows.h"
#else
# include <fcntl.h>
# include <stdlib.h>
# include <unistd.h>
# include <sys/stat.h>
# include <sys/mman.h>
#endif

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
size_t winClusterSize()
{
    SYSTEM_INFO info;
    ::GetSystemInfo(&info); 
    return info.dwAllocationGranularity;
}

#endif

namespace {

#ifdef _WIN32
size_t mmap_alignment = winClusterSize();
#else
size_t mmap_alignment = sysconf(_SC_PAGE_SIZE);
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

#if 0
template <unsigned int N, class T>
class ChunkedArrayChunk
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
#ifndef _WIN32
    typedef int HANDLE;
#endif    
    
    ChunkedArrayChunk()
    : pointer_()
    {}
    
    ~ChunkedArrayChunk()
    {
        unmap();
    }
    
    void reshape(shape_type const & shape)
    {
        file_ = 0;
        offset_ = 0;
        size_t size = prod(shape)*sizeof(T);
        alloc_size_ = (size + alloc_mask_) & ~alloc_mask_;
        shape_ = shape;
        strides_ = detail::defaultStride(shape_);
    }
    
    void setFile(HANDLE file, size_t offset)
    {
        file_ = file;
        offset_ = offset;
    }
    
    size_t size() const
    {
        return prod(shape_);
    }
    
    size_t alloc_size() const
    {
        return alloc_size_;
    }
    
    void map()
    {
        if(!pointer_)
        {
#ifdef _WIN32
            static const size_t bits = sizeof(DWORD)*8,
                                mask = (size_t(1) << bits) - 1;
            pointer_ = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                              offset_ >> bits, offset_ & mask, alloc_size_);
            if(!pointer_)
                winErrorToException("ChunkedArrayChunk::map(): ");
#else
            pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                     file_, offset_);
            if(!pointer_)
                throw std::runtime_error("ChunkedArrayChunk::map(): mmap() failed.");
#endif
       }
    }
    
    void unmap()
    {
#ifdef _WIN32
        if(pointer_)
            ::UnmapViewOfFile((void*)pointer_);
#else
        if(pointer_)
            munmap((void*)pointer_, alloc_size_);
#endif
        pointer_ = 0;
    }
    
    pointer pointer_;
    size_t offset_, alloc_size_;
    shape_type shape_, strides_;
    HANDLE file_;
    
    static size_t alloc_mask_;
};

#ifdef _WIN32
    template <unsigned int N, class T>
    size_t ChunkedArrayChunk<N, T>::alloc_mask_ = winClusterSize() - 1;
#else
    template <unsigned int N, class T>
    size_t ChunkedArrayChunk<N, T>::alloc_mask_ = sysconf(_SC_PAGE_SIZE) - 1;
#endif

#endif

// template <unsigned int N, class T>
// class ChunkedArray
// {
  // public:
    // typedef MultiArray<N, MultiArray<N, T> > ChunkStorage;
    // typedef typename ChunkStorage::difference_type  shape_type;
    // typedef T value_type;
    // typedef value_type * pointer;
    // typedef value_type & reference;
    
    // ChunkedArray(shape_type const & shape)
    // : shape_(shape),
      // default_chunk_shape_(ChunkShape<N>::value),
      // outer_array_(outerShape(shape))
    // {
        // auto i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
            // i->reshape(min(default_chunk_shape_, shape - i.point()));
    // }
    
    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    // pointer ptr(shape_type const & p, shape_type & strides, shape_type & border)
    // {
        // // if(!isInside(p))
            // // return 0;
        
        // shape_type chunkIndex(SkipInitialization);
        // ChunkShape<N>::chunkIndex(p, chunkIndex);
        // MultiArray<N, T> & chunk = outer_array_[chunkIndex];
        // strides = chunk.stride();
        // border = (chunkIndex + shape_type(1)) * default_chunk_shape_;
        // std::size_t offset = ChunkShape<N>::offset(p, strides);
        // return chunk.data() + offset;
    // }
    
    // shape_type const & shape() const
    // {
        // return shape_;
    // }
    
    // bool isInside (shape_type const & p) const
    // {
        // for(int d=0; d<N; ++d)
            // if(p[d] < 0 || p[d] >= shape_[d])
                // return false;
        // return true;
    // }
    
    // shape_type shape_, default_chunk_shape_;
    // ChunkStorage outer_array_;
// };

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
    
    ChunkBase(ChunkBase const & rhs)
    : pointer_(0),
      shape_(rhs.shape_),
      strides_(rhs.strides_),
      refcount_(0)
    {}
    
    threading::atomic<pointer> pointer_;
    shape_type shape_, strides_;
    mutable threading::atomic<int> refcount_;
};


/*
The present implementation uses a memory-mapped sparse file to store the chunks.
A sparse file is created on Linux using the O_TRUNC flag (possibly, it is even 
the default file behavior on Linux => check), and on Windows by
calling DeviceIoControl(file_handle, FSCTL_SET_SPARSE,...) after file creation.

We can automatically delete the file upon closing. On Windows, this happens
if the file was opened with FILE_FLAG_DELETE_ON_CLOSE (it may be useful to
combine this with the flag FILE_ATTRIBUTE_TEMPORARY, which tells the OS to
avoid writing the file to disk if possible). On Linux, you can call
fileno(tmpfile()) for the same purpose.

Alternatives are:
* Keep the array in memory, but compress unused chunks.
* Don't create a file explicitly, but use the swap file instead. This is 
  achieved on Linux by mmap(..., MAP_PRIVATE | MAP_ANONYMOUS, -1, ...), 
  on Windows by calling CreateFileMapping(INVALID_HANDLE_VALUE, ...).
* Back chunks by HDF5 chunks, possibly using on-the-fly compression. This
  is in particular useful for existing HDF5 files.
* Back chunks by HDF5 datasets. This can be combined with compression 
  (both explicit and on-the-fly).
*/
template <unsigned int N, class T>
class ChunkedArray
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
        
    ChunkedArray(shape_type const & shape)
    : shape_(shape),
      default_chunk_shape_(ChunkShape<N>::value)
    {}
    
    virtual ~ChunkedArray()
    {}
    
    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & p, shape_type & strides, shape_type & border, ChunkBase<N, T> ** chunk) = 0;
    
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
    
    shape_type shape_, default_chunk_shape_;
};

template <unsigned int N, class T>
class ChunkedArrayFile
: public ChunkedArray<N, T>
{
  public:
#ifndef _WIN32
    typedef int HANDLE;
#endif    
    
    class Chunk
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Chunk()
        : pointer_()
        {}
        
        ~Chunk()
        {
            unmap();
        }
        
        void reshape(shape_type const & shape, size_t alignment)
        {
            file_ = 0;
            offset_ = 0;
            size_t size = prod(shape)*sizeof(T);
            size_t mask = alignment - 1;
            alloc_size_ = (size + mask) & ~mask;
            shape_ = shape;
            strides_ = detail::defaultStride(shape_);
        }
        
        void setFile(HANDLE file, size_t offset)
        {
            file_ = file;
            offset_ = offset;
        }
        
        size_t size() const
        {
            return prod(shape_);
        }
        
        size_t alloc_size() const
        {
            return alloc_size_;
        }
        
        void map()
        {
            if(!pointer_)
            {
    #ifdef _WIN32
                static const size_t bits = sizeof(DWORD)*8,
                                    mask = (size_t(1) << bits) - 1;
                pointer_ = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                                  offset_ >> bits, offset_ & mask, alloc_size_);
                if(!pointer_)
                    winErrorToException("ChunkedArrayChunk::map(): ");
    #else
                pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                         file_, offset_);
                if(!pointer_)
                    throw std::runtime_error("ChunkedArrayChunk::map(): mmap() failed.");
    #endif
           }
        }
        
        void unmap()
        {
    #ifdef _WIN32
            if(pointer_)
                ::UnmapViewOfFile((void*)pointer_);
    #else
            if(pointer_)
                munmap((void*)pointer_, alloc_size_);
    #endif
            pointer_ = 0;
        }
        
        pointer pointer_;
        size_t offset_, alloc_size_;
        shape_type shape_, strides_;
        HANDLE file_;
    };

    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayFile(shape_type const & shape, char const * name)
    : ChunkedArray<N, T>(shape),
      outer_array_(outerShape(shape))
    {
        // set shape of the chunks
        typename ChunkStorage::iterator i = outer_array_.begin(), 
             end = outer_array_.end();
        size_t size = 0;
        for(; i != end; ++i)
        {
            i->reshape(min(this->default_chunk_shape_, shape - i.point()*this->default_chunk_shape_),
                       mmap_alignment);
            size += i->alloc_size();
        }
        
        std::cerr << "file size: " << size << "\n";

#ifdef _WIN32
        // create file or open it when it already exists
        file_ = ::CreateFile(name, GENERIC_READ | GENERIC_WRITE,
                             0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (file_ == INVALID_HANDLE_VALUE) 
            winErrorToException("ChunkedArrayFile(): ");
        
        bool isNewFile = ::GetLastError() != ERROR_ALREADY_EXISTS;
        if(isNewFile)
        {
            std::cerr << "creating new file\n";
            
            // make it a sparse file
            DWORD dwTemp;
            if(!::DeviceIoControl(file_, FSCTL_SET_SPARSE, NULL, 0, NULL, 0, &dwTemp, NULL))
                winErrorToException("ChunkedArrayFile(): ");

            // set the file size
            static const size_t bits = sizeof(LONG)*8,
                                mask = (size_t(1) << bits) - 1;
            LONG fileSizeHigh = size >> bits;
            if(::SetFilePointer(file_, size & mask, &fileSizeHigh, FILE_BEGIN) == INVALID_SET_FILE_POINTER)
                winErrorToException("ChunkedArrayFile(): ");
            if(!::SetEndOfFile(file_)) 
                winErrorToException("ChunkedArrayFile(): ");
        }
        
        // memory-map the file
        mappedFile_ = CreateFileMapping(file_, NULL, PAGE_READWRITE, 0, 0, NULL);
        if(!mappedFile_)
            winErrorToException("ChunkedArrayFile(): ");
#else
        mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
        mappedFile_ = file_ = open(name, O_RDWR | O_CREAT | O_TRUNC, mode);
        if(file_ == -1)
            throw std::runtime_error("ChunkedArrayFile(): unable to open file.");
        lseek(file_, size-1, SEEK_SET);
        if(write(file_, "0", 1) == -1)
            throw std::runtime_error("ChunkedArrayFile(): unable to resize file.");
#endif
        
        // tell the chunks about the memory mapping
        i = outer_array_.begin();
        int offset = 0;
        for(; i != end; ++i)
        {
            i->setFile(mappedFile_, offset);
            offset += i->alloc_size();
        }
    }
    
    ~ChunkedArrayFile()
    {
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
    
    virtual pointer ptr(shape_type const & p, shape_type & strides, shape_type & border)
    {
        if(!this->isInside(p))
            return 0;
        
        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(p, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
        chunk.map();
        strides = chunk.strides_;
        border = (chunkIndex + shape_type(1)) * this->default_chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(p, strides);
        return chunk.pointer_ + offset;
    }
    
    HANDLE file_, mappedFile_;
    ChunkStorage outer_array_;
};

template <unsigned int N, class T>
class ChunkedArrayTmpFile
: public ChunkedArray<N, T>
{
  public:
#ifndef _WIN32
    typedef int HANDLE;
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
        
        void reshape(shape_type const & shape, size_t alignment)
        {
            file_ = 0;
            offset_ = 0;
            size_t size = prod(shape)*sizeof(T);
            size_t mask = alignment - 1;
            alloc_size_ = (size + mask) & ~mask;
            this->shape_ = shape;
            this->strides_ = detail::defaultStride(shape_);
        }
        
        void setFile(HANDLE file, ptrdiff_t offset)
        {
            file_ = file;
            offset_ = offset;
        }
        
        size_t size() const
        {
            return prod(this->shape_);
        }
        
        size_t alloc_size() const
        {
            return alloc_size_;
        }
        
        pointer map(ptrdiff_t offset)
        {
            pointer p = this->pointer_.load();
            if(!p)
            {
                if(offset_ < 0) 
                    offset_ = offset;
            #ifdef _WIN32
                static const size_t bits = sizeof(DWORD)*8,
                                    mask = (size_t(1) << bits) - 1;
                p = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                           size_t(offset_) >> bits, size_t(offset_) & mask, alloc_size_);
                if(!p)
                    winErrorToException("ChunkedArrayChunk::map(): ");
            #else
                p = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                  file_, size_t(offset_));
                if(!p)
                    throw std::runtime_error("ChunkedArrayChunk::map(): mmap() failed.");
            #endif
                this->pointer_.store(p);
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
        
        ptrdiff_t offset_;
        size_t alloc_size_;
        HANDLE file_;
        Chunk * cache_next_;
    };

    typedef MultiArray<N, Chunk> ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    ChunkedArrayTmpFile(shape_type const & shape, std::string const & path = "")
    : ChunkedArray<N, T>(shape),
      outer_array_(outerShape(shape)),
      file_size_(),
      file_capacity_(),
      cache_first_(0), cache_last_(0),
      cache_size_(0),
      // cache_max_size_(max(outer_array_.shape())*2)
      cache_max_size_(8)
    {
        // set shape of the chunks
        typename ChunkStorage::iterator i = outer_array_.begin(), 
             end = outer_array_.end();
        size_t size = 0;
        for(; i != end; ++i)
        {
            i->reshape(min(this->default_chunk_shape_, shape - i.point()*this->default_chunk_shape_),
                       mmap_alignment);
            size += i->alloc_size();
        }
        
        std::cerr << "file size: " << size << "\n";

    #ifdef VIGRA_NO_SPARSE_FILE
        file_capacity_ = 4*prod(this->default_chunk_shape_)*sizeof(T);
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

        // resize and memory-map the file
        static const size_t bits = sizeof(LONG)*8, mask = (size_t(1) << bits) - 1;
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
        std::cerr << "final cache size: " << cache_size_ << "\n";
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
    
    virtual pointer ptr(shape_type const & point, shape_type & strides, shape_type & border, ChunkBase<N, T> ** chunk_ptr)
    {
        if(*chunk_ptr)
        {
            (*chunk_ptr)->refcount_.fetch_sub(1, threading::memory_order_seq_cst);
            *chunk_ptr = 0;
        }
        
        if(!this->isInside(point))
            return 0;

        shape_type chunkIndex(SkipInitialization);
        ChunkShape<N>::chunkIndex(point, chunkIndex);
        Chunk & chunk = outer_array_[chunkIndex];
        
        T * p = 0;
        while(true)
        {
            int rc = chunk.refcount_.load(threading::memory_order_seq_cst);
            if(rc < 0)
            {
                // cache management in progress => try again later
                // (we can effort this simplistic collision handling 
                //  method because collisions will occur very rarely)
                threading::this_thread::yield();
            }
            else if(chunk.refcount_.compare_exchange_strong(rc, rc+1, threading::memory_order_seq_cst))
            {
                // we successfully obtained a reference
                p = chunk.pointer_.load(threading::memory_order_seq_cst);
                if(p == 0)
                {
                    // apply the double-checked locking pattern
                    threading::lock_guard<threading::mutex> guard(cache_lock_);
                    p = chunk.pointer_.load(threading::memory_order_seq_cst);
                    if(p == 0)
                    {
                        ptrdiff_t offset = file_size_;
                    #ifdef VIGRA_NO_SPARSE_FILE
                        size_t chunk_size = chunk.alloc_size();
                        if(chunk.offset_ < 0 && offset + chunk_size > file_capacity_)
                        {
                            file_capacity_ *= 2;
                            lseek(file_, file_capacity_-1, SEEK_SET);
                            if(write(file_, "0", 1) == -1)
                                throw std::runtime_error("ChunkedArrayTmpFile(): unable to resize file.");
                        }
                        file_size_ += chunk_size;
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
                            rc = 0;
                            if(c->refcount_.compare_exchange_strong(rc, -1, std::memory_order_seq_cst))
                            {
                                // refcount was zero => can unmap
                                c->unmap();
                                c->cache_next_ = 0;
                                --cache_size_;
                                c->refcount_.store(0, std::memory_order_seq_cst);
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
                break;
            }
        }
        
        strides = chunk.strides_;
        border = (chunkIndex + shape_type(1)) * this->default_chunk_shape_;
        std::size_t offset = ChunkShape<N>::offset(point, strides);
        *chunk_ptr = &chunk;
        return p + offset;
    }
    
    ChunkStorage outer_array_;  // the array of chunks

    HANDLE file_, mappedFile_;  // the file back-end
    size_t file_size_, file_capacity_;
    
    Chunk * cache_first_, * cache_last_;  // cache of loaded chunks
    std::size_t cache_max_size_, cache_size_;
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

    typedef ChunkedArray<dimensions, T>           array_type;
    typedef ChunkBase<dimensions, T>              chunk_type;
    typedef ChunkShape<dimensions>                chunk_shape;
    typedef T                                     value_type;
    typedef value_type *                          pointer;
    typedef value_type const *                    const_pointer;
    typedef value_type &                          reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;
    
    typedef T* (*SetPointerFct)(shape_type const & p, array_type * array, shape_type & strides, shape_type & border);
    
    static SetPointerFct setPointer;
    
    static T* setPointerFct(shape_type const & p, array_type * array, shape_type & strides, shape_type & border)
    {
        return array->ptr(p, strides, border);
    }

    CoupledHandle()
    : base_type(),
      pointer_(), 
      array_()
    {}

    CoupledHandle(array_type & array, NEXT const & next)
    : base_type(next),
      pointer_(), 
      array_(&array),
      chunk_()
    {
        pointer_ = array_->ptr(point(), strides_, border_, &chunk_);
    }
    
    using base_type::point;
    using base_type::shape;

    inline void incDim(int dim) 
    {
        base_type::incDim(dim);
        // if((point()[dim] & chunk_shape::mask) == 0)
        // {
            // if(point()[dim] < shape()[dim])
                // pointer_ = (*setPointer)(point(), array_, strides_, border_);
        // }
        // else
        // {
            // pointer_ += strides_[dim];
        // }
        pointer_ += strides_[dim];
        if(point()[dim] == border_[dim])
        {
            if(point()[dim] < shape()[dim])
                // pointer_ = (*setPointer)(point(), array_, strides_, border_);
                pointer_ = array_->ptr(point(), strides_, border_, &chunk_);
        }
    }

    // inline void decDim(int dim) 
    // {
        // view_.unsafePtr() -= strides_[dim];
        // base_type::decDim(dim);
    // }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        base_type::addDim(dim, d);
        if(point()[dim] < shape()[dim])
            pointer_ = array_->ptr(point(), strides_, border_, &chunk_);
        // if(dim == 0)
            // std::cerr << point() << " " << pointer_ << " " << strides_ << " " << border_ << "\n";
    }

    // inline void add(shape_type const & d) 
    // {
        // view_.unsafePtr() += dot(d, strides_);
        // base_type::add(d);
    // }
    
    // template<int DIMENSION>
    // inline void increment() 
    // {
        // view_.unsafePtr() += strides_[DIMENSION];
        // base_type::template increment<DIMENSION>();
    // }
    
    // template<int DIMENSION>
    // inline void decrement() 
    // {
        // view_.unsafePtr() -= strides_[DIMENSION];
        // base_type::template decrement<DIMENSION>();
    // }
    
    // // TODO: test if making the above a default case of the this hurts performance
    // template<int DIMENSION>
    // inline void increment(MultiArrayIndex offset) 
    // {
        // view_.unsafePtr() += offset*strides_[DIMENSION];
        // base_type::template increment<DIMENSION>(offset);
    // }
    
    // template<int DIMENSION>
    // inline void decrement(MultiArrayIndex offset) 
    // {
        // view_.unsafePtr() -= offset*strides_[DIMENSION];
        // base_type::template decrement<DIMENSION>(offset);
    // }
    
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
    shape_type strides_, border_;
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

    // typedef ChunkedArray<dimensions, T>           array_type;
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


#if 0
template <class T, class NEXT>
class CoupledHandle
: public NEXT
{
public:
    typedef NEXT                            base_type;
    typedef CoupledHandle<T, NEXT>          self_type;
    
    static const int index =                NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =  NEXT::dimensions;

    typedef T                               value_type;
    typedef T *                             pointer;
    typedef T const *                       const_pointer;
    typedef T &                             reference;
    typedef T const &                       const_reference;
    typedef typename base_type::shape_type  shape_type;

    CoupledHandle()
    : base_type(),
      pointer_(), 
      strides_()
    {}

    CoupledHandle(const_pointer p, shape_type const & strides, NEXT const & next)
    : base_type(next),
      pointer_(const_cast<pointer>(p)), 
      strides_(strides)
    {}

    template <class Stride>
    CoupledHandle(MultiArrayView<dimensions, T, Stride> const & v, NEXT const & next)
    : base_type(next),
      pointer_(const_cast<pointer>(v.data())), 
      strides_(v.stride())
    {
        vigra_precondition(v.shape() == this->shape(), "createCoupledIterator(): shape mismatch.");
    }

    inline void incDim(int dim) 
    {
        pointer_ += strides_[dim];
        base_type::incDim(dim);
    }

    inline void decDim(int dim) 
    {
        pointer_ -= strides_[dim];
        base_type::decDim(dim);
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        pointer_ += d*strides_[dim];
        base_type::addDim(dim, d);
    }

    inline void add(shape_type const & d) 
    {
        pointer_ += dot(d, strides_);
        base_type::add(d);
    }
    
    template<int DIMENSION>
    inline void increment() 
    {
        pointer_ += strides_[DIMENSION];
        base_type::template increment<DIMENSION>();
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        pointer_ -= strides_[DIMENSION];
        base_type::template decrement<DIMENSION>();
    }
    
    // TODO: test if making the above a default case of the this hurts performance
    template<int DIMENSION>
    inline void increment(MultiArrayIndex offset) 
    {
        pointer_ += offset*strides_[DIMENSION];
        base_type::template increment<DIMENSION>(offset);
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex offset) 
    {
        pointer_ -= offset*strides_[DIMENSION];
        base_type::template decrement<DIMENSION>(offset);
    }
    
    void restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        pointer_ += dot(start, strides_);
        base_type::restrictToSubarray(start, end);
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

    shape_type const & strides() const
    {
        return strides_;
    }

    pointer pointer_;
    shape_type strides_;
};

template <unsigned int N, class T, class TAIL>
struct ComposeCoupledHandle<N, TypeList<T, TAIL> >
{
    typedef typename ComposeCoupledHandle<N, TAIL>::type  BaseType;
    typedef typename MultiArrayShape<N>::type             shape_type;
    typedef CoupledHandle<T, BaseType>                    type;
    
    template <class S>
    type exec(MultiArrayView<N, T, S> const & m, 
              shape_type const & start, shape_type const & end,
              BaseType const & base)
    {
        return type(m.subarray(start, end).data(), m.stride(), base);
    }
    
    template <class S>
    type exec(MultiArrayView<N, T, S> const & m, BaseType const & base)
    {
        return type(m.data(), m.stride(), base);
    }
};

template <unsigned int N>
struct ComposeCoupledHandle<N, void>
{
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef CoupledHandle<shape_type, void>    type;
    
    type exec(shape_type const & shape)
    {
        return type(shape);
    }
    
    type exec(shape_type const & start, shape_type const & end)
    {
        return type(end-start);
    }
};

template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct CoupledHandleType
{
    // reverse the order to get the desired index order
    typedef typename MakeTypeList<T5, T4, T3, T2, T1>::type TypeList;
    typedef typename ComposeCoupledHandle<N, TypeList>::type type;
};

template <unsigned int N, class T1, class T2, class T3, class T4, class T5>
struct CoupledHandleType<N, Multiband<T1>, T2, T3, T4, T5>
{
    // reverse the order to get the desired index order
    typedef typename MakeTypeList<T5, T4, T3, T2, Multiband<T1> >::type TypeList;
    typedef typename ComposeCoupledHandle<N-1, TypeList>::type type;
};

#endif

} // namespace vigra

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HXX */
