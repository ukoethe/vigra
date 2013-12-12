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

#ifndef VIGRA_MULTI_ARRAY_BLOCKED_HXX
#define VIGRA_MULTI_ARRAY_BLOCKED_HXX

#include "metaprogramming.hxx"
#include "multi_array.hxx"

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
struct BlockedMemory;

template <unsigned int N>
struct BlockShape;

template <>
struct BlockShape<1>
{
    static const unsigned int bits = 18;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct BlockShape<2>
{
    static const unsigned int bits = 9;
    static const unsigned int value = 1 << bits;
    static const unsigned int mask = value -1;
};

template <>
struct BlockShape<3>
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
    
    static void blockIndex(Shape3 const & p, Shape3 & b)
    {
        typedef std::size_t UI;
        b[0] = (UI)p[0] >> bits;
        b[1] = (UI)p[1] >> bits;
        b[2] = (UI)p[2] >> bits;
    }
    
    static std::size_t blockOffset(Shape3 const & p, Shape3 const & s)
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
        shape[k] = (shape[k] + BlockShape<N>::mask) >> BlockShape<N>::bits;
    return shape;
}

#if 0
template <unsigned int N, class T>
class BlockedArrayChunk
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
#ifndef _WIN32
    typedef int HANDLE;
#endif    
    
    BlockedArrayChunk()
    : pointer_()
    {}
    
    ~BlockedArrayChunk()
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
                winErrorToException("BlockedArrayChunk::map(): ");
#else
            pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                     file_, offset_);
            if(!pointer_)
                throw std::runtime_error("BlockedArrayChunk::map(): mmap() failed.");
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
    size_t BlockedArrayChunk<N, T>::alloc_mask_ = winClusterSize() - 1;
#else
    template <unsigned int N, class T>
    size_t BlockedArrayChunk<N, T>::alloc_mask_ = sysconf(_SC_PAGE_SIZE) - 1;
#endif

#endif

// template <unsigned int N, class T>
// class BlockedArray
// {
  // public:
    // typedef MultiArray<N, MultiArray<N, T> > BlockStorage;
    // typedef typename BlockStorage::difference_type  shape_type;
    // typedef T value_type;
    // typedef value_type * pointer;
    // typedef value_type & reference;
    
    // BlockedArray(shape_type const & shape)
    // : shape_(shape),
      // default_block_shape_(BlockShape<N>::value),
      // outer_array_(outerShape(shape))
    // {
        // auto i = outer_array_.begin(), 
             // end = outer_array_.end();
        // for(; i != end; ++i)
            // i->reshape(min(default_block_shape_, shape - i.point()));
    // }
    
    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    // pointer ptr(shape_type const & p, shape_type & strides, shape_type & border)
    // {
        // // if(!isInside(p))
            // // return 0;
        
        // shape_type blockIndex(SkipInitialization);
        // BlockShape<N>::blockIndex(p, blockIndex);
        // MultiArray<N, T> & block = outer_array_[blockIndex];
        // strides = block.stride();
        // border = (blockIndex + shape_type(1)) * default_block_shape_;
        // std::size_t offset = BlockShape<N>::offset(p, strides);
        // return block.data() + offset;
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
    
    // shape_type shape_, default_block_shape_;
    // BlockStorage outer_array_;
// };


/*
The present implementation uses a memory-mapped sparse file to store the blocks.
A sparse file is created on Linux using the O_TRUNC flag (possibly, it is even 
the default file behavior on Linux => check), and on Windows by
calling DeviceIoControl(file_handle, FSCTL_SET_SPARSE,...) after file creation.

We can automatically delete the file upon closing. On Windows, this happens
if the file was opened with FILE_FLAG_DELETE_ON_CLOSE (it may be useful to
combine this with the flag FILE_ATTRIBUTE_TEMPORARY, which tells the OS to
avoid writing the file to disk if possible). On Linux, you can call
fileno(tmpfile()) for the same purpose.

Alternatives are:
* Keep the array in memory, but compress unused blocks.
* Don't create a file explicitly, but use the swap file instead. This is 
  achieved on Linux by mmap(..., MAP_PRIVATE | MAP_ANONYMOUS, -1, ...), 
  on Windows by calling CreateFileMapping(INVALID_HANDLE_VALUE, ...).
* Back blocks by HDF5 chunks, possibly using on-the-fly compression. This
  is in particular useful for existing HDF5 files.
* Back blocks by HDF5 datasets. This can be combined with compression 
  (both explicit and on-the-fly).
*/
template <unsigned int N, class T>
class BlockedArray
{
  public:
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
        
    BlockedArray(shape_type const & shape)
    : shape_(shape),
      default_block_shape_(BlockShape<N>::value)
    {}
    
    virtual ~BlockedArray()
    {}
    
    // reference operator[](shape_type const & p)
    // {
        // return *ptr(p);
    // }
    
    virtual pointer ptr(shape_type const & p, shape_type & strides, shape_type & border) = 0;
    
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
    
    shape_type shape_, default_block_shape_;
};

template <unsigned int N, class T>
class BlockedArrayFile
: public BlockedArray<N, T>
{
  public:
#ifndef _WIN32
    typedef int HANDLE;
#endif    
    
    class Block
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Block()
        : pointer_()
        {}
        
        ~Block()
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
                    winErrorToException("BlockedArrayChunk::map(): ");
    #else
                pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                         file_, offset_);
                if(!pointer_)
                    throw std::runtime_error("BlockedArrayChunk::map(): mmap() failed.");
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

    typedef MultiArray<N, Block> BlockStorage;
    typedef typename BlockStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    BlockedArrayFile(shape_type const & shape, char const * name)
    : BlockedArray<N, T>(shape),
      outer_array_(outerShape(shape))
    {
        // set shape of the blocks
        typename BlockStorage::iterator i = outer_array_.begin(), 
             end = outer_array_.end();
        size_t size = 0;
        for(; i != end; ++i)
        {
            i->reshape(min(this->default_block_shape_, shape - i.point()*this->default_block_shape_),
                       mmap_alignment);
            size += i->alloc_size();
        }
        
        std::cerr << "file size: " << size << "\n";

#ifdef _WIN32
        // create file or open it when it already exists
        file_ = ::CreateFile(name, GENERIC_READ | GENERIC_WRITE,
                             0, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (file_ == INVALID_HANDLE_VALUE) 
            winErrorToException("BlockedArrayFile(): ");
        
        bool isNewFile = ::GetLastError() != ERROR_ALREADY_EXISTS;
        if(isNewFile)
        {
            std::cerr << "creating new file\n";
            
            // make it a sparse file
            DWORD dwTemp;
            if(!::DeviceIoControl(file_, FSCTL_SET_SPARSE, NULL, 0, NULL, 0, &dwTemp, NULL))
                winErrorToException("BlockedArrayFile(): ");

            // set the file size
            static const size_t bits = sizeof(LONG)*8,
                                mask = (size_t(1) << bits) - 1;
            LONG fileSizeHigh = size >> bits;
            if(::SetFilePointer(file_, size & mask, &fileSizeHigh, FILE_BEGIN) == INVALID_SET_FILE_POINTER)
                winErrorToException("BlockedArrayFile(): ");
            if(!::SetEndOfFile(file_)) 
                winErrorToException("BlockedArrayFile(): ");
        }
        
        // memory-map the file
        mappedFile_ = CreateFileMapping(file_, NULL, PAGE_READWRITE, 0, 0, NULL);
        if(!mappedFile_)
            winErrorToException("BlockedArrayFile(): ");
#else
        mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
        mappedFile_ = file_ = open(name, O_RDWR | O_CREAT | O_TRUNC, mode);
        if(file_ == -1)
            throw std::runtime_error("BlockedArrayFile(): unable to open file.");
        lseek(file_, size-1, SEEK_SET);
        if(write(file_, "0", 1) == -1)
            throw std::runtime_error("BlockedArrayFile(): unable to resize file.");
#endif
        
        // tell the blocks about the memory mapping
        i = outer_array_.begin();
        int offset = 0;
        for(; i != end; ++i)
        {
            i->setFile(mappedFile_, offset);
            offset += i->alloc_size();
        }
    }
    
    ~BlockedArrayFile()
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
        typename BlockStorage::iterator  i = outer_array_.begin(), 
             end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->map();
        }
    }
    
    void unmap()
    {
        typename BlockStorage::iterator  i = outer_array_.begin(), 
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
        
        shape_type blockIndex(SkipInitialization);
        BlockShape<N>::blockIndex(p, blockIndex);
        Block & block = outer_array_[blockIndex];
        block.map();
        strides = block.strides_;
        border = (blockIndex + shape_type(1)) * this->default_block_shape_;
        std::size_t offset = BlockShape<N>::offset(p, strides);
        return block.pointer_ + offset;
    }
    
    HANDLE file_, mappedFile_;
    BlockStorage outer_array_;
};


template <unsigned int N, class T>
class BlockedArrayTmpFile
: public BlockedArray<N, T>
{
  public:
#ifndef _WIN32
    typedef int HANDLE;
#endif    
    
    class Block
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;
        
        Block()
        : pointer_()
        {}
        
        ~Block()
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
        
        void map(size_t offset)
        {
            if(!pointer_)
            {
            #ifdef VIGRA_NO_SPARSE_FILE
                offset_ = offset;
            #endif
            #ifdef _WIN32
                static const size_t bits = sizeof(DWORD)*8,
                                    mask = (size_t(1) << bits) - 1;
                pointer_ = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                                  offset_ >> bits, offset_ & mask, alloc_size_);
                if(!pointer_)
                    winErrorToException("BlockedArrayChunk::map(): ");
            #else
                pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                         file_, offset_);
                if(!pointer_)
                    throw std::runtime_error("BlockedArrayChunk::map(): mmap() failed.");
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

    typedef MultiArray<N, Block> BlockStorage;
    typedef typename BlockStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
    BlockedArrayTmpFile(shape_type const & shape, std::string const & path = "")
    : BlockedArray<N, T>(shape),
      outer_array_(outerShape(shape)),
      file_size_(),
      file_capacity_()
    {
        // set shape of the blocks
        typename BlockStorage::iterator i = outer_array_.begin(), 
             end = outer_array_.end();
        size_t size = 0;
        for(; i != end; ++i)
        {
            i->reshape(min(this->default_block_shape_, shape - i.point()*this->default_block_shape_),
                       mmap_alignment);
            size += i->alloc_size();
        }
        
        std::cerr << "file size: " << size << "\n";

#ifdef _WIN32
        file_capacity_ = size;
        
        // create a temp file
        file_ = ::CreateFile(winTempFileName(path).c_str(), GENERIC_READ | GENERIC_WRITE,
                             0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE, NULL);
        if (file_ == INVALID_HANDLE_VALUE) 
            winErrorToException("BlockedArrayTmpFile(): ");
            
        // make it a sparse file
        DWORD dwTemp;
        if(!::DeviceIoControl(file_, FSCTL_SET_SPARSE, NULL, 0, NULL, 0, &dwTemp, NULL))
            winErrorToException("BlockedArrayTmpFile(): ");

        // resize and memory-map the file
        static const size_t bits = sizeof(LONG)*8, mask = (size_t(1) << bits) - 1;
        mappedFile_ = CreateFileMapping(file_, NULL, PAGE_READWRITE, 
                                        file_capacity_ >> bits, file_capacity_ & mask, NULL);
        if(!mappedFile_)
            winErrorToException("BlockedArrayTmpFile(): ");
#else
        mappedFile_ = file_ = fileno(tmpfile());
        if(file_ == -1)
            throw std::runtime_error("BlockedArrayTmpFile(): unable to open file.");
    #ifdef VIGRA_NO_SPARSE_FILE
        file_capacity_ = 4*prod(this->default_block_shape_)*sizeof(T);
    #else
        file_capacity_ = size;
    #endif
        lseek(file_, file_capacity_-1, SEEK_SET);
        if(write(file_, "0", 1) == -1)
            throw std::runtime_error("BlockedArrayTmpFile(): unable to resize file.");
#endif
        
        // tell the blocks about the memory mapping
        i = outer_array_.begin();
        int offset = 0;
        for(; i != end; ++i)
        {
            i->setFile(mappedFile_, offset);
            offset += i->alloc_size();
        }
    }
    
    ~BlockedArrayTmpFile()
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
        typename BlockStorage::iterator  i = outer_array_.begin(), 
             end = outer_array_.end();
        for(; i != end; ++i)
        {
            i->map();
        }
    }
    
    void unmap()
    {
        typename BlockStorage::iterator  i = outer_array_.begin(), 
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
        
        shape_type blockIndex(SkipInitialization);
        BlockShape<N>::blockIndex(p, blockIndex);
        Block & block = outer_array_[blockIndex];
        if(block.pointer_ == 0)
        {
            size_t block_size = block.alloc_size();
        #ifdef VIGRA_NO_SPARSE_FILE
            if(file_size_ + block_size > file_capacity_)
            {
                file_capacity_ *= 2;
                lseek(file_, file_capacity_-1, SEEK_SET);
                if(write(file_, "0", 1) == -1)
                    throw std::runtime_error("BlockedArrayTmpFile(): unable to resize file.");
            }
        #endif
            block.map(file_size_);
            file_size_ += block_size;
        }
        strides = block.strides_;
        border = (blockIndex + shape_type(1)) * this->default_block_shape_;
        std::size_t offset = BlockShape<N>::offset(p, strides);
        return block.pointer_ + offset;
    }
    
    HANDLE file_, mappedFile_;
    BlockStorage outer_array_;
    size_t file_size_, file_capacity_;
};



template <class T, class NEXT>
class CoupledHandle<BlockedMemory<T>, NEXT>
: public NEXT
{
public:
    typedef NEXT                                  base_type;
    typedef CoupledHandle<BlockedMemory<T>, NEXT> self_type;
    
    static const int index =                      NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef BlockedArray<dimensions, T>           array_type;
    typedef BlockShape<dimensions>                block_shape;
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
      pointer_((*setPointer)(shape_type(), &array, strides_, border_)), 
      array_(&array)
    {}
    
    using base_type::point;
    using base_type::shape;

    inline void incDim(int dim) 
    {
        base_type::incDim(dim);
        // if((point()[dim] & block_shape::mask) == 0)
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
                pointer_ = array_->ptr(point(), strides_, border_);
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
            pointer_ = (*setPointer)(point(), array_, strides_, border_);
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
};

template <class T, class NEXT>
typename CoupledHandle<BlockedMemory<T>, NEXT>::SetPointerFct 
CoupledHandle<BlockedMemory<T>, NEXT>::setPointer = &CoupledHandle<BlockedMemory<T>, NEXT>::setPointerFct;


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

#endif /* VIGRA_MULTI_ARRAY_BLOCKED_HXX */
