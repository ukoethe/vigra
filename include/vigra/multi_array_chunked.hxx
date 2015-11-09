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
    (These numbers refer to nested loop iteration. Scan-order iteration
     is unfortunately 3.5 times slower on the Mac. On the other hand,
     two-level indexing as faster on a Mac than on Linux and Windows --
     the speed penalty is only a factor of 2 rather than 3.)

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

#include <queue>
#include <string>

#include "multi_fwd.hxx"
#include "multi_handle.hxx"
#include "multi_array.hxx"
#include "memory.hxx"
#include "metaprogramming.hxx"
#include "threading.hxx"
#include "compression.hxx"

// // FIXME: why is this needed when compiling the Python bindng,
// //        but not when compiling test_multiarray_chunked?
// #if defined(__GNUC__)
// #  define memory_order_release memory_order_seq_cst
// #  define memory_order_acquire memory_order_seq_cst
// #endif

#ifdef _WIN32
# include "windows.h"
#else
# include <fcntl.h>
# include <stdlib.h>
# include <unistd.h>
# include <sys/stat.h>
# include <sys/mman.h>
# include <cstdio>
#endif

// Bounds checking Macro used if VIGRA_CHECK_BOUNDS is defined.
#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(this->isInside(diff), "Index out of bounds")
#else
#define VIGRA_ASSERT_INSIDE(diff)
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

template <unsigned int N, class T>
class IteratorChunkHandle;

namespace detail {

template <unsigned int N>
struct ChunkIndexing
{
    template <class T, int M>
    static void chunkIndex(TinyVector<T, M> const & p,
                           TinyVector<T, M> const & bits,
                           TinyVector<T, M> & index)
    {
        typedef std::size_t UI;
        ChunkIndexing<N-1>::chunkIndex(p, bits, index);
        index[N-1] = (UI)p[N-1] >> bits[N-1];
    }

    template <class T, int M>
    static std::size_t chunkOffset(TinyVector<T, M> const & p,
                                   TinyVector<T, M> const & bits,
                                   TinyVector<T, M> const & strides)
    {
        typedef std::size_t UI;
        return ChunkIndexing<N-1>::chunkOffset(p, bits, strides) +
               ((UI)p[N-1] >> bits[N-1]) * strides[N-1];
    }

    template <class T, int M>
    static std::size_t offsetInChunk(TinyVector<T, M> const & p,
                                     TinyVector<T, M> const & mask,
                                     TinyVector<T, M> const & strides)
    {
        typedef std::size_t UI;
        return ChunkIndexing<N-1>::offsetInChunk(p, mask, strides) +
               ((UI)p[N-1] & (UI)mask[N-1]) * strides[N-1];
    }
};

template <>
struct ChunkIndexing<1>
{
    template <class T, int M>
    static void chunkIndex(TinyVector<T, M> const & p,
                           TinyVector<T, M> const & bits,
                           TinyVector<T, M> & index)
    {
        typedef std::size_t UI;
        index[0] = (UI)p[0] >> bits[0];
    }

    template <class T, int M>
    static std::size_t chunkOffset(TinyVector<T, M> const & p,
                                   TinyVector<T, M> const & bits,
                                   TinyVector<T, M> const & strides)
    {
        typedef std::size_t UI;
        return ((UI)p[0] >> bits[0]) * strides[0];
    }

    template <class T, int M>
    static std::size_t offsetInChunk(TinyVector<T, M> const & p,
                                     TinyVector<T, M> const & mask,
                                     TinyVector<T, M> const & strides)
    {
        typedef std::size_t UI;
        return ((UI)p[0] & (UI)mask[0]) * strides[0];
    }
};

template <class T, int M>
inline TinyVector<T, M>
computeChunkArrayShape(TinyVector<T, M> shape,
                       TinyVector<T, M> const & bits,
                       TinyVector<T, M> const & mask)
{
    for(int k=0; k<M; ++k)
        shape[k] = (shape[k] + mask[k]) >> bits[k];
    return shape;
}

template <class T, int M>
inline T
defaultCacheSize(TinyVector<T, M> const & shape)
{
    T res = max(shape);
    for(int k=0; k<M-1; ++k)
        for(int j=k+1; j<M; ++j)
            res = std::max(res, shape[k]*shape[j]);
    return res + 1;
}

} // namespace detail

template <unsigned int N, class T>
class ChunkBase
{
  public:
    typedef typename MultiArrayShape<N>::type shape_type;
    typedef T value_type;
    typedef T* pointer;

    ChunkBase()
    : strides_()
    , pointer_()
    {}

    ChunkBase(shape_type const & strides, pointer p = 0)
    : strides_(strides)
    , pointer_(p)
    {}

    typename MultiArrayShape<N>::type strides_;
    T * pointer_;
};

template <unsigned int N, class T>
class SharedChunkHandle
{
  public:
    typedef typename MultiArrayShape<N>::type shape_type;

    static const long chunk_asleep = -2;
    static const long chunk_uninitialized = -3;
    static const long chunk_locked = -4;
    static const long chunk_failed = -5;

    SharedChunkHandle()
    : pointer_(0)
    , chunk_state_()
    {
        chunk_state_ = chunk_uninitialized;
    }

    SharedChunkHandle(SharedChunkHandle const & rhs)
    : pointer_(rhs.pointer_)
    , chunk_state_()
    {
        chunk_state_ = chunk_uninitialized;
    }

    shape_type const & strides() const
    {
        return pointer_->strides_;
    }

    ChunkBase<N, T> * pointer_;
    mutable threading::atomic_long chunk_state_;

  private:
    SharedChunkHandle & operator=(SharedChunkHandle const & rhs);
};

template <unsigned int N, class T>
class ChunkedArrayBase
{
  public:
    enum ActualDimension{ actual_dimension = (N == 0) ? 1 : N };
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    typedef ChunkBase<N, T> Chunk;

    ChunkedArrayBase()
    : shape_()
    , chunk_shape_()
    {}

    ChunkedArrayBase(shape_type const & shape, shape_type const & chunk_shape)
    : shape_(shape)
    , chunk_shape_(prod(chunk_shape) > 0 ? chunk_shape : detail::ChunkShape<N, T>::defaultShape())
    {}

    virtual ~ChunkedArrayBase()
    {}

    virtual void unrefChunk(IteratorChunkHandle<N, T> * h) const = 0;

    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h) = 0;

    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h) const = 0;

    virtual std::string backend() const = 0;

    virtual shape_type chunkArrayShape() const = 0;

    virtual bool isReadOnly() const
    {
        return false;
    }

    MultiArrayIndex size() const
    {
        return prod(shape_);
    }

    shape_type const & shape() const
    {
        return shape_;
    }

    MultiArrayIndex shape(MultiArrayIndex d) const
    {
        return shape_[d];
    }

    shape_type const & chunkShape() const
    {
        return chunk_shape_;
    }

    MultiArrayIndex chunkShape(MultiArrayIndex d) const
    {
        return chunk_shape_[d];
    }

    bool isInside(shape_type const & p) const
    {
        for(int d=0; d<N; ++d)
            if(p[d] < 0 || p[d] >= shape_[d])
                return false;
        return true;
    }

    shape_type shape_, chunk_shape_;
};

template <unsigned int N, class T>
class ChunkedArray;

struct ChunkUnrefProxyBase
{
    virtual ~ChunkUnrefProxyBase() {}
};

template <unsigned int N, class T_MaybeConst>
class MultiArrayView<N, T_MaybeConst, ChunkedArrayTag>
: public ChunkedArrayBase<N, typename UnqualifiedType<T_MaybeConst>::type>
{
  public:
    enum ActualDimension { actual_dimension = (N==0) ? 1 : N };
    typedef typename UnqualifiedType<T_MaybeConst>::type     T;
    typedef T value_type;   // FIXME: allow Multiband<T> ???
    typedef T_MaybeConst & reference;
    typedef const value_type &const_reference;
    typedef T_MaybeConst * pointer;
    typedef const value_type *const_pointer;
    typedef typename MultiArrayShape<actual_dimension>::type difference_type;
    typedef difference_type key_type;
    typedef difference_type size_type;
    typedef difference_type shape_type;
    typedef MultiArrayIndex difference_type_1;
    typedef ChunkIterator<actual_dimension, T_MaybeConst>         chunk_iterator;
    typedef ChunkIterator<actual_dimension, T const>   chunk_const_iterator;
    typedef StridedScanOrderIterator<actual_dimension, ChunkedMemory<T_MaybeConst>, T_MaybeConst&, T_MaybeConst*> iterator;
    typedef StridedScanOrderIterator<actual_dimension, ChunkedMemory<T const>, T const &, T const *> const_iterator;
    typedef MultiArrayView<N, T_MaybeConst, ChunkedArrayTag> view_type;
    typedef MultiArrayView<N, T const, ChunkedArrayTag> const_view_type;
    typedef ChunkedArrayTag StrideTag;
    typedef ChunkBase<N, T> Chunk;

    typedef MultiArray<N, Chunk> ChunkHolder;

    struct UnrefProxy
    : public ChunkUnrefProxyBase
    {
        UnrefProxy(int size, ChunkedArray<N, T> * array)
        : chunks_(size)
        , array_(array)
        {}

        ~UnrefProxy()
        {
            if(array_)
                array_->unrefChunks(chunks_);
        }

        ArrayVector<SharedChunkHandle<N, T> *> chunks_;
        ChunkedArray<N, T> * array_;
    };

    virtual shape_type chunkArrayShape() const
    {
        return chunks_.shape();
    }

    shape_type chunkStart(shape_type const & global_start) const
    {
        shape_type chunk_start(SkipInitialization);
        detail::ChunkIndexing<N>::chunkIndex(global_start, bits_, chunk_start);
        return chunk_start;
    }

    shape_type chunkStop(shape_type global_stop) const
    {
        global_stop -= shape_type(1);
        shape_type chunk_stop(SkipInitialization);
        detail::ChunkIndexing<N>::chunkIndex(global_stop, bits_, chunk_stop);
        chunk_stop += shape_type(1);
        return chunk_stop;
    }

    virtual void unrefChunk(IteratorChunkHandle<N, T> *) const {}

    virtual T* chunkForIterator(shape_type const & point,
                                shape_type & strides, shape_type & upper_bound,
                                IteratorChunkHandle<N, T> * h)
    {
        return const_cast<MultiArrayView const *>(this)->chunkForIterator(point, strides, upper_bound, h);
    }

    virtual T* chunkForIterator(shape_type const & point,
                                shape_type & strides, shape_type & upper_bound,
                                IteratorChunkHandle<N, T> * h) const
    {
        shape_type global_point = point + h->offset_;

        if(!this->isInside(global_point))
        {
            upper_bound = point + this->chunk_shape_;
            return 0;
        }

        global_point += offset_;
        shape_type coffset = offset_ + h->offset_;

        shape_type chunkIndex = chunkStart(global_point);
        Chunk const * chunk = &chunks_[chunkIndex];
        strides = chunk->strides_;
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_ - coffset;
        std::size_t offset = detail::ChunkIndexing<N>::offsetInChunk(global_point, mask_, strides);
        return const_cast<T*>(chunk->pointer_ + offset);
    }

    virtual std::string backend() const
    {
        return "MultiArrayView<ChunkedArrayTag>";
    }

    MultiArrayView()
    : ChunkedArrayBase<N, T>()
    {}

    MultiArrayView(shape_type const & shape, shape_type const & chunk_shape)
    : ChunkedArrayBase<N, T>(shape, chunk_shape)
    {}

    MultiArrayView & operator=(MultiArrayView const & rhs)
    {
        if(this != &rhs)
        {
            if(!hasData())
            {
                ChunkedArrayBase<N, T>::operator=(rhs);
                chunks_ = rhs.chunks_;
                offset_ = rhs.offset_;
                bits_ = rhs.bits_;
                mask_ = rhs.mask_;
                unref_ = rhs.unref_;
            }
            else
            {
                vigra_precondition(this->shape() == rhs.shape(),
                                   "MultiArrayView::operator=(): shape mismatch.");
                iterator i = begin(), ie = end();
                const_iterator j = rhs.begin();
                for(; i != ie; ++i, ++j)
                    *i = *j;
            }
        }
        return *this;
    }

    #define VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(op) \
    template<class U, class C1> \
    MultiArrayView & operator op(MultiArrayView<N, U, C1> const & rhs) \
    { \
        vigra_precondition(this->shape() == rhs.shape(), \
                           "MultiArrayView::operator" #op "(): shape mismatch."); \
        iterator i = begin(), ie = end(); \
        typename MultiArrayView<N, U, C1>::const_iterator j = rhs.begin(); \
        for(; i != ie; ++i, ++j) \
            *i op detail::RequiresExplicitCast<value_type>::cast(*j); \
        return *this; \
    } \
     \
    MultiArrayView & operator op(value_type const & v) \
    { \
        if(hasData()) \
        { \
            iterator i = begin(), ie = end(); \
            for(; i != ie; ++i) \
                *i op v; \
        } \
        return *this; \
    }

    VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(=)
    VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(+=)
    VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(-=)
    VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(*=)
    VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN(/=)

    #undef VIGRA_CHUNKED_ARRAY_VIEW_ASSIGN

    // template<class Expression>
    // MultiArrayView & operator=(multi_math::MultiMathOperand<Expression> const & rhs)
    // {
        // multi_math::math_detail::assign(*this, rhs);
        // return *this;
    // }

        // /** Add-assignment of an array expression. Fails with
            // <tt>PreconditionViolation</tt> exception when the shapes do not match.
         // */
    // template<class Expression>
    // MultiArrayView & operator+=(multi_math::MultiMathOperand<Expression> const & rhs)
    // {
        // multi_math::math_detail::plusAssign(*this, rhs);
        // return *this;
    // }

        // /** Subtract-assignment of an array expression. Fails with
            // <tt>PreconditionViolation</tt> exception when the shapes do not match.
         // */
    // template<class Expression>
    // MultiArrayView & operator-=(multi_math::MultiMathOperand<Expression> const & rhs)
    // {
        // multi_math::math_detail::minusAssign(*this, rhs);
        // return *this;
    // }

        // /** Multiply-assignment of an array expression. Fails with
            // <tt>PreconditionViolation</tt> exception when the shapes do not match.
         // */
    // template<class Expression>
    // MultiArrayView & operator*=(multi_math::MultiMathOperand<Expression> const & rhs)
    // {
        // multi_math::math_detail::multiplyAssign(*this, rhs);
        // return *this;
    // }

        // /** Divide-assignment of an array expression. Fails with
            // <tt>PreconditionViolation</tt> exception when the shapes do not match.
         // */
    // template<class Expression>
    // MultiArrayView & operator/=(multi_math::MultiMathOperand<Expression> const & rhs)
    // {
        // multi_math::math_detail::divideAssign(*this, rhs);
        // return *this;
    // }

    reference operator[](shape_type point)
    {
        VIGRA_ASSERT_INSIDE(point);
        point += offset_;
        Chunk * chunk = chunks_.data() +
                        detail::ChunkIndexing<N>::chunkOffset(point, bits_, chunks_.stride());
        return *(chunk->pointer_ +
                 detail::ChunkIndexing<N>::offsetInChunk(point, mask_, chunk->strides_));
    }

    const_reference operator[](shape_type const & point) const
    {
        return const_cast<MultiArrayView *>(this)->operator[](point);
    }

    template <int M>
    MultiArrayView <N-M, T, ChunkedArrayTag>
    operator[](const TinyVector<MultiArrayIndex, M> &d) const
    {
        return bindInner(d);
    }

    reference operator[](difference_type_1 d)
    {
        return operator[](scanOrderIndexToCoordinate(d));
    }

    const_reference operator[](difference_type_1 d) const
    {
        return operator[](scanOrderIndexToCoordinate(d));
    }

    difference_type scanOrderIndexToCoordinate(difference_type_1 d) const
    {
        difference_type coord(SkipInitialization);
        detail::ScanOrderToCoordinate<actual_dimension>::exec(d, this->shape_, coord);
        return coord;
    }

        /** convert coordinate to scan-order index.
         */
    difference_type_1 coordinateToScanOrderIndex(const difference_type &d) const
    {
        return detail::CoordinateToScanOrder<actual_dimension>::exec(this->shape_, d);
    }

        // /** 1D array access. Use only if N == 1.
         // */
    // reference operator() (difference_type_1 x)
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x));
        // return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x)];
    // }

        // /** 2D array access. Use only if N == 2.
         // */
    // reference operator() (difference_type_1 x, difference_type_1 y)
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x, y));
        // return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x, y)];
    // }

        // /** 3D array access. Use only if N == 3.
         // */
    // reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z)
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x, y, z));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    // }

        // /** 4D array access. Use only if N == 4.
         // */
    // reference operator() (difference_type_1 x, difference_type_1 y,
                          // difference_type_1 z, difference_type_1 u)
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x, y, z, u));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    // }

        // /** 5D array access. Use only if N == 5.
         // */
    // reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z,
                          // difference_type_1 u, difference_type_1 v)
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x, y,z, u,v));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    // }

        // /** 1D const array access. Use only if N == 1.
         // */
    // const_reference operator() (difference_type_1 x) const
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x));
        // return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x)];
    // }

        // /** 2D const array access. Use only if N == 2.
         // */
    // const_reference operator() (difference_type_1 x, difference_type_1 y) const
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x, y));
        // return m_ptr [detail::CoordinatesToOffest<StrideTag>::exec(m_stride, x, y)];
    // }

        // /** 3D const array access. Use only if N == 3.
         // */
    // const_reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z) const
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x,y,z));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z];
    // }

        // /** 4D const array access. Use only if N == 4.
         // */
    // const_reference operator() (difference_type_1 x, difference_type_1 y,
                                // difference_type_1 z, difference_type_1 u) const
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x,y,z,u));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u];
    // }

        // /** 5D const array access. Use only if N == 5.
         // */
    // const_reference operator() (difference_type_1 x, difference_type_1 y, difference_type_1 z,
                                // difference_type_1 u, difference_type_1 v) const
    // {
        // VIGRA_ASSERT_INSIDE(difference_type(x,y,z,u,v));
        // return m_ptr [m_stride[0]*x + m_stride[1]*y + m_stride[2]*z + m_stride[3]*u + m_stride[4]*v];
    // }

    template <class U>
    MultiArrayView & init(const U & init)
    {
        return operator=(init);
    }

    template <class U, class CN>
    void copy(const MultiArrayView <N, U, CN>& rhs)
    {
        operator=(rhs);
    }

    template <class T2, class C2>
    void swapData(MultiArrayView <N, T2, C2> rhs)
    {
        if(this == &rhs)
            return;
        vigra_precondition(this->shape() == rhs.shape(),
                           "MultiArrayView::swapData(): shape mismatch.");
        iterator i = begin(), ie = end();
        typename MultiArrayView<N, T2, C2>::iterator j = rhs.begin();
        for(; i != ie; ++i, ++j)
            std::swap(*i, *j);
    }

    bool isUnstrided(unsigned int dimension = N-1) const
    {
        if(chunks_.size() > 1)
            return false;
        difference_type s = vigra::detail::defaultStride<actual_dimension>(this->shape());
        for(unsigned int k = 0; k <= dimension; ++k)
            if(chunks_.data()->strides_[k] != s[k])
                return false;
        return true;
    }

    MultiArrayView<N-1, value_type, ChunkedArrayTag>
    bindAt(MultiArrayIndex m, MultiArrayIndex d) const
    {
        MultiArrayView<N-1, value_type, ChunkedArrayTag> res(this->shape_.dropIndex(m), this->chunk_shape_.dropIndex(m));
        res.offset_ = offset_.dropIndex(m);
        res.bits_   = bits_.dropIndex(m);
        res.mask_   = mask_.dropIndex(m);
        res.chunks_.reshape(chunks_.shape().dropIndex(m));
        res.unref_ = unref_;

        typedef std::size_t UI;
        UI start = offset_[m] + d;
        UI chunk_start = start >> bits_[m];
        UI startInChunk = start - chunk_start * this->chunk_shape_[m];

        MultiArrayView<N-1, Chunk> view(chunks_.bindAt(m, chunk_start));
        MultiCoordinateIterator<N-1> i(view.shape()),
                                     end(i.getEndIterator());
        for(; i != end; ++i)
        {
            res.chunks_[*i].pointer_ = view[*i].pointer_ + startInChunk*view[*i].strides_[m];
            res.chunks_[*i].strides_ = view[*i].strides_.dropIndex(m);
        }

        return res;
    }

    template <unsigned int M>
    MultiArrayView <N-1, value_type, ChunkedArrayTag>
    bind (difference_type_1 d) const
    {
        return bindAt(M, d);
    }

    MultiArrayView <N-1, value_type, ChunkedArrayTag>
    bindOuter (difference_type_1 d) const
    {
        return bindAt(N-1, d);
    }

    template <int M, class Index>
    MultiArrayView <N-M, value_type, ChunkedArrayTag>
    bindOuter(const TinyVector <Index, M> &d) const
    {
        return bindAt(N-1, d[M-1]).bindOuter(d.dropIndex(M-1));
    }

    template <class Index>
    MultiArrayView <N-1, value_type, ChunkedArrayTag>
    bindOuter(const TinyVector <Index, 1> &d) const
    {
        return bindAt(N-1, d[0]);
    }

    MultiArrayView <N-1, value_type, ChunkedArrayTag>
    bindInner (difference_type_1 d) const
    {
        return bindAt(0, d);
    }

    template <int M, class Index>
    MultiArrayView <N-M, value_type, ChunkedArrayTag>
    bindInner(const TinyVector <Index, M> &d) const
    {
        return bindAt(0, d[0]).bindInner(d.dropIndex(0));
    }

    template <class Index>
    MultiArrayView <N-1, value_type, ChunkedArrayTag>
    bindInner(const TinyVector <Index, 1> &d) const
    {
        return bindAt(0, d[0]);
    }

    // MultiArrayView <N, typename ExpandElementResult<T>::type, StridedArrayTag>
    // bindElementChannel(difference_type_1 i) const
    // {
        // vigra_precondition(0 <= i && i < ExpandElementResult<T>::size,
              // "MultiArrayView::bindElementChannel(i): 'i' out of range.");
        // return expandElements(0).bindInner(i);
    // }

    // MultiArrayView <N+1, typename ExpandElementResult<T>::type, StridedArrayTag>
    // expandElements(difference_type_1 d) const;

    // MultiArrayView <N+1, T, StrideTag>
    // insertSingletonDimension (difference_type_1 i) const;

    // MultiArrayView<N, Multiband<value_type>, StrideTag> multiband() const
    // {
        // return MultiArrayView<N, Multiband<value_type>, StrideTag>(*this);
    // }

    // MultiArrayView<1, T, StridedArrayTag> diagonal() const
    // {
        // return MultiArrayView<1, T, StridedArrayTag>(Shape1(vigra::min(m_shape)),
                                                     // Shape1(vigra::sum(m_stride)), m_ptr);
    // }

    inline void
    checkSubarrayBounds(shape_type const & start, shape_type const & stop,
                        std::string message) const
    {
        message += ": subarray out of bounds.";
        vigra_precondition(allLessEqual(shape_type(), start) &&
                           allLess(start, stop) &&
                           allLessEqual(stop, this->shape_),
                           message);
    }

    MultiArrayView<N, value_type, ChunkedArrayTag>
    subarray(shape_type start, shape_type stop)
    {
        checkSubarrayBounds(start, stop, "MultiArrayView<N-1, T, ChunkedArrayTag>::subarray()");
        start += offset_;
        stop  += offset_;
        shape_type chunk_start(chunkStart(start));

        MultiArrayView<N, value_type, ChunkedArrayTag> view(stop-start, this->chunk_shape_);
        view.chunks_ = chunks_.subarray(chunk_start, chunkStop(stop));
        view.offset_ = start - chunk_start * this->chunk_shape_;
        view.bits_   = bits_;
        view.mask_   = mask_;
        view.unref_ = unref_;
        return view;
    }

        // /** apply an additional striding to the image, thereby reducing
            // the shape of the array.
            // for example, multiplying the stride of dimension one by three
            // turns an appropriately laid out (interleaved) rgb image into
            // a single band image.
        // */
    // MultiArrayView <N, T, StridedArrayTag>
    // stridearray (const difference_type &s) const
    // {
        // difference_type shape = m_shape;
        // for (unsigned int i = 0; i < actual_dimension; ++i)
            // shape [i] /= s [i];
        // return MultiArrayView <N, T, StridedArrayTag>(shape, m_stride * s, m_ptr);
    // }

    MultiArrayView <N, value_type, ChunkedArrayTag>
    transpose () const
    {
        return transpose(difference_type::linearSequence(N-1, -1));
    }

    MultiArrayView <N, value_type, ChunkedArrayTag>
    transpose(const difference_type &permutation) const
    {
        MultiArrayView<N, value_type, ChunkedArrayTag>
            view(vigra::transpose(this->shape_, permutation), vigra::transpose(this->chunk_shape_, permutation));
        view.chunks_        = chunks_.transpose(permutation); // also checks if permutation is valid
        view.offset_        = vigra::transpose(offset_, permutation);
        view.bits_          = vigra::transpose(bits_, permutation);
        view.mask_          = vigra::transpose(mask_, permutation);
        view.unref_         = unref_;
        typename MultiArray<N, Chunk>::iterator i = view.chunks_.begin(),
                                                iend = view.chunks_.end();
        for(; i != iend; ++i)
            i->strides_ = vigra::transpose(i->strides_, permutation);
        return view;
    }

    // MultiArrayView <N, T, StridedArrayTag>
    // permuteDimensions (const difference_type &s) const;

        // /** Permute the dimensions of the array so that the strides are in ascending order.
            // Determines the appropriate permutation and then calls permuteDimensions().
        // */
    // MultiArrayView <N, T, StridedArrayTag>
    // permuteStridesAscending() const;

        // /** Permute the dimensions of the array so that the strides are in descending order.
            // Determines the appropriate permutation and then calls permuteDimensions().
        // */
    // MultiArrayView <N, T, StridedArrayTag>
    // permuteStridesDescending() const;

        // /** Compute the ordering of the strides in this array.
            // The result is describes the current permutation of the axes relative
            // to the standard ascending stride order.
        // */
    // difference_type strideOrdering() const
    // {
        // return strideOrdering(m_stride);
    // }

        // /** Compute the ordering of the given strides.
            // The result is describes the current permutation of the axes relative
            // to the standard ascending stride order.
        // */
    // static difference_type strideOrdering(difference_type strides);

    template <class U, class C1>
    bool operator==(MultiArrayView<N, U, C1> const & rhs) const
    {
        if(this->shape() != rhs.shape())
            return false;
        const_iterator i = begin(), ie = end();
        typename MultiArrayView<N, U, C1>::const_iterator j = rhs.begin();
        for(; i != ie; ++i, ++j)
            if(*i != *j)
                return false;
        return true;
    }

    template <class U, class C1>
    bool operator!=(MultiArrayView<N, U, C1> const & rhs) const
    {
        return !operator==(rhs);
    }

    // bool all() const
    // {
        // bool res = true;
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::AllTrueReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // return res;
    // }

    // bool any() const
    // {
        // bool res = false;
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::AnyTrueReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // return res;
    // }

    // void minmax(T * minimum, T * maximum) const
    // {
        // std::pair<T, T> res(NumericTraits<T>::max(), NumericTraits<T>::min());
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::MinmaxReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // *minimum = res.first;
        // *maximum = res.second;
    // }

    // template <class U>
    // void meanVariance(U * mean, U * variance) const
    // {
        // typedef typename NumericTraits<U>::RealPromote R;
        // R zero = R();
        // triple<double, R, R> res(0.0, zero, zero);
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::MeanVarianceReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // *mean     = res.second;
        // *variance = res.third / res.first;
    // }

    // template <class U>
    // U sum() const
    // {
        // U res = NumericTraits<U>::zero();
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::SumReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // return res;
    // }

    // template <class U, class S>
    // void sum(MultiArrayView<N, U, S> sums) const
    // {
        // transformMultiArray(srcMultiArrayRange(*this),
                            // destMultiArrayRange(sums),
                            // FindSum<U>());
    // }

    // template <class U>
    // U product() const
    // {
        // U res = NumericTraits<U>::one();
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::ProdReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // return res;
    // }

    // typename NormTraits<MultiArrayView>::SquaredNormType
    // squaredNorm() const
    // {
        // typedef typename NormTraits<MultiArrayView>::SquaredNormType SquaredNormType;
        // SquaredNormType res = NumericTraits<SquaredNormType>::zero();
        // detail::reduceOverMultiArray(traverser_begin(), shape(),
                                     // res,
                                     // detail::SquaredL2NormReduceFunctor(),
                                     // MetaInt<actual_dimension-1>());
        // return res;
    // }

    // typename NormTraits<MultiArrayView>::NormType
    // norm(int type = 2, bool useSquaredNorm = true) const;

    bool hasData () const
    {
        return chunks_.hasData();
    }

    iterator begin()
    {
        return createCoupledIterator(*this);
    }

    iterator end()
    {
        return begin().getEndIterator();
    }

    const_iterator cbegin() const
    {
        return createCoupledIterator(const_cast<MultiArrayView const &>(*this));
    }

    const_iterator cend() const
    {
        return cbegin().getEndIterator();
    }

    const_iterator begin() const
    {
        return createCoupledIterator(*this);
    }

    const_iterator end() const
    {
        return begin().getEndIterator();
    }

    chunk_iterator chunk_begin(shape_type const & start, shape_type const & stop)
    {
        checkSubarrayBounds(start, stop, "MultiArrayView<N-1, T, ChunkedArrayTag>::chunk_begin()");
        return chunk_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    chunk_iterator chunk_end(shape_type const & start, shape_type const & stop)
    {
        return chunk_begin(start, stop).getEndIterator();
    }

    chunk_const_iterator chunk_begin(shape_type const & start, shape_type const & stop) const
    {
        checkSubarrayBounds(start, stop, "MultiArrayView<N-1, T, ChunkedArrayTag>::chunk_begin()");
        return chunk_const_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    chunk_const_iterator chunk_end(shape_type const & start, shape_type const & stop) const
    {
        return chunk_begin(start, stop).getEndIterator();
    }

    chunk_const_iterator chunk_cbegin(shape_type const & start, shape_type const & stop) const
    {
        checkSubarrayBounds(start, stop, "MultiArrayView<N-1, T, ChunkedArrayTag>::chunk_cbegin()");
        return chunk_const_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    chunk_const_iterator chunk_cend(shape_type const & start, shape_type const & stop) const
    {
        return chunk_cbegin(start, stop).getEndIterator();
    }

    view_type view ()
    {
        return *this;
    }

    MultiArray<N, Chunk> chunks_;
    shape_type offset_, bits_, mask_;
    VIGRA_SHARED_PTR<ChunkUnrefProxyBase> unref_;
};

template <unsigned int N, class T>
typename MultiArrayView<N, T, ChunkedArrayTag>::iterator
createCoupledIterator(MultiArrayView<N, T, ChunkedArrayTag> & m)
{
    typedef typename MultiArrayView<N, T, ChunkedArrayTag>::iterator    IteratorType;
    typedef typename IteratorType::handle_type           P1;
    typedef typename P1::base_type                       P0;

    return IteratorType(P1(m,
                        P0(m.shape())));
}

template <unsigned int N, class T>
typename MultiArrayView<N, T, ChunkedArrayTag>::const_iterator
createCoupledIterator(MultiArrayView<N, T, ChunkedArrayTag> const & m)
{
    typedef typename MultiArrayView<N, T, ChunkedArrayTag>::const_iterator    IteratorType;
    typedef typename IteratorType::handle_type           P1;
    typedef typename P1::base_type                       P0;

    return IteratorType(P1(m,
                        P0(m.shape())));
}

/** \addtogroup ChunkedArrayClasses Chunked arrays

    Store big data (potentially larger than RAM) as a collection of rectangular blocks.
*/
//@{

/** \brief Option object for \ref ChunkedArray construction.
*/
class ChunkedArrayOptions
{
  public:
    /** \brief Initialize options with defaults.
    */
    ChunkedArrayOptions()
    : fill_value(0.0)
    , cache_max(-1)
    , compression_method(DEFAULT_COMPRESSION)
    {}

    /** \brief Element value for read-only access of uninitialized chunks.

        Default: 0
    */
    ChunkedArrayOptions & fillValue(double v)
    {
        fill_value = v;
        return *this;
    }

    ChunkedArrayOptions fillValue(double v) const
    {
        return ChunkedArrayOptions(*this).fillValue(v);
    }

    /** \brief Maximum number of chunks in the cache.

        Default: -1 ( = use a heuristic depending on array shape)
    */
    ChunkedArrayOptions & cacheMax(int v)
    {
        cache_max = v;
        return *this;
    }

    ChunkedArrayOptions cacheMax(int v) const
    {
        return ChunkedArrayOptions(*this).cacheMax(v);
    }

    /** \brief Compress inactive chunks with the given method.

        Default: DEFAULT_COMPRESSION (depends on backend)
    */
    ChunkedArrayOptions & compression(CompressionMethod v)
    {
        compression_method = v;
        return *this;
    }

    ChunkedArrayOptions compression(CompressionMethod v) const
    {
        return ChunkedArrayOptions(*this).compression(v);
    }

    double fill_value;
    int cache_max;
    CompressionMethod compression_method;
};

/** \brief Interface and base class for chunked arrays.

Very big data arrays (possibly bigger than the available RAM) can
only be processed in smaller pieces. To support quick access to
these pieces, it is advantegeous to store big arrays in chunks,
i.e. as a collection of small rectagular subarrays. The class
ChunkedArray encapsulates storage and handling of these chunks and
provides various APIs to easily access the data.

<b>\#include</b> \<vigra/multi_array_chunked.hxx\> <br/>
Namespace: vigra

@tparam N the array dimension
@tparam T the type of the array elements

(these are the same as in \ref MultiArrayView). The actual way of chunk storage is determined by the derived class the program uses:

<ul>
    <li>ChunkedArrayFull: Provides the chunked array API for a standard
    \ref MultiArray (i.e. there is only one chunk for the entire array).

    <li>ChunkedArrayLazy: All chunks reside in memory, but are only
    allocated upon first access.

    <li>ChunkedArrayCompressed: Like ChunkedArrayLazy, but temporarily
    unused chunks are compressed in memory to save space.

    <li>ChunkedArrayTmpFile: Chunks are stored in a memory-mapped file.
    Temporarily unused chunks are written to the hard-drive and deleted from
    memory.

    <li>ChunkedArrayHDF5: Chunks are stored in a HDF5 dataset by means of
    HDF5's native chunked storage capabilities. Temporarily unused chunks are
    written to the hard-drive in compressed form and deleted from memory.
</ul>
You must use these derived classes to construct a chunked array because
ChunkedArray itself is an abstract class.

Chunks can be in one of the following states:
<ul>
    <li>uninitialized: Chunks are only initialized (i.e. allocated) upon the first
    write access. If an uninitialized chunk is accessed in a read-only manner, the
    system returns a pseudo-chunk whose elements have a user-provided fill value.

    <li>asleep: The chunk is currently unused and has been compressed and/or
    swapped out to the hard drive.

    <li>inactive: The chunk is currently unused, but still resides in memory.

    <li>active: The chunk resides in memory and is currently in use.

    <li>locked: Chunks are briefly in this state during transitions
    between the other states (e.g. while loading and/or decompression is
    in progress).

    <li>failed: An unexpected error occured, e.g. the system is out of memory
    or a write to the hard drive failed.
</ul>
In-memory chunks (active and inactive) are placed in a cache. If a chunk
transitions from the 'asleep' to the 'active' state, it is added to the cache,
and an 'inactive' chunk is removed and sent 'asleep'. If there is no 'inactive'
chunk in the cache, the cache size is temporarily increased. All state
transitions are thread-safe.

In order to optimize performance, the user should adjust the cache size (via
\ref setCacheMaxSize() or \ref ChunkedArrayOptions) so that it can hold all
chunks that are frequently needed (e.g. all chunks forming a row of the full
array).

Another performance critical parameter is the chunk shape. While the system
uses sensible defaults (512<sup>2</sup> for 2D arrays, 64<sup>3</sup> for 3D,
64x64x16x4 for 4D, and 64x64x16x4x4 for 5D), the shape may need to be adjusted
via the array's constructor to match the access patterns of the algorithms to
be used. For speed reasons, chunk shapes must be powers of 2.

The data in the array can be accessed in several ways. The simplest is
via calls to <tt>checkoutSubarray()</tt> and <tt>commitSubarray()</tt>: These
functions copy an arbitrary subregion of a chunked array (possibly straddling
many chunks) into a standard \ref MultiArrayView for processing, and write
results back into the chunked array:
\code
    ChunkedArray<3, float> & chunked_array = ...;

    Shape3 roi_start(1000, 500, 500);
    MultiArray<3, float> work_array(Shape3(100, 100, 100));

    // copy data from region (1000,500,500)...(1100,600,600)
    chunked_array.checkoutSubarray(roi_start, work_array);

    ... // work phase: process data in work_array as usual

    // write results back into chunked_array
    chunked_array.commitSubarray(roi_start, work_array);
\endcode
The required chunks in <tt>chunked_array</tt> will only be active while the
checkout and commit calls are executing. During the work phase, other threads
can use the chunked array's cache to checkout or commit different subregions.

Alternatively, one can work directly on the chunk storage. This is most easily
achieved by means of chunk iterators:
\code
    ChunkedArray<3, float> & chunked_array = ...;

    // define the ROI to be processed
    Shape3 roi_start(100, 200, 300), roi_end(1000, 2000, 600);

    // get a pair of chunk iterators ( = iterators over chunks)
    auto chunk = chunked_array.chunk_begin(roi_start, roi_end),
         end   = chunked_array.chunk_end(roi_start, roi_end);

    // iterate over the chunks in the ROI
    for(; chunk != end; ++chunk)
    {
        // get a view to the current chunk's data
        // Note: The view actually refers to the intersection of the
        //       current chunk with the ROI. Thus, chunks which are
        //       partially outside the ROI are appropriately trimmed.
        MultiArrayView<3, float> chunk_view = *chunk;

        ... // work phase: process data in chunk_view as usual
    }
\endcode
No memory is duplicated in this approach, and only the current chunk needs
to be active, so that a small chunk cache is sufficient. The iteration
over chunks can be distributed over several threads that process the array
data in parallel. The programmer must make sure that write operations to
individual elements are synchronized between threads. This is usually
achieved by ensuring that the threads are responsible for non-overlapping
regions of the output array.

An even simpler method is direct element access via indexing. However, the
chunked array has no control over the access order in this case, so it must
potentially activate the present chunk upon each access. This is rather
expensive and should only be used for debugging:
\code
    ChunkedArray<3, float> & chunked_array = ...;

    Shape3 index(100, 200, 300);
    // access data at coordinate 'index'
    chunked_array.setItem(index, chunked_array.getItem(index) + 2.0);
\endcode

Two additional APIs provide access in a way compatible with an ordinary
\ref MultiArrayView. These APIs should be used in functions that are
supposed to work unchanged on both ordinary and chunked arrays. The first
possibility is the chunked scan-order iterator:
\code
    ChunkedArray<3, float> & chunked_array = ...;

    // get a pair of scan-order iterators ( = iterators over elements)
    auto iter = chunked_array.begin(),
         end  = chunked_array.end();

    // iterate over all array elements
    for(; iter != end; ++iter)
    {
        // access current element
        *iter = *iter + 2.0;
    }
\endcode
A new chunk must potentially be activated whenever the iterator crosses
a chunk boundary. Since the overhead of the activation operation can be
amortized over many within-chunk steps, the iteration (excluding the
workload within the loop) takes only twice as long as the iteration over an
unstrided array using an ordinary \ref StridedScanOrderIterator.

The final possibility is the creation of a MultiArrayView that accesses
an arbitrary ROI directly:
\code
    ChunkedArray<3, float> & chunked_array = ...;

    // define the ROI to be processed
    Shape3 roi_start(100, 200, 300), roi_end(1000, 2000, 600);

    // create view for ROI
    MultiArrayView<3, float, ChunkedArrayTag> view =
                    chunked_array.subarray(roi_start, roi_stop);

    ... // work phase: process view like any ordinary MultiArrayView
\endcode
Similarly, a lower-dimensional view can be created with one of the
<tt>bind</tt> functions. This approach has the advantage that 'view'
can be passed to any function which is implemented in terms of
MultiArrayViews. However, there are two disadvantages: First, data access
in the view requires two steps (first find the chunk, then find the
appropriate element in the chunk), which causes the chunked view to
be slower than an ordinary MultiArrayView. Second, all chunks intersected
by the view must remain active throughout the view's lifetime, which
may require a big chunk cache and thus keeps many chunks in memory.
*/
template <unsigned int N, class T>
class ChunkedArray
: public ChunkedArrayBase<N, T>
{
    /*
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
    * provide a ChunkIterator that iterates over all chunks in a given ROI and returns a
      MultiArrayView for the present chunk (which remains locked in cache until the
      iterator is advanced).
    * implement proper copy constructors and assignment for all backends
    * test HDF5 constructor from existing dataset
    * put HDF5 into header of its own
    * is the full chunkForIterator() function slow? Check this with a simplified one
      in a ChunkedArrayLazy where all chunlks are already implemented, so that
      we can simply can skip the check
    * add support for Multiband and TinyVector pixels

    */

  public:
    typedef ChunkedArrayBase<N, T> base_type;
    typedef typename MultiArrayShape<N>::type  shape_type;
    typedef typename shape_type::value_type  difference_type_1;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef ChunkIterator<N, T>         chunk_iterator;
    typedef ChunkIterator<N, T const>   chunk_const_iterator;
    typedef StridedScanOrderIterator<N, ChunkedMemory<T>, reference, pointer>   iterator;
    typedef StridedScanOrderIterator<N, ChunkedMemory<T const>, const_reference, const_pointer>   const_iterator;
    typedef SharedChunkHandle<N, T> Handle;
    typedef ChunkBase<N, T> Chunk;
    typedef MultiArrayView<N, T, ChunkedArrayTag>                   view_type;
    typedef MultiArrayView<N, T const, ChunkedArrayTag>             const_view_type;
    typedef std::queue<Handle*> CacheType;

    static const long chunk_asleep = Handle::chunk_asleep;
    static const long chunk_uninitialized = Handle::chunk_uninitialized;
    static const long chunk_locked = Handle::chunk_locked;
    static const long chunk_failed = Handle::chunk_failed;

    // constructor only called by derived classes (ChunkedArray is abstract)
    explicit ChunkedArray(shape_type const & shape,
                          shape_type const & chunk_shape = shape_type(),
                          ChunkedArrayOptions const & options = ChunkedArrayOptions())
    : ChunkedArrayBase<N, T>(shape, chunk_shape)
    , bits_(initBitMask(this->chunk_shape_))
    , mask_(this->chunk_shape_ -shape_type(1))
    , cache_max_size_(options.cache_max)
    , chunk_lock_(new threading::mutex())
    , fill_value_(T(options.fill_value))
    , fill_scalar_(options.fill_value)
    , handle_array_(detail::computeChunkArrayShape(shape, bits_, mask_))
    , data_bytes_()
    , overhead_bytes_(handle_array_.size()*sizeof(Handle))
    {
        fill_value_chunk_.pointer_ = &fill_value_;
        fill_value_handle_.pointer_ = &fill_value_chunk_;
        fill_value_handle_.chunk_state_.store(1);
    }

    // compute masks needed for fast index access
    static shape_type initBitMask(shape_type const & chunk_shape)
    {
        shape_type res;
        for(unsigned int k=0; k<N; ++k)
        {
            UInt32 bits = log2i(chunk_shape[k]);
            vigra_precondition(chunk_shape[k] == MultiArrayIndex(1 << bits),
                               "ChunkedArray: chunk_shape elements must be powers of 2.");
            res[k] = bits;
         }
         return res;
    }

    virtual ~ChunkedArray()
    {
        // std::cerr << "    final cache size: " << cacheSize() << " (max: " << cacheMaxSize() << ")\n";
    }

    /** \brief Number of chunks currently fitting into the cache.
    */
    int cacheSize() const
    {
        return cache_.size();
    }

    /** \brief Bytes of main memory occupied by the array's data.

        Compressed chunks are only counted with their compressed size.
        Chunks swapped out to the hard drive are not counted.
    */
    std::size_t dataBytes() const
    {
        return data_bytes_;
    }

    /** \brief Bytes of main memory needed to manage the chunked storage.
    */
    std::size_t overheadBytes() const
    {
        return overhead_bytes_;
    }

    /** \brief Number of chunks along each coordinate direction.
    */
    virtual shape_type chunkArrayShape() const
    {
        return handle_array_.shape();
    }

    virtual std::size_t dataBytes(Chunk * c) const = 0;

    /** \brief Number of data bytes in an uncompressed chunk.
    */
    std::size_t dataBytesPerChunk() const
    {
        return prod(this->chunk_shape_)*sizeof(T);
    }

    /** \brief Bytes of main memory needed to manage a single chunk.
    */
    virtual std::size_t overheadBytesPerChunk() const = 0;

    /** \brief Find the chunk that contains array element 'global_start'.
    */
    shape_type chunkStart(shape_type const & global_start) const
    {
        shape_type chunk_start(SkipInitialization);
        detail::ChunkIndexing<N>::chunkIndex(global_start, bits_, chunk_start);
        return chunk_start;
    }

    /** \brief Find the chunk that is beyond array element 'global_stop'.

        Specifically, this computes
        \code
        chunkStart(global_stop - shape_type(1)) + shape_type(1)
        \endcode
    */
    shape_type chunkStop(shape_type global_stop) const
    {
        global_stop -= shape_type(1);
        shape_type chunk_stop(SkipInitialization);
        detail::ChunkIndexing<N>::chunkIndex(global_stop, bits_, chunk_stop);
        chunk_stop += shape_type(1);
        return chunk_stop;
    }

    /** \brief Find the shape of the chunk indexed by 'chunk_index'.

         This may differ from the global chunk shape because chunks at the
         right/lower border of the array may be smaller than usual.
    */
    shape_type chunkShape(shape_type const & chunk_index) const
    {
        return min(this->chunk_shape_,
                   this->shape_ - chunk_index*this->chunk_shape_);
    }

    using base_type::chunkShape;

#ifdef DOXYGEN
    /** \brief Return the global chunk shape.

        This is the shape of all chunks that are completely contained
        in the array's domain.
    */
    shape_type const & chunkShape() const;

    /** \brief Return the shape in this array.
    */
    shape_type const & shape() const;

    /** \brief Return the number of elements in this array.
    */
    MultiArrayIndex size() const;

    /** \brief Check if the given point is in the array domain.
    */
    bool isInside(shape_type const & p) const;

    /** \brief Return the class that implements this ChunkedArray.
    */
    std::string backend() const;

#endif

    inline void
    checkSubarrayBounds(shape_type const & start, shape_type const & stop,
                        std::string message) const
    {
        message += ": subarray out of bounds.";
        vigra_precondition(allLessEqual(shape_type(), start) &&
                           allLess(start, stop) &&
                           allLessEqual(stop, this->shape_),
                           message);
    }

    /** \brief Check if two arrays are elementwise equal.
    */
    template <class U, class C1>
    bool operator==(MultiArrayView<N, U, C1> const & rhs) const
    {
        if(this->shape() != rhs.shape())
            return false;
        const_iterator i = begin(), ie = end();
        typename MultiArrayView<N, U, C1>::const_iterator j = rhs.begin();
        for(; i != ie; ++i, ++j)
            if(*i != *j)
                return false;
        return true;
    }

    /** \brief Check if two arrays differ in at least one element.
    */
    template <class U, class C1>
    bool operator!=(MultiArrayView<N, U, C1> const & rhs) const
    {
        return !operator==(rhs);
    }

    // internal function to activate a chunk
    virtual pointer loadChunk(Chunk ** chunk, shape_type const & chunk_index) = 0;

    // internal function to send a chunk asleep or delete it
    // entirely (when destroy = true).
    // returns true if the chunk was deleted, false otherwise
    virtual bool unloadHandle(Handle * handle, bool destroy = false)
    {
        if(handle == &fill_value_handle_)
            return false;
        return unloadChunk(handle->pointer_, destroy);
    }

    virtual bool unloadChunk(Chunk * chunk, bool destroy = false) = 0;

    Handle * lookupHandle(shape_type const & index)
    {
        return &handle_array_[index];
    }

    // Decrease the reference counter of the given chunk.
    // Will inactivate the chunk when reference counter reaches zero.
    virtual void unrefChunk(IteratorChunkHandle<N, T> * h) const
    {
        unrefChunk(h->chunk_);
        h->chunk_ = 0;
    }

    // Likewise
    void unrefChunk(Handle * chunk) const
    {
        if(chunk)
        {
            long rc = chunk->chunk_state_.fetch_sub(1);
          #ifdef VIGRA_CHECK_BOUNDS
            vigra_invariant(rc >= 0,
                            "ChunkedArray::unrefChunk(): chunk refcount got negative!");
          #endif
        }
    }

    // Decrease the reference counter of several chunks simultaneously.
    void unrefChunks(ArrayVector<Handle*> const & chunks)
    {
        for(unsigned int k=0; k<chunks.size(); ++k)
            unrefChunk(chunks[k]);

        if(cacheMaxSize() > 0)
        {
            threading::lock_guard<threading::mutex> guard(*chunk_lock_);
            cleanCache(cache_.size());
        }
    }

    // Increase the reference counter of the given chunk.
    // If the chunk was asleep, the function first awakens it.
    long acquireRef(Handle * handle) const
    {
        // Obtain a reference to the current chunk handle.
        // We use a simple spin-lock here because it is very fast in case of success,
        // and failures (i.e. collisions with another thread) are presumably
        // very rare.
        //
        // the function returns the old value of chunk_state_
        long rc = handle->chunk_state_.load(threading::memory_order_acquire);
        while(true)
        {
            if(rc >= 0)
            {
                if(handle->chunk_state_.compare_exchange_weak(rc, rc+1, threading::memory_order_seq_cst))
                {
                    return rc;
                }
            }
            else
            {
                if(rc == chunk_failed)
                {
                    vigra_precondition(false,
                     "ChunkedArray::acquireRef() attempt to access failed chunk.");
                }
                else if(rc == chunk_locked)
                {
                    // cache management in progress => try again later
                    threading::this_thread::yield();
                    rc = handle->chunk_state_.load(threading::memory_order_acquire);
                }
                else if(handle->chunk_state_.compare_exchange_weak(rc, chunk_locked, threading::memory_order_seq_cst))
                {
                    return rc;
                }
            }
        }
    }

    pointer
    getChunk(Handle * handle, bool isConst, bool insertInCache, shape_type const & chunk_index) const
    {
        ChunkedArray * self = const_cast<ChunkedArray *>(this);

        long rc = acquireRef(handle);
        if(rc >= 0)
            return handle->pointer_->pointer_;

        threading::lock_guard<threading::mutex> guard(*chunk_lock_);
        try
        {
            T * p = self->loadChunk(&handle->pointer_, chunk_index);
            Chunk * chunk = handle->pointer_;
            if(!isConst && rc == chunk_uninitialized)
                std::fill(p, p + prod(chunkShape(chunk_index)), this->fill_value_);

            self->data_bytes_ += dataBytes(chunk);

            if(cacheMaxSize() > 0 && insertInCache)
            {
                // insert in queue of mapped chunks
                self->cache_.push(handle);

                // do cache management if cache is full
                // (note that we still hold the chunk_lock_)
                self->cleanCache(2);
            }
            handle->chunk_state_.store(1, threading::memory_order_release);
            return p;
        }
        catch(...)
        {
            handle->chunk_state_.store(chunk_failed);
            throw;
        }
    }

    // helper function for chunkForIterator()
    inline pointer
    chunkForIteratorImpl(shape_type const & point,
                         shape_type & strides, shape_type & upper_bound,
                         IteratorChunkHandle<N, T> * h,
                         bool isConst) const
    {
        ChunkedArray * self = const_cast<ChunkedArray *>(this);

        unrefChunk(h->chunk_);
        h->chunk_ = 0;

        shape_type global_point = point + h->offset_;

        if(!this->isInside(global_point))
        {
            upper_bound = point + this->chunk_shape_;
            return 0;
        }

        shape_type chunkIndex(chunkStart(global_point));

        bool insertInCache = true;
        Handle * handle = self->lookupHandle(chunkIndex);
        if(isConst && handle->chunk_state_.load() == chunk_uninitialized)
        {
            handle = &self->fill_value_handle_;
            insertInCache = false;
        }

        pointer p = getChunk(handle, isConst, insertInCache, chunkIndex);
        strides = handle->strides();
        upper_bound = (chunkIndex + shape_type(1)) * this->chunk_shape_ - h->offset_;
        std::size_t offset = detail::ChunkIndexing<N>::offsetInChunk(global_point, mask_, strides);
        h->chunk_ = handle;
        return p + offset;
    }

    // called by chunked scan-order iterator to obtain the new data pointer
    // when the iterator enters a new chunk
    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h)
    {
        return chunkForIteratorImpl(point, strides, upper_bound, h, false);
    }

    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h) const
    {
        return chunkForIteratorImpl(point, strides, upper_bound, h, true);
    }

    // NOTE: This function must only be called while we hold the chunk_lock_.
    //       This implies refcount != chunk_locked, so that race conditions are avoided.
    long releaseChunk(Handle * handle, bool destroy = false)
    {
        long rc = 0;
        bool mayUnload = handle->chunk_state_.compare_exchange_strong(rc, chunk_locked);
        if(!mayUnload && destroy)
        {
            rc = chunk_asleep;
            mayUnload = handle->chunk_state_.compare_exchange_strong(rc, chunk_locked);
        }
        if(mayUnload)
        {
            // refcount was zero or chunk_asleep => can unload
            try
            {
                vigra_invariant(handle != &fill_value_handle_,
                   "ChunkedArray::releaseChunk(): attempt to release fill_value_handle_.");
                Chunk * chunk = handle->pointer_;
                this->data_bytes_ -= dataBytes(chunk);
                int didDestroy = unloadChunk(chunk, destroy);
                this->data_bytes_ += dataBytes(chunk);
                if(didDestroy)
                    handle->chunk_state_.store(chunk_uninitialized);
                else
                    handle->chunk_state_.store(chunk_asleep);
            }
            catch(...)
            {
                handle->chunk_state_.store(chunk_failed);
                throw;
            }
        }
        return rc;
    }

    // NOTE: this function must only be called while we hold the chunk_lock_
    void cleanCache(int how_many = -1)
    {
        if(how_many == -1)
            how_many = cache_.size();
        for(; cache_.size() > cacheMaxSize() && how_many > 0; --how_many)
        {
            Handle * handle = cache_.front();
            cache_.pop();
            long rc = releaseChunk(handle);
            if(rc > 0) // refcount was positive => chunk is still needed
                cache_.push(handle);
        }
    }

    /** Sends all chunks asleep which are completely inside the given ROI.
        If destroy == true and the backend supports destruction (currently:
        ChunkedArrayLazy and ChunkedArrayCompressed), chunks will be deleted
        entirely. The chunk's contents after releaseChunks() are undefined.
        Currently, chunks retain their values when sent asleep, and assume the
        array's fill_value when deleted, but applications should not rely on this
        behavior.
    */
    void releaseChunks(shape_type const & start, shape_type const & stop, bool destroy = false)
    {
        checkSubarrayBounds(start, stop, "ChunkedArray::releaseChunks()");

        MultiCoordinateIterator<N> i(chunkStart(start), chunkStop(stop)),
                                   end(i.getEndIterator());
        for(; i != end; ++i)
        {
            shape_type chunkOffset = *i * this->chunk_shape_;
            if(!allLessEqual(start, chunkOffset) ||
               !allLessEqual(min(chunkOffset+this->chunk_shape_, this->shape()), stop))
            {
                // chunk is only partially covered by the ROI
                continue;
            }

            Handle * handle = this->lookupHandle(*i);
            threading::lock_guard<threading::mutex> guard(*chunk_lock_);
            releaseChunk(handle, destroy);
        }

        // remove all chunks from the cache that are asleep or unitialized
        threading::lock_guard<threading::mutex> guard(*chunk_lock_);
        int cache_size = cache_.size();
        for(int k=0; k < cache_size; ++k)
        {
            Handle * handle = cache_.front();
            cache_.pop();
            if(handle->chunk_state_.load() >= 0)
                cache_.push(handle);
        }
    }

    /** \brief Copy an ROI of the chunked array into an ordinary MultiArrayView.

        The ROI's lower bound is given by 'start', its upper bound (in 'beyond' sense)
        is 'start + subarray.shape()'. Chunks in the ROI are only activated while
        the read is in progress.
    */
    template <class U, class Stride>
    void
    checkoutSubarray(shape_type const & start,
                     MultiArrayView<N, U, Stride> & subarray) const
    {
        shape_type stop   = start + subarray.shape();

        checkSubarrayBounds(start, stop, "ChunkedArray::checkoutSubarray()");

        chunk_const_iterator i = chunk_cbegin(start, stop);
        for(; i.isValid(); ++i)
        {
            subarray.subarray(i.chunkStart()-start, i.chunkStop()-start) = *i;
        }
    }

    /** \brief Copy an ordinary MultiArrayView into an ROI of the chunked array.

        The ROI's lower bound is given by 'start', its upper bound (in 'beyond' sense)
        is 'start + subarray.shape()'. Chunks in the ROI are only activated while
        the write is in progress.
    */
    template <class U, class Stride>
    void
    commitSubarray(shape_type const & start,
                   MultiArrayView<N, U, Stride> const & subarray)
    {
        shape_type stop   = start + subarray.shape();

        vigra_precondition(!this->isReadOnly(),
                           "ChunkedArray::commitSubarray(): array is read-only.");
        checkSubarrayBounds(start, stop, "ChunkedArray::commitSubarray()");

        chunk_iterator i = chunk_begin(start, stop);
        for(; i.isValid(); ++i)
        {
            *i = subarray.subarray(i.chunkStart()-start, i.chunkStop()-start);
        }
    }

    // helper function for subarray()
    template <class View>
    void subarrayImpl(shape_type const & start, shape_type const & stop,
                      View & view,
                      bool isConst) const
    {
        vigra_precondition(isConst || !this->isReadOnly(),
                           "ChunkedArray::subarray(): array is read-only.");
        checkSubarrayBounds(start, stop, "ChunkedArray::subarray()");
        shape_type chunk_start(chunkStart(start)), chunk_stop(chunkStop(stop));

        view.shape_ = stop-start;
        view.chunk_shape_ = this->chunk_shape_;
        view.chunks_.reshape(chunk_stop-chunk_start);
        view.offset_ = start - chunk_start * this->chunk_shape_;
        view.bits_   = bits_;
        view.mask_   = mask_;

        typedef typename View::UnrefProxy Unref;
        ChunkedArray* self = const_cast<ChunkedArray*>(this);
        Unref * unref = new Unref(view.chunks_.size(), self);
        view.unref_ = VIGRA_SHARED_PTR<Unref>(unref);

        MultiCoordinateIterator<N> i(chunk_start, chunk_stop),
                                   end(i.getEndIterator());
        for(; i != end; ++i)
        {
            Handle * handle = self->lookupHandle(*i);

            if(isConst && handle->chunk_state_.load() == chunk_uninitialized)
                handle = &self->fill_value_handle_;

            // This potentially acquires the chunk_lock_ in each iteration.
            // Would it be better to acquire it once before the loop?
            pointer p = getChunk(handle, isConst, true, *i);

            ChunkBase<N, T> * mini_chunk = &view.chunks_[*i - chunk_start];
            mini_chunk->pointer_ = p;
            mini_chunk->strides_ = handle->strides();
            unref->chunks_[i.scanOrderIndex()] = handle;
        }
    }

    /** \brief Create a view to the specified ROI.

        The view can be used like an ordinary \ref MultiArrayView, but is
        a but slower. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    view_type
    subarray(shape_type const & start, shape_type const & stop)
    {
        view_type view;
        subarrayImpl(start, stop, view, false);
        return view;
    }

    /** \brief Create a read-only view to the specified ROI.

        The view can be used like an ordinary \ref MultiArrayView, but is
        a but slower. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    const_view_type
    subarray(shape_type const & start, shape_type const & stop) const
    {
        const_view_type view;
        subarrayImpl(start, stop, view, true);
        return view;
    }

    /** \brief Create a read-only view to the specified ROI.

        The view can be used like an ordinary \ref MultiArrayView, but is
        a but slower. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    const_view_type
    const_subarray(shape_type const & start, shape_type const & stop) const
    {
        const_view_type view;
        subarrayImpl(start, stop, view, true);
        return view;
    }

    /** \brief Read the array element at index 'point'.

        Since the corresponding chunk must potentially be activated
        first, this function may be slow and should mainly be used in
        debugging.
    */
    value_type getItem(shape_type const & point) const
    {
        vigra_precondition(this->isInside(point),
            "ChunkedArray::getItem(): index out of bounds.");

        ChunkedArray * self = const_cast<ChunkedArray*>(this);
        shape_type chunk_index(chunkStart(point));
        Handle * handle = self->lookupHandle(chunk_index);
        if(handle->chunk_state_.load() == chunk_uninitialized)
            return fill_value_;
        pointer p = self->getChunk(handle, true, false, chunk_index);
        value_type res = *(p +
                           detail::ChunkIndexing<N>::offsetInChunk(point, mask_, handle->strides()));
        self->unrefChunk(handle);
        return res;
    }

    /** \brief Write the array element at index 'point'.

        Since the corresponding chunk must potentially be activated
        first, this function may be slow and should mainly be used in
        debugging.
    */
    void setItem(shape_type const & point, value_type const & v)
    {
        vigra_precondition(!this->isReadOnly(),
            "ChunkedArray::setItem(): array is read-only.");
        vigra_precondition(this->isInside(point),
            "ChunkedArray::setItem(): index out of bounds.");

        shape_type chunk_index(chunkStart(point));
        Handle * handle = lookupHandle(chunk_index);
        pointer p = getChunk(handle, false, false, chunk_index);
        *(p + detail::ChunkIndexing<N>::offsetInChunk(point, mask_, handle->strides())) = v;
        unrefChunk(handle);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        Dimension 'dim' is bound at 'index', all other dimensions remain
        unchanged. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    MultiArrayView<N-1, T, ChunkedArrayTag>
    bindAt(MultiArrayIndex dim, MultiArrayIndex index) const
    {
        shape_type start, stop(this->shape());
        start[dim] = index;
        stop[dim] = index+1;
        return subarray(start, stop).bindAt(dim, 0);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        Dimension 'M' (given as a template parameter) is bound at 'index',
        all other dimensions remain unchanged. All chunks intersecting the
        view remain active throughout the view's lifetime.
    */
    template <unsigned int M>
    MultiArrayView <N-1, T, ChunkedArrayTag>
    bind (difference_type_1 index) const
    {
        return bindAt(M, index);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        Dimension 'N-1' is bound at 'index', all other dimensions remain
        unchanged. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    MultiArrayView <N-1, T, ChunkedArrayTag>
    bindOuter (difference_type_1 index) const
    {
        return bindAt(N-1, index);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        The M rightmost dimensions are bound to the indices given in 'd'.
        All chunks intersecting the view remain active throughout the view's lifetime.
    */
    template <int M, class Index>
    MultiArrayView <N-M, T, ChunkedArrayTag>
    bindOuter(const TinyVector <Index, M> & d) const
    {
        return bindAt(N-1, d[M-1]).bindOuter(d.dropIndex(M-1));
    }

    // terminate the recursion of the above function
    template <class Index>
    MultiArrayView <N-1, T, ChunkedArrayTag>
    bindOuter(const TinyVector <Index, 1> & d) const
    {
        return bindAt(N-1, d[0]);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        Dimension '0' is bound at 'index', all other dimensions remain
        unchanged. All chunks intersecting the view remain active
        throughout the view's lifetime.
    */
    MultiArrayView <N-1, T, ChunkedArrayTag>
    bindInner (difference_type_1 index) const
    {
        return bindAt(0, index);
    }

    /** \brief Create a lower dimensional view to the chunked array.

        The M leftmost dimensions are bound to the indices given in 'd'.
        All chunks intersecting the view remain active throughout the view's lifetime.
    */
    template <int M, class Index>
    MultiArrayView <N-M, T, ChunkedArrayTag>
    bindInner(const TinyVector <Index, M> & d) const
    {
        return bindAt(0, d[0]).bindInner(d.dropIndex(0));
    }

    // terminate the recursion of the above function
    template <class Index>
    MultiArrayView <N-1, T, ChunkedArrayTag>
    bindInner(const TinyVector <Index, 1> & d) const
    {
        return bindAt(0, d[0]);
    }

    /** \brief Get the number of chunks the cache will hold.

        If there are any inactive chunks in the cache, these will be
        sent asleep until the max cahce size is reached. The max cache
        size may be temporarily overridden when more chunks need
        to be active simultaneously.
    */
    std::size_t cacheMaxSize() const
    {
        if(cache_max_size_ < 0)
            const_cast<int &>(cache_max_size_) = detail::defaultCacheSize(this->chunkArrayShape());
        return cache_max_size_;
    }

    /** \brief Set the number of chunks the cache will hold.

        This should be big enough to hold all chunks that are frequently needed
        and must therefore be adopted to the application's access pattern.
    */
    void setCacheMaxSize(std::size_t c)
    {
        cache_max_size_ = c;
        if(c < cache_.size())
        {
            threading::lock_guard<threading::mutex> guard(*chunk_lock_);
            cleanCache();
        }
    }

    /** \brief Create a scan-order iterator for the entire chunked array.
    */
    iterator begin()
    {
        return createCoupledIterator(*this);
    }

    /** \brief Create the end iterator for scan-order iteration over
        the entire chunked array.
    */
    iterator end()
    {
        return begin().getEndIterator();
    }

    /** \brief Create a read-only scan-order iterator for the entire
         chunked array.
    */
    const_iterator cbegin() const
    {
        return createCoupledIterator(const_cast<ChunkedArray const &>(*this));
    }

    /** \brief Create the end iterator for read-only scan-order iteration over
        the entire chunked array.
    */
    const_iterator cend() const
    {
        return cbegin().getEndIterator();
    }

    /** \brief Create a read-only scan-order iterator for the entire
         chunked array.
    */
    const_iterator begin() const
    {
        return createCoupledIterator(*this);
    }

    /** \brief Create the end iterator for read-only scan-order iteration over
        the entire chunked array.
    */
    const_iterator end() const
    {
        return begin().getEndIterator();
    }

    /** \brief Create an iterator over all chunks intersected by the given ROI.
    */
    chunk_iterator chunk_begin(shape_type const & start, shape_type const & stop)
    {
        checkSubarrayBounds(start, stop, "ChunkedArray::chunk_begin()");
        return chunk_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    /** \brief Create the end iterator for iteration over all chunks
        intersected by the given ROI.
    */
    chunk_iterator chunk_end(shape_type const & start, shape_type const & stop)
    {
        return chunk_begin(start, stop).getEndIterator();
    }

    /** \brief Create a read-only iterator over all chunks intersected
        by the given ROI.
    */
    chunk_const_iterator chunk_begin(shape_type const & start, shape_type const & stop) const
    {
        checkSubarrayBounds(start, stop, "ChunkedArray::chunk_begin()");
        return chunk_const_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    /** \brief Create the end iterator for read-only iteration over all chunks
        intersected by the given ROI.
    */
    chunk_const_iterator chunk_end(shape_type const & start, shape_type const & stop) const
    {
        return chunk_begin(start, stop).getEndIterator();
    }

    /** \brief Create a read-only iterator over all chunks intersected
        by the given ROI.
    */
    chunk_const_iterator chunk_cbegin(shape_type const & start, shape_type const & stop) const
    {
        checkSubarrayBounds(start, stop, "ChunkedArray::chunk_cbegin()");
        return chunk_const_iterator(this, start, stop, chunkStart(start), chunkStop(stop), this->chunk_shape_);
    }

    /** \brief Create the end iterator for read-only iteration over all chunks
        intersected by the given ROI.
    */
    chunk_const_iterator chunk_cend(shape_type const & start, shape_type const & stop) const
    {
        return chunk_cbegin(start, stop).getEndIterator();
    }

    shape_type bits_, mask_;
    int cache_max_size_;
    VIGRA_SHARED_PTR<threading::mutex> chunk_lock_;
    CacheType cache_;
    Chunk fill_value_chunk_;
    Handle fill_value_handle_;
    value_type fill_value_;
    double fill_scalar_;
    MultiArray<N, Handle> handle_array_;
    std::size_t data_bytes_, overhead_bytes_;
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

template <unsigned int N, class T>
typename ChunkedArray<N, T>::const_iterator
createCoupledIterator(ChunkedArray<N, T> const & m)
{
    typedef typename ChunkedArray<N, T>::const_iterator  IteratorType;
    typedef typename IteratorType::handle_type           P1;
    typedef typename P1::base_type                       P0;

    return IteratorType(P1(m,
                        P0(m.shape())));
}

/** Implement ChunkedArray as an ordinary MultiArray with a single chunk.

    <b>\#include</b> \<vigra/multi_array_chunked.hxx\> <br/>
    Namespace: vigra
*/
template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayFull
: public ChunkedArray<N, T>,
  public MultiArray<N, T, Alloc>
{
  public:

    typedef MultiArray<N, T, Alloc>             Storage;
    typedef typename Storage::value_type        value_type;
    typedef typename Storage::pointer           pointer;
    typedef typename Storage::const_pointer     const_pointer;
    typedef typename Storage::reference         reference;
    typedef typename Storage::const_reference   const_reference;
    typedef typename Storage::difference_type   difference_type;
    typedef typename Storage::difference_type   shape_type;
    typedef typename Storage::key_type          key_type;
    typedef typename Storage::size_type         size_type;
    typedef typename Storage::difference_type_1 difference_type_1;
    typedef typename Storage::iterator          iterator;
    typedef typename Storage::const_iterator    const_iterator;
    typedef typename Storage::view_type         view_type;

    typedef typename ChunkedArray<N, T>::Chunk       Chunk;

    static shape_type computeChunkShape(shape_type s)
    {
        for(int k=0; k<N; ++k)
            s[k] = ceilPower2(s[k]);
        return s;
    }

    using Storage::subarray;
    using Storage::bindOuter;
    using Storage::bindInner;
    using Storage::bind;
    using Storage::bindAt;
    using Storage::isInside;
    using Storage::shape;
    using Storage::size;
    using Storage::begin;
    using Storage::end;

#ifndef DOXYGEN  // doxygen doesn't understand this
    using Storage::operator==;
    using Storage::operator!=;
#endif

    /** \brief Construct with given 'shape' and 'options', using the allocator
        'alloc' to manage the memory.
    */
    explicit ChunkedArrayFull(shape_type const & shape,
                              ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                              Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape, computeChunkShape(shape), options.cacheMax(0)),
      Storage(shape, this->fill_value_, alloc),
      upper_bound_(shape),
      chunk_(detail::defaultStride(shape), this->data())
    {
        this->handle_array_[0].pointer_ = &chunk_;
        this->handle_array_[0].chunk_state_.store(1);
        this->data_bytes_ = size()*sizeof(T);
        this->overhead_bytes_ = overheadBytesPerChunk();
    }

    ChunkedArrayFull(ChunkedArrayFull const & rhs)
    : ChunkedArray<N, T>(rhs),
      Storage(rhs),
      upper_bound_(rhs.upper_bound_),
      chunk_(detail::defaultStride(shape), this->data())
    {
        this->handle_array_[0].pointer_ = &chunk_;
        this->handle_array_[0].chunk_state_.store(1);
    }

    ChunkedArrayFull & operator=(ChunkedArrayFull const & rhs)
    {
        if(this != &rhs)
        {
            ChunkedArray<N, T>::operator=(rhs);
            Storage::operator=(rhs);
            upper_bound_ = rhs.upper_bound_;
        }
        return *this;
    }

    ~ChunkedArrayFull()
    {}

    virtual shape_type chunkArrayShape() const
    {
        return shape_type(1);
    }

    virtual pointer loadChunk(ChunkBase<N, T> **, shape_type const &)
    {
        return this->data();
    }

    virtual bool unloadChunk(ChunkBase<N, T> *, bool /* destroy */)
    {
        return false; // never destroys the data
    }

    virtual std::size_t dataBytes(Chunk * c) const
    {
        return prod(this->shape());
    }

    virtual std::size_t overheadBytesPerChunk() const
    {
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>);
    }

    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h) const
    {
        shape_type global_point = point + h->offset_;

        if(!this->isInside(global_point))
        {
            upper_bound = point + this->chunk_shape_;
            return 0;
        }

        strides = this->stride();
        upper_bound = upper_bound_;
        return const_cast<pointer>(&Storage::operator[](global_point));
    }

    virtual pointer chunkForIterator(shape_type const & point,
                                     shape_type & strides, shape_type & upper_bound,
                                     IteratorChunkHandle<N, T> * h)
    {
        shape_type global_point = point + h->offset_;

        if(!this->isInside(global_point))
        {
            upper_bound = point + this->chunk_shape_;
            return 0;
        }

        strides = this->stride();
        upper_bound = upper_bound_;
        return &Storage::operator[](global_point);
    }

    virtual std::string backend() const
    {
        return "ChunkedArrayFull";
    }

    shape_type upper_bound_;
    Chunk chunk_;    // a dummy chunk to fulfill the API
};

/** Implement ChunkedArray as a collection of in-memory chunks.

    This optimizes over an ordinary MultiArray by allocating chunks only
    upon the first write. This is especially useful when only a small
    part of the entire array is actually needed, e.g. in a data viewer.

    <b>\#include</b> \<vigra/multi_array_chunked.hxx\> <br/>
    Namespace: vigra
*/
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

        Chunk(shape_type const & shape, Alloc const & alloc = Alloc())
        : ChunkBase<N, T>(detail::defaultStride(shape))
        , size_(prod(shape))
        , alloc_(alloc)
        {}

        ~Chunk()
        {
            deallocate();
        }

        pointer allocate()
        {
            if(this->pointer_ == 0)
                this->pointer_ = detail::alloc_initialize_n<T>(size_, T(), alloc_);
            return this->pointer_;
        }

        void deallocate()
        {
            detail::destroy_dealloc_n(this->pointer_, size_, alloc_);
            this->pointer_ = 0;
        }

        MultiArrayIndex size_;
        Alloc alloc_;

      private:
        Chunk & operator=(Chunk const &);
    };

    typedef MultiArray<N, SharedChunkHandle<N, T> > ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;

    /** \brief Construct with given 'shape', 'chunk_shape' and 'options',
        using the allocator 'alloc' to manage the memory.
    */
    explicit ChunkedArrayLazy(shape_type const & shape,
                              shape_type const & chunk_shape=shape_type(),
                              ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                              Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape, chunk_shape, options.cacheMax(0))
    , alloc_(alloc)
    {}

    ~ChunkedArrayLazy()
    {
        typename ChunkStorage::iterator i   = this->handle_array_.begin(),
                                        end = this->handle_array_.end();
        for(; i != end; ++i)
        {
            if(i->pointer_)
                delete static_cast<Chunk*>(i->pointer_);
            i->pointer_ = 0;
        }
    }

    virtual pointer loadChunk(ChunkBase<N, T> ** p, shape_type const & index)
    {
        if(*p == 0)
        {
            *p = new Chunk(this->chunkShape(index));
            this->overhead_bytes_ += sizeof(Chunk);
        }
        return static_cast<Chunk *>(*p)->allocate();
    }

    virtual bool unloadChunk(ChunkBase<N, T> * chunk, bool destroy)
    {
        if(destroy)
            static_cast<Chunk *>(chunk)->deallocate();
        return destroy;
    }

    virtual std::string backend() const
    {
        return "ChunkedArrayLazy";
    }

    virtual std::size_t dataBytes(ChunkBase<N,T> * c) const
    {
        return c->pointer_ == 0
                 ? 0
                 : static_cast<Chunk*>(c)->size_*sizeof(T);
    }

    virtual std::size_t overheadBytesPerChunk() const
    {
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>);
    }

    Alloc alloc_;
};

/** Implement ChunkedArray as a collection of potentially compressed
    in-memory chunks.

    This works like \ref ChunkedArrayLazy, but inactive chunks are compressed
    when sent asleep. This is especially appropriate for highly compressible
    data such as label images.

    <b>\#include</b> \<vigra/multi_array_chunked.hxx\> <br/>
    Namespace: vigra
*/
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

        Chunk(shape_type const & shape)
        : ChunkBase<N, T>(detail::defaultStride(shape))
        , compressed_()
        , size_(prod(shape))
        {}

        ~Chunk()
        {
            deallocate();
        }

        pointer allocate()
        {
            if(this->pointer_ == 0)
                this->pointer_ = detail::alloc_initialize_n<T>(size_, T(), alloc_);
            return this->pointer_;
        }

        void deallocate()
        {
            detail::destroy_dealloc_n(this->pointer_, size_, alloc_);
            this->pointer_ = 0;
            compressed_.clear();
        }

        void compress(CompressionMethod method)
        {
            if(this->pointer_ != 0)
            {
                vigra_invariant(compressed_.size() == 0,
                    "ChunkedArrayCompressed::Chunk::compress(): compressed and uncompressed pointer are both non-zero.");

                ::vigra::compress((char const *)this->pointer_, size_*sizeof(T), compressed_, method);

                // std::cerr << "compression ratio: " << double(compressed_.size())/(this->size()*sizeof(T)) << "\n";
                detail::destroy_dealloc_n(this->pointer_, size_, alloc_);
                this->pointer_ = 0;
            }
        }

        pointer uncompress(CompressionMethod method)
        {
            if(this->pointer_ == 0)
            {
                if(compressed_.size())
                {
                    this->pointer_ = alloc_.allocate((typename Alloc::size_type)size_);

                    ::vigra::uncompress(compressed_.data(), compressed_.size(),
                                        (char*)this->pointer_, size_*sizeof(T), method);
                    compressed_.clear();
                }
                else
                {
                    this->pointer_ = allocate();
                }
            }
            else
            {
                vigra_invariant(compressed_.size() == 0,
                    "ChunkedArrayCompressed::Chunk::uncompress(): compressed and uncompressed pointer are both non-zero.");
            }
            return this->pointer_;
        }

        ArrayVector<char> compressed_;
        MultiArrayIndex size_;
        Alloc alloc_;

      private:
        Chunk & operator=(Chunk const &);
    };

    typedef MultiArray<N, SharedChunkHandle<N, T> > ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;

    /** \brief Construct with given 'shape', 'chunk_shape' and 'options'.

        The most important option concerns the compression algorithm. Supported
        algorithms are:
        <ul>
        <li>LZ4: Very fast algorithm that achieves decent compression ratios.
        <li>ZLIB_FAST: Fast compression using 'zlib' (slower than LZ4, but higher compression).
        <li>ZLIB_BEST: Best compression using 'zlib', slow.
        <li>ZLIB_NONE: Use 'zlib' format without compression.
        <li>DEFAULT_COMPRESSION: Same as LZ4.
        </ul>
    */
    explicit ChunkedArrayCompressed(shape_type const & shape,
                                    shape_type const & chunk_shape=shape_type(),
                                    ChunkedArrayOptions const & options = ChunkedArrayOptions())
    : ChunkedArray<N, T>(shape, chunk_shape, options),
       compression_method_(options.compression_method)
    {
        if(compression_method_ == DEFAULT_COMPRESSION)
            compression_method_ = LZ4;
    }

    ~ChunkedArrayCompressed()
    {
        typename ChunkStorage::iterator i   = this->handle_array_.begin(),
                                        end = this->handle_array_.end();
        for(; i != end; ++i)
        {
            if(i->pointer_)
                delete static_cast<Chunk*>(i->pointer_);
            i->pointer_ = 0;
        }
    }

    virtual pointer loadChunk(ChunkBase<N, T> ** p, shape_type const & index)
    {
        if(*p == 0)
        {
            *p = new Chunk(this->chunkShape(index));
            this->overhead_bytes_ += sizeof(Chunk);
        }
        return static_cast<Chunk *>(*p)->uncompress(compression_method_);
    }

    virtual bool unloadChunk(ChunkBase<N, T> * chunk, bool destroy)
    {
        if(destroy)
            static_cast<Chunk *>(chunk)->deallocate();
        else
            static_cast<Chunk *>(chunk)->compress(compression_method_);
        return destroy;
    }

    virtual std::string backend() const
    {
        switch(compression_method_)
        {
          case ZLIB:
            return "ChunkedArrayCompressed<ZLIB>";
          case ZLIB_NONE:
            return "ChunkedArrayCompressed<ZLIB_NONE>";
          case ZLIB_FAST:
            return "ChunkedArrayCompressed<ZLIB_FAST>";
          case ZLIB_BEST:
            return "ChunkedArrayCompressed<ZLIB_BEST>";
          case LZ4:
            return "ChunkedArrayCompressed<LZ4>";
          default:
            return "unknown";
        }
    }

    virtual std::size_t dataBytes(ChunkBase<N,T> * c) const
    {
        return c->pointer_ == 0
                 ? static_cast<Chunk*>(c)->compressed_.size()
                 : static_cast<Chunk*>(c)->size_*sizeof(T);
    }

    virtual std::size_t overheadBytesPerChunk() const
    {
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>);
    }

    CompressionMethod compression_method_;
};

/** Implement ChunkedArray as a collection of chunks that can be
    swapped out into a temporary file when asleep.

    <b>\#include</b> \<vigra/multi_array_chunked.hxx\> <br/>
    Namespace: vigra

    The present implementation uses a memory-mapped sparse file to store the chunks.
    A sparse file is created on Linux using the O_TRUNC flag (this seems to be
    the default file behavior on Linux anyway), and on Windows by
    calling DeviceIoControl(file_handle, FSCTL_SET_SPARSE,...) after file creation.

    The file is automatically deleted upon closing. On Windows, this happens
    because the file was opened with FILE_FLAG_DELETE_ON_CLOSE in combination
    with the flag FILE_ATTRIBUTE_TEMPORARY, which tells the OS to avoid writing
    the file to disk if possible. (However, judging from the timings,
    something is still written, or cleanup takes considerable time.)
    On Linux, automated deletion is achieved via <tt>fileno(tmpfile())</tt>.
*/
template <unsigned int N, class T>
class ChunkedArrayTmpFile
: public ChunkedArray<N, T>
{
    /* REMARKS

    Alternatives are:
    * Don't create a file explicitly, but use the swap file instead. This is
      achieved on Linux by mmap(..., MAP_PRIVATE | MAP_ANONYMOUS, -1, ...),
      on Windows by calling CreateFileMapping(INVALID_HANDLE_VALUE, ...).
       * On Linux, the memory must not be unmapped because this
         looses the data. In fact, anonymous mmap() is very similar to
         malloc(), and there is probably no good reason to use anonymous mmap().
       * On Windows, this is much faster, because the OS will never try to
         actually write something to disk (unless swapping is necessary).
    */
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

        Chunk(shape_type const & shape,
              std::size_t offset, size_t alloc_size,
              FileHandle file)
        : ChunkBase<N, T>(detail::defaultStride(shape))
        , offset_(offset)
        , alloc_size_(alloc_size)
        , file_(file)
        {}

        ~Chunk()
        {
            unmap();
        }

        pointer map()
        {
            if(this->pointer_ == 0)
            {
            #ifdef _WIN32
                static const std::size_t bits = sizeof(DWORD)*8,
                                         mask = (std::size_t(1) << bits) - 1;
                this->pointer_ = (pointer)MapViewOfFile(file_, FILE_MAP_ALL_ACCESS,
                                           std::size_t(offset_) >> bits, offset_ & mask, alloc_size_);
                if(this->pointer_ == 0)
                    winErrorToException("ChunkedArrayChunk::map(): ");
            #else
                this->pointer_ = (pointer)mmap(0, alloc_size_, PROT_READ | PROT_WRITE, MAP_SHARED,
                                  file_, offset_);
                if(this->pointer_ == 0)
                    throw std::runtime_error("ChunkedArrayChunk::map(): mmap() failed.");
            #endif
            }
            return this->pointer_;
        }

        void unmap()
        {
            if(this->pointer_ != 0)
            {
        #ifdef _WIN32
                ::UnmapViewOfFile(this->pointer_);
        #else
                munmap(this->pointer_, alloc_size_);
        #endif
                this->pointer_ = 0;
            }
        }

        std::size_t offset_, alloc_size_;
        FileHandle file_;

      private:
        Chunk & operator=(Chunk const &);
    };

    typedef MultiArray<N, SharedChunkHandle<N, T>  > ChunkStorage;
    typedef MultiArray<N, std::size_t>               OffsetStorage;
    typedef typename ChunkStorage::difference_type   shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;

    static std::size_t computeAllocSize(shape_type const & shape)
    {
        std::size_t size = prod(shape)*sizeof(T);
        std::size_t mask = mmap_alignment - 1;
        return (size + mask) & ~mask;
    }

    /** \brief Construct with given 'shape', 'chunk_shape' and 'options'.

        If the optional 'path' is given, the file is created in this directory.
        Otherwise (default), the path specified by the $TMP or $TEMP environment
        variables (in that order) is used.
    */
    explicit ChunkedArrayTmpFile(shape_type const & shape,
                                 shape_type const & chunk_shape=shape_type(),
                                 ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                                 std::string const & path = "")
    : ChunkedArray<N, T>(shape, chunk_shape, options)
    #ifndef VIGRA_NO_SPARSE_FILE
    , offset_array_(this->chunkArrayShape())
    #endif
    , file_size_()
    , file_capacity_()
    {
    #ifdef VIGRA_NO_SPARSE_FILE
        file_capacity_ = 4*prod(this->chunk_shape_)*sizeof(T);
    #else
        // compute offset in file
        typename OffsetStorage::iterator i = offset_array_.begin(),
                                         end = offset_array_.end();
        std::size_t size = 0;
        for(; i != end; ++i)
        {
            *i = size;
            size += computeAllocSize(this->chunkShape(i.point()));
        }
        file_capacity_ = size;
        this->overhead_bytes_ += offset_array_.size()*sizeof(std::size_t);
        // std::cerr << "    file size: " << size << "\n";
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
    }

    ~ChunkedArrayTmpFile()
    {
        typename ChunkStorage::iterator  i = this->handle_array_.begin(),
                                         end = this->handle_array_.end();
        for(; i != end; ++i)
        {
            if(i->pointer_)
                delete static_cast<Chunk*>(i->pointer_);
            i->pointer_ = 0;
        }
    #ifdef _WIN32
        ::CloseHandle(mappedFile_);
        ::CloseHandle(file_);
    #else
        ::close(file_);
    #endif
    }

    virtual pointer loadChunk(ChunkBase<N, T> ** p, shape_type const & index)
    {
        if(*p == 0)
        {
            shape_type shape = this->chunkShape(index);
            std::size_t chunk_size = computeAllocSize(shape);
        #ifdef VIGRA_NO_SPARSE_FILE
            std::size_t offset = file_size_;
            if(offset + chunk_size > file_capacity_)
            {
                file_capacity_ = max<std::size_t>(offset+chunk_size, file_capacity_ * 120 / 100); // extend file by 20%
                if(lseek(file_, file_capacity_-1, SEEK_SET) == -1)
                    throw std::runtime_error("ChunkedArrayTmpFile(): unable to reset file size.");
                if(write(file_, "0", 1) == -1)
                    throw std::runtime_error("ChunkedArrayTmpFile(): unable to resize file.");
            }
            file_size_ += chunk_size;
        #else
            std::size_t offset = offset_array_[index];
        #endif
            *p = new Chunk(shape, offset, chunk_size, mappedFile_);
            this->overhead_bytes_ += sizeof(Chunk);
        }
        return static_cast<Chunk*>(*p)->map();
    }

    virtual bool unloadChunk(ChunkBase<N, T> * chunk, bool /* destroy*/)
    {
        static_cast<Chunk *>(chunk)->unmap();
        return false; // never destroys the data
    }

    virtual std::string backend() const
    {
        return "ChunkedArrayTmpFile";
    }

    virtual std::size_t dataBytes(ChunkBase<N,T> * c) const
    {
        return c->pointer_ == 0
                 ? 0
                 : static_cast<Chunk*>(c)->alloc_size_;
    }

    virtual std::size_t overheadBytesPerChunk() const
    {
      #ifdef VIGRA_NO_SPARSE_FILE
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>);
      #else
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>) + sizeof(std::size_t);
      #endif
    }

  #ifndef VIGRA_NO_SPARSE_FILE
    OffsetStorage offset_array_;  // the array of chunks
  #endif
    FileHandle file_, mappedFile_;  // the file back-end
    std::size_t file_size_, file_capacity_;
};

template<unsigned int N, class U>
class ChunkIterator
: public MultiCoordinateIterator<N>
, private MultiArrayView<N, typename UnqualifiedType<U>::type>
{
  public:
    typedef typename UnqualifiedType<U>::type      T;
    typedef MultiCoordinateIterator<N>             base_type;
    typedef MultiArrayView<N, T>                   base_type2;

    typedef typename base_type::shape_type         shape_type;
    typedef typename base_type::difference_type    difference_type;
    typedef ChunkIterator                          iterator;
    typedef std::random_access_iterator_tag        iterator_category;

    typedef MultiArrayView<N, T>                   value_type;
    typedef MultiArrayView<N, T> &                 reference;
    typedef MultiArrayView<N, T> const &           const_reference;
    typedef MultiArrayView<N, T> *                 pointer;
    typedef MultiArrayView<N, T> const *           const_pointer;

    typedef typename IfBool<UnqualifiedType<U>::isConst,
                          ChunkedArrayBase<N, T> const,
                          ChunkedArrayBase<N, T> >::type array_type;
    typedef IteratorChunkHandle<N, T>        Chunk;


    ChunkIterator()
    : base_type()
    , base_type2()
    {}

    ChunkIterator(array_type * array,
                  shape_type const & start, shape_type const & end,
                  shape_type const & chunk_start, shape_type const & chunk_end,
                  shape_type const & chunk_shape)
    : base_type(chunk_start, chunk_end)
    , array_(array)
    , chunk_(chunk_start * chunk_shape)
    , start_(start - chunk_.offset_)
    , stop_(end - chunk_.offset_)
    , chunk_shape_(chunk_shape)
    {
        getChunk();
    }

    ChunkIterator(ChunkIterator const & rhs)
    : base_type(rhs)
    , base_type2(rhs)
    , array_(rhs.array_)
    , chunk_(rhs.chunk_)
    , start_(rhs.start_)
    , stop_(rhs.stop_)
    , chunk_shape_(rhs.chunk_shape_)
    {
        getChunk();
    }

    ChunkIterator & operator=(ChunkIterator const & rhs)
    {
        if(this != &rhs)
        {
            base_type::operator=(rhs);
            array_ = rhs.array_;
            chunk_ = rhs.chunk_;
            start_ = rhs.start_;
            stop_ = rhs.stop_;
            chunk_shape_ = rhs.chunk_shape_;
            getChunk();
        }
        return *this;
    }

    reference operator*()
    {
        return *this;
    }

    const_reference operator*() const
    {
        return *this;
    }

    pointer operator->()
    {
        return this;
    }

    const_pointer operator->() const
    {
        return this;
    }

    value_type operator[](MultiArrayIndex i) const
    {
        return *(ChunkIterator(*this) += i);
    }

    value_type operator[](const shape_type &coordOffset) const
    {
        return *(ChunkIterator(*this) += coordOffset);
    }

    void getChunk()
    {
        if(array_)
        {
            shape_type array_point = max(start_, this->point()*chunk_shape_),
                       upper_bound(SkipInitialization);
            this->m_ptr = array_->chunkForIterator(array_point, this->m_stride, upper_bound, &chunk_);
            this->m_shape = min(upper_bound, stop_) - array_point;
        }
    }

    shape_type chunkStart() const
    {
        return max(start_, this->point()*chunk_shape_) + chunk_.offset_;
    }

    shape_type chunkStop() const
    {
        return chunkStart() + this->m_shape;
    }

    ChunkIterator & operator++()
    {
        base_type::operator++();
        getChunk();
        return *this;
    }

    ChunkIterator operator++(int)
    {
        ChunkIterator res(*this);
        ++*this;
        return res;
    }

    ChunkIterator & operator+=(MultiArrayIndex i)
    {
        base_type::operator+=(i);
        getChunk();
        return *this;
    }

    ChunkIterator & operator+=(const shape_type &coordOffset)
    {
        base_type::operator+=(coordOffset);
        getChunk();
        return *this;
    }

    ChunkIterator & operator--()
    {
        base_type::operator--();
        getChunk();
        return *this;
    }

    ChunkIterator operator--(int)
    {
        ChunkIterator res(*this);
        --*this;
        return res;
    }

    ChunkIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    ChunkIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

    ChunkIterator getEndIterator() const
    {
        ChunkIterator res(*this);
        static_cast<base_type &>(res) = base_type::getEndIterator();
        res.getChunk();
        return res;
    }

    ChunkIterator operator+(MultiArrayIndex d) const
    {
        return ChunkIterator(*this) += d;
    }

    ChunkIterator operator-(MultiArrayIndex d) const
    {
        return ChunkIterator(*this) -= d;
    }

    ChunkIterator operator+(const shape_type &coordOffset) const
    {
        return ChunkIterator(*this) += coordOffset;
    }

    ChunkIterator operator-(const shape_type &coordOffset) const
    {
        return ChunkIterator(*this) -= coordOffset;
    }

    MultiArrayIndex operator-(const ChunkIterator & other) const
    {
        return base_type::operator-(other);
    }

#ifndef DOXYGEN  // doxygen doesn't understand this
    using base_type::operator==;
    using base_type::operator!=;
#endif
    using base_type::shape;

    array_type * array_;
    Chunk chunk_;
    shape_type start_, stop_, chunk_shape_, array_point_;
};

//@}

} // namespace vigra

#undef VIGRA_ASSERT_INSIDE

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HXX */
