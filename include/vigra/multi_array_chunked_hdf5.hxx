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

#ifndef VIGRA_MULTI_ARRAY_CHUNKED_HDF5_HXX
#define VIGRA_MULTI_ARRAY_CHUNKED_HDF5_HXX

#include <queue>

#include "multi_array_chunked.hxx"
#include "hdf5impex.hxx"

// Bounds checking Macro used if VIGRA_CHECK_BOUNDS is defined.
#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(this->isInside(diff), "Index out of bounds")
#else
#define VIGRA_ASSERT_INSIDE(diff)
#endif

namespace vigra {

/** \addtogroup ChunkedArrayClasses
*/
//@{

/** Implement ChunkedArray as a chunked dataset in an HDF5 file.

    <b>\#include</b> \<vigra/multi_array_chunked_hdf5.hxx\> <br/>
    Namespace: vigra

    This uses the native chunking and compression functionality provided by the
    HDF5 library. Note: This file must only be included when the HDF5 headers
    and libraries are installed on the system.
*/
template <unsigned int N, class T, class Alloc = std::allocator<T> >
class ChunkedArrayHDF5
: public ChunkedArray<N, T>
{
    /* REMARKS
    Alternatives are:
    * Back chunks by HDF5 chunks, possibly using on-the-fly compression. This
      is in particular useful for existing HDF5 files.
    * Back chunks by HDF5 datasets. This can be combined with compression
      (both explicit and on-the-fly) or with memory mapping (using the
      function H5Dget_offset() to get the offset from the beginning of the file).
    */

  public:

    class Chunk
    : public ChunkBase<N, T>
    {
      public:
        typedef typename MultiArrayShape<N>::type  shape_type;
        typedef T value_type;
        typedef value_type * pointer;
        typedef value_type & reference;

        Chunk(shape_type const & shape, shape_type const & start,
              ChunkedArrayHDF5 * array, Alloc const & alloc)
        : ChunkBase<N, T>(detail::defaultStride(shape))
        , shape_(shape)
        , start_(start)
        , array_(array)
        , alloc_(alloc)
        {}

        ~Chunk()
        {
            write();
        }

        std::size_t size() const
        {
            return prod(shape_);
        }

        void write(bool deallocate = true)
        {
            if(this->pointer_ != 0)
            {
                if(!array_->file_.isReadOnly())
                {
                    herr_t status = array_->file_.writeBlock(array_->dataset_, start_,
                                          MultiArrayView<N, T>(shape_, this->strides_, this->pointer_));
                    vigra_postcondition(status >= 0,
                        "ChunkedArrayHDF5: write to dataset failed.");
                }
                if(deallocate)
                {
                    alloc_.deallocate(this->pointer_, this->size());
                    this->pointer_ = 0;
                }
            }
        }

        pointer read()
        {
            if(this->pointer_ == 0)
            {
                this->pointer_ = alloc_.allocate(this->size());
                herr_t status = array_->file_.readBlock(array_->dataset_, start_, shape_,
                                     MultiArrayView<N, T>(shape_, this->strides_, this->pointer_));
                vigra_postcondition(status >= 0,
                    "ChunkedArrayHDF5: read from dataset failed.");
            }
            return this->pointer_;
        }

        shape_type shape_, start_;
        ChunkedArrayHDF5 * array_;
        Alloc alloc_;

      private:
        Chunk & operator=(Chunk const &);
    };

    typedef ChunkedArray<N, T> base_type;
    typedef MultiArray<N, SharedChunkHandle<N, T> > ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;

    /** \brief Construct with given 'shape', 'chunk_shape' and 'options',
        using 'alloc' to manage the in-memory version of the data..

        The data are placed in 'file' at the internal path 'dataset'. Argument
        'mode' must be one of the following:
        <ul>
        <li>HDF5File::New: Create new dataset, possibly deleting any existing content.
                           It is an error to request this mode when the entire
                           'file' is read-only.
        <li>HDF5File::Replace: Same as New.
        <li>HDF5File::ReadWrite: Open the dataset for reading and writing. Create
                                 the datset if it doesn't exist. It is an error
                                 to request this mode when 'file' is read-only.
        <li>HDF5File::ReadOnly: Open the dataset for reading. It is an error to
                                request this mode when the dataset doesn't exist.
        <li>HDF5File::Default: Resolves to ReadOnly when the dataset exists, and
                               to New otherwise.
        </ul>
        The supported compression algorithms are:
        <ul>
        <li>ZLIB_FAST: Fast compression using 'zlib' (slower than LZ4, but higher compression).
        <li>ZLIB_BEST: Best compression using 'zlib', slow.
        <li>ZLIB_NONE: Use 'zlib' format without compression.
        <li>DEFAULT_COMPRESSION: Same as ZLIB_FAST.
        </ul>
    */
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     HDF5File::OpenMode mode,
                     shape_type const & shape,
                     shape_type const & chunk_shape=shape_type(),
                     ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                     Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape, chunk_shape, options),
      file_(file),
      dataset_name_(dataset),
      dataset_(),
      compression_(options.compression_method),
      alloc_(alloc)
    {
        init(mode);
    }

    /** \brief Construct for an already existing dataset with given 'options',
        using 'alloc' to manage the in-memory version of the data.

        The data must be located in 'file' at the internal path 'dataset'. The
        array's shape and chunk_shape are read from the file. It is an error
        to use this constructor when 'dataset' doesn't exist.

        Argument 'mode' must be one of the following:
        <ul>
        <li>HDF5File::ReadWrite: Open the dataset for reading and writing. It is an error
                                 to request this mode when 'file' is read-only.
        <li>HDF5File::ReadOnly: Open the dataset for reading (default).
        <li>HDF5File::Default: Same as ReadOnly.
        </ul>
        The supported compression algorithms are:
        <ul>
        <li>ZLIB_FAST: Fast compression using 'zlib' (slower than LZ4, but higher compression).
        <li>ZLIB_BEST: Best compression using 'zlib', slow.
        <li>ZLIB_NONE: Use 'zlib' format without compression.
        <li>DEFAULT_COMPRESSION: Same as ZLIB_FAST.
        </ul>
    */
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     HDF5File::OpenMode mode = HDF5File::ReadOnly,
                     ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                     Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape_type(), shape_type(), options),
      file_(file),
      dataset_name_(dataset),
      dataset_(),
      compression_(options.compression_method),
      alloc_(alloc)
    {
        init(mode);
    }

    void init(HDF5File::OpenMode mode)
    {
        bool exists = file_.existsDataset(dataset_name_);

        if(mode == HDF5File::Replace)
        {
            mode = HDF5File::New;
        }
        else if(mode == HDF5File::Default)
        {
            if(exists)
                mode = HDF5File::ReadOnly;
            else
                mode = HDF5File::New;
        }

        if(mode == HDF5File::ReadOnly)
            file_.setReadOnly();
        else
            vigra_precondition(!file_.isReadOnly(),
                "ChunkedArrayHDF5(): 'mode' is incompatible with read-only file.");

        vigra_precondition(exists || !file_.isReadOnly(),
            "ChunkedArrayHDF5(): dataset does not exist, but file is read-only.");

        if(!exists || mode == HDF5File::New)
        {
            // FIXME: set rdcc_nbytes to 0 (disable cache, because we don't
            //        need two caches
            // H5Pset_chunk_cache (dapl, rdcc_nslots, rdcc_nbytes, rdcc_w0);
            // Chunk cache size (rdcc_nbytes) should be large
            // enough to hold all the chunks in a selection
            // • If this is not possible, it may be best to disable chunk
            // caching altogether (set rdcc_nbytes to 0)
            // • rdcc_slots should be a prime number that is at
            // least 10 to 100 times the number of chunks that can fit
            // into rdcc_nbytes
            // • rdcc_w0 should be set to 1 if chunks that have been
            // fully read/written will never be read/written again
            //
            // the above may be WRONG in general - it may only apply if the
            // chunk size in the file matches the chunk size in the CachedArray.
            // Otherwise, make sure that the file cache can hold at least as many
            // chunks as are needed for a single array chunk.
            if(compression_ == DEFAULT_COMPRESSION)
                compression_ = ZLIB_FAST;
            vigra_precondition(compression_ != LZ4,
                "ChunkedArrayHDF5(): HDF5 does not support LZ4 compression.");

            vigra_precondition(this->size() > 0,
                "ChunkedArrayHDF5(): invalid shape.");
            typename detail::HDF5TypeTraits<T>::value_type init(this->fill_scalar_);
            dataset_ = file_.createDataset<N, T>(dataset_name_,
                                                 this->shape_,
                                                 init,
                                                 this->chunk_shape_,
                                                 compression_);
        }
        else
        {
            dataset_ = file_.getDatasetHandleShared(dataset_name_);

            // check shape
            ArrayVector<hsize_t> fileShape(file_.getDatasetShape(dataset_name_));
            typedef detail::HDF5TypeTraits<T> TypeTraits;
            if(TypeTraits::numberOfBands() > 1)
            {
                vigra_precondition(fileShape.size() == N+1,
                    "ChunkedArrayHDF5(file, dataset): dataset has wrong dimension.");
                vigra_precondition(fileShape[0] == TypeTraits::numberOfBands(),
                    "ChunkedArrayHDF5(file, dataset): dataset has wrong number of bands.");
                shape_type shape(fileShape.begin()+1);
                if(this->size() > 0)
                {
                    vigra_precondition(shape == this->shape_,
                        "ChunkedArrayHDF5(file, dataset, shape): shape mismatch between dataset and shape argument.");
                }
                else
                {
                    this->shape_ = shape;
                }
            }
            else
            {
                vigra_precondition(fileShape.size() == N,
                    "ChunkedArrayHDF5(file, dataset): dataset has wrong dimension.");
                shape_type shape(fileShape.begin());
                if(this->size() > 0)
                {
                    vigra_precondition(shape == this->shape_,
                        "ChunkedArrayHDF5(file, dataset, shape): shape mismatch between dataset and shape argument.");
                }
                else
                {
                    this->shape_ = shape;
                    ChunkStorage(detail::computeChunkArrayShape(shape, this->bits_, this->mask_)).swap(this->handle_array_);
                }
            }
            typename ChunkStorage::iterator i   = this->handle_array_.begin(),
                                            end = this->handle_array_.end();
            for(; i != end; ++i)
            {
                i->chunk_state_.store(base_type::chunk_asleep);
            }
        }
    }

    ~ChunkedArrayHDF5()
    {
        closeImpl(true);
    }

    void close()
    {
        closeImpl(false);
    }

    void closeImpl(bool force_destroy)
    {
        flushToDiskImpl(true, force_destroy);
        file_.close();
    }

    void flushToDisk()
    {
        flushToDiskImpl(false, false);
    }

    void flushToDiskImpl(bool destroy, bool force_destroy)
    {
        if(file_.isReadOnly())
            return;

        threading::lock_guard<threading::mutex> guard(*this->chunk_lock_);
        typename ChunkStorage::iterator i   = this->handle_array_.begin(),
                                        end = this->handle_array_.end();
        if(destroy && !force_destroy)
        {
            for(; i != end; ++i)
            {
                vigra_precondition(i->chunk_state_.load() <= 0,
                    "ChunkedArrayHDF5::close(): cannot close file because there are active chunks.");
            }
            i   = this->handle_array_.begin();
        }
        for(; i != end; ++i)
        {
            Chunk * chunk = static_cast<Chunk*>(i->pointer_);
            if(!chunk)
                continue;
            if(destroy)
            {
                delete chunk;
                i->pointer_ = 0;
            }
            else
            {
                chunk->write(false);
            }
        }
        file_.flushToDisk();
    }

    virtual bool isReadOnly() const
    {
        return file_.isReadOnly();
    }

    virtual pointer loadChunk(ChunkBase<N, T> ** p, shape_type const & index)
    {
        vigra_precondition(file_.isOpen(),
            "ChunkedArrayHDF5::loadChunk(): file was already closed.");
        if(*p == 0)
        {
            *p = new Chunk(this->chunkShape(index), index*this->chunk_shape_, this, alloc_);
            this->overhead_bytes_ += sizeof(Chunk);
        }
        return static_cast<Chunk *>(*p)->read();
    }

    virtual bool unloadChunk(ChunkBase<N, T> * chunk, bool /* destroy */)
    {
        if(!file_.isOpen())
            return true;
        static_cast<Chunk *>(chunk)->write();
        return false;
    }

    virtual std::string backend() const
    {
        return "ChunkedArrayHDF5<'" + file_.filename() + "/" + dataset_name_ + "'>";
    }

    virtual std::size_t dataBytes(ChunkBase<N,T> * c) const
    {
        return c->pointer_ == 0
                 ? 0
                 : static_cast<Chunk*>(c)->size()*sizeof(T);
    }

    virtual std::size_t overheadBytesPerChunk() const
    {
        return sizeof(Chunk) + sizeof(SharedChunkHandle<N, T>);
    }

    std::string fileName() const
    {
        return file_.filename();
    }

    std::string datasetName() const
    {
        return dataset_name_;
    }

    HDF5File file_;
    std::string dataset_name_;
    HDF5HandleShared dataset_;
    CompressionMethod compression_;
    Alloc alloc_;
};

//@}

} // namespace vigra

#undef VIGRA_ASSERT_INSIDE

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HDF5_HXX */
