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
    
    typedef MultiArray<N, SharedChunkHandle<N, T> > ChunkStorage;
    typedef typename ChunkStorage::difference_type  shape_type;
    typedef T value_type;
    typedef value_type * pointer;
    typedef value_type & reference;
    
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
    
    ChunkedArrayHDF5(HDF5File const & file, std::string const & dataset,
                     ChunkedArrayOptions const & options = ChunkedArrayOptions(),
                     Alloc const & alloc = Alloc())
    : ChunkedArray<N, T>(shape_type(), shape_type(), options),
      file_(file),
      dataset_name_(dataset),
      dataset_(),
      compression_(options.compression_method),
      alloc_(alloc)
    {
        init(HDF5File::OpenReadOnly);
    }
    
    void init(HDF5File::OpenMode mode)
    {
        if(mode == HDF5File::OpenReadOnly)
            file_.setReadOnly();
        else
            vigra_precondition(!file_.isReadOnly(),
                "ChunkedArrayHDF5(): 'mode' is incompatible with read-only file.");
                
        bool exists = file_.existsDataset(dataset_name_);
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
            if(compression_ == DEFAULT_COMPRESSION)
                compression_ = ZLIB_FAST;
            vigra_precondition(compression_ != LZ4,
                "ChunkedArrayHDF5(): HDF5 does not support LZ4 compression.");
            
            vigra_precondition(this->size() > 0,
                "ChunkedArrayHDF5(): invalid shape.");
            typename detail::HDF5TypeTraits<T>::value_type init(fill_scalar_);
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
                }
            }
        }
    }
    
    ~ChunkedArrayHDF5()
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
    
    void flushToDisk()
    {
        if(file_.isReadOnly())
            return;
            
        typename ChunkStorage::iterator i   = this->handle_array_.begin(), 
                                        end = this->handle_array_.end();
        for(; i != end; ++i)
        {
            Chunk * chunk = static_cast<Chunk*>(i->pointer_);
            if(chunk)
                chunk->write(false);
        }
        file_.flushToDisk();
    }
    
    virtual bool isReadOnly() const
    {
        return file_.isReadOnly();
    }
    
    virtual pointer loadChunk(ChunkBase<N, T> ** p, shape_type const & index)
    {
        if(*p == 0)
        {
            *p = new Chunk(this->chunkShape(index), index*this->chunk_shape_, this, alloc_);
            this->overhead_bytes_ += sizeof(Chunk);
        }
        return static_cast<Chunk *>(*p)->read();
    }
    
    virtual bool unloadChunk(ChunkBase<N, T> * chunk, bool /* destroy */)
    {
        static_cast<Chunk *>(chunk)->write();
        return false; // never destroys the data
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
    
    std::string const & datasetName() const
    {
        return dataset_name_;
    }
    
    HDF5File file_;
    std::string dataset_name_;
    HDF5HandleShared dataset_;
    CompressionMethod compression_;
    Alloc alloc_;
};

} // namespace vigra

#undef VIGRA_ASSERT_INSIDE

#endif /* VIGRA_MULTI_ARRAY_CHUNKED_HDF5_HXX */
