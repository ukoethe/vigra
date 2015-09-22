/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2014 by Ullrich Koethe                            */
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

#ifndef MULTI_HANDLE_HXX
#define MULTI_HANDLE_HXX

#include "multi_fwd.hxx"
#include "metaprogramming.hxx"
#include "multi_shape.hxx"

namespace vigra {

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX=Handle::index>
struct CoupledHandleCast;

#ifndef _MSC_VER  // Visual Studio doesn't like these forward declarations
template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::reference
get(Handle & handle);

template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::const_reference
get(Handle const & handle);
#endif

/** \addtogroup MultiIteratorGroup
*/
//@{

  /**
     Handle class, used by CoupledScanOrderIterator as the value type to simultaneously itearate over multiple images.
  */
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

    template <class NEXT1>
    CoupledHandle(CoupledHandle<T, NEXT1> const & h, NEXT const & next)
    : base_type(next),
      pointer_(h.pointer_), 
      strides_(h.strides_)
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
    
    MultiArrayView<dimensions, T>
    arrayView() const
    {
        return MultiArrayView<dimensions, T>(this->shape(), strides(), ptr() - dot(this->point(), strides()));
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::reference
    get()
    {
        return vigra::get<TARGET_INDEX>(*this);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(*this);
    }


    void updatePtrAdresse(T & val){
        pointer_ = & val;
    }

    pointer pointer_;
    shape_type strides_;
};

    //  CoupledHandle holding the current coordinate
    // (always the end of a CoupledHandle chain)
template <int N>
class CoupledHandle<TinyVector<MultiArrayIndex, N>, void>
{
public:
    static const unsigned int index      = 0; // index of this member of the chain
    static const unsigned int dimensions = N;

    typedef typename MultiArrayShape<N>::type   value_type;
    typedef value_type const *                  pointer;
    typedef value_type const *                  const_pointer;
    typedef value_type const &                  reference;
    typedef value_type const &                  const_reference;
    typedef value_type                          shape_type;
    typedef CoupledHandle<value_type, void>     self_type;

    CoupledHandle()
    : point_(),
      shape_(),
      scanOrderIndex_()
    {}

    CoupledHandle(value_type const & shape)
    : point_(),
      shape_(shape),
      scanOrderIndex_()
    {}

    CoupledHandle(typename MultiArrayShape<N+1>::type const & shape)
    : point_(),
      shape_(shape.begin()),
      scanOrderIndex_()
    {}

    inline void incDim(int dim) 
    {
        ++point_[dim];
    }

    inline void decDim(int dim) 
    {
        --point_[dim];
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        point_[dim] += d;
    }

    inline void add(shape_type const & d) 
    {
        point_ += d;
    }
    
    template<int DIMENSION>
    inline void increment() 
    {
        ++point_[DIMENSION];
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        --point_[DIMENSION];
    }
    
    // TODO: test if making the above a default case of the this hurts performance
    template<int DIMENSION>
    inline void increment(MultiArrayIndex offset) 
    {
        point_[DIMENSION] += offset;
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex offset) 
    {
        point_[DIMENSION] -= offset;
    }

    void restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        point_ = shape_type();
        shape_ = end - start;
        scanOrderIndex_ = 0;
    }

    inline void incrementIndex() 
    {
        ++scanOrderIndex_;
    }
    
    inline void decrementIndex() 
    {
        --scanOrderIndex_;
    }
    
    inline void incrementIndex(MultiArrayIndex offset) 
    {
        scanOrderIndex_ += offset;
    }
    
    inline void decrementIndex(MultiArrayIndex offset) 
    {
        scanOrderIndex_ -= offset;
    }

    // access
    MultiArrayIndex scanOrderIndex() const
    {
        return scanOrderIndex_;
    }
    
    // access
    const_reference point() const
    {
        return point_;
    }
    
    // access
    const_reference shape() const
    {
        return shape_;
    }

    const_reference operator*() const
    {
        return point_;
    }

    const_pointer operator->() const
    {
        return &point_;
    }

    const_pointer ptr() const
    {
        return &point_;
    }

    unsigned int borderType() const
    {
        return detail::BorderTypeImpl<N>::exec(point_, shape_);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::reference
    get()
    {
        return vigra::get<TARGET_INDEX>(*this);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(*this);
    }

    value_type point_, shape_;
    MultiArrayIndex scanOrderIndex_;
};

    //  CoupledHandle for multi-band data
template <class T, class NEXT>
class CoupledHandle<Multiband<T>, NEXT>
: public NEXT
{
public:
    typedef NEXT                                  base_type;
    typedef CoupledHandle<Multiband<T>, NEXT>     self_type;
    
    static const unsigned int index =             NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef MultiArrayView<1, T, StridedArrayTag> value_type;
    typedef value_type *                          pointer;
    typedef value_type const *                    const_pointer;
    typedef value_type &                          reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;

    CoupledHandle()
    : base_type(),
      view_(), 
      strides_()
    {}
    
    template <class NEXT1>
    CoupledHandle(CoupledHandle<Multiband<T>, NEXT1> const & h, NEXT const & next)
    : base_type(next),
      view_(h.view_), 
      strides_(h.strides_)
    {}

    CoupledHandle(const_reference p, shape_type const & strides, NEXT const & next)
    : base_type(next),
      view_(p), 
      strides_(strides)
    {}

    template <class Stride>
    CoupledHandle(MultiArrayView<dimensions+1, Multiband<T>, Stride> const & v, NEXT const & next)
    : base_type(next),
      view_(v.bindInner(shape_type())), 
      strides_(v.bindOuter(0).stride())
    {
        vigra_precondition(v.bindOuter(0).shape() == this->shape(), "createCoupledIterator(): shape mismatch.");
    }

    inline void incDim(int dim) 
    {
        view_.unsafePtr() += strides_[dim];
        base_type::incDim(dim);
    }

    inline void decDim(int dim) 
    {
        view_.unsafePtr() -= strides_[dim];
        base_type::decDim(dim);
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        view_.unsafePtr() += d*strides_[dim];
        base_type::addDim(dim, d);
    }

    inline void add(shape_type const & d) 
    {
        view_.unsafePtr() += dot(d, strides_);
        base_type::add(d);
    }
    
    template<int DIMENSION>
    inline void increment() 
    {
        view_.unsafePtr() += strides_[DIMENSION];
        base_type::template increment<DIMENSION>();
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        view_.unsafePtr() -= strides_[DIMENSION];
        base_type::template decrement<DIMENSION>();
    }
    
    // TODO: test if making the above a default case of the this hurts performance
    template<int DIMENSION>
    inline void increment(MultiArrayIndex offset) 
    {
        view_.unsafePtr() += offset*strides_[DIMENSION];
        base_type::template increment<DIMENSION>(offset);
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex offset) 
    {
        view_.unsafePtr() -= offset*strides_[DIMENSION];
        base_type::template decrement<DIMENSION>(offset);
    }

    void restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        view_.unsafePtr() += dot(start, strides_);
        base_type::restrictToSubarray(start, end);
    }

    // ptr access
    reference operator*()
    {
        return view_;
    }

    const_reference operator*() const
    {
        return view_;
    }

    pointer operator->()
    {
        return &view_;
    }

    const_pointer operator->() const
    {
        return &view_;
    }

    pointer ptr()
    {
        return &view_;
    }

    const_pointer ptr() const
    {
        return &view_;
    }

    shape_type const & strides() const
    {
        return strides_;
    }
    
    MultiArrayView<dimensions+1, Multiband<T> >
    arrayView() const
    {
        typedef MultiArrayView<dimensions+1, T> View;
        typename View::difference_type vshape(SkipInitialization), vstride(SkipInitialization);
        vshape.template subarray<0, dimensions>() = this->shape();
        vstride.template subarray<0, dimensions>() = strides();
        vshape[dimensions] = view_.shape(0);
        vstride[dimensions] = view_.stride(0);
        return View(vshape, vstride, view_.data() - dot(this->point(), strides())).multiband();
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::reference
    get()
    {
        return vigra::get<TARGET_INDEX>(*this);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(*this);
    }

    value_type view_;
    shape_type strides_;
};

    //  helper class for CoupledHandle for CunkedArray
template <unsigned int N, class T>
class IteratorChunkHandle
{
  public:
    typedef ChunkedArray<N, T>             array_type;
    typedef typename MultiArrayShape<N>::type  shape_type;
    
    IteratorChunkHandle()
    : offset_(),
      chunk_(0)
    {}
    
    IteratorChunkHandle(shape_type const & offset)
    : offset_(offset),
      chunk_(0)
    {}
    
    IteratorChunkHandle(IteratorChunkHandle const & other)
    : offset_(other.offset_),
      chunk_(0)
    {}
    
    IteratorChunkHandle & operator=(IteratorChunkHandle const & other)
    {
        offset_ = other.offset_;
        chunk_ = 0;
        return *this;
    }
    
    shape_type offset_;
    SharedChunkHandle<N, T> * chunk_;
};

    /*  CoupledHandle for CunkedArray
    
        The handle must store a pointer to a chunk because the chunk knows 
        about memory menagement, and to an array view because it knows about
        subarrays and slices.
        
        Perhaps we can reduce this to a single pointer or otherwise reduce 
        the handle memory to make it faster?
    */
template <class U, class NEXT>
class CoupledHandle<ChunkedMemory<U>, NEXT>
: public NEXT,
  public IteratorChunkHandle<NEXT::dimensions, typename UnqualifiedType<U>::type>
{
public:
    typedef typename UnqualifiedType<U>::type     T;
    typedef NEXT                                  base_type;
    typedef IteratorChunkHandle<NEXT::dimensions, T>    base_type2;
    typedef CoupledHandle<ChunkedMemory<U>, NEXT> self_type;
    
    static const unsigned int index =             NEXT::index + 1;    // index of this member of the chain
    static const unsigned int dimensions =        NEXT::dimensions;

    typedef typename IfBool<UnqualifiedType<U>::isConst,
          ChunkedArrayBase<dimensions, T> const,
          ChunkedArrayBase<dimensions, T> >::type array_type;
    typedef detail::ChunkShape<dimensions, T>     chunk_shape;
    typedef T                                     value_type;
    typedef U *                                   pointer;
    typedef value_type const *                    const_pointer;
    typedef U &                                   reference;
    typedef value_type const &                    const_reference;
    typedef typename base_type::shape_type        shape_type;
    
    CoupledHandle()
    : base_type(),
      base_type2(),
      pointer_(), 
      strides_(),
      upper_bound_(),
      array_()
    {}

    CoupledHandle(CoupledHandle const & other)
    : base_type(other),
      base_type2(other),
      pointer_(other.pointer_),
      strides_(other.strides_),
      upper_bound_(other.upper_bound_),
      array_(other.array_)
    {
        if(array_)
            pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
    }

    CoupledHandle(array_type const & array, NEXT const & next)
    : base_type(next),
      base_type2(),
      pointer_(), 
      array_(const_cast<array_type*>(&array))
    {
        if(array_)
            pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
    }

    ~CoupledHandle()
    {
        // deref the present chunk
        if(array_)
            array_->unrefChunk(this);
    }

    CoupledHandle & operator=(CoupledHandle const & other)
    {
        if(this != &other)
        {
            // deref the present chunk
            if(array_)
                array_->unrefChunk(this);
            base_type::operator=(other);
            base_type2::operator=(other);
            array_ = other.array_;
            if(array_)
            {
                pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
            }
            else
            {
                pointer_ = other.pointer_;
                strides_ = other.strides_;
                upper_bound_ = other.upper_bound_;
            }
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
                pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
        }
    }

    inline void decDim(int dim) 
    {
        base_type::decDim(dim);
        pointer_ -= strides_[dim];
        if(point()[dim] < upper_bound_[dim] - array_->chunk_shape_[dim])
        {
            // if(point()[dim] >= 0)
                pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
        }
    }

    inline void addDim(int dim, MultiArrayIndex d) 
    {
        base_type::addDim(dim, d);
        if(point()[dim] < shape()[dim] && point()[dim] >= 0)
            pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
    }

    inline void add(shape_type const & d) 
    {
        base_type::add(d);
        pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
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
                pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
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
                pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
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
        pointer_ = array_->chunkForIterator(point(), strides_, upper_bound_, this);
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
    
    array_type const &
    arrayView() const
    {
        return *array_;
    }

    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::reference
    get()
    {
        return vigra::get<TARGET_INDEX>(*this);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, CoupledHandle, index>::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(*this);
    }

    pointer pointer_;
    shape_type strides_, upper_bound_;
    array_type * array_;
};

    // meta-programming helper classes to implement 'get<INDEX>(CoupledHandle)'
template <unsigned TARGET_INDEX>
struct Error__CoupledHandle_index_out_of_range;

namespace detail {

template <unsigned TARGET_INDEX, class Handle, bool isValid, unsigned int INDEX=Handle::index>
struct CoupledHandleCastImpl
{
    typedef typename CoupledHandleCastImpl<TARGET_INDEX, typename Handle::base_type, isValid>::type type;
    typedef typename type::value_type value_type;
    typedef typename type::reference reference;
    typedef typename type::const_reference const_reference;
};

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX>
struct CoupledHandleCastImpl<TARGET_INDEX, Handle, false, INDEX>
{
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> type;
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> value_type;
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> reference;
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> const_reference;
};

template <unsigned TARGET_INDEX, class Handle>
struct CoupledHandleCastImpl<TARGET_INDEX, Handle, true, TARGET_INDEX>
{
    typedef Handle type;
    typedef typename type::value_type value_type;
    typedef typename type::reference reference;
    typedef typename type::const_reference const_reference;
};

} // namespace detail

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX>
struct CoupledHandleCast
: public detail::CoupledHandleCastImpl<TARGET_INDEX, Handle, (TARGET_INDEX <= INDEX), INDEX>
{};

template <unsigned int TARGET_INDEX, class Handle>
inline
typename CoupledHandleCast<TARGET_INDEX, Handle>::type &
cast(Handle & handle)
{
    return handle;
}

template <unsigned int TARGET_INDEX, class Handle>
inline
typename CoupledHandleCast<TARGET_INDEX, Handle>::type const &
cast(Handle const & handle)
{
    return handle;
}

  /** Returns reference to the element in the band of the handle with index TARGET_INDEX.
   */
template <unsigned int TARGET_INDEX, class Handle>
inline
typename CoupledHandleCast<TARGET_INDEX, Handle>::reference
get(Handle & handle)
{
    return *cast<TARGET_INDEX>(handle);
}

  /** Returns a constant reference to the element in the band of the handle with index TARGET_INDEX.
   */
template <unsigned int TARGET_INDEX, class Handle>
inline
typename CoupledHandleCast<TARGET_INDEX, Handle>::const_reference
get(Handle const & handle)
{
    return *cast<TARGET_INDEX>(handle);
}

    // meta-programming helper classes to infer the type of 
    // a CoupledHandle for a set of arrays
template <unsigned int N, class List>
struct ComposeCoupledHandle;

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

    // meta-programming helper classes to implement 'zip(iterator1, iterator2)'
template <class A, class B>
struct ZipCoupledHandles;

template <class A, class Head, class Tail>
struct ZipCoupledHandles<A, CoupledHandle<Head, Tail> >
{
    typedef typename ZipCoupledHandles<A, Tail>::type Next;
    typedef CoupledHandle<Head, Next> type;
    
    static type construct(A const & a, CoupledHandle<Head, Tail> const & h)
    {
        return type(h, ZipCoupledHandles<A, Tail>::construct(a, (Tail const &)h));
    }
};

template <class A, class Shape>
struct ZipCoupledHandles<A, CoupledHandle<Shape, void> >
{
    typedef A type;
    
    static type construct(A const & a, CoupledHandle<Shape, void> const &)
    {
        return a;
    }
};

    // allow an iterator that uses CoupledHandle to specialize its 
    // dereferencing functions, such that
    //    '*iter' returns a referenc to the current point if 
    //            the handle is just a coordinate handle
    //    '*iter' returns a reference to the current data element 
    //            if the handle referes to just one array
    //    '*iter' returns a reference to the handle itself if it refers to
    //            several arrays simultaneously (i.e. is actualy a coupled handle)
template <class Handle, unsigned int INDEX=Handle::index>
struct CoupledHandleTraits
{
    typedef Handle value_type;
    typedef Handle & reference;
    typedef Handle const & const_reference;
    typedef Handle * pointer;
    typedef Handle const * const_pointer;
    
    static reference dereference(Handle & h)
    {
        return h;
    }
    
    static const_reference dereference(Handle const & h)
    {
        return h;
    }
};

template <class Handle>
struct CoupledHandleTraits<Handle, 0>
{
    typedef typename Handle::value_type value_type;
    typedef typename Handle::reference reference;
    typedef typename Handle::const_reference const_reference;
    typedef typename Handle::pointer pointer;
    typedef typename Handle::const_pointer const_pointer;
    
    static reference dereference(Handle & h)
    {
        return *h;
    }
    
    static const_reference dereference(Handle const & h)
    {
        return *h;
    }
};

template <class Handle>
struct CoupledHandleTraits<Handle, 1>
{
    typedef typename Handle::value_type value_type;
    typedef typename Handle::reference reference;
    typedef typename Handle::const_reference const_reference;
    typedef typename Handle::pointer pointer;
    typedef typename Handle::const_pointer const_pointer;
    
    static reference dereference(Handle & h)
    {
        return *h;
    }
    
    static const_reference dereference(Handle const & h)
    {
        return *h;
    }
};


//@}

} // namespace vigra

#endif /* MULTI_HANDLE_HXX */
