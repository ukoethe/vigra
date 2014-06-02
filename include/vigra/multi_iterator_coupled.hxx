/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 by Stefan Schmidt and Ullrich Koethe         */
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

#ifndef MULTI_ITERATOR_COUPLED_HXX
#define MULTI_ITERATOR_COUPLED_HXX

#include "metaprogramming.hxx"
#include "multi_shape.hxx"

namespace vigra {

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


template <int N>
class CoupledHandle<TinyVector<MultiArrayIndex, N>, void>
{
public:
    static const int index = 0;                   // index of this member of the chain
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

    value_type point_, shape_;
    MultiArrayIndex scanOrderIndex_;
};

template <unsigned int N, class T, class StrideTag>
class MultiArrayView<N, Multiband<T>, StrideTag>
: public MultiArrayView<N, T, StrideTag>
{
  public:
    MultiArrayView(MultiArrayView<N, T, StrideTag> const & v)
    : MultiArrayView<N, T, StrideTag>(v)
    {}
};

template <class T, class NEXT>
class CoupledHandle<Multiband<T>, NEXT>
: public NEXT
{
public:
    typedef NEXT                                  base_type;
    typedef CoupledHandle<Multiband<T>, NEXT>     self_type;
    
    static const int index =                      NEXT::index + 1;    // index of this member of the chain
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

    value_type view_;
    shape_type strides_;
};

template <unsigned TARGET_INDEX>
struct Error__CoupledHandle_index_out_of_range;

namespace detail {

template <unsigned TARGET_INDEX, class Handle, bool isValid, unsigned int INDEX=Handle::index>
struct CoupledHandleCastImpl
{
    typedef typename CoupledHandleCastImpl<TARGET_INDEX, typename Handle::base_type, isValid>::type type;
    typedef typename type::reference reference;
    typedef typename type::const_reference const_reference;
};

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX>
struct CoupledHandleCastImpl<TARGET_INDEX, Handle, false, INDEX>
{
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> type;
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> reference;
    typedef Error__CoupledHandle_index_out_of_range<TARGET_INDEX> const_reference;
};

template <unsigned TARGET_INDEX, class Handle>
struct CoupledHandleCastImpl<TARGET_INDEX, Handle, true, TARGET_INDEX>
{
    typedef Handle type;
    typedef typename type::reference reference;
    typedef typename type::const_reference const_reference;
};

} // namespace detail

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX=Handle::index>
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

/********************************************************/
/*                                                      */
/*               CoupledScanOrderIterator<N>            */
/*                                                      */
/********************************************************/

template <class Iterator>
class CoupledDimensionProxy
: public Iterator
{
  public:
    typedef typename Iterator::value_type          value_type;
    typedef typename Iterator::difference_type     difference_type;
    typedef typename Iterator::reference           reference;
    typedef typename Iterator::const_reference     const_reference;
    typedef typename Iterator::pointer             pointer;
    typedef CoupledDimensionProxy                  iterator;
    typedef std::random_access_iterator_tag        iterator_category;
    
    static const int dimension = Iterator::dimension;
 
    CoupledDimensionProxy & operator++()
    {
        this->incDim(dimension);
        return *this;
    }
    
    CoupledDimensionProxy operator++(int)
    {
        CoupledDimensionProxy ret(*this);
        this->incDim(dimension);
        return ret;
    }
    
    CoupledDimensionProxy & operator--()
    {
        this->decDim(dimension);
        return *this;
    }
    
    CoupledDimensionProxy operator--(int)
    {
        CoupledDimensionProxy ret(*this);
        this->decDim(dimension);
        return ret;
    }
    
    CoupledDimensionProxy & operator+=(MultiArrayIndex d)
    {
        this->addDim(dimension, d);
        return *this;
    }
    
    CoupledDimensionProxy & operator-=(MultiArrayIndex d)
    {
        this->addDim(dimension, -d);
        return *this;
    }
    
    value_type operator[](MultiArrayIndex d) const
    {
        *(CoupledDimensionProxy(*this) += d);
    }
    
    CoupledDimensionProxy & operator=(MultiArrayIndex d)
    {
        this->setDim(dimension, d);
        return *this;
    }
    
    bool operator==(MultiArrayIndex d) const
    {
        return this->point(dimension) == d;
    }
    
    bool operator!=(MultiArrayIndex d) const
    {
        return this->point(dimension) != d;
    }
    
    bool operator<(MultiArrayIndex d) const
    {
        return this->point(dimension) < d;
    }
    
    bool operator<=(MultiArrayIndex d) const
    {
        return this->point(dimension) <= d;
    }
    
    bool operator>(MultiArrayIndex d) const
    {
        return this->point(dimension) > d;
    }
    
    bool operator>=(MultiArrayIndex d) const
    {
        return this->point(dimension) >= d;
    }
};

/** \brief Iterate over multiple images simultaneously in scan order. 

    The value type of this iterator is an instance of the handle class CoupledHandle. This allows to iterate over multiple arrays simultaneously. The coordinates can be accessed as a special band (index 0) in the handle. The scan-order is defined such that dimensions are iterated from front to back (first to last).
    
    Instances of this class are usually constructed by calling createCoupledIterator() .

    To get the type of a CoupledScanOrderIterator for arrays of a certain dimension and element types use CoupledIteratorType::type.

    The iterator supports all functions listed in the STL documentation for
        <a href="http://www.sgi.com/tech/stl/RandomAccessIterator.html">Random Access Iterators</a>.

    Example of use:
    \code
    using namespace vigra;
    MultiArray<2, double> image1(Shape2(5, 5));
    MultiArray<2, double> image2(Shape2(5, 5));
    // fill image with data ...
    
    typedef CoupledIteratorType<2, double, double>::type Iterator; // the type of the CoupledScanOrderIterator
    
    Iterator start = createCoupledIterator(image1, image2); // create coupled iterator for simultaneous iteration over image1, image2 and their coordinates
    Iterator end = start.getEndIterator();
    
    for (Iterator it = start; it < end; ++it) {
      std::cout << "coordinates: " << it.get<0>() << std::endl;
      std::cout << "image1: " << it.get<1>() << std::endl;
      std::cout << "image2: " << it.get<2>() << std::endl;
    }
    
    //random access:
    Iterator::value_type handle = start[15];
    std::cout << "image1: " << get<1>(handle) << std::endl;
    \endcode
    
    <b>\#include</b> \<vigra/multi_iterator_coupled.hxx\> <br/>
    Namespace: vigra
*/

template <unsigned int N,
          class HANDLES=CoupledHandle<TinyVector<MultiArrayIndex, N>, void>,
          int DIMENSION = N-1>
class CoupledScanOrderIterator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: public CoupledScanOrderIterator<N, HANDLES, DIMENSION-1>
#endif
{
    typedef CoupledScanOrderIterator<N, HANDLES, DIMENSION-1> base_type;

  public:
     static const int dimension = DIMENSION;
 
    typedef MultiArrayIndex                   difference_type;
    typedef CoupledScanOrderIterator          iterator;
    typedef std::random_access_iterator_tag   iterator_category;

    typedef typename base_type::value_type      value_type;

#ifdef DOXYGEN
  /** The type of the CoupledHandle.
   */
    typedef HANDLES value_type;
#endif

    typedef typename base_type::shape_type      shape_type;
    typedef typename base_type::reference       reference;
    typedef typename base_type::const_reference const_reference; // FIXME: do we need both?
    typedef typename base_type::pointer         pointer;
    typedef CoupledDimensionProxy<iterator>     dimension_proxy;

    explicit CoupledScanOrderIterator(value_type const & handles = value_type())
    : base_type(handles)
    {}

    CoupledScanOrderIterator & operator++()
    {
        base_type::operator++();
        if(this->point()[dimension-1] == this->shape()[dimension-1])
        {
            resetAndIncrement();
        }
        return *this;
    }

    CoupledScanOrderIterator operator++(int)
    {
        CoupledScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    CoupledScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        base_type::operator+=(i);
        return *this;
    }

    CoupledScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        base_type::operator+=(coordOffset);
        return *this;
    }

    CoupledScanOrderIterator & operator--()
    {
        base_type::operator--();
        if(this->point()[dimension-1] == -1)
        {
            resetAndDecrement();
        }
        return *this;
    }

    CoupledScanOrderIterator operator--(int)
    {
        CoupledScanOrderIterator res(*this);
        --*this;
        return res;
    }

    CoupledScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    CoupledScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

        /** Returns CoupledScanOrderIterator pointing beyond the last element.
        */
    CoupledScanOrderIterator getEndIterator() const
    {
        return operator+(prod(this->shape()) - this->scanOrderIndex());
    }
    
    CoupledScanOrderIterator operator+(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) += d;
    }

    CoupledScanOrderIterator operator-(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) -= d;
    }

    CoupledScanOrderIterator operator+(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) += coordOffset;
    }

    CoupledScanOrderIterator operator-(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) -= coordOffset;
    }

    MultiArrayIndex operator-(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator-(r);
    }

    CoupledScanOrderIterator & 
    restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        base_type::restrictToSubarray(start, end);
        return *this;
    }

#ifdef DOXYGEN
  
        /** Returns reference to the element in the band with index TARGET_INDEX.
        */
    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::type::reference
    get();

        /** Returns constant reference to the element in the band with index TARGET_INDEX.
        */
    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::type::const_reference
    get() const;
    
#endif

  protected:
        // placing these functions out-of-line prevents MSVC
        // from stupid optimizations
    void resetAndIncrement();
    void resetAndDecrement();

    void reset()
    {
        this->handles_.template decrement<dimension>(this->shape()[dimension]);
    }

    void inverseReset()
    {
        this->handles_.template increment<dimension>(this->shape()[dimension]);
    }
};

template <unsigned int N, class HANDLES, int DIMENSION>
void CoupledScanOrderIterator<N, HANDLES, DIMENSION>::resetAndIncrement()
{
    base_type::reset();
    this->handles_.template increment<dimension>();
}

template <unsigned int N, class HANDLES, int DIMENSION>
void CoupledScanOrderIterator<N, HANDLES, DIMENSION>::resetAndDecrement()
{
    base_type::inverseReset();
    this->handles_.template decrement<dimension>();
}

template <unsigned int N, class HANDLES>
class CoupledScanOrderIterator<N, HANDLES, 0>
{
  public:

    static const int dimension = 0;

    typedef CoupledScanOrderIterator<N, HANDLES, 0>  self_type;
    typedef HANDLES                                  value_type;
    typedef MultiArrayIndex                          difference_type;
    typedef value_type &                             reference;
    typedef value_type const &                       const_reference; 
    typedef value_type *                             pointer;
    typedef typename MultiArrayShape<N>::type        shape_type;
    typedef CoupledScanOrderIterator                 iterator;
    typedef std::random_access_iterator_tag          iterator_category;
    typedef CoupledDimensionProxy<iterator>          dimension_proxy;

    explicit CoupledScanOrderIterator(value_type const & handles = value_type())
    : handles_(handles),
      strides_(detail::defaultStride(handles_.shape()))
    {}
    
    template <unsigned int DIM>
    typename CoupledScanOrderIterator<N, HANDLES, DIM>::dimension_proxy & 
    dim()
    {
        typedef CoupledScanOrderIterator<N, HANDLES, DIM> Iter;
        typedef typename Iter::dimension_proxy Proxy;
        return static_cast<Proxy &>(static_cast<Iter &>(*this));
    }
    
    template <unsigned int DIM>
    typename CoupledScanOrderIterator<N, HANDLES, DIM>::dimension_proxy const & 
    dim() const
    {
        typedef CoupledScanOrderIterator<N, HANDLES, DIM> Iter;
        typedef typename Iter::dimension_proxy Proxy;
        return static_cast<Proxy const &>(static_cast<Iter const &>(*this));
    }

    void incDim(int dim)
    {
        handles_.incDim(dim);
        handles_.incrementIndex(strides_[dim]);
    }

    void decDim(int dim)
    {
        handles_.decDim(dim);
        handles_.decrementIndex(strides_[dim]);
    }

    void addDim(int dim, MultiArrayIndex d)
    {
        handles_.addDim(dim, d);
        handles_.incrementIndex(d*strides_[dim]);
    }

    void setDim(int dim, MultiArrayIndex d)
    {
        d -= point(dim);
        handles_.addDim(dim, d);
        handles_.incrementIndex(d*strides_[dim]);
    }

    void resetDim(int dim)
    {
        MultiArrayIndex d = -point(dim);
        handles_.addDim(dim, d);
        handles_.incrementIndex(d*strides_[dim]);
    }

    CoupledScanOrderIterator & operator++()
    {
        handles_.template increment<dimension>();
        handles_.incrementIndex();
        return *this;
    }

    CoupledScanOrderIterator operator++(int)
    {
        CoupledScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    CoupledScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        // FIXME: this looks very expensive
        shape_type coordOffset;
        detail::ScanOrderToCoordinate<N>::exec(i+scanOrderIndex(), shape(), coordOffset);
        coordOffset -= point();
        handles_.add(coordOffset);
        handles_.scanOrderIndex_ += i;
        return *this;
    }

    CoupledScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        handles_.add(coordOffset);
        handles_.scanOrderIndex_ += detail::CoordinateToScanOrder<N>::exec(shape(), coordOffset);
        return *this;
    }

    CoupledScanOrderIterator & operator--()
    {
        handles_.template decrement<dimension>();
        handles_.decrementIndex();
        return *this;
    }

    CoupledScanOrderIterator operator--(int)
    {
        CoupledScanOrderIterator res(*this);
        --this;
        return res;
    }

    CoupledScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    CoupledScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
    }

    value_type operator[](MultiArrayIndex i) const
    {
        return *(CoupledScanOrderIterator(*this) += i);
    }

    CoupledScanOrderIterator
    operator+(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) += d;
    }

    CoupledScanOrderIterator
    operator-(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) -= d;
    }

    CoupledScanOrderIterator operator+(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) += coordOffset;
    }
    
    CoupledScanOrderIterator operator-(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) -= coordOffset;
    }

    MultiArrayIndex
    operator-(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() - r.scanOrderIndex();
    }

    bool operator==(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() == r.scanOrderIndex();
    }

    bool operator!=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() != r.scanOrderIndex();
    }

    bool operator<(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() < r.scanOrderIndex();
    }

    bool operator<=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() <= r.scanOrderIndex();
    }

    bool operator>(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() > r.scanOrderIndex();
    }

    bool operator>=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() >= r.scanOrderIndex();
    }

    bool isValid() const
    {
        return handles_.scanOrderIndex() < prod(shape());
    }

    bool atEnd() const
    {
        return handles_.scanOrderIndex() >= prod(shape());
    }

    MultiArrayIndex scanOrderIndex() const
    {
        return handles_.scanOrderIndex();
    }

    shape_type const & coord() const
    {
        return handles_.point();
    }

    MultiArrayIndex coord(unsigned int dim) const
    {
        return coord()[dim];
    }

    shape_type const & point() const
    {
        return handles_.point();
    }

    MultiArrayIndex point(unsigned int dim) const
    {
        return point()[dim];
    }

    shape_type const & shape() const
    {
        return handles_.shape();
    }

    MultiArrayIndex shape(unsigned int dim) const
    {
        return handles_.shape()[dim];
    }

    reference operator*()
    {
        return handles_;
    }

    const_reference operator*() const
    {
        return handles_;
    }

    CoupledScanOrderIterator & 
    restrictToSubarray(shape_type const & start, shape_type const & end)
    {
        operator+=(-point());
        handles_.restrictToSubarray(start, end);
        strides_ = detail::defaultStride(shape());
        return *this;
    }

    CoupledScanOrderIterator getEndIterator() const
    {
        return operator+(prod(shape())-scanOrderIndex());
    }

    bool atBorder() const
    {
        return (handles_.borderType() != 0);
    }

    unsigned int borderType() const
    {
        return handles_.borderType();
    }

    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::reference
    get() 
    {
        return vigra::get<TARGET_INDEX>(handles_);
    }

    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(handles_);
    }
    
    reference handles()
    {
        return handles_;
    }
    
    const_reference handles() const
    {
        return handles_;
    }

  protected:
    void reset()
    {
        handles_.template decrement<dimension>(shape()[dimension]);
    }

    void inverseReset()
    {
        handles_.template increment<dimension>(shape()[dimension]);
    }
    
    value_type handles_;
    shape_type strides_;
};


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

/** Helper class to easliy get the type of a CoupledScanOrderIterator (and corresponding CoupledHandle) for up to five arrays of dimension N with element types T1,...,T5.
 */
template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct CoupledIteratorType
{
    /** Type of the CoupledHandle.*/
    typedef typename CoupledHandleType<N, T1, T2, T3, T4, T5>::type HandleType;
  
    /** Type of the CoupledScanOrderIterator.*/
    typedef CoupledScanOrderIterator<HandleType::dimensions, HandleType> IteratorType;
    typedef IteratorType                                                 type;
};

/** Alias for \ref vigra::CoupledIteratorType (maybe easier to remember).
 */
template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct CoupledArrays
: public CoupledIteratorType<N, T1, T2, T3, T4, T5>
{};

/** Returns a CoupledScanOrderIterator from shape to iterate over coordinates. 
 */
template <int N>
typename CoupledIteratorType<N>::type
createCoupledIterator(TinyVector<MultiArrayIndex, N> const & shape)
{
    typedef typename CoupledHandleType<N>::type   P0;
    typedef CoupledScanOrderIterator<N, P0> IteratorType;
    
    return IteratorType(P0(shape));
}

/** Returns a CoupledScanOrderIterator to simultaneously iterate over image m1 and its coordinates. 
 */
template <unsigned int N1, class T1, class S1>
typename CoupledIteratorType<N1, T1>::type
createCoupledIterator(MultiArrayView<N1, T1, S1> const & m1)
{
    typedef typename CoupledHandleType<N1, T1>::type             P1;
    typedef typename P1::base_type                               P0;
    typedef CoupledScanOrderIterator<P1::dimensions, P1>         IteratorType;
    
    return IteratorType(P1(m1, 
                        P0(m1.shape())));
}

/** Returns a CoupledScanOrderIterator to simultaneously iterate over images m1, m2 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2>
typename CoupledIteratorType<N1, T1, T2>::type
createCoupledIterator(MultiArrayView<N1, T1, S1> const & m1,
                      MultiArrayView<N2, T2, S2> const & m2)
{
    typedef typename CoupledHandleType<N1, T1, T2>::type         P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef CoupledScanOrderIterator<P2::dimensions, P2> IteratorType;
    
    return IteratorType(P2(m2, 
                        P1(m1, 
                        P0(m1.shape()))));
}

/** Returns a CoupledScanOrderIterator to simultaneously iterate over images m1, m2, m3 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3>
typename CoupledIteratorType<N1, T1, T2, T3>::type
createCoupledIterator(MultiArrayView<N1, T1, S1> const & m1,
                      MultiArrayView<N2, T2, S2> const & m2,
                      MultiArrayView<N3, T3, S3> const & m3)
{
    typedef typename CoupledHandleType<N1, T1, T2, T3>::type     P3;
    typedef typename P3::base_type                               P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef CoupledScanOrderIterator<P3::dimensions, P3> IteratorType;
    
    return IteratorType(P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape())))));
}

/** Returns a CoupledScanOrderIterator to simultaneously iterate over images m1, m2, m3, m4 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3,
          unsigned int N4, class T4, class S4>
typename CoupledIteratorType<N1, T1, T2, T3, T4>::type
createCoupledIterator(MultiArrayView<N1, T1, S1> const & m1,
                      MultiArrayView<N2, T2, S2> const & m2,
                      MultiArrayView<N3, T3, S3> const & m3,
                      MultiArrayView<N4, T4, S4> const & m4)
{
    typedef typename CoupledHandleType<N1, T1, T2, T3, T4>::type P4;
    typedef typename P4::base_type                               P3;
    typedef typename P3::base_type                               P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef CoupledScanOrderIterator<P4::dimensions, P4> IteratorType;
    
    return IteratorType(P4(m4, 
                        P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape()))))));
}

/** Returns a CoupledScanOrderIterator to simultaneously iterate over images m1, m2, m3, m4, m5 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3,
          unsigned int N4, class T4, class S4,
          unsigned int N5, class T5, class S5>
typename CoupledIteratorType<N1, T1, T2, T3, T4, T5>::type
createCoupledIterator(MultiArrayView<N1, T1, S1> const & m1,
                      MultiArrayView<N2, T2, S2> const & m2,
                      MultiArrayView<N3, T3, S3> const & m3,
                      MultiArrayView<N4, T4, S4> const & m4,
                      MultiArrayView<N5, T5, S5> const & m5)
{
    typedef typename CoupledHandleType<N1, T1, T2, T3, T4, T5>::type P5;
    typedef typename P5::base_type                                   P4;
    typedef typename P4::base_type                                   P3;
    typedef typename P3::base_type                                   P2;
    typedef typename P2::base_type                                   P1;
    typedef typename P1::base_type                                   P0;
    typedef CoupledScanOrderIterator<P1::dimensions, P5> IteratorType;
    
    return IteratorType(P5(m5, 
                        P4(m4, 
                        P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape())))))));
}

//@}

} // namespace vigra

namespace std {

template <unsigned int N, class HANDLES, int DIMENSION>
ostream & operator<<(ostream & o, vigra::CoupledScanOrderIterator<N, HANDLES, DIMENSION> const & i)
{
    o << i.point();
    return o;
}

} // namespace std

#endif /* MULTI_ITERATOR_COUPLED_HXX */
