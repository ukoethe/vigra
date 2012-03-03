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

#ifndef MULTI_ITERATOR_COUPLED_HXX_
#define MULTI_ITERATOR_COUPLED_HXX_

#include "metaprogramming.hxx"
#include "multi_iterator.hxx"

namespace vigra {

template <unsigned int N, unsigned int DIMENSION=N-1>
struct NeighborhoodTypeImpl
{
    typedef typename MultiArrayShape<N>::type shape_type;
    
    static unsigned int exec(shape_type const & point, shape_type const & shape)
    {
        unsigned int res = NeighborhoodTypeImpl<N, DIMENSION-1>::exec(point, shape);
        if(point[DIMENSION] == 0)
            res |= (1 << 2*DIMENSION);
        if(point[DIMENSION] == shape[DIMENSION]-1)
            res |= (2 << 2*DIMENSION);
        return res;
    }
};

template <unsigned int N>
struct NeighborhoodTypeImpl<N, 0>
{
    typedef typename MultiArrayShape<N>::type shape_type;
    static const unsigned int DIMENSION = 0;
    
    static unsigned int exec(shape_type const & point, shape_type const & shape)
    {
        unsigned int res = 0;
        if(point[DIMENSION] == 0)
            res |= (1 << 2*DIMENSION);
        if(point[DIMENSION] == shape[DIMENSION]-1)
            res |= (2 << 2*DIMENSION);
        return res;
    }
};

template <class T, class NEXT>
class CoupledHandle
: public NEXT
{
public:
    typedef NEXT                            base_type;
    typedef CoupledHandle<T, NEXT>          self_type;
    
    static const int index =                NEXT::index + 1;    // index of this member of the chain

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

    template<unsigned int DIMENSION>
    inline void increment() 
    {
        pointer_ += strides_[DIMENSION];
        base_type::increment<DIMENSION>();
    }
    
    template<int DIMENSION>
    inline void decrement() 
    {
        pointer_ -= strides_[DIMENSION];
        base_type::decrement<DIMENSION>();
    }
    
    // TODO: test if making the above a default case of the this hurts performance
    template<int DIMENSION>
    inline void increment(MultiArrayIndex offset) 
    {
        pointer_ += offset*strides_[DIMENSION];
        base_type::increment<DIMENSION>(offset);
    }
    
    template<int DIMENSION>
    inline void decrement(MultiArrayIndex offset) 
    {
        pointer_ -= offset*strides_[DIMENSION];
        base_type::decrement<DIMENSION>(offset);
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


template <unsigned int N>
class CoupledHandle<TinyVector<MultiArrayIndex, N>, void>
{
public:
    static const int index = 0;                   // index of this member of the chain

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

    template<unsigned int DIMENSION>
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
    
    unsigned int neighborhoodType() const
    {
        return NeighborhoodTypeImpl<N>::exec(point_, shape_);
    }

    value_type point_, shape_;
    MultiArrayIndex scanOrderIndex_;
};

template <unsigned TARGET_INDEX, class Handle, unsigned int INDEX=Handle::index>
struct CoupledHandleCast
{
    typedef typename CoupledHandleCast<TARGET_INDEX, typename Handle::base_type>::type type;
};

template <unsigned TARGET_INDEX, class Handle>
struct CoupledHandleCast<TARGET_INDEX, Handle, TARGET_INDEX>
{
    typedef Handle type;
};

template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::type &
cast(Handle & handle)
{
    return handle;
};

template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::type const &
cast(Handle const & handle)
{
    return handle;
};

template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::type::reference
get(Handle & handle)
{
    return *cast<TARGET_INDEX>(handle);
};

template <unsigned int TARGET_INDEX, class Handle>
typename CoupledHandleCast<TARGET_INDEX, Handle>::type::const_reference
get(Handle const & handle)
{
    return *cast<TARGET_INDEX>(handle);
};

/********************************************************/
/*                                                      */
/*               CoupledScanOrderIterator<N>            */
/*                                                      */
/********************************************************/

/** \brief Iterate over multiple images simultaneously in scan order. Zero images is a special case;
    The coordinates can be accessed as a special band.

<b>\#include</b> \<vigra/multi_iterator_coupled.hxx\>

Namespace: vigra
*/

template <unsigned int N,
          class HANDLES,
          unsigned int DIMENSION = N-1>
class CoupledScanOrderIterator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: protected CoupledScanOrderIterator<N, HANDLES, DIMENSION-1>
#endif
{
    typedef CoupledScanOrderIterator<N, HANDLES, DIMENSION-1> base_type;
    enum { dimension = DIMENSION };

  public:

    typedef typename MultiArrayShape<dimension+1>::type        shape_type;
    typedef MultiArrayIndex                   difference_type;
    typedef CoupledScanOrderIterator          iterator;
    typedef std::random_access_iterator_tag   iterator_category;

    typedef typename base_type::value_type      value_type;
    typedef typename base_type::reference       reference;
    typedef typename base_type::const_reference const_reference; // FIXME: do we need both?
    typedef typename base_type::pointer         pointer;

    CoupledScanOrderIterator(value_type const & handles = value_type())
    : base_type(handles)
    {}

    reference operator[](MultiArrayIndex i)
    {
        return *(CoupledScanOrderIterator(*this) += i);
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        return *(CoupledScanOrderIterator(*this) += i);
    }

    CoupledScanOrderIterator & operator++()
    {
        base_type::operator++();
        if(this->point()[dimension-1] == this->shape()[dimension-1])
        {
            base_type::reset();
            this->handles_.template increment<dimension>();
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
        // FIXME: this looks very expensive
        shape_type coordOffset;
        detail::ScanOrderToCoordinate<N>::exec(i+scanOrderIndex(), this->shape(), coordOffset);
        coordOffset -= point();
        moveRelative(coordOffset);
        this->handles_.scanOrderIndex_ += i;
        return *this;
    }

    CoupledScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        moveRelative(coordOffset);
        this->handles_.scanOrderIndex_ += detail::CoordinateToScanOrder<N>::exec(this->shape(), coordOffset);
        return *this;
    }

    CoupledScanOrderIterator & operator--()
    {
        base_type::operator--();
        if(this->point()[dimension-1] == -1)
        {
            base_type::inverseReset();
            this->handles_.template decrement<dimension>();
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

    CoupledScanOrderIterator getEndIterator() const
    {
        return operator+(prod(this->shape()));
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

    bool operator==(CoupledScanOrderIterator const & r)
    {
        return base_type::operator==(r);
    }

    bool operator!=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator!=(r);
    }

    bool operator<(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator<(r);
    }

    bool operator<=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator<=(r);
    }

    bool operator>(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator>(r);
    }

    bool operator>=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator>=(r);
    }

    using base_type::operator*;
    using base_type::point;
    using base_type::shape;
    using base_type::scanOrderIndex;
    using base_type::atBorder;
    using base_type::neighborhoodType;
    using base_type::get;

    // using base_type::top;
    // using base_type::get;

  protected:
    void reset()
    {
        this->handles_.template decrement<dimension>(this->shape()[dimension]);
    }

    void inverseReset()
    {
        this->handles_.template increment<dimension>(this->shape()[dimension]);
    }

    void moveRelative(typename value_type::shape_type const & coordOffset)
    {
        base_type::moveRelative(coordOffset);
        this->handles_.template increment<dimension>(coordOffset[dimension]);
    }
};



template <unsigned int N, class HANDLES>
class CoupledScanOrderIterator<N, HANDLES, 0>
{
    enum { dimension = 0 };

  public:

    typedef CoupledScanOrderIterator<N, HANDLES, 0>  self_type;
    typedef HANDLES                                  value_type;
    typedef MultiArrayIndex                          difference_type;
    typedef value_type &                             reference;
    typedef value_type const &                       const_reference; 
    typedef value_type *                             pointer;
    typedef typename MultiArrayShape<1>::type        shape_type;
    typedef CoupledScanOrderIterator                 iterator;
    typedef std::random_access_iterator_tag          iterator_category;

    CoupledScanOrderIterator(value_type const & handles = value_type())
    : handles_(handles)
    {}

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
        shape_type coordOffset;
        detail::ScanOrderToCoordinate<N>::exec(i, shape(), coordOffset);
        moveRelative(coordOffset);
        handles_.scanOrderIndex_ += i;
        return *this;
    }

    CoupledScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        moveRelative(coordOffset);
        handles_.scanOrderIndex_ += detail::CoordinateToScanOrder<N>::exec(shape(), coordOffset);
        return *this;
    }

    CoupledScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
        return operator+=(-coordOffset);
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

    reference operator[](MultiArrayIndex i)
    {
        return *(CoupledScanOrderIterator(*this) += i);
    }

    const_reference operator[](MultiArrayIndex i) const
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

    bool
    operator==(CoupledScanOrderIterator const & r)
    {
        return scanOrderIndex() == r.scanOrderIndex();
    }

    bool
    operator!=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() != r.scanOrderIndex();
    }

    bool
    operator<(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() < r.scanOrderIndex();
    }

    bool
    operator<=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() <= r.scanOrderIndex();
    }

    bool
    operator>(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() > r.scanOrderIndex();
    }

    bool
    operator>=(CoupledScanOrderIterator const & r) const
    {
        return scanOrderIndex() >= r.scanOrderIndex();
    }

    MultiArrayIndex scanOrderIndex() const
    {
        return handles_.scanOrderIndex();
    }

    typename value_type::shape_type const & point() const
    {
        return handles_.point();
    }

    typename value_type::shape_type const & shape() const
    {
        return handles_.shape();
    }

    reference operator*()
    {
        return handles_;
    }

    const_reference operator*() const
    {
        return handles_;
    }

    void restrictToSubarray(shape_type const & start, shape_type const & end) const
    {
        operator+=(-point());
        handles_.restricToSubarray(start, end);
    }

    CoupledScanOrderIterator getEndIterator() const
    {
        return operator+(prod(shape()));
    }

    bool atBorder() const
    {
        return (handles_.neighborhoodType() != 0);
    }

    unsigned int neighborhoodType() const
    {
        return handles_.neighborhoodType();
    }

    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::type::reference
    get() 
    {
        return vigra::get<TARGET_INDEX>(handles_);
    }

    template<unsigned int TARGET_INDEX> 
    typename CoupledHandleCast<TARGET_INDEX, value_type>::type::const_reference
    get() const
    {
        return vigra::get<TARGET_INDEX>(handles_);
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

    void moveRelative(typename value_type::shape_type const & coordOffset)
    {
        handles_.template increment<dimension>(coordOffset[dimension]);
    }   
    
    value_type handles_;
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

template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct CoupledIteratorType
{
    typedef typename CoupledHandleType<N, T1, T2, T3, T4, T5>::type HandleType;
    typedef CoupledScanOrderIterator<N, HandleType> type;
};

template <unsigned int N>
typename CoupledIteratorType<N>::type
createCoupledIterator(TinyVector<MultiArrayIndex, N> const & shape)
{
    typedef typename CoupledHandleType<N>::type   P0;
    typedef CoupledScanOrderIterator<N, P0> IteratorType;
    
    return IteratorType(P0(shape));
}

template <unsigned int N, class T1, class S1>
typename CoupledIteratorType<N, T1>::type
createCoupledIterator(MultiArrayView<N, T1, S1> const & m1)
{
    typedef typename CoupledHandleType<N, T1>::type             P1;
    typedef typename P1::base_type                              P0;
    typedef CoupledScanOrderIterator<N, P1> IteratorType;
    
    return IteratorType(P1(m1.data(), m1.stride(), 
                        P0(m1.shape())));
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2>
typename CoupledIteratorType<N, T1, T2>::type
createCoupledIterator(MultiArrayView<N, T1, S1> const & m1,
                      MultiArrayView<N, T2, S2> const & m2)
{
    typedef typename CoupledHandleType<N, T1, T2>::type         P2;
    typedef typename P2::base_type                              P1;
    typedef typename P1::base_type                              P0;
    typedef CoupledScanOrderIterator<N, P2> IteratorType;
    
    return IteratorType(P2(m2.data(), m2.stride(), 
                        P1(m1.data(), m1.stride(), 
                        P0(m1.shape()))));
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3>
typename CoupledIteratorType<N, T1, T2, T3>::type
createCoupledIterator(MultiArrayView<N, T1, S1> const & m1,
                      MultiArrayView<N, T2, S2> const & m2,
                      MultiArrayView<N, T3, S3> const & m3)
{
    typedef typename CoupledHandleType<N, T1, T2, T3>::type     P3;
    typedef typename P3::base_type                              P2;
    typedef typename P2::base_type                              P1;
    typedef typename P1::base_type                              P0;
    typedef CoupledScanOrderIterator<N, P3> IteratorType;
    
    return IteratorType(P3(m3.data(), m3.stride(), 
                        P2(m2.data(), m2.stride(), 
                        P1(m1.data(), m1.stride(), 
                        P0(m1.shape())))));
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3,
                          class T4, class S4>
typename CoupledIteratorType<N, T1, T2, T3, T4>::type
createCoupledIterator(MultiArrayView<N, T1, S1> const & m1,
                      MultiArrayView<N, T2, S2> const & m2,
                      MultiArrayView<N, T3, S3> const & m3,
                      MultiArrayView<N, T4, S4> const & m4)
{
    typedef typename CoupledHandleType<N, T1, T2, T3, T4>::type P4;
    typedef typename P4::base_type                              P3;
    typedef typename P3::base_type                              P2;
    typedef typename P2::base_type                              P1;
    typedef typename P1::base_type                              P0;
    typedef CoupledScanOrderIterator<N, P4> IteratorType;
    
    return IteratorType(P4(m4.data(), m4.stride(), 
                        P3(m3.data(), m3.stride(), 
                        P2(m2.data(), m2.stride(), 
                        P1(m1.data(), m1.stride(), 
                        P0(m1.shape()))))));
}

template <unsigned int N, class T1, class S1,
                          class T2, class S2,
                          class T3, class S3,
                          class T4, class S4,
                          class T5, class S5>
typename CoupledIteratorType<N, T1, T2, T3, T4, T5>::type
createCoupledIterator(MultiArrayView<N, T1, S1> const & m1,
                      MultiArrayView<N, T2, S2> const & m2,
                      MultiArrayView<N, T3, S3> const & m3,
                      MultiArrayView<N, T4, S4> const & m4,
                      MultiArrayView<N, T5, S5> const & m5)
{
    typedef typename CoupledHandleType<N, T1, T2, T3, T4, T5>::type P5;
    typedef typename P5::base_type                                  P4;
    typedef typename P4::base_type                                  P3;
    typedef typename P3::base_type                                  P2;
    typedef typename P2::base_type                                  P1;
    typedef typename P1::base_type                                  P0;
    typedef CoupledScanOrderIterator<N, P5> IteratorType;
    
    return IteratorType(P5(m5.data(), m5.stride(), 
                        P4(m4.data(), m4.stride(), 
                        P3(m3.data(), m3.stride(), 
                        P2(m2.data(), m2.stride(), 
                        P1(m1.data(), m1.stride(), 
                        P0(m1.shape())))))));
}

}; // namespace vigra

#endif /* MULTI_ITERATOR_COUPLED_HXX_ */
