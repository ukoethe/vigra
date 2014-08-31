/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2012 by Gunnar Kedenburg and Ullrich Koethe       */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.3.0, Sep 10 2004 )                                    */
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

#ifndef VIGRA_MULTI_HIERARCHICAL_ITERATOR_HXX
#define VIGRA_MULTI_HIERARCHICAL_ITERATOR_HXX

#include <sys/types.h>
#include "multi_fwd.hxx"
#include "iteratortags.hxx"
#include "multi_handle.hxx"

namespace vigra {

/** \addtogroup MultiIteratorGroup
*/
//@{

/********************************************************/
/*                                                      */
/*                 MultiHierarchicalIterator            */
/*                                                      */
/********************************************************/

template <unsigned int N,
          class HANDLES,
          int DIMENSION = N-1>
class HierarchicalIterator
: public HierarchicalIterator<N, HANDLES, DIMENSION-1>
{
  public:

    typedef HierarchicalIterator<N, HANDLES, DIMENSION-1> base_type;
    static const int level = DIMENSION;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::const_pointer const_pointer;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::shape_type shape_type;

    explicit HierarchicalIterator(HANDLES const & handles = HANDLES())
    : base_type(handles)
    {}

    void operator++ ()
    {
        this->handles_.template increment<level>();
    }

    void operator-- ()
    {
        this->handles_.template decrement<level>();
    }

    HierarchicalIterator operator++ (int)
    {
        HierarchicalIterator ret = *this;
        ++(*this);
        return ret;
    }

    HierarchicalIterator operator-- (int)
    {
        HierarchicalIterator ret = *this;
        --(*this);
        return ret;
    }

    HierarchicalIterator & operator+= (difference_type n)
    {
        this->handles_.addDim(level, n);
        return *this;
    }

    HierarchicalIterator & operator+= (shape_type const & d)
    {
        base_type::operator+=(d);
        return *this;
    }

    HierarchicalIterator  &operator-= (difference_type n)
    {
        this->handles_.addDim(level, -n);
        return *this;
    }

    // HierarchicalIterator & operator-= (multi_difference_type const & d)
    // {
        // this->m_ptr -= total_stride(d.begin());
        // return *this;
    // }

    HierarchicalIterator operator+ (difference_type n) const
    {
        return HierarchicalIterator(*this) += n;
    }

    // HierarchicalIterator operator+ (multi_difference_type const & d) const
    // {
        // StridedMultiIterator ret = *this;
        // ret += d;
        // return ret;
    // }

    difference_type operator- (HierarchicalIterator const & d) const
    {
        return this->point()[level] - d.point()[level];
    }

    HierarchicalIterator operator- (difference_type n) const
    {
        return HierarchicalIterator(*this) -= n;
    }

    // HierarchicalIterator operator- (multi_difference_type const & d) const
    // {
        // HierarchicalIterator ret = *this;
        // ret -= d;
        // return ret;
    // }

    bool operator== (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] == rhs.point()[level];
    }

    bool operator!= (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] != rhs.point()[level];
    }

    bool operator< (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] < rhs.point()[level];
    }

    bool operator<= (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] <= rhs.point()[level];
    }

    bool operator> (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] > rhs.point()[level];
    }

    bool operator>= (const HierarchicalIterator &rhs) const
    {
        return this->point()[level] >= rhs.point()[level];
    }

    base_type begin () const
    {
        return *this;
    }

    base_type end () const
    {
        return base_type(*this) += this->shape()[level-1] - this->point()[level-1];
    }
    
    HierarchicalIterator getEndIterator() const
    {
        return HierarchicalIterator(*this) += this->shape() - this->point();
    }

    // iterator iteratorForDimension(unsigned int d) const
    // {
        // vigra_precondition(d <= level,
            // "StridedMultiIterator<N>::iteratorForDimension(d): d < N required");
        // return iterator(this->m_ptr, stride_traits::shift(m_stride, d), 0);
    // }

    template <int K>
    HierarchicalIterator<N, HANDLES, K> &
    dim()
    {
        return *this;
    }

    template <int K>
    HierarchicalIterator<N, HANDLES, K> const &
    dim() const
    {
        return *this;
    }
};

/********************************************************/

template <unsigned int N,
          class HANDLES>
class HierarchicalIterator<N, HANDLES, 0>
{
  public:
    static const int level = 0;
    
    typedef CoupledHandleTraits<HANDLES> HandleTraits;
    typedef typename HandleTraits::value_type value_type;
    typedef typename HandleTraits::reference reference;
    typedef typename HandleTraits::const_reference const_reference;
    typedef typename HandleTraits::pointer pointer;
    typedef typename HandleTraits::const_pointer const_pointer;
    typedef MultiArrayIndex difference_type;
    typedef typename MultiArrayShape<N>::type shape_type;
    typedef HierarchicalIterator<N, HANDLES, 0> iterator;
    typedef std::random_access_iterator_tag iterator_category;

  protected:
    HANDLES handles_;

  public:
    explicit HierarchicalIterator(HANDLES const & handles = HANDLES())
    : handles_(handles)
    {}

    void operator++ ()
    {
        handles_.template increment<level>();
    }

    void operator-- ()
    {
        handles_.template decrement<level>();
    }

    HierarchicalIterator operator++ (int)
    {
        HierarchicalIterator ret = *this;
        ++(*this);
        return ret;
    }

    HierarchicalIterator operator-- (int)
    {
        HierarchicalIterator ret = *this;
        --(*this);
        return ret;
    }

    HierarchicalIterator & operator+= (difference_type n)
    {
        handles_.addDim(level, n);
        return *this;
    }

    HierarchicalIterator & operator+= (shape_type const & d)
    {
        handles_.add(d);
        handles_.scanOrderIndex_ += detail::CoordinateToScanOrder<N>::exec(shape(), d);
        return *this;
    }

    HierarchicalIterator & operator-= (difference_type n)
    {
        handles_.addDim(level, -n);
        return *this;
    }

    // HierarchicalIterator & operator-= (multi_difference_type const & d)
    // {
        // m_ptr -= d[level];
        // return *this;
    // }

    HierarchicalIterator operator+ (difference_type n) const
    {
        return HierarchicalIterator(*this) += n;
    }

    // HierarchicalIterator operator+ (multi_difference_type const & d) const
    // {
        // HierarchicalIterator ret = *this;
        // ret += d;
        // return ret;
    // }

    difference_type operator- (HierarchicalIterator const & d) const
    {
        return point()[level] - d.point()[level];
    }

    HierarchicalIterator operator- (difference_type n) const
    {
        return HierarchicalIterator(*this) -= n;
    }

    // HierarchicalIterator operator- (multi_difference_type const & d) const
    // {
        // HierarchicalIterator ret = *this;
        // ret -= d;
        // return ret;
    // }

    // reference operator[] (difference_type n) const
    // {
        // return m_ptr [n];
    // }

    // reference operator[] (multi_difference_type const & d) const
    // {
        // return m_ptr [d[level]];
    // }
    
    reference operator* ()
    {
        return HandleTraits::dereference(handles_);
    }

    const_reference operator* () const
    {
        return HandleTraits::dereference(handles_);
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, HANDLES>::reference
    get()
    {
        return handles_.get<TARGET_INDEX>();
    }
    
    template <unsigned int TARGET_INDEX>
    typename CoupledHandleCast<TARGET_INDEX, HANDLES>::const_reference
    get() const
    {
        return handles_.get<TARGET_INDEX>();
    }

    pointer operator->()
    {
        return &HandleTraits::dereference(handles_);
    }

    const_pointer operator->() const
    {
        return &HandleTraits::dereference(handles_);
    }

    bool operator== (const HierarchicalIterator &rhs) const
    {
        return point()[level] == rhs.point()[level];
    }

    bool operator!= (const HierarchicalIterator &rhs) const
    {
        return point()[level] != rhs.point()[level];
    }

    bool operator< (const HierarchicalIterator &rhs) const
    {
        return point()[level] < rhs.point()[level];
    }

    bool operator<= (const HierarchicalIterator &rhs) const
    {
        return point()[level] <= rhs.point()[level];
    }

    bool operator> (const HierarchicalIterator &rhs) const
    {
        return point()[level] > rhs.point()[level];
    }

    bool operator>= (const HierarchicalIterator &rhs) const
    {
        return point()[level] >= rhs.point()[level];
    }

    // iterator iteratorForDimension(unsigned int d) const
    // {
        // vigra_precondition(d == 0,
            // "MultiIterator<1>::iteratorForDimension(d): d == 0 required");
        // const difference_type stride = 1;
        // return iterator(m_ptr, &stride, 0);
    // }
    
    HierarchicalIterator getEndIterator() const
    {
        return HierarchicalIterator(*this) += shape() - point();
    }

    shape_type const & point() const
    {
        return handles_.point();
    }

    shape_type const & shape() const
    {
        return handles_.shape();
    }

    difference_type scanOrderIndex() const
    {
        return handles_.scanOrderIndex();
    }
    
    HANDLES const & handles() const
    {
        return handles_;
    }

    template <int K>
    HierarchicalIterator<N, HANDLES, 0> &
    dim()
    {
        return *this;
    }
    
    template <int K>
    HierarchicalIterator<N, HANDLES, 0> const &
    dim() const
    {
        return *this;
    }

  // protected:

    // difference_type 
    // total_stride(typename multi_difference_type::const_iterator d) const
    // {
        // return d[level];
    // }
};

/** Helper class to easliy get the type of a CoupledScanOrderIterator (and corresponding CoupledHandle) for up to five arrays of dimension N with element types T1,...,T5.
 */
template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
struct HierarchicalIteratorType
{
    /** Type of the CoupledHandle.*/
    typedef typename CoupledHandleType<N, T1, T2, T3, T4, T5>::type HandleType;
  
    /** Type of the CoupledScanOrderIterator.*/
    typedef HierarchicalIterator<HandleType::dimensions, HandleType> IteratorType;
    typedef IteratorType                                             type;
};

// /** Alias for \ref vigra::CoupledIteratorType (maybe easier to remember).
 // */
// template <unsigned int N, class T1=void, class T2=void, class T3=void, class T4=void, class T5=void>
// struct CoupledArrays
// : public CoupledIteratorType<N, T1, T2, T3, T4, T5>
// {};

/** Returns a HierarchicalIterator from shape to iterate over coordinates. 
 */
template <int N>
typename HierarchicalIteratorType<N>::type
createHierarchicalIterator(TinyVector<MultiArrayIndex, N> const & shape)
{
    typedef typename CoupledHandleType<N>::type   P0;
    typedef HierarchicalIterator<N, P0> IteratorType;
    
    return IteratorType(P0(shape));
}

/** Returns a HierarchicalIterator to simultaneously iterate over image m1 and its coordinates. 
 */
template <unsigned int N1, class T1, class S1>
typename HierarchicalIteratorType<N1, T1>::type
createHierarchicalIterator(MultiArrayView<N1, T1, S1> const & m1)
{
    typedef typename CoupledHandleType<N1, T1>::type             P1;
    typedef typename P1::base_type                               P0;
    typedef HierarchicalIterator<P1::dimensions, P1>         IteratorType;
    
    return IteratorType(P1(m1, 
                        P0(m1.shape())));
}

/** Returns a HierarchicalIterator to simultaneously iterate over images m1, m2 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2>
typename HierarchicalIteratorType<N1, T1, T2>::type
createHierarchicalIterator(MultiArrayView<N1, T1, S1> const & m1,
                           MultiArrayView<N2, T2, S2> const & m2)
{
    typedef typename CoupledHandleType<N1, T1, T2>::type         P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef HierarchicalIterator<P2::dimensions, P2> IteratorType;
    
    return IteratorType(P2(m2, 
                        P1(m1, 
                        P0(m1.shape()))));
}

/** Returns a HierarchicalIterator to simultaneously iterate over images m1, m2, m3 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3>
typename HierarchicalIteratorType<N1, T1, T2, T3>::type
createHierarchicalIterator(MultiArrayView<N1, T1, S1> const & m1,
                           MultiArrayView<N2, T2, S2> const & m2,
                           MultiArrayView<N3, T3, S3> const & m3)
{
    typedef typename CoupledHandleType<N1, T1, T2, T3>::type     P3;
    typedef typename P3::base_type                               P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef HierarchicalIterator<P3::dimensions, P3> IteratorType;
    
    return IteratorType(P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape())))));
}

/** Returns a HierarchicalIterator to simultaneously iterate over images m1, m2, m3, m4 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3,
          unsigned int N4, class T4, class S4>
typename HierarchicalIteratorType<N1, T1, T2, T3, T4>::type
createHierarchicalIterator(MultiArrayView<N1, T1, S1> const & m1,
                           MultiArrayView<N2, T2, S2> const & m2,
                           MultiArrayView<N3, T3, S3> const & m3,
                           MultiArrayView<N4, T4, S4> const & m4)
{
    typedef typename CoupledHandleType<N1, T1, T2, T3, T4>::type P4;
    typedef typename P4::base_type                               P3;
    typedef typename P3::base_type                               P2;
    typedef typename P2::base_type                               P1;
    typedef typename P1::base_type                               P0;
    typedef HierarchicalIterator<P4::dimensions, P4> IteratorType;
    
    return IteratorType(P4(m4, 
                        P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape()))))));
}

/** Returns a HierarchicalIterator to simultaneously iterate over images m1, m2, m3, m4, m5 and their coordinates. 
 */
template <unsigned int N1, class T1, class S1,
          unsigned int N2, class T2, class S2,
          unsigned int N3, class T3, class S3,
          unsigned int N4, class T4, class S4,
          unsigned int N5, class T5, class S5>
typename HierarchicalIteratorType<N1, T1, T2, T3, T4, T5>::type
createHierarchicalIterator(MultiArrayView<N1, T1, S1> const & m1,
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
    typedef HierarchicalIterator<P5::dimensions, P5> IteratorType;
    
    return IteratorType(P5(m5, 
                        P4(m4, 
                        P3(m3, 
                        P2(m2, 
                        P1(m1, 
                        P0(m1.shape())))))));
}


template <unsigned int N, class A, class B>
HierarchicalIterator<N, typename ZipCoupledHandles<A, B>::type>
zip(HierarchicalIterator<N, A> const & a, HierarchicalIterator<N, B> const & b)
{
    vigra_precondition(a.shape() == b.shape() && a.scanOrderIndex() == b.scanOrderIndex(),
         "zip(HierarchicalIterator): iterators must have identical shape and position.");
         
    typedef typename ZipCoupledHandles<A, B>::type Handle;
    typedef HierarchicalIterator<N, Handle> IteratorType;
    return IteratorType(ZipCoupledHandles<A, B>::construct(a.handles(), b.handles()));
}

//@}

} // namespace vigra

// namespace std {

// template <unsigned int N, class T, class REFERENCE, class POINTER>
// ostream & operator<<(ostream & o, vigra::StridedScanOrderIterator<N, T, REFERENCE, POINTER> const & i)
// {
    // o << *i;
    // return o;
// }

// } // namespace std

#endif // VIGRA_MULTI_HIERARCHICAL_ITERATOR_HXX
