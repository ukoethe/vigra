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

#include "multi_fwd.hxx"
#include "multi_shape.hxx"
#include "multi_handle.hxx"
#include "metaprogramming.hxx"

namespace vigra {

/** \addtogroup MultiIteratorGroup
*/
//@{

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
          class HANDLES,
          int DIMENSION>  // NOTE: default template arguments are defined in multi_fwd.hxx
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

template <unsigned int N, class A, class B>
CoupledScanOrderIterator<N, typename ZipCoupledHandles<A, B>::type>
zip(CoupledScanOrderIterator<N, A> const & a, CoupledScanOrderIterator<N, B> const & b)
{
    vigra_precondition(a.shape() == b.shape() && a.scanOrderIndex() == b.scanOrderIndex(),
         "zip(CoupledScanOrderIterator): iterators must have identical shape and position.");
         
    typedef typename ZipCoupledHandles<A, B>::type Handle;
    typedef CoupledScanOrderIterator<N, Handle> IteratorType;
    return IteratorType(ZipCoupledHandles<A, B>::construct(*a, *b));
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
