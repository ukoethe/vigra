/************************************************************************/
/*                                                                      */
/*    Copyright 2011-2012 by Markus Nullmeier and Ullrich Koethe        */
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

#ifndef VIGRA_COORDINATE_ITERATOR_HXX
#define VIGRA_COORDINATE_ITERATOR_HXX

#include <complex>

#include "tuple.hxx"
#include "accessor.hxx"
#include "tinyvector.hxx"
#include "numerictraits.hxx"
#include "multi_iterator.hxx"
#include "multi_array.hxx"

namespace vigra {

template<unsigned N>
struct StridePair
{
    typedef typename MultiArrayShape<N>::type  index_type;
    typedef          TinyVector<double, N>     coord_type;
    typedef          coord_type                deref_type;
    typedef          StridePair                type;
    typedef          StridePair                stride_type;
    typedef          TinyVector<type, N>       stride_array_type;
    typedef          TinyVector<index_type, N> shape_array_type;
    typedef          shape_array_type          shape_type;

    index_type index;
    coord_type coord;

    StridePair(const index_type & i) : index(i), coord(i) {}
    StridePair(const coord_type & c) : index(),  coord(c) {}
    StridePair(const index_type & i, const coord_type & c)
        :            index       (i),             coord(c) {}
    StridePair(  MultiArrayIndex  i, const coord_type & c)
        :        index(index_type(i)),            coord(c) {}
    StridePair() {}

    // use just the coordinates for further processing ...
    const coord_type & operator*() const
    {
        return this->coord;
    }

    void operator+=(const StridePair & x)
    {
        index += x.index;
        coord += x.coord;
    }
    void operator-=(const StridePair & x)
    {
        index -= x.index;
        coord -= x.coord;
    }
    StridePair operator+(const StridePair & x)
    {
        StridePair ret = *this;
        ret += x;
        return ret;
    }
    StridePair operator-(const StridePair & x)
    {
        StridePair ret = *this;
        ret -= x;
        return ret;
    }
    StridePair operator*(const StridePair & x)
    {
        StridePair ret = *this;
        ret.index *= x.index;
        ret.coord *= x.coord;
        return ret;
    }
    StridePair operator/(const StridePair & x)
    {
        StridePair ret = *this;
        ret.index /= x.index;
        ret.coord /= x.coord;
        return ret;
    }

    MultiArrayIndex & idx0()
    {
        return index[0];
    }
    const index_type & idx() const
    {
        return index;
    }
    
    double & dim0()
    {
        return coord[0];
    }
    double dim0() const
    {
        return coord[0];
    }
};

template<unsigned M>
struct NumericTraits<StridePair<M> >
    : public NumericTraits<typename StridePair<M>::index_type>
{};

template<unsigned N>
struct StridePairCoord : public TinyVector<double, N>
{
    typedef TinyVector<double, N> entry_type;

    StridePairCoord(const entry_type & c) : entry_type(c) {}
    StridePairCoord() {}

    double & dim0()
    {
        return (*this)[0];
    }
    double dim0() const
    {
        return (*this)[0];
    }
};
template<unsigned M>
struct NumericTraits<StridePairCoord<M> >
    : public NumericTraits<typename StridePairCoord<M>::entry_type>
{};

template<unsigned N>
struct StridePairDiff : public StridePairCoord<N>
{
    MultiArrayIndex c;

    typedef StridePairCoord<N> base_type;
    StridePairDiff(MultiArrayIndex c_, const base_type & x)
        : base_type(x), c(c_) {}
    StridePairDiff(const base_type & x)
        : base_type(x), c(0) {}
    StridePairDiff(const TinyVector<double, N> & x)
        : base_type(x), c(0) {}
    StridePairDiff(const TinyVector<MultiArrayIndex, N> & x)
        : base_type(x), c(0) {}
    StridePairDiff() : c(0) {}

    const base_type & base() const
    {
        return *this;
    }
    StridePairDiff operator*(const StridePairDiff & x)
    {
        StridePairDiff ret = base() * x.base();
        ret.c = c * x.c;
        return ret;
    }
};

template<unsigned M>
struct NumericTraits<StridePairDiff<M> >
    : public NumericTraits<StridePairCoord<M> >
{};

template<unsigned N, class T>
struct StridePairPointer : public StridePairCoord<N>
{
    typedef const T*                          index_type;
    typedef StridePairCoord<N>                coord_type;
    typedef typename coord_type::entry_type   coord_num_type;
    typedef StridePairPointer                 type;
    typedef type                              deref_type;
    typedef StridePairDiff<N>                 stride_type;
    typedef TinyVector<stride_type, N>        stride_array_type;
    typedef typename MultiArrayShape<N>::type shape_array_type;
    typedef          shape_array_type         shape_type;

    index_type index;

    StridePairPointer(const index_type & i, const coord_type & c)
        : coord_type(c), index(i) {}

    const type & operator*() const
    {
        return *this;
    }
    const T & value() const
    {
        return *index;
    }
    const coord_type & coord() const
    {
        return *this;
    }

    index_type & idx0()
    {
        return index;
    }
    const index_type & idx() const
    {
        return index;
    }
    
    void operator+=(stride_type x)
    {
        index += x.c;
        coord_type::operator+=(x);
    }
    void operator-=(stride_type x)
    {
        index -= x.c;
        coord_type::operator-=(x);
    }
};

template<unsigned M, class T>
struct NumericTraits<StridePairPointer<M, T> >
    : public NumericTraits<typename StridePairPointer<M, T>::coord_type>
{};

namespace detail {

template<class T, bool is_complex = NumericTraits<T>::isComplex::value,
                  bool is_vector = !NumericTraits<T>::isScalar::value>
struct weighted_abs
{
    static double get(const T & x)
    {
        return x;
    }
};

template<class T>
struct weighted_abs<T, true, false>
{
    static double get(const T & x)
    {
        using std::abs;
        return abs(x);
    }
};

template<class T, bool is_complex>
struct weighted_abs<T, is_complex, true>
{
    static double get(const T & x)
    {
        return x.magnitude();
    }
};

template<class T>
struct accumulable_coord_access;
template<class T>
struct accumulable_value_access;
template<class T>
struct accumulable_weighted_access;

template<unsigned N, class T>
struct accumulable_coord_access<StridePairPointer<N, T> >
{
    typedef StridePairPointer<N, T> accumulable_type;
    typedef typename accumulable_type::coord_num_type type;
    static const type & get(const accumulable_type & v) { return v.coord(); }
};

template<unsigned N, class T>
struct accumulable_value_access<StridePairPointer<N, T> >
{
    typedef StridePairPointer<N, T> accumulable_type;
    typedef T type;
    static const type & get(const accumulable_type & v) { return v.value(); }
};

template<unsigned N, class T>
struct accumulable_weighted_access<StridePairPointer<N, T> >
{
    typedef StridePairPointer<N, T> accumulable_type;
    typedef typename accumulable_type::coord_num_type type;
    static type get(const accumulable_type & v)
    {
        return weighted_abs<T>::get(v.value()) * v.coord();
    }
};

template<class X>
void dismember(X & r, const X & x, unsigned i)
{
    r[i] = x[i];
}
template<unsigned N>
void dismember(StridePair<N> & r, const StridePair<N> & x, unsigned i)
{
    r.index[i] = x.index[i];
    r.coord[i] = x.coord[i];
}
template<unsigned N>
void dismember(StridePairDiff<N> & r, const StridePairDiff<N> & x, unsigned i)
{
    r.c = static_cast<MultiArrayIndex>(r[i] = x[i]);
}

template<unsigned N, class X>
TinyVector<X, N>
dismember(const X & x)
{
    TinyVector<X, N> ret;
    for (unsigned i = 0; i != N; ++i)
        dismember(ret[i], x, i);
    return ret;
}
template<unsigned N>
TinyVector<StridePairDiff<N>, N>
dismember(const TinyVector<MultiArrayIndex, N> & x,
          const StridePairCoord<N> &             y)
{
    typedef StridePairDiff<N> type;
    TinyVector<type, N> ret;
    for (unsigned i = 0; i != N; ++i)
    {
        ret[i].c  = x[i];
        ret[i][i] = y[i];
    }
    return ret;
}

} // namespace detail

// A fake "pointer" for MultiIterator containing coordinates.
// Indices (or a pointer) cannot be circumvented in coordiante iterators,
// since floating point addition is not associative and
// iterator comparison is done via via '<' or '!='. Class CoordinateStride
// thus forwards iterator comparison to the index or pointer part
// of its template parameter S.
template<unsigned N, class S = StridePair<N> >
class CoordinateStride : protected S
{
  public:
    typedef          MultiArrayIndex     difference_type;
    typedef typename S::stride_type      stride_type;
    typedef typename S::deref_type       deref_type;
    typedef          CoordinateStride<N> type;
    typedef typename S::coord_type       coord_type;
    typedef typename S::index_type       index_type;
    typedef typename S::shape_array_type shape_array_type;

  protected:
    double stride_0;

    CoordinateStride(void*) {} // used MultiIterator ctor, unused.

  public:
    CoordinateStride(const S & x, double s0)
        : S(x), stride_0(s0) {}

#ifndef DOXYGEN
    using S::operator*;
    using S::idx0;
    using S::idx;
    using S::dim0;
    using S::operator+=;
    using S::operator-=;
#endif

    void operator++()
    {
        ++idx0();
        dim0() += stride_0;
    }
    void operator--()
    {
        --idx0();
        dim0() -= stride_0;
    }
    void operator+=(difference_type n)
    {
        idx0() += n;
        dim0() += n * stride_0;
    }
    void operator-=(difference_type n)
    {
        idx0() -= n;
        dim0() -= n * stride_0;
    }

    stride_type operator[](difference_type n) const
    {
        type ret = *this;
        ret[0] += n;
        return ret;
    }

    stride_type operator[](stride_type x) const
    {
        return *this + x;
    }

    // ... but use the idx() for comparisons:
    bool operator!=(const CoordinateStride & y) const
    {
        return idx() != y.idx();
    }
    bool operator==(const CoordinateStride & y) const
    {
        if (stride_0 != y.stride_0)
            return false;
        return idx() == y.idx();
    }
    bool operator<(const CoordinateStride & y) const
    {
        return idx() < y.idx();
    }

    bool operator<=(const CoordinateStride & y) const
    {
        if (stride_0 == y.stride_0)
            return true;
        return *this < y;
    }
    bool operator>(const CoordinateStride & y) const
    {
        return y < *this;
    }
    bool operator>=(const CoordinateStride & y) const
    {
        if (stride_0 == y.stride_0)
            return true;
        return operator>(y);
    }

    friend std::ostream &
    operator<<(std::ostream & os, const CoordinateStride & x)
    {
        os << "{" << x.stride_0 << ": " << static_cast<const S &>(x) << "}";
        return os;
    }

    typedef MultiIterator<N, deref_type, const deref_type &, CoordinateStride>
        iterator_type;
};

template <unsigned N, class S>
struct MultiIteratorStrideTraits<CoordinateStride<N, S> >
{
    typedef typename S::stride_type           stride_type;
    typedef typename S::stride_array_type     stride_array_type;
    typedef typename S::shape_array_type      shape_array_type;
    static stride_array_type shift(const stride_array_type & s, unsigned d)
    {
        stride_array_type ret;
        for (unsigned i = d; i != N; ++i)
            ret[i - d] = s[i];
        return ret;
    }
};

template <unsigned N>
struct CoordinateMultiIterator : public CoordinateStride<N>::iterator_type
{
    typedef CoordinateStride<N>                       ptr_type;
    typedef typename ptr_type::iterator_type          base_type;
    typedef typename ptr_type::stride_type            stride_type;
    typedef typename ptr_type::shape_array_type       shape_array_type;
    typedef typename ptr_type::coord_type             coord_type;
    typedef typename ptr_type::index_type             index_type;

    CoordinateMultiIterator(const stride_type & origin,
                            const stride_type & stride,
                            const index_type & shape)

        : base_type(ptr_type(origin, stride.dim0()),
                    detail::dismember<N>(stride),
                    detail::dismember<N>(shape)) {}
                                                 
    CoordinateMultiIterator(const base_type & x) : base_type(x) {}
};

namespace detail {

template<unsigned N>
struct CoordinateMultiRangeReturns
{
    typedef          CoordinateMultiIterator<N>   iterator_type;
    typedef typename iterator_type::coord_type    coord_type;
    typedef          StridePair<N>                pair_type;
    typedef typename pair_type::type              stride_type;
    typedef typename pair_type::stride_array_type stride_array_type;

    typedef typename AccessorTraits<coord_type>::default_const_accessor
        access_type;
    typedef triple<iterator_type, stride_array_type, access_type> type;
};

} // namespace detail

template <unsigned N>
typename detail::CoordinateMultiRangeReturns<N>::type
coordinateMultiRange(const typename MultiArrayShape<N>::type & shape,
                     const TinyVector<double, N> & stride
                                                   = TinyVector<double, N>(1.0),
                     const TinyVector<double, N> & origin
                                                   = TinyVector<double, N>(0.0))
{
    typedef typename
        detail::CoordinateMultiRangeReturns<N>::stride_type stride_type;
    typedef typename
        detail::CoordinateMultiRangeReturns<N>::access_type access_type;

    return typename detail::CoordinateMultiRangeReturns<N>::type
        (CoordinateMultiIterator<N>(stride_type(0, origin),
                                    stride_type(1, stride),
                                    shape),
         detail::dismember<N>(stride_type(shape)),
         access_type());
}

template <unsigned N, class T>
struct CombinedMultiIterator
    : public CoordinateStride<N, StridePairPointer<N, T> >::iterator_type
{
    typedef StridePairPointer<N, T>                   pair_type;
    typedef CoordinateStride<N, pair_type>            ptr_type;
    typedef typename ptr_type::iterator_type          base_type;
    typedef typename ptr_type::stride_type            stride_type;
    typedef typename ptr_type::coord_type             coord_type;
    typedef typename pair_type::shape_array_type      shape_array_type;
    
    CombinedMultiIterator(const T*                               raw_pointer,
                          const stride_type &                    origin,
                          const TinyVector<MultiArrayIndex, N> & pointer_stride,
                          const stride_type &                    stride,
                          const shape_array_type &               shape)

        : base_type(ptr_type(pair_type(raw_pointer, origin), stride.dim0()),
                    detail::dismember<N>(pointer_stride, stride),
                    shape) {}
                                                 
    CombinedMultiIterator(const base_type & x) : base_type(x) {}
};

template<unsigned N, class T>
struct SrcCoordinateMultiArrayRangeReturns
{
    typedef          CombinedMultiIterator<N, T>  iterator_type;
    typedef typename iterator_type::coord_type    coord_type;
    typedef typename iterator_type::pair_type     pair_type;
    typedef typename iterator_type::ptr_type      ptr_type;
    typedef typename ptr_type::deref_type         deref_type;
    typedef typename iterator_type::stride_type   stride_type;
    typedef typename pair_type::stride_array_type stride_array_type;
    typedef typename pair_type::shape_array_type  shape_array_type;

    typedef typename AccessorTraits<deref_type>::default_const_accessor
        access_type;
    typedef triple<iterator_type, stride_array_type, access_type> type;
};

// work around GCC 4.4.3 template argument deduction bug:
template<unsigned N>
struct CoordinateSteps
{
    typedef const TinyVector<double, N> & type;
};

template <unsigned int N, class T, class StrideTag>
inline typename SrcCoordinateMultiArrayRangeReturns<N, T>::type
srcCoordinateMultiArrayRange(const MultiArrayView<N, T, StrideTag> & array,
                             typename CoordinateSteps<N>::type stride
                                                   = TinyVector<double, N>(1.0),
                             typename CoordinateSteps<N>::type origin
                                                   = TinyVector<double, N>(0.0))
{
    typedef SrcCoordinateMultiArrayRangeReturns<N, T> returns;
    typedef typename returns::type                    type;
    typedef typename returns::stride_type             stride_type;
    typedef typename returns::access_type             access_type;
    typedef typename returns::iterator_type           iterator_type;
    typedef typename returns::shape_array_type        shape_array_type;
    
    shape_array_type shape = array.shape();
    return type(iterator_type(array.traverser_begin().get(),
                              stride_type(origin),
                              array.stride(),
                              stride_type(stride),
                              shape),
                detail::dismember<N>(stride_type(shape)),
                access_type());
}

template <class VALUETYPE, class COORD>
struct AccessorCoordinatePair
{
    typedef VALUETYPE                  value_type;
    typedef COORD                      coord_type;
    typedef AccessorCoordinatePair     type;

          value_type   v;
    const coord_type & c;

    AccessorCoordinatePair(const value_type & v_, const coord_type & c_)
        : v(v_), c(c_) {}

    const value_type & value() const
    {
        return v;
    }
    const coord_type & coord() const
    {
        return c;
    }
};

/** \brief Forward accessor to the value() part of the values an iterator
           points to.

    CoordinateConstValueAccessor is a accessor that forwards
    the underlying accessor's operator() read functions.
    It passes its arguments <em>by value</em>.

    <b>\#include</b> \<vigra/coordinate_iterator.hxx\><br>
    Namespace: vigra
*/
template <class Accessor, class COORD>
class CoordinateConstValueAccessor
{
  public:
    typedef typename Accessor::value_type               forward_type;
    typedef AccessorCoordinatePair<forward_type, COORD> value_type;
    Accessor a;
    CoordinateConstValueAccessor(const Accessor & a_) : a(a_) {}
        /** Read the current data item.
        */
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const
    {
        const typename ITERATOR::value_type & x = *i;
        return value_type(a(&x.value()), x.coord()); 
    }
        /** Read the data item at an offset.
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE const & diff) const
    {
        const typename ITERATOR::value_type & x = i[diff];
        return value_type(a(&x.value()), x.coord());
    }
};

template<unsigned N, class T, class Accessor>
struct SrcCoordinateMultiArrayRangeAccessorReturns
{
    typedef          CombinedMultiIterator<N, T>  iterator_type;
    typedef typename iterator_type::coord_type    coord_type;
    typedef typename iterator_type::pair_type     pair_type;
    typedef typename iterator_type::ptr_type      ptr_type;
    typedef typename ptr_type::deref_type         deref_type;
    typedef typename iterator_type::stride_type   stride_type;
    typedef typename pair_type::stride_array_type stride_array_type;
    typedef typename pair_type::shape_array_type  shape_array_type;

    typedef CoordinateConstValueAccessor<Accessor, coord_type> access_type;
    typedef triple<iterator_type, stride_array_type, access_type> type;
};

template <unsigned int N, class T, class StrideTag, class Access>
inline typename SrcCoordinateMultiArrayRangeAccessorReturns<N, T, Access>::type
srcCoordinateMultiArrayRangeAccessor(const MultiArrayView<N, T, StrideTag> &
                                                                          array,
                                     Access a,
                                     typename CoordinateSteps<N>::type stride
                                                   = TinyVector<double, N>(1.0),
                                     typename CoordinateSteps<N>::type origin
                                                   = TinyVector<double, N>(0.0))
{
    typedef SrcCoordinateMultiArrayRangeAccessorReturns<N, T, Access> returns;
    typedef typename returns::type             type;
    typedef typename returns::stride_type      stride_type;
    typedef typename returns::access_type      access_type;
    typedef typename returns::iterator_type    iterator_type;
    typedef typename returns::shape_array_type shape_array_type;
    
    shape_array_type shape = array.shape();
    return type(iterator_type(array.traverser_begin().get(),
                              stride_type(origin),
                              array.stride(),
                              stride_type(stride),
                              shape),
                detail::dismember<N>(stride_type(shape)),
                access_type(a));
}

} // namespace vigra

namespace std {

template <unsigned N>
ostream &
operator<<(ostream & os, const vigra::StridePair<N> & x)
{
    os << "[" << x.index << ", " << x.coord << "]";
    return os;
}

template <unsigned N>
ostream &
operator<<(ostream & os, const vigra::StridePairDiff<N> & x)
{
    os << "<" << x.c << "; "
       << static_cast<vigra::StridePairCoord<N> >(x) << ">";
    return os;
}

template <unsigned N, class T>
ostream &
operator<<(ostream & os, const vigra::StridePairPointer<N, T> & x)
{
    os << "[" << x.value() << ", " << x.coord() << "]";
    return os;
}

template <class VALUETYPE, class COORD>
ostream &
operator<<(ostream & os,
           const vigra::AccessorCoordinatePair<VALUETYPE, COORD> & x)
{
    os << "[" << x.value() << ", " << x.coord() << "]";
    return os;
}

} // namespace std

#endif // VIGRA_COORDINATE_ITERATOR_HXX
