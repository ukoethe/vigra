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

#ifndef VIGRA_MULTI_SHAPE_HXX
#define VIGRA_MULTI_SHAPE_HXX

#include <sys/types.h>
#include "tinyvector.hxx"
#include "array_vector.hxx"
#include "numerictraits.hxx"

namespace vigra {

/** \addtogroup MultiIteratorGroup  Multi-dimensional Shapes and Array Iterators

    \brief Shape objects and general iterators for arrays of arbitrary dimension.
*/
//@{

    /** Index type for a single dimension of a MultiArrayView or
        MultiArray.
    */
typedef std::ptrdiff_t MultiArrayIndex;

/********************************************************/
/*                                                      */
/*              Singleband and Multiband                */
/*                                                      */
/********************************************************/

template <class T>
struct Singleband  // the resulting MultiArray has no explicit channel axis 
                   // (i.e. the number of channels is implicitly one)
{
    typedef T value_type;
};

template <class T>
struct Multiband  // the last axis is explicitly designated as channel axis
{
    typedef T value_type;
};

template<class T>
struct NumericTraits<Singleband<T> >
: public NumericTraits<T>
{};

template<class T>
struct NumericTraits<Multiband<T> >
{
    typedef Multiband<T> Type;
/*
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
*/
    typedef Type ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef typename NumericTraits<T>::isSigned isOrdered;
    typedef typename NumericTraits<T>::isSigned isComplex;
/*
    static signed char zero() { return 0; }
    static signed char one() { return 1; }
    static signed char nonZero() { return 1; }
    static signed char min() { return SCHAR_MIN; }
    static signed char max() { return SCHAR_MAX; }

#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = SCHAR_MIN, maxConst = SCHAR_MIN };
#else
    static const signed char minConst = SCHAR_MIN;
    static const signed char maxConst = SCHAR_MIN;
#endif

    static Promote toPromote(signed char v) { return v; }
    static RealPromote toRealPromote(signed char v) { return v; }
    static signed char fromPromote(Promote v) {
        return ((v < SCHAR_MIN) ? SCHAR_MIN : (v > SCHAR_MAX) ? SCHAR_MAX : v);
    }
    static signed char fromRealPromote(RealPromote v) {
        return ((v < 0.0)
                   ? ((v < (RealPromote)SCHAR_MIN)
                       ? SCHAR_MIN
                       : static_cast<signed char>(v - 0.5))
                   : (v > (RealPromote)SCHAR_MAX)
                       ? SCHAR_MAX
                       : static_cast<signed char>(v + 0.5));
    }
*/
};

namespace detail {

template <class T>
struct ResolveMultiband
{
    typedef T type;
    typedef StridedArrayTag Stride;
    static const bool value = false;

    template <int N>
    static TinyVector <MultiArrayIndex, N>
    defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
    {
        return vigra::detail::defaultStride(shape);
    }
};

template <class T>
struct ResolveMultiband<Singleband<T> >
{
    typedef T type;
    typedef StridedArrayTag Stride;
    static const bool value = false;

    template <int N>
    static TinyVector <MultiArrayIndex, N>
    defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
    {
        return vigra::detail::defaultStride(shape);
    }
};

template <class T>
struct ResolveMultiband<Multiband<T> >
{
    typedef T type;
    typedef StridedArrayTag Stride;
    static const bool value = true;

    template <int N>
    static TinyVector <MultiArrayIndex, N>
    defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
    {
        return vigra::detail::defaultMultibandStride(shape);
    }
};

} // namespace detail

template <unsigned int N, class T, class C = StridedArrayTag>
class MultiArrayView;

template <unsigned int N, class T, 
          class A = std::allocator<typename vigra::detail::ResolveMultiband<T>::type> >
class MultiArray;


    /** Traits class for the difference type of all MultiIterator, MultiArrayView, and
        MultiArray variants.
    */
template <unsigned int N>
class MultiArrayShape
{
  public:
        /** The difference type of all MultiIterator, MultiArrayView, and
            MultiArray variants.
        */
    typedef TinyVector<MultiArrayIndex, N> type;
};

typedef MultiArrayShape<1>::type Shape1; ///< shape type for MultiArray<1, T>
typedef MultiArrayShape<2>::type Shape2; ///< shape type for MultiArray<2, T>
typedef MultiArrayShape<3>::type Shape3; ///< shape type for MultiArray<3, T>
typedef MultiArrayShape<4>::type Shape4; ///< shape type for MultiArray<4, T>
typedef MultiArrayShape<5>::type Shape5; ///< shape type for MultiArray<5, T>

    /** \brief Choose the neighborhood system in a dimension-independent way.  
    
        DirectNeighborhood corresponds to 4-neighborhood in 2D and 6-neighborhood in 3D, whereas
        IndirectNeighborhood means 8-neighborhood in 2D and 26-neighborhood in 3D. The general
        formula for N dimensions are 2*N for direct neighborhood and 3^N-1 for indirect neighborhood. 
    */
enum NeighborhoodType { 
        DirectNeighborhood=0,   ///< use only direct neighbors
        IndirectNeighborhood=1  ///< use direct and indirect neighbors
};

// Helper functions

namespace detail {

/********************************************************/
/*                                                      */
/*                    defaultStride                     */
/*                                                      */
/********************************************************/

    /* generates the stride for a gapless shape.
    */
template <int N>
inline TinyVector <MultiArrayIndex, N>
defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
{
    TinyVector <MultiArrayIndex, N> ret;
    ret [0] = 1;
    for (int i = 1; i < (int)N; ++i)
        ret [i] = ret [i-1] * shape [i-1];
    return ret;
}

    /* generates the stride for a gapless shape.
    */
template <int N>
inline TinyVector <MultiArrayIndex, N>
defaultMultibandStride(const TinyVector <MultiArrayIndex, N> &shape)
{
    TinyVector <MultiArrayIndex, N> ret;
    ret [N-1] = 1;
    for (int i = 0; i < (int)N-1; ++i)
    {
        int j = (i + int(N - 1)) % N;
        ret [i] = ret [j] * shape [j];
    }
    return ret;
}

/********************************************************/
/*                                                      */
/*                CoordinateToScanOrder                 */
/*                                                      */
/********************************************************/

    /* Convert multi-dimensional index (i.e. a grid coordinate) to scan-order index.
    */
template <int K>
struct CoordinateToScanOrder
{
    template <int N, class D1, class D2, class D3, class D4>
    static MultiArrayIndex
    exec(const TinyVectorBase <MultiArrayIndex, N, D1, D2> &shape,
         const TinyVectorBase <MultiArrayIndex, N, D3, D4> & coordinate)
    {
        return coordinate[N-K] + shape[N-K] * CoordinateToScanOrder<K-1>::exec(shape, coordinate);
    }
};

template <>
struct CoordinateToScanOrder<1>
{
    template <int N, class D1, class D2, class D3, class D4>
    static MultiArrayIndex
    exec(const TinyVectorBase <MultiArrayIndex, N, D1, D2> & /*shape*/,
         const TinyVectorBase <MultiArrayIndex, N, D3, D4> & coordinate)
    {
        return coordinate[N-1];
    }
};

/********************************************************/
/*                                                      */
/*                ScanOrderToCoordinate                 */
/*                                                      */
/********************************************************/

    /* Convert scan-order index to multi-dimensional index (i.e. a grid coordinate).
    */
template <int K>
struct ScanOrderToCoordinate
{
    template <int N>
    static void
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> &shape,
         TinyVector <MultiArrayIndex, N> & result)
    {
        result[N-K] = (d % shape[N-K]);
        ScanOrderToCoordinate<K-1>::exec(d / shape[N-K], shape, result);
    }
};

template <>
struct ScanOrderToCoordinate<1>
{
    template <int N>
    static void
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> & /*shape*/,
         TinyVector <MultiArrayIndex, N> & result)
    {
        result[N-1] = d;
    }
};

/********************************************************/
/*                                                      */
/*                 ScanOrderToOffset                    */
/*                                                      */
/********************************************************/

    /* transforms an index in scan order sense to a pointer offset in a possibly
       strided, multi-dimensional array.
    */
template <int K>
struct ScanOrderToOffset
{
    template <int N>
    static MultiArrayIndex
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> &shape,
         const TinyVector <MultiArrayIndex, N> & stride)
    {
        return stride[N-K] * (d % shape[N-K]) +
               ScanOrderToOffset<K-1>::exec(d / shape[N-K], shape, stride);
    }
};

template <>
struct ScanOrderToOffset<1>
{
    template <int N>
    static MultiArrayIndex
    exec(MultiArrayIndex d, const TinyVector <MultiArrayIndex, N> & /*shape*/,
         const TinyVector <MultiArrayIndex, N> & stride)
    {
        return stride[N-1] * d;
    }
};

/********************************************************/
/*                                                      */
/*                 ScanOrderToOffset                    */
/*                                                      */
/********************************************************/

    /* transforms a multi-dimensional index (grid coordinate) to a pointer offset in a possibly
       strided, multi-dimensional array.
    */
template <class C>
struct CoordinatesToOffest
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x)
    {
        return stride[0] * x;
    }
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x, MultiArrayIndex y)
    {
        return stride[0] * x + stride[1] * y;
    }
};

template <>
struct CoordinatesToOffest<UnstridedArrayTag>
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & /*stride*/, MultiArrayIndex x)
    {
        return x;
    }
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & stride, MultiArrayIndex x, MultiArrayIndex y)
    {
        return x + stride[1] * y;
    }
};

/********************************************************/
/*                                                      */
/*            RelativeToAbsoluteCoordinate              */
/*                                                      */
/********************************************************/

    /* transforms a coordinate object with negative indices into the corresponding 
       'shape - abs(index)'.
    */
template <int M>
struct RelativeToAbsoluteCoordinate
{
    template <int N>
    static void
    exec(const TinyVector<MultiArrayIndex, N> & shape, TinyVector<MultiArrayIndex, N> & coord)
    {
        RelativeToAbsoluteCoordinate<M-1>::exec(shape, coord);
        if(coord[M] < 0)
            coord[M] += shape[M];
    }
};

template <>
struct RelativeToAbsoluteCoordinate<0>
{
    template <int N>
    static void
    exec(const TinyVector<MultiArrayIndex, N> & shape, TinyVector<MultiArrayIndex, N> & coord)
    {
        if(coord[0] < 0)
            coord[0] += shape[0];
    }
};

/********************************************************/
/*                                                      */
/*                   BorderTypeImpl                     */
/*                                                      */
/********************************************************/

// a border type is a compact bit-wise encoding of the fact that a 
// given coordinate is at the border of the ROI. Each border corresponds
// to one bit in the encoding, e.g. the left, right, top, bottom borders
// of a 2D image are represented by bits 0 to 3 respectively. 
// If a bit is set, the point in question is at the corresponding border.
// A code of all zeros therefore means that the point is in the interior
// of the ROI
template <unsigned int N, unsigned int DIMENSION=N-1>
struct BorderTypeImpl
{
    typedef typename MultiArrayShape<N>::type shape_type;
    
    static unsigned int exec(shape_type const & point, shape_type const & shape)
    {
        unsigned int res = BorderTypeImpl<N, DIMENSION-1>::exec(point, shape);
        if(point[DIMENSION] == 0)
            res |= (1 << 2*DIMENSION);
        if(point[DIMENSION] == shape[DIMENSION]-1)
            res |= (2 << 2*DIMENSION);
        return res;
    }
};

template <unsigned int N>
struct BorderTypeImpl<N, 0>
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

/********************************************************/
/*                                                      */
/*                makeArrayNeighborhood                 */
/*                                                      */
/********************************************************/

// Create the offsets to all direct neighbors, starting from the given Level (=dimension)
// and append them to the given array. The algorithm is designed so that the offsets are 
// sorted by ascending strides. This has two important consequences:
//  * The first half of the array contains the causal neighbors (negative strides),
//    the second half the anti-causal ones (positive strides), where 'causal' refers
//    to all scan-order predecessors of the center pixel, and 'anticausal' to its successors.
//  * For any neighbor k, its opposite (=point-reflected) neighbor is located at index
//    'N-1-k', where N is the total number of neighbors.
// The function 'exists' returns an array of flags that contains 'true' when the corresponding 
// neighbor is inside the ROI for the given borderType, 'false' otherwise. 
template <unsigned int Level>
struct MakeDirectArrayNeighborhood
{
    template <class Array>
    static void offsets(Array & a)
    {
        typedef typename Array::value_type Shape;
        
        Shape point;
        point[Level] = -1;
        a.push_back(point);
        MakeDirectArrayNeighborhood<Level-1>::offsets(a);
        point[Level] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int borderType)
    {
        a.push_back((borderType & (1 << 2*Level)) == 0);
        MakeDirectArrayNeighborhood<Level-1>::exists(a, borderType);
        a.push_back((borderType & (2 << 2*Level)) == 0);
    }
};

template <>
struct MakeDirectArrayNeighborhood<0>
{
    template <class Array>
    static void offsets(Array & a)
    {
        typedef typename Array::value_type Shape;
        
        Shape point;
        point[0] = -1;
        a.push_back(point);
        point[0] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int borderType)
    {
        a.push_back((borderType & 1) == 0);
        a.push_back((borderType & 2) == 0);
    }
};

// Likewise, create the offsets to all indirect neighbors according to the same rules.
template <unsigned int Level>
struct MakeIndirectArrayNeighborhood
{
    template <class Array, class Shape>
    static void offsets(Array & a, Shape point, bool isCenter = true)
    {
        point[Level] = -1;
        MakeIndirectArrayNeighborhood<Level-1>::offsets(a, point, false);
        point[Level] = 0;
        MakeIndirectArrayNeighborhood<Level-1>::offsets(a, point, isCenter);
        point[Level] = 1;
        MakeIndirectArrayNeighborhood<Level-1>::offsets(a, point, false);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int borderType, bool isCenter = true)
    {
        if((borderType & (1 << 2*Level)) == 0)
            MakeIndirectArrayNeighborhood<Level-1>::exists(a, borderType, false);
        else 
            MakeIndirectArrayNeighborhood<Level-1>::markOutside(a);

        MakeIndirectArrayNeighborhood<Level-1>::exists(a, borderType, isCenter);

        if((borderType & (2 << 2*Level)) == 0)
            MakeIndirectArrayNeighborhood<Level-1>::exists(a, borderType, false);
        else
            MakeIndirectArrayNeighborhood<Level-1>::markOutside(a);
    }

    template <class Array>
    static void markOutside(Array & a)
    {
        // Call markOutside() three times, for each possible offset at (Level-1)
        MakeIndirectArrayNeighborhood<Level-1>::markOutside(a);
        MakeIndirectArrayNeighborhood<Level-1>::markOutside(a);
        MakeIndirectArrayNeighborhood<Level-1>::markOutside(a);
    }

};

template <>
struct MakeIndirectArrayNeighborhood<0>
{
    template <class Array, class Shape>
    static void offsets(Array & a, Shape point, bool isCenter = true)
    {
        point[0] = -1;
        a.push_back(point);
        if(!isCenter) // the center point is not a neighbor, it's just convenient to do the enumeration this way...
        {
            point[0] = 0;
            a.push_back(point);
        }
        point[0] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int borderType, bool isCenter = true)
    {
        a.push_back((borderType & 1) == 0);
        if(!isCenter)
        {
            a.push_back(true);
        }
        a.push_back((borderType & 2) == 0);
    }

    template <class Array>
    static void markOutside(Array & a)
    {
        // Push 'false' three times, for each possible offset at level 0, whenever the point was 
        // outside the ROI in one of the higher levels.
        a.push_back(false);
        a.push_back(false);
        a.push_back(false);
    }
};

// Create the list of neighbor offsets for the given neighborhood type 
// and dimension (the dimension is implicitly defined by the Shape type)
// an return it in 'neighborOffsets'. Moreover, create a list of flags
// for each BorderType that is 'true' when the corresponding neighbor exists
// in this border situation and return the result in 'neighborExists'.
template <class Shape>
void
makeArrayNeighborhood(ArrayVector<Shape> & neighborOffsets, 
                      ArrayVector<ArrayVector<bool> > & neighborExists,
                      NeighborhoodType neighborhoodType = DirectNeighborhood)
{
    enum { N = Shape::static_size };
    
    neighborOffsets.clear();
    if(neighborhoodType == DirectNeighborhood)
    {
        MakeDirectArrayNeighborhood<N-1>::offsets(neighborOffsets);
    }
    else
    {
        Shape point; // represents the center
        MakeIndirectArrayNeighborhood<N-1>::offsets(neighborOffsets, point);
    }
    
    unsigned int borderTypeCount = 1 << 2*N;
    neighborExists.resize(borderTypeCount);

    for(unsigned int k=0; k<borderTypeCount; ++k)
    {
        neighborExists[k].clear();
        if(neighborhoodType == DirectNeighborhood)
        {
            MakeDirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        else
        {
            MakeIndirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
    }
}

} // namespace detail

//@}

} // namespace vigra

#endif // VIGRA_MULTI_SHAPE_HXX
