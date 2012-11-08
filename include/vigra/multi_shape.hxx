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


namespace vigra {

/** \addtogroup MultiIteratorGroup  Multi-dimensional Shapes and Array Iterators

    \brief Shape objects and general iterators for arrays of arbitrary dimension.
*/
//@{

    /** Index type for a single dimension of a MultiArrayView or
        MultiArray.
    */
typedef std::ptrdiff_t MultiArrayIndex;

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

enum NeighborhoodType { DirectNeighborhood=0, IndirectNeighborhood=1 };

// Helper functions

namespace detail {

/********************************************************/
/*                                                      */
/*                    defaultStride                     */
/*                                                      */
/********************************************************/

    /* generates the stride for a gapless shape.
    */
template <unsigned int N>
inline TinyVector <MultiArrayIndex, N>
defaultStride(const TinyVector <MultiArrayIndex, N> &shape)
{
    TinyVector <MultiArrayIndex, N> ret;
    ret [0] = 1;
    for (int i = 1; i < (int)N; ++i)
        ret [i] = ret [i-1] * shape [i-1];
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
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> &shape,
         const TinyVector <MultiArrayIndex, N> & coordinate)
    {
        return coordinate[N-K] + shape[N-K] * CoordinateToScanOrder<K-1>::exec(shape, coordinate);
    }
};

template <>
struct CoordinateToScanOrder<1>
{
    template <int N>
    static MultiArrayIndex
    exec(const TinyVector <MultiArrayIndex, N> & /*shape*/,
         const TinyVector <MultiArrayIndex, N> & coordinate)
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

#if 0
/********************************************************/
/*                                                      */
/*                MoveToScanOrderIndex                  */
/*                                                      */
/********************************************************/

    /* Move a Coupled.
    */
template <unsigned int M>
struct MoveToScanOrderIndex
{
    template <class Shape, class Ptr>
    static void 
    exec(MultiArrayIndex newIndex, Shape const & shape,
         Shape & point, Ptr & p, Shape const & strides)
    {
        enum { N = Shape::static_size };
        MultiArrayIndex newPos = newIndex % shape[N-1-M];
        p += (newPos - point[N-1-M]) * strides[N-1-M];
        point[N-1-M] = newPos;
        MoveToScanOrderIndex<M-1>::exec(newIndex / shape[N-1-M], shape, point, p, strides);
    }
    
    template <class Shape, class Ptr1, class Ptr2>
    static void 
    exec(MultiArrayIndex newIndex, Shape const & shape, Shape & point, 
         Ptr1 & p1, Shape const & strides1, Ptr2 & p2, Shape const & strides2)
    {
        enum { N = Shape::static_size };
        MultiArrayIndex newPos = newIndex % shape[N-1-M];
        p1 += (newPos - point[N-1-M]) * strides1[N-1-M];
        p2 += (newPos - point[N-1-M]) * strides2[N-1-M];
        point[N-1-M] = newPos;
        MoveToScanOrderIndex<M-1>::exec(newIndex / shape[N-1-M], shape, point, 
                                         p1, strides1, p2, strides2);
    }
};

template <>
struct MoveToScanOrderIndex<0>
{
    template <class Shape, class Ptr>
    static void 
    exec(MultiArrayIndex newIndex, Shape const & shape,
         Shape & point, Ptr & p, Shape const & strides)
    {
        enum { N = Shape::static_size }; 
        MultiArrayIndex newPos = std::min(newIndex, shape[N-1]);
        p += (newPos - point[N-1]) * strides[N-1];
        point[N-1] = newPos;
    }
    
    template <class Shape, class Ptr1, class Ptr2>
    static void 
    exec(MultiArrayIndex newIndex, Shape const & shape, Shape & point, 
         Ptr1 & p1, Shape const & strides1, Ptr2 & p2, Shape const & strides2)
    {
        enum { N = Shape::static_size }; 
        MultiArrayIndex newPos = std::min(newIndex, shape[N-1]);
        p1 += (newPos - point[N-1]) * strides1[N-1];
        p2 += (newPos - point[N-1]) * strides2[N-1];
        point[N-1] = newPos;
    }
};
#endif

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


template <class Shape>
void
makeArrayNeighborhood(ArrayVector<ArrayVector<Shape> > & neighborOffsets, 
                      ArrayVector<ArrayVector<bool> > & neighborExists,
                      ArrayVector<ArrayVector<bool> > & causalNeighborExists,
                      ArrayVector<ArrayVector<bool> > & anticausalNeighborExists,
                      ArrayVector<ArrayVector<int> > & neighborIndexLookup,
                      NeighborhoodType neighborhoodType = DirectNeighborhood)
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    Shape strides = cumprod(Shape(MultiArrayIndex(3))) / 3; 
    
    neighborOffsets.resize(size);
    neighborOffsets[0].clear(); // [0] is the standard case of all neighbors present
    if(neighborhoodType == DirectNeighborhood)
    {
        MakeDirectArrayNeighborhood<N-1>::offsets(neighborOffsets[0]);
    }
    else
    {
        Shape point; // represents the center
        MakeIndirectArrayNeighborhood<N-1>::offsets(neighborOffsets[0], point);
    }
    
    unsigned int neighborCount = neighborOffsets[0].size(); // maximal number of neighbors

#ifdef VERBOSE    
    std::cerr << " size " << neighborCount << ": " << neighborOffsets[0] << "\n strides ";
    for(unsigned int l=0; l<neighborCount; ++l)
        std::cerr << dot(neighborOffsets[0][l], strides) << ", ";
    std::cerr << "\n\n";
#endif
    
    neighborExists.resize(size);
    causalNeighborExists.resize(size);
    anticausalNeighborExists.resize(size);
    neighborIndexLookup.resize(size);

    for(unsigned int k=0; k<size; ++k) // iterate all k neighborhood codes
    {
        if (k>0) 
            neighborOffsets[k].clear();
        neighborExists[k].clear();
        if(neighborhoodType == DirectNeighborhood)
        {
            MakeDirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        else
        {
            MakeIndirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        
        causalNeighborExists[k].resize(neighborCount);
        anticausalNeighborExists[k].resize(neighborCount);
        
        for(unsigned int l = 0; l<neighborCount; ++l)
        {
            MultiArrayIndex stride = dot(neighborOffsets[0][l], strides);
            if(stride < 0)
            {
                causalNeighborExists[k][l] = neighborExists[k][l];
                anticausalNeighborExists[k][l] = false;
            }
            else
            {
                causalNeighborExists[k][l] = false;
                anticausalNeighborExists[k][l] = neighborExists[k][l];
            }
            if (neighborExists[k][l])
                neighborIndexLookup[k].push_back(l);
            if (k>0)
                if (neighborExists[k][l])
                    neighborOffsets[k].push_back(neighborOffsets[0][l]);
        }
    }

}

template <class Shape>
void
makeArraySubNeighborhood(const ArrayVector<Shape> & allNeighborOffsets, 
             const ArrayVector<ArrayVector<bool> > & neighborExists,
             const Shape strides,
             ArrayVector<ArrayVector<MultiArrayIndex> > & neighborIndices
             )
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    
    neighborIndices.resize(size);
    const unsigned int neighborCount = allNeighborOffsets.size(); // maximal number of neighbors
    
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
    for(unsigned int l=0; l<neighborCount; ++l) 
        if (neighborExists[k][l])
        neighborIndices[k].push_back(dot(allNeighborOffsets[l], strides));
#if 0
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
    {
    std::cerr << " NB-type " << k << ": ";
    for(unsigned int l=0; l<neighborCount; ++l) 
        if (neighborExists[k][l])
        {
        std::cerr << neighborIndices[k].back() << ", ";
        }
    std::cerr << std::endl;
    }
#endif
}

} // namespace detail

//@}

} // namespace vigra

#endif // VIGRA_MULTI_SHAPE_HXX
