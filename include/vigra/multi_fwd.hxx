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

#ifndef VIGRA_MULTI_FWD_HXX
#define VIGRA_MULTI_FWD_HXX

#include <memory>
#include "metaprogramming.hxx"
#include "tinyvector.hxx"

namespace vigra {

/** \addtogroup MultiIteratorGroup  Multi-dimensional Shapes and Array Iterators

    \brief Shape objects and general iterators for arrays of arbitrary dimension.
*/
//@{

/********************************************************/
/*                                                      */
/*                        shape                         */
/*                                                      */
/********************************************************/

    /** Index type for a single dimension of a MultiArrayView or
        MultiArray.
    */
typedef std::ptrdiff_t MultiArrayIndex;

template <unsigned int N>
class MultiArrayShape;

/********************************************************/
/*                                                      */
/*                    memory layout                     */
/*                                                      */
/********************************************************/

// the resulting MultiArray has no explicit channel axis
// (i.e. the number of channels is implicitly one)
template <class T>
struct Singleband;

// the last axis is explicitly designated as channel axis
template <class T>
struct Multiband;

// the array is organised in chunks
template <class T>
struct ChunkedMemory;

// helper classes to map memory layout to actual type
namespace detail {

template <class T>
struct ResolveMultiband;

template <class T>
struct ResolveChunkedMemory;

} // namespace detail

/********************************************************/
/*                                                      */
/*                   element handles                    */
/*                                                      */
/********************************************************/

template <class T, class NEXT>
class CoupledHandle;

template <unsigned int N, class T>
class SharedChunkHandle;

/********************************************************/
/*                                                      */
/*                     array types                      */
/*                                                      */
/********************************************************/

template <unsigned int N, class T, class C = StridedArrayTag>
class MultiArrayView;

template <unsigned int N, class T, 
          class A = std::allocator<typename detail::ResolveMultiband<T>::type> >
class MultiArray;

template <unsigned int N, class T>
class ChunkedArrayBase;

template <unsigned int N, class T>
class ChunkedArray;

/********************************************************/
/*                                                      */
/*                iterators / traversers                */
/*                                                      */
/********************************************************/

// DEPRECATED!
// hierarchical iterator whose innermost dimension is unstrided
template <unsigned int N, class T, class REFERENCE = T &, class POINTER = T *>
class MultiIterator;

// DEPRECATED!
// hierarchical iterator whose innermost dimension may be strided
template <unsigned int N, class T, class REFERENCE = T &, class POINTER = T *>
class StridedMultiIterator;

// multi-dimensional scan-order iteration returning the coordinate of the current element
template<unsigned int N>
class MultiCoordinateIterator;

// scan-order iteration returning the current element of a multi-dimensional array
template <unsigned int N, class V, class REFERENCE, class POINTER>
class StridedScanOrderIterator;

// simultaneous scan-order iteration over several multi-dimensional arrays
template <unsigned int N,
          class HANDLES=CoupledHandle<TinyVector<MultiArrayIndex, N>, void>,
          int DIMENSION = N-1>
class CoupledScanOrderIterator;

// scan-order iteration over the chunks of a chunked array
// returns a MultiArrayView for the current chunk
template<unsigned int N, class T>
class ChunkIterator;

/********************************************************/
/*                                                      */
/*                 neighborhood types                   */
/*                                                      */
/********************************************************/

    /** \brief Choose the neighborhood system in a dimension-independent way.  
    
        DirectNeighborhood corresponds to 4-neighborhood in 2D and 6-neighborhood in 3D, whereas
        IndirectNeighborhood means 8-neighborhood in 2D and 26-neighborhood in 3D. The general
        formula for N dimensions are 2*N for direct neighborhood and 3^N-1 for indirect neighborhood. 
    */
enum NeighborhoodType { 
        DirectNeighborhood=0,   ///< use only direct neighbors
        IndirectNeighborhood=1  ///< use direct and indirect neighbors
};


/********************************************************/
/*                                                      */
/*                     grid graph                       */
/*                                                      */
/********************************************************/

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphNeighborIterator;

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphEdgeIterator;

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphOutEdgeIterator;

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphArcIterator;

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphOutArcIterator;

template<unsigned int N, bool BackEdgesOnly=false>
class GridGraphInArcIterator;

template<unsigned int N, class DirectedTag>
class GridGraph;

//@}

} // namespace vigra

#endif // VIGRA_MULTI_FWD_HXX
