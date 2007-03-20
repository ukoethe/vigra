/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2007 by Kasim Terzic, Christian-Dennis Rahn       */
/*                        and Ullrich Koethe                            */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_MULTI_MORPHOLOGY_HXX
#define VIGRA_MULTI_MORPHOLOGY_HXX

#include <vector>
#include <cmath>
#include "multi_distance.hxx"
#include "array_vector.hxx"
#include "multi_array.hxx"
#include "accessor.hxx"
#include "numerictraits.hxx"
#include "navigator.hxx"
#include "metaprogramming.hxx"
#include "multi_pointoperators.hxx"
#include "functorexpression.hxx"

namespace vigra
{

/** \addtogroup MultiArrayMorphology Morphological operators for multi-dimensional arrays.

    These functions perform morphological operations on an arbitrary
    dimensional array that is specified by iterators (compatible to \ref MultiIteratorPage)
    and shape objects. It can therefore be applied to a wide range of data structures
    (\ref vigra::MultiArrayView, \ref vigra::MultiArray etc.).
*/
//@{

/********************************************************/
/*                                                      */
/*             multiBinaryErosion                       */
/*                                                      */
/********************************************************/
/** \brief Binary erosion on multi-dimensional arrays.

    This function applies a flat circular erosion operator with a given radius. The
    operation is isotropic.
    The input is a binary multi-dimensional array where non-zero pixels represent 
    foreground and zero pixels represent background.
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.
           
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryErosion(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    int radius);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__morphology_8hxx-source.html">vigra/multi_morphology.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned char> dest(shape);
    ...

    // perform isotropic binary erosion
    multiBinaryErosion(srcMultiArrayRange(source), destMultiArray(dest), 3);
    \endcode

    \see vigra::discErosion()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiBinaryErosion( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, float radius)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MaxValue = NumericTraits<DestType>::max();
    float radius2 = (float) radius * radius;
    enum { N = 1 + SrcIterator::level };
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
   
    using namespace vigra::functor;

    // Get the distance squared transform of the image
    if(N*MaxDim*MaxDim > MaxValue)
    {
        // Allocate a new temporary array if the distances squared wouldn't fit
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        //detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
        //    typename AccessorTraits<TmpType>::default_accessor()/*, false*/ );
        
        separableMultiDistSquared(s, shape, src, tmpArray.traverser_begin(),
                typename AccessorTraits<TmpType>::default_accessor(), false );
        
        // threshold everything less than radius away from the edge
        // std::cerr << "Thresholding!!!!!" << std::endl;
        transformMultiArray( tmpArray.traverser_begin(), shape, 
            typename AccessorTraits<TmpType>::default_accessor(), d, dest, 
            ifThenElse( Arg1() > Param(radius2),
                Param(MaxValue), Param(0) ) );
    }
    else    // work directly on the destination array
    {
        //detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest/*, false*/ );
        separableMultiDistSquared( s, shape, src, d, dest, false );
        
        // threshold everything less than radius away from the edge
        transformMultiArray( d, shape, dest, d, dest, 
            ifThenElse( Arg1() > Param(radius2),
                Param(MaxValue), Param(0) ) );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void multiBinaryErosion(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, int radius)
{
    multiBinaryErosion( source.first, source.second, source.third,
                                 dest.first, dest.second, radius );
}


/********************************************************/
/*                                                      */
/*             multiBinaryDilation                      */
/*                                                      */
/********************************************************/

/** \brief Binary dilation on multi-dimensional arrays.

    This function applies a flat circular dilation operator with a given radius. The
    operation is isotropic.
    The input is a binary multi-dimensional array where non-zero pixels represent
    foreground and zero pixels represent background.
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.
           
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryDilation(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    int radius);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__morphology_8hxx-source.html">vigra/multi_morphology.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned char> dest(shape);
    ...

    // perform isotropic binary erosion
    multiBinaryDilation(srcMultiArrayRange(source), destMultiArray(dest), 3);
    \endcode

    \see vigra::discDilation()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiBinaryDilation( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, float radius)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MaxValue = NumericTraits<DestType>::max();
    float radius2 = (float) radius * radius;
    enum { N = 1 + SrcIterator::level };
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
   
    using namespace vigra::functor;

    // Get the distance squared transform of the image
    if(N*MaxDim*MaxDim > MaxValue)
    {
        // Allocate a new temporary array if the distances squared wouldn't fit
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        //detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
        //    typename AccessorTraits<TmpType>::default_accessor(), true );
        
        separableMultiDistSquared(s, shape, src, tmpArray.traverser_begin(),
            typename AccessorTraits<TmpType>::default_accessor(), true );
      
        // threshold everything less than radius away from the edge
        transformMultiArray( tmpArray.traverser_begin(), shape, 
            typename AccessorTraits<TmpType>::default_accessor(), d, dest, 
            ifThenElse( Arg1() > Param(radius2),
                Param(0), Param(MaxValue) ) );
    }
    else    // work directly on the destination array
    {
        //detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest, true );
        separableMultiDistSquared( s, shape, src, d, dest, true );
        
        // threshold everything less than radius away from the edge
        transformMultiArray( d, shape, dest, d, dest, 
            ifThenElse( Arg1() > Param(radius2),
                Param(0), Param(MaxValue) ) );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void multiBinaryDilation(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, int radius)
{
    multiBinaryDilation( source.first, source.second, source.third,
                                 dest.first, dest.second, radius );
}

/********************************************************/
/*                                                      */
/*             multiGrayscaleErosion                    */
/*                                                      */
/********************************************************/
/** \brief Parabolic grayscale erosion on multi-dimensional arrays.

    This function applies a parabolic erosion operator with a given spread (sigma) on
    a grayscale array. The operation is isotropic.
    The input is a grayscale multi-dimensional array.
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.
           
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, float sigma);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleErosion(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    float sigma);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__morphology_8hxx-source.html">vigra/multi_morphology.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned char> dest(shape);
    ...

    // perform isotropic grayscale erosion
    multiGrayscaleErosion(srcMultiArrayRange(source), destMultiArray(dest), 3.0);
    \endcode

    \see vigra::discErosion()
*/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiGrayscaleErosion( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                       DestIterator d, DestAccessor dest, float sigma)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MaxValue = NumericTraits<DestType>::max();
    enum { N = 1 + SrcIterator::level };
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );
        
    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
    
    using namespace vigra::functor;
    
    // Allocate a new temporary array if the distances squared wouldn't fit
    if(N*MaxDim*MaxDim > MaxValue)
    {
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);

        detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
            typename AccessorTraits<TmpType>::default_accessor(), sigma );
        
        transformMultiArray( tmpArray.traverser_begin(), shape,
                typename AccessorTraits<TmpType>::default_accessor(), d, dest,
                ifThenElse( Arg1() > Param(MaxValue), Param(MaxValue), Arg1() ) );
        //copyMultiArray( tmpArray.traverser_begin(), shape,
        //        typename AccessorTraits<TmpType>::default_accessor(), d, dest );
    }
    else
    {
        detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest, sigma );
    }

}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void multiGrayscaleErosion(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, float sigma)
{
    multiGrayscaleErosion( source.first, source.second, source.third, 
            dest.first, dest.second, sigma);
}

/********************************************************/
/*                                                      */
/*             multiGrayscaleDilation                   */
/*                                                      */
/********************************************************/
/** \brief Parabolic grayscale dilation on multi-dimensional arrays.

    This function applies a parabolic dilation operator with a given spread (sigma) on
    a grayscale array. The operation is isotropic.
    The input is a grayscale multi-dimensional array.
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.
           
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, float sigma);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleDilation(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    float sigma);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__morphology_8hxx-source.html">vigra/multi_morphology.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned char> dest(shape);
    ...

    // perform isotropic grayscale erosion
    multiGrayscaleDilation(srcMultiArrayRange(source), destMultiArray(dest), 3.0);
    \endcode

    \see vigra::discErosion()
*/


template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void multiGrayscaleDilation( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, float sigma)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MinValue = NumericTraits<DestType>::min();
    DestType MaxValue = NumericTraits<DestType>::max();
    enum { N = 1 + SrcIterator::level };
        
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );
        
    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
    
    using namespace vigra::functor;

    // Allocate a new temporary array if the distances squared wouldn't fit
    if(-N*MaxDim*MaxDim < MinValue || N*MaxDim*MaxDim > MaxValue)
    {
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);

        detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
            typename AccessorTraits<TmpType>::default_accessor(), sigma, true );
        
        transformMultiArray( tmpArray.traverser_begin(), shape,
                typename AccessorTraits<TmpType>::default_accessor(), d, dest,
                ifThenElse( Arg1() > Param(MaxValue), Param(MaxValue), 
                    ifThenElse( Arg1() < Param(MinValue), Param(MinValue), Arg1() ) ) );
    }
    else
    {
        detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest, sigma, true );
    }

}


template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void multiGrayscaleDilation(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, float sigma)
{
    multiGrayscaleDilation( source.first, source.second, source.third, 
            dest.first, dest.second, sigma);
}


//@}

} //-- namespace vigra


#endif        //-- VIGRA_MULTI_MORPHOLOGY_HXX
