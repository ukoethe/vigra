/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2007 by Kasim Terzic, Christian-Dennis Rahn       */
/*                        and Ullrich Koethe                            */
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

namespace detail {

// this class simplifies the design, but more importantly, it makes sure
// that the in-place code doesn't get compiled for boolean arrays 
// (were it would never executed anyway -- see the specializations below)
template <class DestType, class TmpType>
struct MultiBinaryMorphologyImpl
{
    template <class SrcIterator, class SrcShape, class SrcAccessor,
              class DestIterator, class DestAccessor>
    static void
    exec( SrcIterator s, SrcShape const & shape, SrcAccessor src,
          DestIterator d, DestAccessor dest, 
          double radius, bool dilation)
    {
        using namespace vigra::functor;
        
        // Allocate a new temporary array if the distances squared wouldn't fit
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
            
        separableMultiDistSquared(s, shape, src, 
                                  tmpArray.traverser_begin(), typename AccessorTraits<TmpType>::default_accessor(), dilation );
            
        // threshold everything less than radius away from the edge
        double radius2 = radius * radius;
        DestType foreground = dilation 
                                 ? NumericTraits<DestType>::zero()
                                 : NumericTraits<DestType>::one(),
                 background = dilation 
                                 ? NumericTraits<DestType>::one()
                                 : NumericTraits<DestType>::zero();
        transformMultiArray( tmpArray.traverser_begin(), shape, StandardValueAccessor<double>(), 
                             d, dest, 
                             ifThenElse( Arg1() >= Param(radius2),
                                         Param(foreground), Param(background) ) );
    }
};

template <class DestType>
struct MultiBinaryMorphologyImpl<DestType, DestType>
{
    template <class SrcIterator, class SrcShape, class SrcAccessor,
              class DestIterator, class DestAccessor>
    static void
    exec( SrcIterator s, SrcShape const & shape, SrcAccessor src,
          DestIterator d, DestAccessor dest, 
          double radius, bool dilation)
    {
        using namespace vigra::functor;

        separableMultiDistSquared( s, shape, src, d, dest, dilation );
        
        // threshold everything less than radius away from the edge
        DestType radius2 = detail::RequiresExplicitCast<DestType>::cast(radius * radius);
        DestType foreground = dilation 
                                 ? NumericTraits<DestType>::zero()
                                 : NumericTraits<DestType>::one(),
                 background = dilation 
                                 ? NumericTraits<DestType>::one()
                                 : NumericTraits<DestType>::zero();
        transformMultiArray( d, shape, dest, d, dest, 
                             ifThenElse( Arg1() > Param(radius2),
                                         Param(foreground), Param(background) ) );
    }
};

template <>
struct MultiBinaryMorphologyImpl<bool, bool>
{
    template <class SrcIterator, class SrcShape, class SrcAccessor,
              class DestIterator, class DestAccessor>
    static void
    exec( SrcIterator s, SrcShape const & shape, SrcAccessor src,
          DestIterator d, DestAccessor dest, double radius, bool dilation)
    {
        vigra_fail("multiBinaryMorphology(): Internal error (this function should never be called).");
    }
};

} // namespace detail

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
    operation is isotropic. The input is interpreted as a binary multi-dimensional 
    array where non-zero pixels represent foreground and zero pixels represent 
    background. In the output, foreground is always represented by ones 
    (i.e. NumericTrais<typename DestAccessor::value_type>::one()).
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A temporary internal array is only allocated if working on the destination
    array directly would cause overflow errors (that is if
    <tt> NumericTraits<typename DestAccessor::value_type>::max() < squaredNorm(shape)</tt>, 
    i.e. the squared length of the image diagonal doesn't fit into the destination type).
           
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
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

    <b>\#include</b> \<vigra/multi_morphology.hxx\>

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
doxygen_overloaded_function(template <...> void multiBinaryErosion)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiBinaryErosion( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, double radius)
{
    typedef typename DestAccessor::value_type DestType;
    typedef Int32 TmpType;
    
    double dmax = squaredNorm(shape);

    // Get the distance squared transform of the image
    if(dmax > NumericTraits<DestType>::toRealPromote(NumericTraits<DestType>::max()))
    {
        detail::MultiBinaryMorphologyImpl<DestType, TmpType>::exec(s, shape, src, d, dest, radius, false);
    }
    else    // work directly on the destination array
    {
        detail::MultiBinaryMorphologyImpl<DestType, DestType>::exec(s, shape, src, d, dest, radius, false);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
multiBinaryErosion(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                   pair<DestIterator, DestAccessor> const & dest, double radius)
{
    multiBinaryErosion( source.first, source.second, source.third,
                        dest.first, dest.second, radius );
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void
multiBinaryErosion(MultiArrayView<N, T1, S1> const & source,
                   MultiArrayView<N, T2, S2> dest, double radius)
{
    multiBinaryErosion( srcMultiArrayRange(source),
                        destMultiArray(dest), radius );
}

/********************************************************/
/*                                                      */
/*             multiBinaryDilation                      */
/*                                                      */
/********************************************************/

/** \brief Binary dilation on multi-dimensional arrays.

    This function applies a flat circular dilation operator with a given radius. The
    operation is isotropic. The input is interpreted as a binary multi-dimensional 
    array where non-zero pixels represent foreground and zero pixels represent 
    background. In the output, foreground is always represented by ones 
    (i.e. NumericTrais<typename DestAccessor::value_type>::one()).
    
    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A temporary internal array is only allocated if working on the destination
    array directly would cause overflow errors (that is if
    <tt> NumericTraits<typename DestAccessor::value_type>::max() < squaredNorm(shape)</tt>, 
    i.e. the squared length of the image diagonal doesn't fit into the destination type).
           
    <b> Declarations:</b>

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiBinaryDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, int radius);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
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

    <b>\#include</b> \<vigra/multi_morphology.hxx\>

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
doxygen_overloaded_function(template <...> void multiBinaryDilation)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiBinaryDilation( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, double radius)
{
    typedef typename DestAccessor::value_type DestType;
    typedef Int32 TmpType;
    
    double dmax = squaredNorm(shape);

    // Get the distance squared transform of the image
    if(dmax > NumericTraits<DestType>::toRealPromote(NumericTraits<DestType>::max()))
    {
        detail::MultiBinaryMorphologyImpl<DestType, TmpType>::exec(s, shape, src, d, dest, radius, true);
    }
    else    // work directly on the destination array
    {
        detail::MultiBinaryMorphologyImpl<DestType, DestType>::exec(s, shape, src, d, dest, radius, true);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
multiBinaryDilation(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                    pair<DestIterator, DestAccessor> const & dest, double radius)
{
    multiBinaryDilation( source.first, source.second, source.third,
                         dest.first, dest.second, radius );
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void 
multiBinaryDilation(MultiArrayView<N, T1, S1> const & source,
                    MultiArrayView<N, T2, S2> dest, double radius)
{
    multiBinaryDilation( srcMultiArrayRange(source),
                         destMultiArray(dest), radius );
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

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, double sigma);

    }
    \endcode

    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleErosion(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, double sigma);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleErosion(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    double sigma);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_morphology.hxx\>

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
doxygen_overloaded_function(template <...> void multiGrayscaleErosion)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
multiGrayscaleErosion( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                       DestIterator d, DestAccessor dest, double sigma)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MaxValue = NumericTraits<DestType>::max();
    enum { N = 1 + SrcIterator::level };
    
    // temporary array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );
        
    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
    
    using namespace vigra::functor;
    
    ArrayVector<double> sigmas(shape.size(), sigma);
    
    // Allocate a new temporary array if the distances squared wouldn't fit
    if(N*MaxDim*MaxDim > MaxValue)
    {
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);

        detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
            typename AccessorTraits<TmpType>::default_accessor(), sigmas );
        
        transformMultiArray( tmpArray.traverser_begin(), shape,
                typename AccessorTraits<TmpType>::default_accessor(), d, dest,
                ifThenElse( Arg1() > Param(MaxValue), Param(MaxValue), Arg1() ) );
        //copyMultiArray( tmpArray.traverser_begin(), shape,
        //        typename AccessorTraits<TmpType>::default_accessor(), d, dest );
    }
    else
    {
        detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest, sigmas );
    }

}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
multiGrayscaleErosion(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                      pair<DestIterator, DestAccessor> const & dest, double sigma)
{
    multiGrayscaleErosion( source.first, source.second, source.third, 
                           dest.first, dest.second, sigma);
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void
multiGrayscaleErosion(MultiArrayView<N, T1, S1> const & source,
                      MultiArrayView<N, T2, S2> dest, double sigma)
{
    multiGrayscaleErosion( srcMultiArrayRange(source), 
                           destMultiArray(dest), sigma);
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

    pass arbitrary-dimensional array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, double sigma);

    }
    \endcode

    pass \ref MultiIteratorPage "MultiIterators" and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleDilation(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, double sigma);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        multiGrayscaleDilation(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, 
                                    double sigma);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_morphology.hxx\>

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
doxygen_overloaded_function(template <...> void multiGrayscaleDilation)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void multiGrayscaleDilation( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, double sigma)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MinValue = NumericTraits<DestType>::min();
    DestType MaxValue = NumericTraits<DestType>::max();
    enum { N = 1 + SrcIterator::level };
        
    // temporary array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );
        
    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;
    
    int MaxDim = 0; 
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
    
    using namespace vigra::functor;

    ArrayVector<double> sigmas(shape.size(), sigma);

    // Allocate a new temporary array if the distances squared wouldn't fit
    if(-N*MaxDim*MaxDim < MinValue || N*MaxDim*MaxDim > MaxValue)
    {
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);

        detail::internalSeparableMultiArrayDistTmp( s, shape, src, tmpArray.traverser_begin(),
            typename AccessorTraits<TmpType>::default_accessor(), sigmas, true );
        
        transformMultiArray( tmpArray.traverser_begin(), shape,
                typename AccessorTraits<TmpType>::default_accessor(), d, dest,
                ifThenElse( Arg1() > Param(MaxValue), Param(MaxValue), 
                    ifThenElse( Arg1() < Param(MinValue), Param(MinValue), Arg1() ) ) );
    }
    else
    {
        detail::internalSeparableMultiArrayDistTmp( s, shape, src, d, dest, sigmas, true );
    }

}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
multiGrayscaleDilation(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                       pair<DestIterator, DestAccessor> const & dest, double sigma)
{
    multiGrayscaleDilation( source.first, source.second, source.third, 
                            dest.first, dest.second, sigma);
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void
multiGrayscaleDilation(MultiArrayView<N, T1, S1> const & source,
                       MultiArrayView<N, T2, S2> dest, double sigma)
{
    multiGrayscaleDilation( srcMultiArrayRange(source), 
                            destMultiArray(dest), sigma);
}

//@}

} //-- namespace vigra


#endif        //-- VIGRA_MULTI_MORPHOLOGY_HXX
