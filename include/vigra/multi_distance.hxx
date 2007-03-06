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

#ifndef VIGRA_MULTI_DISTANCE_HXX
#define VIGRA_MULTI_DISTANCE_HXX

#include <vector>
#include <functional>
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


namespace detail
{

/********************************************************/
/*                                                      */
/*                distParabola                          */
/*                                                      */
/*  Version with sigma (parabola spread) for morphology */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor >
void distParabola(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                              DestIterator id, DestAccessor da, float sigma )
{
    // We assume that the data in the input is distance squared and treat it as such
    int w = std::distance( is, iend );
    
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType ValueType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote SumType;
    
    // Define the stack we use to determine the nearest background row 
    // (from previous dimension), the items on the stack will separate this column into
    // separate regions of influence. Each region of influence is closest to the same 
    // background row from the previous dimension.
    typedef triple<int, ValueType, int> influence;
    std::vector<influence> _stack;

    SrcIterator ibegin = is;
    _stack.push_back(influence(0, *is, w));
    
    ++is;
    int current = 1;
    
    int y0, y1, y2, y_dash, delta_y;
    sigma = sigma * sigma;
    bool nosigma = closeAtTolerance( sigma, 1.0 );
    
    y0 = 0;   // The beginning of the influence of row y1
   
    while( is != iend && current < w )
    {
        y1 = _stack.back().first;
        y2 = current;
        delta_y = y2 - y1;

        // If sigma is 1 (common case) avoid float multiplication here.
        if(nosigma)
            y_dash = (*is - _stack.back().second) - delta_y*delta_y;
        else
            y_dash = (int)(sigma * (*is - _stack.back().second)) - delta_y*delta_y;
        y_dash = y_dash / (delta_y + delta_y);
        y_dash += y2;

        if( y_dash > y0)      
        {
            if( y_dash <= w )   // CASE 2 -- A new region of influence
            {
                y0 = y_dash;
                
                _stack.back().third = y_dash; 
                
                _stack.push_back(influence(current, *is, w));
            }

            // CASE 1 -- This parabola is never active
            ++is;
            ++current;
            continue;
        } 
        else    // CASE 3 -- Parabola shadows the previous one completely
        {
            _stack.pop_back();

            if(_stack.size() < 2) 
                y0=0;
            else    
                y0=_stack[_stack.size()-2].third;
            
            if(_stack.empty())  // This row influences all previous rows.
            {
                _stack.push_back(influence(current, *is, w));

                ++is;
                ++current;
                continue;
            }
        }
    }

    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.

    typename std::vector<influence>::iterator it = _stack.begin();

    int distance = 0;   // The distance squared
    current = 0;
    delta_y = 0;
    is = ibegin;

    for(; is != iend; ++current, ++id, ++is)
    {
        // FIXME FIXME Bound checking incorrect here? vvv
        if( current >= (*it).third && it != _stack.end()) ++it; 
       
        // FIXME FIXME The following check speeds things up for distance, but completely
        // messes up the grayscale morphology. Use an extra flag???
  /*      if( *is == 0 ) // Skip background pixels
        {
            *id = 0;
            continue;
        }
  */      
        delta_y = current - (*it).first;
        distance = delta_y * delta_y + (*it).second;
        *id = distance;
    }

}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void distParabola(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest, float sigma )
{
    distParabola(src.first, src.second, src.third,
                 dest.first, dest.second, sigma);
}

/********************************************************/
/*                                                      */
/*        internalSeparableMultiArrayDistTmp            */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void internalSeparableMultiArrayDistTmp(
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, float sigma, bool invert)
{
    // Sigma is the spread of the parabolas and is only used for ND morphology. When
    // calculating the distance transform, it is set to 1
    enum { N = 1 + SrcIterator::level };

    // we need the Promote type here if we want to invert the image (dilation)
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    
    // temporary array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;
    
    
    // only operate on first dimension here
    SNavigator snav( si, shape, 0 );
    DNavigator dnav( di, shape, 0 );

    using namespace vigra::functor;

    for( ; snav.hasMore(); snav++, dnav++ )
    {
            // first copy source to temp for maximum cache efficiency
            // Invert the values if necessary. Only needed for grayscale morphology
            if(invert)
                transformLine( snav.begin(), snav.end(), src, tmp.begin(),
                               typename AccessorTraits<TmpType>::default_accessor(), 
                               Param(NumericTraits<TmpType>::zero())-Arg1());
            else
                copyLine( snav.begin(), snav.end(), src, tmp.begin(),
                          typename AccessorTraits<TmpType>::default_accessor() );

            detail::distParabola( srcIterRange(tmp.begin(), tmp.end(),
                          typename AccessorTraits<TmpType>::default_const_accessor()),
                          destIter( dnav.begin(), dest ), sigma );
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d )
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( dnav.begin(), dnav.end(), dest,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             detail::distParabola( srcIterRange(tmp.begin(), tmp.end(),
                           typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ), sigma );
        }
    }
    if(invert) transformMultiArray( di, shape, dest, di, dest, -Arg1());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void internalSeparableMultiArrayDistTmp( SrcIterator si, SrcShape const & shape, SrcAccessor src,
                                                DestIterator di, DestAccessor dest, float sigma)
{
    internalSeparableMultiArrayDistTmp( si, shape, src, di, dest, sigma, false );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void internalSeparableMultiArrayDistTmp( SrcIterator si, SrcShape const & shape, SrcAccessor src,
                                                DestIterator di, DestAccessor dest)
{
    internalSeparableMultiArrayDistTmp( si, shape, src, di, dest, 1, false );
}

} // namespace detail

/** \addtogroup MultiArrayDistanceTransform Euclidean distance transform for multi-dimensional arrays.

    These functions perform the Euclidean distance transform an arbitrary dimensional
    array that is specified by iterators (compatible to \ref MultiIteratorPage)
    and shape objects. It can therefore be applied to a wide range of data structures
    (\ref vigra::MultiArrayView, \ref vigra::MultiArray etc.).
*/
//@{

/********************************************************/
/*                                                      */
/*             separableMultiDistSquared                */
/*                                                      */
/********************************************************/

/** \brief Euclidean distance squared on multi-dimensional arrays.

    This function performs a Euclidean distance squared transform on the given
    multi-dimensional array. Both source and destination
    arrays are represented by iterators, shape objects and accessors.
    The destination array is required to already have the correct size.

    This function expects a mask as its source, where background pixels are 
    marked as zero, and non-background pixels as non-zero. If the parameter 
    <i>background</i> is true, then the squared distance of all background
    pixels to the nearest object is calculated. Otherwise, the distance of all
    object pixels to the nearest background pixel is calculated.

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistSquared(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, bool background);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistance(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    bool background);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__distance_8hxx-source.html">vigra/multi_distance.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned int> dest(shape);
    ...

    // Calculate Euclidean distance squared for all background pixels 
    separableMultiDistSquared(srcMultiArrayRange(source), destMultiArray(dest), true);
    \endcode

    \see vigra::distanceTransform()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void separableMultiDistSquared( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, bool background)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::ValueType DestType;
    typedef typename NumericTraits<typename DestAccessor::value_type>::Promote TmpType;
    DestType MaxValue = NumericTraits<DestType>::max();
    enum { N = 1 + SrcIterator::level };

    int MaxDim = 0;
    for( int i=0; i<N; i++)
        if(MaxDim < shape[i]) MaxDim = shape[i];
    int MaxDist = MaxDim*MaxDim;

    using namespace vigra::functor;
   
    if(N*MaxDim*MaxDim > MaxValue) // need a temporary array to avoid overflows
    {
        // Threshold the values so all objects have infinity value in the beginning
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        if(background == true)
            transformMultiArray( s, shape, src, tmpArray.traverser_begin(),
                                 typename AccessorTraits<TmpType>::default_accessor(),
                                 ifThenElse( Arg1() == Param(0), Param(MaxDist), Param(0) ));
        else
            transformMultiArray( s, shape, src, tmpArray.traverser_begin(),
                                 typename AccessorTraits<TmpType>::default_accessor(),
                                 ifThenElse( Arg1() != Param(0), Param(MaxDist), Param(0) ));
        
        detail::internalSeparableMultiArrayDistTmp( tmpArray.traverser_begin(), 
                shape, typename AccessorTraits<TmpType>::default_accessor(),
                tmpArray.traverser_begin(), 
                typename AccessorTraits<TmpType>::default_accessor());
        
        //copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
        transformMultiArray( tmpArray.traverser_begin(), shape,
                             typename AccessorTraits<TmpType>::default_accessor(), d, dest,
                             ifThenElse( Arg1() > Param(MaxValue), Param(MaxValue), Arg1() ) );
              
    }
    else        // work directly on the destination array    
    {
        // Threshold the values so all objects have infinity value in the beginning
        if(background == true)
            transformMultiArray( s, shape, src, d, dest,
                                 ifThenElse( Arg1() == Param(0), Param(MaxDist), Param(0) ));
        else
            transformMultiArray( s, shape, src, d, dest, 
                                 ifThenElse( Arg1() != Param(0), Param(MaxDist), Param(0) ));
     
        detail::internalSeparableMultiArrayDistTmp( d, shape, dest, d, dest);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void separableMultiDistSquared( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest, bool background)
{
    separableMultiDistSquared( source.first, source.second, source.third,
                               dest.first, dest.second, background );
}

/********************************************************/
/*                                                      */
/*             separableMultiDistance                   */
/*                                                      */
/********************************************************/

/** \brief Euclidean distance on multi-dimensional arrays.

    This function performs a Euclidean distance transform on the given
    multi-dimensional array. Both source and destination
    arrays are represented by iterators, shape objects and accessors.
    The destination array is required to already have the correct size.

    This function expects a mask as its source, where background pixels are 
    marked as zero, and non-background pixels as non-zero. If the parameter 
    <i>background</i> is true, then the squared distance of all background
    pixels to the nearest object is calculated. Otherwise, the distance of all
    object pixels to the nearest background pixel is calculated.

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> typeid(typename DestAccessor::value_type) < N * M*M</tt>, where M is the
    size of the largest dimension of the array.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistance(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest, bool background);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistance(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    bool background);

    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__distance_8hxx-source.html">vigra/multi_distance.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned float> dest(shape);
    ...

    // Calculate Euclidean distance squared for all background pixels 
    separableMultiDistance(srcMultiArrayRange(source), destMultiArray(dest), true);
    \endcode

    \see vigra::distanceTransform()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void separableMultiDistance( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, bool background)
{
    separableMultiDistSquared( s, shape, src, d, dest, background);
    
    // Finally, calculate the square root of the distances
    transformMultiArray( d, shape, dest, d, dest, (double(*)(double))&std::sqrt );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void separableMultiDistance( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, bool background)
{
    separableMultiDistance( source.first, source.second, source.third,
                            dest.first, dest.second, background );
}

//@}

} //-- namespace vigra


#endif        //-- VIGRA_MULTI_DISTANCE_HXX
