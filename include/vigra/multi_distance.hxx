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

template <class Value>
struct DistParabolaStackEntry
{
    double left, center, right;
    Value prevVal;
    
    DistParabolaStackEntry(Value const & p, double l, double c, double r)
    : left(l), center(c), right(r), prevVal(p)
    {}
};

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
                  DestIterator id, DestAccessor da, double sigma )
{
    // We assume that the data in the input is distance squared and treat it as such
    double w = iend - is;
    double sigma2 = sigma * sigma;
    double sigma22 = 2.0 * sigma2;
    
    typedef typename SrcAccessor::value_type SrcType;
    typedef DistParabolaStackEntry<SrcType> Influence;
    std::vector<Influence> _stack;
    _stack.push_back(Influence(sa(is), 0.0, 0.0, w));
    
    ++is;
    double current = 1.0;
    while(current < w )
    {
        Influence & s = _stack.back();
        double diff = current - s.center;
        double intersection = current + (sa(is) - s.prevVal - sigma2*sq(diff)) / (sigma22 * diff);
        
        if( intersection < s.left) // previous point has no influence
        {
            _stack.pop_back();
            if(_stack.empty())
                _stack.push_back(Influence(sa(is), 0.0, current, w));
            else
                continue; // try new top of stack without advancing current
        }
        else if(intersection < s.right)
        {
            s.right = intersection;
            _stack.push_back(Influence(sa(is), intersection, current, w));
        }
        ++is;
        ++current;
    }

    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(current = 0.0; current < w; ++current, ++id)
    {
        if( current >= it->right) 
            ++it; 
        da.set(sigma2 * sq(current - it->center) + it->prevVal, id);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void distParabola(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest, double sigma )
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
          class DestIterator, class DestAccessor, class Array>
void internalSeparableMultiArrayDistTmp(
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, Array const & sigmas, bool invert)
{
    // Sigma is the spread of the parabolas. It determines the structuring element size
    // for ND morphology. When calculating the distance transforms, sigma is usually set to 1,
    // unless one wants to account for anisotropic pixel pitch
    enum { N =  SrcShape::static_size};

    // we need the Promote type here if we want to invert the image (dilation)
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
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
                          destIter( dnav.begin(), dest ), sigmas[0] );
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
                           destIter( dnav.begin(), dest ), sigmas[d] );
        }
    }
    if(invert) transformMultiArray( di, shape, dest, di, dest, -Arg1());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void internalSeparableMultiArrayDistTmp( SrcIterator si, SrcShape const & shape, SrcAccessor src,
                                                DestIterator di, DestAccessor dest, Array const & sigmas)
{
    internalSeparableMultiArrayDistTmp( si, shape, src, di, dest, sigmas, false );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void internalSeparableMultiArrayDistTmp( SrcIterator si, SrcShape const & shape, SrcAccessor src,
                                                DestIterator di, DestAccessor dest)
{
    ArrayVector<double> sigmas(shape.size(), 1.0);
    internalSeparableMultiArrayDistTmp( si, shape, src, di, dest, sigmas, false );
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

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // explicitly specify pixel pitch for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class Array>
        void 
        separableMultiDistSquared( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                   DestIterator d, DestAccessor dest, 
                                   bool background,
                                   Array const & pixelPitch);
                                        
        // use default pixel pitch = 1.0 for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistSquared(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                  DestIterator diter, DestAccessor dest, 
                                  bool background);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // explicitly specify pixel pitch for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class Array>
        void 
        separableMultiDistSquared( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                   pair<DestIterator, DestAccessor> const & dest, 
                                   bool background,
                                   Array const & pixelPitch);
                                               
        // use default pixel pitch = 1.0 for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistSquared(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                  pair<DestIterator, DestAccessor> const & dest,
                                  bool background);

    }
    \endcode

    This function performs a squared Euclidean distance squared transform on the given
    multi-dimensional array. Both source and destination
    arrays are represented by iterators, shape objects and accessors.
    The destination array is required to already have the correct size.

    This function expects a mask as its source, where background pixels are 
    marked as zero, and non-background pixels as non-zero. If the parameter 
    <i>background</i> is true, then the squared distance of all background
    pixels to the nearest object is calculated. Otherwise, the distance of all
    object pixels to the nearest background pixel is calculated.
    
    Optionally, one can pass an array that specifies the pixel pitch in each direction. 
    This is necessary when the data have non-uniform resolution (as is common in confocal
    micriscopy, for example). 

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause overflow errors (i.e. if
    <tt> NumericTraits<typename DestAccessor::value_type>::max() < N * M*M</tt>, where M is the
    size of the largest dimension of the array.

    <b> Usage:</b>

    <b>\#include</b> \<<a href="multi__distance_8hxx-source.html">vigra/multi_distance.hxx</a>\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, unsigned int> dest(shape);
    ...

    // Calculate Euclidean distance squared for all background pixels 
    separableMultiDistSquared(srcMultiArrayRange(source), destMultiArray(dest), true);
    \endcode

    \see vigra::distanceTransform(), vigra::separableMultiDistance()
*/
doxygen_overloaded_function(template <...> void separableMultiDistSquared)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void separableMultiDistSquared( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, bool background,
                                Array const & pixelPitch)
{
    int N = shape.size();

    typedef typename DestAccessor::value_type DestType;
    typedef typename NumericTraits<DestType>::RealPromote Real;
    
    double dmax = 0.0;
    bool pixelPitchIsReal = false;
    for( int k=0; k<N; ++k)
    {
        if(int(pixelPitch[k]) != pixelPitch[k])
            pixelPitchIsReal = true;
        dmax += sq(pixelPitch[k]*shape[k]);
    }
            
    using namespace vigra::functor;
   
    if(dmax > NumericTraits<DestType>::toPromote(NumericTraits<DestType>::max()) || pixelPitchIsReal) // need a temporary array to avoid overflows
    {
        // Threshold the values so all objects have infinity value in the beginning
        double maxDist = dmax;
        MultiArray<SrcShape::static_size, Real> tmpArray(shape);
        if(background == true)
            transformMultiArray( s, shape, src, 
                                 tmpArray.traverser_begin(), typename AccessorTraits<Real>::default_accessor(),
                                 ifThenElse( Arg1() == Param(0), Param(maxDist), Param(0) ));
        else
            transformMultiArray( s, shape, src, 
                                 tmpArray.traverser_begin(), typename AccessorTraits<Real>::default_accessor(),
                                 ifThenElse( Arg1() != Param(0), Param(maxDist), Param(0) ));
        
        detail::internalSeparableMultiArrayDistTmp( tmpArray.traverser_begin(), 
                shape, typename AccessorTraits<Real>::default_accessor(),
                tmpArray.traverser_begin(), 
                typename AccessorTraits<Real>::default_accessor(), pixelPitch);
        
        copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
    }
    else        // work directly on the destination array    
    {
        // Threshold the values so all objects have infinity value in the beginning
        DestType maxDist = DestType(std::ceil(dmax));
        if(background == true)
            transformMultiArray( s, shape, src, d, dest,
                                 ifThenElse( Arg1() == Param(0), Param(maxDist), Param(0) ));
        else
            transformMultiArray( s, shape, src, d, dest, 
                                 ifThenElse( Arg1() != Param(0), Param(maxDist), Param(0) ));
     
        detail::internalSeparableMultiArrayDistTmp( d, shape, dest, d, dest, pixelPitch);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void separableMultiDistSquared( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest, bool background,
                                       Array const & pixelPitch)
{
    separableMultiDistSquared( source.first, source.second, source.third,
                               dest.first, dest.second, background, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void separableMultiDistSquared( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, bool background)
{
    ArrayVector<double> pixelPitch(shape.size(), 1.0);
    separableMultiDistSquared( s, shape, src, d, dest, background, pixelPitch );
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

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // explicitly specify pixel pitch for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class Array>
        void 
        separableMultiDistance( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, 
                                bool background,
                                Array const & pixelPitch);
                                        
        // use default pixel pitch = 1.0 for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistance(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                               DestIterator diter, DestAccessor dest, 
                               bool background);

    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // explicitly specify pixel pitch for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class Array>
        void 
        separableMultiDistance( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                pair<DestIterator, DestAccessor> const & dest, 
                                bool background,
                                Array const & pixelPitch);
                                               
        // use default pixel pitch = 1.0 for each coordinate
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        separableMultiDistance(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                               pair<DestIterator, DestAccessor> const & dest,
                               bool background);

    }
    \endcode

    This function performs a Euclidean distance transform on the given
    multi-dimensional array. It simply calls \ref separableMultiDistSquared()
    and takes the pixel-wise square root of the result. See \ref separableMultiDistSquared()
    for more documentation.
    
    <b> Usage:</b>

    <b>\#include</b> \<<a href="multi__distance_8hxx-source.html">vigra/multi_distance.hxx</a>\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...

    // Calculate Euclidean distance squared for all background pixels 
    separableMultiDistance(srcMultiArrayRange(source), destMultiArray(dest), true);
    \endcode

    \see vigra::distanceTransform(), vigra::separableMultiDistSquared()
*/
doxygen_overloaded_function(template <...> void separableMultiDistance)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void separableMultiDistance( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, bool background,
                             Array const & pixelPitch)
{
    separableMultiDistSquared( s, shape, src, d, dest, background, pixelPitch);
    
    // Finally, calculate the square root of the distances
    transformMultiArray( d, shape, dest, d, dest, (double(*)(double))&std::sqrt );
}

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
          class DestIterator, class DestAccessor, class Array>
inline void separableMultiDistance( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest, bool background,
                                    Array const & pixelPitch)
{
    separableMultiDistance( source.first, source.second, source.third,
                            dest.first, dest.second, background, pixelPitch );
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
