/************************************************************************/
/*                                                                      */
/*    Copyright 2003-2013 by Thorben Kroeger, Kasim Terzic,             */
/*        Christian-Dennis Rahn and Ullrich Koethe                      */
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

#ifndef VIGRA_VECTORIAL_DISTANCE_HXX
#define VIGRA_VECTORIAL_DISTANCE_HXX

#include <vector>
#include <set>
#include <functional>
#include "array_vector.hxx"
#include "multi_array.hxx"
#include "accessor.hxx"
#include "numerictraits.hxx"
#include "navigator.hxx"
#include "metaprogramming.hxx"
#include "multi_pointoperators.hxx"
#include "functorexpression.hxx"

#undef VECTORIAL_DIST_DEBUG 

namespace vigra
{

namespace detail
{

template <class Vector, class Value>
struct VectorialDistParabolaStackEntry
{
    double left, center, right;
    Value prevVal;
    Vector prevVector;
    
    VectorialDistParabolaStackEntry(const Vector& vec, Value prev, double l, double c, double r)
    : left(l), center(c), right(r), prevVal(prev), prevVector(vec)
    {}
};

#ifdef VECTORIAL_DIST_DEBUG
template <class Vector, class Value>
std::ostream& operator<<(std::ostream&o, const VectorialDistParabolaStackEntry<Vector, Value>& e) {
    o << "l=" << e.left << ", c=" << e.center << ", r=" << e.right << ", pV=" << e.prevVal << ", pVec=" << e.prevVector;
    return o;
}
#endif

/********************************************************/
/*                                                      */
/*               vectorialDistParabola                  */
/*                                                      */
/********************************************************/

template <class VEC>
inline double partialSquaredMagnitude(const VEC& vec, MultiArrayIndex dim)
{
    //computes the squared magnitude of vec
    //considering only the first dim dimensions
    double sqMag = 0.0; 
    for(MultiArrayIndex i=0; i<dim; ++i)
    {
        sqMag += vec[i]*vec[i];
    }
    return sqMag;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor >
void vectorialDistParabola(MultiArrayIndex dimension, double dmax,
                  SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da, double sigma )
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename SrcType::value_type SrcPixelType;
    typedef VectorialDistParabolaStackEntry<DestType, SrcPixelType> Influence;
    
    double sigma2 = sq(sigma);
    double w = iend - is; //width of the scanline
   
    #ifdef VECTORIAL_DIST_DEBUG
    {
        std::cout << "* dim=" << dimension << ", w=" << w << std::endl;
        
        std::cout << "  source line: ";
        SrcIterator is_copy(is);
        for(; is != iend; ++is) {
            std::cout << sa(is) << " ";
        }
        std::cout << std::endl;
        is = is_copy; 
    }
    #endif
    
    std::vector<Influence> _stack; //stack of influence parabolas
   
    
   
    double psm = partialSquaredMagnitude(sa(is), dimension+1);
    _stack.push_back(Influence(sa(is), psm, 0.0, 0.0, w));
    #ifdef VECTORIAL_DIST_DEBUG
    std::cout << "  building stack:" << std::endl;
    std::cout << "    +stack @[0] " << _stack.back() << std::endl;
    #endif
    
    ++is;
    psm = partialSquaredMagnitude(sa(is), dimension+1);
    
    double current = 1.0;
    while(current < w)
    {
        Influence & s = _stack.back();
        double diff = current - s.center;
        //Bailey 2004, eq. (14)
        //Compute the intersection of the two parabolas
        double intersection = current + ( psm - s.prevVal - sq(sigma*diff) ) / (2.0*sigma2 * diff);
        
        if( intersection < s.left) // previous point has no influence
        {
            #ifdef VECTORIAL_DIST_DEBUG
            std::cout << "    -stack @[" << current << "] " << _stack.back() << std::endl;
            #endif
            _stack.pop_back();
            if(_stack.empty())
            {
                _stack.push_back(Influence(sa(is), psm, 0.0, current, w));
                #ifdef VECTORIAL_DIST_DEBUG
                std::cout << "    +stack @[" << current << "] " << _stack.back() << std::endl;
                #endif
            }
            else
            {
                continue; // try new top of stack without advancing current
            }
        }
        else if(intersection < s.right)
        {
            s.right = intersection;
            _stack.push_back(Influence(sa(is), psm, intersection, current, w));
            #ifdef VECTORIAL_DIST_DEBUG
            std::cout << "    +stack @[" << current << "] " << _stack.back() << std::endl;
            #endif
        }
        else {
        }
        ++is;
        psm = partialSquaredMagnitude(sa(is), dimension+1);
        ++current;
    }
   
    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.
    #ifdef VECTORIAL_DIST_DEBUG
    std::cout << "  stack:" << std::endl;
    for(int i=0; i<_stack.size(); ++i) {
        std::cout << "    " << _stack[i] << std::endl;
    }
    #endif
    
    #ifdef VECTORIAL_DIST_DEBUG
    std::cout << "  writing line: " << std::endl;
    #endif
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(current = 0.0; current < w; ++current, ++id)
    {
        while( current >= it->right) {
            #ifdef VECTORIAL_DIST_DEBUG
            std::cout << "    -stack: " << *it << " (current=" << current << " >= it->right=" << it->right << ")" << std::endl;
            #endif
            ++it; 
        }
        
        da.set(it->prevVector, id);
        
        if( it->prevVector[dimension] == 0 ) {
        //if(it->prevVector != DestType(dmax)) {
            da.setComponent(sigma * (it->center - current) , id, dimension);
            #ifdef VECTORIAL_DIST_DEBUG
            std::cout << "  [" << current << "] = " << da(id) << " (prevVector=" << it->prevVector << ", set @[" << dimension << "]=" << sigma*(it->center - current) << std::endl;
            #endif
        }
        else {
            #ifdef VECTORIAL_DIST_DEBUG
            std::cout << "  [" << current << "] = " << da(id) << " (old value)" << std::endl;
            #endif
        }
    }
    #ifdef VECTORIAL_DIST_DEBUG
    std::cout << std::endl;
    #endif
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void vectorialDistParabola(MultiArrayIndex dimension, double dmax,
                         triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest, double sigma)
{
    vectorialDistParabola(dimension, dmax,
                 src.first, src.second, src.third,
                 dest.first, dest.second, sigma);
}

/********************************************************/
/*                                                      */
/*      internalSeparableMultiVectorialDistTmp          */
/*                                                      */
/********************************************************/

template<int T, class U>
struct error_wrong_pixel_type : vigra::staticAssert::AssertBool<
        (T == U::static_size)
> {};

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void internalSeparableMultiVectorialDistTmp(
                      double dmax,
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, Array const & sigmas, bool invert)
{
    // Sigma is the spread of the parabolas. It determines the structuring element size
    // for ND morphology. When calculating the distance transforms, sigma is usually set to 1,
    // unless one wants to account for anisotropic pixel pitch
    enum { N =  SrcShape::static_size};
    
    vigra_precondition(sigmas.size() == SrcShape::static_size, "sigmas has wrong length");
    VIGRA_STATIC_ASSERT((error_wrong_pixel_type<SrcShape::static_size, typename DestAccessor::value_type>));

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

            detail::vectorialDistParabola(0 /*dimension*/, dmax,
                          srcIterRange(tmp.begin(), tmp.end(),
                          typename AccessorTraits<TmpType>::default_const_accessor()),
                          destIter( dnav.begin(), dest ),
                           
                          sigmas[0] );
    }
    
    #if 1
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

             detail::vectorialDistParabola(d, dmax,
                           srcIterRange(tmp.begin(), tmp.end(),
                           typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ), sigmas[d] );
        }
    }
    #endif
    
    if(invert) transformMultiArray( di, shape, dest, di, dest, -Arg1());

}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void internalSeparableMultiVectorialDistTmp(
    double dmax,
    SrcIterator si, SrcShape const & shape, SrcAccessor src,
    DestIterator di, DestAccessor dest, Array const & sigmas)
{
    internalSeparableMultiVectorialDistTmp( dmax, si, shape, src, di, dest, sigmas, false );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void internalSeparableMultiVectorialDistTmp(
    double dmax,
    SrcIterator si, SrcShape const & shape, SrcAccessor src,
    DestIterator di, DestAccessor dest)
{
    ArrayVector<double> sigmas(shape.size(), 1.0);
    internalSeparableMultiVectorialDistTmp( dmax, si, shape, src, di, dest, sigmas, false );
}

} // namespace detail

/********************************************************/
/*                                                      */
/*           separableMultiVectorialDist                */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void separableMultiVectorialDist( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, bool background,
                                Array const & pixelPitch)
{
    int N = shape.size();

    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type DestPixelType;
    typedef typename NumericTraits<DestType>::RealPromote Real;

    SrcType zero = NumericTraits<SrcType>::zero();

    using namespace vigra::functor;
    
    double dmax = 0.0;
    for( int k=0; k<N; ++k)
    {
        dmax += pixelPitch[k]*shape[k];
    }
    
    // Threshold the values so all objects have infinity value in the beginning
    DestType maxDist = DestType(std::ceil(dmax)), rzero = DestType(static_cast<DestPixelType>(0));
    if(background == true)
        transformMultiArray( s, shape, src, d, dest,
                                ifThenElse( Arg1() == Param(zero), Param(maxDist), Param(rzero) ));
    else
        transformMultiArray( s, shape, src, d, dest, 
                                ifThenElse( Arg1() != Param(zero), Param(maxDist), Param(rzero) ));
    
    detail::internalSeparableMultiVectorialDistTmp( dmax, d, shape, dest, d, dest, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void separableMultiVectorialDist( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest, bool background,
                                       Array const & pixelPitch)
{
    separableMultiVectorialDist( source.first, source.second, source.third,
                               dest.first, dest.second, background, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline 
void separableMultiVectorialDist( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest, bool background)
{
    ArrayVector<double> pixelPitch(shape.size(), 1.0);
    separableMultiVectorialDist( s, shape, src, d, dest, background, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void separableMultiVectorialDist( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest, bool background)
{
    separableMultiVectorialDist( source.first, source.second, source.third,
                               dest.first, dest.second, background );
}

} //-- namespace vigra

#endif        //-- VIGRA_VECTORIAL_DISTANCE_HXX
