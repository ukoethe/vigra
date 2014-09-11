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

template <class VEC, class ARRAY>
inline double 
partialSquaredMagnitude(const VEC& vec, MultiArrayIndex dim, ARRAY const & pixel_pitch)
{
    //computes the squared magnitude of vec
    //considering only the first dim dimensions
    double sqMag = 0.0; 
    for(MultiArrayIndex i=0; i<dim; ++i)
    {
        sqMag += sq(pixel_pitch[i]*vec[i]);
    }
    return sqMag;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Array>
void 
vectorialDistParabola(MultiArrayIndex dimension,
                      SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestAccessor da, 
                      Array const & pixel_pitch )
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename SrcType::value_type SrcPixelType;
    typedef VectorialDistParabolaStackEntry<DestType, SrcPixelType> Influence;
    
    double sigma = pixel_pitch[dimension],
           sigma2 = sq(sigma);
    double w = iend - is; //width of the scanline
   
    std::vector<Influence> _stack; //stack of influence parabolas
    double psm = partialSquaredMagnitude(sa(is), dimension+1, pixel_pitch);
    _stack.push_back(Influence(sa(is), psm, 0.0, 0.0, w));
    ++is;
    psm = partialSquaredMagnitude(sa(is), dimension+1, pixel_pitch);
    double current = 1.0;
    while(current < w)
    {
        Influence & s = _stack.back();
        double diff = current - s.center;
        double intersection = current + ( psm - s.prevVal - sq(sigma*diff) ) / (2.0*sigma2 * diff);
        
        if( intersection < s.left) // previous point has no influence
        {
            _stack.pop_back();
            if(_stack.empty())
            {
                _stack.push_back(Influence(sa(is), psm, 0.0, current, w));
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
        }
        else {
        }
        ++is;
        psm = partialSquaredMagnitude(sa(is), dimension+1, pixel_pitch);
        ++current;
    }
   
    // Now we have the stack indicating which rows are influenced by (and therefore
    // closest to) which row. We can go through the stack and calculate the
    // distance squared for each element of the column.
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(current = 0.0; current < w; ++current, ++id)
    {
        while( current >= it->right) {
            ++it; 
        }
        
        da.set(it->prevVector, id);
        da.setComponent((it->center - current) , id, dimension);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class Array>
inline void vectorialDistParabola(MultiArrayIndex dimension,
                         triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         Array const & pixel_pitch)
{
    vectorialDistParabola(dimension,
                 src.first, src.second, src.third,
                 dest.first, dest.second, pixel_pitch);
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

            detail::vectorialDistParabola(0 /*dimension*/, 
                          srcIterRange(tmp.begin(), tmp.end(),
                          typename AccessorTraits<TmpType>::default_const_accessor()),
                          destIter( dnav.begin(), dest ),
                          sigmas);
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

             detail::vectorialDistParabola(d, 
                           srcIterRange(tmp.begin(), tmp.end(),
                           typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ), sigmas );
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
/*               separableVectorDistance                */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void 
separableVectorDistance(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                        DestIterator d, DestAccessor dest, 
                        bool background,
                        Array const & pixelPitch)
{
    int N = shape.size();

    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type DestPixelType;
    //typedef typename NumericTraits<DestType>::RealPromote Real;

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
inline void 
separableVectorDistance(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                        pair<DestIterator, DestAccessor> const & dest, 
                        bool background,
                        Array const & pixelPitch)
{
    separableVectorDistance( source.first, source.second, source.third,
                               dest.first, dest.second, background, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
separableVectorDistance(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                        DestIterator d, DestAccessor dest, 
                        bool background)
{
    ArrayVector<double> pixelPitch(shape.size(), 1.0);
    separableVectorDistance( s, shape, src, d, dest, background, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
separableVectorDistance(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                        pair<DestIterator, DestAccessor> const & dest, 
                        bool background)
{
    separableVectorDistance( source.first, source.second, source.third,
                               dest.first, dest.second, background );
}


template <unsigned int N, class T1, class S1,
          class T2, class S2, class Array>
inline void 
separableVectorDistance(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest, 
                        bool background,
                        Array const & pixelPitch)
{
    vigra_precondition(source.shape() == dest.shape(),
        "separableVectorDistance(): shape mismatch between input and output.");
    separableVectorDistance(srcMultiArrayRange(source),
                            destMultiArray(dest), background, pixelPitch);
}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void 
separableVectorDistance(MultiArrayView<N, T1, S1> const & source,
                        MultiArrayView<N, T2, S2> dest, 
                        bool background)
{
    vigra_precondition(source.shape() == dest.shape(),
        "separableVectorDistance(): shape mismatch between input and output.");
    separableVectorDistance(srcMultiArrayRange(source),
                            destMultiArray(dest), background );
}

} //-- namespace vigra

#endif        //-- VIGRA_VECTORIAL_DISTANCE_HXX
