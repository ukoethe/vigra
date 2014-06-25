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

#ifndef VIGRA_VECTORIAL_BOUNDARY_DISTANCE_HXX
#define VIGRA_VECTORIAL_BOUNDARY_DISTANCE_HXX


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

#undef VECTORIAL_BOUNDARY_DIST_DEBUG

namespace vigra
{

namespace detail
{

/********************************************************/
/*                                                      */
/*               boundaryVectorialDistParabola          */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor >
void boundaryVectorialDistParabola(MultiArrayIndex dimension, double dmax,
                  SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da, double sigma )
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename SrcType::value_type SrcPixelType;
    typedef VectorialDistParabolaStackEntry<DestType, SrcPixelType> Influence;

    double sigma2 = sq(sigma);
    double w = iend - is; //width of the scanline
    if(w <= 0)
        return;

    std::vector<Influence> _stack; //stack of influence parabolas
    double psm;
    SrcType value;
    SrcType label_check = sa(is);
    bool label_check2;
    if (sa(is) != sa(++is)) {
        psm = partialSquaredMagnitude(SrcType(0.0), dimension+1);
        _stack.push_back(Influence(SrcType(0.0), psm, 0.0, 0.0, w));
        label_check2 = true;
    }
    else  {
        psm = partialSquaredMagnitude(SrcType(dmax), dimension+1);
        _stack.push_back(Influence(SrcType(dmax), psm, 0.0, 0.0, w));
        label_check2 = false;
    }


    bool label_check3 = false;
    double begin = 0.0, end;

    double current = 1.0;
    while(current < w)
    {
        if (label_check3 == true) (label_check3 = false);
        else {
            if (label_check2 == true){
                current -= 1;
                end = current;
                // Now we have the stack indicating which rows are influenced by (and therefore
                // closest to) which row. We can go through the stack and calculate the
                // distance squared for each element of the column.
                typename std::vector<Influence>::iterator it = _stack.begin();
                for(float i = begin ; i < end; ++i, ++id)
                {
                    while( i >= it->right)
                        ++it;
                    da.set(it->prevVector, id);
                    if( it->prevVector[dimension] == 0 ) {
                    //if(it->prevVector != DestType(dmax)) {
                    da.setComponent(sigma * (it->center - current) , id, dimension);
                    }
                }
                while (_stack.empty() == false) (_stack.pop_back());
                begin = current;
                psm = partialSquaredMagnitude(SrcType(0.0), dimension+1);
                _stack.push_back(Influence(SrcType(0.0), psm, begin-1, begin-1, w));
                --is;
                label_check2 = false;
                value = dmax;
            }
            else if (label_check == sa(is)){
                value = dmax;
            }
            else if (label_check != sa(is)){
                label_check2 = true;
                value = 0.0;
            }
            label_check = sa(is);
            psm = partialSquaredMagnitude(value, dimension+1);
        }
        Influence & s = _stack.back();
        double diff = current - s.center;
        //Bailey 2004, eq. (14)
        //Compute the intersection of the two parabolas
        double intersection = current + ( psm - s.prevVal - sq(sigma*diff) ) / (2.0*sigma2 * diff);

        if( intersection < s.left) // previous point has no influence
        {
            _stack.pop_back();
            if(_stack.empty())
            {
                _stack.push_back(Influence(value, psm, 0.0, current, w));

            }
            else
            {
                label_check3 = true;
                continue; // try new top of stack without advancing current
            }
        }
        else if(intersection < s.right)
        {
            s.right = intersection;
            _stack.push_back(Influence(value, psm, intersection, current, w));
        }
        ++is;
        ++current;
    }
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(float i = 0.0; i < w; ++i, ++id)
    {
        while( i >= it->right) {
            ++it;
        }
        da.set(it->prevVector, id);
        if( it->prevVector[dimension] == 0 ) {
        //if(it->prevVector != DestType(dmax)) {
        da.setComponent(sigma * (it->center - current) , id, dimension);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void boundaryVectorialDistParabola(MultiArrayIndex dimension, double dmax,
                         triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest, double sigma)
{
    boundaryVectorialDistParabola(dimension, dmax,
                 src.first, src.second, src.third,
                 dest.first, dest.second, sigma);
}


/********************************************************/
/*                                                      */
/*      internalBoundaryMultiVectorialDistTmp           */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void internalBoundaryMultiVectorialDistTmp(
                      double dmax,
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, Array const & sigmas)
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

            copyLine( snav.begin(), snav.end(), src, tmp.begin(),
                          typename AccessorTraits<TmpType>::default_accessor() );

            detail::boundaryVectorialDistParabola(0 /*dimension*/, dmax,
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

}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void internalBoundaryMultiVectorialDistTmp(
    double dmax,
    SrcIterator si, SrcShape const & shape, SrcAccessor src,
    DestIterator di, DestAccessor dest)
{
    ArrayVector<double> sigmas(shape.size(), 1.0);
    internalBoundaryMultiVectorialDistTmp( dmax, si, shape, src, di, dest, sigmas);
}

} // namespace detail

/********************************************************/
/*                                                      */
/*           boundaryMultiVectorialDist                 */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
void boundaryMultiVectorialDist( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest,
                                Array const & pixelPitch)
{
    int N = shape.size();

    typedef typename DestAccessor::value_type DestType;
    typedef typename NumericTraits<DestType>::RealPromote Real;

    using namespace vigra::functor;

    double dmax = 0.0;
    bool pixelPitchIsReal = false;
    for( int k=0; k<N; ++k)
    {
        if(int(pixelPitch[k]) != pixelPitch[k])
            pixelPitchIsReal = true;
        dmax += sq(pixelPitch[k]*shape[k]);
    }
//    if(dmax > NumericTraits<DestType>::toRealPromote(NumericTraits<DestType>::max())
//       || pixelPitchIsReal) // need a temporary array to avoid overflows
//    {
//        MultiArray<SrcShape::static_size, Real> tmpArray(shape);
//        transformMultiArray( s, shape, src,
//                             tmpArray.traverser_begin(), typename AccessorTraits<Real>::default_accessor(),
//                             Arg1());
//        detail::internalBoundaryMultiVectorialDistTmp( dmax, tmpArray.traverser_begin(),
//                shape, typename AccessorTraits<Real>::default_accessor(),
//                tmpArray.traverser_begin(),
//                typename AccessorTraits<Real>::default_accessor(), pixelPitch);

//        copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
//    }
//    else        // work directly on the destination array
//    {
        detail::internalBoundaryMultiVectorialDistTmp( dmax, s, shape, src, d, dest, pixelPitch);
//    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void boundaryMultiVectorialDist( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest,
                                       Array const & pixelPitch)
{
    boundaryMultiVectorialDist( source.first, source.second, source.third,
                               dest.first, dest.second, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void boundaryMultiVectorialDist( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                                DestIterator d, DestAccessor dest)
{
    ArrayVector<double> pixelPitch(shape.size(), 1.0);
    boundaryMultiVectorialDist( s, shape, src, d, dest, pixelPitch );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void boundaryMultiVectorialDist( triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest)
{
    boundaryMultiVectorialDist( source.first, source.second, source.third,
                               dest.first, dest.second);
}

} //-- namespace vigra

#endif        //-- VIGRA_VECTORIAL_DISTANCE_HXX
