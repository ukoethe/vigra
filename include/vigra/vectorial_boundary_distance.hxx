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

#include "multi_gridgraph.hxx" 	//boundaryMultiDistance
#include "union_find.hxx"		//boundaryMultiDistance

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
          class BufIterator, class BufAccessor,
          class DestIterator, class DestAccessor >
void boundaryVectorialDistParabola(MultiArrayIndex dimension, double dmax,
                  SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  BufIterator bis, BufIterator biend, BufAccessor ba,
                  DestIterator id, DestAccessor da, double sigma )
{
    typedef typename SrcAccessor::value_type SrcType;
    typedef typename BufAccessor::value_type BufType;
    typedef typename DestAccessor::value_type DestType;
    typedef typename BufType::value_type BufPixelType;
    typedef VectorialDistParabolaStackEntry<DestType, BufPixelType> Influence;

    double sigma2 = sq(sigma);
    double w = iend - is; //width of the scanline
    if(w <= 0)
        return;

    std::vector<Influence> _stack; //stack of influence parabolas
    bool label_check2 = false;
    SrcType label_check = sa(is);
    BufType value;
    double psm= partialSquaredMagnitude(ba(bis), dimension+1);
    _stack.push_back(Influence(ba(bis), psm, 0.0, 0.0, w));
    ++is;
    ++bis;
    bool label_check3 = false;
    double begin = 0.0, current = 1.0;
    while(current < w)
    {
        if (label_check3 == true) (label_check3 = false);
        else {
            if (label_check2 == true){
                // Now we have the stack indicating which rows are influenced by (and therefore
                // closest to) which row. We can go through the stack and calculate the
                // distance squared for each element of the column.
                typename std::vector<Influence>::iterator it = _stack.begin();
                for(float i = begin ; i < current; ++i, ++id)
                {
                    while( i >= it->right)
                        ++it;
                    da.set(it->prevVector, id);
                    if( it->prevVector[dimension] == 0 ) {
                    //if(it->prevVector != DestType(dmax)) {
                    da.setComponent(sigma * (it->center - i) , id, dimension);
                    }
                }
                while (_stack.empty() == false) (_stack.pop_back());
                begin = current;
                psm = partialSquaredMagnitude(BufType(0.0), dimension+1);
                _stack.push_back(Influence(BufType(0.0), psm, begin-1, begin-1, w));
                label_check2 = false;
                value = ba(bis);
            }
            else if (label_check == sa(is)){
                value = ba(bis);
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
                _stack.push_back(Influence(value, psm, begin, current, w));

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
        if (label_check2 == false)
        {
            ++is;
            ++current;
            ++bis;
        }
    }
    typename std::vector<Influence>::iterator it = _stack.begin();
    for(float current = begin; current < w; ++current, ++id)
    {
        while( current >= it->right)
            ++it;
        da.set(it->prevVector, id);
        if( it->prevVector[dimension] == 0 ) {
        //if(it->prevVector != DestType(dmax)) {
        da.setComponent(sigma * (it->center - current) , id, dimension);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class BufIterator, class BufAccessor,
          class DestIterator, class DestAccessor>
inline void boundaryVectorialDistParabola(MultiArrayIndex dimension, double dmax,
                         triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         triple<BufIterator, BufIterator, BufAccessor> buffer,
                         pair<DestIterator, DestAccessor> dest, double sigma)
{
    boundaryVectorialDistParabola(dimension, dmax,
                 src.first, src.second, src.third, buffer.first, buffer.second,
                 buffer.third, dest.first, dest.second, sigma);
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
    tmp.init(dmax);
    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    // only operate on first dimension here
    SNavigator snav( si, shape, 0 );
    DNavigator dnav( di, shape, 0 );

    using namespace vigra::functor;

    for( ; snav.hasMore(); snav++, dnav++ )
    {
            // first copy source to temp for maximum cache efficiency

            detail::boundaryVectorialDistParabola(0 /*dimension*/, dmax,
                          srcIterRange(snav.begin(), snav.end(), src),
                          srcIterRange(tmp.begin(), tmp.end(),
                          typename AccessorTraits<TmpType>::default_accessor()),
                          destIter( dnav.begin(), dest ), sigmas[0] );
    }

    #if 1
    // operate on further dimensions
    for( int d = 1; d < N; ++d )
    {
        DNavigator dnav( di, shape, d );
        SNavigator snav( si, shape, d );
        tmp.resize( shape[d] );


        for( ; dnav.hasMore(); dnav++, snav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( dnav.begin(), dnav.end(), dest,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             detail::boundaryVectorialDistParabola(d, dmax,
                           srcIterRange(snav.begin(), snav.end(), src),
                           srcIterRange(tmp.begin(), tmp.end(),
                           typename AccessorTraits<TmpType>::default_accessor()),
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

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Array>
inline void internalBoundaryMultiVectorialDistTmp(
                                       double dmax,
                                       triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest,
                                       Array const & pixelPitch)
{
    internalBoundaryMultiVectorialDistTmp( dmax, source.first, source.second, source.third,
                               dest.first, dest.second, pixelPitch );
}

} // namespace detail

/********************************************************/
/*                                                      */
/*           boundaryMultiVectorialDist                 */
/*                                                      */
/********************************************************/

template <unsigned int N, class T1, class S1,
          class T2, class S2, class Array>
void boundaryMultiVectorialDist( MultiArrayView<N, T1, S1> const & source,
                                 MultiArrayView<N, T2, S2> dest,
                                 Array const & pixelPitch)
{
    using namespace vigra::functor;

    double dmax = 0.0;
    for( int k=0; k<N; ++k)
    {
        dmax += sq(pixelPitch[k]*source.shape(k));
    }
        detail::internalBoundaryMultiVectorialDistTmp( dmax, srcMultiArrayRange(source), destMultiArray(dest), pixelPitch);

    typedef typename GridGraph<N, undirected_tag>::NodeIt        graph_scanner;
    typedef typename GridGraph<N, undirected_tag>::OutArcIt  neighbor_iterator;
    GridGraph<N, undirected_tag> g(source.shape());
    double min_mag;
    typename MultiArrayView<N, T2, S2>::value_type min_pos, min_vec, vec_to_pix;
    for (graph_scanner node(g); node != lemon_graph::INVALID; ++node)
    {
        vec_to_pix = dest[*node];
        min_mag = detail::partialSquaredMagnitude(vec_to_pix, N);
        min_pos = *node;

        //go to adjacent neighbour with different label of target pixel with smallest distance to origin pixel
        for (neighbor_iterator arc(g, *node+vec_to_pix); arc != lemon_graph::INVALID; ++arc)
        {
                            if(source[*node+vec_to_pix] != source[g.target(*arc)])
                            {
                                if (min_mag > detail::partialSquaredMagnitude(g.target(*arc)-*node,N))
                                {
                                    min_mag = detail::partialSquaredMagnitude(g.target(*arc)-*node,N);
                                    min_pos = g.target(*arc);
                                    //std::cout << vec_to_pix - g.target(*arc)-*node << std::endl;
                                    //std::cout << *node << " found smaller position with " << min_mag << " at " << g.target(*arc) << std::endl;
                                 }
                            }
        }
        //from this pixel look for the vector which points to the nearest interpixel between two label
        min_mag = detail::partialSquaredMagnitude(vec_to_pix, N);
        for (neighbor_iterator arc(g, min_pos); arc != lemon_graph::INVALID; ++arc)
        {
            if(source[min_pos] != source[g.target(*arc)])
            {
                if (min_mag > detail::partialSquaredMagnitude(vec_to_pix - (g.target(*arc) - min_pos)*0.5, N))
                {
                    min_vec = vec_to_pix - (g.target(*arc) - min_pos)*0.5;
                    min_mag = detail::partialSquaredMagnitude(min_vec, N);
                    //std::cout << *node << " found smaller vec with " << min_mag << " at " << g.target(*arc) << std::endl;
                }
            }
        }
    dest[*node] = min_vec;
    }


}

template <unsigned int N, class T1, class S1,
          class T2, class S2>
inline void
boundaryMultiVectorialDist(MultiArrayView<N, T1, S1> const & source,
                       MultiArrayView<N, T2, S2> dest)
{
    ArrayVector<double> pixelPitch(source.shape().size(), 1.0);
    vigra_precondition(source.shape() == dest.shape(),
        "boundaryMultiDistance(): shape mismatch between input and output.");
    boundaryMultiVectorialDist( source, dest, pixelPitch );
}

} //-- namespace vigra

#endif        //-- VIGRA_VECTORIAL_DISTANCE_HXX
