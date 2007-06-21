/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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


#ifndef VIGRA_MULTI_RESIZE_HXX
#define VIGRA_MULTI_RESIZE_HXX

#include <vector>
#include "resizeimage.hxx"

namespace vigra {

namespace detail {

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Kernel>
void
internalResizeMultiArrayOneDimension(
                      SrcIterator si, Shape const & sshape, SrcAccessor src,
                      DestIterator di, Shape const & dshape, DestAccessor dest, 
                      Kernel const & spline, unsigned int d)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    SNavigator snav( si, sshape, d );
    DNavigator dnav( di, dshape, d );
    
    int ssize = sshape[d];
    int dsize = dshape[d];

    vigra_precondition(ssize > 1,
                 "resizeMultiArraySplineInterpolation(): "
                 "Source array too small.\n");

    Rational<int> ratio(dsize - 1, ssize - 1);
    Rational<int> offset(0);
    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(ratio, offset);
    int period = lcm(ratio.numerator(), ratio.denominator());
    
    ArrayVector<double> const & prefilterCoeffs = spline.prefilterCoefficients();
    ArrayVector<Kernel1D<double> > kernels(period);
    createResamplingKernels(spline, mapCoordinate, kernels);

    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( ssize );
    ArrayVector<TmpType>::iterator t = tmp.begin(), tend = tmp.end();
    typename AccessorTraits<TmpType>::default_accessor ta;
    
    for( ; snav.hasMore(); snav++, dnav++ )
    {
        // first copy source to temp for maximum cache efficiency
        copyLine( snav.begin(), snav.end(), src, t, ta);

        for(unsigned int b = 0; b < prefilterCoeffs.size(); ++b)
        {
            recursiveFilterLine(t, tend, ta, t, ta,
                                prefilterCoeffs[b], BORDER_TREATMENT_REFLECT);
        }
        resamplingConvolveLine(t, tend, ta,
                               dnav.begin(), dnav.begin() + dsize, dest,
                               kernels, mapCoordinate);
    }
}

} // namespace detail

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Kernel>
void
resizeMultiArraySplineInterpolation(
                      SrcIterator si, Shape const & sshape, SrcAccessor src,
                      DestIterator di, Shape const & dshape, DestAccessor dest, 
                      Kernel const & spline)
{
    enum { N = 1 + SrcIterator::level };
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    typedef MultiArray<N, TmpType> TmpArray;
    typedef typename AccessorTraits<TmpType>::default_accessor TmpAccessor;
        
    if(N==1)
    {
        detail::internalResizeMultiArrayOneDimension(si, sshape, src, 
                      di, dshape, dest, spline, 0);
    }
    else
    {
        unsigned int d = 0;
        Shape tmpShape(sshape);
        tmpShape[d] = dshape[d];
        MultiArray<N, TmpType> tmp(tmpShape);
        TmpAccessor ta;
        
        detail::internalResizeMultiArrayOneDimension(si, sshape, src, 
                             tmp.traverser_begin(), tmpShape, ta, spline, d);
        d = 1;
        for(; d<N-1; ++d)
        {
            tmpShape[d] = dshape[d];
            MultiArray<N, TmpType> dtmp(tmpShape);
            
            detail::internalResizeMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                                  dtmp.traverser_begin(), tmpShape, ta, spline, d);
            dtmp.swap(tmp);
        }
        detail::internalResizeMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                                        di, dshape, dest, spline, d);
    }
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor, 
          class Kernel>
inline void
resizeMultiArraySplineInterpolation(triple<SrcIterator, Shape, SrcAccessor> src,
                      triple<DestIterator, Shape, DestAccessor> dest,
                      Kernel const & spline)
{
    resizeMultiArraySplineInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third, spline);
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeMultiArraySplineInterpolation(
                      SrcIterator si, Shape const & sshape, SrcAccessor src,
                      DestIterator di, Shape const & dshape, DestAccessor dest)
{
    resizeMultiArraySplineInterpolation(si, sshape, src, di, dshape, dest, BSpline<3, double>());
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
resizeMultiArraySplineInterpolation(triple<SrcIterator, Shape, SrcAccessor> src,
                      triple<DestIterator, Shape, DestAccessor> dest)
{
    resizeMultiArraySplineInterpolation(src.first, src.second, src.third,
                                   dest.first, dest.second, dest.third);
}

//@}

} // namespace vigra

#endif // VIGRA_MULTI_RESIZE_HXX
