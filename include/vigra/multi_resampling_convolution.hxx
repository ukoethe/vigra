/************************************************************************/
/*                                                                      */
/*         Copyright 1998-2016 by Sven Peter and Ullrich Koethe         */
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

#ifndef VIGRA_MULTIRESAMPLINGCONVOLUTION_HXX
#define VIGRA_MULTIRESAMPLINGCONVOLUTION_HXX

#include "utilities.hxx"
#include "array_vector.hxx"
#include "copyimage.hxx"
#include "multi_array.hxx"

namespace vigra {

namespace detail {
template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor, class Kernel>
void
internalResamplingSeparableConvolveMultiArrayOneDimension(
                      SrcIterator si, Shape const & sshape, SrcAccessor src,
                      DestIterator di, Shape const & dshape, DestAccessor dest, 
                      Kernel const & kernel, unsigned int d)
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
    ArrayVector<Kernel > kernels(1);
    kernels[0] = kernel;

    // temporary array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( ssize );
    typename ArrayVector<TmpType>::iterator t = tmp.begin(), tend = tmp.end();
    typename AccessorTraits<TmpType>::default_accessor ta;
    
    for( ; snav.hasMore(); snav++, dnav++ )
    {
        // first copy source to temp for maximum cache efficiency
        copyLine( snav.begin(), snav.end(), src, t, ta);
        resamplingConvolveLine(t, tend, ta,
                               dnav.begin(), dnav.begin() + dsize, dest,
                               kernels, mapCoordinate);
    }
}

};


template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class KernelIterator>
void resamplingSeparableConvolveMultiArray(SrcIterator si, Shape const & shape, SrcAccessor src,
                             DestIterator di, Shape const & dshape, DestAccessor dest,
                             KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    typedef typename AccessorTraits<TmpType>::default_accessor TmpAccessor;

    if (N == 1) {
        detail::internalResamplingSeparableConvolveMultiArrayOneDimension(si, shape, src, 
                      di, dshape, dest, kit[0], 0);
        return;
    }

    unsigned int d = 0;
    Shape tmpShape(shape);
    tmpShape[d] = dshape[d];
    MultiArray<N, TmpType> tmp(tmpShape);
    TmpAccessor ta;

    detail::internalResamplingSeparableConvolveMultiArrayOneDimension(si, shape, src, 
                         tmp.traverser_begin(), tmpShape, ta, kit[d], d);

    for (d = 1; d < N-1; ++d) {
        tmpShape[d] = dshape[d];
        MultiArray<N, TmpType> dtmp(tmpShape);
        detail::internalResamplingSeparableConvolveMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                             dtmp.traverser_begin(), tmpShape, ta, kit[d], d);
    }
    detail::internalResamplingSeparableConvolveMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                         di, dshape, dest, kit[d], d);
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void resamplingSeparableConvolveMultiArray(SrcIterator s, Shape const & shape, SrcAccessor src,
                             DestIterator d, Shape const & dshape, DestAccessor dest,
                             Kernel1D<T> const & kernel)
{
    ArrayVector<Kernel1D<T> > kernels(shape.size(), kernel);

    resamplingSeparableConvolveMultiArray( s, shape, src, d, dshape, dest, kernels.begin());
}

} // namespace vigra

#endif