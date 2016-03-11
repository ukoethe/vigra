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
#include "navigator.hxx"

#include <type_traits>

namespace vigra {

namespace detail {
template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
void
internalResamplingSeparableConvolveMultiArrayOneDimension(
                      SrcIterator si, Shape const & sshape, SrcAccessor src,
                      DestIterator di, Shape const & dshape, DestAccessor dest, 
                      KernelIterator const & kernels, unsigned int d,
                      Rational<int> ratio = Rational<int>(0,1))
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
                 "internalResamplingSeparableConvolveMultiArrayOneDimension(): "
                 "Source array too small.\n");

    if (!ratio.numerator())
      ratio = Rational<int>(dsize - 1, ssize - 1);

    Rational<int> offset(0);

    resampling_detail::MapTargetToSourceCoordinate mapCoordinate(ratio, offset);

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
                             KernelIterator const & kit, ArrayVector<Rational<int> > & ratios = ArrayVector<Rational<int> >())
{
    enum { N = 1 + SrcIterator::level };
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    typedef typename AccessorTraits<TmpType>::default_accessor TmpAccessor;

    if (ratios.size() == 0)
      ratios = ArrayVector<Rational<int> >(N, Rational<int>(0, 1));

    vigra_precondition(ratios.size() == N,
                 "resamplingSeparableConvolveMultiArray(): "
                 "Number of ratios does not match number of dimensions.\n");

    if (N == 1) {
        detail::internalResamplingSeparableConvolveMultiArrayOneDimension(si, shape, src, 
                      di, dshape, dest, kit[0], 0, ratios[0]);
        return;
    }

    unsigned int d = 0;
    Shape tmpShape(shape);
    tmpShape[d] = dshape[d];
    MultiArray<N, TmpType> tmp(tmpShape);
    TmpAccessor ta;

    detail::internalResamplingSeparableConvolveMultiArrayOneDimension(si, shape, src, 
                         tmp.traverser_begin(), tmpShape, ta, kit[d], d, ratios[d]);

    for (d = 1; d < N-1; ++d) {
        tmpShape[d] = dshape[d];
        MultiArray<N, TmpType> dtmp(tmpShape);
        detail::internalResamplingSeparableConvolveMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                             dtmp.traverser_begin(), tmpShape, ta, kit[d], d, ratios[d]);
        dtmp.swap(tmp);
    }
    detail::internalResamplingSeparableConvolveMultiArrayOneDimension(tmp.traverser_begin(), tmp.shape(), ta, 
                         di, dshape, dest, kit[d], d, ratios[d]);
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void resamplingSeparableConvolveMultiArray(SrcIterator s, Shape const & shape, SrcAccessor src,
                             DestIterator d, Shape const & dshape, DestAccessor dest,
                             ArrayVector<Kernel1D<T> > const & kernel, ArrayVector<Rational<int> > & ratios = ArrayVector<Rational<int> >())
{
    ArrayVector<ArrayVector<Kernel1D<T> > > kernels(shape.size(), kernel);

    resamplingSeparableConvolveMultiArray( s, shape, src, d, dshape, dest, kernels.begin(), ratios);
}


template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void resamplingSeparableConvolveMultiArray(SrcIterator s, Shape const & shape, SrcAccessor src,
                             DestIterator d, Shape const & dshape, DestAccessor dest,
                             Kernel1D<T> const & kernel, ArrayVector<Rational<int> > & ratios = ArrayVector<Rational<int> >())
{
    ArrayVector<Kernel1D<T> > kernels(1, kernel);

    resamplingSeparableConvolveMultiArray( s, shape, src, d, dshape, dest, kernels, ratios);
}

template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void resamplingSeparableConvolveMultiArray(SrcIterator s, Shape const & shape, SrcAccessor src,
                             DestIterator d, Shape const & dshape, DestAccessor dest,
                             ArrayVector<Kernel1D<T> > const & kernel, Rational<int> ratio)
{
    ArrayVector<ArrayVector<Kernel1D<T> > > kernels(shape.size(), kernel);
    ArrayVector<Rational<int> > ratios(shape.size(), ratio);

    resamplingSeparableConvolveMultiArray( s, shape, src, d, dshape, dest, kernels.begin(), ratios);
}


template <class SrcIterator, class Shape, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void resamplingSeparableConvolveMultiArray(SrcIterator s, Shape const & shape, SrcAccessor src,
                             DestIterator d, Shape const & dshape, DestAccessor dest,
                             Kernel1D<T> const & kernel, Rational<int> ratio)
{
    ArrayVector<Kernel1D<T> > kernels(1, kernel);
    ArrayVector<Rational<int> > ratios(shape.size(), ratio);

    resamplingSeparableConvolveMultiArray( s, shape, src, d, dshape, dest, kernels, ratios);
}

} // namespace vigra

#endif