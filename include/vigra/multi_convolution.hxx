//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Christian-Dennis Rahn                */
/*                        and Ullrich Koethe                            */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_MULTI_CONVOLUTION_H
#define VIGRA_MULTI_CONVOLUTION_H

#include <vigra/separableconvolution.hxx>
#include <vigra/array_vector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/accessor.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/navigator.hxx>
#include <vigra/metaprogramming.hxx>


namespace vigra
{


namespace detail
{

/********************************************************/
/*                                                      */
/*        internalSeparableConvolveMultiArray           */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
void
internalSeparableConvolveMultiArrayTmp(
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;

    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );

        for( ; snav.hasMore(); snav++, dnav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( snav.begin(), snav.end(), src,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             convolveLine( srcIterRange(tmp.begin(), tmp.end(),
                                        typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ),
                           kernel1d( *kit ) );
        }
        ++kit;
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d, ++kit )
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( dnav.begin(), dnav.end(), dest,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             convolveLine( srcIterRange(tmp.begin(), tmp.end(),
                                        typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ),
                           kernel1d( *kit ) );
        }
    }
}


} // namespace detail

/** \addtogroup MultiArrayConvolutionFilters Convolution filters for multi-dimensional arrays.

    These functions realize a separable convolution on an arbitrary dimensional
    array that is specified by iterators (compatible to \ref MultiIteratorPage)
    and shape objects. It can therefore be applied to a wide range of data structures
    (\ref vigra::MultiArrayView, \ref vigra::MultiArray etc.).
*/
//@{

/********************************************************/
/*                                                      */
/*             separableConvolveMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Separated convolution on multi-dimensional arrays.

    This function computes a separated convolution on all dimensions
    of the given multi-dimensional array. Both source and destination
    arrays are represented by iterators, shape objects and accessors.
    The destination array is required to already have the correct size.

    There are two variants of this functions: one takes a single kernel
    of type \ref vigra::Kernel1D which is then applied to all dimensions,
    whereas the other requires an iterator referencing a sequence of
    \ref vigra::Kernel1D objects, one for every dimension of the data.
    Then the first kernel in this sequence is applied to the innermost
    dimension (e.g. the x-dimension of an image), while the last is applied to the
    outermost dimension (e.g. the z-dimension in a 3D image).

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause round-off errors (i.e. if
    <tt>typeid(typename NumericTraits<typename DestAccessor::value_type>::RealPromote)
    != typeid(typename DestAccessor::value_type)</tt>.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        separableConvolveMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    Kernel1D<T> const & kernel);

        // apply each kernel from the sequence `kernels´ in turn
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class KernelIterator>
        void
        separableConvolveMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    KernelIterator kernels);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    Kernel1D<T> const & kernel);

        // apply each kernel from the sequence `kernels´ in turn
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class KernelIterator>
        void
        separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    KernelIterator kernels);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__convolution_8hxx-source.html">vigra/multi_convolution.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    Kernel1D<float> gauss;
    gauss.initGaussian(sigma);

    // perform Gaussian smoothing on all dimensions
    separableConvolveMultiArray(srcMultiArrayRange(source), destMultiArray(dest), gauss);
    \endcode

    \see vigra::Kernel1D, convolveLine()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
void
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, KernelIterator kernels )
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;

    if(!IsSameType<TmpType, typename DestAccessor::value_type>::boolResult)
    {
        // need a temporary array to avoid rounding errors
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src,
             tmpArray.traverser_begin(), typename AccessorTraits<TmpType>::default_accessor(), kernels );
        copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
    }
    else
    {
        // work directly on the destination array
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src, d, dest, kernels );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
inline
void separableConvolveMultiArray(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, KernelIterator kit )
{
    separableConvolveMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, kit );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest,
                             Kernel1D<T> const & kernel )
{
    ArrayVector<Kernel1D<T> > kernels(shape.size(), kernel);

    separableConvolveMultiArray( s, shape, src, d, dest, kernels.begin() );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest,
                            Kernel1D<T> const & kernel )
{
    ArrayVector<Kernel1D<T> > kernels(source.second.size(), kernel);

    separableConvolveMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, kernels.begin() );
}

/********************************************************/
/*                                                      */
/*            convolveMultiArrayOneDimension            */
/*                                                      */
/********************************************************/

/** \brief Convolution along a single dimension of a multi-dimensional arrays.

    This function computes a convolution along one dimension (specified by
    the parameter <tt>dim</tt> of the given multi-dimensional array with the given
    <tt>kernel</tt>. Both source and destination arrays are represented by
    iterators, shape objects and accessors. The destination array is required to
    already have the correct size.

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        convolveMultiArrayOneDimension(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                       DestIterator diter, DestAccessor dest,
                                       unsigned int dim, vigra::Kernel1D<T> const & kernel);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        convolveMultiArrayOneDimension(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest,
                                       unsigned int dim, vigra::Kernel1D<T> const & kernel);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__convolution_8hxx-source.html">vigra/multi_convolution.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    Kernel1D<float> gauss;
    gauss.initGaussian(sigma);

    // perform Gaussian smoothing along dimensions 1 (height)
    convolveMultiArrayOneDimension(srcMultiArrayRange(source), destMultiArray(dest), 1, gauss);
    \endcode

    \see separableConvolveMultiArray()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
void
convolveMultiArrayOneDimension(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                               DestIterator d, DestAccessor dest,
                               unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    enum { N = 1 + SrcIterator::level };
    vigra_precondition( dim < N,
                        "convolveMultiArrayOneDimension(): The dimension number to convolve must be smaller "
                        "than the data dimensionality" );

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    ArrayVector<TmpType> tmp( shape[dim] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    SNavigator snav( s, shape, dim );
    DNavigator dnav( d, shape, dim );

    for( ; snav.hasMore(); snav++, dnav++ )
    {
         // first copy source to temp for maximum cache efficiency
         copyLine( snav.begin(), snav.end(), src,
		   tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

         convolveLine( srcIterRange( tmp.begin(), tmp.end(), typename AccessorTraits<TmpType>::default_const_accessor()),
                       destIter( dnav.begin(), dest ),
                       kernel1d( kernel ) );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
convolveMultiArrayOneDimension(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                               pair<DestIterator, DestAccessor> const & dest,
                               unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    convolveMultiArrayOneDimension( source.first, source.second, source.third,
                                   dest.first, dest.second, dim, kernel );
}

/********************************************************/
/*                                                      */
/*             gaussianSmoothMultiArray                 */
/*                                                      */
/********************************************************/

/** \brief Isotropic Gaussian smoothing of a multi-dimensional arrays.

    This function computes an isotropic convolution of the given multi-dimensional
    array with a Gaussian filter at the given standard deviation <tt>sigma</tt>.
    Both source and destination arrays are represented by
    iterators, shape objects and accessors. The destination array is required to
    already have the correct size. This function may work in-place, which means
    that <tt>siter == diter</tt> is allowed. It is implemented by a call to
    \ref separableConvolveMultiArray() with the appropriate kernel.
    If the data are anisotropic (different pixel size along different dimensions)
    you should call \ref separableConvolveMultiArray() directly with the appropriate
    anisotropic Gaussians.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianSmoothMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                 DestIterator diter, DestAccessor dest,
                                 double sigma);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianSmoothMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                 pair<DestIterator, DestAccessor> const & dest,
                                 double sigma);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__convolution_8hxx-source.html">vigra/multi_convolution.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    // perform isotropic Gaussian smoothing at scale `sigma´
    gaussianSmoothMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma);
    \endcode

    \see separableConvolveMultiArray()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
gaussianSmoothMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                   DestIterator d, DestAccessor dest, double sigma )
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote kernel_type;
    Kernel1D<kernel_type> gauss;
    gauss.initGaussian( sigma );

    separableConvolveMultiArray( s, shape, src, d, dest, gauss);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                  pair<DestIterator, DestAccessor> const & dest,
                  double sigma )
{
    gaussianSmoothMultiArray( source.first, source.second, source.third,
                              dest.first, dest.second, sigma );
}

/********************************************************/
/*                                                      */
/*             gaussianGradientMultiArray               */
/*                                                      */
/********************************************************/

/** \brief Calculate Gaussian gradient of a multi-dimensional arrays.

    This function computes the Gaussian gradient of the given multi-dimensional
    array with a sequence of first-derivative-of-Gaussian filters at the given
    standard deviation <tt>sigma</tt> (differentiation is applied to each dimension
    in turn, starting with the innermost dimension). Both source and destination arrays
    are represented by iterators, shape objects and accessors. The destination array is
    required to have a vector valued pixel type with as many elements as the number of
    dimensions. This function is implemented by calls to
    \ref separableConvolveMultiArray() with the appropriate kernels.
    If the data are anisotropic (different pixel size along different dimensions)
    you should call \ref separableConvolveMultiArray() directly with the appropriate
    anisotropic Gaussian derivatives.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianGradientMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                   DestIterator diter, DestAccessor dest,
                                   double sigma);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                   pair<DestIterator, DestAccessor> const & dest,
                                   double sigma);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__convolution_8hxx-source.html">vigra/multi_convolution.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    ...
    // compute Gaussian gradient at scale sigma
    gaussianGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma);
    \endcode

    <b> Required Interface:</b>

    see \ref convolveImage(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see separableConvolveMultiArray()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
gaussianGradientMultiArray( SrcIterator si, SrcShape const & shape, SrcAccessor src,
			    DestIterator di, DestAccessor dest, double sigma )
{
    typedef typename DestAccessor::value_type DestType;
    typedef typename NumericTraits<typename DestType::value_type>::RealPromote kernel_type;

    Kernel1D<kernel_type> gauss, derivative;
    gauss.initGaussian(sigma);
    derivative.initGaussianDerivative(sigma, 1);

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;

    // compute gradient components
    for(unsigned int d = 0; d < shape.size(); ++d )
    {
        ArrayVector<Kernel1D<kernel_type> > kernels(shape.size(), gauss);
        kernels[d] = derivative;
        separableConvolveMultiArray( si, shape, src, di, ElementAccessor(d, dest), kernels.begin());
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                  pair<DestIterator, DestAccessor> const & dest, double sigma )
{
    gaussianGradientMultiArray( source.first, source.second, source.third,
                              dest.first, dest.second, sigma );
}

/********************************************************/
/*                                                      */
/*             symmetricGradientMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate gradient of a multi-dimensional arrays using symmetric difference filters.

    This function computes the gradient of the given multi-dimensional
    array with a sequence of symmetric difference filters a (differentiation is applied
    to each dimension in turn, starting with the innermost dimension). Both source and
    destination arrays are represented by iterators, shape objects and accessors.
    The destination array is required to have a vector valued pixel type with as many
    elements as the number of dimensions. This function is implemented by calls to
    \ref convolveMultiArrayOneDimension() with the symmetric difference kernel.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        symmetricGradientMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        symmetricGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="multi__convolution_8hxx-source.html">vigra/multi_convolution.hxx</a>"

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    ...
    // compute gradient
    symmetricGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
    \endcode

    <b> Required Interface:</b>

    see \ref convolveImage(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see convolveMultiArrayOneDimension()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
symmetricGradientMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                   DestIterator di, DestAccessor dest)
{
    typedef typename DestAccessor::value_type DestType;
    typedef typename NumericTraits<typename DestType::value_type>::RealPromote kernel_type;

    Kernel1D<kernel_type> filter;
    filter.initSymmetricGradient();

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;

    // compute gradient components
    for(unsigned int d = 0; d < shape.size(); ++d )
    {
        convolveMultiArrayOneDimension(si, shape, src,
				                       di, ElementAccessor(d, dest),
				                       d, filter);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
symmetricGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest )
{
    symmetricGradientMultiArray(source.first, source.second, source.third,
                                dest.first, dest.second);
}

//@}

} //-- namespace vigra


#endif        //-- VIGRA_MULTI_CONVOLUTION_H
