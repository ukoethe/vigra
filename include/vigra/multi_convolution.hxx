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
inline void 
internalSeparableConvolveMultiArrayTmp( 
                      SrcIterator si, SrcShape const & shape, SrcAccessor src, 
                      DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp;

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );        

        tmp.resize( shape[0] );

        for( ; snav.hasMore(); snav++, dnav++ ) 
        {
             convolveLine( srcIterRange( snav.begin(), snav.end(), src ),
                           destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
                           kernel1d( *kit ) );
             // copy temp result to target object
             copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
                       dnav.begin(), dest );
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
             convolveLine( srcIterRange( dnav.begin(), dnav.end(), dest ),
                           destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
                           kernel1d( *kit ) );

             // copy temp result to target object
             copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
                       dnav.begin(), dest);
        }
    }
}


} // namespace detail


/********************************************************/
/*                                                      */
/*             separableConvolveMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Separated convolution on MultiArray

      This function computes a separated convolution on all dimensions
      of the given MultiArray \a source. The result will be stored in
      \a dest. \a dest is required to already have the correct size.

      The argument \a kernels is an iterator referencing 1D kernels, each
      for every dimension of the data. The first kernel in this list
      is applied to the innermost dimension (e.g. the x-dimension of
      an image), whilst the last is applied to the outermost dimension
      (e.g. the z-dimension in a 3D image).

      In order to reduce memory consumption this uses an internal
      buffer only of the size of the largest dimension of the source
      object. However, this also results in a slightly reduced
      execution speed.

      \note This function may work in-place, which means that \a
      source == \a dest is allowed.

      \note This function is optimized in terms of memory consumption
      and not for speed.  It allocs temporary memory only of minimum
      size.

      \see Kernel1D, convolveLine()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class KernelIterator>
inline void 
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src, 
                             DestIterator d, DestAccessor dest, KernelIterator kit )
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    if(typeid(TmpType) != typeid(typename DestAccessor::value_type))
    {
        // need a temporary array to avoid rounding errors
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src, 
             tmpArray.traverser_begin(), typename AccessorTraits<TmpType>::default_accessor(), kit );
        copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
    }
    else
    {
        // work directly on the destination array
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src, d, dest, kit );
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



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** \overload

      This overloaded function computes an isotropic convolution where
      the given kernel \a kernel is applied to each dimension
      individually.
   */
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
    ArrayVector<Kernel1D<T> > kernels(shape.size(), kernel);
  
    separableConvolveMultiArray( source.first, source.second, source.third, 
                                 dest.first, dest.second, kernels.begin() );
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** \overload

      This overloaded function computes a convolution with the given
      kernel \a kernel only on the dimension \a dim of the MultiArray
      \a source.

      This may be used to e.g. smooth data only in one dimension.
   */
template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class T>
inline void 
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src, 
                             DestIterator d, DestAccessor dest, 
                             unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    enum { N = 1 + SrcIterator::level };
    vigra_precondition( dim < N,
                        "separableConvolveMultiArray(): The dimension number to convolve must be smaller "
                        "than the data dimensionality" );

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    
    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp(shape[dim]);

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    SNavigator snav( s, shape, dim );
    DNavigator dnav( d, shape, dim );        

    for( ; snav.hasMore(); snav++, dnav++ ) 
    {
         convolveLine( srcIterRange( snav.begin(), snav.end(), src ),
                       destIter( tmp.begin(), StandardValueAccessor<TmpType>() ),
                       kernel1d( kernel ) );
         // copy temp result to target object
         copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
                   dnav.begin(), dest );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor, class T>
inline void 
separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source, 
                            pair<DestIterator, DestAccessor> const & dest, 
                                  unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    separableConvolveMultiArray( source.first, source.second, source.third, 
                                 dest.first, dest.second, dim, kernel );
}


/********************************************************/
/*                                                      */
/*                  gaussianSmoothing                   */
/*                                                      */
/********************************************************/


/** \brief Isotropic Gaussian smoothing of a MultiArrayView

      This performs an isotropic smoothing using a Gaussian filter
      kernel on the given MultiArrayView data \a source. The result
      will be stored in \a dest.

      Both arrays, \a source and \a dest, must be equal in dimensions
      and size.

      An anisotropic smoothing can be realized using the
      separableConvolveMultiArray() function given a set of the
      respective filter kernels as a parameter.

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



/** \brief Gaussian smoothed gradient

    This function computes the gradient of the data
    \a source and stores the gradient in the vector-valued \a
    dest.

    The data is differentiated by symmetric differences for each
    dimension.

    In order to apply this function to anisotropic data, one should
    consider first applying an anisotropic smoothing before calling
    this function. Another option is calling
    separableConvolveMultiArray() several times with appropriate
    derivative filters.

    The production of vector-valued output eases the operation on
    n-dimensional input data instead of n separate component data
    sets. Also, typically the whole vector is of interest, not only
    one component at a time. If this would be the case, simple
    separated convolution is preferablly used.

    As vector-valued element type a TinyVector of suiting dimension
    can be used. However, at least a VectorComponentAccessor must be
    applicable to the destination element data type.

    \see separableConvolveMultiArray()
*/
template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void 
gaussianGradientMultiArray( SrcIterator si, SrcShape const & shape, SrcAccessor src, 
                   DestIterator di, DestAccessor dest, double sigma )
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote kernel_type;

    Kernel1D<kernel_type> gauss, derivative;
    gauss.initGaussian(sigma);
    derivative.initGaussianDerivative(sigma, 1);

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;
    
    // compute gradient components
    for( int d = 0; d < shape.size(); ++d ) 
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

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void 
symmetricGradientMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src, 
                   DestIterator di, DestAccessor dest)
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote kernel_type;

    Kernel1D<kernel_type> filter;
    filter.initSymmetricGradient();

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;
    
    // compute gradient components
    for( int d = 0; d < shape.size(); ++d ) 
    {
        separableConvolveMultiArray( si, shape, src, di, ElementAccessor(d, dest), d, filter);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
symmetricGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source, 
                  pair<DestIterator, DestAccessor> const & dest )
{
    symmetricGradientMultiArray( source.first, source.second, source.third, 
                              dest.first, dest.second );
}



} //-- namespace vigra


#endif        //-- VIGRA_MULTI_CONVOLUTION_H
