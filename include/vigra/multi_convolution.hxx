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

 //#define MC_SHOW_CONVOLUTION       //-- should not be uncommented

#include <vector>
#include <list>
#include <algorithm>
#include <iterator>

#include <vigra/separableconvolution.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/accessor.hxx>
#include <vigra/navigator.hxx>


namespace vigra
{


//---------------------------------------------------------------------------

/********************************************************/
/*                                                      */
/*           Some internal slave functions              */
/*                                                      */
/********************************************************/

namespace detail
{

    /** A small helper
     */
    template <typename Iterator, typename Accessor>
    inline triple<Iterator, Iterator, Accessor>
    srcRange( Iterator b, Iterator e, Accessor a)
    {
      return triple<Iterator, Iterator, Accessor>( b, e, a);
    }


    /** A small helper
     */
    template <typename Iterator, typename Accessor>
    pair<Iterator,Accessor>
    destRange( Iterator it, Accessor a )
    {
      return pair<Iterator,Accessor>( it, a );
    }


    //-----------------------------------------------------------------------


    /********************************************************/
    /*                                                      */
    /*        internalSeparableConvolveMultiarray           */
    /*                                                      */
    /********************************************************/

    template <class S, class D, typename KernelIterator>
    inline
    void internalSeparableConvolveMultiarrayTmp( S source, D dest,
						 KernelIterator kit )
    {
      typedef typename S::first_type SMIT;
      typedef typename D::first_type DMIT;
      enum { N = 1 + SMIT::level };

      typedef vigra::MultiArrayNavigator<SMIT, N> SNavigator;
      typedef vigra::MultiArrayNavigator<DMIT, N> DNavigator;

      { //-- only operate on first dimension here
#if defined MC_SHOW_CONVOLUTION
	std::cerr << "*Dim: " << 0 << std::endl;
#endif
        SNavigator snav( source.first, source.second, 0 );
	DNavigator dnav( dest.first, source.second, 0 );	

	for( ; snav.hasMore(); snav++, dnav++ ) {
 	  vigra::convolveLine( srcRange( snav.begin(), snav.end(), source.third ),
 			       destRange( dnav.begin(), dest.second ),
 			       kernel1d( *kit ) );
	}
	++kit;
      }

      //-- operate on further dimensions
      typedef typename NumericTraits<typename D::second_type::value_type>::RealPromote TmpType;
      std::vector<TmpType> tmp;

      for( int d = 1; d < N; ++d, ++kit ) {
	DNavigator dnav( dest.first, source.second, d );

	tmp.resize( source.second[d] );
#if defined MC_SHOW_CONVOLUTION
 	std::cerr << "*Dim: " << d << " = " << source.second[d]
		  << ", " << dnav.end() - dnav.begin()
		  << ", " << std::distance( dnav.begin(), dnav.end() )
		  << std::endl;
#endif

	for( ; dnav.hasMore(); dnav++ ) {
#if defined MC_SHOW_CONVOLUTION
	  std::cerr << ".";
#endif
 	  vigra::convolveLine( srcRange( dnav.begin(), dnav.end(), dest.second ),
 		 	       destRange( tmp.begin(), StandardValueAccessor<TmpType>() ),
 			       kernel1d( *kit ) );

	  //-- copy temp result to target object
	  vigra::copyLine( tmp.begin(), tmp.end(), StandardConstValueAccessor<TmpType>(),
			   dnav.begin(), dest.second );
	}
#if defined MC_SHOW_CONVOLUTION
	std::cerr << std::endl;
#endif
      }
    }


};	//-- namespace detail




// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/********************************************************/
/*                                                      */
/*             separableConvolveMultiarray              */
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
template <class S, class D, class KernelIterator>
inline
void separableConvolveMultiarray( S source, D dest, KernelIterator kit )
{
    vigra::detail::internalSeparableConvolveMultiarrayTmp( source, dest, kit );
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** \overload

      This overloaded function computes an isotropic convolution where
      the given kernel \a kernel is applied to each dimension
      individually.
   */
template <class S, class D, class A>
inline
void separableConvolveMultiarray( S source, D dest,
				    const vigra::Kernel1D<A> &kernel )
{
    typedef typename S::first_type SF;
    const int N = 1 + SF::level;
    
    std::list<const vigra::Kernel1D<A> *> kernelList( N );
    for( int i = 0; i < N; ++i )
      kernelList[i] = kernel;
  
    vigra::internalSeparableConvolveMultiarray( source, dest,
						kernelList.begin() );
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** \overload

      This overloaded function computes a convolution with the given
      kernel \a kernel only on the dimension \a dim of the MultiArray
      \a source.

      This may be used to e.g. smooth data only in one dimension.
   */
template <class S, class D, class A>
inline
void separableConvolveMultiarray( S source, D dest,
				    unsigned int dim,
				    const vigra::Kernel1D<A> &kernel )
{
    typedef typename S::first_type SF;
    const int N = 1 + SF::level;
    vigra_precondition( dim < N,
			"The dimension number to convolve must be smaller "
			"than the data dimensionality" );

#if 0
    //-- for the moment being this is a suboptimal
    //-- implementation in terms of speed.
    
    std::vector<vigra::Kernel1D<A> > kernels( N );
    kernels[dim] = kernel;
    vigra::separableConvolveMultiarray( source, dest,
					kernels.begin() );

#else
    typedef typename S::first_type SMIT;
    typedef typename D::first_type DMIT;

    typedef vigra::MultiArrayNavigator<SMIT, N> SNavigator;
    typedef vigra::MultiArrayNavigator<DMIT, N> DNavigator;

    SNavigator snav( source.first, source.second, dim );
    DNavigator dnav( dest.first, source.second, dim );

    for( ; snav.hasMore(); snav++, dnav++ ) {
      vigra::convolveLine( detail::srcRange( snav.begin(), snav.end(), source.third ),
			   detail::destRange( dnav.begin(), dest.second ),
			   kernel1d(kernel) );
    }
#endif
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
      separableConvolveMultiarray() function given a set of the
      respective filter kernels as a parameter.

      \see separableConvolveMultiarray()
*/
template <class S, class D>
inline
void gaussianSmoothing( S source, D dest, double sigma )
{
    typedef typename S::first_type SF;
    const unsigned int N = 1 + SF::level;
    
    vigra_precondition( sigma > 0,
			"The kernel standard deviation sigma "
			"must be greater than zero" );

    std::vector<vigra::Kernel1D<double> > kernels( N );
    for( unsigned int i = 0; i < N; ++i )
      kernels[i].initGaussian( sigma );

    vigra::separableConvolveMultiarray( source, dest, kernels.begin() );
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
    separableConvolveMultiarray() several times with appropriate
    derivative filters.

    The production of vector-valued output eases the operation on
    n-dimensional input data instead of n separate component data
    sets. Also, typically the whole vector is of interest, not only
    one component at a time. If this would be the case, simple
    separated convolution is preferablly used.

    As vector-valued element type a TinyVector of suiting dimension
    can be used. However, at least a VectorComponentAccessor must be
    applicable to the destination element data type.

    \see separableConvolveMultiarray()
*/
template <typename S, typename D>
inline
void symmetricGradient( S source, D dest )
{
    typedef typename S::first_type SF;
    const int N = 1 + SF::level;

    //  vigra_precondition( D::ActualDimension == N,
    //  vigra_precondition( D::value_type == N,

    //   vigra_precondition( (*dest).size() == N,
    // 		      "The pixel dimension of the destination data "
    // 		      "must match the data dimensions of the source data" );

    typedef typename NumericTraits<typename S::first_type::value_type>::RealPromote kernel_type;

    Kernel1D<kernel_type> filter;
    filter.initSymmetricGradient();

    typedef typename D::first_type DestType;
    typedef typename DestType::value_type VectorType;
    typedef VectorComponentAccessor<VectorType> VAccessor;
    
    //-- compute gradient components
    for( int d = 0; d < N; ++d ) {
      separableConvolveMultiarray( source,
				   std::pair<DestType, VAccessor>( dest.first, VAccessor(d) ),
				   d, filter );
    }
}



};	//-- namespace vigra


#endif	//-- VIGRA_MULTI_CONVOLUTION_H
