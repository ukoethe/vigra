//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Christian-Dennis Rahn                */
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

#include <vector>
#include <list>
#include <algorithm>

#include <vigra/separableconvolution.hxx>
#include <vigra/multi_array.hxx>

namespace vigra
{

  /********************************************************/
  /*                                                      */
  /*        internalSeparableConvolveMultiarray           */
  /*                                                      */
  /********************************************************/

  template <class T, class C, class A>
  inline
  void internalSeparableConvolveMultiarray( const MultiArrayView<1,T,C> &src,
					    MultiArrayView<1,T,C> &dst,
					    std::list<const vigra::Kernel1D<A> *> &kernels,
					    std::vector<T> &tmp )
  {
    const int len = src.shape(0);

    tuple5< typename Kernel1D<A>::const_iterator,
      typename Kernel1D<A>::ConstAccessor,
      int, int, BorderTreatmentMode> kernel = vigra::kernel1d( * kernels.front() );


    vigra::convolveLine( src.traverser_begin(), src.traverser_end(),
			 StandardConstAccessor<T>(),
			 tmp.begin(), StandardAccessor<T>(),
			 kernel.first, kernel.second, kernel.third,
			 kernel.fourth, kernel.fifth );

    std::copy( tmp.begin(), tmp.begin() + len,
	       dst.traverser_begin() );

//     vigra::convolveLine( srcIterRange( src.traverser_begin(), src.traverser_end(),
// 				       vigra::StandardConstAccessor<T>() ),
// 			 destIter( dst.traverser_begin() ),
// 			 kernel );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <unsigned int N, class T, class C, class A>
  inline
  void internalSeparableConvolveMultiarray( const MultiArrayView<N,T,C> &source,
					    MultiArrayView<N,T,C> &dest,
					    std::list<const vigra::Kernel1D<A> *> &kernels,
					    std::vector<T> &tmp )
  {
    const int dimensions = N;
    const MultiArrayView<N,T,C> *src = &source;
    MultiArrayView<N,T,C> *dst = &dest;

    typedef typename std::list<const vigra::Kernel1D<A> *>::iterator KIterator;

    KIterator kit = kernels.begin();

    for( int d = 0; d < dimensions; ++d, ++kit ) {
      const Kernel1D<A> *ktmp = (*kit);
      kit = kernels.erase( kit );

      const int dsize = src->shape( d );      
      for( int i = 0; i < dsize; ++i ) {
	  
	MultiArrayView<N-1,T,StridedArrayTag> hyposrc = src->bindAt( d, i );
	MultiArrayView<N-1,T,StridedArrayTag> hypodst = dst->bindAt( d, i );
	
	vigra::internalSeparableConvolveMultiarray( hyposrc, hypodst, kernels, tmp );
      }
      src = dst;
      kit = kernels.insert( kit, ktmp );
    }
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /** \brief Separated convolution on MultiArray

      This function computes a separated convolution on all dimensions
      of the given MultiArray \a source. The result will be stored in
      \a dest.

      The parameter \a kernels is a vector containing 1D kernels, each
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
  template <unsigned int N, class T, class C, class A>
  inline
  void separableConvolveMultiarray( const MultiArrayView<N,T,C> &source,
				    MultiArrayView<N,T,C> &dest,
				    const std::vector<vigra::Kernel1D<A> > &kernels )
{
  vigra_precondition( kernels.size() >= N, "The kernels vector must have a length of"
		      " at least the data dimensionality" );

  const size_t tsize = *std::max_element( source.shape().begin(),
					  source.shape().end() );
  std::vector<T> tmp( tsize );

  std::list<const vigra::Kernel1D<A> *> kernelList;
  for( int i = 0; i < N; ++i )
    kernelList.push_back( &kernels[i] );

  vigra::internalSeparableConvolveMultiarray( source, dest, kernelList, tmp );
}



  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /** \overload
   */
  template <unsigned int N, class T, class C, class A>
  inline
  void separableConvolveMultiarray( const MultiArrayView<N,T,C> &source,
				    MultiArrayView<N,T,C> &dest,
				    const vigra::Kernel1D<A> &kernel )
{
  const size_t tsize = *std::max_element( source.shape().begin(),
					  source.shape().end() );
  std::vector<T> tmp( tsize );
  
  std::list<const vigra::Kernel1D<A> *> kernelList;
  for( int i = 0; i < N; ++i )
    kernelList.push_back( &kernel );
  
  vigra::internalSeparableConvolveMultiarray( source, dest, kernelList, tmp );
}



};	//-- namespace vigra


#endif	//-- VIGRA_MULTI_CONVOLUTION_H
