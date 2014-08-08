#ifndef OPENMP_MULTI_POINTOPERATORS_H_
#define OPENMP_MULTI_POINTOPERATORS_H_

#include "multi_pointoperators.hxx"
#include "openmp_def.h"

namespace vigra
{
namespace omp
{
#ifdef OPENMP

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(SrcIterator s,
               SrcShape const & shape, SrcAccessor src,
               DestIterator d, DestAccessor dest)
{
#pragma omp parallel
	{
		//Compute size
		shape.size() //?
#pragma omp for scheduled(guided) nowait
		for ()
		{
			//For each thread: compute its reponsible ROI
			vigra::copyMultiArray(s+begin, shape, src, d+begin, dest);
		}
	}//omp parallel
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void
copyMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest)
{

}

#else
template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(SrcIterator s,
               SrcShape const & shape, SrcAccessor src,
               DestIterator d, DestAccessor dest)
{
    vigra::copyMultiArray(s, shape, src, d, dest);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
void
copyMultiArray(SrcIterator s, SrcShape const & sshape, SrcAccessor src,
               DestIterator d, DestShape const & dshape, DestAccessor dest)
{
    vigra::copyMultiArray(s, sshape, src, d, dshape, dest);
}
#endif //OPENMP

//
// Argument Object Factory versions
//

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               pair<DestIterator, DestAccessor> const & dest)
{

	vigra::omp::copyMultiArray(src.first, src.second, src.third, dest.first, dest.second);
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestShape, class DestAccessor>
inline void
copyMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & src,
               triple<DestIterator, DestShape, DestAccessor> const & dest)
{

	vigra::omp::copyMultiArray(src.first, src.second, src.third, dest.first, dest.second, dest.third);
}

} //namespace omp
} //namespace vigra

#endif //OPENMP_MULTI_POINTOPERATORS_H_
