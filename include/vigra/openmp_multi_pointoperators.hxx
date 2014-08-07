#ifndef OPENMP_MULTI_POINTOPERATORS_H_
#define OPENMP_MULTI_POINTOPERATORS_H_

#include "multi_pointoperators.hxx"
#include "openmp_def.h"

namespace vigra
{
namespace omp
{
#ifdef OPENMP

template <unsigned int N, class T1, class S1,
class T2, class S2>
inline void
copyMultiArray(MultiArrayView<N, T1, S1> const & source,
		MultiArrayView<N, T2, S2> dest)
{
#pragma omp parallel
	{
#pragma omp for schedule(guided) nowait
		for(int y = 0; y < source.size(0); y++ )
		{
			vigra::copyMultiArray(source.bind<0>(y), dest-iter);
		}
	} // omp parallel
}
#else
template<unsigned int N, class T1, class S1, class T2, class S2>
inline void copyMultiArray(MultiArrayView<N, T1, S1> const & source,
		MultiArrayView<N, T2, S2> dest) {
	vigra::copyMultiArray(source, dest);
}
}
#endif //OPENMP
} //namespace omp
} //namespace vigra

#endif //OPENMP_MULTI_POINTOPERATORS_H_
