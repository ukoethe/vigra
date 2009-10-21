#ifndef VIGRANUMPYCONVOLUTION_HXX
#define VIGRANUMPYCONVOLUTION_HXX

#include <vigra/stdconvolution.hxx>
#include <boost/tuple/tuple.hpp>

namespace vigra
{

typedef double KernelValueType;
typedef Kernel2D< KernelValueType > TwoDKernel;
typedef Kernel1D< KernelValueType > Kernel;

template<class T>
struct kTriplet
{
	typedef boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > Type;
};

void exportVigraNumpyConvolution();

} // namespace vigra

#endif // VIGRANUMPYCONVOLUTION_HXX
