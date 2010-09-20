#ifndef VIGRANUMPYKERNEL_HXX
#define VIGRANUMPYKERNEL_HXX

#include <vigra/separableconvolution.hxx>
#include <vigra/stdconvolution.hxx>

namespace vigra
{

typedef double KernelValueType;
typedef Kernel2D< KernelValueType > TwoDKernel;
typedef Kernel1D< KernelValueType > Kernel;

} // namespace vigra

#endif // VIGRANUMPYKERNEL_HXX
