#ifndef VIGRA_MATLAB_TENSORS_HXX
#define VIGRA_MATLAB_TENSORS_HXX

#include <vigra/mathutil.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/matrix.hxx>
#include <vigra/eigensystem.hxx>
#include <vigra/array_vector.hxx>
#include <vigra/static_assert.hxx>

namespace vigra {

template <int N, int N1>
struct Tensors__Dimension_mismatch_between_input_and_output
: staticAssert::AssertBool<(N1 == N + 1)>
{};

template <class T>
struct GradMagFunctor
{
  mutable T sum;
  GradMagFunctor()
    : sum(0)
  {}

  void operator()(T t) const { sum += t*t; }

  T operator()() const { return std::sqrt(sum); }
};

template <class T>
struct FunctorTraits<GradMagFunctor<T> > 
{
  typedef VigraTrueType isInitializer;
  typedef VigraTrueType isUnaryAnalyser;
};

template <int N, class T, class C1, class U, class C2>
void gradientMagnitude(MultiArrayView<N, T, C1> const & in,
                       MultiArrayView<N, U, C2> & mag, double scale)
{
	for(int k=0; k<N; ++k)
		vigra_precondition(in.shape(k) == mag.shape(k),
	    	"gradientMagnitude(): Shape mismatch between input and output.");
	vigra_precondition(scale > 0.0,
	    "gradientMagnitude(): Scale must be positive.");

  Kernel1D<double> gauss;
	gauss.initGaussian(scale);

  typename MultiArrayShape<N+1>::type g_shape, g_stride;
	for(int k=0; k<N; ++k)
  {
     g_shape[k] = in.shape(k);
     g_stride[k] = mag.stride(k);
  }
  g_shape[N] = N;
  MultiArray<N+1, T> grad(g_shape); 

	for(int b=0; b<N; ++b)
	{
			MultiArrayView<N, T, UnstridedArrayTag> gband = grad.bindOuter(b);
			ArrayVector<Kernel1D<double> > kernels(N, gauss);
			kernels[b].initGaussianDerivative(scale, 1);
			separableConvolveMultiArray(srcMultiArrayRange(in), destMultiArray(gband),
										              kernels.begin());
	}

  g_shape[N] = 1;
  g_stride[N] = 1;
  MultiArrayView<N+1, U, UnstridedArrayTag> magv(g_shape, g_stride, mag.data());
  transformMultiArray(srcMultiArrayRange(grad), destMultiArrayRange(magv), GradMagFunctor<U>());
}

template <unsigned int N, class T, class C1, unsigned int N1, class U, class C2>
void hessianOfGaussian(MultiArrayView<N, T, C1> const & in,
                       MultiArrayView<N1, U, C2> & hesse, double scale)
{
	const int matSize = N*(N+1)/2;

	VIGRA_STATIC_ASSERT((Tensors__Dimension_mismatch_between_input_and_output<N, N1>));

	for(int k=0; k<N; ++k)
		vigra_precondition(in.shape(k) == hesse.shape(k),
	    	"hessianOfGaussian(): Shape mismatch between input and output.");
	vigra_precondition(hesse.shape(N) == matSize,
	    "hessianOfGaussian(): Wrong number of bands in output array.");
	vigra_precondition(scale > 0.0,
	    "hessianOfGaussian(): Scale must be positive.");

  Kernel1D<double> gauss;
	gauss.initGaussian(scale);

	for(int b=0, i=0; i<N; ++i)
	{
		for(int j=i; j<N; ++j, ++b)
		{
			MultiArrayView<N, U, C2> hband = hesse.bindOuter(b);
			ArrayVector<Kernel1D<double> > kernels(N, gauss);
			if(i == j)
			{
				kernels[i].initGaussianDerivative(scale, 2);
			}
			else
			{
				kernels[i].initGaussianDerivative(scale, 1);
				kernels[j].initGaussianDerivative(scale, 1);
			}
			separableConvolveMultiArray(srcMultiArrayRange(in), destMultiArray(hband),
										kernels.begin());
		}
	}
}

template <unsigned int N, class T, class C1, unsigned int N1, class U, class C2>
void structureTensor(MultiArrayView<N, T, C1> const & in,
                     MultiArrayView<N1, U, C2> & st, double innerScale, double outerScale)
{
	const int matSize = N*(N+1)/2;

	VIGRA_STATIC_ASSERT((Tensors__Dimension_mismatch_between_input_and_output<N, N1>));

	for(int k=0; k<N; ++k)
		vigra_precondition(in.shape(k) == st.shape(k),
	    	"structureTensor(): Shape mismatch between input and output.");
	vigra_precondition(st.shape(N) == matSize,
	    "structureTensor(): Wrong number of bands in output array.");
	vigra_precondition(innerScale > 0.0 && outerScale >= 0.0,
	    "structureTensor(): Scales must be positive.");
	    
	typename MultiArrayShape<N+1>::type gradShape(st.shape(0), st.shape(1), N);
	MultiArray<N+1, double> gradient(gradShape);

  Kernel1D<double> gauss;
	gauss.initGaussian(innerScale);

	for(int b=0; b<N; ++b)
	{
    MultiArrayView<N, double, UnstridedArrayTag> gband = gradient.bindOuter(b);
    ArrayVector<Kernel1D<double> > kernels(N, gauss);
    kernels[b].initGaussianDerivative(innerScale, 1);
		separableConvolveMultiArray(srcMultiArrayRange(in), destMultiArray(gband),
                                    kernels.begin());
	}



	for(int y=0; y<st.shape(1); ++y)
	{
		for(int x=0; x<st.shape(0); ++x)
		{
			for(int b=0, i=0; i<N; ++i)
			{
				for(int j=i; j<N; ++j, ++b)
				{
					st(x, y, b) = detail::RequiresExplicitCast<T>::cast(gradient(x, y, i)*gradient(x, y, j));
				}
			}
		}
	}

	
	for(int b=0; b<matSize; ++b)
	{
		MultiArrayView<N, U, C2> stband = st.bindOuter(b);
		
		gaussianSmoothMultiArray(srcMultiArrayRange(stband), destMultiArray(stband),
                              	 outerScale);
  }
}

template <class T, class C1, class C2>
void eigenValuesPerPixel(MultiArrayView<3, T, C1> const & tensors,
                         MultiArrayView<3, T, C2> & eigenValues)
{
	const int N = 2;
	const int matSize = N*(N+1)/2;

	for(int k=0; k<N; ++k)
		vigra_precondition(tensors.shape(k) == eigenValues.shape(k),
	    	"eigenValuesPerPixel(): Shape mismatch between input and output.");
	vigra_precondition(tensors.shape(N) == matSize,
	    "eigenValuesPerPixel(): Wrong number of bands in input array.");
	vigra_precondition(eigenValues.shape(N) == N,
	    "eigenValuesPerPixel(): Wrong number of bands in output array.");
	    
	Matrix<double> tensor(N, N), ev(N, 1);
	for(int y=0; y<tensors.shape(1); ++y)
	{
		for(int x=0; x<tensors.shape(0); ++x)
		{
			symmetric2x2Eigenvalues(tensors(x,y,0), tensors(x,y,1), tensors(x,y,2), 
                              &eigenValues(x, y, 0), &eigenValues(x, y, 1));
		}
	}
}

} // namespace vigra

#endif // VIGRA_MATLAB_TENSORS_HXX
