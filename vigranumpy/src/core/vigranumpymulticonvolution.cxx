/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyconvolution_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/functorexpression.hxx>
#include "vigranumpykernel.hxx"

namespace python = boost::python;

namespace vigra
{

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonConvolveOneDimensionND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                           unsigned int dim,
                                           Kernel const & kernel,
                                           NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "convolveOneDimensionND(): Output array has wrong shape.");
    for(int k=0;k<volume.shape(ndim-1);++k)
	{
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	convolveMultiArrayOneDimension(srcMultiArrayRange(bvolume), destMultiArray(bres), dim, kernel);
	}
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonSeparableConvolveND_1Kernel(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                                Kernel const & kernel,
                                                NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "separableConvolveND(): Output array has wrong shape.");
    for(int k=0;k<volume.shape(ndim-1);++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	separableConvolveMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), kernel);
    }
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonSeparableConvolveND_NKernels(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                                 python::tuple pykernels,
                                                 NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    if(python::len(pykernels) == 1)
        return pythonSeparableConvolveND_1Kernel(volume, python::extract<Kernel1D<KernelValueType> const &>(pykernels[0]), res);
        
    vigra_precondition(python::len(pykernels) == ndim-1,
       "separableConvolveND(): Number of kernels must be 1 or equal to the number of spatial dimensions.");
       
    res.reshapeIfEmpty(volume.shape(), "separableConvolveND(): Output array has wrong shape.");
    
    ArrayVector<Kernel1D<KernelValueType> > kernels;
    for(unsigned int k=0; k < ndim-1; ++k)
        kernels.push_back(python::extract<Kernel1D<KernelValueType> const &>(pykernels[k]));

    for(int k=0; k < volume.shape(ndim-1); ++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	separableConvolveMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), kernels.begin());
    }
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonGaussianSmoothing(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                        python::tuple sigmas,
                                        NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    unsigned int sigmaCount = python::len(sigmas);
    vigra_precondition(sigmaCount == 1 || sigmaCount == ndim-1,
       "gaussianSmoothing(): Number of kernels must be 1 or equal to the number of spatial dimensions.");
       
    ArrayVector<Kernel1D<KernelValueType> > kernels;
    for(unsigned int k=0; k < sigmaCount; ++k)
    {
        KernelValueType sigma = python::extract<KernelValueType>(sigmas[k]);        
        kernels.push_back(Kernel1D<KernelValueType>());
        kernels.back().initGaussian(sigma);
    }
    for(unsigned int k=sigmaCount; k < ndim-1; ++k)
    {
        kernels.push_back(Kernel1D<KernelValueType>(kernels.back()));
    }
    
    res.reshapeIfEmpty(volume.shape(), "gaussianSmoothing(): Output array has wrong shape.");

    for(int k=0;k<volume.shape(ndim-1);++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	separableConvolveMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), kernels.begin());
    }
    return res;

}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonGaussianSmoothingIsotropic(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                               double sigma,
                                               NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    return pythonGaussianSmoothing(volume, python::make_tuple(sigma), res);
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonGaussianGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                       double sigma,
                                       NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "gaussianGradientND(): Output array has wrong shape.");
    gaussianGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res), sigma);
    return res;
}


template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                  double sigma, 
                                  NumpyArray<ndim-1, Singleband<VoxelType> > res=python::object())
{
    using namespace vigra::functor;
    
    typename MultiArrayShape<ndim-1>::type tmpShape(volume.shape().begin());
    res.reshapeIfEmpty(tmpShape, "gaussianGradientMagnitude(): Output array has wrong shape.");
    res.init(NumericTraits<VoxelType>::zero());
    MultiArray<ndim-1, TinyVector<VoxelType, (int)(ndim-1)> > grad(tmpShape);
    for(int k=0; k<volume.shape(ndim-1); ++k)
    {
        MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(grad), sigma);
        combineTwoMultiArrays(srcMultiArrayRange(grad), srcMultiArray(res), destMultiArray(res), 
                              squaredNorm(Arg1())+Arg2());
    }
    transformMultiArray(srcMultiArrayRange(res), destMultiArray(res), sqrt(Arg1()));
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                  double sigma, 
                                  NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    using namespace vigra::functor;
    
    res.reshapeIfEmpty(volume.shape(), "gaussianGradientMagnitude(): Output array has wrong shape.");
    
    typename MultiArrayShape<ndim-1>::type tmpShape(volume.shape().begin());
    MultiArray<ndim-1, TinyVector<VoxelType, (int)(ndim-1)> > grad(tmpShape);
    for(int k=0; k<volume.shape(ndim-1); ++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(grad), sigma);
        transformMultiArray(srcMultiArrayRange(grad), destMultiArray(bres), norm(Arg1()));
    }
    return res;
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray 
pythonGaussianGradientMagnitude(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                double sigma, bool accumulate, 
                                NumpyAnyArray res=python::object())
{
    return accumulate
              ? pythonGaussianGradientMagnitudeND(volume, sigma, NumpyArray<ndim-1, Singleband<VoxelType> >(res))
              : pythonGaussianGradientMagnitudeND(volume, sigma, NumpyArray<ndim, Multiband<VoxelType> >(res));
}

template < class VoxelType, unsigned int ndim >
NumpyAnyArray pythonSymmetricGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                        double sigma,
                                        NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "symmetricGradientND(): Output array has wrong shape.");
    symmetricGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res));
    return res;
}

void defineMultiConvolutionFunctions()
{
    using namespace python;

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothingIsotropic<float,3>),               
        (arg("image"), arg("sigma"), arg("out")=python::object()),
        "Perform Gaussian smoothing of a 2D or 3D scalar or multiband image.\n\n"
        "Each channel of the array is smoothed independently. "
        "If 'sigma' is a single value, an isotropic Gaussian filter at this scale is "
        "applied (i.e. each dimension is smoothed in the same way). "
        "If 'sigma' is a tuple of values, the amount of smoothing will be different "
        "for each spatial dimension. The length of the tuple must be equal to the "
        "number of spatial dimensions.\n\n"        
        "For details see gaussianSmoothing_ in the vigra C++ documentation.\n");

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothing<float,3>),               
        (arg("volume"), arg("sigma"), arg("out")=python::object()));

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothingIsotropic<float,4>),               
        (arg("volume"), arg("sigma"), arg("out")=python::object()));

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothing<float,4>),               
        (arg("volume"), arg("sigma"), arg("out")=python::object()));

    def("gaussianGradient",
    	registerConverters(&pythonGaussianGradientND<float,2>),
    	(arg("image"), arg("sigma"), arg("out")=python::object()),
        "Calculate the gradient vector by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D or 3D scalar image.\n"
        "\n"
        "For details see gaussianGradientMultiArray_ in the vigra C++ documentation.\n");

    def("gaussianGradient",
    	registerConverters(&pythonGaussianGradientND<float,3>),
    	(arg("volume"), arg("sigma"), arg("out")=python::object()));

    def("gaussianGradientMagnitude",
    	registerConverters(&pythonGaussianGradientMagnitude<float,3>),
    	(arg("image"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object()),
        "Calculate the gradient magnitude by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D or 3D scalar or multiband image.\n"
        "If 'accumulate' is True (the default), the gradients are accumulated (in the "
        "L2-norm sense) over all  channels of a multi-channel array. Otherwise, "
        "a separate gradient magnitude is computed for each channel.\n"
        "\n"
        "For details see gaussianGradientMultiArray_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
    	registerConverters(&pythonGaussianGradientMagnitude<float,4>),
    	(arg("volume"), arg("sigma"), arg("accumulate")=true, arg("out")=python::object()));

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,2>),
        (arg("image"), arg("out")=python::object()),
        "Calculate gradient of a scalar 2D image or 3D volume using symmetric difference filters."
        "\n"
        "For details see symmetricGradientMultiArray_ in the vigra C++ documentation.\n");

    def("symmetricGradient",
        registerConverters(&pythonSymmetricGradientND<float,3>), 
        (arg("volume"), arg("out")=python::object()));

    def("convolveOneDimension",
    	registerConverters(&pythonConvolveOneDimensionND<float,3>),
    	(arg("image"), arg("dim"), arg("kernel"), arg("out")=python::object()),
        "Convolution along a single dimension of a 2D or 3D scalar or multiband image. "
        "'kernel' must be an instance of Kernel1D.\n"
        "\n"
        "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("convolveOneDimension",
    	registerConverters(&pythonConvolveOneDimensionND<float,4>),
    	(arg("volume"), arg("dim"), arg("kernel"), arg("out")=python::object()));

    def("convolve", registerConverters(&pythonSeparableConvolveND_1Kernel<float,3>),
        (arg("image"), arg("kernel"), arg("out")=python::object()),
        "Convolve an image or volume with the given 'kernel' (or kernels).\n"
        "If the input has multiple channels, the filter is applied to each channel\n"
        "independently. The function can be used in 3 different ways:\n"
        "\n"
        "* When 'kernel' is a single object of type :class:`Kernel1D`, this kernel\n"
        "  is applied along all spatial dimensions of the data (separable filtering).\n"
        "* When 'kernel' is a tuple of :class:`Kernel1D` objects, one different kernel\n"
        "  is used for each spatial dimension (separable filtering). The number of\n"
        "  kernels must equal the number of dimensions).\n"
        "* When 'kernel' is an instance of :class:`Kernel2D`, a 2-dimensional convolution\n"
        "  is performed (non-separable filtering). This is only applicable to 2D images.\n"
        "\n"
        "For details see separableConvolveMultiArray_ and "
        "|StandardConvolution.convolveImage|_ in the vigra C++ documentation.\n\n");

    def("convolve", registerConverters(&pythonSeparableConvolveND_1Kernel<float,4>),
        (arg("volume"), arg("kernel"), arg("out")=python::object()));

    def("convolve", registerConverters(&pythonSeparableConvolveND_NKernels<float,3>),
        (arg("image"), arg("kernels"), arg("out")=python::object()));

    def("convolve", registerConverters(&pythonSeparableConvolveND_NKernels<float,4>),
        (arg("volume"), arg("kernels"), arg("out")=python::object()));
}

} // namespace vigra
