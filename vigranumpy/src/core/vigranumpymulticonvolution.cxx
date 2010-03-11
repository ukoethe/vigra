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
NumpyAnyArray pythonGaussianSmoothingND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                        double sigma,
                                        NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "gaussianSmoothingND(): Output array has wrong shape.");

    for(int k=0;k<volume.shape(ndim-1);++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	gaussianSmoothMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), sigma);
    }
    return res;

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
NumpyAnyArray pythonGaussianGradientMagnitudeND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                                double sigma,
                                                NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    using namespace vigra::functor;
    
    res.reshapeIfEmpty(volume.shape(), "gaussianGradientMagnitudeND(): Output array has wrong shape.");
    
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
NumpyAnyArray pythonSymmetricGradientND(NumpyArray<ndim, Singleband<VoxelType> > volume,
                                        double sigma,
                                        NumpyArray<ndim, TinyVector<VoxelType, (int)ndim> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "symmetricGradientND(): Output array has wrong shape.");
    symmetricGradientMultiArray(srcMultiArrayRange(volume), destMultiArray(res));
    return res;
}

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

void defineMultiConvolutionFunctions()
{
    using namespace python;

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothingND<float,3>),               
        (arg("image"), arg("sigma"), arg("out")=python::object()),
        "Perform isotropic Gaussian smoothing at the given scale for a 2D or 3D scalar or multiband image.\n"
        "\n"
        "For details see gaussianSmoothing_ in the vigra C++ documentation.\n");

    def("gaussianSmoothing",
        registerConverters(&pythonGaussianSmoothingND<float,4>),               
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
    	registerConverters(&pythonGaussianGradientMagnitudeND<float,3>),
    	(arg("image"), arg("sigma"), arg("out")=python::object()),
        "Calculate the gradient magnitude by means of a 1st derivative of "
        "Gaussian filter at the given scale for a 2D or 3D scalar or multiband image.\n"
        "\n"
        "For details see gaussianGradientMultiArray_ in the vigra C++ documentation.\n");

    def("gaussianGradientMagnitude",
    	registerConverters(&pythonGaussianGradientMagnitudeND<float,4>),
    	(arg("volume"), arg("sigma"), arg("out")=python::object()));

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

    def("separableConvolve", registerConverters(&pythonSeparableConvolveND_1Kernel<float,3>),
        (arg("image"), arg("kernel"), arg("out")=python::object()),
        "Apply the same 1D kernel to all dimensions of a 2D or 3D scalar or multiband image. "
        "'kernel' must be a single instance of Kernel1D or a tuple of Kernel1D objects, "
        "one kernel for each spatial dimension. If only a single kernel is given, "
        "it is applied to all dimensions.\n"
        "\n"
        "For details see separableConvolveMultiArray_ in the vigra C++ documentation.\n");

    def("separableConvolve", registerConverters(&pythonSeparableConvolveND_1Kernel<float,4>),
        (arg("volume"), arg("kernel"), arg("out")=python::object()));

    def("separableConvolve", registerConverters(&pythonSeparableConvolveND_NKernels<float,3>),
        (arg("image"), arg("kernels"), arg("out")=python::object()));

    def("separableConvolve", registerConverters(&pythonSeparableConvolveND_NKernels<float,4>),
        (arg("volume"), arg("kernels"), arg("out")=python::object()));
}

} // namespace vigra
