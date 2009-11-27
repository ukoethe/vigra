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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycmodule_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/multi_impex.hxx>
#include "vigranumpykernel.hxx"

namespace python = boost::python;

namespace vigra
{

template < class VoxelType, int ndim >
NumpyAnyArray pythonGaussianSmoothingND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                          double sigma,
                                          NumpyArray<ndim, Multiband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "gaussianSmoothingND(): Output array has wrong shape.");
    //gaussianSmoothMultiArray(srcMultiArrayRange(MultiArrayView<ndim, float, StridedArrayTag>(volume)),
    		                 //destMultiArray((MultiArrayView<ndim, float, StridedArrayTag>&)res)), sigma);
    for(int k=0;k<volume.shape(ndim-1);++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	gaussianSmoothMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), sigma);
    }
    return res;

}

template < class VoxelType, int ndim >
NumpyAnyArray pythonGaussianGradientND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                          double sigma,
                                          NumpyArray<ndim, Multiband< TinyVector<VoxelType, ndim> > > res=python::object())
{
	vigra_precondition(typename NumericTraits< VoxelType >::isScalar().asBool,
	        "Only single band arrays are supported.");
    res.reshapeIfEmpty(volume.shape(), "gaussianGradientND(): Output array has wrong shape.");
    for(int k=0;k<volume.shape(ndim-1);++k)
    {
        MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
        MultiArrayView<ndim-1, TinyVector<VoxelType, ndim>, StridedArrayTag> bres = res.bindOuter(k);
        gaussianGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), sigma);
    }
    return res;
}

template < class VoxelType, int ndim >
NumpyAnyArray pythonSymmetricGradientND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                          double sigma,
                                          NumpyArray<ndim, Multiband< TinyVector<VoxelType, ndim> > > res=python::object())
{
    vigra_precondition(typename NumericTraits< VoxelType >::isScalar().asBool,
        "Only single band arrays are supported.");
    res.reshapeIfEmpty(volume.shape(), "symmetricGradientND(): Output array has wrong shape.");
    for(int k=0;k<volume.shape(ndim-1);++k)
    {
    	MultiArrayView<ndim-1, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<ndim-1, TinyVector<VoxelType, ndim>, StridedArrayTag> bres = res.bindOuter(k);
    	symmetricGradientMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres));
    }
    return res;
}

template < class VoxelType, int ndim >
NumpyAnyArray pythonConvolveOneDimensionND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                          int dim,
                                          Kernel kernel,
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

template < class VoxelType, int ndim >
NumpyAnyArray pythonSeparableConvolveND(NumpyArray<ndim, Multiband<VoxelType> > volume,
                                          Kernel kernel,
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

// Todo: more cases for separableConvolve.
// what happens, if separableConvolveMultiArray gets called for data with dim > nkernels?

template < class VoxelType >
NumpyAnyArray pythonSeparableConvolve_3_3D(NumpyArray<4, Multiband<VoxelType> > volume,
                                          Kernel kernel0,
                                          Kernel kernel1,
                                          Kernel kernel2,
                                          NumpyArray<4, Multiband<VoxelType> > res=python::object())
{
    res.reshapeIfEmpty(volume.shape(), "separableConvolve_3_3D(): Output array has wrong shape.");
    ArrayVector< Kernel > kernels(3);
    kernels.push_back(kernel0);
    kernels.push_back(kernel1);
    kernels.push_back(kernel2);
    for(int k=0;k<volume.shape(3);++k)
	{
    	MultiArrayView<3, VoxelType, StridedArrayTag> bvolume = volume.bindOuter(k);
    	MultiArrayView<3, VoxelType, StridedArrayTag> bres = res.bindOuter(k);
    	separableConvolveMultiArray(srcMultiArrayRange(bvolume), destMultiArray(bres), kernels.begin());
	}
    return res;
}

#if 0

template < class VoxelType, int pyArrayTypeConstant, class Tag >
PyObject* separableConvolve2(MultiArrayView< 2, VoxelType, Tag > const
    & volume, Kernel const & kernel0,
    Kernel const & kernel1)
{
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< VoxelType >::isScalar());
    MultiArrayView< 2, VoxelType, UnstridedArrayTag >* arrayView =
        createView< VoxelType, 2 >(array,
        typename NumericTraits< VoxelType >::isScalar());
    ArrayVector< Kernel > kernels(2);
    kernels.push_back(kernel0);
    kernels.push_back(kernel1);
    separableConvolveMultiArray(srcMultiArrayRange(volume),
        destMultiArray(*arrayView), kernels.begin());
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

void pyKernel_initExplicitly(Kernel & self, int left, int right,
    python::object const & args)
{
    vigra_precondition(left <= 0, "left should be <= 0");
    vigra_precondition(right >= 0, "right should be >= 0");

    if(! PySequence_Check(args.ptr()))
    {
        KernelValueType value = python::extract<KernelValueType>(args);
        self.initExplicitly(left, right) = value;
    }
    else
    {
        KernelValueType value = python::extract<KernelValueType>(args[0]);
        Kernel::InitProxy ip = self.initExplicitly(left, right) = value;
        if(python::len(args) != self.size())
        {
            std::stringstream str;
            str << "Wrong number of init values. The number must be ";
            str << self.size();
            PyErr_SetString(PyExc_ValueError, str.str().c_str());
            python::throw_error_already_set();
        }
        else
        {
            int size = self.size();
            for(int i=1; i<size; ++i)
            {
                ip,(python::extract<KernelValueType>(args[i]));
            }
        }
    }
}

void pyKernel_copy(Kernel & self, Kernel const & kernel)
{
    self=kernel;
}

KernelValueType pyKernel_getitem(Kernel const & self, int position)
{
    if(self.left() <= position && self.right() >= position)
    {
        return self[position];
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.left() << " <= position <= " << self.right();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
        return 0;
    }
}

void pyKernel_setitem(Kernel & self, int position, KernelValueType value)
{
    if(self.left() <= position && self.right() >= position)
    {
        self[position] = value;
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.left() << " <= position <= " << self.right();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
}

#endif

namespace { void DUMMY_FUNCTION(int, int, int, int, int) {} }

void defineMultiConvolutionFunctions()
{
    using namespace python;

    def("gaussianSmoothing3D",
        registerConverters(&pythonGaussianSmoothingND<float,4>),               // also multiband
        (arg("volume"), arg("sigma"), arg("out")=python::object()) );

    def("gaussianGradient3D",
    	registerConverters(&pythonGaussianGradientND<float,4>),               // also multiband?
    	(arg("volume"), arg("sigma"), arg("out")=python::object()) );

    def("symmetricGradient3D",
        registerConverters(&pythonSymmetricGradientND<float,4>),                // also multiband?
        (arg("volume"), arg("out")=python::object()) );

    def("convolveOneDimension3D",
    	registerConverters(&pythonConvolveOneDimensionND<float,4>),				// also multiband
    	(arg("volume"), arg("dim"), arg("kernel1D"), arg("out")=python::object()) );

    def("separableConvolve3D",
        // also multiband
        registerConverters(&pythonSeparableConvolveND<float,4>),
        (arg("volume"), arg("kernel1D"), arg("out")=python::object()), "apply one 1D kernel to all dimensions" );

    def("separableConvolve3D",
        // also multiband
        registerConverters(&pythonSeparableConvolve_3_3D<float>),
        (arg("volume"), arg("kernel1D"), arg("kernel1D"), arg("kernel1D"), arg("out")=python::object()), "apply one kernel for every dimension of the data. The first kernel in this sequence is applied to the innermost dimension (e.g. the x-dimension of an image), while the last is applied to the outermost dimension (e.g. the z-dimension in a 3D image)."
       );
    }

} // namespace vigra
