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

namespace python = boost::python;

namespace vigra
{

#if 0
template < class VoxelType, int pyArrayTypeConstant, class Tag, int ndim >
PyObject* gaussianSmoothing(MultiArrayView< ndim, VoxelType, Tag > const
    & volume, double sigma)
{
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< VoxelType >::isScalar());
    MultiArrayView< ndim, VoxelType, UnstridedArrayTag >* arrayView =
        createView< VoxelType, ndim >(array,
        typename NumericTraits< VoxelType >::isScalar());
    gaussianSmoothMultiArray(srcMultiArrayRange(volume),
        destMultiArray(*arrayView), sigma);
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class VoxelType, int pyArrayTypeConstant, class Tag, int ndim >
PyObject* gaussianGradient(MultiArrayView< ndim, VoxelType, Tag > const
    & volume, double sigma)
{
    vigra_precondition(typename NumericTraits< VoxelType >::isScalar().asBool,
        "Only single band arrays are supported.");
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< TinyVector< VoxelType, 3 > >::isScalar());
    MultiArrayView< ndim,  TinyVector< VoxelType, 3 >, UnstridedArrayTag >*
        arrayView = createView< TinyVector< VoxelType, 3 >, ndim >(array,
        typename NumericTraits< TinyVector< VoxelType, 3 > >::isScalar());
    gaussianGradientMultiArray(srcMultiArrayRange(volume),
        destMultiArray(*arrayView), sigma);
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class VoxelType, int pyArrayTypeConstant, class Tag, int ndim >
PyObject* symmetricGradient(MultiArrayView< ndim, VoxelType, Tag > const
    & volume)
{
    vigra_precondition(typename NumericTraits< VoxelType >::isScalar().asBool,
        "Only single band arrays are supported.");
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< TinyVector< VoxelType, 3 > >::isScalar());
    MultiArrayView< ndim,  TinyVector< VoxelType, 3 >, UnstridedArrayTag >*
        arrayView = createView< TinyVector< VoxelType, 3 >, ndim >(array,
        typename NumericTraits< TinyVector< VoxelType, 3 > >::isScalar());
    symmetricGradientMultiArray(srcMultiArrayRange(volume),
        destMultiArray(*arrayView));
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class VoxelType, int pyArrayTypeConstant, class Tag, int ndim >
PyObject* convolveOneDimension(MultiArrayView< ndim, VoxelType, Tag > const
    & volume, int dim, Kernel const & kernel)
{
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< VoxelType >::isScalar());
    MultiArrayView< ndim, VoxelType, UnstridedArrayTag >* arrayView =
        createView< VoxelType, ndim >(array,
        typename NumericTraits< VoxelType >::isScalar());
    convolveMultiArrayOneDimension(srcMultiArrayRange(volume),
        destMultiArray(*arrayView), dim, kernel);
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

template < class VoxelType, int pyArrayTypeConstant, class Tag, int ndim >
PyObject* separableConvolve(MultiArrayView< ndim, VoxelType, Tag > const
    & volume, Kernel const & kernel)
{
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< VoxelType >::isScalar());
    MultiArrayView< ndim, VoxelType, UnstridedArrayTag >* arrayView =
        createView< VoxelType, ndim >(array,
        typename NumericTraits< VoxelType >::isScalar());
    separableConvolveMultiArray(srcMultiArrayRange(volume),
        destMultiArray(*arrayView), kernel);
    delete arrayView;
    return PyArray_Return((PyArrayObject*) array);
}

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

template < class VoxelType, int pyArrayTypeConstant, class Tag >
PyObject* separableConvolve3(MultiArrayView< 3, VoxelType, Tag > const
    & volume, Kernel const & kernel0,
    Kernel const & kernel1, Kernel const & kernel2)
{
    PyObject* array = createNumpyArray(volume.shape(), pyArrayTypeConstant,
        typename NumericTraits< VoxelType >::isScalar());
    MultiArrayView< 3, VoxelType, UnstridedArrayTag >* arrayView =
        createView< VoxelType, 3 >(array,
        typename NumericTraits< VoxelType >::isScalar());
    ArrayVector< Kernel > kernels(3);
    kernels.push_back(kernel0);
    kernels.push_back(kernel1);
    kernels.push_back(kernel2);
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
        &DUMMY_FUNCTION,               // also multiband
        (arg("volume"), arg("sigma")));

    def("gaussianGradient3D",
        &DUMMY_FUNCTION,               // also multiband?
        (arg("volume"), arg("sigma")));

    def("symmetricGradient3D",
        &DUMMY_FUNCTION,               // also multiband?
        (arg("volume")));

    def("convolveOneDimension3D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("volume"), arg("dim"), arg("kernel1D")));

    def("separableConvolve3D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("volume"), arg("kernel1D")));

    def("separableConvolve3D",
        &DUMMY_FUNCTION,               // also multiband
        (arg("volume"), arg("kernel1Ddim0"), arg("kernel1Ddim1"),
        arg("kernel1Ddim2")));
}

} // namespace vigra
