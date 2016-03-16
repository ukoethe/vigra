/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Ullrich Koethe and Hans Meine                */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycore_PyArray_API
#define NO_IMPORT_ARRAY

#include <Python.h>
#include <vigra/matrix.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <boost/python.hpp>
#include <boost/python/to_python_converter.hpp>
#include <numpy/arrayscalars.h>

namespace python = boost::python;

namespace vigra {

#define VIGRA_NUMPY_TYPECHECKER(type) \
    if(python::object((python::detail::new_reference)PyArray_TypeObjectFromType(type)).ptr() == obj) \
        return obj;

#define VIGRA_NUMPY_TYPECONVERTER(type) \
    if(python::object((python::detail::new_reference)PyArray_TypeObjectFromType(type)).ptr() == obj) \
        typeID = type;

struct NumpyTypenumConverter
{
    NumpyTypenumConverter()
    {
        python::converter::registry::insert(&convertible, &construct,
                                            python::type_id<NPY_TYPES>());
        python::to_python_converter<NPY_TYPES, NumpyTypenumConverter>();
    }

    static void* convertible(PyObject* obj)
    {
        // FIXME: there should be a more elegant way to program this...
        if(obj == 0)
            return 0;
        if(obj->ob_type == &PyArrayDescr_Type)
            return obj;
        if(!PyType_Check(obj))
            return 0;
        VIGRA_NUMPY_TYPECHECKER(NPY_BOOL)
        VIGRA_NUMPY_TYPECHECKER(NPY_INT8)
        VIGRA_NUMPY_TYPECHECKER(NPY_UINT8)
        VIGRA_NUMPY_TYPECHECKER(NPY_INT16)
        VIGRA_NUMPY_TYPECHECKER(NPY_UINT16)
        VIGRA_NUMPY_TYPECHECKER(NPY_INT32)
        VIGRA_NUMPY_TYPECHECKER(NPY_UINT32)
        VIGRA_NUMPY_TYPECHECKER(NPY_INT)
        VIGRA_NUMPY_TYPECHECKER(NPY_UINT)
        VIGRA_NUMPY_TYPECHECKER(NPY_INT64)
        VIGRA_NUMPY_TYPECHECKER(NPY_UINT64)
        VIGRA_NUMPY_TYPECHECKER(NPY_FLOAT32)
        VIGRA_NUMPY_TYPECHECKER(NPY_FLOAT64)
        VIGRA_NUMPY_TYPECHECKER(NPY_LONGDOUBLE)
        VIGRA_NUMPY_TYPECHECKER(NPY_CFLOAT)
        VIGRA_NUMPY_TYPECHECKER(NPY_CDOUBLE)
        VIGRA_NUMPY_TYPECHECKER(NPY_CLONGDOUBLE)
        return 0;
    }

    // from Python
    static void construct(PyObject* obj,
        python::converter::rvalue_from_python_stage1_data* data)
    {
        void* const storage =
            ((python::converter::rvalue_from_python_storage<NumpyAnyArray>* ) data)->storage.bytes;

        // FIXME: there should be a more elegant way to program this...
        int typeID = -1;
        if(obj->ob_type == &PyArrayDescr_Type)
            typeID = (NPY_TYPES)((PyArray_Descr*)obj)->type_num;
        VIGRA_NUMPY_TYPECONVERTER(NPY_BOOL)
        VIGRA_NUMPY_TYPECONVERTER(NPY_INT8)
        VIGRA_NUMPY_TYPECONVERTER(NPY_UINT8)
        VIGRA_NUMPY_TYPECONVERTER(NPY_INT16)
        VIGRA_NUMPY_TYPECONVERTER(NPY_UINT16)
        VIGRA_NUMPY_TYPECONVERTER(NPY_INT32)
        VIGRA_NUMPY_TYPECONVERTER(NPY_UINT32)
        VIGRA_NUMPY_TYPECONVERTER(NPY_INT)
        VIGRA_NUMPY_TYPECONVERTER(NPY_UINT)
        VIGRA_NUMPY_TYPECONVERTER(NPY_INT64)
        VIGRA_NUMPY_TYPECONVERTER(NPY_UINT64)
        VIGRA_NUMPY_TYPECONVERTER(NPY_FLOAT32)
        VIGRA_NUMPY_TYPECONVERTER(NPY_FLOAT64)
        VIGRA_NUMPY_TYPECONVERTER(NPY_LONGDOUBLE)
        VIGRA_NUMPY_TYPECONVERTER(NPY_CFLOAT)
        VIGRA_NUMPY_TYPECONVERTER(NPY_CDOUBLE)
        VIGRA_NUMPY_TYPECONVERTER(NPY_CLONGDOUBLE)

        new (storage) NPY_TYPES((NPY_TYPES)typeID);

        data->convertible = storage;
    }

    // to Python
    static PyObject* convert(NPY_TYPES typeID)
    {
        return PyArray_TypeObjectFromType(typeID);
    }
};

#undef VIGRA_NUMPY_TYPECHECKER
#undef VIGRA_NUMPY_TYPECONVERTER

struct NumpyAnyArrayConverter
{
    NumpyAnyArrayConverter()
    {
        python::converter::registry::insert(&convertible, &construct,
                                            python::type_id<NumpyAnyArray>());
        python::to_python_converter<NumpyAnyArray, NumpyAnyArrayConverter>();
    }

    static void* convertible(PyObject* obj)
    {
        return obj && (obj == Py_None || PyArray_Check(obj))
                 ? obj
                 : 0;
    }

    // from Python
    static void construct(PyObject* obj,
        python::converter::rvalue_from_python_stage1_data* data)
    {
        void* const storage =
            ((python::converter::rvalue_from_python_storage<NumpyAnyArray>* ) data)->storage.bytes;

        if(obj == Py_None)
            obj = 0;

        new (storage) NumpyAnyArray(obj);

        data->convertible = storage;
    }

    static PyObject* convert(NumpyAnyArray const& a)
    {
        return returnNumpyArray(a);
    }
};

namespace detail {

template <int N, class T>
struct MultiArrayShapeConverterTraits
{
    typedef TinyVector<T, N> ShapeType;

    static void construct(void* const storage, PyObject * obj)
    {
        ShapeType * shape = new (storage) ShapeType();
        for(int i=0; i<PySequence_Length(obj); ++i)
            (*shape)[i] = python::extract<T>(PySequence_ITEM(obj, i));
    }
};

template <class T>
struct MultiArrayShapeConverterTraits<0, T>
{
    typedef ArrayVector<T> ShapeType;

    static void construct(void* const storage, PyObject * obj)
    {
        int len = (obj == Py_None)
                           ? 0
                           : PySequence_Length(obj);
        ShapeType * shape = new (storage) ShapeType(len);
        for(int i=0; i<len; ++i)
            (*shape)[i] = python::extract<T>(PySequence_ITEM(obj, i));
    }
};

} // namespace detail


template <int M, class T>
struct MultiArrayShapeConverter
{

    typedef typename detail::MultiArrayShapeConverterTraits<M, T>::ShapeType ShapeType;

    MultiArrayShapeConverter()
    {
        python::converter::registry::insert(&convertible, &construct,
                                            python::type_id<ShapeType>());
        python::to_python_converter<ShapeType, MultiArrayShapeConverter>();
    }

    static void* convertible(PyObject* obj)
    {
        if(obj == 0)
            return 0;
        if(M == 0 && obj == Py_None)
            return obj;
        if(!PySequence_Check(obj) || (M != 0 && PySequence_Length(obj) != M))
            return 0;
        for(int i=0; i<PySequence_Length(obj); ++i)
            if(!PyNumber_Check(PySequence_ITEM(obj, i)))
                return 0;
        return obj;
    }

    // from Python
    static void construct(PyObject* obj,
        python::converter::rvalue_from_python_stage1_data* data)
    {
        void* const storage =
            ((python::converter::rvalue_from_python_storage<ShapeType>* ) data)->storage.bytes;

        detail::MultiArrayShapeConverterTraits<M, T>::construct(storage, obj);
        data->convertible = storage;
    }

    // to Python
    static PyObject* convert(ShapeType const& shape)
    {
        return shapeToPythonTuple(shape).release();
    }
};

python_ptr point2DToPythonTuple(Point2D const & point)
{
    python_ptr tuple(PyTuple_New(2), python_ptr::keep_count);
    pythonToCppException(tuple);
    PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), 0 ,pythonFromData(point.x).release());
    PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), 1 ,pythonFromData(point.y).release());
    return tuple;
}

struct Point2DConverter
{
    Point2DConverter()
    {
        python::converter::registry::insert(&convertible, &construct,
                                            python::type_id<Point2D>());
        python::to_python_converter<Point2D, Point2DConverter>();
    }

    static void* convertible(PyObject* obj)
    {
        if(obj == 0 || !PySequence_Check(obj) || (PySequence_Length(obj) !=2))
            return 0;
        if(!PyNumber_Check(PySequence_Fast_GET_ITEM(obj,0)))
            return 0;
        if(!PyNumber_Check(PySequence_Fast_GET_ITEM(obj,0)))
            return 0;
        return obj;
    }

    //from python
    static void construct(PyObject* obj, python::converter::rvalue_from_python_stage1_data* data)
    {
        void* const storage =
            ((python::converter::rvalue_from_python_storage<Point2D>*) data)->storage.bytes;
        new (storage) Point2D(python::extract<int>(PySequence_Fast_GET_ITEM(obj,0)),
                              python::extract<int>(PySequence_Fast_GET_ITEM(obj,1)));
        data->convertible = storage;
    }

    //to python
    static PyObject* convert(Point2D const& p)
    {
        return point2DToPythonTuple(p).release();
    }

};

void registerNumpyPoint2DConverter()
{
    Point2DConverter();
}

template <class T>
void registerNumpyShapeConvertersOneType()
{
    MultiArrayShapeConverter<0, T>();
    MultiArrayShapeConverter<1, T>();
    MultiArrayShapeConverter<2, T>();
    MultiArrayShapeConverter<3, T>();
    MultiArrayShapeConverter<4, T>();
    MultiArrayShapeConverter<5, T>();
    MultiArrayShapeConverter<6, T>();
    MultiArrayShapeConverter<7, T>();
    MultiArrayShapeConverter<8, T>();
    MultiArrayShapeConverter<9, T>();
    MultiArrayShapeConverter<10, T>();
}

void registerNumpyShapeConvertersAllTypes()
{
    registerNumpyShapeConvertersOneType<MultiArrayIndex>();
    registerNumpyShapeConvertersOneType<float>();
    registerNumpyShapeConvertersOneType<double>();
    registerNumpyShapeConvertersOneType<short>();
    if(typeid(int) != typeid(MultiArrayIndex))
        registerNumpyShapeConvertersOneType<int>();
    if(typeid(npy_intp) != typeid(MultiArrayIndex))
        MultiArrayShapeConverter<0, npy_intp>();
}

#if 0 // FIXME: reimplement to replace the Python versions for consistence?
PyObject *
constructNumpyArrayFromShape(python::object type, ArrayVector<npy_intp> const & shape,
                       unsigned int spatialDimensions, unsigned int channels,
                       NPY_TYPES typeCode, std::string order, bool init)
{
    PyObject * obj = type.ptr();
    vigra_precondition(obj && PyType_Check(obj) && PyType_IsSubtype((PyTypeObject *)obj, &PyArray_Type),
           "constructNumpyArray(type, ...): first argument was not an array type.");
    return detail::constructNumpyArrayImpl((PyTypeObject *)obj, shape, spatialDimensions, channels, typeCode, order, init).release();
}

PyObject *
constructNumpyArrayFromArray(python::object type, NumpyAnyArray array,
                       unsigned int spatialDimensions, unsigned int channels,
                       NPY_TYPES typeCode, std::string order, bool init)
{
    PyObject * obj = type.ptr();
    vigra_precondition(obj && PyType_Check(obj) && PyType_IsSubtype((PyTypeObject *)obj, &PyArray_Type),
           "constructNumpyArray(type, ...): first argument was not an array type.");
    PyObject * res = detail::constructNumpyArrayImpl((PyTypeObject *)obj, array.shape(), spatialDimensions, channels,
                                                     typeCode, order, false, array.strideOrdering()).release();
    if(init)
    {
        NumpyAnyArray lhs(res);
        lhs = array;
    }
    return res;
}
#endif

PyObject *
constructArrayFromAxistags(python::object type, ArrayVector<npy_intp> const & shape,
                           NPY_TYPES typeCode, AxisTags const & axistags, bool init)
{
    PyAxisTags pyaxistags(python_ptr(python::object(axistags).ptr()));

    ArrayVector<npy_intp> norm_shape(shape);
    if(pyaxistags.size() > 0)
    {
        ArrayVector<npy_intp> permutation(pyaxistags.permutationToNormalOrder());
        applyPermutation(permutation.begin(), permutation.end(), shape.begin(), norm_shape.begin());
    }

    TaggedShape tagged_shape(norm_shape, pyaxistags);

	// FIXME: check that type is an array class?
	return constructArray(tagged_shape, typeCode, init, python_ptr(type.ptr()));
}

template <class T>
struct MatrixConverter
{
    typedef linalg::Matrix<T> ArrayType;

    MatrixConverter();

    // to Python
    static PyObject* convert(ArrayType const& a)
    {
        return returnNumpyArray(NumpyArray<2, T>(a));
    }
};

template <class T>
MatrixConverter<T>::MatrixConverter()
{
    using namespace boost::python;

    converter::registration const * reg = converter::registry::query(type_id<ArrayType>());

    // register the to_python_converter only once
    if(!reg || !reg->rvalue_chain)
    {
        to_python_converter<ArrayType, MatrixConverter>();
    }
}

template <typename ScalarType>
struct NumpyScalarConverter
{
    NumpyScalarConverter()
    {
        using namespace boost::python;
        converter::registry::push_back( &convertible, &construct, type_id<ScalarType>());
    }

    // Determine if obj_ptr is a supported numpy.number
    static void* convertible(PyObject* obj_ptr)
    {
        if (PyArray_IsScalar(obj_ptr, Float32) ||
            PyArray_IsScalar(obj_ptr, Float64) ||
            PyArray_IsScalar(obj_ptr, Int8)    ||
            PyArray_IsScalar(obj_ptr, Int16)   ||
            PyArray_IsScalar(obj_ptr, Int32)   ||
            PyArray_IsScalar(obj_ptr, Int64)   ||
            PyArray_IsScalar(obj_ptr, UInt8)   ||
            PyArray_IsScalar(obj_ptr, UInt16)  ||
            PyArray_IsScalar(obj_ptr, UInt32)  ||
            PyArray_IsScalar(obj_ptr, UInt64))
        {
            return obj_ptr;
        }
        return 0;
    }

    static void construct( PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;

        // Grab pointer to memory into which to construct the C++ scalar
        void* storage = ((converter::rvalue_from_python_storage<ScalarType>*) data)->storage.bytes;

        // in-place construct the new scalar value
        ScalarType * scalar = new (storage) ScalarType;

        if (PyArray_IsScalar(obj_ptr, Float32))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Float32);
        else if (PyArray_IsScalar(obj_ptr, Float64))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Float64);
        else if (PyArray_IsScalar(obj_ptr, Int8))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Int8);
        else if (PyArray_IsScalar(obj_ptr, Int16))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Int16);
        else if (PyArray_IsScalar(obj_ptr, Int32))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Int32);
        else if (PyArray_IsScalar(obj_ptr, Int64))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, Int64);
        else if (PyArray_IsScalar(obj_ptr, UInt8))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, UInt8);
        else if (PyArray_IsScalar(obj_ptr, UInt16))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, UInt16);
        else if (PyArray_IsScalar(obj_ptr, UInt32))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, UInt32);
        else if (PyArray_IsScalar(obj_ptr, UInt64))
            (*scalar) = PyArrayScalar_VAL(obj_ptr, UInt64);

        // Stash the memory chunk pointer for later use by boost.python
        data->convertible = storage;
    }
};


void registerNumpyArrayConverters()
{
    NumpyTypenumConverter();
    registerNumpyShapeConvertersAllTypes();
    registerNumpyPoint2DConverter();
    NumpyAnyArrayConverter();
    MatrixConverter<float>();
    MatrixConverter<double>();

    NumpyScalarConverter<signed char>();
    NumpyScalarConverter<short>();
    NumpyScalarConverter<int>();
    NumpyScalarConverter<long>();
    NumpyScalarConverter<long long>();
    NumpyScalarConverter<unsigned char>();
    NumpyScalarConverter<unsigned short>();
    NumpyScalarConverter<unsigned int>();
    NumpyScalarConverter<unsigned long>();
    NumpyScalarConverter<unsigned long long>();
    NumpyScalarConverter<float>();
    NumpyScalarConverter<double>();

    python::docstring_options doc_options(true, true, false);

    doc_options.disable_all();
    python::def("constructArrayFromAxistags", &constructArrayFromAxistags);
    // python::def("constructNumpyArray", &constructNumpyArrayFromArray);
}

} // namespace vigra

