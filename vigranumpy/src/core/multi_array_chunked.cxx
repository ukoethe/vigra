/************************************************************************/
/*                                                                      */
/*               Copyright 2013-2014 by Ullrich Koethe                  */
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

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/axistags.hxx>
#include <vigra/multi_array_chunked.hxx>
#ifdef HasHDF5
#include <vigra/multi_array_chunked_hdf5.hxx>
#endif
#include <vigra/compression.hxx>
#include <vigra/python_utility.hxx>
#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include <sstream>

namespace python = boost::python;

namespace vigra {

#ifdef HasHDF5
// workaround for compiler bug in VS 2015 (compiler fails to match the template
// function get_pointer() at line 20 of boost/get_pointer.hpp)
inline ChunkedArrayHDF5<1, float, class std::allocator<float> > const volatile *
get_pointer(ChunkedArrayHDF5<1, float, class std::allocator<float> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<1, unsigned long, class std::allocator<unsigned long> > const volatile *
get_pointer(ChunkedArrayHDF5<1, unsigned long, class std::allocator<unsigned long> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<1, unsigned char, class std::allocator<unsigned char> > const volatile *
get_pointer(ChunkedArrayHDF5<1, unsigned char, class std::allocator<unsigned char> > const volatile * p) { return p; }


inline ChunkedArrayHDF5<2, float, class std::allocator<float> > const volatile *
get_pointer(ChunkedArrayHDF5<2, float, class std::allocator<float> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<2, unsigned char, class std::allocator<unsigned char> > const volatile *
get_pointer(ChunkedArrayHDF5<2, unsigned char, class std::allocator<unsigned char> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<2, unsigned long, class std::allocator<unsigned long> > const volatile *
get_pointer(ChunkedArrayHDF5<2, unsigned long, class std::allocator<unsigned long> > const volatile * p) { return p; }


inline ChunkedArrayHDF5<3, float, class std::allocator<float> > const volatile *
get_pointer(ChunkedArrayHDF5<3, float, class std::allocator<float> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<3, unsigned char, class std::allocator<unsigned char> > const volatile *
get_pointer(ChunkedArrayHDF5<3, unsigned char, class std::allocator<unsigned char> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<3, unsigned long, class std::allocator<unsigned long> > const volatile *
get_pointer(ChunkedArrayHDF5<3, unsigned long, class std::allocator<unsigned long> > const volatile * p) { return p; }


inline ChunkedArrayHDF5<4, float, class std::allocator<float> > const volatile *
get_pointer(ChunkedArrayHDF5<4, float, class std::allocator<float> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<4, unsigned char, class std::allocator<unsigned char> > const volatile *
get_pointer(ChunkedArrayHDF5<4, unsigned char, class std::allocator<unsigned char> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<4, unsigned long, class std::allocator<unsigned long> > const volatile *
get_pointer(ChunkedArrayHDF5<4, unsigned long, class std::allocator<unsigned long> > const volatile * p) { return p; }


inline ChunkedArrayHDF5<5, float, class std::allocator<float> > const volatile *
get_pointer(ChunkedArrayHDF5<5, float, class std::allocator<float> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<5, unsigned char, class std::allocator<unsigned char> > const volatile *
get_pointer(ChunkedArrayHDF5<5, unsigned char, class std::allocator<unsigned char> > const volatile * p) { return p; }

inline ChunkedArrayHDF5<5, unsigned long, class std::allocator<unsigned long> > const volatile *
get_pointer(ChunkedArrayHDF5<5, unsigned long, class std::allocator<unsigned long> > const volatile * p) { return p; }
#endif // HasHDF5

inline ChunkedArray<2, float> const volatile *
get_pointer(ChunkedArray<2, float> const volatile * p) { return p; }

inline ChunkedArray<2, unsigned long> const volatile *
get_pointer(ChunkedArray<2, unsigned long> const volatile * p) { return p; }

inline ChunkedArray<2, unsigned char> const volatile *
get_pointer(ChunkedArray<2, unsigned char> const volatile * p) { return p; }


inline ChunkedArray<3, float> const volatile *
get_pointer(ChunkedArray<3, float> const volatile * p) { return p; }

inline ChunkedArray<3, unsigned long> const volatile *
get_pointer(ChunkedArray<3, unsigned long> const volatile * p) { return p; }

inline ChunkedArray<3, unsigned char> const volatile *
get_pointer(ChunkedArray<3, unsigned char> const volatile * p) { return p; }


inline ChunkedArray<4, float> const volatile *
get_pointer(ChunkedArray<4, float> const volatile * p) { return p; }

inline ChunkedArray<4, unsigned long> const volatile *
get_pointer(ChunkedArray<4, unsigned long> const volatile * p) { return p; }

inline ChunkedArray<4, unsigned char> const volatile *
get_pointer(ChunkedArray<4, unsigned char> const volatile * p) { return p; }


inline ChunkedArray<5, float> const volatile *
get_pointer(ChunkedArray<5, float> const volatile * p) { return p; }

inline ChunkedArray<5, unsigned long> const volatile *
get_pointer(ChunkedArray<5, unsigned long> const volatile * p) { return p; }

inline ChunkedArray<5, unsigned char> const volatile *
get_pointer(ChunkedArray<5, unsigned char> const volatile * p) { return p; }

template <unsigned int N, class T>
TinyVector<MultiArrayIndex, N>
ChunkedArray_shape(ChunkedArray<N, T> const & array)
{
    return array.shape();
}

template <unsigned int N, class T>
TinyVector<MultiArrayIndex, N>
ChunkedArray_chunkShape(ChunkedArray<N, T> const & array)
{
    return array.chunkShape();
}

template <unsigned int N, class T>
TinyVector<MultiArrayIndex, N>
ChunkedArray_chunkArrayShape(ChunkedArray<N, T> const & array)
{
    return array.chunkArrayShape();
}

template <unsigned int N, class T>
std::string ChunkedArray_repr(ChunkedArray<N, T> const & array)
{
    std::stringstream s;
    s << array.backend() << "( shape=" << array.shape() <<
           ", dtype=" << NumpyArrayValuetypeTraits<T>::typeName() << ")";
    return s.str();
}

template <unsigned int N, class T>
std::string ChunkedArray_str(ChunkedArray<N, T> const & array)
{
    return ChunkedArray_repr(array);
}

template <unsigned int N, class T>
PyObject * ChunkedArray_dtype(ChunkedArray<N, T> const &)
{
    return NumpyArrayValuetypeTraits<T>::typeObject();
}

template <unsigned int N, class T>
unsigned int ChunkedArray_ndim(ChunkedArray<N, T> const &)
{
    return N;
}

template <unsigned int N, class T>
NumpyAnyArray
ChunkedArray_checkoutSubarray(python::object array,
                              TinyVector<MultiArrayIndex, N> const & start,
                              TinyVector<MultiArrayIndex, N> const & stop,
                              NumpyArray<N, T> res = NumpyArray<N, T>())
{
    ChunkedArray<N, T> const & self = python::extract<ChunkedArray<N, T> const &>(array)();

    python_ptr pytags;
    if(PyObject_HasAttrString(array.ptr(), "axistags"))
    {
        pytags = python_ptr(PyObject_GetAttrString(array.ptr(), "axistags"), python_ptr::keep_count);
    }
    PyAxisTags tags(pytags, true);
    TaggedShape shape(stop-start, tags);
    res.reshapeIfEmpty(shape,
        "ChunkedArray::checkoutSubarray(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        self.checkoutSubarray(start, res);
    }
    return res;
}

template <unsigned int N, class T>
void
ChunkedArray_commitSubarray(ChunkedArray<N, T> & self,
                            TinyVector<MultiArrayIndex, N> const & start,
                            NumpyArray<N, T> array)
{
    PyAllowThreads _pythread;
    self.commitSubarray(start, array);
}

template <class Shape>
python::object
bindNumpyArray(NumpyAnyArray self, Shape const & stop)
{
    if(stop == Shape())
        return python::object(self);

    python_ptr func(pythonFromData("__getitem__"));
    pythonToCppException(func);
    python_ptr index(PyTuple_New(stop.size()), python_ptr::keep_count);
    pythonToCppException(index);
    for(unsigned int k=0; k<stop.size(); ++k)
    {
        PyObject * item = stop[k] == 0
                            ? pythonFromData(0)
                            : PySlice_New(0,0,0);
        pythonToCppException(item);
        PyTuple_SET_ITEM((PyTupleObject *)index.ptr(), k, item);
    }

    return python::object(python::detail::new_non_null_reference(PyObject_CallMethodObjArgs(self.pyObject(), func.ptr(), index.ptr(), NULL)));
}

template <unsigned int N, class T>
python::object
ChunkedArray_getitem(python::object array, python::object index)
{
    typedef typename ChunkedArray<N, T>::shape_type Shape;

    ChunkedArray<N, T> const & self = python::extract<ChunkedArray<N, T> const &>(array)();
    Shape start, stop;
    numpyParseSlicing(self.shape(), index.ptr(), start, stop);
    if(start == stop)
    {
        // return a single point
        return python::object(self.getItem(start));
    }
    else if(allLessEqual(start, stop))
    {
        // return a slice
        NumpyAnyArray subarray = ChunkedArray_checkoutSubarray<N,T>(array, start, max(start + Shape(1), stop));
        return python::object(subarray.getitem(Shape(), stop-start));
    }
    else
    {
        vigra_precondition(false,
            "ChunkedArray.__getitem__(): index out of bounds.");
        return python::object();
    }
}

template <unsigned int N, class T>
void
ChunkedArray_setitem(ChunkedArray<N, T> & self, python::object index, T value)
{
    typedef typename ChunkedArray<N, T>::shape_type Shape;
    Shape start, stop;
    numpyParseSlicing(self.shape(), index.ptr(), start, stop);
    if(start == stop)
    {
        self.setItem(start, value);
    }
    else
    {
        PyAllowThreads _pythread;
        stop = max(start + Shape(1), stop);
        typename ChunkedArray<N, T>::iterator i(self.begin().restrictToSubarray(start, stop)),
                                              end(i.getEndIterator());
        for(; i != end; ++i)
            *i = value;
    }
}

template <unsigned int N, class T>
void
ChunkedArray_setitem2(ChunkedArray<N, T> & self, python::object index, NumpyArray<N, T> array)
{
    typedef typename ChunkedArray<N, T>::shape_type Shape;
    Shape start, stop;
    numpyParseSlicing(self.shape(), index.ptr(), start, stop);
    stop = max(start + Shape(1), stop);
    vigra_precondition(array.shape() == stop - start,
        "ChunkedArray.__setitem__(): shape mismatch");

    PyAllowThreads _pythread;
    self.commitSubarray(start, array);
}

python::object defaultDtype()
{
    PyObject * dtype = NumpyArrayValuetypeTraits<npy_float32>::typeObject();
    return python::object(python::detail::new_reference(dtype));
}

int numpyScalarTypeNumber(python::object obj)
{
    PyArray_Descr* dtype;
    if(!PyArray_DescrConverter(obj.ptr(), &dtype))
        return NPY_NOTYPE;
    int typeNum = dtype->type_num;
    Py_DECREF(dtype);
    return typeNum;
}

template <class Array>
PyObject *
ptr_to_python(Array * array, python::object axistags)
{
    python_ptr py_array(python::to_python_indirect<Array*,
                                      python::detail::make_owning_holder>()(array),
                        python_ptr::new_nonzero_reference);
    if(axistags != python::object())
    {
        AxisTags at;
#if PY_MAJOR_VERSION < 3
        if(PyString_Check(axistags.ptr()))
#else
		if(PyUnicode_Check(axistags.ptr()))
#endif
			at = AxisTags(python::extract<std::string>(axistags)());
        else
            at = AxisTags(python::extract<AxisTags const &>(axistags)());
        int N = Array::shape_type::static_size;
        vigra_precondition(at.size() == 0 || at.size() == unsigned(N),
            "ChunkedArray(): axistags have invalid length.");
        if(at.size() == unsigned(N))
        {
            int res = PyObject_SetAttrString(py_array, "axistags", python::object(at).ptr());
            pythonToCppException(res != 0);
        }
    }
    return py_array.release();
}

template <class T, int N>
ChunkedArray<N, T> *
construct_ChunkedArrayFullImpl(TinyVector<MultiArrayIndex, N> const & shape,
                               double fill_value)
{
    return new ChunkedArrayFull<N, T>(shape,
                                      ChunkedArrayOptions().fillValue(fill_value));
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayFull(TinyVector<MultiArrayIndex, N> const & shape,
                           python::object dtype, double fill_value,
                           python::object axistags)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayFullImpl<npy_uint8>(shape, fill_value), axistags);
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayFullImpl<npy_uint32>(shape, fill_value), axistags);
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayFullImpl<npy_float32>(shape, fill_value), axistags);
      default:
        vigra_precondition(false, "ChunkedArrayFull(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> *
construct_ChunkedArrayLazyImpl(TinyVector<MultiArrayIndex, N> const & shape,
                               TinyVector<MultiArrayIndex, N> const & chunk_shape,
                               double fill_value)
{
    return new ChunkedArrayLazy<N, T>(shape, chunk_shape,
                                      ChunkedArrayOptions().fillValue(fill_value));
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayLazy(TinyVector<MultiArrayIndex, N> const & shape,
                           python::object dtype,
                           TinyVector<MultiArrayIndex, N> const & chunk_shape,
                           double fill_value,
                           python::object axistags)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayLazyImpl<npy_uint8>(shape, chunk_shape, fill_value), axistags);
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayLazyImpl<npy_uint32>(shape, chunk_shape, fill_value), axistags);
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayLazyImpl<npy_float32>(shape, chunk_shape, fill_value), axistags);
      default:
        vigra_precondition(false, "ChunkedArrayLazy(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> *
construct_ChunkedArrayCompressedImpl(TinyVector<MultiArrayIndex, N> const & shape,
                                     CompressionMethod method,
                                     TinyVector<MultiArrayIndex, N> const & chunk_shape,
                                     int cache_max,
                                     double fill_value)
{
    return new ChunkedArrayCompressed<N, T>(shape, chunk_shape,
            ChunkedArrayOptions().compression(method).cacheMax(cache_max).fillValue(fill_value));
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayCompressed(TinyVector<MultiArrayIndex, N> const & shape,
                                 CompressionMethod method,
                                 python::object dtype,
                                 TinyVector<MultiArrayIndex, N> const & chunk_shape,
                                 int cache_max,
                                 double fill_value,
                                 python::object axistags)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayCompressedImpl<npy_uint8>(shape, method, chunk_shape,
                             cache_max, fill_value), axistags);
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayCompressedImpl<npy_uint32>(shape, method, chunk_shape,
                             cache_max, fill_value), axistags);
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayCompressedImpl<npy_float32>(shape, method, chunk_shape,
                             cache_max, fill_value), axistags);
      default:
        vigra_precondition(false, "ChunkedArrayCompressed(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> *
construct_ChunkedArrayTmpFileImpl(TinyVector<MultiArrayIndex, N> const & shape,
                                  TinyVector<MultiArrayIndex, N> const & chunk_shape,
                                  int cache_max,
                                  std::string path,
                                  double fill_value)
{
    return new ChunkedArrayTmpFile<N, T>(shape, chunk_shape,
                                ChunkedArrayOptions().cacheMax(cache_max).fillValue(fill_value), path);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayTmpFile(TinyVector<MultiArrayIndex, N> const & shape,
                              python::object dtype,
                              TinyVector<MultiArrayIndex, N> const & chunk_shape,
                              int cache_max,
                              std::string path,
                              double fill_value,
                              python::object axistags)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayTmpFileImpl<npy_uint8>(shape, chunk_shape, cache_max,
                             path, fill_value), axistags);
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayTmpFileImpl<npy_uint32>(shape, chunk_shape, cache_max,
                             path, fill_value), axistags);
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayTmpFileImpl<npy_float32>(shape, chunk_shape, cache_max,
                             path, fill_value), axistags);
      default:
        vigra_precondition(false, "ChunkedArrayTmpFile(): unsupported dtype.");
    }
    return 0;
}

#ifdef HasHDF5

template <class T, int N>
ChunkedArrayHDF5<N, T> *
construct_ChunkedArrayHDF5Impl(HDF5File const & file,
                               std::string datasetName,
                               TinyVector<MultiArrayIndex, N> const & shape,
                               HDF5File::OpenMode mode,
                               CompressionMethod method,
                               TinyVector<MultiArrayIndex, N> const & chunk_shape,
                               int cache_max,
                               double fill_value)
{
    return new ChunkedArrayHDF5<N, T>(file, datasetName, mode,
                                      shape, chunk_shape,
                                      ChunkedArrayOptions().compression(method).cacheMax(cache_max).fillValue(fill_value));
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayHDF5Impl(HDF5File const & file,
                               std::string datasetName,
                               TinyVector<MultiArrayIndex, N> const & shape,
                               python::object dtype,
                               HDF5File::OpenMode mode,
                               CompressionMethod compression,
                               TinyVector<MultiArrayIndex, N> const & chunk_shape,
                               int cache_max,
                               double fill_value,
                               python::object axistags)
{
    int dtype_code = NPY_FLOAT32;

    if(dtype != python::object())
    {
        dtype_code = numpyScalarTypeNumber(dtype);
    }
    else if(file.existsDataset(datasetName))
    {
        std::string type = file.getDatasetType(datasetName);
        if(type == "UINT8")
            dtype_code = NPY_UINT8;
        else if(type == "UINT32")
            dtype_code = NPY_UINT32;
    }
    switch(dtype_code)
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayHDF5Impl<npy_uint8>(file, datasetName, shape,
                             mode, compression, chunk_shape, cache_max, fill_value), axistags);
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayHDF5Impl<npy_uint32>(file, datasetName, shape,
                             mode, compression, chunk_shape, cache_max, fill_value), axistags);
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayHDF5Impl<npy_float32>(file, datasetName, shape,
                             mode, compression, chunk_shape, cache_max, fill_value), axistags);
      default:
        vigra_precondition(false, "ChunkedArrayHDF5(): unsupported dtype.");
    }
    return 0;
}

PyObject *
construct_ChunkedArrayHDF5Impl(HDF5File const & file,
                               std::string datasetName,
                               python::object py_shape,
                               python::object dtype,
                               HDF5File::OpenMode mode,
                               CompressionMethod compression,
                               python::object py_chunk_shape,
                               int cache_max,
                               double fill_value,
                               python::object axistags)
{
    int ndim = 0;
    bool has_shape = PySequence_Check(py_shape.ptr());
    bool use_existing_dataset = file.existsDataset(datasetName) &&
                                mode != HDF5File::New;
    if(use_existing_dataset)
    {
        ndim = file.getDatasetDimensions(datasetName);
        vigra_precondition(!has_shape || ndim == python::len(py_shape),
            "ChunkedArrayHDF5(): dimension mismatch between dataset and requested shape.");
    }
    else
    {
        vigra_precondition(has_shape,
            "ChunkedArrayHDF5(): cannot create dataset because no shape is given.");
        ndim = python::len(py_shape);
    }

    bool has_chunk_shape = false;
    if(PySequence_Check(py_chunk_shape.ptr()))
    {
        vigra_precondition(python::len(py_chunk_shape) == ndim,
            "ChunkedArrayHDF5(): chunk_shape has wrong dimension.");
        has_chunk_shape = true;
    }

    switch(ndim)
    {
      case 1:
      {
        typedef Shape1 shape_type;

        shape_type shape = has_shape
                                ? python::extract<shape_type>(py_shape)()
                                : shape_type(),
                   chunk_shape = has_chunk_shape
                                ? python::extract<shape_type>(py_chunk_shape)()
                                : shape_type();
        return construct_ChunkedArrayHDF5Impl<1>(file, datasetName, shape, dtype,
                             mode, compression, chunk_shape, cache_max, fill_value, axistags);
      }
      case 2:
      {
        typedef Shape2 shape_type;

        shape_type shape = has_shape
                                ? python::extract<shape_type>(py_shape)()
                                : shape_type(),
                   chunk_shape = has_chunk_shape
                                ? python::extract<shape_type>(py_chunk_shape)()
                                : shape_type();
        return construct_ChunkedArrayHDF5Impl<2>(file, datasetName, shape, dtype,
                             mode, compression, chunk_shape, cache_max, fill_value, axistags);
      }
      case 3:
      {
        typedef Shape3 shape_type;

        shape_type shape = has_shape
                                ? python::extract<shape_type>(py_shape)()
                                : shape_type(),
                   chunk_shape = has_chunk_shape
                                ? python::extract<shape_type>(py_chunk_shape)()
                                : shape_type();
        return construct_ChunkedArrayHDF5Impl<3>(file, datasetName, shape, dtype,
                             mode, compression, chunk_shape, cache_max, fill_value, axistags);
      }
      case 4:
      {
        typedef Shape4 shape_type;

        shape_type shape = has_shape
                                ? python::extract<shape_type>(py_shape)()
                                : shape_type(),
                   chunk_shape = has_chunk_shape
                                ? python::extract<shape_type>(py_chunk_shape)()
                                : shape_type();
        return construct_ChunkedArrayHDF5Impl<4>(file, datasetName, shape, dtype,
                             mode, compression, chunk_shape, cache_max, fill_value, axistags);
      }
      case 5:
      {
        typedef Shape5 shape_type;

        shape_type shape = has_shape
                                ? python::extract<shape_type>(py_shape)()
                                : shape_type(),
                   chunk_shape = has_chunk_shape
                                ? python::extract<shape_type>(py_chunk_shape)()
                                : shape_type();
        return construct_ChunkedArrayHDF5Impl<5>(file, datasetName, shape, dtype,
                             mode, compression, chunk_shape, cache_max, fill_value, axistags);
      }
      default:
        vigra_precondition(false, "ChunkedArrayHDF5(): unsupported array dimension (1 <= ndim <= 5 required).");
    }
    return 0;
}

PyObject *
construct_ChunkedArrayHDF5(std::string filename,
                           std::string datasetName,
                           python::object shape,
                           python::object dtype,
                           HDF5File::OpenMode mode,
                           CompressionMethod compression,
                           python::object chunk_shape,
                           int cache_max,
                           double fill_value,
                           python::object axistags)
{
    bool file_exists = isHDF5(filename.c_str());
    if(mode == HDF5File::Default)
    {
        if(!file_exists)
            mode = HDF5File::New;
        else if(HDF5File(filename, HDF5File::ReadOnly).existsDataset(datasetName))
            mode = HDF5File::ReadOnly;
        else
            mode = HDF5File::Replace;
    }
    HDF5File::OpenMode filemode = mode;
    if(mode == HDF5File::Replace)
    {
        mode = HDF5File::New;
        if(file_exists)
            filemode = HDF5File::ReadWrite;
        else
            filemode = HDF5File::New;
    }
    HDF5File file(filename, filemode);
    return construct_ChunkedArrayHDF5Impl(file, datasetName, shape, dtype, mode, compression,
                               chunk_shape, cache_max, fill_value, axistags);
}

PyObject *
construct_ChunkedArrayHDF5id(hid_t file_id,
                             std::string datasetName,
                             python::object shape,
                             python::object dtype,
                             HDF5File::OpenMode mode,
                             CompressionMethod compression,
                             python::object chunk_shape,
                             int cache_max,
                             double fill_value,
                             python::object axistags)
{
    HDF5HandleShared handle(file_id, 0, "");
    HDF5File file(handle);
    return construct_ChunkedArrayHDF5Impl(file, datasetName, shape, dtype, mode, compression,
                               chunk_shape, cache_max, fill_value, axistags);
}

#endif

template <unsigned int N, class T>
void defineChunkedArrayImpl()
{
    using namespace boost::python;

    docstring_options doc_options(true, false, false);

    typedef ChunkedArray<N, T> Array;
    class_<Array, boost::noncopyable>("ChunkedArrayBase",
         "\n"
         "Base class for chunked arrays, can only be created via factory functions\n"
         "like :func:`~vigra.ChunkedArrayCompressed` or :func:`~vigra.ChunkedArrayHDF5`.\n\n",
         no_init)
        .add_property("shape", &ChunkedArray_shape<N, T>,
             "\nshape of the array.\n")
        .add_property("chunk_shape", &ChunkedArray_chunkShape<N, T>,
             "\nshape of (interior) chunks.\n")
        .add_property("chunk_array_shape", &ChunkedArray_chunkArrayShape<N, T>,
             "\nshape of internal array of chunks.\n")
        .add_property("size", &Array::size,
             "\nnumber of elements of the array.\n")
        .add_property("overhead_bytes", &Array::overheadBytes,
             "\nsize of the overhead caused by chunked storage.\n")
        .add_property("data_bytes", (std::size_t (Array::*)() const)&Array::dataBytes,
             "\nsize of the currently allocated part of the data.\n")
        .add_property("overhead_bytes_per_chunk", &Array::overheadBytesPerChunk,
             "\nsize of the overhead caused by chunked storage for a single chunk.\n")
        .add_property("data_bytes_per_chunk", &Array::dataBytesPerChunk,
             "\nsize of the data of a single chunk.\n")
        .add_property("backend", &Array::backend,
             "\nthe backend driver of this array.\n")
        .add_property("read_only", &Array::isReadOnly,
             "\n'True' if array values cannot be changed.\n")
        .add_property("cache_max_size",
             &Array::cacheMaxSize, &Array::setCacheMaxSize,
             "\nget/set the size of the chunk cache.\n")
        .add_property("dtype", &ChunkedArray_dtype<N, T>,
             "\nthe array's value type\n")
        .add_property("ndim", &ChunkedArray_ndim<N, T>,
             "\nthe array's dimension\n")
        .def("__repr__", &ChunkedArray_repr<N, T>)
        .def("__str__", &ChunkedArray_str<N, T>)
        .def("checkoutSubarray",
             registerConverters(&ChunkedArray_checkoutSubarray<N, T>),
             (arg("start"), arg("stop"), arg("out")=python::object()),
             "\n    checkoutSubarray(start, stop, res=None) => array\n\n"
             "Obtain a copy of the subarray in the ROI '[start, stop)'.\n"
             "If 'res' is given, it must have matching shape and will be used\n"
             "to store the data instead of allocating new storage for 'array'.\n\n"
             "The index operator provides a shorthand for this function, e.g.\n"
             "for a 2-dimensional array you can equivalently write::\n\n"
             "    roi = chunked_array.checkoutSubarray((5,10), (12,19))\n"
             "    roi = chunked_array[5:12, 10:19]\n\n"
             "to read the ROI from 'start=(5,10)' to 'stop=(12,19)' (exclusive).\n"
             "Note that 'roi' is a copy, so overwriting it has no effect on the\n"
             "chunked array. Use 'commitSubarray()' to overwrite data.\n")
        .def("commitSubarray",
             registerConverters(&ChunkedArray_commitSubarray<N, T>),
             (arg("start"), arg("array")),
             "\n    commitSubarray(start, array)\n\n"
             "Write the given 'array' at offset 'start'.\n"
             "The index operator provides a shorthand for this function, e.g.\n"
             "for a 2-dimensional array you can equivalently write::\n\n"
             "    chunked_array.commitSubarray((5,10), roi)\n"
             "    chunked_array[5:12, 10:19] = roi\n\n"
             "to write an ROI with shape (5,7) starting at 'start=(5,10)'.\n")
        .def("releaseChunks",
             &Array::releaseChunks,
             (arg("start"), arg("stop"),arg("destroy")=false),
             "\n    releaseChunks(start, stop, destroy=False)\n\n"
             "\nrelease or destroy all chunks that are completely contained in [start, stop).\n")
        .def("__getitem__", &ChunkedArray_getitem<N, T>,
             "\nRead data from a chunked array with the usual index or slicing syntax::\n\n"
             "    value = chunked_array[5, 20]\n"
             "    roi   = chunked_array[5:12, 10:19]\n\n"
             "Note that the roi is not a slice view of the original array\n"
             "(as in numpy.ndarray), but a copy of the data.\n")
        .def("__setitem__", &ChunkedArray_setitem<N, T>)
        .def("__setitem__", &ChunkedArray_setitem2<N, T>,
             "\nWrite data to a chunked array with the usual index or slicing syntax::\n\n"
             "    chunked_array[5, 20] = value\n"
             "    chunked_array[5:12, 10:19] = roi\n")
        ;

#ifdef HasHDF5
    typedef ChunkedArrayHDF5<N, T> ArrayHDF5;
    class_<ChunkedArrayHDF5<N, T>, bases<Array>, boost::noncopyable>(
         "ChunkedArrayHDF5Base",
         "\n"
         "Base class for HDF5-based chunked arrays, can only be created via\n"
         "the factory function :func:`~vigra.ChunkedArrayHDF5`.\n\n",
         no_init)
        .def("close", &ArrayHDF5::close,
             "\nFlush data to disk and close the underlying HDF5 file.\n")
        .def("flush", &ArrayHDF5::flushToDisk,
             "\nFlush data to disk.\n")
        .add_property("filename", &ArrayHDF5::fileName,
             "\nName of the file backend of this array.\n")
        .add_property("dataset_name", &ArrayHDF5::datasetName,
             "\nName of the dataset backend of this array.\n")
        .add_property("readonly", &ArrayHDF5::isReadOnly,
             "\nTrue if this array is read-only.\n")
    ;
#endif
}

template <unsigned int N>
void defineChunkedArrayFactories(bool export_docu=false)
{
    using namespace boost::python;
    typedef typename MultiArrayShape<N>::type shape_type;

    docstring_options doc_options(true, false, false);

    def("ChunkedArrayFull", &construct_ChunkedArrayFull<N>,
        (arg("shape"), arg("dtype")=defaultDtype(),
         arg("fill_value")=0.0, arg("axistags")=python::object()),
        !export_docu ? "" :
        "Create a chunked array (type :class:`~vigra.vigranumpycore.ChunkedArrayBase`)\n"
        "backed by a plain (consecutive) array::\n\n"
        "  ChunkedArrayFull(shape, dtype=float32, fill_value=0, axistags=None)\n\n"
        "'shape' can be up to 5-dimensional.\n\n"
        "'dtype' can currently be ``uint8``, ``uint32``, and ``float32``.\n\n"
        "'fill_value' is returned for all array elements that have never been written.\n\n"
        "For more details see ChunkedArray_ in the vigra C++ documentation.\n");
    def("ChunkedArrayLazy", &construct_ChunkedArrayLazy<N>,
        (arg("shape"), arg("dtype")=defaultDtype(), arg("chunk_shape")=shape_type(),
         arg("fill_value")=0.0, arg("axistags")=python::object()),
        !export_docu ? "" :
        "Create a chunked array (type :class:`~vigra.vigranumpycore.ChunkedArrayBase`)\n"
        "backed by one plain array for each chunk (rectangular data block)::\n\n"
        "  ChunkedArrayLazy(shape, dtype=float32, chunk_shape=None, fill_value=0, axistags=None)\n\n"
        "The individual chunks are allocated lazily upon first write. Reads before the\n"
        "first write will simply return the 'fill_value' without allocating memory.\n"
        "All allocated chunks reside in memory.\n\n"
        "'shape' can be up to 5-dimensional.\n\n"
        "'chunk_shape' must have the same dimension as 'shape', and its elements must\n"
        "be powers of 2.\n\n"
        "'dtype' can currently be ``uint8``, ``uint32``, and ``float32``.\n\n"
        "'fill_value' is returned for all array elements that have never been written.\n\n"
        "For more details see ChunkedArray_ in the vigra C++ documentation.\n");
    def("ChunkedArrayCompressed", &construct_ChunkedArrayCompressed<N>,
        (arg("shape"), arg("compression")=LZ4, arg("dtype")=defaultDtype(),
         arg("chunk_shape")=shape_type(), arg("cache_max")=-1, arg("fill_value")=0.0,
         arg("axistags")=python::object()),
        !export_docu ? "" :
        "Create a chunked array (type :class:`~vigra.vigranumpycore.ChunkedArrayBase`)\n"
        "backed by one plain array for each chunk (rectangular data block)::\n\n"
        "  ChunkedArrayCompressed(shape, compression=LZ4, dtype=float32, chunk_shape=None, \n"
        "                         cache_max=-1, fill_value=0, axistags=None)\n\n"
        "The individual chunks are allocated lazily upon first write. Reads before the\n"
        "first write will simply return the 'fill_value' without allocating memory.\n"
        "All allocated chunks reside in memory, but may be compressed when not in use.\n"
        "This is especially beneficial for highly compressible data like label images.\n\n"
        "'shape' can be up to 5-dimensional.\n\n"
        "'chunk_shape' must have the same dimension as 'shape', and its elements must\n"
        "be powers of 2.\n\n"
        "'dtype' can currently be ``uint8``, ``uint32``, and ``float32``.\n\n"
        "'fill_value' is returned for all array elements that have never been written.\n\n"
        "'compression' can be any of the flags defined in the :class:`~vigra.Compression` enum.\n\n"
        "'cache_max' specifies how many uncompressed chunks may reside in memory at the\n"
        "same time. If it is '-1', vigra will choose a sensible default, but other values\n"
        "may better fit your data access patterns. This is a soft limit, i.e. may be\n"
        "exceeded temporarily if more chunks are needed simultaneously in a single\n"
        "operation.\n\n"
        "For more details see ChunkedArray_ in the vigra C++ documentation.\n");
    def("ChunkedArrayTmpFile", &construct_ChunkedArrayTmpFile<N>,
        (arg("shape"), arg("dtype")=defaultDtype(), arg("chunk_shape")=shape_type(),
         arg("cache_max")=-1, arg("path")="", arg("fill_value")=0.0,
         arg("axistags")=python::object()),
        !export_docu ? "" :
        "Create a chunked array (type :class:`~vigra.vigranumpycore.ChunkedArrayBase`)\n"
        "backed by a temporary file::\n\n"
        "  ChunkedArrayTmpFile(shape, dtype=float32, chunk_shape=None, cache_max=-1,\n"
        "                      path="", fill_value=0, axistags=None)\n\n"
        "The individual chunks are allocated lazily upon first write. Reads before the\n"
        "first write will simply return the 'fill_value' without allocating memory.\n"
        "Unused chunks will be moved to the file to free their memory. The file is\n"
        "automatically deleted when the object is deleted. Use :func:`~vigra.ChunkedArrayHDF5` if\n"
        "you need persistent storage.\n\n"
        "'shape' can be up to 5-dimensional.\n\n"
        "'chunk_shape' must have the same dimension as 'shape', and its elements must\n"
        "be powers of 2.\n\n"
        "'dtype' can currently be ``uint8``, ``uint32``, and ``float32``.\n\n"
        "'cache_max' specifies how many uncompressed chunks may reside in memory at the\n"
        "same time. If it is '-1', vigra will choose a sensible default, but other values\n"
        "may better fit your data access patterns. This is a soft limit, i.e. may be\n"
        "exceeded temporarily if more chunks are needed simultaneously in a single\n"
        "operation.\n\n"
        "'fill_value' is returned for all array elements that have never been written.\n\n"
        "'path' is the directory where the file is located (default: the system's TMP\n"
        "directory).\n\n"
        "For more details see ChunkedArray_ in the vigra C++ documentation.\n");
}

void defineChunkedArray()
{
    using namespace boost::python;

    docstring_options doc_options(true, false, false);

    enum_<CompressionMethod>("Compression",
         "\nEnum to encode the type of compression for\n"
         ":func:`~vigra.ChunkedArrayCompressed` and :func:`~vigra.ChunkedArrayHDF5`:\n\n"
         "   ``Compression.ZLIB:``\n      ZLIB default compression\n"
         "   ``Compression.ZLIB_NONE:``\n      ZLIB no compression (level = 0)\n"
         "   ``Compression.ZLIB_FAST:``\n      ZLIB fast compression (level = 1)\n"
         "   ``Compression.ZLIB_BEST:``\n      ZLIB best compression (level = 9)\n"
         "   ``Compression.LZ4:``\n      LZ4 compression (very fast)\n\n")
        .value("ZLIB", vigra::ZLIB)
        .value("ZLIB_NONE", vigra::ZLIB_NONE)
        .value("ZLIB_FAST", vigra::ZLIB_FAST)
        .value("ZLIB_BEST", vigra::ZLIB_BEST)
        .value("LZ4", vigra::LZ4)
    ;

#ifdef HasHDF5
    enum_<HDF5File::OpenMode>("HDF5Mode",
         "\nEnum to encode open mode for :func:`~vigra.ChunkedArrayHDF5`:\n\n"
         "   ``HDF5Mode.Default:``\n      Use the default strategy (ReadOnly when file and dataset exist, New otherwise)\n"
         "   ``HDF5Mode.New:``\n      Create new file (existing file will be deleted)\n"
         "   ``HDF5Mode.ReadWrite:``\n      Open file (create when not existing) and allow creation of new datasets.\n"
         "                                  Contents of existing datasets may be changed, but not their shape.\n"
         "   ``HDF5Mode.ReadOnly:``\n     Open files and datasets read-only, fail when not existing.\n"
         "   ``HDF5Mode.Replace:``\n      Like ReadWrite, but always replace existing datasets.\n\n")
        .value("New", HDF5File::New)
        .value("ReadWrite", HDF5File::ReadWrite)
        .value("ReadOnly", HDF5File::ReadOnly)
        .value("Replace", HDF5File::Replace)
        .value("Default", HDF5File::Default)
    ;
#endif

    defineChunkedArrayImpl<2, npy_uint8>();
    defineChunkedArrayImpl<3, npy_uint8>();
    defineChunkedArrayImpl<4, npy_uint8>();
    defineChunkedArrayImpl<5, npy_uint8>();

    defineChunkedArrayImpl<2, npy_uint32>();
    defineChunkedArrayImpl<3, npy_uint32>();
    defineChunkedArrayImpl<4, npy_uint32>();
    defineChunkedArrayImpl<5, npy_uint32>();

    defineChunkedArrayImpl<2, npy_float32>();
    defineChunkedArrayImpl<3, npy_float32>();
    defineChunkedArrayImpl<4, npy_float32>();
    defineChunkedArrayImpl<5, npy_float32>();

    defineChunkedArrayFactories<2>();
    defineChunkedArrayFactories<3>();
    defineChunkedArrayFactories<4>();
    defineChunkedArrayFactories<5>(true);

#ifdef HasHDF5
    def("ChunkedArrayHDF5", &construct_ChunkedArrayHDF5id,
        (arg("file_id"), arg("dataset_name"), arg("shape")=python::object(),
         arg("dtype")=python::object(), arg("mode")=HDF5File::ReadOnly, arg("compression")=ZLIB_FAST,
         arg("chunk_shape")=python::object(), arg("cache_max")=-1, arg("fill_value")=0.0,
         arg("axistags")=python::object()),
        "");
    def("ChunkedArrayHDF5", &construct_ChunkedArrayHDF5,
        (arg("file_name"), arg("dataset_name"), arg("shape")=python::object(),
         arg("dtype")=python::object(), arg("mode")=HDF5File::Default, arg("compression")=ZLIB_FAST,
         arg("chunk_shape")=python::object(), arg("cache_max")=-1, arg("fill_value")=0.0,
         arg("axistags")=python::object()),
        "Create a chunked array (type :class:`~vigra.vigranumpycore.ChunkedArrayHDF5Base`)\n"
        "backed by a HDF5 file::\n\n"
        "  ChunkedArrayHDF5(file, dataset_name, shape=None, dtype=None,\n"
        "                   mode=HDF5Mode.Default, compression=Compression.ZLIB_FAST, \n"
        "                   chunk_shape=None, cache_max=-1, fill_value=0, axistags=None)\n\n"
        "Parameters 'shape', 'dtype', 'compression', 'chunk_shape', 'fill_value', and\n"
        "'axistags' may only be provided when a new dataset is created.\n\n"
        "'file' can be either a file name or a file ID as returned by ``h5py.File.id.id``.\n\n"
        "'shape' can be up to 5-dimensional.\n\n"
        "'chunk_shape' must have the same dimension as 'shape', and its elements must\n"
        "be powers of 2.\n\n"
        "'dtype' can currently be ``uint8``, ``uint32``, and ``float32``.\n\n"
        "'fill_value' is returned for all array elements that have never been written.\n\n"
        "'compression' can be any of the flags defined in the :class:`~vigra.Compression` enum\n"
        "except for `LZ4`.\n\n"
        "'cache_max' specifies how many chunks may reside in memory at the same time.\n"
        "If it is '-1', vigra will choose a sensible default, but other values may\n"
        "better fit your data access patterns. This is a soft limit, i.e. may be exceeded\n"
        "temporarily if more chunks are needed simultaneously in a single operation.\n\n"
        "'mode' defines the access rights to the file and may be any of the flags defined\n"
        "in the :class:`~vigra.HDF5Mode` enum. By default, you get read permission for an existing\n"
        "dataset and read/write permission for a new dataset.\n\n"
        "For more details see ChunkedArray_ in the vigra C++ documentation.\n");
#endif
}

} // namespace vigra

