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
#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include <sstream>

namespace python = boost::python;

namespace vigra {

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
ChunkedArray_checkoutSubarray(ChunkedArray<N, T> const & self,
                              TinyVector<MultiArrayIndex, N> const & start,
                              TinyVector<MultiArrayIndex, N> const & stop,
                              NumpyArray<N, T> res = NumpyArray<N, T>())
{
    res.reshapeIfEmpty(stop-start, 
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

template <unsigned int N, class T, class Shape>
python::object
bindNumpyArray(NumpyArray<N, T> const & self, Shape const & start, Shape const & stop)
{
    if(allLess(start, stop))
        return python::object(self);
        
    int k=0;
    for(; k<N; ++k)
        if(start[k] == stop[k])
            break;
    
    // FIXME: we should slice the array instead of copying
    NumpyArray<N-1, T> reduced(self.bindAt(k, 0));
    return bindNumpyArray(reduced, start.dropIndex(k), stop.dropIndex(k));
}

template <class T, class Shape>
python::object
bindNumpyArray(NumpyArray<1, T> const & self, Shape const & start, Shape const & stop)
{
    if(start[0] < stop[0])
        return python::object(self);
    else
        return python::object(self[0]);
}

template <unsigned int N, class T>
python::object
ChunkedArray_getitem(ChunkedArray<N, T> const & self, python::object index)
{
    typedef typename MultiArrayShape<N>::type Shape;
    python::object ellipsis = python::object(python::detail::borrowed_reference(Py_Ellipsis));
    if(!PySequence_Check(index.ptr()))
    {
        index = python::make_tuple(index);
    }
    Shape start, stop(self.shape());
    int lindex = len(index);
    int kindex = 0;
    for(; kindex<lindex; ++kindex)
    {
        python::object item = index[kindex];
        if(item == ellipsis)
            break;
    }
    if(kindex == lindex && lindex < N)
    {
        index += python::make_tuple(ellipsis);
        ++lindex;
    }
    kindex = 0;
    for(int k=0; k < N; ++k)
    {
        python::object item = index[kindex];
        if(PyInt_Check(item.ptr()))
        {
            MultiArrayIndex i = python::extract<MultiArrayIndex>(item)();
            start[k] = i;
            if(start[k] < 0)
                start[k] += self.shape(k);
            stop[k] = start[k];
            ++kindex;
        }
        else if(python::extract<python::slice>(item).check())
        {
            python::extract<python::slice> slice(item);
            python::extract<MultiArrayIndex> step(slice().step());
            vigra_precondition(!step.check() || step() == 1,
                "ChunkedArray.__getitem__(): only unit steps are supported.");
            python::extract<MultiArrayIndex> sstart(slice().start());
            if(sstart.check())
            {
                start[k] = sstart();
                if(start[k] < 0)
                    start[k] += self.shape(k);
            }
            python::extract<MultiArrayIndex> sstop(slice().stop());
            if(sstop.check())
            {
                stop[k] = sstop();
                if(stop[k] < 0)
                    stop[k] += self.shape(k);
            }
            ++kindex;
        }
        else if(item == ellipsis)
        {
            if(lindex == N)
                ++kindex;
            else
                ++lindex;
        }
        else
        {
            vigra_precondition(false,
                "ChunkedArray.__getitem__(): unsupported index object.");
        }
    }
    if(start == stop)
    {
        // return a single point
        return python::object(self.getItem(start));
    }
    else if(allLessEqual(start, stop))
    {
        // return a slice
        Shape subStop = max(start + Shape(1), stop);
        NumpyArray<N, T> subarray(subStop - start);
        self.checkoutSubarray(start, subarray);
        return bindNumpyArray(subarray, start, stop);
    }
    else
    {
        vigra_precondition(false,
            "ChunkedArray.__getitem__(): invalid index bounds.");
        return python::object();
    }
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

template <unsigned int N, class T>
PyObject *
ptr_to_python(ChunkedArray<N, T> * a)
{
    return python::to_python_indirect<ChunkedArray<N, T>*, 
                                      python::detail::make_owning_holder>()(a);
}

template <class T, int N>
ChunkedArray<N, T> *
construct_ChunkedArrayFull(TinyVector<MultiArrayIndex, N> const & shape)
{
    return new ChunkedArrayFull<N, T>(shape);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayFull(TinyVector<MultiArrayIndex, N> const & shape,
                           python::object dtype)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayFull<npy_uint8>(shape));
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayFull<npy_uint32>(shape));
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayFull<npy_float32>(shape));
      default:
        vigra_precondition(false, "ChunkedArrayFull(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> * 
construct_ChunkedArrayLazy(TinyVector<MultiArrayIndex, N> const & shape,
                           TinyVector<MultiArrayIndex, N> const & chunk_shape)
{
    return new ChunkedArrayLazy<N, T>(shape, chunk_shape);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayLazy(TinyVector<MultiArrayIndex, N> const & shape,
                           python::object dtype,
                           TinyVector<MultiArrayIndex, N> const & chunk_shape)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayLazy<npy_uint8>(shape, chunk_shape));
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayLazy<npy_uint32>(shape, chunk_shape));
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayLazy<npy_float32>(shape, chunk_shape));
      default:
        vigra_precondition(false, "ChunkedArrayLazy(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> * 
construct_ChunkedArrayCompressed(TinyVector<MultiArrayIndex, N> const & shape,
                                 CompressionMethod method,
                                 TinyVector<MultiArrayIndex, N> const & chunk_shape,
                                 int cache_max)
{
    return new ChunkedArrayCompressed<N, T>(method, shape, chunk_shape, cache_max);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayCompressed(TinyVector<MultiArrayIndex, N> const & shape,
                                 CompressionMethod method,
                                 python::object dtype,
                                 TinyVector<MultiArrayIndex, N> const & chunk_shape,
                                 int cache_max)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayCompressed<npy_uint8>(shape, method, chunk_shape, cache_max));
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayCompressed<npy_uint32>(shape, method, chunk_shape, cache_max));
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayCompressed<npy_float32>(shape, method, chunk_shape, cache_max));
      default:
        vigra_precondition(false, "ChunkedArrayCompressed(): unsupported dtype.");
    }
    return 0;
}

template <class T, int N>
ChunkedArray<N, T> * 
construct_ChunkedArrayTmpFile(TinyVector<MultiArrayIndex, N> const & shape,
                              TinyVector<MultiArrayIndex, N> const & chunk_shape,
                              int cache_max,
                              std::string path)
{
    return new ChunkedArrayTmpFile<N, T>(shape, chunk_shape, cache_max, path);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayTmpFile(TinyVector<MultiArrayIndex, N> const & shape,
                              python::object dtype,
                              TinyVector<MultiArrayIndex, N> const & chunk_shape,
                              int cache_max,
                              std::string path)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayTmpFile<npy_uint8>(shape, chunk_shape, cache_max, path));
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayTmpFile<npy_uint32>(shape, chunk_shape, cache_max, path));
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayTmpFile<npy_float32>(shape, chunk_shape, cache_max, path));
      default:
        vigra_precondition(false, "ChunkedArrayTmpFile(): unsupported dtype.");
    }
    return 0;
}

#ifdef HasHDF5
template <class T, int N>
ChunkedArray<N, T> * 
construct_ChunkedArrayHDF5(std::string filename,
                           std::string datasetName,
                           TinyVector<MultiArrayIndex, N> const & shape,
                           HDF5File::OpenMode mode,
                           CompressionMethod compression,
                           TinyVector<MultiArrayIndex, N> const & chunk_shape,
                           int cache_max)
{
    HDF5File file(filename, mode);
    return new ChunkedArrayHDF5<N, T>(file, datasetName, mode, compression, 
                                      shape, chunk_shape, cache_max);
}

template <unsigned int N>
PyObject *
construct_ChunkedArrayHDF5(std::string filename,
                           std::string datasetName,
                           TinyVector<MultiArrayIndex, N> const & shape,
                           python::object dtype,
                           HDF5File::OpenMode mode,
                           CompressionMethod compression,
                           TinyVector<MultiArrayIndex, N> const & chunk_shape,
                           int cache_max)
{
    switch(numpyScalarTypeNumber(dtype))
    {
      case NPY_UINT8:
        return ptr_to_python(construct_ChunkedArrayHDF5<npy_uint8>(filename, datasetName, shape, mode, compression, chunk_shape, cache_max));
      case NPY_UINT32:
        return ptr_to_python(construct_ChunkedArrayHDF5<npy_uint32>(filename, datasetName, shape, mode, compression, chunk_shape, cache_max));
      case NPY_FLOAT32:
        return ptr_to_python(construct_ChunkedArrayHDF5<npy_float32>(filename, datasetName, shape, mode, compression, chunk_shape, cache_max));
      default:
        vigra_precondition(false, "ChunkedArrayHDF5(): unsupported dtype.");
    }
    return 0;
}
#endif

template <unsigned int N, class T>
void defineChunkedArrayImpl()
{
    using namespace boost::python;

    docstring_options doc_options(true, false, false);

    typedef ChunkedArray<N, T> Array;
    class_<Array, Array*>("ChunkedArray", 
         "\n"
         "Base class for chunked arrays.\n\n", 
         no_init)
        .add_property("shape", &ChunkedArray_shape<N, T>,
             "\nshape of the array.\n")
        .add_property("chunk_shape", &ChunkedArray_chunkShape<N, T>,
             "\nshape of (interior) chunks.\n")
        .add_property("chunk_array_shape", &ChunkedArray_chunkArrayShape<N, T>,
             "\nshape of internal array of chunks.\n")
        .add_property("size", &Array::size,
             "\nnumber of elements of the array.\n")
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
        // .def("axistags", &ChunkedArray_axistags<N, T>,
             // "\n(read-only property, type :class:`~vigra.AxisTags`) the axis tags.\n")
        .def("__repr__", &ChunkedArray_repr<N, T>)
        .def("__str__", &ChunkedArray_str<N, T>)
        .def("checkoutSubarray", 
             registerConverters(&ChunkedArray_checkoutSubarray<N, T>), 
             (arg("start"), arg("stop"), arg("out")=python::object()),
             "\nobtain a copy of the specified subarray.\n")
        .def("commitSubarray", 
             registerConverters(&ChunkedArray_commitSubarray<N, T>),
             (arg("start"), arg("array")),
             "\nwrite the given array at offset 'start'.\n")
        .def("__getitem__", &ChunkedArray_getitem<N, T>)
        // .def("__setitem__", ChunkedArray_setitem)
        ;
        
}

template <unsigned int N>
void defineChunkedArrayFactories()
{
    using namespace boost::python;
    typedef typename MultiArrayShape<N>::type shape_type;

    docstring_options doc_options(true, false, false);
    
    PyObject * dtype = NumpyArrayValuetypeTraits<npy_float32>::typeObject();
    object default_dtype = object(python::detail::new_reference(dtype));

    def("ChunkedArrayFull", &construct_ChunkedArrayFull<N>,
        (arg("shape"), arg("dtype")=default_dtype));
    def("ChunkedArrayLazy", &construct_ChunkedArrayLazy<N>,
        (arg("shape"), arg("dtype")=default_dtype, arg("chunk_shape")=shape_type()));
    def("ChunkedArrayCompressed", &construct_ChunkedArrayCompressed<N>,
        (arg("shape"), arg("method")=LZ4, arg("dtype")=default_dtype, 
        arg("chunk_shape")=shape_type(), arg("cache_max")=-1));
    def("ChunkedArrayTmpFile", &construct_ChunkedArrayTmpFile<N>,
        (arg("shape"), arg("dtype")=default_dtype, arg("chunk_shape")=shape_type(), 
         arg("cache_max")=-1, arg("path")=""));
#ifdef HasHDF5
    def("ChunkedArrayHDF5", &construct_ChunkedArrayHDF5<N>,
        (arg("file_name"), arg("dataset_name"), arg("shape")=shape_type(),
         arg("dtype")=default_dtype, arg("mode")=HDF5File::New, arg("compression")=ZLIB_FAST, 
         arg("chunk_shape")=shape_type(), arg("cache_max")=-1));
#endif
}

void defineChunkedArray()
{
    using namespace boost::python;

    docstring_options doc_options(true, false, false);

    enum_<CompressionMethod>("Compression",
         "\nEnum to encode the type of compression for\n"
         "ChunkedArrayCompressed and ChunkedArrayHDF5:\n\n"
         "   ``Compression.ZLIB:``\n      ZLIB default compression\n"
         "   ``Compression.ZLIB_NONE:``\n      ZLIB no compression (level = 0)\n"
         "   ``Compression.ZLIB_FAST:``\n      ZLIB fast compression (level = 1)\n"
         "   ``Compression.ZLIB_BEST:``\n      ZLIB best compression (level = 9)\n"
         "   ``Compression.LZ4:``\n      LZ4 compression (very fast)\n\n")
        .value("ZLIB", CompressionMethod::ZLIB)
        .value("ZLIB_NONE", CompressionMethod::ZLIB_NONE)
        .value("ZLIB_FAST", CompressionMethod::ZLIB_FAST)
        .value("ZLIB_BEST", CompressionMethod::ZLIB_BEST)
        .value("LZ4", CompressionMethod::LZ4)
    ;

#ifdef HasHDF5
    enum_<HDF5File::OpenMode>("HDF5Mode",
         "\nEnum to encode open mode for ChunkedArrayHDF5:\n\n"
         "   ``HDF5Mode.New:``\n      Create new dataset (existing dataset will be deleted)\n"
         "   ``HDF5Mode.Open:``\n      Open dataset, create when not existing\n"
         "   ``HDF5Mode.OpenReadOnly:``\n      Open dataset read-only, fail when not existing\n\n")
        .value("New", HDF5File::New)
        .value("Open", HDF5File::Open)
        .value("OpenReadOnly", HDF5File::OpenReadOnly)
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
    defineChunkedArrayFactories<5>();
}

} // namespace vigra

