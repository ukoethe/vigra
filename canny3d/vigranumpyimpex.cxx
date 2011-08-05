#include "vigranumpyimpex.hxx"

namespace python = boost::python;

namespace vigra
{
  template< class ValueType, int dim >
  PyObject* createNumpyArray(TinyVector< ValueType, dim > const & shape, 
                 int pyArrayTypeConstant)
  {
//     npy_intp dimensions[dim];
//     for(int i = dim - 1; i >= 0; i--)
//       dimensions[i] = shape[i];

    npy_intp dims[dim];
    for(int i = 0; i < dim; i++)
      dims[i] = shape[i];

    //PyObject* array = PyArray_SimpleNew(dim, dimensions, pyArrayTypeConstant);
    //PyObject* arrayFortranOrder = PyArray_Transpose((PyArrayObject*) array, 
    //						    NULL);
    PyObject* arrayFortranOrder =  PyArray_New(&PyArray_Type, 
                           dim, 
                           dims, 
                           pyArrayTypeConstant, 
                           NULL, 
                           NULL, 
                           0, 
                           1, //fortran
                           NULL);
    //Py_XDECREF(array);
    return PyArray_Return((PyArrayObject*) arrayFortranOrder);
  }

  PyObject* createNumpyArray(int shape, int pyArrayTypeConstant)
  {
    npy_intp dims[1];
    dims[0] = shape;

    PyObject* arrayFortranOrder =  PyArray_New(&PyArray_Type, 
                           1, 
                           dims, 
                           pyArrayTypeConstant, 
                           NULL, 
                           NULL, 
                           0, 
                           1, //fortran
                           NULL);
    //Py_XDECREF(array);
    return PyArray_Return((PyArrayObject*) arrayFortranOrder);
  }

  template
  PyObject* createNumpyArray<int,2>(TinyVector<int,2> const & shape, 
                    int pyArrayTypeConstant);
  template
  PyObject* createNumpyArray<int,3>(TinyVector<int,3> const & shape, 
                    int pyArrayTypeConstant);
  template
  PyObject* createNumpyArray<int,4>(TinyVector<int,4> const & shape, 
                    int pyArrayTypeConstant);


//   template
//   PyObject* createNumpyArray<double,3>(TinyVector<long,3> const & shape, 
//  				       int pyArrayTypeConstant);


  template< class ValueType, int dim, class Tag, int numpyTypeConstant >
  struct MultiArrayViewFromNumpyNdarray
  {
    MultiArrayViewFromNumpyNdarray()
    {
    }
    
    static void* convertible(PyObject* obj)
    {
      //std::cerr << "convertible called." << std::endl;
      return 0;
    }
    
    static void construct(PyObject* obj, python::converter::rvalue_from_python_stage1_data* data)
    {
    }
  };

  template< class ValueType, int dim, int numpyTypeConstant >
  struct MultiArrayViewFromNumpyNdarray< ValueType, 
                     dim,
                     UnstridedArrayTag, 
                     numpyTypeConstant >
{
  MultiArrayViewFromNumpyNdarray()
    {
      python::converter::registry::insert(&convertible, &construct, python::type_id< MultiArrayView< dim, ValueType, UnstridedArrayTag> >());
    }
  
  static void* convertible(PyObject* obj)
    {

      bool iscon =  (PyArray_Check(obj) && // PyArrayObject*?
             PyArray_ISBEHAVED_RO(obj) && 
             // aligned and machine byte order and readable or writeable?
             (PyArray_NDIM(obj) == dim) && // three dimensions?
             (PyArray_TYPE(obj) == numpyTypeConstant) && // UInt8?
             (PyArray_ISFORTRAN(obj) || dim == 1) && // Fortran-contiguous?
             (sizeof(ValueType) == PyArray_ITEMSIZE(obj)));

      if(iscon)
    return obj;
    
//       std::cerr << "convertible called." << std::endl;
//       std::cerr << "check:" << PyArray_Check(obj) << std::endl;
//       std::cerr << "behd :" << PyArray_ISBEHAVED_RO(obj) << std::endl;
//       std::cerr << "dim  :" << PyArray_NDIM(obj) << "!=" << dim  << std::endl;
//       std::cerr << "type :" << (PyArray_TYPE(obj) == numpyTypeConstant) << std::endl;
//       std::cerr << "fort :" << (PyArray_ISFORTRAN(obj) || dim == 1) << std::endl;
//       std::cerr << "size :" << (sizeof(ValueType) == PyArray_ITEMSIZE(obj)) << std::endl;

      return 0;
    }
  
  static void construct(PyObject* obj, python::converter::rvalue_from_python_stage1_data* data)
    {
      void* const storage =   
    ((python::converter::rvalue_from_python_storage<MultiArrayView< 
      dim, ValueType, UnstridedArrayTag > >*) 
     data)->storage.bytes;

      TinyVector<int, dim> shape;
      for(int i = 0; i < dim; i++)
    shape[i] = PyArray_DIMS(obj)[i];
      
      new (storage) MultiArrayView< dim, ValueType, UnstridedArrayTag>
    (shape, (ValueType*) PyArray_DATA(obj));
      
      data->convertible = storage;
//       std::cerr << "constructed." << std::endl;
    }
};

}

using namespace vigra;
using namespace python;

void exportVigraNumpyImpex() 
{
  import_array();

  MultiArrayViewFromNumpyNdarray< UInt64, 
    4, UnstridedArrayTag, NPY_UINT64 >();  
  MultiArrayViewFromNumpyNdarray< UInt64, 
    3, UnstridedArrayTag, NPY_UINT64 >();  
  MultiArrayViewFromNumpyNdarray< UInt64, 
    2, UnstridedArrayTag, NPY_UINT64 >();  
  MultiArrayViewFromNumpyNdarray< UInt64, 
    1, UnstridedArrayTag, NPY_UINT64 >();  

  MultiArrayViewFromNumpyNdarray< UInt32, 
    4, UnstridedArrayTag, NPY_UINT32 >();  
  MultiArrayViewFromNumpyNdarray< UInt32, 
    3, UnstridedArrayTag, NPY_UINT32 >();  
  MultiArrayViewFromNumpyNdarray< UInt32, 
    2, UnstridedArrayTag, NPY_UINT32 >();  
  MultiArrayViewFromNumpyNdarray< UInt32, 
    1, UnstridedArrayTag, NPY_UINT32 >();  

  MultiArrayViewFromNumpyNdarray< Int32, 
    4, UnstridedArrayTag, NPY_INT32 >();  
  MultiArrayViewFromNumpyNdarray< Int32, 
    2, UnstridedArrayTag, NPY_INT32 >();  
  MultiArrayViewFromNumpyNdarray< Int32, 
    1, UnstridedArrayTag, NPY_INT32 >();  

  MultiArrayViewFromNumpyNdarray< Int64, 
    4, UnstridedArrayTag, NPY_INT64 >();  
  MultiArrayViewFromNumpyNdarray< Int64, 
    2, UnstridedArrayTag, NPY_INT64 >();  
  MultiArrayViewFromNumpyNdarray< Int64, 
    1, UnstridedArrayTag, NPY_INT64 >();  

  MultiArrayViewFromNumpyNdarray< UInt16, 
    4, UnstridedArrayTag, NPY_UINT16 >();  
  MultiArrayViewFromNumpyNdarray< UInt16, 
    3, UnstridedArrayTag, NPY_UINT16 >();  
  MultiArrayViewFromNumpyNdarray< UInt16, 
    2, UnstridedArrayTag, NPY_UINT16 >();  
  MultiArrayViewFromNumpyNdarray< UInt16, 
    1, UnstridedArrayTag, NPY_UINT16 >();  

  MultiArrayViewFromNumpyNdarray< bool, 
    1, UnstridedArrayTag, NPY_BOOL >();  

  MultiArrayViewFromNumpyNdarray< double, 
    4, UnstridedArrayTag, NPY_FLOAT64 >();  
  MultiArrayViewFromNumpyNdarray< double, 
    3, UnstridedArrayTag, NPY_FLOAT64 >();  
  MultiArrayViewFromNumpyNdarray< double, 
    2, UnstridedArrayTag, NPY_FLOAT64 >();  
  MultiArrayViewFromNumpyNdarray< double, 
    1, UnstridedArrayTag, NPY_FLOAT64 >();  

  MultiArrayViewFromNumpyNdarray< float, 
    4, UnstridedArrayTag, NPY_FLOAT32 >();  
  MultiArrayViewFromNumpyNdarray< float, 
    2, UnstridedArrayTag, NPY_FLOAT32 >();  
  MultiArrayViewFromNumpyNdarray< float, 
    1, UnstridedArrayTag, NPY_FLOAT32 >();  
}


