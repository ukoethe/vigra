#ifndef CGP2D_CGP2D_PYTHON
#define CGP2D_CGP2D_PYTHON


#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>

#include <sstream>
#include <Python.h>
#include <numpy/arrayobject.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>





namespace vigra{

	typedef vigra::UInt32 CgpLabelType;
	typedef vigra::UInt32 CgpCoordinateType;





	template <typename T>
	inline int typeEnumFromType(void) {
	PyErr_SetString(PyExc_ValueError, "no mapping available for this type");
	boost::python::throw_error_already_set();
	return NPY_VOID;
	}

	//NPY_BOOL
	template <> inline int typeEnumFromType<bool>(void) {
	return NPY_BOOL;
	}


	template <> inline int typeEnumFromType<vigra::UInt8>(void) {
	return NPY_UINT8;
	}

	template <> inline int typeEnumFromType<vigra::UInt16>(void) {
	return NPY_UINT16;
	}

	template <> inline int typeEnumFromType<vigra::UInt32>(void) {
	return NPY_UINT32;
	}

	template <> inline int typeEnumFromType<vigra::UInt64>(void) {
	return NPY_UINT64;
	}

	template <> inline int typeEnumFromType<vigra::Int8>(void) {
	return NPY_INT8;
	}

	template <> inline int typeEnumFromType<vigra::Int16>(void) {
	return NPY_INT16;
	}

	template <> inline int typeEnumFromType<vigra::Int32>(void) {
	return NPY_INT32;
	}

	template <> inline int typeEnumFromType<vigra::Int64>(void) {
	return NPY_INT64;
	}

	template <> inline int typeEnumFromType<float>(void) {
	return NPY_FLOAT32;
	}

	template <> inline int typeEnumFromType<double>(void) {
	return NPY_FLOAT64;
	}



	template<class T>
	static boost::python::object  numpyView1d(T * ptr, const size_t size,const bool isMutable)
	{

	    npy_intp dims[]={static_cast<npy_intp>(size)};
	    void * vptr = static_cast< void *>(ptr);
	    int ndim  = 1;
	    PyObject * pyobj = PyArray_SimpleNewFromData(ndim, dims, NPY_UINT64, vptr );

	    if(!isMutable){
	    

	    	//boost::python::numeric::array boostArray=boost::python::extract<boost::python::numeric::array > (    boost::python::object(boost::python::handle<>(pyobj))   );

	    	//PyArray_CLEARFLAGS( reinterpret_cast<PyArrayObject *>(pyobj), NPY_ARRAY_WRITEABLE);
		}

	    //return pyobj;
	    return boost::python::object(boost::python::handle<>(pyobj));
	}

}

#endif