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

#ifndef VIGRA_NUMPY_ARRAY_UTILITIES_HXX
#define VIGRA_NUMPY_ARRAY_UTILITIES_HXX

#include <Python.h>
#include <algorithm>
#include <string>
#include <numpy/arrayobject.h>
#include "python_utility.hxx"

int _import_array();

namespace vigra {

inline void import_vigranumpy()
{
    if(_import_array() < 0)
        pythonToCppException(0);
    python_ptr module(PyImport_ImportModule("vigra.vigranumpycore"), python_ptr::keep_count);
    pythonToCppException(module);
}

/********************************************************/
/*                                                      */
/*                a few Python utilities                */
/*                                                      */
/********************************************************/

namespace detail {

inline long getAttrLong(PyObject * obj, PyObject * key, long defaultValue)
{
    if(!obj)
        return defaultValue;
        
    python_ptr pres(PyObject_GetAttr(obj, key), python_ptr::keep_count);

    if(pres)
        return PyInt_Check(pres)
                     ? PyInt_AsLong(pres)
                     : defaultValue;

    PyErr_Clear();
    return defaultValue;
}

inline long spatialDimensions(PyObject * obj)
{
    static python_ptr key(PyString_FromString("spatialDimensions"), python_ptr::keep_count);
    return getAttrLong(obj, key, -1);
}

inline long channelIndex(PyObject * obj, long defaultVal)
{
    static python_ptr key(PyString_FromString("channelIndex"), python_ptr::keep_count);
    return getAttrLong(obj, key, defaultVal);
}

inline long channelIndex(PyArrayObject * array, long defaultVal)
{
    return channelIndex((PyObject*)array, defaultVal);
}

inline long majorNonchannelIndex(PyArrayObject * array, long defaultVal)
{
    static python_ptr key(PyString_FromString("majorNonchannelIndex"), python_ptr::keep_count);
    return getAttrLong((PyObject*)array, key, defaultVal);
}

inline
python_ptr getArrayTypeObject()
{
    python_ptr vigra(PyImport_ImportModule("vigra"));
    if(vigra)
    {
        static python_ptr key(PyString_FromString("standardArrayType"), python_ptr::keep_count);
        python_ptr arraytype(PyObject_GetAttr(vigra, key), python_ptr::keep_count);
        if(arraytype)
            return arraytype;
    }
    PyErr_Clear();
    return python_ptr((PyObject*)&PyArray_Type);
}

inline std::string defaultOrder(std::string defaultValue = "C")
{
    python_ptr arraytype = getArrayTypeObject();
    static python_ptr key(PyString_FromString("defaultOrder"), python_ptr::keep_count);
    python_ptr pres(PyObject_GetAttr(arraytype, key), python_ptr::keep_count);
    if(pres)
        return PyString_Check(pres)
                     ? PyString_AsString(pres)
                     : defaultValue;
    PyErr_Clear();
    return defaultValue;
}

inline 
python_ptr defaultAxistags(int ndim, std::string order = "")
{
    if(order == "")
        order = defaultOrder();
    python_ptr arraytype = getArrayTypeObject();
    python_ptr func(PyString_FromString("defaultAxistags"), python_ptr::keep_count);
    python_ptr d(PyInt_FromLong(ndim), python_ptr::keep_count);
    python_ptr o(PyString_FromString(order.c_str()), python_ptr::keep_count);
    python_ptr axistags(PyObject_CallMethodObjArgs(arraytype, func.get(), d.get(), o.get(), NULL),
                        python_ptr::keep_count);
    if(axistags)
        return axistags;
    PyErr_Clear();
    return python_ptr();
}

inline 
python_ptr copyAxistags(python_ptr axistags)
{
    if(axistags)
    {
        python_ptr func(PyString_FromString("__copy__"), python_ptr::keep_count);
        axistags = python_ptr(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                              python_ptr::keep_count);
    }
    return axistags;
}

inline void setChannelDescription(python_ptr axistags, std::string const & d)
{
    python_ptr pyd(PyString_FromString(d.c_str()), python_ptr::keep_count);
    python_ptr func(PyString_FromString("setChannelDescription"), python_ptr::keep_count);
    python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), pyd.get(), NULL), 
                   python_ptr::keep_count);
    pythonToCppException(res);
}

inline void scaleAxisResolution(python_ptr axistags, long index, double factor)
{
    python_ptr func(PyString_FromString("scaleAxisResolution"), python_ptr::keep_count);
    python_ptr k(PyInt_FromLong(index), python_ptr::keep_count);
    python_ptr f(PyFloat_FromDouble(factor), python_ptr::keep_count);
    python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), k.get(), f.get(), NULL), 
                   python_ptr::keep_count);
    pythonToCppException(res);
}

inline 
ArrayVector<npy_intp> 
getAxisPermutationImpl(python_ptr object, const char * name, bool ignoreErrors)
{
    python_ptr func(PyString_FromString(name), python_ptr::keep_count);
    python_ptr permutation(PyObject_CallMethodObjArgs(object, func.get(), NULL), 
                           python_ptr::keep_count);
    if(!permutation && ignoreErrors)
    {
        PyErr_Clear();
        return ArrayVector<npy_intp>();
    }
    pythonToCppException(permutation);
    
    if(!PySequence_Check(permutation))
    {
        if(ignoreErrors)
            return ArrayVector<npy_intp>();
        std::string message = std::string(name) + "() did not return a sequence.";
        PyErr_SetString(PyExc_ValueError, message.c_str());
        pythonToCppException(false);
    }
        
    ArrayVector<npy_intp> res(PySequence_Length(permutation));
    for(int k=0; k<(int)res.size(); ++k)
    {
        python_ptr i(PySequence_GetItem(permutation, k), python_ptr::keep_count);
        if(!PyInt_Check(i))
        {
            if(ignoreErrors)
                return ArrayVector<npy_intp>();
            std::string message = std::string(name) + "() did not return a sequence of int.";
            PyErr_SetString(PyExc_ValueError, message.c_str());
            pythonToCppException(false);
        }
        res[k] = PyInt_AsLong(i);
    }
    return res;
}

inline 
ArrayVector<npy_intp> permutationToNormalOrder(python_ptr object, bool ignoreErrors = false)
{
    return getAxisPermutationImpl(object, "permutationToNormalOrder", ignoreErrors);
}

inline 
ArrayVector<npy_intp> permutationFromNormalOrder(python_ptr object, bool ignoreErrors = false)
{
    return getAxisPermutationImpl(object, "permutationFromNormalOrder", ignoreErrors);
}

#if 0

// FIXME: reimplement in terms of TaggedShape?
template <class TINY_VECTOR>
inline
python_ptr constructNumpyArrayFromData(
    std::string const & typeKeyFull,
    std::string const & typeKey,
    TINY_VECTOR const & shape, npy_intp *strides,
    NPY_TYPES typeCode, void *data)
{
    ArrayVector<npy_intp> pyShape(shape.begin(), shape.end());

    python_ptr type = detail::getArrayTypeObject(typeKeyFull);
    if(type == 0)
        type = detail::getArrayTypeObject(typeKey, &PyArray_Type);

    python_ptr array(PyArray_New((PyTypeObject *)type.ptr(), shape.size(), pyShape.begin(), typeCode, strides, data, 0, NPY_WRITEABLE, 0),
                     python_ptr::keep_count);
    pythonToCppException(array);

    return array;
}

#endif

} // namespace detail

} // namespace vigra

#endif // VIGRA_NUMPY_ARRAY_UTILITIES_HXX
