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

#ifndef VIGRA_PYTHON_UTILITY_HXX
#define VIGRA_PYTHON_UTILITY_HXX

#include <Python.h>
#include <algorithm>
#include <string>
#include "vigra/error.hxx"
#include "vigra/tinyvector.hxx"

namespace vigra {

template <class PYOBJECT_PTR>
void pythonToCppException(PYOBJECT_PTR obj)
{
    if(obj != 0)
        return;
    PyObject * type, * value, * trace;
    PyErr_Fetch(&type, &value, &trace);
    if(type == 0)
        return;
    std::string message(((PyTypeObject *)type)->tp_name);
    if(PyString_Check(value))
    {
        message += std::string(": ") + PyString_AS_STRING(value);
    }

    Py_XDECREF(type);
    Py_XDECREF(value);
    Py_XDECREF(trace);
    throw std::runtime_error(message.c_str());
}

/********************************************************/
/*                                                      */
/*                       python_ptr                     */
/*                                                      */
/********************************************************/

class python_ptr
{
  private:
    PyObject * ptr_;

  public:

    typedef PyObject element_type;
    typedef PyObject value_type;
    typedef PyObject * pointer;
    typedef PyObject & reference;

    enum refcount_policy { increment_count, borrowed_reference = increment_count,
                           keep_count, new_reference = keep_count };

    explicit python_ptr(pointer p = 0, refcount_policy rp = increment_count)
    : ptr_( p )
    {
        if(rp == increment_count)
        {
            Py_XINCREF(ptr_);
        }
    }

    python_ptr(python_ptr const & p)
    : ptr_(p.ptr_)
    {
        Py_XINCREF(ptr_);
    }

    python_ptr & operator=(pointer p)
    {
        reset(p);
        return *this;
    }

    python_ptr & operator=(python_ptr const & r)
    {
        reset(r.ptr_);
        return *this;
    }

    ~python_ptr()
    {
        reset();
    }

    void reset(pointer p = 0, refcount_policy rp = increment_count)
    {
        if(p == ptr_)
            return;
        if(rp == increment_count)
        {
            Py_XINCREF(p);
        }
        Py_XDECREF(ptr_);
        ptr_ = p;
    }

    pointer release(bool return_borrowed_reference = false)
    {
        pointer p = ptr_;
        ptr_ = 0;
        if(return_borrowed_reference)
        {
            Py_XDECREF(p);
        }
        return p;
    }

    reference operator* () const
    {
        vigra_precondition(ptr_ != 0, "python_ptr::operator*(): Cannot dereference NULL pointer.");
        return *ptr_;
    }

    pointer operator-> () const
    {
        vigra_precondition(ptr_ != 0, "python_ptr::operator->(): Cannot dereference NULL pointer.");
        return ptr_;
    }

    pointer ptr() const
    {
        return ptr_;
    }

    pointer get() const
    {
        return ptr_;
    }

    operator pointer() const
    {
        return ptr_;
    }

    bool operator! () const
    {
        return ptr_ == 0;
    }

    bool unique() const
    {
        return ptr_ && ptr_->ob_refcnt == 1;
    }

    void swap(python_ptr & other)
    {
        std::swap(ptr_, other.ptr_);
    }

    bool operator==(python_ptr const & p) const
    {
        return ptr_ == p.ptr_;
    }

    bool operator==(pointer p) const
    {
        return ptr_ == p;
    }

    bool operator!=(python_ptr const & p) const
    {
        return ptr_ != p.ptr_;
    }

    bool operator!=(pointer p) const
    {
        return ptr_ != p;
    }
};

inline void swap(python_ptr & a, python_ptr & b)
{
    a.swap(b);
}

/****************************************************************/

inline python_ptr 
makePythonDictionary(char const * k1 = 0, PyObject * a1 = 0,
                    char const * k2 = 0, PyObject * a2 = 0,
                    char const * k3 = 0, PyObject * a3 = 0)
{
    python_ptr dict(PyDict_New(), python_ptr::keep_count);
    pythonToCppException(dict);
    if(k1 && a1)
        PyDict_SetItemString(dict, k1, a1);
    if(k2 && a2)
        PyDict_SetItemString(dict, k2, a2);
    if(k3 && a3)
        PyDict_SetItemString(dict, k3, a3);
    return dict;
}

/****************************************************************/

inline python_ptr pythonFromData(bool t)
{
    python_ptr res(PyBool_FromLong(t ? 1 : 0), python_ptr::keep_count);
    pythonToCppException(res);
    return res;
}

inline python_ptr pythonFromData(std::string const & s)
{
    python_ptr res(PyString_FromString(s.c_str()), python_ptr::keep_count);
    pythonToCppException(res);
    return res;
}

inline python_ptr pythonFromData(long long t)
{
    python_ptr res;
    if(t > (long long)NumericTraits<long>::max() || t < (long long)NumericTraits<long>::min())
        res = python_ptr(PyLong_FromLongLong(t), python_ptr::keep_count);
    else
        res = python_ptr(PyInt_FromLong((long)t), python_ptr::keep_count);
    pythonToCppException(res);
    return res;
}

inline python_ptr pythonFromData(unsigned long long t)
{
    python_ptr res;
    if(t > (unsigned long long)NumericTraits<long>::max())
        res = python_ptr(PyLong_FromUnsignedLongLong(t), python_ptr::keep_count);
    else
        res = python_ptr(PyInt_FromLong((long)t), python_ptr::keep_count);
    pythonToCppException(res);
    return res;
}

#define VIGRA_PYTHON_FROM_DATA(type, fct, cast_type) \
inline python_ptr pythonFromData(type t) \
{ \
    python_ptr res(fct((cast_type)t), python_ptr::keep_count); \
    pythonToCppException(res); \
    return res; \
}

VIGRA_PYTHON_FROM_DATA(signed char, PyInt_FromLong, long)
VIGRA_PYTHON_FROM_DATA(unsigned char, PyInt_FromLong, long)
VIGRA_PYTHON_FROM_DATA(short, PyInt_FromLong, long)
VIGRA_PYTHON_FROM_DATA(unsigned short, PyInt_FromLong, long)
VIGRA_PYTHON_FROM_DATA(long, PyInt_FromLong, long)
VIGRA_PYTHON_FROM_DATA(unsigned long, PyInt_FromSize_t, size_t)
VIGRA_PYTHON_FROM_DATA(int, PyInt_FromSsize_t, Py_ssize_t)
VIGRA_PYTHON_FROM_DATA(unsigned int, PyInt_FromSize_t, size_t)
VIGRA_PYTHON_FROM_DATA(float, PyFloat_FromDouble, double)
VIGRA_PYTHON_FROM_DATA(double, PyFloat_FromDouble, double)
VIGRA_PYTHON_FROM_DATA(char const *, PyString_FromString, char const *)

#undef VIGRA_PYTHON_FROM_DATA

/****************************************************************/

#define VIGRA_DATA_FROM_PYTHON(type, check, extract) \
inline type dataFromPython(PyObject * data, type const & defaultVal) \
{ \
    return data && check(data) \
             ? (type)extract(data) \
             : defaultVal; \
}

VIGRA_DATA_FROM_PYTHON(signed char, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(unsigned char, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(short, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(unsigned short, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(long, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(unsigned long, PyInt_Check, PyInt_AsUnsignedLongMask)
VIGRA_DATA_FROM_PYTHON(int, PyInt_Check, PyInt_AsLong)
VIGRA_DATA_FROM_PYTHON(unsigned int, PyInt_Check, PyInt_AsUnsignedLongMask)
VIGRA_DATA_FROM_PYTHON(long long, PyInt_Check, PyInt_AsSsize_t)
VIGRA_DATA_FROM_PYTHON(unsigned long long, PyInt_Check, PyInt_AsUnsignedLongLongMask)
VIGRA_DATA_FROM_PYTHON(float, PyFloat_Check, PyFloat_AsDouble)
VIGRA_DATA_FROM_PYTHON(double, PyFloat_Check, PyFloat_AsDouble)

inline std::string dataFromPython(PyObject * data, const char * defaultVal) 
{ 
    return data && PyString_Check(data) 
             ? std::string(PyString_AsString(data)) 
             : std::string(defaultVal); 
}

inline std::string dataFromPython(PyObject * data, std::string const & defaultVal) 
{ 
    return data && PyString_Check(data) 
             ? std::string(PyString_AsString(data)) 
             : defaultVal; 
}

inline python_ptr dataFromPython(PyObject * data, python_ptr defaultVal) 
{ 
    return data
             ? python_ptr(data) 
             : defaultVal; 
}

#undef VIGRA_DATA_FROM_PYTHON

/****************************************************************/

template <class T>
T pythonGetAttr(PyObject * obj, const char * key, T defaultValue)
{
    if(!obj)
        return defaultValue;
        
    python_ptr k(PyString_FromString(key), python_ptr::keep_count);
    pythonToCppException(k);
    python_ptr pres(PyObject_GetAttr(obj, k), python_ptr::keep_count);
    if(!pres)
        PyErr_Clear();
    return dataFromPython(pres, defaultValue);
}

inline std::string 
pythonGetAttr(PyObject * obj, const char * key, const char * defaultValue)
{
    if(!obj)
        return std::string(defaultValue);
        
    python_ptr k(PyString_FromString(key), python_ptr::keep_count);
    pythonToCppException(k);
    python_ptr pres(PyObject_GetAttr(obj, k), python_ptr::keep_count);
    if(!pres)
        PyErr_Clear();
    return dataFromPython(pres, defaultValue);
}

/****************************************************************/

template <class T, int N>
python_ptr shapeToPythonTuple(TinyVector<T, N> const & shape)
{
    python_ptr tuple(PyTuple_New(N), python_ptr::keep_count);
    pythonToCppException(tuple);
    for(unsigned int k=0; k<N; ++k)
    {
        PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), k, pythonFromData(shape[k]).release());
    }
    return tuple;
}

template <class T>
python_ptr shapeToPythonTuple(ArrayVectorView<T> const & shape)
{
    python_ptr tuple(PyTuple_New(shape.size()), python_ptr::keep_count);
    pythonToCppException(tuple);
    for(unsigned int k=0; k<shape.size(); ++k)
    {
        PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), k, pythonFromData(shape[k]).release());
    }
    return tuple;
}

/****************************************************************/

class PyAllowThreads
{
    PyThreadState * save_;
    
    // make it non-copyable
    PyAllowThreads(PyAllowThreads const &);
    PyAllowThreads & operator=(PyAllowThreads const &);
  
  public:
    PyAllowThreads()
    : save_(PyEval_SaveThread())
    {}
    
    ~PyAllowThreads()
    {
        PyEval_RestoreThread(save_);
    }
};

} // namespace vigra

#endif  // VIGRA_PYTHON_UTILITY_HXX
