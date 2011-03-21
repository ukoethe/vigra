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

inline python_ptr pythonFromNumber(bool t)
{
    python_ptr res(PyBool_FromLong(t ? 1 : 0), python_ptr::keep_count);
    pythonToCppException(res);
    return res;
}

#define VIGRA_PYTHON_NUMBER_CONVERSION(type, condition, fct1, fct2) \
inline python_ptr pythonFromNumber(type t) \
{ \
    python_ptr res; \
    if(condition) \
        res = python_ptr(fct1(t), python_ptr::keep_count); \
    else \
        res = python_ptr(fct2(t), python_ptr::keep_count); \
    pythonToCppException(res); \
    return res; \
}

VIGRA_PYTHON_NUMBER_CONVERSION(signed char, true, PyInt_FromLong, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(unsigned char, true, PyInt_FromLong, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(short, true, PyInt_FromLong, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(unsigned short, true, PyInt_FromLong, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(long, true, PyInt_FromLong, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(unsigned long, sizeof(unsigned long) < sizeof(Py_ssize_t), PyInt_FromSsize_t, PyLong_FromUnsignedLongLong)
VIGRA_PYTHON_NUMBER_CONVERSION(int, sizeof(long) < sizeof(Py_ssize_t), PyInt_FromSsize_t, PyInt_FromLong)
VIGRA_PYTHON_NUMBER_CONVERSION(unsigned int, sizeof(unsigned int) < sizeof(Py_ssize_t), PyInt_FromSsize_t, PyLong_FromUnsignedLongLong)
VIGRA_PYTHON_NUMBER_CONVERSION(long long, true, PyLong_FromLongLong, PyLong_FromLongLong)
VIGRA_PYTHON_NUMBER_CONVERSION(unsigned long long, true, PyLong_FromUnsignedLongLong, PyLong_FromUnsignedLongLong)
VIGRA_PYTHON_NUMBER_CONVERSION(float, true, PyFloat_FromDouble, PyFloat_FromDouble)
VIGRA_PYTHON_NUMBER_CONVERSION(double, true, PyFloat_FromDouble, PyFloat_FromDouble)

#undef VIGRA_PYTHON_NUMBER_CONVERSION

template <class T, int N>
python_ptr shapeToPythonTuple(TinyVector<T, N> const & shape)
{
    python_ptr tuple(PyTuple_New(N), python_ptr::keep_count);
    pythonToCppException(tuple);
    for(unsigned int k=0; k<N; ++k)
    {
        PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), k, pythonFromNumber(shape[k]).release());
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
        PyTuple_SET_ITEM((PyTupleObject *)tuple.get(), k, pythonFromNumber(shape[k]).release());
    }
    return tuple;
}

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
