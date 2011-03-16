/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/functorexpression.hxx>
#include <vigra/mathutil.hxx>
#include <vigra/utilities.hxx>
#include <vector>

//NUMPY_ARRAY_DUMMY_REGISTRY
NUMPY_ARRAY_INITIALIZE_REGISTRY

namespace python = boost::python;

namespace vigra {

UInt32 pychecksum(python::str const & s)
{
    unsigned int size = len(s);
    return checksum(PyString_AsString(s.ptr()), size);
}

python::object arraytype(std::string const & cmd)
{
    return python::eval(cmd.c_str());
}

std::string simple(std::string const & cmd)
{
    std::string defaultValue("C");
    
    PyObject *g = PyEval_GetGlobals();
    python_ptr pres(PyRun_String(cmd.c_str(), Py_eval_input, g, g), python_ptr::keep_count);
    if(pres)
        return PyString_Check(pres)
                     ? std::string(PyString_AsString(pres))
                     : defaultValue;
    PyErr_Clear();
    return defaultValue;
    
}

python::list defaultAxisPermutation(int ndim)
{
    PyObject *g = PyEval_GetGlobals();
    std::string command = std::string("vigra.arraytypes.defaultAxisPermutation(") + asString(ndim) + ")";
    python_ptr pres(PyRun_String(command.c_str(), Py_eval_input, g, g), 
                    python_ptr::keep_count);
    python::list res;
    if(pres && PySequence_Check(pres))
    {
        for(int k=0; k<ndim; ++k)
        {
            python_ptr i(PySequence_GetItem(pres, k), python_ptr::keep_count);
            res.append(PyInt_AsLong(i));
        }
    }
    else
    {
        for(int k=0; k<ndim; ++k)
            res.append(k);
    }
    return res;
}

PyObject * defaultAxistags(int ndim)
{
    PyObject *g = PyEval_GetGlobals();
    std::string command = std::string("vigra.arraytypes.defaultAxistags(") + asString(ndim) + ")";
    python_ptr pres(PyRun_String(command.c_str(), Py_eval_input, g, g), 
                    python_ptr::keep_count);
    return pres.release();
}

void setChannelDescription(PyObject * axistags, std::string const & d)
{
    python_ptr pyd(PyString_FromString(d.c_str()), python_ptr::keep_count);
    python_ptr func(PyString_FromString("setChannelDescription"), python_ptr::keep_count);
//    std::cerr << "is function: " << PyCallable_Check(func) << "\n";
    if(!PyObject_CallMethodObjArgs(axistags, func.get(), pyd.get(), NULL))
        python::throw_error_already_set();
}


void registerNumpyArrayConverters();
void defineAxisTags();

} // namespace vigra

using namespace boost::python;
using namespace vigra;

BOOST_PYTHON_MODULE_INIT(vigranumpycore)
{
    import_array();
    registerNumpyArrayConverters();
    defineAxisTags();
    
    def("checksum", &pychecksum, args("data"));
    def("arraytype", &arraytype, (arg("cmd")="vigra.tagged_array.TaggedArray"));
    def("defaultAxisPermutation", &defaultAxisPermutation);
    def("defaultAxistags", &defaultAxistags);
    def("setChannelDescription", &setChannelDescription);
}
