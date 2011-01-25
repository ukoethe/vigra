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
#include <vector>

NUMPY_ARRAY_INITIALIZE_REGISTRY

namespace python = boost::python;

namespace vigra {

enum AxisType { Unknown = 0, Space = 1, Time = 2, Channels = 4, Frequency = 8, Angle = 16 };

class AxisInfo
{
  public:
  
    AxisInfo(std::string key = "?", AxisType flags = Unknown, 
             double resolution = 0.0, std::string description = "")
    : key_(key),
      description_(description),
      resolution_(resolution),
      flags_(flags)
    {}
    
    std::string key() const
    {
        return key_;
    }
    
    std::string description() const
    {
        return description_;
    }
    
    double resolution() const
    {
        return resolution_;
    }
    
    AxisType flags() const
    {
        return flags_;
    }
    
    bool isUnknown() const
    {
        return flags_ == Unknown;
    }
    
    bool isSpatial() const
    {
        return (flags_ & Space) != 0;
    }
    
    bool isTemporal() const
    {
        return (flags_ & Time) != 0;
    }
    
    bool isChannel() const
    {
        return (flags_ & Channels) != 0;
    }
    
    bool isFrequency() const
    {
        return (flags_ & Frequency) != 0;
    }
    
    bool isAngular() const
    {
        return (flags_ & Angle) != 0;
    }
    
    std::string repr() const
    {
        std::string res("AxisInfo: '");
        res += key_ + "' (type:";
        if(flags_ == Unknown)
        {
            res += " none";
        }
        else
        {
            if(isSpatial())
                res += " Space";
            if(isTemporal())
                res += " Time";
            if(isChannel())
                res += " Channels";
            if(isFrequency())
                res += " Frequency";
            if(isAngular())
                res += " Angle";
        }
        if(resolution_ > 0.0)
        {
            res += ", resolution=";
            res += asString(resolution_);
        }
        res += ")";
        if(description_ != "")
        {
            res += " ";
            res += description_;
        }
		return res;
    }
    
    std::string key_, description_;
    double resolution_;
    AxisType flags_;
};

// split up AxisTags into a purely C++ part and a Python-dependent part (PyAxisTags below)
template <class T = AxisInfo, class GetFunctor = functor::Identity>
class AxisTags
{
  public:
   
    AxisTags()
    {}
    
    AxisTags(AxisInfo const & i1)
    {
        push_back(i1);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2)
    {
        push_back(i1);
        push_back(i2);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2, AxisInfo const & i3)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2,
             AxisInfo const & i3, AxisInfo const & i4)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
        push_back(i4);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2,
             AxisInfo const & i3, AxisInfo const & i4, AxisInfo const & i5)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
        push_back(i4);
        push_back(i5);
    }
    
    T & operator[](int k)
    {
        return axes_[k];
    }
    
    T const & operator[](int k) const
    {
        return axes_[k];
    }
    
    AxisInfo const & get(int k) const
    {
        vigra_precondition(checkIndex(k),
           "AxisTags::get(): Invalid index or key.");
        return k < 0
                  ? getFunctor(axes_[size()+k])
                  : getFunctor(axes_[k]);
    }

    AxisInfo const & get(std::string const & key) const
    {
        return get(findKey(key));
    }

	unsigned int size() const
	{
		return axes_.size();
	}
    
    std::string repr() const
    {
        std::string res;
        if(size() > 0)
            res += get(0).key();
        for(unsigned int k=1; k<size(); ++k)
        {
            res += " ";
            res += get(k).key();
        }
		return res;
    }
    
    void push_back(T const & i)
    {
        std::string key = getFunctor(i).key();
        vigra_precondition(key == "?" || findKey(key) == size(), 
            std::string("AxisTags::push_back(): duplicate key '" + key + "'."));
        axes_.push_back(i);
    }
    
    void dropAxis(int k)
    {
        vigra_precondition(checkIndex(k),
           "AxisTags::dropAxis(): Invalid index or key.");
        typename std::vector<T>::iterator i = k < 0
                                                 ? axes_.end() + k
                                                 : axes_.begin() + k;
        axes_.erase(i, i+1);
    }
    
    void dropAxis(std::string const & key)
    {
        dropAxis(findKey(key));
    }
    
    bool checkIndex(int k) const
    {
        return k < (int)size() && k >= -(int)size();
    }
    
    int findKey(std::string const & key) const
    {
        for(unsigned int k=0; k<size(); ++k)
            if(get(k).key() == key)
                return k;
        return (int)size();
    }
    
  protected:
	std::vector<T> axes_;
    GetFunctor getFunctor;
};

struct PyGetFunctor
{
	AxisInfo const & operator()(python::object const & o) const
	{
		return python::extract<AxisInfo const &>(o)();
	}
};

class PyAxisTags
: public AxisTags<python::object, PyGetFunctor>
{
  public:
    PyAxisTags()
	{}

    PyAxisTags(python::object i1, python::object i2,
             python::object i3, python::object i4, python::object i5)
    {
        if(PySequence_Check(i1.ptr()))
        {
			int size = len(i1);
            for(int k=0; k<size; ++k)
                if(python::extract<AxisInfo const &>(i1[k]).check())
                    push_back(i1[k]);
        }
        else if(PyInt_Check(i1.ptr()))
        {
            int size = python::extract<int>(i1)();
            for(int k=0; k<size; ++k)
                push_back(python::object(AxisInfo()));
        }
        else
        {
            if(python::extract<AxisInfo const &>(i1).check())
                push_back(i1);
            if(python::extract<AxisInfo const &>(i2).check())
                push_back(i2);
            if(python::extract<AxisInfo const &>(i3).check())
                push_back(i3);
            if(python::extract<AxisInfo const &>(i4).check())
                push_back(i4);
            if(python::extract<AxisInfo const &>(i5).check())
                push_back(i5);
        }
    }
    
    python::object getitem(int k)
    {
        if(!checkIndex(k))
		{
			PyErr_SetString(PyExc_IndexError, "AxisInfo::getitem(): Invalid index or key.");
			python::throw_error_already_set();
		}
		return k < 0 
			      ? this->axes_[size()+k]
				  : this->axes_[k];
    }
    
    python::object getitem(std::string const & key)
    {
        return getitem(findKey(key));
    }

	python::list spatialAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isSpatial())
				res.append(k);
		return res;
	}

	python::list temporalAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isTemporal())
				res.append(k);
		return res;
	}

	python::list channelAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isChannel())
				res.append(k);
		return res;
	}

	python::list frequencyAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isFrequency())
				res.append(k);
		return res;
	}

	python::list angularAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isAngular())
				res.append(k);
		return res;
	}

	python::list untaggedAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<size(); ++k)
			if(get(k).isUnknown())
				res.append(k);
		return res;
	}
    
    PyAxisTags transform(python::object index, int lnew) const;
};

PyAxisTags 
PyAxisTags::transform(python::object index, int lnew) const
{
    PyAxisTags newTags;
	python::object ellipsis = python::object(python::detail::borrowed_reference(Py_Ellipsis));
	int lold = axes_.size();
    if(!PySequence_Check(index.ptr()))
    {
        index = make_tuple(index);
    }
    int lindex = len(index);
    int lnone = 0, lellipsis = 0;
    for(int k=0; k<lindex; ++k)
    {
        python::object item(index[k]);
        if(item == python::object())
            ++lnone;
        else if(item == ellipsis)
            ++lellipsis;
    }
    lindex -= lnone;
    if(lindex < lold && lellipsis == 0)
    {
        index += make_tuple(ellipsis);
        ++lindex;
    }
    lellipsis = lold - lindex;
    int knew = 0, kold = 0, kindex = 0;
    while(knew < lnew)
    {
        python::object item = index[kindex];
        if(PyInt_Check(item.ptr()))
        {
            ++kold;
            ++kindex;
        }
        else 
        {
            if(item == python::object())
            {
                newTags.push_back(python::object(AxisInfo()));
            }
            else
            {
                newTags.push_back(this->axes_[kold]);
                ++kold;
            }
            ++knew;
            if(lellipsis > 0 && item == ellipsis)
                --lellipsis;
            else
                ++kindex;
        }
    }
    return newTags;
}

void registerNumpyArrayConverters();

} // namespace vigra

using namespace boost::python;
using namespace vigra;

BOOST_PYTHON_MODULE_INIT(vigranumpycore)
{
    import_array();
    registerNumpyArrayConverters();

    enum_<AxisType>("AxisType")
        .value("Unknown", Unknown)
        .value("Space", Space)
        .value("Time", Time)
        .value("Channels", Channels)
        .value("Frequency", Frequency)
        .value("Angle", Angle)
    ;

    class_<AxisInfo>("AxisInfo", no_init)
        .def(init<std::string, AxisType, double, std::string>(
             (arg("name")="?", arg("flags")=Unknown, 
              arg("resolution")=0.0, arg("description")="")))
        .def(init<AxisInfo const &>())
		.def_readonly("key", &AxisInfo::key_)
		.def_readwrite("description", &AxisInfo::description_)
		.def_readwrite("resolution", &AxisInfo::resolution_)
		.def_readonly("flags", &AxisInfo::flags_)
		.def("isSpatial", &AxisInfo::isSpatial)
		.def("isTemporal", &AxisInfo::isTemporal)
		.def("isChannel", &AxisInfo::isChannel)
		.def("isFrequency", &AxisInfo::isFrequency)
		.def("isAngular", &AxisInfo::isAngular)
		.def("__repr__", &AxisInfo::repr)
    ;

	class_<PyAxisTags>("AxisTags", no_init)
        .def(init<object const &, object const &, object const &, 
                  object const &, object const &>(
                  (arg("i1")=object(), arg("i2")=object(), arg("i3")=object(), 
                   arg("i4")=object(), arg("i5")=object())))
		.def("__repr__", &PyAxisTags::repr)
		.def("__getitem__", (object (PyAxisTags::*)(int))&PyAxisTags::getitem)
		.def("__getitem__", (object (PyAxisTags::*)(std::string const &))&PyAxisTags::getitem)
		.def("__delitem__", (void (PyAxisTags::*)(int))&PyAxisTags::dropAxis)
		.def("__delitem__", (void (PyAxisTags::*)(std::string const &))&PyAxisTags::dropAxis)
		.def("__len__", &PyAxisTags::size)
		.def("spatialAxes", &PyAxisTags::spatialAxes)
		.def("temporalAxes", &PyAxisTags::temporalAxes)
		.def("channelAxes", &PyAxisTags::channelAxes)
		.def("frequencyAxes", &PyAxisTags::frequencyAxes)
		.def("angularAxes", &PyAxisTags::angularAxes)
		.def("untaggedAxes", &PyAxisTags::untaggedAxes)
		.def("transform", &PyAxisTags::transform)
    ;
}
