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
#include <vector>

NUMPY_ARRAY_INITIALIZE_REGISTRY

namespace python = boost::python;

namespace vigra {

// this particular assignment of bits to types is crucial for 
// canonical axis ordering
enum AxisType { Space = 1, Angle = 2, Time = 4, Frequency = 8, 
                UnknownAxisType = 16, Channels = 32 };

class AxisInfo
{
  public:
  
    AxisInfo(std::string key = "?", AxisType typeFlags = UnknownAxisType, 
             double resolution = 0.0, std::string description = "")
    : key_(key),
      description_(description),
      resolution_(resolution),
      flags_(typeFlags)
    {}
    
    std::string key() const
    {
        return key_;
    }
    
    std::string description() const
    {
        return description_;
    }
    
    void setDescription(std::string const & description)
    {
        description_ = description;
    }
    
    double resolution() const
    {
        return resolution_;
    }
    
    void setResolution(double resolution)
    {
        resolution_ = resolution;
    }
    
    AxisType typeFlags() const
    {
        return flags_;
    }
    
    bool isUnknown() const
    {
        return (flags_ == UnknownAxisType) || (flags_ == 0);
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
        if(isUnknown())
        {
            res += " none";
        }
        else
        {
            if(isChannel())
                res += " Channels";
            if(isSpatial())
                res += " Space";
            if(isTemporal())
                res += " Time";
            if(isAngular())
                res += " Angle";
            if(isFrequency())
                res += " Frequency";
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
    
    AxisInfo toFrequencyDomain(unsigned int size = 0) const
    {
        vigra_precondition(!isFrequency(),
            "AxisInfo::toFrequencyDomain(): axis is already in the Fourier domain.");
        AxisInfo res(std::string("f")+key(), AxisType(Frequency | flags_), 0.0, description_);
        if(resolution_ > 0.0 && size > 0u)
            res.resolution_ = 1.0 / (resolution_ * size);
        return res;
    }
    
    AxisInfo fromFrequencyDomain(unsigned int size = 0) const
    {
        vigra_precondition(isFrequency(),
            "AxisInfo::fromFrequencyDomain(): axis is not in the Fourier domain.");
        std::string newKey = key();
        if(newKey[0] == 'f')
            newKey.erase(0,1);
        AxisInfo res(newKey, AxisType(~Frequency & flags_), 0.0, description_);
        if(resolution_ > 0.0 && size > 0u)
            res.resolution_ = 1.0 / (resolution_ * size);
        return res;
    }

	bool operator==(AxisInfo const & other) const
	{
		if(isUnknown() || other.isUnknown())
            return true;
        return key() == other.key() && typeFlags() == other.typeFlags() && 
			   (resolution() == 0.0 || other.resolution() == 0.0 || 
                resolution() == other.resolution());
	}

	bool operator!=(AxisInfo const & other) const
	{
		return !operator==(other);
	}
    
    // factory functions for standard tags
    static AxisInfo x(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("x", Space, resolution, description);
    }
    
    static AxisInfo y(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("y", Space, resolution, description);
    }
    
    static AxisInfo z(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("z", Space, resolution, description);
    }
    
    static AxisInfo t(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("t", Time, resolution, description);
    }
    
    static AxisInfo fx(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("fx", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo fy(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("fy", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo fz(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("fz", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo ft(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("ft", AxisType(Time | Frequency), resolution, description);
    }
    
    static AxisInfo c(std::string const & description = "")
    {
        return AxisInfo("c", Channels, 0.0, description);
    }
    
    std::string key_, description_;
    double resolution_;
    AxisType flags_;
};

AxisInfo AxisInfo__call__(AxisInfo const & i, double resolution, std::string const & description)
{
    return AxisInfo(i.key(), i.typeFlags(), resolution, description);
}

AxisInfo AxisInfo_x()
{
	return AxisInfo::x();
}

AxisInfo AxisInfo_y()
{
    return AxisInfo::y();
}

AxisInfo AxisInfo_z()
{
    return AxisInfo::z();
}

AxisInfo AxisInfo_t()
{
    return AxisInfo::t();
}

AxisInfo AxisInfo_fx()
{
    return AxisInfo::fx();
}

AxisInfo AxisInfo_fy()
{
    return AxisInfo::fy();
}

AxisInfo AxisInfo_fz()
{
    return AxisInfo::fz();
}

AxisInfo AxisInfo_ft()
{
    return AxisInfo::ft();
}

AxisInfo AxisInfo_c()
{
    return AxisInfo::c();
}

// the primary ordering is according to axis type:
//     Space < Angle < Time < Frequency < Unknown < Channels
// the secondary ordering is the lexicographic ordering of the keys
//     "x" < "y" < "z"
template <class GetFunctor>
struct AxisInfoSorter
{
    GetFunctor get;
    
    template <class T>
    bool operator()(T const & l, T const & r) const
    {
		UInt32 lAxisType = get(l).typeFlags();
        if(lAxisType == 0)
            lAxisType = UnknownAxisType;
        UInt32 rAxisType = get(r).typeFlags();
        if(rAxisType == 0)
            rAxisType = UnknownAxisType;
		return (lAxisType == rAxisType)
			        ? get(l).key() < get(r).key()
					: lAxisType < rAxisType;
    }
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
        if(k < 0)
            k += size();
        return getFunctor(axes_[k]);
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
        vigra_precondition(key == "?" || findKey(key) == (int)size(), 
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
    
    std::vector<UInt32> matchOrdering(AxisTags const & other)
    {
        vigra_precondition(size() == other.size(),
            "AxisTags::matchOrdering(): size mismatch.");
        
        std::vector<UInt32> permutation(size());
        for(unsigned int k = 0; k<size(); ++k)
        {
            std::string key = other.get(k).key();
            unsigned int l=0;
            for(; l<size(); ++l)
            {
                if(key == get(l).key())
                    break;
            }
            vigra_precondition(l < size(),
                "AxisTags::matchOrdering(): key mismatch.");
            permutation[k] = l;
        }
        return permutation;
    }
    
    std::vector<UInt32> canonicalOrdering()
    {
        std::vector<UInt32> permutation(size());
        indexSort(axes_.begin(), axes_.end(), permutation.begin(), AxisInfoSorter<GetFunctor>());
        return permutation;
    }
    
    void transpose(std::vector<UInt32> const & permutation)
    {
        vigra_precondition(permutation.size() == size(),
            "AxisTags::transpose(): Permutation has wrong size.");
        std::vector<T> newAxes(size());
        applyPermutation(permutation.begin(), permutation.end(), axes_.begin(), newAxes.begin());
        axes_.swap(newAxes);
    }
    
    void transpose()
    {
        std::reverse(axes_.begin(), axes_.end());
    }

	bool operator==(AxisTags const & other) const
	{
		return std::equal(axes_.begin(), axes_.end(), other.axes_.begin());
	}

	bool operator!=(AxisTags const & other) const
	{
		return !operator==(other);
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
    typedef AxisTags<python::object, PyGetFunctor> BaseType;
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
        if(k < 0)
            k += this->size();
		return this->axes_[k];
    }
    
    python::object getitem(std::string const & key)
    {
        return getitem(this->findKey(key));
    }

    void setitem(int k, python::object i)
    {
        if(!this->checkIndex(k))
		{
			PyErr_SetString(PyExc_IndexError, "AxisInfo::setitem(): Invalid index or key.");
			python::throw_error_already_set();
		}
        if(!python::extract<AxisInfo const &>(i).check())
		{
			PyErr_SetString(PyExc_TypeError, "AxisInfo::setitem(): Item type must be AxisInfo.");
			python::throw_error_already_set();
		}
        
        if(k < 0)
            k += this->size();
		this->axes_[k] = i;
    }
    
    void setitem(std::string const & key, python::object i)
    {
        setitem(this->findKey(key), i);
    }

    void append(python::object i)
    {
        insert(size(), i);
    }

    void insert(int k, python::object i)
    {
        if(k < 0)
            k += this->size();
        if(k < 0)
            k = 0;
        if(k > (int)this->size())
            k = this->size();
        if(!python::extract<AxisInfo const &>(i).check())
		{
			PyErr_SetString(PyExc_TypeError, "AxisInfo::insert(): Item type must be AxisInfo.");
			python::throw_error_already_set();
		}
		this->axes_.insert(this->axes_.begin()+k, i);
    }
    
    void insert(std::string const & key, python::object i)
    {
        insert(this->findKey(key), i);
    }

	python::list axesByFlag(AxisType typeFlags) const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).typeFlags() == typeFlags)
				res.append(k);
		return res;
	}

	python::list spatialAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isSpatial())
				res.append(k);
		return res;
	}

	python::list temporalAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isTemporal())
				res.append(k);
		return res;
	}

	python::list channelAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isChannel())
				res.append(k);
		return res;
	}

	python::list frequencyAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isFrequency())
				res.append(k);
		return res;
	}

	python::list angularAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isAngular())
				res.append(k);
		return res;
	}

	python::list untaggedAxes() const
	{
		python::list res;
		for(unsigned int k=0; k<this->size(); ++k)
			if(this->get(k).isUnknown())
				res.append(k);
		return res;
	}

	template <class U>
    python::list vectorToPython(std::vector<U> const & v) const
	{
		python::list res;
        for(unsigned int k=0; k<v.size(); ++k)
			res.append(v[k]);
		return res;
	}

	python::list canonicalOrdering()
	{
		return vectorToPython(BaseType::canonicalOrdering());
	}

	python::list matchOrdering(PyAxisTags const & other)
	{
		return vectorToPython(BaseType::matchOrdering(other));
	}

	void transpose(python::object const & o)
	{
		unsigned int osize = len(o);
        std::vector<UInt32> permutation(osize);
        
        for(unsigned int k=0; k<this->size(); ++k)
			permutation[k] = python::extract<UInt32>(o[k])();
            
		BaseType::transpose(permutation);
	}

	void transpose()
	{
		BaseType::transpose();
	}
    
    PyAxisTags transform(python::object index, int lnew) const;
};

PyAxisTags 
PyAxisTags::transform(python::object index, int lnew) const
{
    PyAxisTags newTags;
	python::object ellipsis = python::object(python::detail::borrowed_reference(Py_Ellipsis));
	int lold = this->size();
    if(!PySequence_Check(index.ptr()))
    {
        index = python::make_tuple(index);
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
        index += python::make_tuple(ellipsis);
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

UInt32 pychecksum(python::str const & s)
{
    unsigned int size = len(s);
    return checksum(PyString_AsString(s.ptr()), size);
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
        .value("UnknownAxisType", UnknownAxisType)
        .value("Space", Space)
        .value("Time", Time)
        .value("Channels", Channels)
        .value("Frequency", Frequency)
        .value("Angle", Angle)
    ;

    class_<AxisInfo>("AxisInfo", no_init)
        .def(init<std::string, AxisType, double, std::string>(
             (arg("name")="?", arg("typeFlags")=UnknownAxisType, 
              arg("resolution")=0.0, arg("description")="")))
        .def(init<AxisInfo const &>())
		.def_readonly("key", &AxisInfo::key_)
		.def_readwrite("description", &AxisInfo::description_)
		.def_readwrite("resolution", &AxisInfo::resolution_)
		.def_readonly("typeFlags", &AxisInfo::flags_)
		.def("toFrequencyDomain", &AxisInfo::toFrequencyDomain, (arg("size") = 0))
		.def("fromFrequencyDomain", &AxisInfo::fromFrequencyDomain, (arg("size") = 0))
		.def("isSpatial", &AxisInfo::isSpatial)
		.def("isTemporal", &AxisInfo::isTemporal)
		.def("isChannel", &AxisInfo::isChannel)
		.def("isFrequency", &AxisInfo::isFrequency)
		.def("isAngular", &AxisInfo::isAngular)
		.def(self == self)
		.def(self != self)
		.def("__repr__", &AxisInfo::repr)
		.def("__call__", &AxisInfo__call__, (arg("resolution") = 0.0, arg("description") = ""))
        .add_static_property("x", &AxisInfo_x)
        .add_static_property("y", &AxisInfo_y)
        .add_static_property("z", &AxisInfo_z)
        .add_static_property("t", &AxisInfo_t)
        .add_static_property("fx", &AxisInfo_fx)
        .add_static_property("fy", &AxisInfo_fy)
        .add_static_property("fz", &AxisInfo_fz)
        .add_static_property("ft", &AxisInfo_ft)
        .add_static_property("c", &AxisInfo_c)
    ;

	class_<PyAxisTags>("AxisTags", no_init)
        .def(init<object const &, object const &, object const &, 
                  object const &, object const &>(
                  (arg("i1")=object(), arg("i2")=object(), arg("i3")=object(), 
                   arg("i4")=object(), arg("i5")=object())))
		.def("__repr__", &PyAxisTags::repr)
		.def("__getitem__", (object (PyAxisTags::*)(int))&PyAxisTags::getitem)
		.def("__getitem__", (object (PyAxisTags::*)(std::string const &))&PyAxisTags::getitem)
		.def("__setitem__", (void (PyAxisTags::*)(int, object))&PyAxisTags::setitem)
		.def("__setitem__", (void (PyAxisTags::*)(std::string const &, object))&PyAxisTags::setitem)
		.def("insert", (void (PyAxisTags::*)(int, object))&PyAxisTags::insert)
		.def("insert", (void (PyAxisTags::*)(std::string const &, object))&PyAxisTags::insert)
		.def("append", &PyAxisTags::append)
		.def("__delitem__", (void (PyAxisTags::*)(int))&PyAxisTags::dropAxis)
		.def("__delitem__", (void (PyAxisTags::*)(std::string const &))&PyAxisTags::dropAxis)
		.def("__len__", &PyAxisTags::size)
		.def(self == self)
		.def(self != self)
		.def("axesByFlag", &PyAxisTags::axesByFlag)
		.def("spatialAxes", &PyAxisTags::spatialAxes)
		.def("temporalAxes", &PyAxisTags::temporalAxes)
		.def("channelAxes", &PyAxisTags::channelAxes)
		.def("frequencyAxes", &PyAxisTags::frequencyAxes)
		.def("angularAxes", &PyAxisTags::angularAxes)
		.def("untaggedAxes", &PyAxisTags::untaggedAxes)
		.def("canonicalOrdering", &PyAxisTags::canonicalOrdering)
		.def("matchOrdering", &PyAxisTags::matchOrdering)
		.def("transpose", (void (PyAxisTags::*)(object const &))&PyAxisTags::transpose)
		.def("transpose", (void (PyAxisTags::*)())&PyAxisTags::transpose)
		.def("transform", &PyAxisTags::transform)
    ;
    
    def("checksum", &pychecksum, args("data"));
}
