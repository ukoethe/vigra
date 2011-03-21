/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
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
#include <vigra/axistags.hxx>
#include <vigra/numpy_array.hxx>
#include <boost/python.hpp>
#include <boost/python/slice.hpp>

namespace python = boost::python;

namespace vigra {

template<class T>
inline PyObject * managingPyObject(T *p)
{
    return typename python::manage_new_object::apply<T *>::type()(p);
}

template<class Copyable>
python::object
generic__copy__(python::object copyable)
{
    Copyable* newCopyable(new Copyable(python::extract<const Copyable &>(copyable)()));
    python::object result(python::detail::new_reference(managingPyObject(newCopyable)));

    python::extract<python::dict>(result.attr("__dict__"))().update(copyable.attr("__dict__"));

    return result;
}

template<class Copyable>
python::object
generic__deepcopy__(python::object copyable, python::dict memo)
{
    python::object copyMod = python::import("copy");
    python::object deepcopy = copyMod.attr("deepcopy");
    python::object builtin = python::import("__builtin__");
    python::object globals = builtin.attr("__dict__");
    
    Copyable* newCopyable(new Copyable(python::extract<const Copyable &>(copyable)()));
    python::object result(python::detail::new_reference(managingPyObject(newCopyable)));

    python::dict locals;
    locals["copyable"] = copyable;
    int copyableId = python::extract<int>(python::eval("id(copyable)", globals, locals))();
    memo[copyableId] = result;

    python::object dict_copy = deepcopy(python::extract<python::dict>(copyable.attr("__dict__"))(),
                                        memo);    
    python::extract<python::dict>(result.attr("__dict__"))().update(dict_copy);
    return result;
}

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

AxisTags *
AxisTags_create(python::object i1, python::object i2,
                python::object i3, python::object i4, python::object i5)
{
    std::auto_ptr<AxisTags> res(new AxisTags());
    
    python::extract<AxisTags const &> tags(i1);
    if(tags.check())
    {
        res = std::auto_ptr<AxisTags>(new AxisTags(tags()));
    }
    else if(PySequence_Check(i1.ptr()))
    {
        int size = len(i1);
        for(int k=0; k<size; ++k)
        {
            python::extract<AxisInfo const &> info(i1[k]);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
    }
    else if(PyInt_Check(i1.ptr()))
    {
        int size = python::extract<int>(i1)();
        for(int k=0; k<size; ++k)
            res->push_back(AxisInfo());
    }
    else
    {
        if(i1 != python::object())
        {
            python::extract<AxisInfo const &> info(i1);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
        if(i2 != python::object())
        {
            python::extract<AxisInfo const &> info(i2);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
        if(i3 != python::object())
        {
            python::extract<AxisInfo const &> info(i3);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
        if(i4 != python::object())
        {
            python::extract<AxisInfo const &> info(i4);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
        if(i5 != python::object())
        {
            python::extract<AxisInfo const &> info(i5);
            if(!info.check())
            {
                PyErr_SetString(PyExc_TypeError, "AxisTags(): Argument must be a sequence of AxisInfo objects.");
                python::throw_error_already_set();
            }
            res->push_back(info());
        }
    }
    
    return res.release();
}

void AxisTags_insertChannelAxis(AxisTags & axistags)
{
    int k = axistags.channelIndex();
    vigra_precondition(k == (int)axistags.size(),
         "AxisTags::insertChannelAxis(): already has a channel axis.");
    if(detail::defaultOrder() == "F")
        axistags.insert(0, AxisInfo::c());
    else
        axistags.push_back(AxisInfo::c());
}

python::object
AxisTags_permutationToNormalOrder(AxisTags & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToNormalOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationFromNormalOrder(AxisTags & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationFromNormalOrder(permutation);
    return python::object(permutation);
}

AxisTags *
AxisTags_transform(AxisTags const & oldTags, python::object index, int lnew)
{
    std::auto_ptr<AxisTags> newTags(new AxisTags(lnew));
	python::object ellipsis = python::object(python::detail::borrowed_reference(Py_Ellipsis));
	int lold = oldTags.size();
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
            if(item != python::object())
            {
                newTags->get(knew) = oldTags.get(kold);
                // adjust the resolution if item has a valid 'step' member
                python::extract<python::slice> slice(item);
                if(slice.check())
                {
                    python::extract<int> step(slice().step());
                    if(step.check())
                    {
                        newTags->get(knew).resolution_ *= step();
                    }
                }
                ++kold;
            }
            ++knew;
            if(lellipsis > 0 && item == ellipsis)
                --lellipsis;
            else
                ++kindex;
        }
    }
    return newTags.release();
}

// #if 0
// void printAxistags(NumpyAnyArray a)
// {
    // python::object array(python::detail::borrowed_reference(a.pyObject()));
    // python::object tags(getattr(array, "axistags", PyAxisTags()));
	// std::cerr << "Axistags via boost::python:\n";
    // std::cerr << python::extract<PyAxisTags const &>(tags)().repr();

	// std::cerr << "Axistags via C-API:\n";
	// if(PyObject_HasAttrString(a.pyObject(), "axistags"))
	// {
		// python::object tags(python::detail::new_reference(PyObject_GetAttrString(a.pyObject(), "axistags")));
		// std::cerr << python::extract<PyAxisTags const &>(tags)().repr();
	// }
	// else
	// {
		// std::cerr << "attribute 'axistags' missing\n";
	// }
// }
// #endif

void defineAxisTags()
{
    using namespace boost::python;

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
		.def("isType", &AxisInfo::isType)
		.def(self == self)
		.def(self != self)
		.def("__copy__", &generic__copy__<AxisInfo>)
		.def("__deepcopy__", &generic__deepcopy__<AxisInfo>)
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

	class_<AxisTags >("AxisTags", no_init)
		.def("__init__", make_constructor(&AxisTags_create,
            default_call_policies(),
            (arg("i1")=object(), arg("i2")=object(), arg("i3")=object(), 
             arg("i4")=object(), arg("i5")=object())))
		.def("__repr__", &AxisTags::repr)
		.def("__copy__", &generic__copy__<AxisTags>)
		.def("__deepcopy__", &generic__deepcopy__<AxisTags>)
		.def("__len__", &AxisTags::size)
		.def("__getitem__", (AxisInfo const & (AxisTags::*)(int) const)&AxisTags::get,
                             return_value_policy<copy_const_reference>())
		.def("__getitem__", 
            (AxisInfo const & (AxisTags::*)(std::string const &) const)&AxisTags::get,
                             return_value_policy<copy_const_reference>())
		.def("__setitem__", (void (AxisTags::*)(int, AxisInfo const &))&AxisTags::set)
		.def("__setitem__", 
            (void (AxisTags::*)(std::string const &, AxisInfo const &))&AxisTags::set)
		.def("__delitem__", (void (AxisTags::*)(int))&AxisTags::dropAxis)
		.def("__delitem__", (void (AxisTags::*)(std::string const &))&AxisTags::dropAxis)
		.def("insert", &AxisTags::insert)
		.def("append", &AxisTags::push_back)
		.def("dropChannelAxis", &AxisTags::dropChannelAxis)
		.def("insertChannelAxis", &AxisTags_insertChannelAxis)
		.def("swapaxes", &AxisTags::swapaxes)
		.def("transpose", (void (AxisTags::*)())&AxisTags::transpose)
		.def("transpose", (void (AxisTags::*)(ArrayVector<npy_intp> const &))&AxisTags::transpose)
		.def("index", &AxisTags::index)
		.def("resolution", (double (AxisTags::*)(int) const)&AxisTags::resolution)
		.def("resolution", (double (AxisTags::*)(std::string const &) const)&AxisTags::resolution)
		.def("setResolution", (void (AxisTags::*)(int, double))&AxisTags::setResolution)
		.def("setResolution", 
            (void (AxisTags::*)(std::string const &, double))&AxisTags::setResolution)
		.def("scaleAxisResolution", (void (AxisTags::*)(int, double))&AxisTags::scaleAxisResolution)
		.def("scaleAxisResolution", 
               (void (AxisTags::*)(std::string const &, double))&AxisTags::scaleAxisResolution)
		.def("toFrequencyDomain", (void (AxisTags::*)(int, int, int))&AxisTags::toFrequencyDomain,
                (arg("index"), arg("size")=0, arg("sign")=1))
		.def("toFrequencyDomain", 
               (void (AxisTags::*)(std::string const &, int, int))&AxisTags::toFrequencyDomain,
               (arg("key"), arg("size")=0, arg("sign")=1))
		.def("fromFrequencyDomain", (void (AxisTags::*)(int, int))&AxisTags::fromFrequencyDomain,
                (arg("index"), arg("size")=0))
		.def("fromFrequencyDomain", 
               (void (AxisTags::*)(std::string const &, int))&AxisTags::fromFrequencyDomain,
               (arg("key"), arg("size")=0))
		.add_property("channelIndex", &AxisTags::channelIndex)
		.add_property("majorNonchannelIndex", &AxisTags::majorNonchannelIndex)
		.def("axisTypeCount", &AxisTags::axisTypeCount)
		.def("setChannelDescription", &AxisTags::setChannelDescription)
		.def("permutationToNormalOrder", &AxisTags_permutationToNormalOrder)
		.def("permutationFromNormalOrder", &AxisTags_permutationFromNormalOrder)
		.def("transform", &AxisTags_transform,
                             return_value_policy<manage_new_object>())
		.def(self == self)
		.def(self != self)
    ;
}

} // namespace vigra

