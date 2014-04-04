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

#include <vigra/numpy_array.hxx>
#include <vigra/axistags.hxx>
#include <boost/python.hpp>
#include <boost/python/slice.hpp>

// disable axistag construction from json for the moment,
// because these libraries require boost 1.41, whereas
// our Ubuntu longterm support has only boost 1.40
//
// as a workaround, we will do this in Python using the 'json' module
#define VIGRA_DISABLE_FROMJSON

#ifndef VIGRA_DISABLE_FROMJSON  
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#endif

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
    size_t copyableId = python::extract<size_t>(python::eval("id(copyable)", globals, locals))();
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

#ifndef VIGRA_DISABLE_FROMJSON  
AxisTags AxisTags::fromJSON(std::string const & repr)
{
    using boost::property_tree::ptree;
    
    std::istringstream s(repr);
    ptree pt;
    read_json(s, pt);
    
    AxisTags res;
    for(ptree::iterator v = pt.get_child("axes").begin(); 
                         v != pt.get_child("axes").end(); ++v)
    {
        std::string key(v->second.get<std::string>("key"));
        unsigned int typeFlags(v->second.get<unsigned int>("typeFlags"));
        double resolution(v->second.get<double>("resolution"));
        std::string description(v->second.get<std::string>("description"));
        
        res.push_back(AxisInfo(key, (AxisInfo::AxisType)typeFlags, resolution, description));
    }    
    return res;
}

AxisTags * AxisTags_readJSON(std::string const & repr)
{
    return new AxisTags(AxisTags::fromJSON(repr));
}
#endif

AxisTags *
AxisTags_create(python::object i1, python::object i2,
                python::object i3, python::object i4, python::object i5)
{
    VIGRA_UNIQUE_PTR<AxisTags> res(new AxisTags());
    
    python::extract<AxisTags const &> tags(i1);
    if(tags.check())
    {
        res = VIGRA_UNIQUE_PTR<AxisTags>(new AxisTags(tags()));
    }
    else if(PyString_Check(i1.ptr()))
    {
        res = VIGRA_UNIQUE_PTR<AxisTags>(new AxisTags(python::extract<std::string>(i1)()));
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

AxisInfo & AxisTags_getitem(AxisTags & axistags, int index)
{
    if(index < 0)
        index += axistags.size();
        
    if(index >= (int)axistags.size())
    {
        PyErr_SetString(PyExc_IndexError, "AxisTags.__getitem__(): Invalid index or key.");
        python::throw_error_already_set();
    }
    
    return axistags.get(index);
}

std::string AxisTags_str(AxisTags const & axistags)
{
    std::string res;
    for(unsigned int k=0; k<axistags.size(); ++k)
        res += axistags.get(k).repr() + "\n";
    return res;
}

python::object
AxisTags_permutationToNormalOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToNormalOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationToNormalOrder2(AxisTags const & axistags, unsigned int types)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToNormalOrder(permutation, (AxisInfo::AxisType)types);
    return python::object(permutation);
}

python::object
AxisTags_permutationFromNormalOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationFromNormalOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationFromNormalOrder2(AxisTags const & axistags, unsigned int types)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationFromNormalOrder(permutation, (AxisInfo::AxisType)types);
    return python::object(permutation);
}

python::object
AxisTags_permutationToNumpyOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToNumpyOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationFromNumpyOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationFromNumpyOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationToVigraOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToVigraOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationFromVigraOrder(AxisTags const & axistags)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationFromVigraOrder(permutation);
    return python::object(permutation);
}

python::object
AxisTags_permutationToOrder(AxisTags const & axistags, std::string const & order)
{
    ArrayVector<npy_intp> permutation;
    axistags.permutationToOrder(permutation, order);
    return python::object(permutation);
}



AxisTags *
AxisTags_transform(AxisTags const & oldTags, python::object index, int lnew)
{
    VIGRA_UNIQUE_PTR<AxisTags> newTags(new AxisTags());
    python::object ellipsis = python::object(python::detail::borrowed_reference(Py_Ellipsis));
    int lold = oldTags.size();
    if(!PySequence_Check(index.ptr()))
    {
        index = python::make_tuple(index);
    }
    int lindex = len(index);
    int lnewaxis = 0, lellipsis = 0;
    for(int k=0; k<lindex; ++k)
    {
        python::object item(index[k]);
        if(item == python::object() || python::extract<AxisInfo const &>(item).check())
            ++lnewaxis;
        else if(item == ellipsis)
            ++lellipsis;
    }
    lindex -= lnewaxis;
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
                python::extract<AxisInfo const &> newaxis(item);
                
                if(newaxis.check())
                {
                    newTags->push_back(newaxis);
                }
                else
                {
                    newTags->push_back(oldTags.get(kold));
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
            }
            else
            {
                newTags->push_back(AxisInfo());
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

    docstring_options doc_options(true, false, false);

    enum_<AxisInfo::AxisType>("AxisType",
         "\nEnum to encode the type of an axis described by an\n"
         ":class:`~vigra.AxisInfo` object. Possible values:\n\n"
         "   ``AxisType.Channels:``\n      a channel axis\n"
         "   ``AxisType.Space:``\n      a spatial axis\n"
         "   ``AxisType.Angle:``\n      an axis encoding angles (e.g. polar coordinates)\n"
         "   ``AxisType.Time:``\n      a temporal axis\n"
         "   ``AxisType.Frequency:``\n      an axis in the Fourier domain\n"
         "   ``AxisType.UnknownAxisType:``\n      type not specified\n"
         "   ``AxisType.NonChannel:``\n      any type except Channels\n"
         "   ``AxisType.AllAxes:``\n      any type\n\n"
         "Types can be combined by the bitwise 'or' operator. For example,\n"
         "``Space | Frequency`` denotes a spatial axis in the Fourier domain.\n\n")
        .value("UnknownAxisType", AxisInfo::UnknownAxisType)
        .value("Channels", AxisInfo::Channels)
        .value("Space", AxisInfo::Space)
        .value("Angle", AxisInfo::Angle)
        .value("Time", AxisInfo::Time)
        .value("Frequency", AxisInfo::Frequency)
        .value("NonChannel", AxisInfo::NonChannel)
        .value("AllAxes", AxisInfo::AllAxes)
    ;

    class_<AxisInfo>("AxisInfo", 
         "\n"
         "An object describing a single axis.\n\nConstructor:\n\n"
         ".. method:: AxisInfo(key='?', typeFlags=AxisType.UnknownAxisType, resolution=0.0, description='')\n\n"
         "    :param key: the key of the axis,\n"
         "                e.g. 'x' (x-axis), 'c' (channel axis), '?' (unknown)\n" 
         "    :type key: string\n"
         "    :param typeFlags: the type of the axis,\n"
         "                      e.g. AxisType.Space or AxisType.Time\n"
         "    :type typeFlags: :class:`~vigra.AxisType`\n"
         "    :param resolution: the resolution (step size) of the axis\n"
         "                       (e.g. 0.0 means 'unknown')\n"
         "    :param description: an arbitrary string giving additional information \n"
         "                        about the axis.\n\n"
         "In addition, AxisInfo defines the following factories for the most common\n"
         "cases:\n\n"
         "   ``AxisInfo.c`` or ``AxisInfo.c(description='a description')``:\n"
         "        Factory for an axisinfo object describing the 'c' (channel) axis.\n"
         "   ``AxisInfo.x`` or ``AxisInfo.x(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'x' (spatial) axis.\n"
         "   ``AxisInfo.y`` or ``AxisInfo.y(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'y' (spatial) axis.\n"
         "   ``AxisInfo.z`` or ``AxisInfo.z(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'z' (spatial) axis.\n"
         "   ``AxisInfo.t`` or ``AxisInfo.t(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 't' (time) axis.\n"
         "   ``AxisInfo.fx`` or ``AxisInfo.fx(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'x' axis\n"
         "        in the Fourier domain.\n"
         "   ``AxisInfo.fy`` or ``AxisInfo.fy(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'y' axis\n"
         "        in the Fourier domain.\n"
         "   ``AxisInfo.fz`` or ``AxisInfo.fz(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 'z' axis\n"
         "        in the Fourier domain.\n"
         "   ``AxisInfo.ft`` or ``AxisInfo.ft(resolution=0.0, description='')``:\n"
         "        Factory for an axisinfo object describing the 't' axis\n"
         "        in the Fourier domain.\n\n", 
         no_init)
        .def(init<std::string, AxisInfo::AxisType, double, std::string>(
             (arg("key")="?", arg("typeFlags")=AxisInfo::UnknownAxisType, 
              arg("resolution")=0.0, arg("description")="")))
        .def(init<AxisInfo const &>())
        .def_readonly("key", &AxisInfo::key_,
             "\n(read-only property, type 'string') the key of the axis.\n")
        .def_readwrite("description", &AxisInfo::description_,
             "\n(read/write property, type 'string') the string description of the axis.\n")
        .def_readwrite("resolution", &AxisInfo::resolution_,
             "\n(read/write property, type 'float') the resolution of the axis. The resolution\n"
             "will be automatically adjusted whenever the image size changes, e.g. due to a call\n"
             "to :func:`~vigra.sampling.resize` or slicing with non-unit step size::\n\n"
             "    >>> a = vigra.RGBImage((200,100))\n"
             "    >>> a.axistags['x'].resolution = 1.0\n"
             "    >>> a.axistags['y'].resolution = 1.2\n"
             "    >>> print a.axistags\n"
             "    AxisInfo: 'x' (type: Space, resolution=1)\n"
             "    AxisInfo: 'y' (type: Space, resolution=1.2)\n"
             "    AxisInfo: 'c' (type: Channels) RGB\n"
             "    >>> b = a[::2, ::4, :]\n"
             "    >>> print b.axistags\n"
             "    AxisInfo: 'x' (type: Space, resolution=2)\n"
             "    AxisInfo: 'y' (type: Space, resolution=4.8)\n"
             "    AxisInfo: 'c' (type: Channels) RGB\n\n")
        .def_readonly("typeFlags", &AxisInfo::flags_,
             "\n(read-only property, type :class:`~vigra.AxisType`) the type of the axis .\n")
        .def("toFrequencyDomain", &AxisInfo::toFrequencyDomain, (arg("size") = 0, arg("sign") = 1))
        .def("fromFrequencyDomain", &AxisInfo::fromFrequencyDomain, (arg("size") = 0))
        .def("isSpatial", &AxisInfo::isSpatial, 
             "\naxisinfo.isSSpactial() yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains AxisType.Space\n")
        .def("isTemporal", &AxisInfo::isTemporal, 
             "\naxisinfo.isTemporal() yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains AxisType.Time\n")
        .def("isChannel", &AxisInfo::isChannel, 
             "\naxisinfo.isChannel() yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains AxisType.Channels\n")
        .def("isFrequency", &AxisInfo::isFrequency, 
             "\naxisinfo.isFrequency() yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains AxisType.Frequency\n")
        .def("isAngular", &AxisInfo::isAngular, 
             "\naxisinfo.isAngular() yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains AxisType.Angle\n")
        .def("isType", &AxisInfo::isType, 
             "\naxisinfo.isType(axistype) yields True when :attr:`~vigra.AxisInfo.typeFlags` "
             "contains the given axistype.\n")
        .def("compatible", &AxisInfo::compatible, 
             "\naxisinfo1.compatible(axisinfo2) yields True when the two axisinfo objects "
             "have the same keys and types, or either of the two is 'unknown'.\n")
        .def(self == self)
        .def(self != self)
        .def(self < self)
        .def(self <= self)
        .def(self > self)
        .def(self >= self)
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

    class_<AxisTags >("AxisTags", 
            "Object to describe axis properties and axis ordering in a "
            ":class:`~vigra.VigraArray`. \n\nConstructor:\n\n"
            ".. method:: AxisTags(i1=None, i2=None, i3=None, i4=None, i5=None)\n\n"
            "    The parameters 'i1'...'i5' are the :class:`~vigra.AxisInfo` objects\n"
            "    describing the axes. If all are None, an empty AxisTags object is\n"
            "    created. Alternatively, 'i1' can also be a Python sequence of\n"
            "    :class:`~vigra.AxisInfo` objects, or an integer (in which case an\n"
            "    AxisTags object with that many '?' entries is created).\n\n"
            "Most AxisTags methods should not be called directly, but via the\n"
            "corresponding array methods, because this ensures that arrays and\n"
            "their axistags are always kept in sync (rule of thumb: if an axistags\n" 
            "function is not documented, you call it on your own risk).\n\n"
            "The entries of an axistags object (i.e. the individual axisinfo objects)\n"
            "can be accessed via the index operator, where the argument can either be\n"
            "the axis index or the axis key::\n\n"
            "    >>> print array.axistags[0]\n"
            "    AxisInfo: 'x' (type: Space, resolution=1.2)\n"
            "    >>> print array.axistags['x']\n"
            "    AxisInfo: 'x' (type: Space, resolution=1.2)\n"
            "    >>> array.axistags['x'].resolution = 2.0\n"
            "    >>> print array.axistags['x']\n"
            "    AxisInfo: 'x' (type: Space, resolution=2)\n\n",
            no_init)
        .def("__init__", make_constructor(&AxisTags_create,
            default_call_policies(),
            (arg("i1")=object(), arg("i2")=object(), arg("i3")=object(), 
             arg("i4")=object(), arg("i5")=object())))
        .def("__repr__", &AxisTags::repr)
        .def("__str__", &AxisTags_str)
        .def("__copy__", &generic__copy__<AxisTags>)
        .def("__deepcopy__", &generic__deepcopy__<AxisTags>)
        .def("__len__", &AxisTags::size)
        .def("__getitem__", &AxisTags_getitem, return_internal_reference<>())
        .def("__getitem__", 
            (AxisInfo & (AxisTags::*)(std::string const &))&AxisTags::get, return_internal_reference<>())
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
        
        // NOTE: in contrast to arrays, AxisTags::transpose() works
        //       in-place, i.e. changes 'self'
        .def("transpose", (void (AxisTags::*)())&AxisTags::transpose)
        .def("transpose", (void (AxisTags::*)(ArrayVector<npy_intp> const &))&AxisTags::transpose)
        .def("index", &AxisTags::index,
             "Get the axis index of a given axis key::\n\n"
             "    >>> axistags.index('x')\n"
             "    0\n\n"     
             "In this example, the 'x'-axis corresponds to index 0 (i.e. the first index).\n")
        .add_property("channelIndex", &AxisTags::channelIndex,
            "(read-only property, type 'int') the index of the channel axis, or ``len(axistags)``\n"
            "when no channel axis exists (i.e. ``axistags.channelIndex`` is similar to\n"
            "``axistags.index('c')``, but doesn't throw an exception when there\n"
            "is no 'c' axis.)\n\n")
        .add_property("innerNonchannelIndex", &AxisTags::innerNonchannelIndex,
            "(read-only property, type 'int') the index of the innermost non-channel axis.\n")
        .def("axisTypeCount", &AxisTags::axisTypeCount,
            "How many axes of the given type(s) are in this axistags object?::\n\n"
            "    axistags.axisTypeCount(types) -> int\n\n"
            "The 'types' of the query must be single :class:`~vigra.AxisType` instances\n"
            "or a combination of them. Examples::\n\n"
            "    >>> a = vigra.defaultAxistags('txyc')\n"
            "    >>> a.axisTypeCount(vigra.AxisType.Space)\n"
            "    2\n"
            "    >>> a.axisTypeCount(vigra.AxisType.Time)\n"
            "    1\n"
            "    >>> a.axisTypeCount(vigra.AxisType(vigra.AxisType.Space | vigra.AxisType.Time))\n"
            "    3\n"
            "    >>> a.axisTypeCount(vigra.AxisType.NonChannel)\n"
            "    3\n\n")
        
        // IMPORTANT: do not remove are rename the following functions,
        //            they are used by the vigranumpy C++ API
        .def("resolution", (double (AxisTags::*)(int) const)&AxisTags::resolution)
        .def("resolution", (double (AxisTags::*)(std::string const &) const)&AxisTags::resolution)
        .def("setResolution", (void (AxisTags::*)(int, double))&AxisTags::setResolution)
        .def("setResolution", 
            (void (AxisTags::*)(std::string const &, double))&AxisTags::setResolution)
        .def("scaleResolution", (void (AxisTags::*)(int, double))&AxisTags::scaleResolution)
        .def("scaleResolution", 
            (void (AxisTags::*)(std::string const &, double))&AxisTags::scaleResolution)
        .def("description", (std::string (AxisTags::*)(int) const)&AxisTags::description)
        .def("description", 
             (std::string (AxisTags::*)(std::string const &) const)&AxisTags::description)
        .def("setDescription", 
            (void (AxisTags::*)(int, std::string const &))&AxisTags::setDescription)
        .def("setDescription", 
            (void (AxisTags::*)(std::string const &, std::string const &))&AxisTags::setDescription)
        .def("setChannelDescription", &AxisTags::setChannelDescription,
             "Set a description for the channel axis, if one exists::\n\n"
             "    axistags.setChannelDescription('colors are in Lab color space')\n\n"
             "It is similar to::\n\n"
             "    axistags['c'].description = 'colors are in Lab color space'\n\n"
             "except when the axistags contain no channel axis, in which case\n"
             "setChannelDescription() is simply ignored, whereas axistags['c']\n"
             "would cause an exception.\n")
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
        .def("permutationToNormalOrder", &AxisTags_permutationToNormalOrder)
        .def("permutationToNormalOrder", &AxisTags_permutationToNormalOrder2)
        .def("permutationFromNormalOrder", &AxisTags_permutationFromNormalOrder)
        .def("permutationFromNormalOrder", &AxisTags_permutationFromNormalOrder2)
        .def("permutationToNumpyOrder", &AxisTags_permutationToNumpyOrder)
        .def("permutationFromNumpyOrder", &AxisTags_permutationFromNumpyOrder)
        .def("permutationToVigraOrder", &AxisTags_permutationToVigraOrder)
        .def("permutationFromVigraOrder", &AxisTags_permutationFromVigraOrder)
        .def("permutationToOrder", &AxisTags_permutationToOrder)
        .def("transform", &AxisTags_transform,
                             return_value_policy<manage_new_object>())
        .def("compatible", &AxisTags::compatible)
        .def(self == self)
        .def(self != self)
        .def("toJSON", &AxisTags::toJSON,
             "Create a string representation of this axistags in JSON format.")
#ifndef VIGRA_DISABLE_FROMJSON
        .def("fromJSON", &AxisTags_readJSON,
                             return_value_policy<manage_new_object>())
        .staticmethod("fromJSON")
#endif
        ;
}

} // namespace vigra

