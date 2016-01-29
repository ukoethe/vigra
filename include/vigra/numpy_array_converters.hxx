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

#ifndef VIGRA_NUMPY_ARRAY_CONVERTERS_HXX
#define VIGRA_NUMPY_ARRAY_CONVERTERS_HXX

#include "numpy_array.hxx"
#include "metaprogramming.hxx"
#include <boost/python.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/raw_function.hpp>
#include <set>
#include <type_traits>

namespace vigra {

template <class Array>
PyObject * returnNumpyArray(Array const & a)
{
    PyObject * pa = a.pyObject();
    if(pa == 0)
        PyErr_SetString(PyExc_ValueError, "returnNumpyArray(): Conversion to Python failed, array has no data.");
    else
        Py_INCREF(pa);
    return pa;
}

VIGRA_EXPORT std::set<std::string> & exportedArrayKeys();

template <class ArrayType>
struct NumpyArrayConverter {};

template <unsigned int N, class T, class Stride>
struct NumpyArrayConverter<NumpyArray<N, T, Stride> >
{
    typedef NumpyArray<N, T, Stride> ArrayType;
    typedef typename ArrayType::ArrayTraits ArrayTraits;

    NumpyArrayConverter();

    static void* convertible(PyObject* obj);

    // from Python
    static void construct(PyObject* obj,
        boost::python::converter::rvalue_from_python_stage1_data* data);

    // to Python
    static PyObject* convert(ArrayType const& a)
    {
        return returnNumpyArray(a);
    }
};

template <unsigned int N, class T, class Stride>
NumpyArrayConverter<NumpyArray<N, T, Stride> >::NumpyArrayConverter()
{
    using namespace boost::python;

    converter::registration const * reg = converter::registry::query(type_id<ArrayType>());

    // register the to_python_converter only once
    // FIXME: I'm not sure if this is correct.
    if(!reg || !reg->rvalue_chain)
    {
        to_python_converter<ArrayType, NumpyArrayConverter>();
    }
    converter::registry::insert(&convertible, &construct, type_id<ArrayType>());
}

template <unsigned int N, class T, class Stride>
void * NumpyArrayConverter<NumpyArray<N, T, Stride> >::convertible(PyObject* obj)
{
    bool isCompatible = obj == Py_None || ArrayType::isStrictlyCompatible(obj);
    // std::cerr << "compatible for " << typeid(NumpyArray<N, T, Stride>).name() << ": " << isCompatible << "\n";
    return isCompatible
             ? obj
             : 0;
}

// from Python
template <unsigned int N, class T, class Stride>
void NumpyArrayConverter<NumpyArray<N, T, Stride> >::construct(PyObject* obj,
                   boost::python::converter::rvalue_from_python_stage1_data* data)
{
    void* const storage =
        ((boost::python::converter::rvalue_from_python_storage<ArrayType>* ) data)->storage.bytes;

    ArrayType * array = new (storage) ArrayType();
    if(obj != Py_None)
        array->makeReferenceUnchecked(obj);

    data->convertible = storage;
}

template <unsigned int N, class T, class Stride>
struct NumpyArrayConverter<MultiArrayView<N, T, Stride> >
: public NumpyArrayConverter<NumpyArray<N, T, Stride> >
{
    typedef NumpyArrayConverter<NumpyArray<N, T, Stride> > BaseType;
    typedef MultiArrayView<N, T, Stride> ArrayType;

    NumpyArrayConverter()
    {
        using namespace boost::python;
        converter::registry::insert(&BaseType::convertible, &BaseType::construct,
                                    type_id<ArrayType>());
    }
};

template <class Iter, class End>
struct RegisterNumpyArrayConverters
{
    static void exec()
    {
        typedef typename UnqualifiedType<typename boost::mpl::deref<Iter>::type>::type Type;
        NumpyArrayConverter<Type>();
        RegisterNumpyArrayConverters<typename boost::mpl::next<Iter>::type, End>::exec();
    }
};

template <class End>
struct RegisterNumpyArrayConverters<End, End>
{
    static void exec()
    {}
};

template <class Typelist>
void registerNumpyArrayConverters(Typelist)
{
    RegisterNumpyArrayConverters<typename boost::mpl::begin<Typelist>::type,
                                 typename boost::mpl::end<Typelist>::type >::exec();
}

template <class FN>
FN registerConverters(FN f)
{
    registerNumpyArrayConverters(boost::python::detail::get_signature(f));
    return f;
}

namespace detail {

template <class T>
struct TypeName;

template <class T>
struct TypeName<Singleband<T>>
: public TypeName<T>
{};

template <class T>
struct TypeName<Multiband<T>>
: public TypeName<T>
{};

template <class T, int N>
struct TypeName<TinyVector<T, N>>
: public TypeName<T>
{};

template <>
struct TypeName<void>
{
    static std::string name() {
        return std::string("void");
    }
    static std::string sized_name() {
        return std::string("void");
    }
};

#define VIGRA_SIGNED_INT_NAME(type) \
template <> \
struct TypeName<type> \
{ \
    static std::string name() { \
        return std::string(#type); \
    } \
    static std::string sized_name() { \
        return std::string("int") + std::to_string(sizeof(type)*8); \
    } \
};

VIGRA_SIGNED_INT_NAME(signed char)
VIGRA_SIGNED_INT_NAME(short)
VIGRA_SIGNED_INT_NAME(int)
VIGRA_SIGNED_INT_NAME(long)
VIGRA_SIGNED_INT_NAME(long long)

#define VIGRA_UNSIGNED_INT_NAME(type) \
template <> \
struct TypeName<type> \
{ \
    static std::string name() { \
        return std::string(#type); \
    } \
    static std::string sized_name() { \
        return std::string("uint") + std::to_string(sizeof(type)*8); \
    } \
};

VIGRA_UNSIGNED_INT_NAME(unsigned char)
VIGRA_UNSIGNED_INT_NAME(unsigned short)
VIGRA_UNSIGNED_INT_NAME(unsigned int)
VIGRA_UNSIGNED_INT_NAME(unsigned long)
VIGRA_UNSIGNED_INT_NAME(unsigned long long)

#define VIGRA_FLOAT_NAME(type) \
template <> \
struct TypeName<type> \
{ \
    static std::string name() { \
        return std::string(#type); \
    } \
    static std::string sized_name() { \
        return std::string("float") + std::to_string(sizeof(type)*8); \
    } \
};

VIGRA_FLOAT_NAME(float)
VIGRA_FLOAT_NAME(double)
VIGRA_FLOAT_NAME(long double)

#undef VIGRA_SIGNED_INT_NAME
#undef VIGRA_UNSIGNED_INT_NAME
#undef VIGRA_FLOAT_NAME

} // namespace detail

} // namespace vigra

namespace boost { namespace python {

#define VIGRA_PYTHON_MULTITYPE_FUNCTOR(functor_name, function) \
template <class T> \
struct functor_name##Impl \
{ \
    static void def(const char * pythonName) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T>)); \
    } \
     \
    template <class A1> \
    static void def(const char * pythonName, A1 const & a1) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T>), a1); \
    } \
     \
    template <class A1, class A2> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T>), a1, a2); \
    } \
     \
    template <class A1, class A2, class A3> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2, A3 const & a3) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T>), a1, a2, a3); \
    } \
}; \
 \
template <> \
struct functor_name##Impl<void> \
{ \
    static void def(const char *) {} \
     \
    template <class A1> \
    static void def(const char *, A1 const &) {} \
     \
    template <class A1, class A2> \
    static void def(const char *, A1 const &, A2 const &) {} \
     \
    template <class A1, class A2, class A3> \
    static void def(const char *, A1 const &, A2 const &, A3 const &) {} \
}; \
 \
template <class T1, \
          class T2 = void, \
          class T3 = void, \
          class T4 = void, \
          class T5 = void, \
          class T6 = void, \
          class T7 = void, \
          class T8 = void, \
          class T9 = void, \
          class T10 = void, \
          class T11 = void, \
          class T12 = void> \
struct functor_name \
: public boost::python::PythonMultidefFunctor \
{ \
    bool install_fallback_; \
    \
    functor_name(bool install_fallback = false) \
    : install_fallback_(install_fallback) \
    {} \
    \
    typedef boost::python::ArgumentMismatchMessage\
        <T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> Message; \
    typedef functor_name##Impl<T1 > F1; \
    typedef functor_name##Impl<T2 > F2; \
    typedef functor_name##Impl<T3 > F3; \
    typedef functor_name##Impl<T4 > F4; \
    typedef functor_name##Impl<T5 > F5; \
    typedef functor_name##Impl<T6 > F6; \
    typedef functor_name##Impl<T7 > F7; \
    typedef functor_name##Impl<T8 > F8; \
    typedef functor_name##Impl<T9 > F9; \
    typedef functor_name##Impl<T10> F10; \
    typedef functor_name##Impl<T11> F11; \
    typedef functor_name##Impl<T12> F12; \
    \
    void def(const char * pythonName) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName); \
        F1 ::def(pythonName); \
        F2 ::def(pythonName); \
        F3 ::def(pythonName); \
        F4 ::def(pythonName); \
        F5 ::def(pythonName); \
        F6 ::def(pythonName); \
        F7 ::def(pythonName); \
        F8 ::def(pythonName); \
        F9 ::def(pythonName); \
        F10::def(pythonName); \
        F11::def(pythonName); \
        F12::def(pythonName); \
    } \
     \
    template <class A1> \
    void def(const char * pythonName, A1 const & a1) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1); \
        F1 ::def(pythonName, a1); \
        F2 ::def(pythonName, a1); \
        F3 ::def(pythonName, a1); \
        F4 ::def(pythonName, a1); \
        F5 ::def(pythonName, a1); \
        F6 ::def(pythonName, a1); \
        F7 ::def(pythonName, a1); \
        F8 ::def(pythonName, a1); \
        F9 ::def(pythonName, a1); \
        F10::def(pythonName, a1); \
        F11::def(pythonName, a1); \
        F12::def(pythonName, a1); \
    } \
     \
    template <class A1, class A2> \
    void def(const char * pythonName, A1 const & a1, A2 const & a2) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1, a2); \
        F1 ::def(pythonName, a1, a2); \
        F2 ::def(pythonName, a1, a2); \
        F3 ::def(pythonName, a1, a2); \
        F4 ::def(pythonName, a1, a2); \
        F5 ::def(pythonName, a1, a2); \
        F6 ::def(pythonName, a1, a2); \
        F7 ::def(pythonName, a1, a2); \
        F8 ::def(pythonName, a1, a2); \
        F9 ::def(pythonName, a1, a2); \
        F10::def(pythonName, a1, a2); \
        F11::def(pythonName, a1, a2); \
        F12::def(pythonName, a1, a2); \
    } \
     \
    template <class A1, class A2, class A3> \
    void def(const char * pythonName, A1 const & a1, A2 const & a2, A3 const & a3) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1, a2, a3); \
        F1 ::def(pythonName, a1, a2, a3); \
        F2 ::def(pythonName, a1, a2, a3); \
        F3 ::def(pythonName, a1, a2, a3); \
        F4 ::def(pythonName, a1, a2, a3); \
        F5 ::def(pythonName, a1, a2, a3); \
        F6 ::def(pythonName, a1, a2, a3); \
        F7 ::def(pythonName, a1, a2, a3); \
        F8 ::def(pythonName, a1, a2, a3); \
        F9 ::def(pythonName, a1, a2, a3); \
        F10::def(pythonName, a1, a2, a3); \
        F11::def(pythonName, a1, a2, a3); \
        F12::def(pythonName, a1, a2, a3); \
    } \
};

#define VIGRA_PYTHON_MULTITYPE_FUNCTOR_NDIM(functor_name, function) \
template <class T, int FROM, int TO> \
struct functor_name##Impl \
{ \
    typedef functor_name##Impl type; \
     \
    static void def(const char * pythonName) \
    { \
        functor_name##Impl<T, FROM, FROM>::def(pythonName); \
        functor_name##Impl<T, FROM+1, TO>::def(pythonName); \
    } \
     \
    template <class A1> \
    static void def(const char * pythonName, A1 const & a1) \
    { \
        functor_name##Impl<T, FROM, FROM>::def(pythonName, a1); \
        functor_name##Impl<T, FROM+1, TO>::def(pythonName, a1); \
    } \
     \
    template <class A1, class A2> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2) \
    { \
        functor_name##Impl<T, FROM, FROM>::def(pythonName, a1, a2); \
        functor_name##Impl<T, FROM+1, TO>::def(pythonName, a1, a2); \
    } \
     \
    template <class A1, class A2, class A3> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2, A3 const & a3) \
    { \
        functor_name##Impl<T, FROM, FROM>::def(pythonName, a1, a2, a3); \
        functor_name##Impl<T, FROM+1, TO>::def(pythonName, a1, a2, a3); \
    } \
}; \
 \
template <class T, int N> \
struct functor_name##Impl<T, N, N> \
{ \
    typedef functor_name##Impl type; \
     \
    static void def(const char * pythonName) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T, N>)); \
    } \
     \
    template <class A1> \
    static void def(const char * pythonName, A1 const & a1) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T, N>), a1); \
    } \
     \
    template <class A1, class A2> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T, N>), a1, a2); \
    } \
     \
    template <class A1, class A2, class A3> \
    static void def(const char * pythonName, A1 const & a1, A2 const & a2, A3 const & a3) \
    { \
        boost::python::def(pythonName, vigra::registerConverters(&function<T, N>), a1, a2, a3); \
    } \
}; \
 \
template <int FROM, int TO> \
struct functor_name##Impl<void, FROM, TO> \
{ \
    static void def(const char *) {} \
     \
    template <class A1> \
    static void def(const char *, A1 const &) {} \
     \
    template <class A1, class A2> \
    static void def(const char *, A1 const &, A2 const &) {} \
     \
    template <class A1, class A2, class A3> \
    static void def(const char *, A1 const &, A2 const &, A3 const &) {} \
}; \
 \
template <int N> \
struct functor_name##Impl<void, N, N> \
{ \
    static void def(const char *) {} \
     \
    template <class A1> \
    static void def(const char *, A1 const &) {} \
     \
    template <class A1, class A2> \
    static void def(const char *, A1 const &, A2 const &) {} \
     \
    template <class A1, class A2, class A3> \
    static void def(const char *, A1 const &, A2 const &, A3 const &) {} \
}; \
 \
template <int FROM, int TO, \
          class T1, \
          class T2 = void, \
          class T3 = void, \
          class T4 = void, \
          class T5 = void, \
          class T6 = void, \
          class T7 = void, \
          class T8 = void, \
          class T9 = void, \
          class T10 = void, \
          class T11 = void, \
          class T12 = void> \
struct functor_name \
: public boost::python::PythonMultidefFunctor \
{ \
    bool install_fallback_; \
    \
    functor_name(bool install_fallback = false) \
    : install_fallback_(install_fallback) \
    { \
        static_assert(FROM <= TO, #functor_name ": dimension range empty (FROM > TO)"); \
    } \
    \
    typedef boost::python::ArgumentMismatchMessage\
        <T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> Message; \
    typedef functor_name##Impl<T1 , FROM, TO> F1; \
    typedef functor_name##Impl<T2 , FROM, TO> F2; \
    typedef functor_name##Impl<T3 , FROM, TO> F3; \
    typedef functor_name##Impl<T4 , FROM, TO> F4; \
    typedef functor_name##Impl<T5 , FROM, TO> F5; \
    typedef functor_name##Impl<T6 , FROM, TO> F6; \
    typedef functor_name##Impl<T7 , FROM, TO> F7; \
    typedef functor_name##Impl<T8 , FROM, TO> F8; \
    typedef functor_name##Impl<T9 , FROM, TO> F9; \
    typedef functor_name##Impl<T10, FROM, TO> F10; \
    typedef functor_name##Impl<T11, FROM, TO> F11; \
    typedef functor_name##Impl<T12, FROM, TO> F12; \
    \
    void def(const char * pythonName) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName); \
        F1 ::def(pythonName); \
        F2 ::def(pythonName); \
        F3 ::def(pythonName); \
        F4 ::def(pythonName); \
        F5 ::def(pythonName); \
        F6 ::def(pythonName); \
        F7 ::def(pythonName); \
        F8 ::def(pythonName); \
        F9 ::def(pythonName); \
        F10::def(pythonName); \
        F11::def(pythonName); \
        F12::def(pythonName); \
    } \
     \
    template <class A1> \
    void def(const char * pythonName, A1 const & a1) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1); \
        F1 ::def(pythonName, a1); \
        F2 ::def(pythonName, a1); \
        F3 ::def(pythonName, a1); \
        F4 ::def(pythonName, a1); \
        F5 ::def(pythonName, a1); \
        F6 ::def(pythonName, a1); \
        F7 ::def(pythonName, a1); \
        F8 ::def(pythonName, a1); \
        F9 ::def(pythonName, a1); \
        F10::def(pythonName, a1); \
        F11::def(pythonName, a1); \
        F12::def(pythonName, a1); \
    } \
     \
    template <class A1, class A2> \
    void def(const char * pythonName, A1 const & a1, A2 const & a2) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1, a2); \
        F1 ::def(pythonName, a1, a2); \
        F2 ::def(pythonName, a1, a2); \
        F3 ::def(pythonName, a1, a2); \
        F4 ::def(pythonName, a1, a2); \
        F5 ::def(pythonName, a1, a2); \
        F6 ::def(pythonName, a1, a2); \
        F7 ::def(pythonName, a1, a2); \
        F8 ::def(pythonName, a1, a2); \
        F9 ::def(pythonName, a1, a2); \
        F10::def(pythonName, a1, a2); \
        F11::def(pythonName, a1, a2); \
        F12::def(pythonName, a1, a2); \
    } \
     \
    template <class A1, class A2, class A3> \
    void def(const char * pythonName, A1 const & a1, A2 const & a2, A3 const & a3) const \
    { \
        if(install_fallback_) \
            Message::def(pythonName, a1, a2, a3); \
        F1 ::def(pythonName, a1, a2, a3); \
        F2 ::def(pythonName, a1, a2, a3); \
        F3 ::def(pythonName, a1, a2, a3); \
        F4 ::def(pythonName, a1, a2, a3); \
        F5 ::def(pythonName, a1, a2, a3); \
        F6 ::def(pythonName, a1, a2, a3); \
        F7 ::def(pythonName, a1, a2, a3); \
        F8 ::def(pythonName, a1, a2, a3); \
        F9 ::def(pythonName, a1, a2, a3); \
        F10::def(pythonName, a1, a2, a3); \
        F11::def(pythonName, a1, a2, a3); \
        F12::def(pythonName, a1, a2, a3); \
    } \
};

struct PythonMultidefFunctor {};

template <class T1,
          class T2 = void,
          class T3 = void,
          class T4 = void,
          class T5 = void,
          class T6 = void,
          class T7 = void,
          class T8 = void,
          class T9 = void,
          class T10 = void,
          class T11 = void,
          class T12 = void>
struct ArgumentMismatchMessage
{
    static std::string message()
    {
        std::string res(
            "No C++ overload matches the arguments. This can have three reasons:\n\n"
            " * The array arguments may have an unsupported element type. You may need\n"
            "   to convert your array(s) to another element type using 'array.astype(...)'.\n"
            "   The function currently supports the following types:\n\n     ");
        res += vigra::detail::TypeName<T1>::sized_name();

        if(vigra::detail::TypeName<T2>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T2>::sized_name();
        if(vigra::detail::TypeName<T3>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T3>::sized_name();
        if(vigra::detail::TypeName<T4>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T4>::sized_name();
        if(vigra::detail::TypeName<T5>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T5>::sized_name();
        if(vigra::detail::TypeName<T6>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T6>::sized_name();
        if(vigra::detail::TypeName<T7>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T7>::sized_name();
        if(vigra::detail::TypeName<T8>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T8>::sized_name();
        if(vigra::detail::TypeName<T9>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T9>::sized_name();
        if(vigra::detail::TypeName<T10>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T10>::sized_name();
        if(vigra::detail::TypeName<T11>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T11>::sized_name();
        if(vigra::detail::TypeName<T12>::sized_name() != "void")
            res += ", " + vigra::detail::TypeName<T12>::sized_name();

        res +=
            "\n\n"
            " * The dimension of your array(s) is currently unsupported (consult the\n"
            "   function's documentation for information about supported dimensions).\n\n"
            " * You provided an unrecognized argument, or an argument with incorrect type\n"
            "   (consult the documentation for valid function signatures).\n\n"
            "Additional overloads can easily be added in the vigranumpy C++ sources.\n"
            "Please submit an issue at http://github.com/ukoethe/vigra/ to let us know\n"
            "what you need (or a pull request if you solved it on your own :-).\n\n";

        return res;
    }

    static void def(const char * pythonName)
    {
        std::string msg    = message(),
                    module = extract<std::string>(scope().attr("__name__"))() + ".";
        msg += "Type 'help(" + module + pythonName + ")' to get full documentation.\n";
        boost::python::def(pythonName,
            raw_function([msg](tuple, dict) -> object {
                throw std::invalid_argument(msg);
                return object();
            }, 0));
    }

    template <class A1>
    static void def(const char * pythonName, A1 const &)
    {
        def(pythonName);
    }

    template <class A1, class A2>
    static void def(const char * pythonName, A1 const &, A2 const &)
    {
        def(pythonName);
    }

    template <class A1, class A2, class A3>
    static void def(const char * pythonName, A1 const &, A2 const &, A3 const &)
    {
        def(pythonName);
    }
};

// in the sequel, the doc string is only registered with the last
// overload, so that it shows up only once
template <class Functor>
inline typename std::enable_if<std::is_base_of<PythonMultidefFunctor, Functor>::value,
                               void>::type
multidef(char const* functor_name, Functor const & f)
{
    f.def(functor_name);
}

template <class Functor, class Args>
inline typename std::enable_if<std::is_base_of<PythonMultidefFunctor, Functor>::value,
                               void>::type
multidef(char const* functor_name, Functor const & f, Args const& args)
{
    f.def(functor_name, args);
}

template <class Functor>
inline typename std::enable_if<std::is_base_of<PythonMultidefFunctor, Functor>::value,
                               void>::type
multidef(char const* functor_name, Functor const & f, const char * help)
{
    f.def(functor_name, help);
}

template <class Functor, class Args>
inline typename std::enable_if<std::is_base_of<PythonMultidefFunctor, Functor>::value,
                               void>::type
multidef(char const* functor_name, Functor const & f, Args const& args, const char * help)
{
    f.def(functor_name, args, help);
}

}} // namespace boost::python

#endif // VIGRA_NUMPY_ARRAY_CONVERTERS_HXX
