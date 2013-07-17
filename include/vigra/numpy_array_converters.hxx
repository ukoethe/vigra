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
#include <set>

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


} // namespace vigra

namespace boost { namespace python {

#define VIGRA_PYTHON_MULTITYPE_FUNCTOR(functor_name, function) \
template <class T> \
struct functor_name##Impl \
{ \
    typedef functor_name##Impl type; \
     \
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
    typedef void type; \
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
: public boost::python::TypeList<typename functor_name##Impl<T1>::type, \
         boost::python::TypeList<typename functor_name##Impl<T2>::type, \
         boost::python::TypeList<typename functor_name##Impl<T3>::type, \
         boost::python::TypeList<typename functor_name##Impl<T4>::type, \
         boost::python::TypeList<typename functor_name##Impl<T5>::type, \
         boost::python::TypeList<typename functor_name##Impl<T6>::type, \
         boost::python::TypeList<typename functor_name##Impl<T7>::type, \
         boost::python::TypeList<typename functor_name##Impl<T8>::type, \
         boost::python::TypeList<typename functor_name##Impl<T9>::type, \
         boost::python::TypeList<typename functor_name##Impl<T10>::type, \
         boost::python::TypeList<typename functor_name##Impl<T11>::type, \
         boost::python::TypeList<typename functor_name##Impl<T12>::type, \
         boost::python::TypeList<void, void> > > > > > > > > > > > > \
{};

template <class Head, class Tail>
struct TypeList
{
    typedef Head head;
    typedef Tail tail;
};

// in the sequel, the doc string is only registered with the last
// overload, so that it shows up only once
template <class Head, class Tail>
inline void multidef(char const* functor_name, TypeList<Head, Tail>)
{
    Head::def(functor_name);
    multidef(functor_name, Tail());
}

template <class Head, class Tail>
inline void multidef(char const* functor_name, TypeList<Head, Tail>, const char * help)
{
    Head::def(functor_name);
    multidef(functor_name, Tail(), help);
}

template <class Head, class Tail, class Args>
inline void multidef(char const* functor_name, TypeList<Head, Tail>, Args const& args)
{
    Head::def(functor_name, args);
    multidef(functor_name, Tail(), args);
}

template <class Head, class Tail, class Args>
inline void multidef(char const* functor_name, TypeList<Head, Tail>, Args const& args, char const * help)
{
    Head::def(functor_name, args);
    multidef(functor_name, Tail(), args, help);
}

template <class Head, class Tail>
inline void multidef(char const* functor_name, TypeList<Head, TypeList<void, Tail> >)
{
    Head::def(functor_name);
}

template <class Head, class Tail, class Args>
inline void multidef(char const* functor_name, TypeList<Head, TypeList<void, Tail> >, Args const& args)
{
    Head::def(functor_name, args);
}

template <class Head, class Tail>
inline void multidef(char const* functor_name, TypeList<Head, TypeList<void, Tail> >, const char * help)
{
    Head::def(functor_name, help);
}

template <class Head, class Tail, class Args>
inline void multidef(char const* functor_name, TypeList<Head, TypeList<void, Tail> >, Args const& args, const char * help)
{
    Head::def(functor_name, args, help);
}

}} // namespace boost::python

#endif // VIGRA_NUMPY_ARRAY_CONVERTERS_HXX
