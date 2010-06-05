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

#ifndef VIGRA_NUMPY_ARRAY_HXX
#define VIGRA_NUMPY_ARRAY_HXX

#include <Python.h>
#include <iostream>
#include <algorithm>
#include <complex>
#include <string>
#include <sstream>
#include <map>
#include <vigra/multi_array.hxx>
#include <vigra/array_vector.hxx>
#include <vigra/sized_int.hxx>
#include <vigra/python_utility.hxx>
#include <numpy/arrayobject.h>

int _import_array();

namespace vigra {

/********************************************************/
/*                                                      */
/*              Singleband and Multiband                */
/*                                                      */
/********************************************************/

typedef float NumpyValueType;

template <class T>
struct Singleband  // the last array dimension is not to be interpreted as a channel dimension
{
    typedef T value_type;
};

template <class T>
struct Multiband  // the last array dimension is a channel dimension
{
    typedef T value_type;
};

template<class T>
struct NumericTraits<Singleband<T> >
: public NumericTraits<T>
{};

template<class T>
struct NumericTraits<Multiband<T> >
{
    typedef Multiband<T> Type;
/*
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
*/
    typedef Type ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef typename NumericTraits<T>::isSigned isOrdered;
    typedef typename NumericTraits<T>::isSigned isComplex;
/*
    static signed char zero() { return 0; }
    static signed char one() { return 1; }
    static signed char nonZero() { return 1; }
    static signed char min() { return SCHAR_MIN; }
    static signed char max() { return SCHAR_MAX; }

#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = SCHAR_MIN, maxConst = SCHAR_MIN };
#else
    static const signed char minConst = SCHAR_MIN;
    static const signed char maxConst = SCHAR_MIN;
#endif

    static Promote toPromote(signed char v) { return v; }
    static RealPromote toRealPromote(signed char v) { return v; }
    static signed char fromPromote(Promote v) {
        return ((v < SCHAR_MIN) ? SCHAR_MIN : (v > SCHAR_MAX) ? SCHAR_MAX : v);
    }
    static signed char fromRealPromote(RealPromote v) {
        return ((v < 0.0)
                   ? ((v < (RealPromote)SCHAR_MIN)
                       ? SCHAR_MIN
                       : static_cast<signed char>(v - 0.5))
                   : (v > (RealPromote)SCHAR_MAX)
                       ? SCHAR_MAX
                       : static_cast<signed char>(v + 0.5));
    }
*/
};

template <class T>
class MultibandVectorAccessor
{
    MultiArrayIndex size_, stride_;

  public:
    MultibandVectorAccessor(MultiArrayIndex size, MultiArrayIndex stride)
    : size_(size),
      stride_(stride)
    {}


    typedef Multiband<T> value_type;

        /** the vector's value_type
        */
    typedef T component_type;

    typedef VectorElementAccessor<MultibandVectorAccessor<T> > ElementAccessor;

        /** Read the component data at given vector index
            at given iterator position
        */
    template <class ITERATOR>
    component_type const & getComponent(ITERATOR const & i, int idx) const
    {
        return *(&*i+idx*stride_);
    }

        /** Set the component data at given vector index
            at given iterator position. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
            In case of a conversion floating point -> intergral this includes rounding and clipping.
        */
    template <class V, class ITERATOR>
    void setComponent(V const & value, ITERATOR const & i, int idx) const
    {
        *(&*i+idx*stride_) = detail::RequiresExplicitCast<component_type>::cast(value);
    }

        /** Read the component data at given vector index
            at an offset of given iterator position
        */
    template <class ITERATOR, class DIFFERENCE>
    component_type const & getComponent(ITERATOR const & i, DIFFERENCE const & diff, int idx) const
    {
        return *(&i[diff]+idx*stride_);
    }

    /** Set the component data at given vector index
        at an offset of given iterator position. The type <TT>V</TT> of the passed
        in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
            In case of a conversion floating point -> intergral this includes rounding and clipping.
    */
    template <class V, class ITERATOR, class DIFFERENCE>
    void
    setComponent(V const & value, ITERATOR const & i, DIFFERENCE const & diff, int idx) const
    {
        *(&i[diff]+idx*stride_) = detail::RequiresExplicitCast<component_type>::cast(value);
    }

    template <class U>
    MultiArrayIndex size(U) const
    {
        return size_;
    }
};

/********************************************************/
/*                                                      */
/*                a few Python utilities                */
/*                                                      */
/********************************************************/

namespace detail {

inline long spatialDimensions(PyObject * obj)
{
    static python_ptr key(PyString_FromString("spatialDimensions"), python_ptr::keep_count);
    python_ptr pres(PyObject_GetAttr(obj, key), python_ptr::keep_count);
    long res = pres && PyInt_Check(pres)
                 ? PyInt_AsLong(pres)
                 : -1;
    return res;
}

/*
 * The registry is used to optionally map specific C++ types to
 * specific python sub-classes of numpy.ndarray (for example,
 * MultiArray<2, Singleband<int> > to a user-defined Python class 'ScalarImage').
 *
 * One needs to use NUMPY_ARRAY_INITIALIZE_REGISTRY once in a python
 * extension module using this technique, in order to actually provide
 * the registry (this is done by vigranumpycmodule and will then be
 * available for other modules, too).  Alternatively,
 * NUMPY_ARRAY_DUMMY_REGISTRY may be used to disable this feature
 * completely.  In both cases, the macro must not be enclosed by any
 * namespace, so it is best put right at the beginning of the file
 * (e.g. below the #includes).
 */

typedef std::map<std::string, std::pair<python_ptr, python_ptr> > ArrayTypeMap;

VIGRA_EXPORT ArrayTypeMap * getArrayTypeMap();

#define NUMPY_ARRAY_INITIALIZE_REGISTRY                                 \
    namespace vigra { namespace detail {                                \
    ArrayTypeMap * getArrayTypeMap()                                    \
    {                                                                   \
        static ArrayTypeMap arrayTypeMap;                               \
        return &arrayTypeMap;                                           \
    }                                                                   \
    }} // namespace vigra::detail

#define NUMPY_ARRAY_DUMMY_REGISTRY                      \
    namespace vigra { namespace detail {                \
    ArrayTypeMap * getArrayTypeMap()                    \
    {                                                   \
        return NULL;                                    \
    }                                                   \
    }} // namespace vigra::detail

inline
void registerPythonArrayType(std::string const & name, PyObject * obj, PyObject * typecheck)
{
    ArrayTypeMap *types = getArrayTypeMap();
    vigra_precondition(
        types != NULL,
        "registerPythonArrayType(): module was compiled without array type registry.");
    vigra_precondition(
        obj && PyType_Check(obj) && PyType_IsSubtype((PyTypeObject *)obj, &PyArray_Type),
        "registerPythonArrayType(obj): obj is not a subtype of numpy.ndarray.");
    if(typecheck && PyCallable_Check(typecheck))
        (*types)[name] = std::make_pair(python_ptr(obj), python_ptr(typecheck));
    else
        (*types)[name] = std::make_pair(python_ptr(obj), python_ptr());
//    std::cerr << "Registering " << ((PyTypeObject *)obj)->tp_name << " for " << name << "\n";
}

inline
python_ptr getArrayTypeObject(std::string const & name, PyTypeObject * def = 0)
{
    ArrayTypeMap *types = getArrayTypeMap();
    if(!types)
        // dummy registry -> handle like empty registry
        return python_ptr((PyObject *)def);

    python_ptr res;
    ArrayTypeMap::iterator i = types->find(name);
    if(i != types->end())
        res = i->second.first;
    else
        res = python_ptr((PyObject *)def);
//    std::cerr << "Requested " << name << ", got " << ((PyTypeObject *)res.get())->tp_name << "\n";
    return res;
}

// there are two cases for the return:
// * if a typecheck function was registered, it is returned
// * a null pointer is returned if nothing was registered for either key, or if
//   a type was registered without typecheck function
inline python_ptr
getArrayTypecheckFunction(std::string const & keyFull, std::string const & key)
{
    python_ptr res;
    ArrayTypeMap *types = getArrayTypeMap();
    if(types)
    {
        ArrayTypeMap::iterator i = types->find(keyFull);
        if(i == types->end())
            i = types->find(key);
        if(i != types->end())
            res = i->second.second;
    }
    return res;
}

inline bool
performCustomizedArrayTypecheck(PyObject * obj, std::string const & keyFull, std::string const & key)
{
    if(obj == 0 || !PyArray_Check(obj))
        return false;
    python_ptr typecheck = getArrayTypecheckFunction(keyFull, key);
    if(typecheck == 0)
        return true; // no custom test registered
    python_ptr args(PyTuple_Pack(1, obj), python_ptr::keep_count);
    pythonToCppException(args);
    python_ptr res(PyObject_Call(typecheck.get(), args.get(), 0), python_ptr::keep_count);
    pythonToCppException(res);
    vigra_precondition(PyBool_Check(res),
           "NumpyArray conversion: registered typecheck function did not return a boolean.");
    return (void*)res.get() == (void*)Py_True;
}

inline
python_ptr constructNumpyArrayImpl(
    PyTypeObject * type,
    ArrayVector<npy_intp> const & shape, npy_intp *strides,
    NPY_TYPES typeCode, bool init)
{
    python_ptr array;

    if(strides == 0)
    {
        array = python_ptr(PyArray_New(type, shape.size(), (npy_intp *)shape.begin(), typeCode, 0, 0, 0, 1 /* Fortran order */, 0),
                           python_ptr::keep_count);
    }
    else
    {
        int N = shape.size();
        ArrayVector<npy_intp> pshape(N);
        for(int k=0; k<N; ++k)
            pshape[strides[k]] = shape[k];

        array = python_ptr(PyArray_New(type, N, pshape.begin(), typeCode, 0, 0, 0, 1 /* Fortran order */, 0),
                           python_ptr::keep_count);
        pythonToCppException(array);

        PyArray_Dims permute = { strides, N };
        array = python_ptr(PyArray_Transpose((PyArrayObject*)array.get(), &permute), python_ptr::keep_count);
    }
    pythonToCppException(array);

    if(init)
        PyArray_FILLWBYTE((PyArrayObject *)array.get(), 0);

    return array;
}

// strideOrdering will be ignored unless order == "A"
// TODO: this function should receive some refactoring in order to make
//       the rules clear from the code rather than from comments
inline python_ptr
constructNumpyArrayImpl(PyTypeObject * type, ArrayVector<npy_intp> const & shape,
                       unsigned int spatialDimensions, unsigned int channels,
                       NPY_TYPES typeCode, std::string order, bool init,
                       ArrayVector<npy_intp> strideOrdering = ArrayVector<npy_intp>())
{
    // shape must have at least length spatialDimensions, but can also have a channel dimension
    vigra_precondition(shape.size() == spatialDimensions || shape.size() == spatialDimensions + 1,
           "constructNumpyArray(type, shape, ...): shape has wrong length.");

    // if strideOrdering is given, it must have at least length spatialDimensions,
    // but can also have a channel dimension
    vigra_precondition(strideOrdering.size() == 0 || strideOrdering.size() == spatialDimensions ||
                       strideOrdering.size() == spatialDimensions + 1,
           "constructNumpyArray(type, ..., strideOrdering): strideOrdering has wrong length.");

    if(channels == 0) // if the requested number of channels is not given ...
    {
        // ... deduce it
        if(shape.size() == spatialDimensions)
            channels = 1;
        else
            channels = shape.back();
    }
    else
    {
        // otherwise, if the shape object also contains a channel dimension, they must be consistent
        if(shape.size() > spatialDimensions)
            vigra_precondition(channels == (unsigned int)shape[spatialDimensions],
                   "constructNumpyArray(type, ...): shape contradicts requested number of channels.");
    }

    // if we have only one channel, no explicit channel dimension should be in the shape
    unsigned int shapeSize = channels == 1
                                  ? spatialDimensions
                                  : spatialDimensions + 1;

    // create the shape object with optional channel dimension
    ArrayVector<npy_intp> pshape(shapeSize);
    std::copy(shape.begin(), shape.begin()+std::min(shape.size(), pshape.size()), pshape.begin());
    if(shapeSize > spatialDimensions)
        pshape[spatialDimensions] = channels;

    // order "A" means "preserve order" when an array is copied, and
    // defaults to "V" when a new array is created without explicit strideOrdering
    //
    if(order == "A")
    {
        if(strideOrdering.size() == 0)
        {
            order = "V";
        }
        else if(strideOrdering.size() > shapeSize)
        {
            // make sure that strideOrdering length matches shape length
            ArrayVector<npy_intp> pstride(strideOrdering.begin(), strideOrdering.begin()+shapeSize);

            // adjust the ordering when the channel dimension has been dropped because channel == 1
            if(strideOrdering[shapeSize] == 0)
                for(unsigned int k=0; k<shapeSize; ++k)
                    pstride[k] -= 1;
            pstride.swap(strideOrdering);
        }
        else if(strideOrdering.size() < shapeSize)
        {
            // make sure that strideOrdering length matches shape length
            ArrayVector<npy_intp> pstride(shapeSize);

            // adjust the ordering when the channel dimension has been dropped because channel == 1
            for(unsigned int k=0; k<shapeSize-1; ++k)
                pstride[k] = strideOrdering[k] + 1;
            pstride[shapeSize-1] = 0;
            pstride.swap(strideOrdering);
        }
    }

    // create the appropriate strideOrdering objects for the other memory orders
    // (when strideOrdering already contained data, it is ignored because order != "A")
    if(order == "C")
    {
        strideOrdering.resize(shapeSize);
        for(unsigned int k=0; k<shapeSize; ++k)
            strideOrdering[k] = shapeSize-1-k;
    }
    else if(order == "F" || (order == "V" && channels == 1))
    {
        strideOrdering.resize(shapeSize);
        for(unsigned int k=0; k<shapeSize; ++k)
            strideOrdering[k] = k;
    }
    else if(order == "V")
    {
        strideOrdering.resize(shapeSize);
        for(unsigned int k=0; k<shapeSize-1; ++k)
            strideOrdering[k] = k+1;
        strideOrdering[shapeSize-1] = 0;
    }

    return constructNumpyArrayImpl(type, pshape, strideOrdering.begin(), typeCode, init);
}

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


} // namespace detail

/********************************************************/
/*                                                      */
/*               NumpyArrayValuetypeTraits              */
/*                                                      */
/********************************************************/

template<class ValueType>
struct ERROR_NumpyArrayValuetypeTraits_not_specialized_for_ { };

template<class ValueType>
struct NumpyArrayValuetypeTraits
{
    static bool isValuetypeCompatible(PyArrayObject const * obj)
    {
        return ERROR_NumpyArrayValuetypeTraits_not_specialized_for_<ValueType>();
    }

    static ERROR_NumpyArrayValuetypeTraits_not_specialized_for_<ValueType> typeCode;

    static std::string typeName()
    {
        return std::string("ERROR: NumpyArrayValuetypeTraits not specialized for this case");
    }

    static std::string typeNameImpex()
    {
        return std::string("ERROR: NumpyArrayValuetypeTraits not specialized for this case");
    }

    static PyObject * typeObject()
    {
        return (PyObject *)0;
    }
};

template<class ValueType>
ERROR_NumpyArrayValuetypeTraits_not_specialized_for_<ValueType> NumpyArrayValuetypeTraits<ValueType>::typeCode;

#define VIGRA_NUMPY_VALUETYPE_TRAITS(type, typeID, numpyTypeName, impexTypeName) \
template <> \
struct NumpyArrayValuetypeTraits<type > \
{ \
    static bool isValuetypeCompatible(PyArrayObject const * obj) /* obj must not be NULL */ \
    { \
        return PyArray_EquivTypenums(typeID, PyArray_DESCR((PyObject *)obj)->type_num) && \
               PyArray_ITEMSIZE((PyObject *)obj) == sizeof(type); \
    } \
    \
    static NPY_TYPES const typeCode = typeID; \
    \
    static std::string typeName() \
    { \
        return #numpyTypeName; \
    } \
    \
    static std::string typeNameImpex() \
    { \
        return impexTypeName; \
    } \
    \
    static PyObject * typeObject() \
    { \
        return PyArray_TypeObjectFromType(typeID); \
    } \
};

VIGRA_NUMPY_VALUETYPE_TRAITS(bool,           NPY_BOOL, bool, "UINT8")
VIGRA_NUMPY_VALUETYPE_TRAITS(signed char,    NPY_INT8, int8, "INT16")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned char,  NPY_UINT8, uint8, "UINT8")
VIGRA_NUMPY_VALUETYPE_TRAITS(short,          NPY_INT16, int16, "INT16")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned short, NPY_UINT16, uint16, "UINT16")

#if VIGRA_BITSOF_LONG == 32
VIGRA_NUMPY_VALUETYPE_TRAITS(long,           NPY_INT32, int32, "INT32")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned long,  NPY_UINT32, uint32, "UINT32")
#elif VIGRA_BITSOF_LONG == 64
VIGRA_NUMPY_VALUETYPE_TRAITS(long,           NPY_INT64, int64, "DOUBLE")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned long,  NPY_UINT64, uint64, "DOUBLE")
#endif

#if VIGRA_BITSOF_INT == 32
VIGRA_NUMPY_VALUETYPE_TRAITS(int,            NPY_INT32, int32, "INT32")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned int,   NPY_UINT32, uint32, "UINT32")
#elif VIGRA_BITSOF_INT == 64
VIGRA_NUMPY_VALUETYPE_TRAITS(int,            NPY_INT64, int64, "DOUBLE")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned int,   NPY_UINT64, uint64, "DOUBLE")
#endif

#ifdef PY_LONG_LONG
# if VIGRA_BITSOF_LONG_LONG == 32
VIGRA_NUMPY_VALUETYPE_TRAITS(long long,            NPY_INT32, int32, "INT32")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned long long,   NPY_UINT32, uint32, "UINT32")
# elif VIGRA_BITSOF_LONG_LONG == 64
VIGRA_NUMPY_VALUETYPE_TRAITS(long long,          NPY_INT64, int64, "DOUBLE")
VIGRA_NUMPY_VALUETYPE_TRAITS(unsigned long long, NPY_UINT64, uint64, "DOUBLE")
# endif
#endif

VIGRA_NUMPY_VALUETYPE_TRAITS(npy_float32, NPY_FLOAT32, float32, "FLOAT")
VIGRA_NUMPY_VALUETYPE_TRAITS(npy_float64, NPY_FLOAT64, float64, "DOUBLE")
#if NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE
VIGRA_NUMPY_VALUETYPE_TRAITS(npy_longdouble, NPY_LONGDOUBLE, longdouble, "")
#endif
VIGRA_NUMPY_VALUETYPE_TRAITS(npy_cfloat, NPY_CFLOAT, complex64, "")
VIGRA_NUMPY_VALUETYPE_TRAITS(std::complex<npy_float>, NPY_CFLOAT, complex64, "")
VIGRA_NUMPY_VALUETYPE_TRAITS(npy_cdouble, NPY_CDOUBLE, complex128, "")
VIGRA_NUMPY_VALUETYPE_TRAITS(std::complex<npy_double>, NPY_CDOUBLE, complex128, "")
VIGRA_NUMPY_VALUETYPE_TRAITS(npy_clongdouble, NPY_CLONGDOUBLE, clongdouble, "")
#if NPY_SIZEOF_LONGDOUBLE != NPY_SIZEOF_DOUBLE
VIGRA_NUMPY_VALUETYPE_TRAITS(std::complex<npy_longdouble>, NPY_CLONGDOUBLE, clongdouble, "")
#endif

#undef VIGRA_NUMPY_VALUETYPE_TRAITS

/********************************************************/
/*                                                      */
/*                  NumpyArrayTraits                    */
/*                                                      */
/********************************************************/

template <class U, int N>
bool stridesAreAscending(TinyVector<U, N> const & strides)
{
    for(int k=1; k<N; ++k)
        if(strides[k] < strides[k-1])
            return false;
    return true;
}

template<unsigned int N, class T, class Stride>
struct NumpyArrayTraits;

template<unsigned int N, class T>
struct NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef T dtype;
    typedef T value_type;
    typedef NumpyArrayValuetypeTraits<T> ValuetypeTraits;
    static NPY_TYPES const typeCode = ValuetypeTraits::typeCode;

    enum { spatialDimensions = N, channels = 1 };

    static bool isArray(PyObject * obj)
    {
        return obj && PyArray_Check(obj);
    }

    static bool isClassCompatible(PyObject * obj)
    {
        return detail::performCustomizedArrayTypecheck(obj, typeKeyFull(), typeKey());
    }

    static bool isValuetypeCompatible(PyArrayObject * obj)  /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj);
    }

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return PyArray_NDIM((PyObject *)obj) == N-1 ||
               PyArray_NDIM((PyObject *)obj) == N ||
               (PyArray_NDIM((PyObject *)obj) == N+1 && PyArray_DIM((PyObject *)obj, N) == 1);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKey()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", *>";
        return key;
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", " +
                                 ValuetypeTraits::typeName() + ", StridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, T, UnstridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isShapeCompatible(obj) &&
               PyArray_STRIDES((PyObject *)obj)[0] == PyArray_ITEMSIZE((PyObject *)obj);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", " +
                                 ValuetypeTraits::typeName() + ", UnstridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Singleband<T>, StridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isClassCompatible(PyObject * obj)
    {
        return detail::performCustomizedArrayTypecheck(obj, typeKeyFull(), typeKey());
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKey()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Singleband<*> >";
        return key;
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Singleband<" +
                                 ValuetypeTraits::typeName() + ">, StridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Singleband<T>, UnstridedArrayTag>
: public NumpyArrayTraits<N, Singleband<T>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, UnstridedArrayTag> UnstridedTraits;
    typedef NumpyArrayTraits<N, Singleband<T>, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return UnstridedTraits::isShapeCompatible(obj);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return UnstridedTraits::isPropertyCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Singleband<" +
                                 ValuetypeTraits::typeName() + ">, UnstridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Multiband<T>, StridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    enum { spatialDimensions = N-1, channels = 0 };

    static bool isClassCompatible(PyObject * obj)
    {
        return detail::performCustomizedArrayTypecheck(obj, typeKeyFull(), typeKey());
    }

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return PyArray_NDIM(obj) == N || PyArray_NDIM(obj) == N-1;
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKey()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Multiband<*> >";
        return key;
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Multiband<" +
                                 ValuetypeTraits::typeName() + ">, StridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Multiband<T>, UnstridedArrayTag>
: public NumpyArrayTraits<N, Multiband<T>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, Multiband<T>, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isShapeCompatible(obj) &&
               PyArray_STRIDES((PyObject *)obj)[0] == PyArray_ITEMSIZE((PyObject *)obj);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", Multiband<" +
                                 ValuetypeTraits::typeName() + ">, UnstridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, int M, class T>
struct NumpyArrayTraits<N, TinyVector<T, M>, StridedArrayTag>
{
    typedef T dtype;
    typedef TinyVector<T, M> value_type;
    typedef NumpyArrayValuetypeTraits<T> ValuetypeTraits;
    static NPY_TYPES const typeCode = ValuetypeTraits::typeCode;

    enum { spatialDimensions = N, channels = M };

    static bool isArray(PyObject * obj)
    {
        return obj && PyArray_Check(obj);
    }

    static bool isClassCompatible(PyObject * obj)
    {
        return detail::performCustomizedArrayTypecheck(obj, typeKeyFull(), typeKey());
    }

    static bool isValuetypeCompatible(PyArrayObject * obj)  /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj);
    }

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return PyArray_NDIM((PyObject *)obj) == N+1 &&
               PyArray_DIM((PyObject *)obj, N) == M &&
               PyArray_STRIDES((PyObject *)obj)[N] == PyArray_ITEMSIZE((PyObject *)obj);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N+1> npyShape;
        std::copy(shape.begin(), shape.end(), npyShape.begin());
        npyShape[N] = M;

        TinyVector<npy_intp, N+1> npyStride;
        std::transform(
            stride.begin(), stride.end(), npyStride.begin(),
            std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        npyStride[N] = sizeof(T);

        return detail::constructNumpyArrayFromData(
            typeKeyFull(), typeKey(), npyShape,
            npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKey()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", TinyVector<*, " + asString(M) + "> >";
        return key;
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) +
                      ", TinyVector<" + ValuetypeTraits::typeName() + ", " + asString(M) + ">, StridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, int M, class T>
struct NumpyArrayTraits<N, TinyVector<T, M>, UnstridedArrayTag>
: public NumpyArrayTraits<N, TinyVector<T, M>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, TinyVector<T, M>, StridedArrayTag> BaseType;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isShapeCompatible(obj) &&
               PyArray_STRIDES((PyObject *)obj)[0] == sizeof(TinyVector<T, M>);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return BaseType::isValuetypeCompatible(obj) &&
               isShapeCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N+1> npyShape;
        std::copy(shape.begin(), shape.end(), npyShape.begin());
        npyShape[N] = M;

        TinyVector<npy_intp, N+1> npyStride;
        std::transform(
            stride.begin(), stride.end(), npyStride.begin(),
            std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        npyStride[N] = sizeof(T);

        return detail::constructNumpyArrayFromData(
            typeKeyFull(), BaseType::typeKey(), npyShape,
            npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) +
                      ", TinyVector<" + ValuetypeTraits::typeName() + ", " + asString(M) + ">, UnstridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, RGBValue<T>, StridedArrayTag>
: public NumpyArrayTraits<N, TinyVector<T, 3>, StridedArrayTag>
{
    typedef T dtype;
    typedef RGBValue<T> value_type;
    typedef NumpyArrayValuetypeTraits<T> ValuetypeTraits;

    static bool isClassCompatible(PyObject * obj)
    {
        return detail::performCustomizedArrayTypecheck(obj, typeKeyFull(), typeKey());
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N+1> npyShape;
        std::copy(shape.begin(), shape.end(), npyShape.begin());
        npyShape[N] = 3;

        TinyVector<npy_intp, N+1> npyStride;
        std::transform(
            stride.begin(), stride.end(), npyStride.begin(),
            std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        npyStride[N] = sizeof(T);

        return detail::constructNumpyArrayFromData(
            typeKeyFull(), typeKey(), npyShape,
            npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKey()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) + ", RGBValue<*> >";
        return key;
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) +
                      ", RGBValue<" + ValuetypeTraits::typeName() + ">, StridedArrayTag>";
        return key;
    }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, RGBValue<T>, UnstridedArrayTag>
: public NumpyArrayTraits<N, RGBValue<T>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, TinyVector<T, 3>, UnstridedArrayTag> UnstridedTraits;
    typedef NumpyArrayTraits<N, RGBValue<T>, StridedArrayTag> BaseType;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return UnstridedTraits::isShapeCompatible(obj);
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return UnstridedTraits::isPropertyCompatible(obj);
    }

    template <class U>
    static python_ptr constructor(TinyVector<U, N> const & shape,
                                  T *data, TinyVector<U, N> const & stride)
    {
        TinyVector<npy_intp, N+1> npyShape;
        std::copy(shape.begin(), shape.end(), npyShape.begin());
        npyShape[N] = 3;

        TinyVector<npy_intp, N+1> npyStride;
        std::transform(
            stride.begin(), stride.end(), npyStride.begin(),
            std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        npyStride[N] = sizeof(T);

        return detail::constructNumpyArrayFromData(
            typeKeyFull(), BaseType::typeKey(), npyShape,
            npyStride.begin(), ValuetypeTraits::typeCode, data);
    }

    static std::string typeKeyFull()
    {
        static std::string key = std::string("NumpyArray<") + asString(N) +
                      ", RGBValue<" + ValuetypeTraits::typeName() + ">, UnstridedArrayTag>";
        return key;
    }
};

/********************************************************/
/*                                                      */
/*                    NumpyAnyArray                     */
/*                                                      */
/********************************************************/

/** Wrapper class for a Python array.

    This class stores a reference-counted pointer to an Python numpy array object,
    i.e. an object where <tt>PyArray_Check(object)</tt> returns true (in Python, the
    object is then a subclass of <tt>numpy.ndarray</tt>). This class is mainly used
    as a smart pointer to these arrays, but some basic access and conversion functions
    are also provided.

    <b>\#include</b> \<<a href="numpy__array_8hxx-source.html">vigra/numpy_array.hxx</a>\><br>
    Namespace: vigra
*/
class NumpyAnyArray
{
  protected:
    python_ptr pyArray_;

    // We want to apply broadcasting to the channel dimension.
    // Since only leading dimensions can be added during numpy
    // broadcasting, we permute the array accordingly.
    NumpyAnyArray permuteChannelsToFront() const
    {
        MultiArrayIndex M = ndim();
        ArrayVector<npy_intp> permutation(M);
        for(int k=0; k<M; ++k)
            permutation[k] = M-1-k;
        // explicit cast to int is neede here to avoid gcc c++0x compilation
        // error: narrowing conversion of ‘M’ from ‘vigra::MultiArrayIndex’
        //        to ‘int’ inside { }
        // int overflow should not occur here because PyArray_NDIM returns
        // an integer which is converted to long in NumpyAnyArray::ndim()
        PyArray_Dims permute = { permutation.begin(), (int) M };
        python_ptr array(PyArray_Transpose(pyArray(), &permute), python_ptr::keep_count);
        pythonToCppException(array);
        return NumpyAnyArray(array.ptr());
    }

  public:

        /// difference type
    typedef ArrayVector<npy_intp> difference_type;

        /**
         Construct from a Python object. If \a obj is NULL, or is not a subclass
         of numpy.ndarray, the resulting NumpyAnyArray will have no data (i.e.
         hasData() returns false). Otherwise, it creates a new reference to the array
         \a obj, unless \a createCopy is true, where a new array is created by calling
         the C-equivalent of obj->copy().
         */
    explicit NumpyAnyArray(PyObject * obj = 0, bool createCopy = false, PyTypeObject * type = 0)
    {
        if(obj == 0)
            return;
        vigra_precondition(type == 0 || PyType_IsSubtype(type, &PyArray_Type),
             "NumpyAnyArray(obj, createCopy, type): type must be numpy.ndarray or a subclass thereof.");
        if(createCopy)
            makeCopy(obj, type);
        else
            vigra_precondition(makeReference(obj, type), "NumpyAnyArray(obj): obj isn't a numpy array.");
    }

        /**
         Copy constructor. By default, it creates a new reference to the array
         \a other. When \a createCopy is true, a new array is created by calling
         the C-equivalent of other.copy().
         */
    NumpyAnyArray(NumpyAnyArray const & other, bool createCopy = false, PyTypeObject * type = 0)
    {
        if(!other.hasData())
            return;
        vigra_precondition(type == 0 || PyType_IsSubtype(type, &PyArray_Type),
             "NumpyAnyArray(obj, createCopy, type): type must be numpy.ndarray or a subclass thereof.");
        if(createCopy)
            makeCopy(other.pyObject(), type);
        else
            makeReference(other.pyObject(), type);
    }

    // auto-generated destructor is ok

        /**
         * Assignment operator. If this is already a view with data
         * (i.e. hasData() is true) and the shapes match, the RHS
         * array contents are copied via the C-equivalent of
         * 'self[...] = other[...]'. If the shapes don't matched,
         * broadcasting is tried on the trailing (i.e. channel)
         * dimension.
         * If the LHS is an empty view, assignment is identical to
         * makeReference(other.pyObject()).
         */
    NumpyAnyArray & operator=(NumpyAnyArray const & other)
    {
        if(hasData())
        {
            vigra_precondition(other.hasData(),
                "NumpyArray::operator=(): Cannot assign from empty array.");
            if(PyArray_CopyInto(permuteChannelsToFront().pyArray(), other.permuteChannelsToFront().pyArray()) == -1)
                pythonToCppException(0);
        }
        else
        {
            pyArray_ = other.pyArray_;
        }
        return *this;
    }

        /**
         Returns the number of dimensions of this array, or 0 if
         hasData() is false.
         */
    MultiArrayIndex ndim() const
    {
        if(hasData())
            return PyArray_NDIM(pyObject());
        return 0;
    }

        /**
         Returns the number of spatial dimensions of this array, or 0 if
         hasData() is false. If the enclosed Python array does not define
         the attribute spatialDimensions, ndim() is returned.
         */
    MultiArrayIndex spatialDimensions() const
    {
        if(!hasData())
            return 0;
        MultiArrayIndex s = detail::spatialDimensions(pyObject());
        if(s == -1)
            s = ndim();
        return s;
    }

        /**
         Returns the shape of this array. The size of
         the returned shape equals ndim().
         */
    difference_type shape() const
    {
        if(hasData())
            return difference_type(PyArray_DIMS(pyObject()), PyArray_DIMS(pyObject()) + ndim());
        return difference_type();
    }

        /** Compute the ordering of the strides of this array.
            The result is describes the current permutation of the axes relative
            to an ascending stride order.
        */
    difference_type strideOrdering() const
    {
        if(!hasData())
            return difference_type();
        MultiArrayIndex N = ndim();
        difference_type stride(PyArray_STRIDES(pyObject()), PyArray_STRIDES(pyObject()) + N),
                        permutation(N);
        for(MultiArrayIndex k=0; k<N; ++k)
            permutation[k] = k;
        for(MultiArrayIndex k=0; k<N-1; ++k)
        {
            MultiArrayIndex smallest = k;
            for(MultiArrayIndex j=k+1; j<N; ++j)
            {
                if(stride[j] < stride[smallest])
                    smallest = j;
            }
            if(smallest != k)
            {
                std::swap(stride[k], stride[smallest]);
                std::swap(permutation[k], permutation[smallest]);
            }
        }
        difference_type ordering(N);
        for(MultiArrayIndex k=0; k<N; ++k)
            ordering[permutation[k]] = k;
        return ordering;
    }

        /**
         Returns the value type of the elements in this array, or -1
         when hasData() is false.
         */
    int dtype() const
    {
        if(hasData())
            return PyArray_DESCR(pyObject())->type_num;
        return -1;
    }

        /**
         * Return a borrowed reference to the internal PyArrayObject.
         */
    PyArrayObject * pyArray() const
    {
        return (PyArrayObject *)pyArray_.get();
    }

        /**
         * Return a borrowed reference to the internal PyArrayObject
         * (see pyArray()), cast to PyObject for your convenience.
         */
    PyObject * pyObject() const
    {
        return pyArray_.get();
    }

        /**
           Reset the NumpyAnyArray to the given object. If \a obj is a numpy array object,
           a new reference to that array is created, and the function returns
           true. Otherwise, it returns false and the NumpyAnyArray remains unchanged.
           If \a type is given, the new reference will be a view with that type, provided
           that \a type is a numpy ndarray or a subclass thereof. Otherwise, an
           exception is thrown.
         */
    bool makeReference(PyObject * obj, PyTypeObject * type = 0)
    {
        if(obj == 0 || !PyArray_Check(obj))
            return false;
        if(type != 0)
        {
            vigra_precondition(PyType_IsSubtype(type, &PyArray_Type) != 0,
                "NumpyAnyArray::makeReference(obj, type): type must be numpy.ndarray or a subclass thereof.");
            obj = PyArray_View((PyArrayObject*)obj, 0, type);
            pythonToCppException(obj);
        }
        pyArray_.reset(obj);
        return true;
    }

        /**
           Create a copy of the given array object. If \a obj is a numpy array object,
           a copy is created via the C-equivalent of 'obj->copy()'. If
           this call fails, or obj was not an array, an exception is thrown
           and the NumpyAnyArray remains unchanged.
         */
    void makeCopy(PyObject * obj, PyTypeObject * type = 0)
    {
        vigra_precondition(obj && PyArray_Check(obj),
             "NumpyAnyArray::makeCopy(obj): obj is not an array.");
        vigra_precondition(type == 0 || PyType_IsSubtype(type, &PyArray_Type),
             "NumpyAnyArray::makeCopy(obj, type): type must be numpy.ndarray or a subclass thereof.");
        python_ptr array(PyArray_NewCopy((PyArrayObject*)obj, NPY_ANYORDER), python_ptr::keep_count);
        pythonToCppException(array);
        makeReference(array, type);
    }

         /**
           Check whether this NumpyAnyArray actually points to a Python array.
         */
    bool hasData() const
    {
        return pyArray_ != 0;
    }
};

/********************************************************/
/*                                                      */
/*                     NumpyArray                       */
/*                                                      */
/********************************************************/

/** Provide the MultiArrayView interface for a Python array.

    This class inherits from both \ref vigra::MultiArrayView and \ref vigra::NumpyAnyArray
    in order to support easy and save application of VIGRA functions to Python arrays.

    <b>\#include</b> \<<a href="numpy__array_8hxx-source.html">vigra/numpy_array.hxx</a>\><br>
    Namespace: vigra
*/
template <unsigned int N, class T, class Stride = StridedArrayTag>
class NumpyArray
: public MultiArrayView<N, typename NumpyArrayTraits<N, T, Stride>::value_type, Stride>,
  public NumpyAnyArray
{
  public:
    typedef NumpyArrayTraits<N, T, Stride> ArrayTraits;
    typedef typename ArrayTraits::dtype dtype;
    typedef T pseudo_value_type;

    static NPY_TYPES const typeCode = ArrayTraits::typeCode;

        /** the view type associated with this array.
         */
    typedef MultiArrayView<N, typename ArrayTraits::value_type, Stride> view_type;

    enum { actual_dimension = view_type::actual_dimension };

        /** the array's value type
         */
    typedef typename view_type::value_type value_type;

        /** pointer type
         */
    typedef typename view_type::pointer pointer;

        /** const pointer type
         */
    typedef typename view_type::const_pointer const_pointer;

        /** reference type (result of operator[])
         */
    typedef typename view_type::reference reference;

        /** const reference type (result of operator[] const)
         */
    typedef typename view_type::const_reference const_reference;

        /** size type
         */
    typedef typename view_type::size_type size_type;

        /** difference type (used for multi-dimensional offsets and indices)
         */
    typedef typename view_type::difference_type difference_type;

        /** difference and index type for a single dimension
         */
    typedef typename view_type::difference_type_1 difference_type_1;

        /** traverser type
         */
    typedef typename view_type::traverser traverser;

        /** traverser type to const data
         */
    typedef typename view_type::const_traverser const_traverser;

        /** sequential (random access) iterator type
         */
    typedef value_type * iterator;

        /** sequential (random access) const iterator type
         */
    typedef value_type * const_iterator;

    using view_type::shape;   // resolve ambiguity of multiple inheritance
    using view_type::hasData; // resolve ambiguity of multiple inheritance
    using view_type::strideOrdering; // resolve ambiguity of multiple inheritance

  protected:

    // this function assumes that pyArray_ has already been set, and compatibility been checked
    void setupArrayView();

    static python_ptr getArrayTypeObject()
    {
        python_ptr type = detail::getArrayTypeObject(ArrayTraits::typeKeyFull());
        if(type == 0)
            type = detail::getArrayTypeObject(ArrayTraits::typeKey(), &PyArray_Type);
        return type;
    }

    static python_ptr init(difference_type const & shape, bool init = true)
    {
        ArrayVector<npy_intp> pshape(shape.begin(), shape.end());
        return detail::constructNumpyArrayImpl((PyTypeObject *)getArrayTypeObject().ptr(), pshape,
                       ArrayTraits::spatialDimensions, ArrayTraits::channels,
                       typeCode, "V", init);
    }

    static python_ptr init(difference_type const & shape, difference_type const & strideOrdering, bool init = true)
    {
        ArrayVector<npy_intp> pshape(shape.begin(), shape.end()),
                              pstrideOrdering(strideOrdering.begin(), strideOrdering.end());
        return detail::constructNumpyArrayImpl((PyTypeObject *)getArrayTypeObject().ptr(), pshape,
                       ArrayTraits::spatialDimensions, ArrayTraits::channels,
                       typeCode, "A", init, pstrideOrdering);
    }

  public:

    using view_type::init;

        /**
         * Construct from a given PyObject pointer. When the given
         * python object is NULL, the internal python array will be
         * NULL and hasData() will return false.
         *
         * Otherwise, the function attempts to create a
         * new reference to the given Python object, unless
         * copying is forced by setting \a createCopy to true.
         * If either of this fails, the function throws an exception.
         * This will not happen if isStrictlyCompatible(obj) (in case
         * of creating a new reference) or isCopyCompatible(obj)
         * (in case of copying) have returned true beforehand.
         */
    explicit NumpyArray(PyObject *obj = 0, bool createCopy = false)
    {
        if(obj == 0)
            return;
        if(createCopy)
            makeCopy(obj);
        else
            vigra_precondition(makeReference(obj),
                  "NumpyArray(obj): Cannot construct from incompatible array.");
    }

       /**
         * Copy constructor; does not copy the memory, but creates a
         * new reference to the same underlying python object, unless
         * a copy is forced by setting \a createCopy to true.
         * (If the source object has no data, this one will have
         * no data, too.)
         */
    NumpyArray(const NumpyArray &other, bool createCopy = false) :
            MultiArrayView<N, typename NumpyArrayTraits<N, T, Stride>::value_type, Stride>(other),
            NumpyAnyArray(other, createCopy)
    {
        if(!other.hasData())
            return;
        if(createCopy)
            makeCopy(other.pyObject());
        else
            makeReferenceUnchecked(other.pyObject());
    }

       /**
         * Allocate new memory and copy data from a MultiArrayView.
         */
    explicit NumpyArray(const view_type &other)
    {
        if(!other.hasData())
            return;
        vigra_postcondition(makeReference(init(other.shape(), false)),
                  "NumpyArray(view_type): Python constructor did not produce a compatible array.");
        static_cast<view_type &>(*this) = other;
    }

        /**
         * Construct a new array object, allocating an internal python
         * ndarray of the given shape (in fortran order), initialized
         * with zeros.
         *
         * An exception is thrown when construction fails.
         */
    explicit NumpyArray(difference_type const & shape)
    {
        vigra_postcondition(makeReference(init(shape)),
                     "NumpyArray(shape): Python constructor did not produce a compatible array.");
    }

        /**
         * Construct a new array object, allocating an internal python
         * ndarray of the given shape and given stride ordering, initialized
         * with zeros.
         *
         * An exception is thrown when construction fails.
         */
    NumpyArray(difference_type const & shape, difference_type const & strideOrdering)
    {
        vigra_postcondition(makeReference(init(shape, strideOrdering)),
                     "NumpyArray(shape): Python constructor did not produce a compatible array.");
    }

        /**
         * Constructor from NumpyAnyArray.
         * Equivalent to NumpyArray(other.pyObject())
         */
    NumpyArray(const NumpyAnyArray &other, bool createCopy = false)
    {
        if(!other.hasData())
            return;
        if(createCopy)
            makeCopy(other.pyObject());
        else
            vigra_precondition(makeReference(other.pyObject()), //, false),
                   "NumpyArray(NumpyAnyArray): Cannot construct from incompatible or empty array.");
    }

        /**
         * Assignment operator. If this is already a view with data
         * (i.e. hasData() is true) and the shapes match, the RHS
         * array contents are copied.  If this is an empty view,
         * assignment is identical to makeReferenceUnchecked(other.pyObject()).
         * See MultiArrayView::operator= for further information on
         * semantics.
         */
    NumpyArray &operator=(const NumpyArray &other)
    {
        if(hasData())
            view_type::operator=(other);
        else
            makeReferenceUnchecked(other.pyObject());
        return *this;
    }

        /**
         * Assignment operator. If this is already a view with data
         * (i.e. hasData() is true) and the shapes match, the RHS
         * array contents are copied.
         * If this is an empty view, assignment is identical to
         * makeReference(other.pyObject()).
         * Otherwise, an exception is thrown.
         */
    NumpyArray &operator=(const NumpyAnyArray &other)
    {
        if(hasData())
        {
            NumpyAnyArray::operator=(other);
        }
        else if(isStrictlyCompatible(other.pyObject()))
        {
            makeReferenceUnchecked(other.pyObject());
        }
        else
        {
            vigra_precondition(false,
                "NumpyArray::operator=(): Cannot assign from incompatible array.");
        }
        return *this;
    }

        /**
         * Test whether a given python object is a numpy array that can be
         * converted (copied) into an array compatible to this NumpyArray type.
         * This means that the array's shape conforms to the requirements of
         * makeCopy().
         */
    static bool isCopyCompatible(PyObject *obj)
    {
        return ArrayTraits::isArray(obj) &&
               ArrayTraits::isShapeCompatible((PyArrayObject *)obj);
    }

        /**
         * Test whether a given python object is a numpy array with a
         * compatible dtype and the correct shape and strides, so that it
         * can be referenced as a view by this NumpyArray type (i.e.
         * it conforms to the requirements of makeReference()).
         */
    static bool isReferenceCompatible(PyObject *obj)
    {
        return ArrayTraits::isArray(obj) &&
               ArrayTraits::isPropertyCompatible((PyArrayObject *)obj);
    }

        /**
         * Like isReferenceCompatible(obj), but also executes a customized type compatibility
         * check when such a check has been registered for this class via
         * registerPythonArrayType().
         *
         * This facilitates proper overload resolution between
         * NumpyArray<3, Multiband<T> > (a multiband image) and NumpyArray<3, Singleband<T> > (a scalar volume).
         */
    static bool isStrictlyCompatible(PyObject *obj)
    {
#if VIGRA_CONVERTER_DEBUG
        std::cerr << "class " << typeid(NumpyArray).name() << " got " << obj->ob_type->tp_name << "\n";
        bool isClassCompatible=ArrayTraits::isClassCompatible(obj);
        bool isPropertyCompatible((PyArrayObject *)obj);
        std::cerr<<"isClassCompatible: "<<isClassCompatible<<std::endl;
        std::cerr<<"isPropertyCompatible: "<<isPropertyCompatible<<std::endl;
#endif
        return ArrayTraits::isClassCompatible(obj) &&
               ArrayTraits::isPropertyCompatible((PyArrayObject *)obj);
    }

        /**
         * Create a vector representing the standard stride ordering of a NumpyArray.
         * That is, we get a vector representing the range [0,...,N-1], which
         * denotes the stride ordering for Fortran order.
         */
    static difference_type standardStrideOrdering()
    {
        difference_type strideOrdering;
        for(unsigned int k=0; k<N; ++k)
            strideOrdering[k] = k;
        return strideOrdering;
    }

        /**
         * Set up a view to the given object without checking compatibility.
         * This function must not be used unless isReferenceCompatible(obj) returned
         * true on the given object (otherwise, a crash is likely).
         */
    void makeReferenceUnchecked(PyObject *obj)
    {
        NumpyAnyArray::makeReference(obj);
        setupArrayView();
    }

        /**
         * Try to set up a view referencing the given PyObject.
         * Returns false if the python object is not a compatible
         * numpy array (see isReferenceCompatible() or
         * isStrictlyCompatible(), according to the parameter \a
         * strict).
         */
    bool makeReference(PyObject *obj, bool strict = true)
    {
        if(strict)
        {
            if(!isStrictlyCompatible(obj))
                return false;
        }
        else
        {
            if(!isReferenceCompatible(obj))
                return false;
        }
        makeReferenceUnchecked(obj);
        return true;
    }

        /**
         * Try to set up a view referencing the same data as the given
         * NumpyAnyArray.  This overloaded variant simply calls
         * makeReference() on array.pyObject().
         */
    bool makeReference(const NumpyAnyArray &array, bool strict = true)
    {
        return makeReference(array.pyObject(), strict);
    }

        /**
         * Set up an unsafe reference to the given MultiArrayView.
         * ATTENTION: This creates a numpy.ndarray that points to the
         * same data, but does not own it, so it must be ensured by
         * other means that the memory does not get freed before the
         * end of the ndarray's lifetime!  (One elegant way would be
         * to set the 'base' attribute of the resulting ndarray to a
         * python object which directly or indirectly holds the memory
         * of the given MultiArrayView.)
         */
    void makeReference(const view_type &multiArrayView)
    {
        vigra_precondition(!hasData(), "makeReference(): cannot replace existing view with given buffer");

        // construct an ndarray that points to our data (taking strides into account):
        python_ptr array(ArrayTraits::constructor(multiArrayView.shape(), multiArrayView.data(), multiArrayView.stride()));

        view_type::operator=(multiArrayView);
        pyArray_ = array;
    }

        /**
         Try to create a copy of the given PyObject.
         Raises an exception when obj is not a compatible array
         (see isCopyCompatible() or isStrictlyCompatible(), according to the
         parameter \a strict) or the Python constructor call failed.
         */
    void makeCopy(PyObject *obj, bool strict = false)
    {
        vigra_precondition(strict ? isStrictlyCompatible(obj) : isCopyCompatible(obj),
                     "NumpyArray::makeCopy(obj): Cannot copy an incompatible array.");

        int M = PyArray_NDIM(obj);
        TinyVector<npy_intp, N> shape;
        std::copy(PyArray_DIMS(obj), PyArray_DIMS(obj)+M, shape.begin());
        if(M == N-1)
            shape[M] = 1;
        vigra_postcondition(makeReference(init(shape, false)),
                     "NumpyArray::makeCopy(obj): Copy created an incompatible array.");
        NumpyAnyArray::operator=(NumpyAnyArray(obj));
//        if(PyArray_CopyInto(pyArray(), (PyArrayObject*)obj) == -1)
//            pythonToCppException(0);
    }

        /**
            Allocate new memory with the given shape and initialize with zeros.<br>
            If a stride ordering is given, the resulting array will have this stride
            ordering, when it is compatible with the array's memory layout (unstrided
            arrays only permit the standard ascending stride ordering).

            <em>Note:</em> this operation invalidates dependent objects
            (MultiArrayViews and iterators)
         */
    void reshape(difference_type const & shape, difference_type const & strideOrdering = standardStrideOrdering())
    {
        vigra_postcondition(makeReference(init(shape, strideOrdering)),
                     "NumpyArray(shape): Python constructor did not produce a compatible array.");
    }

        /**
            When this array has no data, allocate new memory with the given \a shape and
            initialize with zeros. Otherwise, check if the new shape matches the old shape
            and throw a precondition exception with the given \a message if not.
         */
    void reshapeIfEmpty(difference_type const & shape, std::string message = "")
    {
        reshapeIfEmpty(shape, standardStrideOrdering(), message);
    }

        /**
            When this array has no data, allocate new memory with the given \a shape and
            initialize with zeros. Otherwise, check if the new shape matches the old shape
            and throw a precondition exception with the given \a message if not. If strict
            is true, the given stride ordering must also match that of the existing data.
         */
    void reshapeIfEmpty(difference_type const & shape, difference_type const & strideOrdering,
                        std::string message = "", bool strict = false)
    {
        if(hasData())
        {
            if(strict)
            {
                if(message == "")
                    message = "NumpyArray::reshapeIfEmpty(shape): array was not empty, and shape or stride ordering did not match.";
                vigra_precondition(shape == this->shape() && strideOrdering == this->strideOrdering(), message.c_str());
            }
            else
            {
                if(message == "")
                    message = "NumpyArray::reshapeIfEmpty(shape): array was not empty, and shape did not match.";
                vigra_precondition(shape == this->shape(), message.c_str());
            }
        }
        else
        {
            reshape(shape, strideOrdering);
        }
    }
};

    // this function assumes that pyArray_ has already been set, and compatibility been checked
template <unsigned int N, class T, class Stride>
void NumpyArray<N, T, Stride>::setupArrayView()
{
    if(NumpyAnyArray::hasData())
    {
        unsigned int dimension = std::min<unsigned int>(actual_dimension, pyArray()->nd);
        std::copy(pyArray()->dimensions, pyArray()->dimensions + dimension, this->m_shape.begin());
        std::copy(pyArray()->strides, pyArray()->strides + dimension, this->m_stride.begin());
        if(pyArray()->nd < actual_dimension)
        {
            this->m_shape[dimension] = 1;
            this->m_stride[dimension] = sizeof(value_type);
        }
        this->m_stride /= sizeof(value_type);
        this->m_ptr = reinterpret_cast<pointer>(pyArray()->data);
    }
    else
    {
        this->m_ptr = 0;
    }
}


typedef NumpyArray<2, float >  NumpyFArray2;
typedef NumpyArray<3, float >  NumpyFArray3;
typedef NumpyArray<4, float >  NumpyFArray4;
typedef NumpyArray<2, Singleband<float> >  NumpyFImage;
typedef NumpyArray<3, Singleband<float> >  NumpyFVolume;
typedef NumpyArray<2, RGBValue<float> >  NumpyFRGBImage;
typedef NumpyArray<3, RGBValue<float> >  NumpyFRGBVolume;
typedef NumpyArray<3, Multiband<float> >  NumpyFMultibandImage;
typedef NumpyArray<4, Multiband<float> >  NumpyFMultibandVolume;

inline void import_vigranumpy()
{
    if(_import_array() < 0)
        pythonToCppException(0);
    python_ptr module(PyImport_ImportModule("vigra.vigranumpycore"), python_ptr::keep_count);
    pythonToCppException(module);
}

/********************************************************/
/*                                                      */
/*   NumpyArray Multiband Argument Object Factories     */
/*                                                      */
/********************************************************/

template <class PixelType, class Stride>
inline triple<ConstStridedImageIterator<PixelType>,
              ConstStridedImageIterator<PixelType>,
              MultibandVectorAccessor<PixelType> >
srcImageRange(NumpyArray<3, Multiband<PixelType>, Stride> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return triple<ConstStridedImageIterator<PixelType>,
                  ConstStridedImageIterator<PixelType>,
                  MultibandVectorAccessor<PixelType> >
        (ul, ul + Size2D(img.shape(0), img.shape(1)), MultibandVectorAccessor<PixelType>(img.shape(2), img.stride(2)));
}

template <class PixelType, class Stride>
inline pair< ConstStridedImageIterator<PixelType>,
             MultibandVectorAccessor<PixelType> >
srcImage(NumpyArray<3, Multiband<PixelType>, Stride> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return pair<ConstStridedImageIterator<PixelType>, MultibandVectorAccessor<PixelType> >
        (ul, MultibandVectorAccessor<PixelType>(img.shape(2), img.stride(2)));
}

template <class PixelType, class Stride>
inline triple< StridedImageIterator<PixelType>,
               StridedImageIterator<PixelType>,
               MultibandVectorAccessor<PixelType> >
destImageRange(NumpyArray<3, Multiband<PixelType>, Stride> & img)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return triple<StridedImageIterator<PixelType>,
                  StridedImageIterator<PixelType>,
                  MultibandVectorAccessor<PixelType> >
        (ul, ul + Size2D(img.shape(0), img.shape(1)),
        MultibandVectorAccessor<PixelType>(img.shape(2), img.stride(2)));
}

template <class PixelType, class Stride>
inline pair< StridedImageIterator<PixelType>,
             MultibandVectorAccessor<PixelType> >
destImage(NumpyArray<3, Multiband<PixelType>, Stride> & img)
{
    StridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    return pair<StridedImageIterator<PixelType>, MultibandVectorAccessor<PixelType> >
        (ul, MultibandVectorAccessor<PixelType>(img.shape(2), img.stride(2)));
}

template <class PixelType, class Stride>
inline pair< ConstStridedImageIterator<PixelType>,
             MultibandVectorAccessor<PixelType> >
maskImage(NumpyArray<3, Multiband<PixelType>, Stride> const & img)
{
    ConstStridedImageIterator<PixelType>
        ul(img.data(), 1, img.stride(0), img.stride(1));
    typedef typename AccessorTraits<PixelType>::default_accessor Accessor;
    return pair<ConstStridedImageIterator<PixelType>, MultibandVectorAccessor<PixelType> >
        (ul, MultibandVectorAccessor<PixelType>(img.shape(2), img.stride(2)));
}

} // namespace vigra

#endif // VIGRA_NUMPY_ARRAY_HXX
