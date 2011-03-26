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

#ifndef VIGRA_NUMPY_ARRAY_TRAITS_HXX
#define VIGRA_NUMPY_ARRAY_TRAITS_HXX

#include "numerictraits.hxx"
#include "multi_array.hxx"
#include "numpy_array_taggedshape.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*              Singleband and Multiband                */
/*                                                      */
/********************************************************/

template <class T>
struct Singleband  // the resulting NumpyArray has no explicit channel axis 
                   // (i.e. the number of channels is implicitly one)
{
    typedef T value_type;
};

template <class T>
struct Multiband  // the last axis is explicitly designated as channel axis
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

template<unsigned int N, class T, class Stride>
struct NumpyArrayTraits;

/********************************************************/

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

    static bool isValuetypeCompatible(PyArrayObject * obj)  /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj);
    }

    static bool isShapeCompatible(PyArrayObject * array) /* array must not be NULL */
    {
		PyObject * obj = (PyObject *)array;
		int ndim = PyArray_NDIM(obj);

        return ndim == N;
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && isValuetypeCompatible(obj);
    }
    
    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, PyAxisTags axistags)
    {
        return TaggedShape(shape, axistags);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, std::string const & order = "")
    {
        return TaggedShape(shape, PyAxisTags(detail::emptyAxistags(shape.size())));
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        vigra_precondition(tagged_shape.size() == N,
                  "reshapeIfEmpty(): tagged_shape has wrong size.");
    }
    
    template <class U>
    static void permutationToNormalOrder(python_ptr array, AxisInfo::AxisType types, 
                                         ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       types, true);

        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
    }
    
    template <class U>
    static void permutationToSetupOrder(python_ptr array, ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AllAxes, true);

        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
    }

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, T, UnstridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* obj must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
        int ndim = PyArray_NDIM(obj);
        long channelIndex = pythonGetAttr(obj, "channelIndex", ndim);
        long majorIndex = pythonGetAttr(obj, "majorNonchannelIndex", ndim);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        if(channelIndex < ndim)
        {
            return (ndim == N && strides[channelIndex] == itemsize);
        }
        else if(majorIndex < ndim)
        {
            return (ndim == N && strides[majorIndex] == itemsize);
        }
        else
        {
            return (ndim == N && strides[ndim-1] == itemsize);
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && BaseType::isValuetypeCompatible(obj);
    }

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Singleband<T>, StridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* array must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
		int ndim = PyArray_NDIM(obj);
        long channelIndex = pythonGetAttr(obj, "channelIndex", N);
        
        // When ndim is right, this array must not have an explicit channel axis.
        if(ndim == N)
            return channelIndex == N;
            
        // When we have one extra axis, we allow to drop it, provided it is
        // a channel axis and the number of channels is 1.
        // When no explicit channel axis is known, we use the last axis by
        // default (we automaticaly get 'channelIndex == N' in this case).
        return ndim == N+1 && PyArray_DIM(obj, channelIndex) == 1;
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && BaseType::isValuetypeCompatible(obj);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, PyAxisTags axistags)
    {
        return TaggedShape(shape, axistags).setChannelCount(1);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, std::string const & order = "")
    {
        return TaggedShape(shape, 
                  PyAxisTags(detail::defaultAxistags(shape.size()+1, order))).setChannelCount(1);
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        tagged_shape.setChannelCount(1);
        vigra_precondition(tagged_shape.size() == N+1,
              "reshapeIfEmpty(): tagged_shape has wrong size.");
    }
    
    template <class U>
    static void permutationToNormalOrder(python_ptr array, AxisInfo::AxisType types, 
                                         ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AxisType(types & ~AxisInfo::Channels), true);

        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
    }
    
    template <class U>
    static void permutationToSetupOrder(python_ptr array, ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AllAxes, true);
        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
        else if(permute.size() == N+1)
        {
            permute.erase(permute.begin());
        }
    }
    
    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Singleband<T>, UnstridedArrayTag>
: public NumpyArrayTraits<N, Singleband<T>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, UnstridedArrayTag> UnstridedTraits;
    typedef NumpyArrayTraits<N, Singleband<T>, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* obj must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
        int ndim = PyArray_NDIM(obj);
        long channelIndex = pythonGetAttr(obj, "channelIndex", ndim);
        long majorIndex = pythonGetAttr(obj, "majorNonchannelIndex", ndim);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        if(channelIndex < ndim)
        {
            // When we have a channel axis, it must be a singleton, so that we can drop it.
            // Moreover, the major axis among the remaining ones must be unstrided.
            // (Note that existence of a valid channelIndex implies a valid majorIndex.)
            return ndim == N+1 && PyArray_DIM(obj, channelIndex) == 1 && strides[majorIndex] == itemsize;
        }
        else if(majorIndex < ndim)
        {
            // We have axistags, but no channel axis
            // => ndim must be right, and the majorAxis must be unstrided
            return ndim == N && strides[majorIndex] == itemsize;
        }
        else
        {
            // We have no axistags. When ndim == N, we assume that
            //     there is no channel index and
            //     majorIndex == ndim-1
            // When ndim == N+1, we assume that
            //     channelIndex == ndim-1
            //     majorIndex == ndim-2
            // and require the channel axis to be a singleton.
            return (ndim == N && strides[ndim-1] == itemsize)  ||
                    (ndim == N+1 && PyArray_DIM(obj, ndim-1) == 1 && strides[ndim-2] == itemsize);
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && BaseType::isValuetypeCompatible(obj);
    }

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Multiband<T>, StridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    enum { spatialDimensions = N-1, channels = 0 };

    static bool isShapeCompatible(PyArrayObject * array) /* array must not be NULL */
    {
        PyObject * obj = (PyObject*)array;
        int ndim = PyArray_NDIM(obj);
        long channelIndex = pythonGetAttr(obj, "channelIndex", ndim);
        long majorIndex = pythonGetAttr(obj, "majorNonchannelIndex", ndim);
        
        if(channelIndex < ndim)
        {
            // When we have a channel axis, ndim must be right.
            return ndim == N;
        }
        else if(majorIndex < ndim)
        {
            // When we have axistags, but no channel axis, we must add a singleton axis.
            return ndim == N-1;
        }
        else
        {
            // Otherwise, we may or may not add a singleton channel axis.
            return ndim == N || ndim == N-1;
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && ValuetypeTraits::isValuetypeCompatible(obj);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, PyAxisTags axistags)
    {
        return TaggedShape(shape, axistags).setChannelIndexLast();
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, std::string const & order = "")
    {
        return TaggedShape(shape, 
                    PyAxisTags(detail::defaultAxistags(shape.size(), order))).setChannelIndexLast();
    }

    template <class U>
    static void permutationToNormalOrder(python_ptr array, AxisInfo::AxisType types, 
                                         ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       types, true);

        if(permute.size() == 0)
        {
            permute.resize(N-1);
            linearSequence(permute.begin(), permute.end(), int(N-2), -1);
            if((types & AxisInfo::Channels) != 0)
                permute.push_back(N-1);
        }
        else if((types & AxisInfo::Channels) != 0)
        {
            // if permutation with channels is requested, 
            //   rotate the channel axis to last position
            int channelIndex = permute[0];
            for(int k=1; k<N; ++k)
                permute[k-1] = permute[k];
            permute[N-1] = channelIndex;
        }
    }
    
    template <class U>
    static void permutationToSetupOrder(python_ptr array, ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AllAxes, true);

        if(permute.size() == 0)
        {
            permute.resize(N-1);
            linearSequence(permute.begin(), permute.end(), int(N-2), -1);
            if(PyArray_NDIM((PyArrayObject*)array.get()) == N)
                permute.push_back(N-1);
        }
        else if(permute.size() == N)
        {
            // if we have a channel axis, rotate it to last position
            int channelIndex = permute[0];
            for(int k=1; k<N; ++k)
                permute[k-1] = permute[k];
            permute[N-1] = channelIndex;
        }
    }
    
    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, Multiband<T>, UnstridedArrayTag>
: public NumpyArrayTraits<N, Multiband<T>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, Multiband<T>, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* obj must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
        int ndim = PyArray_NDIM(obj);
        long channelIndex = pythonGetAttr(obj, "channelIndex", ndim);
        long majorIndex = pythonGetAttr(obj, "majorNonchannelIndex", ndim);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);

        if(channelIndex < ndim)
        {
            // When we have a channel axis, ndim must be right, and the major non-channel
            // axis must be unstrided.
            return ndim == N && strides[majorIndex] == itemsize*PyArray_DIM(obj, channelIndex);
        }
        else if(majorIndex < ndim)
        {
            // We have axistags, but no channel axis
            // => We will add a singleton channel axis later, and the major axis must be unstrided.
            return ndim == N-1 && strides[majorIndex] == itemsize;
        }
        else
        {
            // We have no axistags. When ndim == N, we assume that
            //     channelIndex == ndim-1
            //     majorIndex == ndim-2
            // When ndim == N-1, we assume
            //     there is no channel axis
            //     majorIndex == ndim-1
            return (ndim == N && strides[ndim-2] == itemsize*PyArray_DIM(obj, ndim-1))  ||
                    (ndim == N-1 && strides[ndim-1] == itemsize);
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && BaseType::isValuetypeCompatible(obj);
    }

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N> npyStride(stride * sizeof(T));
        // return detail::constructNumpyArrayFromData(typeKeyFull(), BaseType::typeKey(), shape, npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
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

    static bool isValuetypeCompatible(PyArrayObject * obj)  /* obj must not be NULL */
    {
        return ValuetypeTraits::isValuetypeCompatible(obj);
    }

    static bool isShapeCompatible(PyArrayObject * array) /* array must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
        
        if(PyArray_NDIM(obj) != N+1) // We need an extra channel axis.
            return false;
            
        long channelIndex = pythonGetAttr(obj, "channelIndex", N);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        return PyArray_DIM(obj, channelIndex) == M && strides[channelIndex] == itemsize;
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && ValuetypeTraits::isValuetypeCompatible(obj);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, PyAxisTags axistags)
    {
        return TaggedShape(shape, axistags).setChannelCount(M);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, std::string const & order = "")
    {
        return TaggedShape(shape, 
                     PyAxisTags(detail::defaultAxistags(shape.size()+1, order))).setChannelCount(M);
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        tagged_shape.setChannelCount(M);
        vigra_precondition(tagged_shape.size() == N+1,
              "reshapeIfEmpty(): tagged_shape has wrong size.");
    }

    template <class U>
    static void permutationToNormalOrder(python_ptr array, AxisInfo::AxisType types, 
                                         ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AxisType(types & ~AxisInfo::Channels), true);

        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
    }
    
    template <class U>
    static void permutationToSetupOrder(python_ptr array, ArrayVector<U> & permute)
    {
        detail::getAxisPermutationImpl(permute, array, "permutationToNormalOrder", 
                                       AxisInfo::AllAxes, true);
        if(permute.size() == 0)
        {
            permute.resize(N);
            linearSequence(permute.begin(), permute.end(), int(N-1), -1);
        }
        else if(permute.size() == N+1)
        {
            permute.erase(permute.begin());
        }
    }
    
    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N+1> npyShape;
        // std::copy(shape.begin(), shape.end(), npyShape.begin());
        // npyShape[N] = M;

        // TinyVector<npy_intp, N+1> npyStride;
        // std::transform(
            // stride.begin(), stride.end(), npyStride.begin(),
            // std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        // npyStride[N] = sizeof(T);

        // return detail::constructNumpyArrayFromData(
            // typeKeyFull(), typeKey(), npyShape,
            // npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, int M, class T>
struct NumpyArrayTraits<N, TinyVector<T, M>, UnstridedArrayTag>
: public NumpyArrayTraits<N, TinyVector<T, M>, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, TinyVector<T, M>, StridedArrayTag> BaseType;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* obj must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
        int ndim = PyArray_NDIM(obj);
        
        if(PyArray_NDIM(obj) != N+1) // We need an extra channel axis.
            return false;
            
        long channelIndex = pythonGetAttr(obj, "channelIndex", ndim);
        long majorIndex = pythonGetAttr(obj, "majorNonchannelIndex", ndim);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        if(channelIndex == ndim)
        {
            if(majorIndex != ndim) // we have axis tags, but no channel axis 
                return false;     // => cannot be a vector image
                
            return PyArray_DIM(obj, N) == M && 
                   strides[N] == itemsize &&
                   strides[N-1] == sizeof(TinyVector<T, M>);
        }
        else
        {
            return PyArray_DIM(obj, channelIndex) == M && 
                   strides[channelIndex] == itemsize &&
                   strides[majorIndex] == sizeof(TinyVector<T, M>);
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && BaseType::isValuetypeCompatible(obj);
    }

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N+1> npyShape;
        // std::copy(shape.begin(), shape.end(), npyShape.begin());
        // npyShape[N] = M;

        // TinyVector<npy_intp, N+1> npyStride;
        // std::transform(
            // stride.begin(), stride.end(), npyStride.begin(),
            // std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        // npyStride[N] = sizeof(T);

        // return detail::constructNumpyArrayFromData(
            // typeKeyFull(), BaseType::typeKey(), npyShape,
            // npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

/********************************************************/

template<unsigned int N, class T>
struct NumpyArrayTraits<N, RGBValue<T>, StridedArrayTag>
: public NumpyArrayTraits<N, TinyVector<T, 3>, StridedArrayTag>
{
    typedef T dtype;
    typedef RGBValue<T> value_type;
    typedef NumpyArrayValuetypeTraits<T> ValuetypeTraits;

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N+1> npyShape;
        // std::copy(shape.begin(), shape.end(), npyShape.begin());
        // npyShape[N] = 3;

        // TinyVector<npy_intp, N+1> npyStride;
        // std::transform(
            // stride.begin(), stride.end(), npyStride.begin(),
            // std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        // npyStride[N] = sizeof(T);

        // return detail::constructNumpyArrayFromData(
            // typeKeyFull(), typeKey(), npyShape,
            // npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
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

    // template <class U>
    // static python_ptr constructor(TinyVector<U, N> const & shape,
                                  // T *data, TinyVector<U, N> const & stride)
    // {
        // TinyVector<npy_intp, N+1> npyShape;
        // std::copy(shape.begin(), shape.end(), npyShape.begin());
        // npyShape[N] = 3;

        // TinyVector<npy_intp, N+1> npyStride;
        // std::transform(
            // stride.begin(), stride.end(), npyStride.begin(),
            // std::bind2nd(std::multiplies<npy_intp>(), sizeof(value_type)));
        // npyStride[N] = sizeof(T);

        // return detail::constructNumpyArrayFromData(
            // typeKeyFull(), BaseType::typeKey(), npyShape,
            // npyStride.begin(), ValuetypeTraits::typeCode, data);
    // }
};

} // namespace vigra

#endif // VIGRA_NUMPY_ARRAY_TRAITS_HXX
