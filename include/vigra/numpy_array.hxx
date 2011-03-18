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
#include "multi_array.hxx"
#include "array_vector.hxx"
#include "sized_int.hxx"
#include "python_utility.hxx"
#include "algorithm.hxx"
#include "axistags.hxx"
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

inline long getAttrLong(PyObject * obj, PyObject * key, long defaultValue)
{
    python_ptr pres(PyObject_GetAttr(obj, key), python_ptr::keep_count);
    if(pres)
        return PyInt_Check(pres)
                     ? PyInt_AsLong(pres)
                     : defaultValue;
    PyErr_Clear();
    return defaultValue;
}

inline long spatialDimensions(PyObject * obj)
{
    static python_ptr key(PyString_FromString("spatialDimensions"), python_ptr::keep_count);
    return getAttrLong(obj, key, -1);
}

inline long channelIndex(PyArrayObject * array, long defaultVal)
{
    static python_ptr key(PyString_FromString("channelIndex"), python_ptr::keep_count);
    return getAttrLong((PyObject*)array, key, defaultVal);
}

inline long majorNonchannelIndex(PyArrayObject * array, long defaultVal)
{
    static python_ptr key(PyString_FromString("majorNonchannelIndex"), python_ptr::keep_count);
    return getAttrLong((PyObject*)array, key, defaultVal);
}

inline
python_ptr getArrayTypeObject()
{
    PyObject *g = PyEval_GetGlobals();
    python_ptr arraytype(PyRun_String("vigra.defaultArrayType", Py_eval_input, g, g), 
                         python_ptr::keep_count);
    if(!arraytype)
    {
        PyErr_Clear();
        arraytype = (PyObject*)&PyArray_Type;
    }
    return arraytype;
}

inline std::string defaultOrder(std::string defaultValue = "C")
{
    python_ptr arraytype = getArrayTypeObject();
    static python_ptr key(PyString_FromString("defaultOrder"), python_ptr::keep_count);
    python_ptr pres(PyObject_GetAttr(arraytype, key), python_ptr::keep_count);
    if(pres)
        return PyString_Check(pres)
                     ? PyString_AsString(pres)
                     : defaultValue;
    PyErr_Clear();
    return defaultValue;
}

inline 
python_ptr defaultAxistags(int ndim)
{
    PyObject *g = PyEval_GetGlobals();
    std::string command = std::string("vigra.arraytypes.defaultAxistags(") + asString(ndim) + ")";
    python_ptr axistags(PyRun_String(command.c_str(), Py_eval_input, g, g), 
                        python_ptr::keep_count);
    if(!axistags)
        PyErr_Clear(); // if no axistags are found, resume without
    
    return axistags;
}

inline 
python_ptr copyAxistags(python_ptr axistags)
{
    if(axistags)
    {
        python_ptr func(PyString_FromString("__copy__"), python_ptr::keep_count);
        axistags = python_ptr(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                              python_ptr::keep_count);
    }
    return axistags;
}

inline void setChannelDescription(python_ptr axistags, std::string const & d)
{
    python_ptr pyd(PyString_FromString(d.c_str()), python_ptr::keep_count);
    python_ptr func(PyString_FromString("setChannelDescription"), python_ptr::keep_count);
    python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), pyd.get(), NULL), 
                   python_ptr::keep_count);
    pythonToCppException(res);
}

inline 
ArrayVector<npy_intp> 
getAxisPermutationImpl(python_ptr object, const char * name, bool ignoreErrors)
{
    python_ptr func(PyString_FromString(name), python_ptr::keep_count);
    python_ptr permutation(PyObject_CallMethodObjArgs(object, func.get(), NULL), 
                           python_ptr::keep_count);
    if(!permutation && ignoreErrors)
    {
        PyErr_Clear();
        return ArrayVector<npy_intp>();
    }
    pythonToCppException(permutation);
    
    
    if(!PySequence_Check(permutation))
    {
        if(ignoreErrors)
            return ArrayVector<npy_intp>();
        std::string message = std::string(name) + "() did not return a sequence.";
        PyErr_SetString(PyExc_ValueError, message.c_str());
        pythonToCppException(false);
    }
        
    ArrayVector<npy_intp> res(PySequence_Length(permutation));
    for(int k=0; k<(int)res.size(); ++k)
    {
        python_ptr i(PySequence_GetItem(permutation, k), python_ptr::keep_count);
        if(!PyInt_Check(i))
        {
            if(ignoreErrors)
                return ArrayVector<npy_intp>();
            std::string message = std::string(name) + "() did not return a sequence of int.";
            PyErr_SetString(PyExc_ValueError, message.c_str());
            pythonToCppException(false);
        }
        res[k] = PyInt_AsLong(i);
    }
    return res;
}

inline 
ArrayVector<npy_intp> permutationToNormalOrder(python_ptr object, bool ignoreErrors = false)
{
    return getAxisPermutationImpl(object, "permutationToNormalOrder", ignoreErrors);
}

inline 
ArrayVector<npy_intp> permutationFromNormalOrder(python_ptr object, bool ignoreErrors = false)
{
    return getAxisPermutationImpl(object, "permutationFromNormalOrder", ignoreErrors);
}

#if 0

// FIXME: reimplement in terms of TaggedShape?
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

#endif

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
/*                     TaggedShape                      */
/*                                                      */
/********************************************************/

class TaggedShape
{
  public:
    enum ChannelAxis { first, last, none };
    
    ArrayVector<npy_intp> shape;
    python_ptr axistags;
    ChannelAxis channelAxis;
    std::string channelDescription;
    
    TaggedShape(MultiArrayIndex size)
    : shape(size),
      channelAxis(none)
    {}
    
    template <class U, int N>
    TaggedShape(TinyVector<U, N> const & sh, python_ptr tags = python_ptr())
    : shape(sh.begin(), sh.end()),
      axistags(tags),
      channelAxis(none)
    {}
    
    template <class T>
    TaggedShape(ArrayVector<T> const & sh, python_ptr tags = python_ptr())
    : shape(sh.begin(), sh.end()),
      axistags(tags),
      channelAxis(none)
    {}
    
    template <class U, int N>
    TaggedShape & operator=(TinyVector<U, N> const & sh)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
                        
        vigra_precondition(N == stop - start || size() == 0,
             "TaggedShape.operator=(): size mismatch.");
             
        if(size() == 0)
            shape.resize(N);
        
        for(int k=; k<N; ++k)
            shape[k+start] = sh[k];
            
        return *this;
    }
    
    npy_intp & operator[](int i)
    {
        return shape[i];
    }
    
    npy_intp operator[](int i) const
    {
        return shape[i];
    }
    
    unsigned int size() const
    {
        return shape.size();
    }
    
    TaggedShape & operator+=(int v)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
        for(int k=start; k<stop; ++k)
            shape[k] += v;
            
        return *this;
    }
    
    TaggedShape & operator-=(int v)
    {
        return operator+=(-v);
    }
    
    TaggedShape & operator*=(int factor)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
        for(int k=start; k<stop; ++k)
            shape[k] *= factor;
            
        return *this;
    }
    
    TaggedShape & setChannelDescription(std::string const & description)
    {
        // we only remember the description here, and will actually set
        // it in the constructor (after creating a copy of the axistags)
        channelDescription = description;
        return *this;
    }
    
    TaggedShape & setChannelIndexLast()
    {
        // FIXME: add some checks?
        channelAxis = last;
        return *this;
    }
    
    TaggedShape & setChannelCount(int count)
    {
        switch(channelAxis)
        {
          case first:
            shape[0] = count;
            break;
          case last:
            shape[size()-1] = count;
            break;
          case none:
            shape.push_back(count);
            channelAxis = last;
            break;
        }
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*                    constructArray                    */
/*                                                      */
/********************************************************/

inline // FIXME
PyObject * 
constructArray(TaggedShape const & tagged_shape, NPY_TYPES typeCode, bool init,
               python_ptr arraytype = python_ptr())
{
    int ndim = (int)tagged_shape.size();
    ArrayVector<npy_intp> shape(tagged_shape.shape);
    python_ptr axistags(tagged_shape.axistags);
    
    if(tagged_shape.channelAxis == TaggedShape::last)
    {
        for(int k=0; k<ndim; ++k)
            shape[k] = tagged_shape[(k-1+ndim)%ndim]; // rotate to normal order
    }
    
    ArrayVector<npy_intp> inverse_permutation;
    
    // we assume here that the axistag object belongs to the array to be created
    // so that we can freely edit it
    // FIXME: should we rather create a copy first, even if this results in some 
    //        unnecessary copies?
    if(axistags)
    {
        if(!PySequence_Check(axistags))
        {
            PyErr_SetString(PyExc_TypeError, "constructArray(): axistags have wrong type.");
            pythonToCppException(false);
        }
        
        int ntags = PySequence_Length(axistags);
        static python_ptr key(PyString_FromString("channelIndex"), python_ptr::keep_count);
        long channelIndex = detail::getAttrLong(axistags, key, ntags);

#if 0 // debug only
        std::cerr << "ndim: " << ndim << ", ntags: " << ntags << ", channelIndex: " << channelIndex << "\n";
        static python_ptr func(PyString_FromString("__repr__"), 
                               python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
        std::cerr << "axistags: " << PyString_AsString(res) << "\n";
#endif

        if(tagged_shape.channelAxis == TaggedShape::none)
        {
            // shape has no channel axis
            if(channelIndex == ntags)
            {
                // axistags have no channel axis either => sizes should match
                vigra_precondition(ndim == ntags,
                     "constructArray(): size mismatch between shape and axistags.");
            }
            else
            {
                if(ndim+1 == ntags)
                {
                    // axistags have have one additional element => drop the channel tag
                    // FIXME: would it be cleaner to make this an error ?
                    static python_ptr func(PyString_FromString("dropChannelAxis"), 
                                           python_ptr::keep_count);
                    python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                                   python_ptr::keep_count);
                    pythonToCppException(res);
                }
                else
                    vigra_precondition(ndim == ntags,
                         "constructArray(): size mismatch between shape and axistags.");
                
            }
        }
        else
        {
            // shape has a channel axis
            if(channelIndex == ntags)
            {
                // axistags have no channel axis => should be one element shorter
                vigra_precondition(ndim == ntags+1,
                     "constructArray(): size mismatch between shape and axistags.");
                     
                if(shape[0] == 1)
                {
                    // we have a singleband image => drop the channel axis
                    shape.erase(shape.begin());
                    ndim -= 1;
                }
                else
                {
                    // we have a multiband image => add a channel tag
                    static python_ptr func(PyString_FromString("insertChannelAxis"), 
                                           python_ptr::keep_count);
                    python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                                   python_ptr::keep_count);
                    pythonToCppException(res);
                }
            }
            else
            {
                // axistags have channel axis => sizes should match
                vigra_precondition(ndim == ntags,
                     "constructArray(): size mismatch between shape and axistags.");
            }
        }
            
        if(tagged_shape.channelDescription != "")
            detail::setChannelDescription(axistags, tagged_shape.channelDescription);
        
        inverse_permutation = detail::permutationFromNormalOrder(axistags);
        vigra_precondition(ndim == (int)inverse_permutation.size(),
                     "axistags.permutationFromNormalOrder(): permutation has wrong size.");
    }
    else
    {
        inverse_permutation.resize(ndim);
        linearSequence(inverse_permutation.begin(), inverse_permutation.end(), ndim-1, -1);
    }
    
    if(axistags)
    {
        if(!arraytype)
            arraytype = detail::getArrayTypeObject();
    }
    else
    {
        arraytype = python_ptr((PyObject*)&PyArray_Type);
    }
    
    python_ptr array(PyArray_New((PyTypeObject *)arraytype.get(), ndim, shape.begin(), 
                                  typeCode, 0, 0, 0, 1 /* Fortran order */, 0),
                     python_ptr::keep_count);
    pythonToCppException(array);

    PyArray_Dims permute = { inverse_permutation.begin(), ndim };
    array = python_ptr(PyArray_Transpose((PyArrayObject*)array.get(), &permute), 
                       python_ptr::keep_count);
    pythonToCppException(array);
    
    if(arraytype != (PyObject*)&PyArray_Type && axistags)
        pythonToCppException(PyObject_SetAttrString(array, "axistags", axistags));
    
    if(init)
        PyArray_FILLWBYTE((PyArrayObject *)array.get(), 0);
   
    return array.release();
}

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
        long channelIndex = detail::channelIndex(array, ndim);
        long majorIndex = detail::majorNonchannelIndex(array, ndim);

        if(channelIndex < ndim)
        {
            // When we have a channel axis, there are two cases:
            // 1. ndim is right: everything is ok
            // 2. ndim == N+1: we drop the channel axis when it is a singleton
            return (ndim == N) ||
                    (ndim == N+1 && PyArray_DIM(obj, channelIndex) == 1);
        }
        else if(majorIndex < ndim)
        {
            // We have axistags, but no channel axis. There are again two cases:
            // 1. ndim is right: everything is ok
            // 2. ndim == N-1: we add a singleton channel axis later
            return ndim == N-1 || ndim == N;
        }
        else
        {
            // We have no axistags. 
            // When ndim == N, everything is ok
            // When ndim == N-1, we add a singleton axis later
            // When ndim == N+1, we assume that channelIndex == ndim-1, and we may drop the
            // channel axis when it is a singleton.
            return ndim == N-1 || ndim == N ||
                    (ndim == N+1 && PyArray_DIM(obj, ndim-1) == 1);
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && isValuetypeCompatible(obj);
    }
    
    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, python_ptr axistags)
    {
        return TaggedShape(shape, axistags);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape)
    {
        // FIXME: we construct axistags with one entry too many, so that the channel axis
        //        can be dropped in constructArray(). This is a HACK.
        return TaggedShape(shape, detail::defaultAxistags(shape.size()+1));
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        vigra_precondition(tagged_shape.size() == N,
                  "reshapeIfEmpty(): tagged_shape has wrong size.");
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
        long channelIndex = detail::channelIndex(array, ndim);
        long majorIndex = detail::majorNonchannelIndex(array, ndim);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        if(channelIndex < ndim)
        {
            // When we have a channel axis, there are two cases:
            // 1. ndim is right: the channel axis is the major axis and must be unstrided
            // 2. ndim == N+1: we drop the channel axis when it is a singleton, and require the
            //                 major non-channel axis to be unstrided
            return (ndim == N && strides[channelIndex] == itemsize) ||
                    (ndim == N+1 && PyArray_DIM(obj, channelIndex) == 1 && strides[majorIndex] == itemsize);
        }
        else if(majorIndex < ndim)
        {
            // We have axistags, but no channel axis. There are again two cases:
            // 1. ndim is right: the major axis must be unstrided
            // 2. ndim == N-1: we add a singleton channel axis later, which is automatically unstrided,
            //                 so there is no explicit stride requirement
            return ndim == N-1 || (ndim == N && strides[majorIndex] == itemsize);
        }
        else
        {
            // We have no axistags. 
            // When ndim == N or ndim == N+1, we assume that
            //     channelIndex == ndim-1
            //     majorIndex == ndim-2
            // When ndim == N-1, we assume
            //     there is no channel axis
            //     majorIndex == ndim-1
            // and proceed as above
            return ndim == N-1 ||
                    (ndim == N && strides[ndim-1] == itemsize)  ||
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
struct NumpyArrayTraits<N, Singleband<T>, StridedArrayTag>
: public NumpyArrayTraits<N, T, StridedArrayTag>
{
    typedef NumpyArrayTraits<N, T, StridedArrayTag> BaseType;
    typedef typename BaseType::ValuetypeTraits ValuetypeTraits;

    static bool isShapeCompatible(PyArrayObject * array) /* array must not be NULL */
    {
        PyObject * obj = (PyObject *)array;
		int ndim = PyArray_NDIM(obj);
        long channelIndex = detail::channelIndex(array, N);
        
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
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, python_ptr axistags)
    {
        return TaggedShape(shape, axistags).setChannelCount(1);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape)
    {
        return TaggedShape(shape, detail::defaultAxistags(shape.size()+1)).setChannelCount(1);
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        tagged_shape.setChannelCount(1);
        vigra_precondition(tagged_shape.size() == N+1,
              "reshapeIfEmpty(): tagged_shape has wrong size.");
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
        long channelIndex = detail::channelIndex(array, ndim);
        long majorIndex = detail::majorNonchannelIndex(array, ndim);
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
        int ndim = PyArray_NDIM((PyObject *)array);
        long channelIndex = detail::channelIndex(array, ndim);
        
        if(channelIndex < ndim)
        {
            // When we have a channel axis, ndim must be right.
            return ndim == N;
        }
        else
        {
            // Otherwise, we may add a singleton channel axis, if necessary.
            return ndim == N || ndim == N-1;
        }
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && ValuetypeTraits::isValuetypeCompatible(obj);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, python_ptr axistags)
    {
        return TaggedShape(shape, axistags).setChannelIndexLast();
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape)
    {
        return TaggedShape(shape, detail::defaultAxistags(shape.size())).setChannelIndexLast();
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
        long channelIndex = detail::channelIndex(array, ndim);
        long majorIndex = detail::majorNonchannelIndex(array, ndim);
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
            
        long channelIndex = detail::channelIndex(array, N);
        int itemsize = PyArray_ITEMSIZE(obj);
        npy_intp * strides = PyArray_STRIDES(obj);
        
        return PyArray_DIM(obj, channelIndex) == M && strides[channelIndex] == itemsize;
    }

    static bool isPropertyCompatible(PyArrayObject * obj) /* obj must not be NULL */
    {
        return isShapeCompatible(obj) && ValuetypeTraits::isValuetypeCompatible(obj);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape, python_ptr axistags)
    {
        return TaggedShape(shape, axistags).setChannelCount(M);
    }

    template <class U>
    static TaggedShape taggedShape(TinyVector<U, N> const & shape)
    {
        return TaggedShape(shape, detail::defaultAxistags(shape.size()+1)).setChannelCount(M);
    }

    static void finalizeTaggedShape(TaggedShape & tagged_shape)
    {
        tagged_shape.setChannelCount(M);
        vigra_precondition(tagged_shape.size() == N+1,
              "reshapeIfEmpty(): tagged_shape has wrong size.");
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
            
        long channelIndex = detail::channelIndex(array, ndim);
        long majorIndex = detail::majorNonchannelIndex(array, ndim);
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

    <b>\#include</b> \<vigra/numpy_array.hxx\><br>
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
         Returns the the permutation that will transpose this array into 
         canonical ordering (currently: F-order). The size of
         the returned permutation equals ndim().
         */
    difference_type permutationToNormalOrder() const
    {
		if(!hasData())
            return difference_type();
            
        difference_type res(detail::permutationToNormalOrder(pyArray_, true));
        if(res.size() == 0)
        {
            res.resize(ndim());
            linearSequence(res.begin(), res.end(), ndim()-1, MultiArrayIndex(-1));
        }
        return res;
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
         * Return the AxisTags of this array or a NULL pointer when the attribute
           'axistags' is missing in the Python object.
         */
    python_ptr axistags() const
    {
        static python_ptr key(PyString_FromString("axistags"), python_ptr::keep_count);
        python_ptr axistags(PyObject_GetAttr(pyObject(), key), python_ptr::keep_count);
        if(!axistags)
            PyErr_Clear();
        return axistags;
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

    <b>\#include</b> \<vigra/numpy_array.hxx\><br>
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
        return detail::getArrayTypeObject();
    }

    static python_ptr init(difference_type const & shape, bool init = true)
    {
        return python_ptr(constructArray(ArrayTraits::taggedShape(shape), typeCode, init), 
                          python_ptr::keep_count);
    }

#if 0 // FIXME: not sure if this is still needed when we have axistags
    static python_ptr init(difference_type const & shape, difference_type const & strideOrdering, bool init = true)
    {
        // FIXME: what to do with this function?
        ArrayVector<npy_intp> pshape(shape.begin(), shape.end()),
                              pstrideOrdering(strideOrdering.begin(), strideOrdering.end());
        if(pstrideOrdering.size() == ArrayTraits::spatialDimensions)
        {
            for(unsigned int k=0; k < pstrideOrdering.size(); ++k)
                pstrideOrdering[k] += 1;
            pstrideOrdering.push_back(0);
        }
        return detail::constructNumpyArrayImpl((PyTypeObject *)getArrayTypeObject().ptr(), pshape,
                       ArrayTraits::spatialDimensions, ArrayTraits::channels,
                       typeCode, "A", init, pstrideOrdering);
    }
#endif

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
    NumpyArray(const NumpyArray &other, bool createCopy = false) 
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

#if 0
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
#endif

        /**
         * Constructor from NumpyAnyArray.
         * Equivalent to NumpyArray(other.pyObject())
         */
    explicit NumpyArray(const NumpyAnyArray &other, bool createCopy = false)
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
#if VIGRA_CONVERTER_DEBUG
        std::cerr << "class " << typeid(NumpyArray).name() << " got " << obj->ob_type->tp_name << "\n";
        std::cerr << "using traits " << typeid(ArrayTraits).name() << "\n";
        std::cerr<<"isArray: "<< ArrayTraits::isArray(obj)<<std::endl;
        std::cerr<<"isShapeCompatible: "<< ArrayTraits::isShapeCompatible((PyArrayObject *)obj)<<std::endl;
#endif

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
         * Deprecated, use isReferenceCompatible(obj) instead.
         */
    static bool isStrictlyCompatible(PyObject *obj)
    {
        return isReferenceCompatible(obj);
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
         * numpy array (see isReferenceCompatible()).
         *
         * The parameter \a strict is deprecated and will be ignored
         */
    bool makeReference(PyObject *obj, bool strict = false)
    {
        if(!isReferenceCompatible(obj))
            return false;
        makeReferenceUnchecked(obj);
        return true;
    }

        /**
         * Try to set up a view referencing the same data as the given
         * NumpyAnyArray.  This overloaded variant simply calls
         * makeReference() on array.pyObject(). The parameter \a strict
         * is deprecated and will be ignored.
         */
    bool makeReference(const NumpyAnyArray &array, bool strict = false)
    {
        return makeReference(array.pyObject(), strict);
    }

#if 0 // FIXME: implement this in a different way
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
#endif

        /**
         Try to create a copy of the given PyObject.
         Raises an exception when obj is not a compatible array
         (see isCopyCompatible() or isReferenceCompatible(), according to the
         parameter \a strict) or the Python constructor call failed.
         */
    void makeCopy(PyObject *obj, bool strict = false)
    {
        vigra_precondition(strict ? isReferenceCompatible(obj) : isCopyCompatible(obj),
                     "NumpyArray::makeCopy(obj): Cannot copy an incompatible array.");

        NumpyAnyArray copy(obj, true);
        makeReferenceUnchecked(copy.pyObject());
    }

#if 0
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
#endif

        /**
            Allocate new memory with the given shape and initialize with zeros.<br>
            If a stride ordering is given, the resulting array will have this stride
            ordering, when it is compatible with the array's memory layout (unstrided
            arrays only permit the standard ascending stride ordering).

            <em>Note:</em> this operation invalidates dependent objects
            (MultiArrayViews and iterators)
         */
    void reshape(difference_type const & shape)
    {
        vigra_postcondition(makeReference(init(shape)),
                "NumpyArray.reshape(shape): Python constructor did not produce a compatible array.");
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
            and throw a precondition exception with the given \a message if not.
         */
    void reshapeIfEmpty(TaggedShape tagged_shape, std::string message = "")
    {
        if(hasData())
        {
            vigra_fail("reshapeIfEmpty(): already has data, but shape check is not implemented yet, sorry.");
        }
        else
        {
            tagged_shape.axistags = detail::copyAxistags(tagged_shape.axistags);
            ArrayTraits::finalizeTaggedShape(tagged_shape);

            python_ptr array(constructArray(tagged_shape, typeCode, true), 
                             python_ptr::keep_count);
            vigra_postcondition(makeReference(NumpyAnyArray(array.get())),
                  "NumpyArray.reshapeIfEmpty(): Python constructor did not produce a compatible array.");
        }
    }

#if 0
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
#endif
    
    TaggedShape taggedShape() const
    {
        return ArrayTraits::taggedShape(this->shape(), this->axistags());
    }
};

    // this function assumes that pyArray_ has already been set, and compatibility been checked
template <unsigned int N, class T, class Stride>
void NumpyArray<N, T, Stride>::setupArrayView()
{
    bool channelDimensionMustBeLast = typeid(T) == typeid(Multiband<value_type>);
                        
    if(NumpyAnyArray::hasData())
    {
        NumpyAnyArray::difference_type ordering = permutationToNormalOrder();
        
        if(actual_dimension == pyArray()->nd)
        {
            if(channelDimensionMustBeLast)
            {
                // rotate the channel axis to last position
                for(int k=1; k<actual_dimension; ++k)
                {
                    this->m_shape[k-1] = pyArray()->dimensions[ordering[k]];
                    this->m_stride[k-1] = pyArray()->strides[ordering[k]];
                }
                this->m_shape[actual_dimension-1] = pyArray()->dimensions[ordering[0]];
                this->m_stride[actual_dimension-1] = pyArray()->strides[ordering[0]];
            }
            else
            {
                for(int k=0; k<actual_dimension; ++k)
                {
                    this->m_shape[k] = pyArray()->dimensions[ordering[k]];
                    this->m_stride[k] = pyArray()->strides[ordering[k]];
                }
            }
        }
        else if(actual_dimension == pyArray()->nd - 1)
        {
            // skip the channel axis
            for(int k=0; k<actual_dimension; ++k)
            {
                this->m_shape[k] = pyArray()->dimensions[ordering[k+1]];
                this->m_stride[k] = pyArray()->strides[ordering[k+1]];
            }
        }
        else if(actual_dimension == pyArray()->nd + 1)
        {
            if(channelDimensionMustBeLast)
            {
                // insert a singleton channel at the last position
                for(int k=0; k<actual_dimension-1; ++k)
                {
                    this->m_shape[k] = pyArray()->dimensions[ordering[k]];
                    this->m_stride[k] = pyArray()->strides[ordering[k]];
                }
                this->m_shape[actual_dimension-1] = 1;
                this->m_stride[actual_dimension-1] = sizeof(value_type);
            }
            else
            {
                // insert a singleton channel at the first position
                for(int k=0; k<actual_dimension-1; ++k)
                {
                    this->m_shape[k+1] = pyArray()->dimensions[ordering[k]];
                    this->m_stride[k+1] = pyArray()->strides[ordering[k]];
                }
                this->m_shape[0] = 1;
                this->m_stride[0] = sizeof(value_type);
            }
        }
        else
        {
            vigra_precondition(false,
              "NumpyArray::setupArrayView(): got array of incompatible shape (should never happen).");
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

#ifdef import_array
#undef import_array
#endif

inline void import_array()
{
    if(_import_array() < 0)
        pythonToCppException(0);
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
