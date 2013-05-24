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

#ifndef NPY_NO_DEPRECATED_API
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif 

#include <Python.h>
#include <string>
#include <iostream>
#include <numpy/arrayobject.h>
#include "multi_array.hxx"
#include "array_vector.hxx"
#include "python_utility.hxx"
#include "numpy_array_traits.hxx"
#include "numpy_array_taggedshape.hxx"

// NumPy function called by NumPy's import_array() macro (and our import_vigranumpy() below)
int _import_array();

namespace vigra {

static inline void import_vigranumpy()
{
    // roughly equivalent to import_array():
    if(_import_array() < 0)
        pythonToCppException(0);

    // Import vigra to activate the numpy array converters, but ensure that 
    // cyclic imports (from within vigra itself) are avoided.
    char const * load_vigra = 
        "import sys\n"
        "if not sys.modules.has_key('vigra.vigranumpycore'):\n"
        "    import vigra\n";
    pythonToCppException(PyRun_SimpleString(load_vigra) == 0);
}

/********************************************************/
/*                                                      */
/*               MultibandVectorAccessor                */
/*                                                      */
/********************************************************/

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
            In case of a conversion floating point -> integral this includes rounding and clipping.
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
            In case of a conversion floating point -> integral this includes rounding and clipping.
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

template <class TYPECODE> // pseudo-template to avoid inline expansion of the function
                          // will always be NPY_TYPES
PyObject *
constructArray(TaggedShape tagged_shape, TYPECODE typeCode, bool init,
               python_ptr arraytype = python_ptr());

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

  public:

        /// difference type
    typedef ArrayVector<npy_intp> difference_type;

    static python_ptr getArrayTypeObject()
    {
        return detail::getArrayTypeObject();
    }

    static std::string defaultOrder(std::string defaultValue = "C")
    {
        return detail::defaultOrder(defaultValue);
    }

    static python_ptr defaultAxistags(int ndim, std::string order = "")
    {
        return detail::defaultAxistags(ndim, order);
    }

    static python_ptr emptyAxistags(int ndim)
    {
        return detail::emptyAxistags(ndim);
    }

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

            python_ptr arraytype = getArrayTypeObject();
            python_ptr f(PyString_FromString("_copyValuesImpl"), python_ptr::keep_count);
            if(PyObject_HasAttr(arraytype, f))
            {
                python_ptr res(PyObject_CallMethodObjArgs(arraytype, f.get(),
                                                          pyArray_.get(), other.pyArray_.get(), NULL),
                               python_ptr::keep_count);
                vigra_postcondition(res.get() != 0,
                       "NumpyArray::operator=(): VigraArray._copyValuesImpl() failed.");
            }
            else
            {
                PyArrayObject * sarray = (PyArrayObject *)pyArray_.get();
                PyArrayObject * tarray = (PyArrayObject *)other.pyArray_.get();

                if(PyArray_CopyInto(tarray, sarray) == -1)
                    pythonToCppException(0);
            }
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
            return PyArray_NDIM(pyArray());
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
        return pythonGetAttr(pyObject(), "spatialDimensions", ndim());
    }

    bool hasChannelAxis() const
    {
        if(!hasData())
            return false;
        return channelIndex() == ndim();
    }

    MultiArrayIndex channelIndex() const
    {
        if(!hasData())
            return 0;
        return pythonGetAttr(pyObject(), "channelIndex", ndim());
    }

    MultiArrayIndex innerNonchannelIndex() const
    {
        if(!hasData())
            return 0;
        return pythonGetAttr(pyObject(), "innerNonchannelIndex", ndim());
    }

        /**
         Returns the shape of this array. The size of
         the returned shape equals ndim().
         */
    difference_type shape() const
    {
        if(hasData())
            return difference_type(PyArray_DIMS(pyArray()), PyArray_DIMS(pyArray()) + ndim());
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
        difference_type stride(PyArray_STRIDES(pyArray()), PyArray_STRIDES(pyArray()) + N),
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

        // /**
         // Returns the the permutation that will transpose this array into
         // canonical ordering (currently: F-order). The size of
         // the returned permutation equals ndim().
         // */
    // difference_type permutationToNormalOrder() const
    // {
        // if(!hasData())
            // return difference_type();

        // // difference_type res(detail::getAxisPermutationImpl(pyArray_,
                                               // // "permutationToNormalOrder", true));
        // difference_type res;
        // detail::getAxisPermutationImpl(res, pyArray_, "permutationToNormalOrder", true);
        // if(res.size() == 0)
        // {
            // res.resize(ndim());
            // linearSequence(res.begin(), res.end(), ndim()-1, MultiArrayIndex(-1));
        // }
        // return res;
    // }

        /**
         Returns the value type of the elements in this array, or -1
         when hasData() is false.
         */
    int dtype() const
    {
        if(hasData())
            return PyArray_DESCR(pyArray())->type_num;
        return -1;
    }

        /**
         * Return the AxisTags of this array or a NULL pointer when the attribute
           'axistags' is missing in the Python object or this array has no data.
         */
    python_ptr axistags() const
    {
        python_ptr axistags;
        if(pyObject())
        {
            python_ptr key(PyString_FromString("axistags"), python_ptr::keep_count);
            axistags.reset(PyObject_GetAttr(pyObject(), key), python_ptr::keep_count);
            if(!axistags)
                PyErr_Clear();
        }
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
/*                    constructArray                    */
/*                                                      */
/********************************************************/

namespace detail {

inline bool
nontrivialPermutation(ArrayVector<npy_intp> const & p)
{
    for(unsigned int k=0; k<p.size(); ++k)
        if(p[k] != k)
            return true;
    return false;
}

} // namespace detail

template <class TYPECODE> // pseudo-template to avoid inline expansion of the function
                          // will always be NPY_TYPES
PyObject *
constructArray(TaggedShape tagged_shape, TYPECODE typeCode, bool init, python_ptr arraytype)
{
    ArrayVector<npy_intp> shape = finalizeTaggedShape(tagged_shape);
    PyAxisTags axistags(tagged_shape.axistags);

    int ndim = (int)shape.size();
    ArrayVector<npy_intp> inverse_permutation;
    int order = 1; // Fortran order

    if(axistags)
    {
        if(!arraytype)
            arraytype = NumpyAnyArray::getArrayTypeObject();

        inverse_permutation = axistags.permutationFromNormalOrder();
        vigra_precondition(ndim == (int)inverse_permutation.size(),
                     "axistags.permutationFromNormalOrder(): permutation has wrong size.");
    }
    else
    {
        arraytype = python_ptr((PyObject*)&PyArray_Type);
        order = 0; // C order
    }

//    std::cerr << "constructArray: " << shape << "\n" << inverse_permutation << "\n";

    python_ptr array(PyArray_New((PyTypeObject *)arraytype.get(), ndim, shape.begin(),
                                  typeCode, 0, 0, 0, order, 0),
                     python_ptr::keep_count);
    pythonToCppException(array);

    if(detail::nontrivialPermutation(inverse_permutation))
    {
        PyArray_Dims permute = { inverse_permutation.begin(), ndim };
        array = python_ptr(PyArray_Transpose((PyArrayObject*)array.get(), &permute),
                           python_ptr::keep_count);
        pythonToCppException(array);
    }

    if(arraytype != (PyObject*)&PyArray_Type && axistags)
        pythonToCppException(PyObject_SetAttrString(array, "axistags", axistags.axistags) != -1);

    if(init)
        PyArray_FILLWBYTE((PyArrayObject *)array.get(), 0);

    return array.release();
}

// FIXME: reimplement in terms of TaggedShape?
template <class TINY_VECTOR>
inline
python_ptr constructNumpyArrayFromData(
    TINY_VECTOR const & shape, npy_intp *strides,
    NPY_TYPES typeCode, void *data)
{
    ArrayVector<npy_intp> pyShape(shape.begin(), shape.end());

#if NPY_NO_DEPRECATED_API < NPY_1_7_API_VERSION    // old API
    python_ptr array(PyArray_New(&PyArray_Type, shape.size(), pyShape.begin(),
                                 typeCode, strides, data, 0, NPY_WRITEABLE, 0),
                     python_ptr::keep_count);
#else
    python_ptr array(PyArray_New(&PyArray_Type, shape.size(), pyShape.begin(),
                                 typeCode, strides, data, 0, NPY_ARRAY_WRITEABLE, 0),
                     python_ptr::keep_count);
#endif
    pythonToCppException(array);

    return array;
}

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

        /** type of an array specifying an axis permutation
         */
    typedef typename NumpyAnyArray::difference_type permutation_type;

        /** traverser type
         */
    typedef typename view_type::traverser traverser;

        /** traverser type to const data
         */
    typedef typename view_type::const_traverser const_traverser;

        /** sequential (random access) iterator type
         */
    typedef typename view_type::iterator iterator;

        /** sequential (random access) const iterator type
         */
    typedef typename view_type::const_iterator const_iterator;

    using view_type::shape;   // resolve ambiguity of multiple inheritance
    using view_type::hasData; // resolve ambiguity of multiple inheritance
    using view_type::strideOrdering; // resolve ambiguity of multiple inheritance

  protected:

    // this function assumes that pyArray_ has already been set, and compatibility been checked
    void setupArrayView();

    static python_ptr init(difference_type const & shape, bool init = true,
                           std::string const & order = "")
    {
        vigra_precondition(order == "" || order == "C" || order == "F" ||
                           order == "V" || order == "A",
            "NumpyArray.init(): order must be in ['C', 'F', 'V', 'A', ''].");
        return python_ptr(constructArray(ArrayTraits::taggedShape(shape, order), typeCode, init),
                          python_ptr::keep_count);
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
         * This will not happen if isReferenceCompatible(obj) (in case
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
    : view_type(),
      NumpyAnyArray()
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
    template <class U, class S>
    explicit NumpyArray(const MultiArrayView<N, U, S> &other)
    {
        if(!other.hasData())
            return;
        vigra_postcondition(makeReference(init(other.shape(), false)),
                  "NumpyArray(MultiArrayView): Python constructor did not produce a compatible array.");
        view_type::operator=(other);
    }

        /**
         * Construct a new array object, allocating an internal python
         * ndarray of the given shape in the given order (default: VIGRA order), initialized
         * with zeros.
         *
         * An exception is thrown when construction fails.
         */
    explicit NumpyArray(difference_type const & shape, std::string const & order = "")
    {
        vigra_postcondition(makeReference(init(shape, true, order)),
                     "NumpyArray(shape): Python constructor did not produce a compatible array.");
    }

        /**
         * Construct a new array object, allocating an internal python
         * ndarray according to the given tagged shape, initialized with zeros.
         *
         * An exception is thrown when construction fails.
         */
    explicit NumpyArray(TaggedShape const & tagged_shape)
    {
        reshapeIfEmpty(tagged_shape,
           "NumpyArray(tagged_shape): Python constructor did not produce a compatible array.");
    }

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
         * array contents are copied.  If this is an empty view,
         * assignment is identical to makeReferenceUnchecked(other.pyObject()).
         * See MultiArrayView::operator= for further information on
         * semantics.
         */
    template <class U, class S>
    NumpyArray &operator=(const NumpyArray<N, U, S> &other)
    {
        if(hasData())
        {
            vigra_precondition(shape() == other.shape(),
                "NumpyArray::operator=(): shape mismatch.");
            view_type::operator=(other);
        }
        else if(other.hasData())
        {
            NumpyArray copy;
            copy.reshapeIfEmpty(other.taggedShape(),
                "NumpyArray::operator=(): reshape failed unexpectedly.");
            copy = other;
            makeReferenceUnchecked(copy.pyObject());
        }
        return *this;
    }

        /**
         * Assignment operator. If this is already a view with data
         * (i.e. hasData() is true) and the shapes match, the RHS
         * array contents are copied.  If this is an empty view,
         * a new buffer with the RHS shape is allocated before copying.
         */
    template <class U, class S>
    NumpyArray &operator=(const MultiArrayView<N, U, S> &other)
    {
        if(hasData())
        {
            vigra_precondition(shape() == other.shape(),
                "NumpyArray::operator=(): shape mismatch.");
            view_type::operator=(other);
        }
        else if(other.hasData())
        {
            NumpyArray copy;
            copy.reshapeIfEmpty(other.shape(),
                "NumpyArray::operator=(): reshape failed unexpectedly.");
            copy = other;
            makeReferenceUnchecked(copy.pyObject());
        }
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
        else if(isReferenceCompatible(other.pyObject()))
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
         Permute the entries of the given array \a data exactly like the axes of this NumpyArray
         were permuted upon conversion from numpy.
         */
    template <class U>
    ArrayVector<U>
    permuteLikewise(ArrayVector<U> const & data) const
    {
        vigra_precondition(hasData(),
            "NumpyArray::permuteLikewise(): array has no data.");

        ArrayVector<U> res(data.size());
        ArrayTraits::permuteLikewise(this->pyArray_, data, res);
        return res;
    }

        /**
         Permute the entries of the given array \a data exactly like the axes of this NumpyArray
         were permuted upon conversion from numpy.
         */
    template <class U, int K>
    TinyVector<U, K>
    permuteLikewise(TinyVector<U, K> const & data) const
    {
        vigra_precondition(hasData(),
            "NumpyArray::permuteLikewise(): array has no data.");

        TinyVector<U, K> res;
        ArrayTraits::permuteLikewise(this->pyArray_, data, res);
        return res;
    }

        /**
         Get the permutation of the axes of this NumpyArray
         that was performed upon conversion from numpy.
         */
    template <int K>
    TinyVector<npy_intp, K>
    permuteLikewise() const
    {
        vigra_precondition(hasData(),
            "NumpyArray::permuteLikewise(): array has no data.");

        TinyVector<npy_intp, K> data, res;
        linearSequence(data.begin(), data.end());
        ArrayTraits::permuteLikewise(this->pyArray_, data, res);
        return res;
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
         * The second parameter ('strict') is deprecated and will be ignored.
         */
    bool makeReference(PyObject *obj, bool /* strict */ = false)
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
    void makeUnsafeReference(const view_type &multiArrayView)
    {
        vigra_precondition(!hasData(),
            "makeUnsafeReference(): cannot replace existing view with given buffer");

        // construct an ndarray that points to our data (taking strides into account):
        python_ptr array(ArrayTraits::unsafeConstructorFromData(multiArrayView.shape(),
                                  multiArrayView.data(), multiArrayView.stride()));

        view_type::operator=(multiArrayView);
        pyArray_ = array;
    }

        /**
         Try to create a copy of the given PyObject.
         Raises an exception when obj is not a compatible array
         (see isCopyCompatible() or isReferenceCompatible(), according to the
         parameter \a strict) or the Python constructor call failed.
         */
    void makeCopy(PyObject *obj, bool strict = false)
    {
#if VIGRA_CONVERTER_DEBUG
        int ndim = PyArray_NDIM((PyArrayObject *)obj);
        npy_intp * s = PyArray_DIMS((PyArrayObject *)obj);
        std::cerr << "makeCopy: " << ndim << " " <<  ArrayVectorView<npy_intp>(ndim, s) <<
                     ", strides " << ArrayVectorView<npy_intp>(ndim, PyArray_STRIDES((PyArrayObject *)obj)) << "\n";
        std::cerr << "for " << typeid(*this).name() << "\n";
#endif
        vigra_precondition(strict ? isReferenceCompatible(obj) : isCopyCompatible(obj),
                     "NumpyArray::makeCopy(obj): Cannot copy an incompatible array.");

        NumpyAnyArray copy(obj, true);
        makeReferenceUnchecked(copy.pyObject());
    }

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
        // FIXME: is this really a good replacement?
        // reshapeIfEmpty(shape, standardStrideOrdering(), message);
        reshapeIfEmpty(TaggedShape(shape), message);
    }

        /**
            When this array has no data, allocate new memory with the given \a shape and
            initialize with zeros. Otherwise, check if the new shape matches the old shape
            and throw a precondition exception with the given \a message if not.
         */
    void reshapeIfEmpty(TaggedShape tagged_shape, std::string message = "")
    {
        ArrayTraits::finalizeTaggedShape(tagged_shape);

        if(hasData())
        {
            vigra_precondition(tagged_shape.compatible(taggedShape()), message.c_str());
        }
        else
        {
            python_ptr array(constructArray(tagged_shape, typeCode, true),
                             python_ptr::keep_count);
            vigra_postcondition(makeReference(NumpyAnyArray(array.get())),
                  "NumpyArray.reshapeIfEmpty(): Python constructor did not produce a compatible array.");
        }
    }

    TaggedShape taggedShape() const
    {
        return ArrayTraits::taggedShape(this->shape(), PyAxisTags(this->axistags(), true));
    }
};

    // this function assumes that pyArray_ has already been set, and compatibility been checked
template <unsigned int N, class T, class Stride>
void NumpyArray<N, T, Stride>::setupArrayView()
{
    if(NumpyAnyArray::hasData())
    {
        permutation_type permute;
        ArrayTraits::permutationToSetupOrder(this->pyArray_, permute);

        vigra_precondition(abs((int)permute.size() - actual_dimension) <= 1,
            "NumpyArray::setupArrayView(): got array of incompatible shape (should never happen).");

        applyPermutation(permute.begin(), permute.end(),
                         PyArray_DIMS(pyArray()), this->m_shape.begin());
        applyPermutation(permute.begin(), permute.end(),
                         PyArray_STRIDES(pyArray()), this->m_stride.begin());

        if((int)permute.size() == actual_dimension - 1)
        {
            this->m_shape[actual_dimension-1] = 1;
            this->m_stride[actual_dimension-1] = sizeof(value_type);
        }

        this->m_stride /= sizeof(value_type);
        this->m_ptr = reinterpret_cast<pointer>(PyArray_DATA(pyArray()));
        vigra_precondition(this->checkInnerStride(Stride()),
            "NumpyArray<..., UnstridedArrayTag>::setupArrayView(): First dimension of given array is not unstrided (should never happen).");

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
