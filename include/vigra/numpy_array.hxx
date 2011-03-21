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
#include <string>
#include <numpy/arrayobject.h>
#include "multi_array.hxx"
#include "array_vector.hxx"
#include "python_utility.hxx"
#include "numpy_array_utilities.hxx"
#include "numpy_array_traits.hxx"
#include "numpy_array_taggedshape.hxx"

namespace vigra {

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
/*                    constructArray                    */
/*                                                      */
/********************************************************/

// template <class TYPECODE> // pseudo-template to avoid inline expansion of the function
                          // // will always be NPY_TYPES
// PyObject * 
// constructArray(TaggedShape tagged_shape, TYPECODE typeCode, bool init,
               // python_ptr arraytype = python_ptr())
// {
    // ArrayVector<npy_intp> shape = finalizeTaggedShape(tagged_shape);
    // python_ptr axistags(tagged_shape.axistags);
    
    // int ndim = (int)shape.size();
    // ArrayVector<npy_intp> inverse_permutation;
    
    // if(axistags)
    // {
        // if(!arraytype)
            // arraytype = detail::getArrayTypeObject();

        // inverse_permutation = detail::permutationFromNormalOrder(axistags);
        // vigra_precondition(ndim == (int)inverse_permutation.size(),
                     // "axistags.permutationFromNormalOrder(): permutation has wrong size.");
    // }
    // else
    // {
        // arraytype = python_ptr((PyObject*)&PyArray_Type);

        // inverse_permutation.resize(ndim);
        // linearSequence(inverse_permutation.begin(), inverse_permutation.end(), ndim-1, -1);
    // }
    
    // python_ptr array(PyArray_New((PyTypeObject *)arraytype.get(), ndim, shape.begin(), 
                                  // typeCode, 0, 0, 0, 1 /* Fortran order */, 0),
                     // python_ptr::keep_count);
    // pythonToCppException(array);

    // PyArray_Dims permute = { inverse_permutation.begin(), ndim };
    // array = python_ptr(PyArray_Transpose((PyArrayObject*)array.get(), &permute), 
                       // python_ptr::keep_count);
    // pythonToCppException(array);
    
    // if(arraytype != (PyObject*)&PyArray_Type && axistags)
        // pythonToCppException(PyObject_SetAttrString(array, "axistags", axistags) != -1);
    
    // if(init)
        // PyArray_FILLWBYTE((PyArrayObject *)array.get(), 0);
   
    // return array.release();
// }

template <class TYPECODE> // pseudo-template to avoid inline expansion of the function
                          // will always be NPY_TYPES
PyObject * 
constructArray(TaggedShape tagged_shape, TYPECODE typeCode, bool init,
               python_ptr arraytype = python_ptr())
{
    ArrayVector<npy_intp> shape = finalizeTaggedShape(tagged_shape);
    PyAxisTags axistags(tagged_shape.axistags);
    
    int ndim = (int)shape.size();
    ArrayVector<npy_intp> inverse_permutation;
    
    if(axistags)
    {
        if(!arraytype)
            arraytype = detail::getArrayTypeObject();

        inverse_permutation = axistags.permutationFromNormalOrder();
        vigra_precondition(ndim == (int)inverse_permutation.size(),
                     "axistags.permutationFromNormalOrder(): permutation has wrong size.");
    }
    else
    {
        arraytype = python_ptr((PyObject*)&PyArray_Type);

        inverse_permutation.resize(ndim);
        linearSequence(inverse_permutation.begin(), inverse_permutation.end(), ndim-1, -1);
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
        pythonToCppException(PyObject_SetAttrString(array, "axistags", axistags.axistags) != -1);
    
    if(init)
        PyArray_FILLWBYTE((PyArrayObject *)array.get(), 0);
   
    return array.release();
}

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

    static python_ptr init(difference_type const & shape, bool init = true, 
                           std::string const & order = "")
    {
        vigra_precondition(order == "" || order == "C" || order == "F" || 
                           order == "V" || order == "A",
            "NumpyArray.init(): order must be in ['C', 'F', 'V', 'A', ''].");
        return python_ptr(constructArray(ArrayTraits::taggedShape(shape, order), typeCode, init), 
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
    explicit NumpyArray(difference_type const & shape, std::string const & order = "")
    {
        vigra_postcondition(makeReference(init(shape, true, order)),
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
            // FIXME: implement
            vigra_fail("reshapeIfEmpty(): already has data, but shape check is not implemented yet, sorry.");
        }
        else
        {
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
    
    // TaggedShape taggedShape() const
    // {
        // return ArrayTraits::taggedShape(this->shape(), 
                                         // detail::copyAxistags(this->axistags()));
    // }
    
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
