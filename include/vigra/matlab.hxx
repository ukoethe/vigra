/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


#ifndef VIGRA_MATLAB_HXX
#define VIGRA_MATLAB_HXX

#include <string>
#include <mex.h>
#include "array_vector.hxx"
#include "sized_int.hxx"
#include "matrix.hxx"

#include <time.h>

namespace vigra {

namespace matlab {

class InputArgumentArray
{
    int size_;
    const mxArray ** data_;
      
  public:
  
    typedef const mxArray * value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef int size_type;
    typedef int difference_type;

    InputArgumentArray(size_type size, pointer data)
    : size_(size),
      data_(data)
    {}
      
    const_reference operator[]( difference_type i ) const
    {
        if(!isValid(i))
            mexErrMsgTxt("Input parameter index out of range.");
        return data_[i];
    }

    size_type size() const
    {
        return size_;
    }

    bool isValid( difference_type i ) const
    {
        return i >= 0 && i < size_;
    }
      
};

class OutputArgumentArray
{
    int size_;
    mxArray ** data_;
      
  public:
  
    typedef mxArray * value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef int size_type;
    typedef int difference_type;

    OutputArgumentArray(size_type size, pointer data)
    : size_(size),
      data_(data)
    {}
      
    reference operator[]( difference_type i )
    {
        if(!isValid(i))
            mexErrMsgTxt("Output parameter index out of range.");
        return data_[i];
    }

    const_reference operator[]( difference_type i ) const
    {
        if(!isValid(i))
            mexErrMsgTxt("Output parameter index out of range.");
        return data_[i];
    }

    size_type size() const
    {
        return size_;
    }

    bool isValid( difference_type i ) const
    {
        return i >= 0 && i < size_;
    }
      
};

class ConstCellArray
{
  protected:
    mxArray * matPointer_;
    int size_;
      
  public:
  
    struct Proxy
    {
        mxArray * matPointer_;
        int index_;
        
        Proxy(mxArray * matPointer, int index)
        : matPointer_(matPointer),
          index_(index)
        {}
        
        operator const mxArray *() const
        {
            return mxGetCell(matPointer_, index_);
        }
    };
  
    ConstCellArray(const mxArray * matPointer, int size)
    : matPointer_(const_cast<mxArray *>(matPointer)),
      size_(size)
    {}
      
    Proxy operator[](int i) const
    {
        if(!isValid(i))
            mexErrMsgTxt("CellArray::operator[]: Index out of range.");
        return Proxy(matPointer_, i);
    }

    int size() const
    {
        return size_;
    }

    bool isValid( int i ) const
    {
        return i >= 0 && i < size_;
    }
      
};

class CellArray
: public ConstCellArray
{
  public:
  
    struct Proxy
    : public ConstCellArray::Proxy
    {
        Proxy(mxArray * matPointer, int index)
        : ConstCellArray::Proxy(matPointer, index)
        {}
        
        void operator=(mxArray * v)
        {
            mxDestroyArray(mxGetCell(matPointer_, index_));
            mxSetCell(matPointer_, index_, v);
        }
    };
  
    CellArray(const mxArray * matPointer, int size)
    : ConstCellArray(matPointer, size)
    {}
      
    Proxy operator[](int i)
    {
        if(!isValid(i))
            mexErrMsgTxt("CellArray::operator[]: Index out of range.");
        return Proxy(matPointer_, i);
    }
      
    ConstCellArray::Proxy operator[](int i) const
    {
        if(!isValid(i))
            mexErrMsgTxt("CellArray::operator[]: Index out of range.");
        return ConstCellArray::Proxy(matPointer_, i);
    }
};

namespace detail {

template <class T>
struct ValueType;

#define VIGRA_MATLAB_VALUETYPE_UTIL(type, functionName, typeID, matTypeName) \
template <> \
struct ValueType<type> \
{ \
    static bool check(mxArray const * t) \
    { \
        return mxIs##functionName(t); \
    } \
    \
    static mxClassID const classID = typeID; \
    \
    static std::string typeName() \
    { \
        return #matTypeName; \
    } \
};

VIGRA_MATLAB_VALUETYPE_UTIL(double, Double, mxDOUBLE_CLASS, double)
VIGRA_MATLAB_VALUETYPE_UTIL(float, Single, mxSINGLE_CLASS, single)
VIGRA_MATLAB_VALUETYPE_UTIL(Int8,  Int8, mxINT8_CLASS, int8)
VIGRA_MATLAB_VALUETYPE_UTIL(Int16, Int16, mxINT16_CLASS, int16)
VIGRA_MATLAB_VALUETYPE_UTIL(Int32, Int32, mxINT32_CLASS, int32)
VIGRA_MATLAB_VALUETYPE_UTIL(Int64, Int64, mxINT64_CLASS, int64)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt8,  Uint8, mxUINT8_CLASS, uint8)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt16, Uint16, mxUINT16_CLASS, uint16)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt32, Uint32, mxUINT32_CLASS, uint32)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt64, Uint64, mxUINT64_CLASS, uint64)

#undef VIGRA_MATLAB_VALUETYPE_UTIL

} // namespace detail 

// TODO: 
//    * handle complex matrices
//    * handle sparse matrices
//    * convert from/to BasicImageView

template<class T>
T
getScalar(mxArray const * t)
{
    if(mxIsEmpty(t))
        mexErrMsgTxt("getScalar() on empty input.");
    return static_cast<T>(mxGetScalar(t));
}

template<class T>
mxArray *
createScalar(T v)
{
    mxArray * m;
    createMatrix<double>(1, 1, m)(0,0) = static_cast<double>(v);
    return m;
}

template <unsigned int SIZE, class T>
TinyVectorView<T, SIZE>
getTinyVector(mxArray const * t)
{
    if(!detail::ValueType<T>::check(t))
    {
        std::string msg = std::string("Input array must have type ") + 
                          detail::ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getTinyVector(): Input array has wrong number of elements.");
    }
    
    return TinyVectorView<T, SIZE>((T *)mxGetData(t));
}

template <unsigned int SIZE>
TinyVectorView<MultiArrayIndex, SIZE>
getShape(mxArray const * t)
{
    if(!detail::ValueType<MultiArrayIndex>::check(t))
    {
        std::string msg = std::string("Input array must have type ") + 
                          detail::ValueType<MultiArrayIndex>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getShape(): Input array has wrong number of elements.");
    }
    
    return TinyVectorView<MultiArrayIndex, SIZE>((MultiArrayIndex *)mxGetData(t));
}

template <unsigned int DIM, class T>
MultiArrayView<DIM, T>
getMultiArray(mxArray const * t)
{
    typedef typename MultiArrayView<DIM, T>::difference_type Shape;

    if(!detail::ValueType<T>::check(t))
    {
        std::string msg = std::string("getMultiArray(): Input array must have type ") + 
                          detail::ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    
    Shape shape;
    if(DIM > 1)
    {
        if(DIM != mxGetNumberOfDimensions(t))
        {
            mexErrMsgTxt("getMultiArray(): Input array has wrong number of dimensions.");
        }
        const mwSize * matlabShape = mxGetDimensions(t);
        for(unsigned int k=0; k<DIM; ++k)
        {
            shape[k] = static_cast<typename Shape::value_type>(matlabShape[k]);
        }
    }
    else
    {
        shape[0] = static_cast<typename Shape::value_type>(mxGetNumberOfElements(t));
    }        
    return MultiArrayView<DIM, T>(shape, (T *)mxGetData(t));
}

template <unsigned int DIM, class T>
MultiArrayView<DIM, T>
createMultiArray(typename MultiArrayShape<DIM>::type const & shape, mxArray * & t)
{
    mwSize matlabShape[DIM];
    for(int k=0; k<DIM; ++k)
        matlabShape[k] = static_cast<mwSize>(shape[k]);
    t = mxCreateNumericArray(DIM, matlabShape, detail::ValueType<T>::classID, mxREAL);   
    
    return MultiArrayView<DIM, T>(shape, (T *)mxGetData(t));
}

template <unsigned int DIM, class T>
MultiArrayView<DIM, T>
createMultiArray(typename MultiArrayShape<DIM>::type const & shape, CellArray::Proxy t)
{
    mwSize matlabShape[DIM];
    for(int k=0; k<DIM; ++k)
        matlabShape[k] = static_cast<mwSize>(shape[k]);
    t = mxCreateNumericArray(DIM, matlabShape, detail::ValueType<T>::classID, mxREAL);   
    
    return MultiArrayView<DIM, T>(shape, (T *)mxGetData(t));
}

template <class T>
inline MultiArrayView<1, T>
getArray(mxArray const * t)
{
    return getMultiArray<1, T>(t);
}

template <class T>
inline MultiArrayView<1, T>
createArray(MultiArrayIndex size, mxArray * & t)
{
    return createMultiArray<1, T>(MultiArrayShape<1>::type(size), t);
}

template <class T>
inline MultiArrayView<1, T>
createArray(MultiArrayIndex size, CellArray::Proxy t)
{
    return createMultiArray<1, T>(MultiArrayShape<1>::type(size), t);
}

template <class T>
MultiArrayView<2, T>
getMatrix(mxArray const * t)
{
    typedef typename MultiArrayView<2, T>::difference_type Shape;

    if(!detail::ValueType<T>::check(t))
    {
        std::string msg = std::string("getMatrix(): Input matrix must have type ") + 
                          detail::ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }

    if(2 != mxGetNumberOfDimensions(t))
        mexErrMsgTxt("getMatrix(): Input matrix must have 2 dimensions.");
        
    const mwSize * matlabShape = mxGetDimensions(t);
    Shape shape(static_cast<MultiArrayIndex>(matlabShape[0]),
                static_cast<MultiArrayIndex>(matlabShape[1]));
        
    return MultiArrayView<2, T>(shape, (T *)mxGetData(t));
}

template <class T>
MultiArrayView<2, T>
createMatrix(mwSize rowCount, mwSize columnCount, mxArray * & t)
{
    typedef typename MultiArrayView<2, T>::difference_type Shape;

    Shape shape(rowCount, columnCount);
    t = mxCreateNumericMatrix(rowCount, columnCount, detail::ValueType<T>::classID, mxREAL);  
    
    return MultiArrayView<2, T>(shape, (T *)mxGetData(t));
}

template <class T>
MultiArrayView<2, T>
createMatrix(mwSize rowCount, mwSize columnCount, CellArray::Proxy t)
{
    typedef typename MultiArrayView<2, T>::difference_type Shape;

    Shape shape(rowCount, columnCount);
    t = mxCreateNumericMatrix(rowCount, columnCount, detail::ValueType<T>::classID, mxREAL);  
    
    return MultiArrayView<2, T>(shape, (T *)mxGetData(t));
}

template <class T>
BasicImageView<T>
getImage(mxArray const * t)
{
    if(!detail::ValueType<T>::check(t))
    {
        std::string msg = std::string("getImage(): Input matrix must have type ") + 
                          detail::ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }

    if(2 != mxGetNumberOfDimensions(t))
        mexErrMsgTxt("getImage(): Input matrix must have 2 dimensions.");
        
    const mwSize * matlabShape = mxGetDimensions(t);
    return BasicImageView<T>((T *)mxGetData(t), static_cast<int>(matlabShape[0]),
                                                static_cast<int>(matlabShape[1]));
}

template <class T>
BasicImageView<T>
createImage(mwSize width, mwSize height, mxArray * & t)
{
    t = mxCreateNumericMatrix(width, height, detail::ValueType<T>::classID, mxREAL);  
    
    return BasicImageView<T>((T *)mxGetData(t), width, height);
}

template <class T>
BasicImageView<T>
createImage(mwSize width, mwSize height, CellArray::Proxy t)
{
    t = mxCreateNumericMatrix(width, height, detail::ValueType<T>::classID, mxREAL);  
    
    return BasicImageView<T>((T *)mxGetData(t), width, height);
}

inline ConstCellArray
getCellArray(mxArray const * t)
{
    if(!mxIsCell(t))
        mexErrMsgTxt("getCellArray(): Input must have type CellArray.");

    return ConstCellArray(t, mxGetNumberOfElements(t));
}

inline CellArray
createCellArray(mwSize size, mxArray * & t)
{
    mwSize matSize[] = { size };
    t = mxCreateCellArray(1, matSize);  
    
    return CellArray(t, size);
}

inline CellArray
createCellArray(mwSize size, CellArray::Proxy t)
{
    mwSize matSize[] = { size };
    t = mxCreateCellArray(1, matSize);  
    
    return CellArray(t, size);
}

} // namespace matlab 

} // namespace vigra 

#endif // VIGRA_MATLAB_HXX
