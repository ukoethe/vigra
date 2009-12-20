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
#include <map>
#include <time.h>
#include "matlab_FLEXTYPE.hxx"

namespace vigra {

namespace matlab {

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
VIGRA_MATLAB_VALUETYPE_UTIL(UInt8,  Uint8, mxUINT8_CLASS, uint8)
VIGRA_MATLAB_VALUETYPE_UTIL(Int16, Int16, mxINT16_CLASS, int16)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt16, Uint16, mxUINT16_CLASS, uint16)

#if VIGRA_BITSOF_INT == 32
VIGRA_MATLAB_VALUETYPE_UTIL(int, Int32, mxINT32_CLASS, int32)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned int, Uint32, mxUINT32_CLASS, uint32)
#elif VIGRA_BITSOF_INT == 64
VIGRA_MATLAB_VALUETYPE_UTIL(int, Int64, mxINT64_CLASS, int64)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned int, Uint64, mxUINT64_CLASS, uint64)
#endif

#if VIGRA_BITSOF_LONG == 32
VIGRA_MATLAB_VALUETYPE_UTIL(long, Int32, mxINT32_CLASS, int32)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned long, Uint32, mxUINT32_CLASS, uint32)
#elif VIGRA_BITSOF_LONG == 64
VIGRA_MATLAB_VALUETYPE_UTIL(long, Int64, mxINT64_CLASS, int64)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned long, Uint64, mxUINT64_CLASS, uint64)
#endif

#if VIGRA_BITSOF_LONG_LONG == 32
VIGRA_MATLAB_VALUETYPE_UTIL(long long, Int32, mxINT32_CLASS, int32)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned long long, Uint32, mxUINT32_CLASS, uint32)
#elif VIGRA_BITSOF_LONG_LONG == 64
VIGRA_MATLAB_VALUETYPE_UTIL(long long, Int64, mxINT64_CLASS, int64)
VIGRA_MATLAB_VALUETYPE_UTIL(unsigned long long, Uint64, mxUINT64_CLASS, uint64)
#endif

#undef VIGRA_MATLAB_VALUETYPE_UTIL

class ConstStructArray
{
  protected:
    mxArray * matPointer_;

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
            return mxGetFieldByNumber(matPointer_, 0, index_);
        }
    };

    ConstStructArray(const mxArray * matPointer = 0)
    : matPointer_(const_cast<mxArray *>(matPointer))
    {
        if(matPointer != 0 && !mxIsStruct(matPointer))
            mexErrMsgTxt("StructArray(mxArray *): Argument must be a Matlab struct array.");
    }

    Proxy operator[](const char * field_name) const
    {
        if(matPointer_ == 0)
            mexErrMsgTxt("StructArray::operator[]: Cannot access uninitialized struct array.");

        int i = mxGetFieldNumber(matPointer_, field_name);
        if(i == -1)
            mexErrMsgTxt("StructArray::operator[]: Unknown field name.");

        return Proxy(matPointer_, i);
    }

    Proxy operator[](std::string field_name) const
    {
        return operator[](field_name.c_str());
    }

    bool isValid() const
    {
        return matPointer_ != 0;
    }

    bool isValid(const char * field_name) const
    {
        return isValid() && mxGetFieldNumber(matPointer_, field_name) != -1;
    }

    bool isValid(std::string field_name) const
    {
        return isValid(field_name.c_str());
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

    ConstCellArray(const mxArray * matPointer = 0)
    : matPointer_(const_cast<mxArray *>(matPointer)),
      size_(0)
    {
        if(matPointer != 0 && !mxIsCell(matPointer))
            mexErrMsgTxt("CellArray(mxArray *): Argument must be a Matlab cell array.");
        if(matPointer != 0)
            size_ = mxGetNumberOfElements(matPointer);
        else
            size_ = -1;
    }

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

    CellArray(const mxArray * matPointer)
    : ConstCellArray(matPointer)
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





template <class T, unsigned int SIZE>
TinyVectorView<T, SIZE>
getTinyVector(mxArray const * t)
{
    if(!ValueType<T>::check(t))
    {
        std::string msg = std::string("Input array must have type ") +
                          ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getTinyVector(): Input array has wrong number of elements.");
    }

    return TinyVectorView<T, SIZE>((T *)mxGetData(t));
}

template <unsigned int SIZE>
typename MultiArrayShape<SIZE>::type
getShape(mxArray const * t)
{
    if(!ValueType<Int32>::check(t))
    {
        std::string msg = std::string("Input array must have type 'int32'.");
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getShape(): Input array has wrong number of elements.");
    }
    TinyVectorView<Int32, SIZE> res((MultiArrayIndex *)mxGetData(t));
    return typename MultiArrayShape<SIZE>::type(res);
}

template <unsigned int DIM, class T>
MultiArrayView<DIM, T>
getMultiArray(mxArray const * t)
{
    typedef typename MultiArrayView<DIM, T>::difference_type Shape;

    if(!ValueType<T>::check(t))
    {
        std::string msg = std::string("getMultiArray(): Input array must have type ") +
                          ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }

    Shape shape;
    if(DIM > 1)
    {
        int mdim = mxGetNumberOfDimensions(t);
        if(static_cast<int>(DIM) < mdim)
        {
            mexErrMsgTxt("getMultiArray(): Input array has too many dimensions.");
        }
        const mwSize * matlabShape = mxGetDimensions(t);
        for(int k=0; k<mdim; ++k)
        {
            shape[k] = static_cast<typename Shape::value_type>(matlabShape[k]);
        }
        for(int k=mdim; k<static_cast<int>(DIM); ++k)
        {
            shape[k] = 1;
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
    for(int k=0; k<static_cast<int>(DIM); ++k)
        matlabShape[k] = static_cast<mwSize>(shape[k]);
    t = mxCreateNumericArray(DIM, matlabShape, ValueType<T>::classID, mxREAL);

    return MultiArrayView<DIM, T>(shape, (T *)mxGetData(t));
}

template <unsigned int DIM, class T>
MultiArrayView<DIM, T>
createMultiArray(typename MultiArrayShape<DIM>::type const & shape, CellArray::Proxy t)
{
    mwSize matlabShape[DIM];
    for(int k=0; k<DIM; ++k)
        matlabShape[k] = static_cast<mwSize>(shape[k]);
    t = mxCreateNumericArray(DIM, matlabShape, ValueType<T>::classID, mxREAL);

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

    if(!ValueType<T>::check(t))
    {
        std::string msg = std::string("getMatrix(): Input matrix must have type ") +
                          ValueType<T>::typeName() + ".";
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
    t = mxCreateNumericMatrix(rowCount, columnCount, ValueType<T>::classID, mxREAL);

    return MultiArrayView<2, T>(shape, (T *)mxGetData(t));
}

template <class T>
MultiArrayView<2, T>
createMatrix(mwSize rowCount, mwSize columnCount, CellArray::Proxy t)
{
    typedef typename MultiArrayView<2, T>::difference_type Shape;

    Shape shape(rowCount, columnCount);
    t = mxCreateNumericMatrix(rowCount, columnCount, ValueType<T>::classID, mxREAL);

    return MultiArrayView<2, T>(shape, (T *)mxGetData(t));
}

template <class T>
BasicImageView<T>
getImage(mxArray const * t)
{
    if(!ValueType<T>::check(t))
    {
        std::string msg = std::string("getImage(): Input matrix must have type ") +
                          ValueType<T>::typeName() + ".";
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
    t = mxCreateNumericMatrix(width, height, ValueType<T>::classID, mxREAL);

    return BasicImageView<T>((T *)mxGetData(t), width, height);
}

template <class T>
BasicImageView<T>
createImage(mwSize width, mwSize height, CellArray::Proxy t)
{
    t = mxCreateNumericMatrix(width, height, ValueType<T>::classID, mxREAL);

    return BasicImageView<T>((T *)mxGetData(t), width, height);
}

inline ConstCellArray
getCellArray(mxArray const * t)
{
    return ConstCellArray(t);
}

inline CellArray
createCellArray(mwSize size, mxArray * & t)
{
    mwSize matSize[] = { size };
    t = mxCreateCellArray(1, matSize);

    return CellArray(t);
}

inline CellArray
createCellArray(mwSize size, CellArray::Proxy t)
{
    mwSize matSize[] = { size };
    t = mxCreateCellArray(1, matSize);

    return CellArray(t);
}

inline ConstStructArray
getStructArray(mxArray const * t)
{
    return ConstStructArray(t);
}

template<class T>
T
getScalar(mxArray const * t)
{
    if(mxIsEmpty(t))
        mexErrMsgTxt("getScalar() on empty input.");
    if(!mxIsNumeric(t) && !mxIsLogical(t))
        mexErrMsgTxt("getScalar(): argument is not numeric.");
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

inline std::string
getString(mxArray const * t)
{
    if(mxIsEmpty(t))
        mexErrMsgTxt("getString() on empty input.");
    if(!mxIsChar(t))
        mexErrMsgTxt("getString(): argument is not a string.");
    int size = mxGetNumberOfElements(t) + 1;
    ArrayVector<char> buf(size);
    mxGetString(t, buf.begin(), size);
    return std::string(buf.begin());
}



class CompileTimeError;

namespace detail {

class Required
{
  public:
    void argumentWasProvided() const { /* empty because required arguments are always provided */ }
};


template<class T>
class DefaultImpl
{
  public:

    T defaultValue_;
    mutable bool * argumentWasProvided_;

    DefaultImpl(T v, bool * argFlag = 0)
    : defaultValue_(v),
      argumentWasProvided_(argFlag)
    {
        if(argumentWasProvided_ != 0)
            *argumentWasProvided_ = false;
    }

    void argumentWasProvided() const
    {
        if(argumentWasProvided_ != 0)
            *argumentWasProvided_ = true;
    }
};

class OptionalImpl
{
  public:
    mutable bool * argumentWasProvided_;

    OptionalImpl(bool * argFlag = 0)
    : argumentWasProvided_(argFlag)
    {
        if(argumentWasProvided_ != 0)
            *argumentWasProvided_ = false;
    }

    void argumentWasProvided() const
    {
        if(argumentWasProvided_ != 0)
            *argumentWasProvided_ = true;
    }
};

} // namespace detail

inline detail::Required v_required()
{
    return detail::Required();
}

template<class T>
inline detail::DefaultImpl<T> v_default(T in)
{
    return detail::DefaultImpl<T>(in);
}

template<class T>
inline detail::DefaultImpl<T> v_default(T in, bool & argFlag)
{
    return detail::DefaultImpl<T>(in, &argFlag);
}

inline detail::OptionalImpl v_optional()
{
    return detail::OptionalImpl();
}

inline detail::OptionalImpl v_optional(bool& argFlag)
{
    return detail::OptionalImpl(&argFlag);
}

// TODO:
//    * handle rgb images
//    * handle complex matrices
//    * handle sparse matrices

class InputArray
{
    int size_;
    const mxArray ** data_;

    std::string createErrMsg(std::string name)
    {
        std::string s1;
        s1 =  "Required input '" + name + "' not found in option struct!";
        return s1;
    }
    std::string createErrMsg(int pos)
    {
        char tmp[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
        std::string oi(1, tmp[pos%10]);
        std::string s1  = "Required input in signature of function at position: '"+ oi+"' has not been supplied";
        return s1;
    }


  public:
    ConstStructArray options_;

    /* Local Typedefs */
    typedef const mxArray * value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef int size_type;
    typedef int difference_type;

    /*Constructor*/
    InputArray(size_type size, pointer data)
    : size_(size),
      data_(data),
      options_(isValid(size-1) && mxIsStruct(data_[size-1])
                 ? data_[size-1]
                 : 0)
    {}

    /*Operators*/
    const_reference operator[]( difference_type i ) const
    {
        if(!isValid(i))
            mexErrMsgTxt("Too few input arguments.");
        return data_[i];
    }

    value_type operator[]( std::string name) const
    {
        std::string errMsg = "Not Found " + name +" in OptionStruct or OptionStruct not set";
        if(!isValid(name))
            mexErrMsgTxt(errMsg.c_str());
        return options_[name];
    }


    /*Some More Helper Func*/
    size_type size() const
    {
        return size_;
    }

    bool isValid( difference_type i ) const
    {
        return i >= 0 && i < size_;
    }

    bool isValid(std::string name) const
    {
        return options_.isValid(name);
    }

    bool isEmpty(difference_type i) const
    {
        return mxIsEmpty(data_[i]);
    }

    bool isEmpty(std::string name) const
    {
        return mxIsEmpty(options_[name]);
    }

    bool hasData(difference_type i) const
    {
        return isValid(i) && !isEmpty(i);
    }

    bool hasData(std::string name) const
    {
        return isValid(name) && !isEmpty(name);
    }

    template<class Place>
    mxClassID typeOf(Place posOrName)
    {
        return mxGetClassID((*this)[posOrName]);
    }

    /*Action to take if value not set*/
    template <class T, class U, class Place>
    T errorOrDefault(detail::DefaultImpl<U> const & o, Place posOrName)
    {
        return o.defaultValue_;
    }

    template <class T, class Place>
    T errorOrDefault(detail::OptionalImpl, Place posOrName)
    {
        return T();
    }

    template <class T, class Place>
    T errorOrDefault(detail::Required r, Place posOrName)
    {
        std::string a = createErrMsg(posOrName);
        mexErrMsgTxt( a.c_str());
        return T();
    }

    /*getter Func*/
    template <class Place, class ReqType>
    int getEnum(Place posOrName, ReqType req, std::map<std::string, int> const & converter)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<int>(req, posOrName);
        }
        std::string enumAsString = matlab::getString((*this)[posOrName]);
        typename std::map<std::string, int>::const_iterator m = converter.find(enumAsString);
        if(m == converter.end())
        {
            std::string msg = std::string("Unknown option: ") + enumAsString + ".";
            mexErrMsgTxt(msg.c_str());
        }

        req.argumentWasProvided();
        return (*m).second;
    }


    /*String Type*/
    template <class Place, class ReqType>
    std::string getString(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<std::string>(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            return matlab::getString((*this)[posOrName]);
        }
    }

    /*Scalar Type*/
    template <class T,class Place, class ReqType>
    T getScalar(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<T>(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            return matlab::getScalar<T>((*this)[posOrName]);
        }
    }


    template <class T, class Place, class ReqType, class minClass, class maxClass>
    T getScalarMinMax(Place posOrName, ReqType req, minClass min_, maxClass max_)
    {
        T temp = this->getScalar<T>(posOrName, req);
        if (!is_in_range(temp, min_, max_))
            mexErrMsgTxt("Value out of bounds.");

        return temp;
    }

    template <class T, class Place, class ReqType, class iteratorType>
    T getScalarVals(Place posOrName, ReqType req, iteratorType begin_, iteratorType end_)
    {
        T temp = this->getScalar<T>(posOrName, req);
        for(iteratorType iter = begin_; iter != end_; ++iter)
        {
            if((*iter) == temp) return temp;
        }
        mexErrMsgTxt("Value not allowed");
    }



    template <class T, class Place, class ReqType, class iteratorType>
    T getScalarVals2D3D(Place posOrName, ReqType req, iteratorType begin2D_, iteratorType end2D_,
                                                     iteratorType begin3D_, iteratorType end3D_,
                                                     int dimVar)
    {
        T temp = this->getScalar<T>(posOrName, req);
        switch(dimVar)
        {
            case 2:
                for(iteratorType iter = begin2D_; iter != end2D_; ++iter)
                {
                    if((*iter) == temp) return temp;
                }
                break;
            case 3:
                for(iteratorType iter = begin3D_; iter != end3D_; ++iter)
                {
                    if((*iter) == temp) return temp;
                }
                break;
            default:
                mexErrMsgTxt("dimVar specified must be 2 or 3");
        }
        mexErrMsgTxt("Value not allowed");
    }

    template <class Place, class ReqType>
    bool getBool(Place posOrName, ReqType req)
    {
        return this->getScalarMinMax<int>(posOrName, req, 0, 1) != 0;
    }

    /*Array Type*/
    template <unsigned int N, class T, class Place, class ReqType>
    MultiArrayView<N,T> getMultiArray(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault< MultiArrayView<N,T> >(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            value_type temp = (*this)[posOrName];
            return matlab::getMultiArray<N,T>(temp);
        }
    }

    template < class T, class Place, class ReqType>
    BasicImageView<T> getImage(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<BasicImageView<T> >(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            value_type temp = (*this)[posOrName];
            return matlab::getImage<T>(temp);
        }
    }

    template<class T,unsigned int sze, class Place, class ReqType>
    TinyVectorView< T, sze> getTinyVector(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<TinyVectorView< T, sze> >(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            value_type temp = (*this)[posOrName];
            return matlab::getTinyVector< T, sze>(temp);
        }
    }

    template< unsigned int sze, class Place, class ReqType>
    TinyVectorView<MultiArrayIndex, sze> getShape(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<TinyVectorView<MultiArrayIndex, sze> >(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            value_type temp = (*this)[posOrName];
            return matlab::getShape<sze>(temp);
        }
    }


    template< class Place, class ReqType>
    int getDimOfInput(Place posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<int>(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            return mxGetNumberOfDimensions((*this)[posOrName]);
        }
    }

    template<class ReqType>
    ConstCellArray getCellArray(int posOrName, ReqType req)
    {
        if(!hasData(posOrName))
        {
            return errorOrDefault<ConstCellArray>(req, posOrName);
        }
        else
        {
            req.argumentWasProvided();
            value_type temp = (*this)[posOrName];
            return matlab::getCellArray(temp);
        }
    }

    template<class ReqType>
    ConstCellArray getCellArray(std::string posOrName, ReqType req)
    {
        CompileTimeError ERROR__Const_Cell_Array_May_Not_Be_In_Option_Struct;
        return ConstCellArray();
    }

};

class OutputArray
{
    int size_;
    mxArray ** data_;
    std::string createErrMsgOut(int pos)
    {
        char tmp[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
        std::string oi(1, tmp[pos%10]);
        std::string s1 =  "Required Output at position: '" + oi + "' has not been supplied";
        return s1;
    }
  public:

    typedef mxArray * value_type;
    typedef value_type & reference;
    typedef value_type const & const_reference;
    typedef value_type * pointer;
    typedef value_type const * const_pointer;
    typedef int size_type;
    typedef int difference_type;

    OutputArray(size_type size, pointer data)
    : size_(size),
      data_(data)
    {}

    reference operator[]( difference_type i )
    {
        if(!isValid(i))
            mexErrMsgTxt("Too few output arguments.");
        return data_[i];
    }

    const_reference operator[]( difference_type i ) const
    {
        if(!isValid(i))
            mexErrMsgTxt("Too few output arguments.");
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

    bool isEmpty(difference_type i){
        return mxIsEmpty(data_[i]);
    }

    template <class T>
    T errorOrDefault(detail::OptionalImpl const & o, int Pos)
    {
        return T();
    }

    template <class T>
    T errorOrDefault(detail::Required r, int Pos)
    {
        mexErrMsgTxt(createErrMsgOut(Pos).c_str());
        return T();
    }

    /* creating func */
    template <unsigned int DIM, class T, class ReqType>
    MultiArrayView<DIM, T> createMultiArray(int pos,ReqType req,
                                            const TinyVector<int, DIM>  & shape)
    {
        if(!isValid(pos))
            return errorOrDefault<MultiArrayView<DIM, T> >(req, pos);
        req.argumentWasProvided();
        return matlab::createMultiArray<DIM, T>(shape, (*this)[pos]);
    }

    template <class T, class ReqType>
    BasicImageView<T> createImage(int pos, ReqType req,
                                    mwSize width, mwSize height)
    {
        if(!isValid(pos))
            return errorOrDefault<BasicImageView<T> >(req, pos);
        req.argumentWasProvided();
        return matlab::createImage<T>(width, height, (*this)[pos]);
    }

    template <class T, class ReqType>
    BasicImageView<T> createImage(  int pos, ReqType req,
                                    typename MultiArrayShape<2>::type const & shape)
    {
        return createImage<T>(pos, req, shape[1], shape[0]);
    }

    template <class T, class ReqType>
    T* createScalar(int pos, ReqType req)
    {
        if(!isValid(pos))
            return errorOrDefault<T*>(req, pos);
        req.argumentWasProvided();
        BasicImageView<T> temp = matlab::createImage<T>(1, 1, (*this)[pos]);
        return &temp(0,0);
    }

    template <class T, class ReqType>
    void createScalar(int pos, ReqType req, T val)
    {
        if(!isValid(pos))
        {
            errorOrDefault<T>(req, pos);
            return;
        }
        req.argumentWasProvided();
        BasicImageView<T> temp = matlab::createImage<T>(1, 1, (*this)[pos]);
        temp(0,0) = val;
    }

    template <class ReqType>
    ConstCellArray createCellArray(int pos, ReqType req, mwSize sze)
    {
        if(!isValid(pos))
            return errorOrDefault<ConstCellArray>(req, pos);
        return matlab::createCellArray(sze, (*this)[pos]);
    }
};



/***********************************
Rahuls code starts here
************************************+*/
using namespace vigra;


/*++++++++++++++++++++++++++HELPERFUNC+++++++++++++++++++++++++++++++
 * This is used for better readibility of the test cases            .
 * Nothing to be done here.
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int cantorPair(int x, int y){
        return (int)(((x+y)*(x+y+1))/2+y);
}

int cantorPair(int x, int y, int z){
        return cantorPair(cantorPair(x,y),z);
}

template <int x, int y>
struct cP{
    enum { value = (int)(((x+y)*(x+y+1))/2+y)};
};

template <int x, int y, int z>
struct cP3{
    enum { value = cP<cP<x, y>::value, z>::value};
};

template <class T>
inline bool is_in_range(T in, T min, T max)
{
    return (in >= min && in <= max);
}
template<class T>
inline bool is_in_range(T in, std::string min, T max)
{
    return(in <= max);
}

template<class T>
inline bool is_in_range(T in, T min, std::string max)
{
    return (in >= min);
}



//Wrapper classes to STL-Map for use as a sparse array.

//This is used for the ordering of the map. Lexicographical ordering of the index pairs.
struct ShapeCmp {
  bool operator()( TinyVector<int,2> s1, TinyVector<int,2>  s2 ) const {
    if(s1[0] != s2[0]){
        return (s1[0] < s2[0]);
    } else {
        return s1[1] < s2[1];
    }
  }
};

template<class T>
class SparseArray
{

    std::map<TinyVector<int,2>, T,ShapeCmp> data;
    int width, length;

    public:
    void assign(int i = 1, int j = 1){
        width = j;
        length = i;
    }
    SparseArray(int i = 1 , int j = 1){
        width = j;
        length = i;
    }

    //Any better idea? i would like to unify the get and operator() functions.
    // Problem is that  operator() always passes a reference or creates one.
    template<class indexType>
    T& operator()(indexType i_, indexType j_){
        Int32 i = static_cast<Int32>(i_);
        Int32 j = static_cast<Int32>(j_);
        TinyVector<int,2> newShapew(i, j);
        typename std::map<TinyVector<int,2>, T, ShapeCmp>::iterator iter;
        TinyVector<int,2> newShape;
        return data[newShapew];
    }

    template<class indexType>
    const T get(indexType i_, indexType j_){
        Int32 i = static_cast<Int32>(i_);
        Int32 j = static_cast<Int32>(j_);
        TinyVector<int,2> newShape(i, j);
        if(data.find(newShape) == data.end()) return 0;
        else return data.find(newShape)->second;
    }

    //see dokumentation of mxCreateSparse and the mxGet functions to understand this.
    void mapToMxArray(mxArray * & in){

        int len = data.size();
        in = mxCreateSparse(width, length, len, mxREAL);
        int* jc = mxGetJc(in);
        int* ir = mxGetIr(in);
        double* pr = mxGetPr(in);
        if(len == 0){
            jc[0] = 1;
            return;
        }
        typename std::map<TinyVector<int,2>, T, ShapeCmp>::iterator iter;
        TinyVector<int,2> newShape;
        int ii = 0;
        int jj = 0;
        int curjc = -1;
        for( iter = data.begin(); iter != data.end(); ++iter ) {
            newShape = iter->first;
            ir[ii] = newShape[1];
            pr[ii] = iter->second;
            if(newShape[0]  != curjc){
                curjc = newShape[0] ;
                jc[jj] = ii;
                jj++;
            }

            ii++;
        }
        jc[jj] = len;
    }

};

enum DataDimension {IMAGE = 2, VOLUME = 3};

} // namespace matlab

} // namespace vigra

void vigraMexFunction(vigra::matlab::OutputArray, vigra::matlab::InputArray);

#ifndef VIGRA_CUSTOM_MEXFUNCTION

/*
    DO NOT Comment out this function. If you are using a
    custom mexfunction just #define VIGRA_CUSTOM_MEXFUNCTION
    before #including matlab.hxx.
*/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  try
  {
    vigra::matlab::InputArray inputs(nrhs, prhs);
    vigra::matlab::OutputArray outputs(nlhs, plhs);

    vigraMexFunction(outputs, inputs);
  }
  catch(std::exception & e)
  {
    mexErrMsgTxt(e.what());
  }
}

#endif /*CUSTOM_MEXFUNCTION*/


#define VIGRA_CREATE_ENUM_AND_STD_MAP2(mapName, item1, item2) \
    const int item1 = 1;\
    const int item2 = 2;\
    std::map<std::string,int>  mapName;\
    mapName[#item1] = (int)item1;\
    mapName[#item2] = (int)item2;\


#define VIGRA_CREATE_ENUM_AND_STD_MAP3(mapName, item1, item2, item3) \
    const int item1 = 1;\
    const int item2 = 2;\
    const int item3 = 3;\
    std::map<std::string,int>  mapName;\
    mapName[#item1] = (int)item1;\
    mapName[#item2] = (int)item2;\
    mapName[#item3] = (int)item3;\


#define VIGRA_CREATE_ENUM_AND_STD_MAP4(mapName, item1, item2, item3, item4) \
    const int item1 = 1;\
    const int item2 = 2;\
    const int item3 = 3;\
    const int item4 = 4;\
    std::map<std::string,int>  mapName;\
    mapName[#item1] = (int)item1;\
    mapName[#item2] = (int)item2;\
    mapName[#item3] = (int)item3;\
    mapName[#item4] = (int)item4;\

#define VIGRA_CREATE_ENUM_AND_STD_MAP5(mapName, item1, item2, item3, item4, item5) \
    const int item1 = 1;\
    const int item2 = 2;\
    const int item3 = 3;\
    const int item4 = 4;\
    const int item5 = 5;\
    std::map<std::string, int>  mapName;\
    mapName[#item1] = (int)item1;\
    mapName[#item2] = (int)item2;\
    mapName[#item3] = (int)item3;\
    mapName[#item4] = (int)item4;\
    mapName[#item5] = (int)item5;\

#define VIGRA_CREATE_ENUM_AND_STD_MAP6(mapName, item1, item2, item3, item4, item5, item6) \
    const int item1 = 1;\
    const int item2 = 2;\
    const int item3 = 3;\
    const int item4 = 4;\
    const int item5 = 5;\
    const int item6 = 6;\
    std::map<std::string,int>  mapName;\
    mapName[#item1] = (int)item1;\
    mapName[#item2] = (int)item2;\
    mapName[#item3] = (int)item3;\
    mapName[#item4] = (int)item4;\
    mapName[#item5] = (int)item5;\
    mapName[#item6] = (int)item6;\

#endif // VIGRA_MATLAB_HXX
