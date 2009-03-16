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
#include <sstream>
#include <iostream>

namespace vigra {

typedef enum {
        vUNKNOWN,
        vCELL,
        vSTRUCT,
        vLOGICAL,
        vCHAR,
        vDOUBLE,
        vSINGLE,
        vINT8,
        vUINT8,
        vINT16,
        vUINT16,
        vINT32,
        vUINT32,
        vINT64,
        vUINT64,
        vFUNCTION
} vClassID;

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
VIGRA_MATLAB_VALUETYPE_UTIL(Int16, Int16, mxINT16_CLASS, int16)
VIGRA_MATLAB_VALUETYPE_UTIL(Int32, Int32, mxINT32_CLASS, int32)
VIGRA_MATLAB_VALUETYPE_UTIL(Int64, Int64, mxINT64_CLASS, int64)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt8,  Uint8, mxUINT8_CLASS, uint8)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt16, Uint16, mxUINT16_CLASS, uint16)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt32, Uint32, mxUINT32_CLASS, uint32)
VIGRA_MATLAB_VALUETYPE_UTIL(UInt64, Uint64, mxUINT64_CLASS, uint64)

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

    ConstCellArray(const mxArray * matPointer)
    : matPointer_(const_cast<mxArray *>(matPointer)),
      size_(0)
    {
        if(!mxIsCell(matPointer))
            mexErrMsgTxt("CellArray(mxArray *): Argument must be a Matlab cell array.");
        size_ = mxGetNumberOfElements(matPointer);
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
TinyVector<MultiArrayIndex, SIZE>
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
    return TinyVector<MultiArrayIndex, SIZE>(res);
}

// template <class T, unsigned int SIZE>
// TinyVectorView<T, SIZE>
// getVector(mxArray const * t)
// {
    // if(!ValueType<T>::check(t))
    // {
        // std::string msg = std::string("Input array must have type ") +
                          // ValueType<T>::typeName() + ".";
        // mexErrMsgTxt(msg.c_str());
    // }
    // if(SIZE != mxGetNumberOfElements(t))
    // {
        // mexErrMsgTxt("getVector(): Input array has wrong number of elements.");
    // }

    // return TinyVectorView<T, SIZE>((T *)mxGetData(t));
// }


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
        if(DIM < mdim)
        {
            mexErrMsgTxt("getMultiArray(): Input array has too many dimensions.");
        }
        const mwSize * matlabShape = mxGetDimensions(t);
        for(unsigned int k=0; k<mdim; ++k)
        {
            shape[k] = static_cast<typename Shape::value_type>(matlabShape[k]);
        }
        for(unsigned int k=mdim; k<DIM; ++k)
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
    for(int k=0; k<DIM; ++k)
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
    if(!mxIsNumeric(t))
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


class Required
{};

template<class T>
struct OptionalImpl
{
    T val;
    OptionalImpl(T v): val(v) {}

};

struct OptionalImplVoid
{
};

OptionalImplVoid Optional()
{
    return OptionalImplVoid();
}

template<class T>
OptionalImpl<T> Optional(T in)
{
    return OptionalImpl<T>(in);
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
        std::ostringstream s1;
        s1  << "Required Input In Optionstruct: '" << name \
            <<"' has not been supplied!"<< std::endl;
        return s1.str();
    }
    std::string createErrMsg(int pos)
    {
        std::ostringstream s1;
        s1  <<"Required Input In Signature of Function at position: '"\
            << pos << "' has not been supplied"<<std::endl;
        return s1.str();
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
      options_(this->isValid(size-1) && mxIsStruct((*this)[size-1])
                 ? (*this)[size-1]
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

    bool isEmpty(difference_type i){
        return mxIsEmpty(data_[i]);
    }

    bool isEmpty(std::string name){
        return mxIsEmpty(options_[name]);
    }

    template <class T, class U, class Place>
    T errorOrDefault(OptionalImpl<U> const & o, Place NameOrPos = Place()){ return o.val; }

    template <class T, class Place>
    T errorOrDefault(Required r, Place NameOrPos = Place())
    {
        mexErrMsgTxt(createErrMsg(NameOrPos).c_str());
        return T();
    }

    /*getter Func*/

    /*Scalar Type*/
    template <class T,class place, class ReqType>
    T getScalar(place posOrName, ReqType req)
    {
        if(!isValid(posOrName)|| isEmpty(posOrName))
        {
            return errorOrDefault<T>(req, posOrName);
        }
        else
        {
            return matlab::getScalar<T>((*this)[posOrName]);
        }
    }


    template <class T, class place, class reqClass, class minClass, class maxClass>
    T getScalarMinMax(place posOrName, reqClass req, minClass min_, maxClass max_)
    {
        T temp = this->getScalar<T>(posOrName, req);
        if (!is_in_range(temp, min_, max_))
            mexErrMsgTxt("Value Out of bound");
        else
            return temp;

    }

    template <class place, class reqClass>
    bool getBool(place posOrName, reqClass req)
    {
        return this->getScalarMinMax<bool>(posOrName, req, 0, 1);
    }

    /*Array Type*/
    template <unsigned int N, class T, class place, class reqClass>
    MultiArrayView<N,T> getMultiArray(place posOrName, reqClass req)
    {
        if(!isValid(posOrName) || isEmpty(posOrName))
        {
            return errorOrDefault< MultiArrayView<N,T> >(req, posOrName);
        }
        else
        {
            value_type temp = (*this)[posOrName];
            return matlab::getMultiArray<N,T>(temp);
        }
    }

    template < class T, class place, class reqClass>
    BasicImageView<T> getMultiArray(place posOrName, reqClass req)
    {
        if(!isValid(posOrName) || isEmpty(posOrName))
        {
            return errorOrDefault<BasicImageView<T> >(req, posOrName);
        }
        else
        {
            value_type temp = (*this)[posOrName];
            return matlab::getImage<T>(temp);
        }
    }

    template<class T,unsigned int sze, class place, class reqClass>
    TinyVectorView< T, sze> getTinyVector(place posOrName, reqClass req)
    {
        if(!isValid(posOrName) || isEmpty(posOrName))
        {
            return errorOrDefault<TinyVectorView< T, sze> >(req, posOrName);
        }
        else
        {
            value_type temp = (*this)[posOrName];
            return matlab::getTinyVector< T, sze>(temp);
        }
    }

    template< unsigned int sze, class place, class reqClass>
    TinyVectorView<MultiArrayIndex, sze> getShape(place posOrName, reqClass req)
    {
        if(!isValid(posOrName) || isEmpty(posOrName))
        {
            return errorOrDefault<TinyVectorView<MultiArrayIndex, sze> >(req, posOrName);
        }
        else
        {
            value_type temp = (*this)[posOrName];
            return matlab::getShape<sze>(temp);
        }
    }

    template<class place, class reqClass>
    ConstCellArray getCellArray(place posOrName, reqClass req)
    {
        if(!isValid(posOrName) || isEmpty(posOrName))
        {
            return errorOrDefault<ConstCellArray>(req, posOrName);
        }
        else
        {
            value_type temp = (*this)[posOrName];
            return matlab::getCellArray(temp);
        }
    }

};

class OutputArray
{
    int size_;
    mxArray ** data_;
    std::string createErrMsgOut(int pos)
    {
        std::ostringstream s1;
        s1  <<"Required Output at position: '"\
            << pos << "' has not been supplied"<<std::endl;
        return s1.str();
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
    T errorOrDefault(OptionalImplVoid const & o, int Pos){ return T(); }

    template <class T>
    T errorOrDefault(Required r, int Pos)
    {
        mexErrMsgTxt(createErrMsgOut(Pos).c_str());
        return T();
    }

    template <class T>
    T errorOrDefaultP(OptionalImplVoid const & o, int Pos){ return 0; }

    template <class T>
    T errorOrDefaultP(Required r, int Pos)
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
        return matlab::createMultiArray<DIM, T>(shape, (*this)[pos]);
    }

    template <class T, class ReqType>
    BasicImageView<T> createImage(int pos, ReqType req,
                                    mwSize width, mwSize height)
    {
        if(!isValid(pos))
            return errorOrDefault<BasicImageView<T> >(req, pos);
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
            return errorOrDefaultP<T>(req, pos);
        BasicImageView<T> temp = matlab::createImage<T>(1, 1, pos);
        return &temp[0];
    }

    template <class T, class ReqType>
    void createScalar(int pos, ReqType req, T val)
    {
        if(!isValid(pos))
            errorOrDefault<T>(req, pos);
        BasicImageView<T> temp = matlab::createImage<T>(1, 1, pos);
        temp[0] = val;
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


/*++++++++++++++++++++++++++HELPERFUNC+++++++++++++++++++++++++++++++*/
/* This is used for better readibility of the test cases            .
/* Nothing to be done here.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
    T& operator()(int i, int j){
        TinyVector<int,2> newShapew(i, j);
        typename std::map<TinyVector<int,2>, T, ShapeCmp>::iterator iter;
        TinyVector<int,2> newShape;
        return data[newShapew];
    }

    const T get(int i, int j){
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

} // namespace matlab

struct MeshGridAccessor
{
    typedef TinyVector<Diff2D::MoveX, 2> value_type;

    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const
    {
        return value_type(i->x, i->y);
    }
};

inline
triple<Diff2D, Diff2D, MeshGridAccessor>
meshGrid(Diff2D ul, Diff2D lr)
{
    return triple<Diff2D, Diff2D, MeshGridAccessor>(ul, lr, MeshGridAccessor());
}

} // namespace vigra

void vigraMexFunction(vigra::matlab::OutputArray, vigra::matlab::InputArray);

class vigraFunctor;

#ifndef DO_NOT_USE_VIGRA_MEX_FUNCTION

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

#endif /*DO_NOT_USE_VIGRA_MEX_FUNCTION*/

//Character Valued options
#define LOAD_ENUM_OPTION(title, number, name1_default, name2, name3, name4, name5)\
    title##Enum title;\
    {\
        std::map<std::string, int> title##_str;\
        if(number > 0) title##_str[#name1_default] = 1;\
        if(number > 1) title##_str[#name2] = 2;\
        if(number > 2) title##_str[#name3] = 3;\
        if(number > 3) title##_str[#name4] = 4;\
        if(number > 4) title##_str[#name5] = 5;\
        if(!inputs.options_.isValid(#title))\
        {\
            title =  name1_default;\
        }\
        else\
        {\
            const mxArray* name = inputs.options_[#title];\
            if(!mxIsChar(name))\
                mexErrMsgTxt("option '" #title "' must be a string.");\
            std::string namex = matlab::getString(inputs.options_[#title]); \
            if(title##_str.count(namex) == 0)\
                mexErrMsgTxt("option '" #title "' contains invalid string.");\
            title = (title##Enum)title##_str[namex];\
        }\
    };

#define LOAD_ENUM_OPTION2(title, name1_default, name2);\
    enum title##Enum {name1_default = 1, name2 = 2};\
    LOAD_ENUM_OPTION(title, 2, name1_default, name2, name3, name4, name5);

#define LOAD_ENUM_OPTION3(title, name1_default, name2, name3);\
    enum title##Enum {name1_default = 1, name2 = 2, name3 = 3};\
    LOAD_ENUM_OPTION(title, 3, name1_default, name2, name3, name4, name5);

#define LOAD_ENUM_OPTION4(title, name1_default, name2, name3, name4);\
    enum title##Enum {name1_default = 1, name2 = 2, name3 = 3, name4 = 4};\
    LOAD_ENUM_OPTION(title, 4, name1_default, name2, name3, name4, name5);

#define LOAD_ENUM_OPTION5(title, name1_default, name2, name3, name4, name5);\
    enum title##Enum {name1_default = 1, name2 = 2, name3 = 3, name4 = 4, name5 = 5};\
    LOAD_ENUM_OPTION(title, 5, name1_default, name2, name3, name4, name5);

#endif // VIGRA_MATLAB_HXX
