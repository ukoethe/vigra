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

// TODO: 
//    * handle rgb images
//    * handle complex matrices
//    * handle sparse matrices

class InputArray
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

    InputArray(size_type size, pointer data)
    : size_(size),
      data_(data)
    {}
      
    const_reference operator[]( difference_type i ) const
    {
        if(!isValid(i))
            mexErrMsgTxt("Too few input arguments.");
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
};

class OutputArray
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

template <unsigned int SIZE, class T>
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
TinyVectorView<MultiArrayIndex, SIZE>
getShape(mxArray const * t)
{
    if(!ValueType<MultiArrayIndex>::check(t))
    {
        std::string msg = std::string("Input array must have type ") + 
                          ValueType<MultiArrayIndex>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getShape(): Input array has wrong number of elements.");
    }
    
    return TinyVectorView<MultiArrayIndex, SIZE>((MultiArrayIndex *)mxGetData(t));
}

template <class T, unsigned int SIZE>
TinyVectorView<T, SIZE>
getVector(mxArray const * t)
{
    if(!ValueType<T>::check(t))
    {
        std::string msg = std::string("Input array must have type ") + 
                          ValueType<T>::typeName() + ".";
        mexErrMsgTxt(msg.c_str());
    }
    if(SIZE != mxGetNumberOfElements(t))
    {
        mexErrMsgTxt("getVector(): Input array has wrong number of elements.");
    }
    
    return TinyVectorView<T, SIZE>((T *)mxGetData(t));
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

template<class T>
T
getScalar(mxArray const * t)
{
    if(mxIsEmpty(t))
        mexErrMsgTxt("getScalar() on empty input.");
    return static_cast<T>(mxGetScalar(t));
}


std::string getString(mxArray const * t)
{
    if(mxIsEmpty(t))
        mexErrMsgTxt("getString() on empty input.");
	char temp_save[30];
	int sze =  mxGetN(t);

	mxGetString(t, temp_save, sze+1);
	
	std::string out = temp_save;
	return out;
}

template<class T>
mxArray *
createScalar(T v)
{
    mxArray * m;
    createMatrix<double>(1, 1, m)(0,0) = static_cast<double>(v);
    return m;
}

/***********************************
Rahuls code starts here
************************************+*/
using namespace vigra;

template <class Functor>
void callMexFunctor(matlab::OutputArray outputs, matlab::InputArray inputs)
{
	mxClassID inClass = mxGetClassID(inputs[0]);
	switch(inClass){
		case mxDOUBLE_CLASS:
			Functor::template exec<double>(outputs, inputs);	break;
		case mxSINGLE_CLASS:
			Functor::template exec<float>(outputs, inputs);		break;
        case mxINT8_CLASS:
			Functor::template exec<Int8>(outputs, inputs);		break;
		case mxINT16_CLASS:
			Functor::template exec<Int16>(outputs, inputs);		break;
		case mxINT32_CLASS:
			Functor::template exec<Int32>(outputs, inputs);		break;
		case mxINT64_CLASS:
			Functor::template exec<Int64>(outputs, inputs);		break;
        case mxUINT8_CLASS:
			Functor::template exec<UInt8>(outputs, inputs);		break;
		case mxUINT16_CLASS:
			Functor::template exec<UInt16>(outputs, inputs);	break;
		case mxUINT32_CLASS:
			Functor::template exec<UInt32>(outputs, inputs);	break;
		case mxUINT64_CLASS:
			Functor::template exec<UInt64>(outputs, inputs);	break;
		default:
			mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");
	}
}

template <class Functor>
void callMexAllIntFunctor(matlab::OutputArray outputs, matlab::InputArray inputs)
{
	mxClassID inClass = mxGetClassID(inputs[0]);
	switch(inClass){
        case mxINT8_CLASS:
			Functor::template exec<Int8>(outputs, inputs);		break;
		case mxINT16_CLASS:
			Functor::template exec<Int16>(outputs, inputs);		break;
		case mxINT32_CLASS:
			Functor::template exec<Int32>(outputs, inputs);		break;
		case mxINT64_CLASS:
			Functor::template exec<Int64>(outputs, inputs);		break;
        case mxUINT8_CLASS:
			Functor::template exec<UInt8>(outputs, inputs);		break;
		case mxUINT16_CLASS:
			Functor::template exec<UInt16>(outputs, inputs);	break;
		case mxUINT32_CLASS:
			Functor::template exec<UInt32>(outputs, inputs);	break;
		case mxUINT64_CLASS:
			Functor::template exec<UInt64>(outputs, inputs);	break;
		default:
			mexErrMsgTxt("Input image must have type 'uint8'-16-32-64', 'int8-16-32-64' 'single' or 'double'.");
	}
}
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
bool is_in_range(T in, T min, T max)
{
	return (in >= min && in <= max)? true:false;
}
template<class T>
bool is_in_range(T in, std::string min, T max)
{
	return(in <= max)? true:false;
}

template<class T>
bool is_in_range(T in, T min, std::string max)
{
	return (in >= min)? true:false;
}
/*++++++++++++++++++++++++++HARDCORE-MACROS++++++++++++++++++++++++++++*/

//Definition of Membervariables and assignment functions
#define declScalar(type, name, default); \
    type name;\
	type get_##name(matlab::InputArray inputs)\
	{\
		if(inputs.size() == 1){\
			return default;\
		}else{\
			mxArray* name =mxGetField(inputs[1], 0, #name);\
			return (name!=NULL&&mxIsNumeric(name))\
							? matlab::getScalar<type>(name)\
							: default;\
		}\
	}\

#define declScalar2D3D(type, name, default2D, default3D); \
    type name;\
	type get_##name(matlab::InputArray inputs)\
	{\
		if(inputs.size() == 1){\
			return this->numOfDim == 2? default2D : default3D;;\
		}else{\
			mxArray* name =mxGetField(inputs[1], 0, #name);\
			return (name!=NULL&&mxIsNumeric(name))\
							? matlab::getScalar<type>(name)\
							: this->numOfDim == 2? default2D : default3D;\
		}\
	}\
	
#define declScalarMinMax(type, name, default,min, max); \
	type name;\
	type get_##name(matlab::InputArray inputs)\
	{\
		if(inputs.size() == 1){\
			return default;\
		}else{\
			mxArray* name =mxGetField(inputs[1], 0, #name);\
			type out = (name!=NULL&&mxIsNumeric(name))\
								? matlab::getScalar<type>(name)\
							 	: default;\
			if(!is_in_range<type>(out,min,max)){\
				out = default;\
				mexWarnMsgTxt("Value of '"#name"' out of range ["#min" ... "#max"]. using default: " #default );\
			}\
			return out;\
		}\
	}

#define declCharConstr(title, number, name1_default, name2, name3, name4, name5);\
	enum title##Enum {name1_default = 1, name2 = 2, name3 = 3, name4 = 4, name5 = 5};\
	int title;\
	int get_##title(matlab::InputArray inputs)\
	{\
		std::string title##_str[5];\
		title##_str[0]= #name1_default;\
		title##_str[1]= #name2;\
		title##_str[2]= #name3;\
		title##_str[3]= #name4;\
		title##_str[4]= #name5;\
		int returnval = name1_default;\
		if(inputs.size() == 2){\
			mxArray* title = mxGetField(inputs[1], 0 , #title);\
			std::string test;\
			bool flag = false;\
			if(title != NULL && mxIsChar(title)){\
				flag = true;\
				test = matlab::getString(title);\
				for(int ii = 0; ii < number;ii++)\
				{\
					if(test == title##_str[ii])\
					{\
						returnval = ii+1;\
						flag = false;\
					}\
				}\
			}\
			if(flag)mexWarnMsgTxt("Value of '" #title "' not Valid. Using default: '" #name1_default "'");\
		}\
		return returnval;\
	}

#define declImage(name, type);\
		BasicImageView<type>  name;\
		BasicImageView<type> get_##name(mxArray *opt)\
		{\
			mxArray* nameArr =mxGetField(opt, 0, #name);\
			BasicImageView<type> bla;\
			return (name!=NULL&&mxIsNumeric(name))\
		                 ? matlab::getImage<type>(name)\
						 : bla;\
		}

#define declMultiArray(name, type);\
		MultiArrayView<3,type>  name;\
		MultiArrayView<3,type> get_##name(mxArray *opt)\
		{\
			mxArray* nameArr =mxGetField(opt, 0, #name);\
			MultiArrayView<3,type> bla;\
			return (name!=NULL&&mxIsNumeric(name))\
		                 ? matlab::getMultiArray<3, T><type>(name)\
						 : bla;\
		}

#define declOut(type);\
	BasicImageView<type>  out;\
	MultiArrayView<3,type> out3D;
    
// Some Macros for commonly used output declarations
#define mapOut_SAME(outType);\
	if(this->numOfDim == 2){\
		out = matlab::createImage<outType>(this->in.width(), this->in.height(), outputs[0]);\
		out3D = MultiArrayView<3, outType>(this->in3D.shape(), (outType*)out.data());\
	}else{\
		out3D = matlab::createMultiArray<3,outType>(this->in3D.shape(), outputs[0]);\
	}	
#define mapOut_SIZE(outType,w,h,d);\
	MultiArrayShape<3>::type newShape(w, h, d);\
	if(numOfDim == IMAG){\
		out = matlab::createImage<outType>(w, h, outputs[0]);\
		out3D = MultiArrayView<3, outType>(newShape, (outType*)out.data());\
	}else{\
		out3D = matlab::createMultiArray<3,outType>(newShape, outputs[0]);\
	}	
#define mapOut_2D(outType,w,h,d);\
	MultiArrayShape<3>::type newShape(w, h, 1);\
	out = matlab::createImage<outType>(w, h, outputs[0]);\
	out3D = MultiArrayView<3, outType>(newShape, (outType*)out.data());\

#define mapOut_3D(outType,w,h,d)\
	MultiArrayShape<3>::type scrabawaduwada(w, h, d);\
	out3D = matlab::createMultiArray<3,outType>(scrabawaduwada, outputs[0]);\

//Simplify Member Initialisors 
#define map(name) name(get_##name(inputs))

//enumeration - for GCC 
enum DataDimension {IMAG = 2, VOLUME = 3};

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* The Optons struct contains all the necassary working data and 
/* options for the vigraFunc. This is the minimal struct
/* Add fields as necessary
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct base_data{
	int numOfDim;
//	enum {IMAG = 2, VOLUME = 3} dim;
	BasicImageView<T>  in;
	MultiArrayView<3,T> in3D;

	base_data(matlab::InputArray inputs)
	{
		numOfDim = mxGetNumberOfDimensions(inputs[0]);
		if(numOfDim > 3)  
							mexErrMsgTxt("Currently InputArray may only have 2 or 3 dimensions");
		if(inputs.isEmpty(0)) 
							mexErrMsgTxt("Input Image is empty!");
		if(numOfDim == 2){
			in = matlab::getImage<T>(inputs[0]);
			in3D = matlab::getMultiArray<3, T>(inputs[0]);
		}else{
			in3D = matlab::getMultiArray<3, T>(inputs[0]);
		}
	}
	
};


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
class SparseArray{
	
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

#endif // VIGRA_MATLAB_HXX
