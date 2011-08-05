#ifndef VIGRANUMPYIMPEX_HXX
#define VIGRANUMPYIMPEX_HXX

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <numpy/arrayobject.h>

#include <iostream>
#include "vigra/impex.hxx"

typedef double real_value;

using namespace boost::python;

typedef vigra::MultiArrayView< 3, vigra::UInt64, vigra::UnstridedArrayTag> UMArray3Int64;
typedef vigra::MultiArrayView< 3, vigra::UInt32, vigra::UnstridedArrayTag> UMArray3Int32;
typedef vigra::MultiArrayView< 3, vigra::UInt16, vigra::UnstridedArrayTag> UMArray3Int16;

typedef vigra::MultiArrayView< 4, vigra::UInt64, vigra::UnstridedArrayTag> UMArray4Int64;
typedef vigra::MultiArrayView< 4, vigra::UInt32, vigra::UnstridedArrayTag> UMArray4Int32;
typedef vigra::MultiArrayView< 4, vigra::UInt16, vigra::UnstridedArrayTag> UMArray4Int16;
typedef vigra::MultiArrayView< 4, double, vigra::UnstridedArrayTag> UMArray4FLOAT64;

typedef vigra::MultiArrayView< 1, double, vigra::UnstridedArrayTag> UMArray1FLOAT64;
typedef vigra::MultiArrayView< 1, vigra::UInt32, vigra::UnstridedArrayTag> UMArray1Int64;
typedef vigra::MultiArrayView< 1, vigra::UInt64, vigra::UnstridedArrayTag> UMArray1Int32;

typedef vigra::MultiArrayView< 4, double, vigra::UnstridedArrayTag> UMArray4Float64;

typedef vigra::MultiArrayView< 2, double, vigra::UnstridedArrayTag> UMArray2FLOAT64;
typedef vigra::MultiArrayView< 2, vigra::UInt32, vigra::UnstridedArrayTag> UMArray2UInt32;
typedef vigra::MultiArrayView< 2, vigra::UInt64, vigra::UnstridedArrayTag> UMArray2UInt64;

//typedef vigra::MultiArrayView< 3, vigra::TinyVector<3,double>, vigra::UnstridedArrayTag> UMArray3V3FLOAT64;


void exportVigraNumpyImpex();

namespace vigra
{

  PyObject* createNumpyArray2D(Size2D const & shape, 
                   int pyArrayTypeConstant, VigraTrueType /*gray volume*/);
  
  PyObject* createNumpyArray2D(Size2D const & shape, 
                   int pyArrayTypeConstant, VigraFalseType /*color volume*/);
  
  template<typename ValueType, int dim>
  PyObject* createNumpyArray(TinyVector<ValueType,dim> const & shape, 
                 int pyArrayTypeConstant);

  PyObject* createNumpyArray(int shape, int pyArrayTypeConstant);

enum VigraNumpyExportPixelType { ExportDEFAULT, ExportUINT8, ExportINT16, 
    ExportUINT16, ExportINT32, ExportUINT32, ExportFLOAT32, ExportFLOAT64, 
    ExportNUINT8, ExportNINT16, ExportNUINT16, ExportNINT32, ExportNUINT32, 
    ExportNFLOAT32, ExportNFLOAT64};

enum ImportPixelType { UINT8, INT16, UINT16, INT32, UINT32, FLOAT, DOUBLE };

template <class T>
void setRangeMapping(std::string const & pixeltype, 
    FindMinMax<T> const & minmax, ImageExportInfo & info);

// template<>
// struct NumericTraits<npy_int64>
// {
//     typedef npy_int64 Type;
//     typedef npy_int64 Promote;
//     typedef npy_float64 RealPromote;
//     typedef std::complex<RealPromote> ComplexPromote;
//     typedef Type ValueType;

//     typedef VigraTrueType isIntegral;
//     typedef VigraTrueType isScalar;
//     typedef VigraTrueType isSigned;
//     typedef VigraTrueType isOrdered;
//     typedef VigraFalseType isComplex;
    
//     static npy_int64 zero() { return 0; }
//     static npy_int64 one() { return 1; }
//     static npy_int64 nonZero() { return 1; }
//     static npy_int64 min() { return NPY_MIN_INT64; }
//     static npy_int64 max() { return NPY_MAX_INT64; }
    
// #ifdef NO_INLINE_STATIC_CONST_DEFINITION
//     enum { minConst = NPY_MIN_INT64, maxConst = NPY_MAX_INT64 };
// #else
//     static const npy_int64 minConst = NPY_MIN_INT64;
//     static const npy_int64 maxConst = NPY_MAX_INT64;
// #endif

//     static Promote toPromote(npy_int64 v) { return v; }
//     static RealPromote toRealPromote(npy_int64 v) { return v; }
//     static npy_int64 fromPromote(Promote v) { return v; }
//     static npy_int64 fromRealPromote(RealPromote v) {
//         return ((v < 0.0) 
//                  ? ((v < (RealPromote)NPY_MIN_INT64) 
//                      ? NPY_MIN_INT64
//                      : static_cast<npy_int64>(v - 0.5)) 
//                  : ((v > (RealPromote)NPY_MAX_INT64) 
//                      ? NPY_MAX_INT64 
//                      : static_cast<npy_int64>(v + 0.5))); 
//     }
// };

// template<>
// struct NumericTraits<npy_uint64>
// {
//     typedef npy_uint64 Type;
//     typedef npy_uint64 Promote;
//     typedef npy_float64 RealPromote;
//     typedef std::complex<RealPromote> ComplexPromote;
//     typedef Type ValueType;

//     typedef VigraTrueType isIntegral;
//     typedef VigraTrueType isScalar;
//     typedef VigraFalseType isSigned;
//     typedef VigraTrueType isOrdered;
//     typedef VigraFalseType isComplex;
    
//     static npy_uint64 zero() { return 0; }
//     static npy_uint64 one() { return 1; }
//     static npy_uint64 nonZero() { return 1; }
//     static npy_uint64 min() { return 0; }
//     static npy_uint64 max() { return NPY_MAX_UINT64; }
    
// #ifdef NO_INLINE_STATIC_CONST_DEFINITION
//     enum { minConst = 0, maxConst = NPY_MAX_UINT64 };
// #else
//     static const npy_uint64 minConst = 0;
//     static const npy_uint64 maxConst = NPY_MAX_UINT64;
// #endif

//     static Promote toPromote(npy_uint64 v) { return v; }
//     static RealPromote toRealPromote(npy_uint64 v) { return v; }
//     static npy_uint64 fromPromote(Promote v) { return v; }
//     static npy_uint64 fromRealPromote(RealPromote v) {
//             return ((v < 0.0) 
//                      ? 0 
//                      : ((v > (RealPromote)NPY_MAX_UINT64) 
//                          ? NPY_MAX_UINT64
//                          : static_cast<npy_uint64>(v + 0.5)));
//     }
// };
} // namespace vigra
#endif
