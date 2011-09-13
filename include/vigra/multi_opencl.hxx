/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
/*               Copyright 2011-2011 by Michael Tesch                   */
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

#ifndef VIGRA_OPENCL_HXX
#define VIGRA_OPENCL_HXX

#include "numerictraits.hxx"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

namespace vigra {

/********************************************************/
/*                                                      */
/*                      NumericTraits                   */
/*                                                      */
/********************************************************/

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

#define VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(basetype, n)               \
  template<>                                                            \
  struct NumericTraits< basetype##n >                                   \
  {                                                                     \
    typedef basetype##n Type;                                           \
    typedef Type Promote;                                               \
    typedef Type UnsignedPromote;                                       \
    typedef Type RealPromote;                                           \
    typedef std::complex<Type> ComplexPromote;                          \
    typedef basetype ValueType;                                         \
                                                                        \
    typedef VigraFalseType isIntegral;                                  \
    typedef VigraFalseType isScalar;                                    \
    typedef typename NumericTraits<ValueType>::isSigned isSigned;       \
    typedef VigraFalseType isOrdered;                                   \
    typedef typename NumericTraits<ValueType>::isComplex isComplex;     \
                                                                        \
    static Type zero() { Type x; bzero(&x, sizeof(x)); return x; }      \
    static Type one() { Type x = {{1}}; return x; }                     \
    static Type nonZero() { return one(); }                             \
                                                                        \
    static Promote toPromote(Type const & v) { return v; }              \
    static Type fromPromote(Promote const & v) { return v; }            \
    static Type fromRealPromote(RealPromote v) { return v; }            \
  }

#define VIGRA_OPENCL_VECTYPEN_REAL_TRAITS(basetype, n)                  \
  template<>                                                            \
  struct NumericTraits< basetype##n >                                   \
  {                                                                     \
    typedef basetype##n Type;                                           \
    typedef Type Promote;                                               \
    typedef Type UnsignedPromote;                                       \
    typedef Type RealPromote;                                           \
    typedef std::complex<Type> ComplexPromote;                          \
    typedef basetype ValueType;                                         \
                                                                        \
    typedef VigraFalseType isIntegral;                                  \
    typedef VigraFalseType isScalar;                                    \
    typedef typename NumericTraits<ValueType>::isSigned isSigned;       \
    typedef VigraFalseType isOrdered;                                   \
    typedef typename NumericTraits<ValueType>::isComplex isComplex;     \
                                                                        \
    static Type zero() { Type x; bzero(&x, sizeof(x)); return x; }      \
    static Type one() { Type x = {{1}}; return x; }                     \
    static Type nonZero() { return one(); }                             \
    static Type epsilon() { Type x; x.x = NumericTraits<ValueType>::epsilon(); return x; } \
    static Type smallestPositive() { Type x; x.x = NumericTraits<ValueType>::smallestPositive(); return x; } \
                                                                        \
    static Promote toPromote(Type const & v) { return v; }              \
    static Type fromPromote(Promote const & v) { return v; }            \
    static Type fromRealPromote(RealPromote v) { return v; }            \
  }

/// \todo - fix one() - maybe with .hi and .lo accessors?

#define VIGRA_OPENCL_VECN_TRAITS(n)                          \
  VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_char, n);          \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_uchar, n);       \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_short, n);       \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_ushort, n);      \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_int, n);         \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_uint, n);        \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_long, n);        \
    VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS(cl_ulong, n);       \
    VIGRA_OPENCL_VECTYPEN_REAL_TRAITS(cl_float, n);          \
    VIGRA_OPENCL_VECTYPEN_REAL_TRAITS(cl_double, n);

VIGRA_OPENCL_VECN_TRAITS(2);
VIGRA_OPENCL_VECN_TRAITS(3);
//VIGRA_OPENCL_VECN_TRAITS(4); // cl_type4 is the same as cl_type3
VIGRA_OPENCL_VECN_TRAITS(8);
VIGRA_OPENCL_VECN_TRAITS(16);

#undef VIGRA_OPENCL_VECTYPEN_INTEGER_TRAITS
#undef VIGRA_OPENCL_VECTYPEN_REAL_TRAITS
#undef VIGRA_OPENCL_VECN_TRAITS

/** \todo looks like the windows CL/cl_platform.h does signed/unsigned
 *      strangely, so that the signed properties may not have the right
 *      NumericalTraits::isSigned -- not sure if there's a reason for that.
 */


#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                    SquareRootTraits                  */
/*                                                      */
/********************************************************/

/********************************************************/
/*                                                      */
/*                       NormTraits                     */
/*                                                      */
/********************************************************/

#if 0
template<>
struct NormTraits<fftw_complex>
{
    typedef fftw_complex Type;
    typedef fftw_real    SquaredNormType;
    typedef fftw_real    NormType;
};

template<class Real>
struct NormTraits<FFTWComplex<Real> >
{
    typedef FFTWComplex<Real>  Type;
    typedef typename Type::SquaredNormType   SquaredNormType;
    typedef typename Type::NormType   NormType;
};
#endif

/********************************************************/
/*                                                      */
/*                      PromoteTraits                   */
/*                                                      */
/********************************************************/

#if 0
template<class T>
struct CanSkipInitialization<std::complex<T> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};
#endif

/********************************************************/
/*                                                      */
/*                      multi_math                      */
/*                                                      */
/********************************************************/

namespace multi_math {

/// \todo !

/**     OpenCL 1.1 [6.2] - Convert operators */
/**     OpenCL 1.1 [6.3] - Scalar/vector math operators */

/**     OpenCL 1.1 [6.11.2] - Math Built-in Functions */
/**     OpenCL 1.1 [6.11.3] - Integer Built-in Functions */
/**     OpenCL 1.1 [6.11.4] - Common Built-in Functions */
/**     OpenCL 1.1 [6.11.5] - Geometric Built-in Functions */
/**     OpenCL 1.1 [6.11.6] - Relational Built-in Functions */
/**     OpenCL 1.1 [6.11.7] - Vector Data Load/Store Built-in Functions */

/**     OpenCL 1.1 [6.11.12] - Misc Vector Built-in Functions */

/**     OpenCL 1.1 [6.11.12] - Image Read and Write Built-in Functions */


} // namespace multi_math

/********************************************************/
/*                                                      */
/*                  Channel Accessors                   */
/*                                                      */
/********************************************************/

/** \addtogroup DataAccessors
*/
//@{
/** \defgroup OpenCL-Accessors Accessors for OpenCL types

    Encapsulate access to members of OpenCL vector types.


    <b>\#include</b> \<vigra/multi_opencl.hxx\>

    OpenCL 1.1 [6.1.7] - Vector Components

    - cl_TYPE2Accessor_x
    - cl_TYPE2Accessor_y
    - cl_TYPE2Accessor_s0
    - cl_TYPE2Accessor_s1

    - cl_TYPE2WriteAccessor_x
    - cl_TYPE2WriteAccessor_y
    - cl_TYPE2WriteAccessor_s0
    - cl_TYPE2WriteAccessor_s1

    - cl_TYPE3Accessor_x
    - cl_TYPE3Accessor_y
    - cl_TYPE3Accessor_z
    - cl_TYPE3Accessor_s0
    - cl_TYPE3Accessor_s1
    - cl_TYPE3Accessor_s2

    - cl_TYPE3WriteAccessor_x
    - cl_TYPE3WriteAccessor_y
    - cl_TYPE3WriteAccessor_z
    - cl_TYPE3WriteAccessor_s0
    - cl_TYPE3WriteAccessor_s1
    - ...

    where TYPE is one of {char, uchar, short, ushort, int, uint, long, ulong, float, double }

    For example:

    \code

    #include <vigra/multi_opencl.hxx>

    MultiArrayView<2, cl_double3 > dataView = ...;

    vigra::FindMinMax<double> minmax;
    vigra::inspectMultiArray(srcMultiArrayRange(dataView, cl_double3Accessor_z()), minmax);
    std::cout << "range of .z: " << minmax.min << " - " << minmax.max;

    \endcode
*/
//@{
/**
   \class cl_charNAccessor_COMP
   
   access the first component.

   \class cl_TYPE3WriteAccessor_s1

   access the second component.

   \class cl_TYPE3WriteAccessor_s2

   access the third component.
*/

//@}
//@}

#define VIGRA_OPENCL_TYPE_ACCESSOR(basetype, n, NTH) \
  class basetype##n##Accessor_##NTH                  \
  {                                                             \
  public:                                                       \
    /** The accessor's value type. */                           \
    typedef NumericTraits< basetype##n >::ValueType value_type; \
                                                                \
    /** Read component at iterator position. */                 \
    template <class ITERATOR>                                   \
      value_type operator()(ITERATOR const & i) const {         \
      return (*i).NTH;                                          \
    }                                                           \
                                                                \
    /** Read component at offset from iterator position. */             \
    template <class ITERATOR, class DIFFERENCE>                         \
      value_type operator()(ITERATOR const & i, DIFFERENCE d) const {   \
      return i[d].NTH;                                                  \
    }                                                                   \
                                                                        \
    /** Write component at iterator position from a scalar. */          \
    template <class ITERATOR>                                           \
      void set(value_type const & v, ITERATOR const & i) const {        \
      (*i).NTH = v;                                                     \
    }                                                                   \
                                                                        \
    /** Write component at offset from iterator position from a scalar. */ \
    template <class ITERATOR, class DIFFERENCE>                         \
      void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const { \
      i[d].NTH = v;                                                     \
    }                                                                   \
                                                                        \
    /** Write component at iterator position into a scalar. */          \
    template <class R, class ITERATOR>                                  \
      void set(FFTWComplex<R> const & v, ITERATOR const & i) const {    \
      *i = v.NTH;                                                       \
    }                                                                   \
                                                                        \
    /** Write component at offset from iterator position into a scalar. */ \
    template <class R, class ITERATOR, class DIFFERENCE>                \
      void set(FFTWComplex<R> const & v, ITERATOR const & i, DIFFERENCE d) const { \
      i[d] = v.NTH;                                                     \
    }                                                                   \
  };                                                                    \
  class basetype##n##WriteAccessor_##NTH                                \
    : public basetype##n##Accessor_##NTH                                \
  {                                                                     \
  public:                                                               \
    /** The accessor's value type. */                                   \
    typedef NumericTraits< basetype##n >::ValueType value_type;         \
                                                                        \
    /** Write component at iterator position. */                        \
    template <class ITERATOR>                                           \
    void set(value_type const & v, ITERATOR const & i) const {          \
      (*i).NTH = v;                                                     \
    }                                                                   \
                                                                        \
    /** Write component at offset from iterator position. */            \
    template <class ITERATOR, class DIFFERENCE>                         \
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const { \
      i[d].NTH = v;                                                     \
    }                                                                   \
  }

#define VIGRA_OPENCL_TYPE2_ACCESSORS(basetype)  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 2, s0);  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 2, s1);  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 2, x);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 2, y);

#define VIGRA_OPENCL_TYPE3_ACCESSORS(basetype)  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, s0);  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, s1);  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, s2);  \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, x);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, y);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 3, z);

#define VIGRA_OPENCL_TYPE4_ACCESSORS(basetype)   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, s0);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, s1);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, s2);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, s3);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, x);    \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, y);    \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, z);    \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 4, w);

#define VIGRA_OPENCL_TYPE8_ACCESSORS(basetype)   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s0);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s1);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s2);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s3);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s4);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s5);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s6);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s7);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 8, s8);

#define VIGRA_OPENCL_TYPE16_ACCESSORS(basetype)   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s0);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s1);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s2);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s3);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s4);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s5);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s6);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s7);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, s8);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sa);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sb);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sc);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sd);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, se);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sf);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sA);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sB);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sC);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sD);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sE);   \
  VIGRA_OPENCL_TYPE_ACCESSOR(basetype, 16, sF);

/// \todo figure out half (.hi .lo, .even .odd) and other odd-sized accessors

#define VIGRA_OPENCL_ACCESSORS(basetype)  \
  VIGRA_OPENCL_TYPE2_ACCESSORS(basetype); \
  VIGRA_OPENCL_TYPE3_ACCESSORS(basetype); \
  VIGRA_OPENCL_TYPE4_ACCESSORS(basetype); \
  VIGRA_OPENCL_TYPE8_ACCESSORS(basetype); \
  VIGRA_OPENCL_TYPE16_ACCESSORS(basetype);

VIGRA_OPENCL_ACCESSORS(cl_char);
VIGRA_OPENCL_ACCESSORS(cl_uchar);
VIGRA_OPENCL_ACCESSORS(cl_short);
VIGRA_OPENCL_ACCESSORS(cl_ushort);
VIGRA_OPENCL_ACCESSORS(cl_int);
VIGRA_OPENCL_ACCESSORS(cl_uint);
VIGRA_OPENCL_ACCESSORS(cl_long);
VIGRA_OPENCL_ACCESSORS(cl_ulong);
VIGRA_OPENCL_ACCESSORS(cl_float);
VIGRA_OPENCL_ACCESSORS(cl_double);

#undef VIGRA_OPENCL_TYPE_ACCESSOR
#undef VIGRA_OPENCL_TYPE2_ACCESSORS
#undef VIGRA_OPENCL_TYPE3_ACCESSORS
#undef VIGRA_OPENCL_TYPE4_ACCESSORS
#undef VIGRA_OPENCL_TYPE8_ACCESSORS
#undef VIGRA_OPENCL_TYPE16_ACCESSORS
#undef VIGRA_OPENCL_ACCESSORS

} // namespace vigra

#endif // VIGRA_OPENCL_HXX
