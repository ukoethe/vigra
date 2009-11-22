/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
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

#ifndef VIGRA_METAPROGRAMMING_HXX
#define VIGRA_METAPROGRAMMING_HXX

#include "config.hxx"
#include <climits>
#include <limits>

namespace vigra {

template <int N>
class MetaInt
{
  public:
    enum { value = N };
};

template <int N1, int N2>
class MetaMax
{
  public:
    enum { value = N1 < N2 ? N2 : N1 };
};

template <int N1, int N2>
class MetaMin
{
  public:
    enum { value = N1 < N2 ? N1 : N2 };
};

struct VigraTrueType
{
   enum { asBool = true };
};

struct VigraFalseType
{
    enum { asBool = false };
};

/**  \addtogroup MultiArrayTags Multi-dimensional Array Tags
      Meta-programming tags to mark array's as strided or unstrided.
*/

//@{

/********************************************************/
/*                                                      */
/*                   StridedArrayTag                    */
/*                                                      */
/********************************************************/

/** tag for marking a MultiArray strided.

<b>\#include</b>
\<<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>\>

Namespace: vigra
*/
struct StridedArrayTag {};

/********************************************************/
/*                                                      */
/*                  UnstridedArrayTag                   */
/*                                                      */
/********************************************************/

/** tag for marking a MultiArray unstrided.

<b>\#include</b>
\<<a href="multi__array_8hxx-source.html">vigra/multi_array.hxx</a>\>

Namespace: vigra
*/
struct UnstridedArrayTag {};

template<class T>
class TypeTraits
{
  public:
    typedef VigraFalseType isConst;
    typedef VigraFalseType isPOD;
    typedef VigraFalseType isBuiltinType;
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template<class T>
class TypeTraits<T const>
: public TypeTraits<T>
{
  public:
    typedef VigraTrueType isConst;
};

template<class T> 
class TypeTraits<T *>
{
  public:
    typedef VigraFalseType isConst;
    typedef VigraTrueType isPOD;
    typedef VigraTrueType isBuiltinType;
};

template<class T> 
class TypeTraits<T const *>
{
  public:
    typedef VigraFalseType isConst;
    typedef VigraTrueType isPOD;
    typedef VigraTrueType isBuiltinType;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

namespace detail {

template <int size>
struct SizeToType;

} // namespace detail 

#define VIGRA_TYPE_TRAITS(type, size) \
template<> \
class TypeTraits<type> \
{ \
  public: \
    typedef VigraFalseType isConst; \
    typedef VigraTrueType isPOD; \
    typedef VigraTrueType isBuiltinType; \
    typedef char TypeToSize[size]; \
}; \
 \
namespace detail { \
  TypeTraits<type>::TypeToSize * typeToSize(type); \
  \
  template <> \
  struct SizeToType<size> \
  { \
      typedef type result; \
  }; \
} 

VIGRA_TYPE_TRAITS(char, 1)
VIGRA_TYPE_TRAITS(signed char, 2)
VIGRA_TYPE_TRAITS(unsigned char, 3)
VIGRA_TYPE_TRAITS(short, 4)
VIGRA_TYPE_TRAITS(unsigned short, 5)
VIGRA_TYPE_TRAITS(int, 6)
VIGRA_TYPE_TRAITS(unsigned int, 7)
VIGRA_TYPE_TRAITS(long, 8)
VIGRA_TYPE_TRAITS(unsigned long, 9)
VIGRA_TYPE_TRAITS(float, 10)
VIGRA_TYPE_TRAITS(double, 11)
VIGRA_TYPE_TRAITS(long double, 12)
#ifdef LLONG_MAX
VIGRA_TYPE_TRAITS(long long, 13)
VIGRA_TYPE_TRAITS(unsigned long long, 14)
#endif

#undef VIGRA_TYPE_TRAITS

//@}

template <class A>
struct Not;

template <>
struct Not<VigraTrueType>
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct Not<VigraFalseType>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class L, class R>
struct And;

template <>
struct And<VigraFalseType, VigraFalseType>
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct And<VigraFalseType, VigraTrueType>
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct And<VigraTrueType, VigraFalseType>
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct And<VigraTrueType, VigraTrueType>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class L, class R>
struct Or;

template <>
struct Or<VigraFalseType, VigraFalseType>
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <>
struct Or<VigraTrueType, VigraFalseType>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

template <>
struct Or<VigraFalseType, VigraTrueType>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

template <>
struct Or<VigraTrueType, VigraTrueType>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class PREDICATE, class TRUECASE, class FALSECASE>
struct If;

template <class TRUECASE, class FALSECASE>
struct If<VigraTrueType, TRUECASE, FALSECASE>
{
    typedef TRUECASE type;
};

template <class TRUECASE, class FALSECASE>
struct If<VigraFalseType, TRUECASE, FALSECASE>
{
    typedef FALSECASE type;
};

template <bool PREDICATE, class TRUECASE, class FALSECASE>
struct IfBool;

template <class TRUECASE, class FALSECASE>
struct IfBool<true, TRUECASE, FALSECASE>
{
    typedef TRUECASE type;
};

template <class TRUECASE, class FALSECASE>
struct IfBool<false, TRUECASE, FALSECASE>
{
    typedef FALSECASE type;
};

template <class L, class R>
struct IsSameType
{
    typedef VigraFalseType result;        // deprecated
    static const bool boolResult = false; // deprecated
    typedef VigraFalseType type;
    static const bool value = false;
};

template <class T>
struct IsSameType<T, T>
{
    typedef VigraTrueType result;        // deprecated
    static const bool boolResult = true; // deprecated
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class L, class R>
struct IsDifferentType
{
    typedef VigraTrueType type;
    static const bool value = true;
};

template <class T>
struct IsDifferentType<T, T>
{
    typedef VigraFalseType type;
    static const bool value = false;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class DERIVED, class BASE>
struct IsDerivedFrom
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * testIsDerivedFrom(...);
    static trueResult * testIsDerivedFrom(BASE const *);
    
    enum { resultSize = sizeof(*testIsDerivedFrom((DERIVED const *)0)) };
    
    static const bool value = (resultSize == 2);
    typedef typename 
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;

    static const bool boolResult = value; // deprecated
    typedef type result;                  // deprecated
};

template <class T>
struct UnqualifiedType
{
    typedef T type;
};

template <class T>
struct UnqualifiedType<T &>
: public UnqualifiedType<T>
{};

template <class T>
struct UnqualifiedType<T const &>
: public UnqualifiedType<T>
{};

template <class T>
struct UnqualifiedType<T *>
: public UnqualifiedType<T>
{};

template <class T>
struct UnqualifiedType<T const *>
: public UnqualifiedType<T>
{};


template <class T>
struct has_result_type
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * test(...);
    template <class U>
    static trueResult * test(U *, typename U::result_type * = 0);
    
    enum { resultSize = sizeof(*test((T*)0)) };
    
    static const bool value = (resultSize == 2);
    typedef typename 
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;
};

template <class T>
struct has_value_type
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * test(...);
    template <class U>
    static trueResult * test(U *, typename U::value_type * = 0);
    
    enum { resultSize = sizeof(*test((T*)0)) };
    
    static const bool value = (resultSize == 2);
    typedef typename 
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;
};

} // namespace vigra

#endif /* VIGRA_METAPROGRAMMING_HXX */
