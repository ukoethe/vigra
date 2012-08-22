/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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

#ifndef VIGRA_METAPROGRAMMING_HXX
#define VIGRA_METAPROGRAMMING_HXX

#include "config.hxx"
#include <climits>
#include <limits>
#include <algorithm>

namespace vigra {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

template <int N>
class MetaInt
{
  public:
    static const int value = N;
};

template <int N1, int N2>
class MetaMax
{
  public:
    static const int value = N1 < N2 ? N2 : N1;
};

template <int N1, int N2>
class MetaMin
{
  public:
    static const int value = N1 < N2 ? N1 : N2;
};

struct VigraTrueType
{
   static const bool asBool = true, value = true;
};

struct VigraFalseType
{
    static const bool asBool = false, value = false;
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
\<vigra/multi_array.hxx\>

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
\<vigra/multi_array.hxx\>

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

template <class From, class To>
struct IsConvertibleTo
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static From const & check();
    
    static falseResult * testIsConvertible(...);
    static trueResult * testIsConvertible(To const &);
    
    enum { resultSize = sizeof(*testIsConvertible(check())) };
    
    static const bool value = (resultSize == 2);
    typedef typename 
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;
};

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
struct UnqualifiedType<T const>
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

template <bool, class T = void>
struct enable_if {};
template <class T>
struct enable_if<true, T> { typedef T type; };

struct sfinae_void;

template <class T, template<class> class USER>
struct sfinae_test
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * test(...);
    static trueResult * test(USER<sfinae_void>);
    
    enum { resultSize = sizeof(*test((T*)0)) };
    
    static const bool value = (resultSize == 2);
    typedef typename
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;
};

template <class T>
struct has_argument_type : public sfinae_test<T, has_argument_type>
{
    template <class U> has_argument_type(U*, typename U::argument_type* = 0);
};

template <class T>
struct has_result_type : public sfinae_test<T, has_result_type>
{
    template <class U> has_result_type(U*, typename U::result_type* = 0);
};

template <class T>
struct has_value_type : public sfinae_test<T, has_value_type>
{
    template <class U> has_value_type(U*, typename U::value_type* = 0);
};

template <class T>
struct IsIterator : public sfinae_test<T, IsIterator>
{
    template <class U> IsIterator(U*, typename U::iterator_category* = 0);
};

template <class T>
struct IsIterator<T*>
{
    static const bool value = true;
    typedef VigraTrueType type;
};

template <class T>
struct IsIterator<T const *>
{
    static const bool value = true;
    typedef VigraTrueType type;
};

template <class T>
struct IsArray
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * test(...);
    template <class U, unsigned n>
    static trueResult * test(U (*)[n]);
    
    enum { resultSize = sizeof(*test((T*)0)) };
    
    static const bool value = (resultSize == 2);
    typedef typename
        IfBool<value, VigraTrueType, VigraFalseType>::type
        type;
};


template <class D, class B, class Z> inline
D & static_cast_2(Z & z)
{
    return static_cast<D &>(static_cast<B &>(z));
}

template <class A>
class copy_if_same_as
{
    const bool copied;
    const A *const data;
    copy_if_same_as(const copy_if_same_as &);
    void operator=(const copy_if_same_as &);
public:
    copy_if_same_as(const A & x, const A & y)
        : copied(&x == &y), data(copied ? new A(y) : &x) {}
    ~copy_if_same_as()
    {
        if (copied)
            delete data;
    }
    const A & operator()() const { return *data; }
};


template <class>
struct true_test : public VigraTrueType {};

template <class>
struct false_test : VigraFalseType {};

template <class PC, class T, class F>
struct ChooseBool
{
    static const bool value = IfBool<PC::value, T, F>::type::value;
};

template <bool>
struct choose_type
{
    template <class A, class B>
    static const A & at(const A & a, const B &) { return a; }
    template <class A, class B>
    static       A & at(      A & a,       B &) { return a; }
};
template <>
struct choose_type<false>
{
    template <class A, class B>
    static const B & at(const A &, const B & b) { return b; }
    template <class A, class B>
    static       B & at(      A &,       B & b) { return b; }
};

template <class X>
struct HasMetaLog2
{
    static const bool value =   !std::numeric_limits<X>::is_signed
                              && std::numeric_limits<X>::is_integer;
};
template <class X>
struct EnableMetaLog2
    : public enable_if<HasMetaLog2<X>::value> {};
template <class>
class vigra_error_MetaLog2_accepts_only_unsigned_types_and_no_;

// use a conforming template depth here (below 15 for up to 128 bits)
template <class X = unsigned long,
          X n = ~(X(0)), unsigned s = 1, unsigned t = 0, bool q = 1,
          X m = 0, X z = 0, X u = 1, class = void>
class MetaLog2
    : public vigra_error_MetaLog2_accepts_only_unsigned_types_and_no_<X>
{};
template <class X, X n, unsigned s, unsigned t, bool q, X m, X z, X u>
struct MetaLog2 <X, n, s, t, q, m, z, u, typename EnableMetaLog2<X>::type>
{
    static const unsigned value
        = t + MetaLog2<X, (n >> s), s * (1 + q), s, !q, n / 2, z, u>::value;
};
template <class X, unsigned s, unsigned t, bool q, X m, X z, X u>
struct MetaLog2<X, z, s, t, q, m, z, u, typename EnableMetaLog2<X>::type>
{
    static const unsigned value
        = 1 + MetaLog2<X, m / 2, 2, 1, 1, 0, z, u>::value;
};
template <class X, unsigned s, unsigned t, bool q, X z, X u>
struct MetaLog2<X, z, s, t, q, u, z, u, typename EnableMetaLog2<X>::type>
{
    static const unsigned value = 2;
};
template <class X, unsigned s, unsigned t, bool q, X z, X u>
struct MetaLog2<X, z, s, t, q, z, z, u, typename EnableMetaLog2<X>::type>
{
    static const unsigned value = 1;
};
template <class X, X z, X u>
struct MetaLog2<X, z, 1, 0, 1, z, z, u, typename EnableMetaLog2<X>::type>
{
    // A value of zero for MetaLog2<X, 0> is likely to cause most harm,
    // such as division by zero or zero array sizes, this is actually indended.
    static const unsigned value = 0;
};

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra

#endif /* VIGRA_METAPROGRAMMING_HXX */
