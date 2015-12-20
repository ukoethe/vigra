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

#ifdef __APPLE__
#include <AssertMacros.h>
#undef check
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
/*          tags for MultiArray memory layout           */
/*                                                      */
/********************************************************/

/** tag for marking a MultiArray strided.

    <b>\#include</b> \<vigra/metaprogramming.hxx\> <br/>
    Namespace: vigra
*/
struct StridedArrayTag {};

/** tag for marking a MultiArray unstrided.

    <b>\#include</b> \<vigra/metaprogramming.hxx\> <br/>
    Namespace: vigra
*/
struct UnstridedArrayTag {};

/** tag for marking a MultiArray chunked.

    <b>\#include</b> \<vigra/metaprogramming.hxx\> <br/>
    Namespace: vigra
*/
struct ChunkedArrayTag {};

/********************************************************/
/*                                                      */
/*                      TypeTraits                      */
/*                                                      */
/********************************************************/

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
    
    enum { resultSize = sizeof(*testIsDerivedFrom(static_cast<DERIVED const *>(0))) };
    
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
    static const bool isConst = false;
    static const bool isReference = false;
    static const bool isPointer = false;
};

template <class T>
struct UnqualifiedType<T const>
{
    typedef T type;
    static const bool isConst = true;
    static const bool isReference = false;
    static const bool isPointer = false;
};

template <class T>
struct UnqualifiedType<T &>
: public UnqualifiedType<T>
{
    static const bool isReference = true;
};

template <class T>
struct UnqualifiedType<T const &>
: public UnqualifiedType<T const>
{
    static const bool isReference = true;
};

template <class T>
struct UnqualifiedType<T *>
: public UnqualifiedType<T>
{
    static const bool isPointer = true;
};

template <class T>
struct UnqualifiedType<T const *>
: public UnqualifiedType<T const>
{
    static const bool isPointer = true;
};

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
    
    enum { resultSize = sizeof(*test(static_cast<T*>(0))) };
    
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
struct has_real_promote_type : public sfinae_test<T, has_real_promote_type>
{
    template <class U>
    has_real_promote_type(U*, typename U::real_promote_type* = 0);
};

template <class T, bool P = has_real_promote_type<T>::value>
struct get_optional_real_promote
{
    typedef T type;
};
template <class T>
struct get_optional_real_promote<T, true>
{
    typedef typename T::real_promote_type type;
};

template <class T>
struct IsArray
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * test(...);
    template <class U, unsigned n>
    static trueResult * test(U (*)[n]);
    
    enum { resultSize = sizeof(*test(static_cast<T*>(0))) };
    
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

template <int X, unsigned int N>
struct MetaPow
{
    static const long long value = MetaPow<X, N-1>::value * X;
};

template <int X>
struct MetaPow<X, 0>
{
    static const long long value = 1;
};

/****************************************************************************/
/*                                                                          */
/*                        TypeList and its functions                        */
/*                                                                          */
/****************************************************************************/

template<class HEAD, class TAIL=void>
struct TypeList
{
    typedef TypeList<HEAD, TAIL> type;
    typedef HEAD Head;
    typedef TAIL Tail;
};

template <class List, class T>
struct Contains;

template <class Head, class Tail, class T>
struct Contains<TypeList<Head, Tail>, T>
{
    typedef typename Contains<Tail, T>::type type;
};

template <class Head, class Tail>
struct Contains<TypeList<Head, Tail>, Head>
{
    typedef VigraTrueType type;
};

template <class T>
struct Contains<void, T>
{
    typedef VigraFalseType type;
};

template <class List, class T>
struct Remove;

template <class Head, class Tail, class T>
struct Remove<TypeList<Head, Tail>, T>
{
    typedef TypeList<Head, typename Remove<Tail, T>::type> type;
};

template <class Head, class Tail>
struct Remove<TypeList<Head, Tail>, Head>
{
    typedef Tail type;
};

template <class T>
struct Remove<void, T>
{
    typedef void type;
};

template <class A, class Tail=void>
struct Push
{
    typedef TypeList<A, typename Tail::type> type;
};

template <class Head, class Tail, class List>
struct Push<TypeList<Head, Tail>, List>
{
    typedef typename Push<Tail, List>::type Rest;
    typedef TypeList<Head, Rest> type;
};

template <class Head, class Tail>
struct Push<TypeList<Head, Tail>, void>
{
    typedef TypeList<Head, Tail> type;
};

template <class A>
struct Push<A, void>
{
    typedef TypeList<A> type;
};

template <class A>
struct Push<void, A>
{
    typedef A type;
};

template <>
struct Push<void, void>
{
    typedef void type;
};

template <class A, class Tail=void>
struct PushUnique
{
    typedef typename Contains<Tail, A>::type AlreadyInList;
    typedef typename If<AlreadyInList, typename Tail::type, TypeList<A, typename Tail::type> >::type type;
};

template <class Head, class Tail, class List>
struct PushUnique<TypeList<Head, Tail>, List>
{
    typedef typename PushUnique<Tail, List>::type Rest;
    typedef typename Contains<Rest, Head>::type HeadAlreadyInList;
    typedef typename If<HeadAlreadyInList, Rest, TypeList<Head, Rest> >::type type;
};

template <class Head, class Tail>
struct PushUnique<TypeList<Head, Tail>, void>
{
    typedef TypeList<Head, Tail> type;
};

template <class A>
struct PushUnique<A, void>
{
    typedef TypeList<A> type;
};

template <class A>
struct PushUnique<void, A>
{
    typedef A type;
};

template <>
struct PushUnique<void, void>
{
    typedef void type;
};

template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct MakeTypeList
{
    typedef typename Push<T19, T20>::type L19;
    typedef typename Push<T18, L19>::type L18;
    typedef typename Push<T17, L18>::type L17;
    typedef typename Push<T16, L17>::type L16;
    typedef typename Push<T15, L16>::type L15;
    typedef typename Push<T14, L15>::type L14;
    typedef typename Push<T13, L14>::type L13;
    typedef typename Push<T12, L13>::type L12;
    typedef typename Push<T11, L12>::type L11;
    typedef typename Push<T10, L11>::type L10;
    typedef typename Push<T09, L10>::type L09;
    typedef typename Push<T08, L09>::type L08;
    typedef typename Push<T07, L08>::type L07;
    typedef typename Push<T06, L07>::type L06;
    typedef typename Push<T05, L06>::type L05;
    typedef typename Push<T04, L05>::type L04;
    typedef typename Push<T03, L04>::type L03;
    typedef typename Push<T02, L03>::type L02;
    typedef typename Push<T01, L02>::type L01;
    typedef L01 type;
};

template <class T01=void, class T02=void, class T03=void, class T04=void, class T05=void,
          class T06=void, class T07=void, class T08=void, class T09=void, class T10=void,
          class T11=void, class T12=void, class T13=void, class T14=void, class T15=void,
          class T16=void, class T17=void, class T18=void, class T19=void, class T20=void>
struct MakeTypeListUnique
{
    typedef typename PushUnique<T19, T20>::type L19;
    typedef typename PushUnique<T18, L19>::type L18;
    typedef typename PushUnique<T17, L18>::type L17;
    typedef typename PushUnique<T16, L17>::type L16;
    typedef typename PushUnique<T15, L16>::type L15;
    typedef typename PushUnique<T14, L15>::type L14;
    typedef typename PushUnique<T13, L14>::type L13;
    typedef typename PushUnique<T12, L13>::type L12;
    typedef typename PushUnique<T11, L12>::type L11;
    typedef typename PushUnique<T10, L11>::type L10;
    typedef typename PushUnique<T09, L10>::type L09;
    typedef typename PushUnique<T08, L09>::type L08;
    typedef typename PushUnique<T07, L08>::type L07;
    typedef typename PushUnique<T06, L07>::type L06;
    typedef typename PushUnique<T05, L06>::type L05;
    typedef typename PushUnique<T04, L05>::type L04;
    typedef typename PushUnique<T03, L04>::type L03;
    typedef typename PushUnique<T02, L03>::type L02;
    typedef typename PushUnique<T01, L02>::type L01;
    typedef L01 type;
};

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra

#endif /* VIGRA_METAPROGRAMMING_HXX */
