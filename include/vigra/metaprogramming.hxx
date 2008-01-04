/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

namespace vigra {

template <int N>
class MetaInt
{
  public:
    enum { value = N };
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

#define VIGRA_TYPE_TRAITS(type) \
template<> \
class TypeTraits<type> \
{ \
  public: \
    typedef VigraFalseType isConst; \
    typedef VigraTrueType isPOD; \
    typedef VigraTrueType isBuiltinType; \
};

VIGRA_TYPE_TRAITS(char)
VIGRA_TYPE_TRAITS(signed char)
VIGRA_TYPE_TRAITS(unsigned char)
VIGRA_TYPE_TRAITS(short)
VIGRA_TYPE_TRAITS(unsigned short)
VIGRA_TYPE_TRAITS(int)
VIGRA_TYPE_TRAITS(unsigned int)
VIGRA_TYPE_TRAITS(long)
VIGRA_TYPE_TRAITS(unsigned long)
VIGRA_TYPE_TRAITS(float)
VIGRA_TYPE_TRAITS(double)
VIGRA_TYPE_TRAITS(long double)

#undef VIGRA_TYPE_TRAITS

//@}

template <class L, class R>
struct And;

template <>
struct And<VigraFalseType, VigraFalseType>
{
    typedef VigraFalseType result;
    static const bool boolResult = false;
};

template <>
struct And<VigraFalseType, VigraTrueType>
{
    typedef VigraFalseType result;
    static const bool boolResult = false;
};

template <>
struct And<VigraTrueType, VigraFalseType>
{
    typedef VigraFalseType result;
    static const bool boolResult = false;
};

template <>
struct And<VigraTrueType, VigraTrueType>
{
    typedef VigraTrueType result;
    static const bool boolResult = true;
};

template <class L, class R>
struct Or;

template <>
struct Or<VigraFalseType, VigraFalseType>
{
    typedef VigraFalseType result;
    static const bool boolResult = false;
};

template <>
struct Or<VigraTrueType, VigraFalseType>
{
    typedef VigraTrueType result;
    static const bool boolResult = true;
};

template <>
struct Or<VigraFalseType, VigraTrueType>
{
    typedef VigraTrueType result;
    static const bool boolResult = true;
};

template <>
struct Or<VigraTrueType, VigraTrueType>
{
    typedef VigraTrueType result;
    static const bool boolResult = true;
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
    typedef VigraFalseType result;
    static const bool boolResult = false;
};

template <class T>
struct IsSameType<T, T>
{
    typedef VigraTrueType result;
    static const bool boolResult = true;
};

template <class DERIVED, class BASE>
struct IsDerivedFrom
{
    typedef char falseResult[1];
    typedef char trueResult[2];
    
    static falseResult * testIsDerivedFrom(...);
    static trueResult * testIsDerivedFrom(BASE const *);
    
    enum { resultSize = sizeof(*testIsDerivedFrom((DERIVED const *)0)) };
    
    static const bool boolResult = (resultSize == 2);
    typedef typename 
        IfBool<boolResult, VigraTrueType, VigraFalseType>::type
        result;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

} // namespace vigra

#endif /* VIGRA_METAPROGRAMMING_HXX */
