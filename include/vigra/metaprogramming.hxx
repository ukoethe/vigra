/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_METAPROGRAMMING_HXX
#define VIGRA_METAPROGRAMMING_HXX

namespace vigra {

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
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

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
"<a href="multi_array_8hxx-source.html">vigra/multi_array.hxx</a>"

Namespace: vigra
*/
struct UnstridedArrayTag {};

template<class T>
class TypeTraits
{
  public:
    typedef VigraFalseType isPOD;
    typedef VigraFalseType isBuiltinType;
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template<class T> 
class TypeTraits<T *>
{
  public:
    typedef VigraTrueType isPOD;
    typedef VigraTrueType isBuiltinType;
};

template<class T> 
class TypeTraits<T const *>
{
  public:
    typedef VigraTrueType isPOD;
    typedef VigraTrueType isBuiltinType;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

#define VIGRA_TYPE_TRAITS(type) \
template<> \
class TypeTraits<type> \
{ \
  public: \
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

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class PREDICATE, class TRUECASE, class FALSECASE>
struct If
{
    typedef TRUECASE type;
};

template <class TRUECASE, class FALSECASE>
struct If<VigraFalseType, TRUECASE, FALSECASE>
{
    typedef FALSECASE type;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

} // namespace vigra

#endif /* VIGRA_METAPROGRAMMING_HXX */
