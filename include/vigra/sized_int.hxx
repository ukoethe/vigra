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


#ifndef VIGRA_SIZED_INT_HXX
#define VIGRA_SIZED_INT_HXX

#include "metaprogramming.hxx"

namespace vigra {

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

namespace detail {

template<class T, class NEXT>
struct IntTypeList
{
    enum { size = sizeof(T)*8 };
    typedef T type;
    typedef NEXT next;
};

template<int SIZE, class LIST>
struct SelectIntegerType
{
    typedef typename 
       IfBool<(SIZE == LIST::size), 
           typename LIST::type,
           typename SelectIntegerType<SIZE, typename LIST::next>::type >::type
       type;
};

template<int SIZE>
struct SelectIntegerType<SIZE, VigraFalseType>
{
    typedef VigraFalseType type;
};

typedef IntTypeList<signed char, 
        IntTypeList<signed short,
        IntTypeList<signed int,
        IntTypeList<signed long,
        VigraFalseType > > > > SignedIntTypes;
typedef IntTypeList<unsigned char, 
        IntTypeList<unsigned short,
        IntTypeList<unsigned int,
        IntTypeList<unsigned long,
        VigraFalseType > > > > UnsignedIntTypes;

} // namespace detail

typedef detail::SelectIntegerType<8,  detail::SignedIntTypes>::type Int8;
typedef detail::SelectIntegerType<16, detail::SignedIntTypes>::type Int16;
typedef detail::SelectIntegerType<32, detail::SignedIntTypes>::type Int32;
typedef detail::SelectIntegerType<64, detail::SignedIntTypes>::type Int64;
typedef detail::SelectIntegerType<8,  detail::UnsignedIntTypes>::type UInt8;
typedef detail::SelectIntegerType<16, detail::UnsignedIntTypes>::type UInt16;
typedef detail::SelectIntegerType<32, detail::UnsignedIntTypes>::type UInt32;
typedef detail::SelectIntegerType<64, detail::UnsignedIntTypes>::type UInt64;

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

typedef signed char    Int8;
typedef signed short   Int16;
typedef signed int     Int32;
typedef VigraFalseType Int64;
typedef unsigned char  UInt8;
typedef unsigned short UInt16;
typedef unsigned int   UInt32;
typedef VigraFalseType UInt64;

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

} // namespace vigra

#endif /* VIGRA_SIZED_INT_HXX */
