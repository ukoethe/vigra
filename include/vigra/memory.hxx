/************************************************************************/
/*                                                                      */
/*         Copyright 2002-2003 by Ullrich Koethe, Hans Meine            */
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

#ifndef VIGRA_MEMORY_HXX
#define VIGRA_MEMORY_HXX

#include <vigra/metaprogramming.hxx>

namespace vigra { namespace detail {

template <class T>
void destroy_n(T * p, int n, VigraTrueType /* isPOD */)
{
}

template <class T>
void destroy_n(T * p, int n, VigraFalseType /* isPOD */)
{
	T * end = p + n;
	for(; p != end; ++p)
		p->~T();
}

template <class T>
void destroy_n(T * p, int n)
{
    destroy_n(p, n, TypeTraits<T>::isPOD());
}

/********************************************************************/

template <class T>
void destroy(T * p, VigraTrueType /* isPOD */)
{
}

template <class T>
void destroy(T * p, VigraFalseType /* isPOD */)
{
    p->~T();
}

template <class T>
void destroy(T * p)
{
    destroy(p, TypeTraits<T>::isPOD());
}

} } // namespace vigra::detail

#endif // VIGRA_MEMORY_HXX
