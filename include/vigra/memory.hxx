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
void destroy_n(T * /* p */, int /* n */, VigraTrueType /* isPOD */)
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
    destroy_n(p, n, typename TypeTraits<T>::isPOD());
}

/********************************************************************/

// g++ 2.95 has std::destroy() in the STL
#if !defined(__GNUC__) ||  __GNUC__ >= 3

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
    destroy(p, typename TypeTraits<T>::isPOD());
}

#else

} } // namespace vigra::detail

#include <memory>

namespace vigra { namespace detail {

using std::destroy;

#endif

} } // namespace vigra::detail

#endif // VIGRA_MEMORY_HXX
