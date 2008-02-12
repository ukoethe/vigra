/************************************************************************/
/*                                                                      */
/*         Copyright 2002-2003 by Ullrich Koethe, Hans Meine            */
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

#ifndef VIGRA_MEMORY_HXX
#define VIGRA_MEMORY_HXX

#include "metaprogramming.hxx"

namespace vigra { namespace detail {

template <class T>
inline void destroy_n(T * /* p */, int /* n */, VigraTrueType /* isPOD */)
{
}

template <class T>
inline void destroy_n(T * p, int n, VigraFalseType /* isPOD */)
{
	T * end = p + n;
	for(; p != end; ++p)
		p->~T();
}

template <class T>
inline void destroy_n(T * p, int n)
{
    destroy_n(p, n, typename TypeTraits<T>::isPOD());
}

/********************************************************************/

// g++ 2.95 has std::destroy() in the STL
#if !defined(__GNUC__) ||  __GNUC__ >= 3

template <class T>
inline void destroy(T * p, VigraTrueType /* isPOD */)
{
}

template <class T>
inline void destroy(T * p, VigraFalseType /* isPOD */)
{
    p->~T();
}

template <class T>
inline void destroy(T * p)
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
