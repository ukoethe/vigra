/************************************************************************/
/*                                                                      */
/*         Copyright 2002-2003 by Ullrich Koethe, Hans Meine            */
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

#ifndef VIGRA_MEMORY_HXX
#define VIGRA_MEMORY_HXX


#include "config.hxx"

#ifdef VIGRA_SHARED_PTR_IN_TR1
#  include <tr1/memory>
#else
#  include <memory>
#endif 

#include <cstring>
#include "metaprogramming.hxx"

namespace vigra { 

enum SkipInitializationTag { SkipInitialization};

template<class T>
struct CanSkipInitialization
{
    typedef typename TypeTraits<T>::isBuiltinType type;
    static const bool value = type::asBool;
};

namespace detail {

// differs from std::uninitialized_copy by explicit type conversion
template <class Src, class Dest>
Dest uninitializedCopy(Src s, Src end, Dest d)
{
    typedef typename std::iterator_traits<Dest>::value_type T;
    for(; s != end; ++s, ++d)
        new(d) T(static_cast<T const &>(*s));
    return d;
}

template <class T>
struct PlacementNewAllocator
{
    T * allocate(std::size_t n)
    {
        return (T*)::operator new(n*sizeof(T));
    }

    void deallocate(T * p, std::size_t)
    {
        return ::operator delete(p);
    }

    void construct(T * p, T const & initial)
    {
        new(p) T(initial);
    }
    
    void destroy(T * p)
    {
        p->~T();
    }
};

template <class T>
inline void 
destroy_n(T * /* p */, std::size_t /* n */, VigraTrueType /* isPOD */)
{}

template <class T>
inline void 
destroy_n(T * p, std::size_t n, VigraFalseType /* isPOD */)
{
    for(std::size_t i=0; i < n; ++i, ++p)
        p->~T();
}

template <class T>
inline void
destroy_n(T * p, std::size_t n)
{
    destroy_n(p, n, typename TypeTraits<T>::isPOD());
}

template <class T, class Alloc>
inline T * 
alloc_initialize_n(std::size_t n, T const & initial, Alloc & alloc)
{
    T * p = alloc.allocate(n);
    bool useMemset = TypeTraits<T>::isPOD::value &&
                     (initial == T());
    if(useMemset)
    {
        std::memset(p, 0, n*sizeof(T));
    }
    else
    {
        std::size_t i=0;
        try 
        {
            for (; i < n; ++i)
                alloc.construct(p+i, initial);
        }
        catch (...) 
        {
            for (std::size_t j=0; j < i; ++j)
                alloc.destroy(p+j);
            alloc.deallocate(p, n);
            throw;
        }
    }
    return p;
}

template <class T>
inline T * 
alloc_initialize_n(std::size_t n, T const & initial)
{
    PlacementNewAllocator<T> alloc;
    return alloc_initialize_n<T>(n, initial, alloc);
}

template <class T>
inline T * 
alloc_initialize_n(std::size_t n)
{
    PlacementNewAllocator<T> alloc;
    return alloc_initialize_n<T>(n, T(), alloc);
}

template <class T, class Alloc>
inline void 
destroy_dealloc_impl(T * p, std::size_t n, Alloc & alloc, VigraTrueType /* isPOD */)
{
    alloc.deallocate(p, n);
}

template <class T, class Alloc>
inline void 
destroy_dealloc_impl(T * p, std::size_t n, Alloc & alloc, VigraFalseType /* isPOD */)
{
    for (std::size_t i=0; i < n; ++i)
        alloc.destroy(p + i);
    alloc.deallocate(p, n);
}

template <class T, class Alloc>
inline T * 
destroy_dealloc_n(T * p, std::size_t n, Alloc & alloc)
{
    if(p != 0)
        destroy_dealloc_impl(p, n, alloc, typename TypeTraits<T>::isPOD());
    return 0;
}

template <class T>
inline T * 
destroy_dealloc_n(T * p, std::size_t n)
{
    if(p != 0)
    {
        PlacementNewAllocator<T> alloc;
        destroy_dealloc_impl(p, n, alloc, typename TypeTraits<T>::isPOD());
    }
    return 0;
}

/********************************************************************/

// g++ 2.95 has std::destroy() in the STL
#if !defined(__GNUC__) ||  __GNUC__ >= 3

template <class T>
inline void destroy(T * /* p */, VigraTrueType /* isPOD */)
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
