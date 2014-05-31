/*
 * Copyright (C) 2009-2014 Christoph L. Spiel
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef OPENMP_DEF_H_INCLUDED_
#define OPENMP_DEF_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <limits>

#include <signal.h>             // sig_atomic_t

#ifdef HAVE_TCMALLOC_H
#include <tcmalloc.h>
#endif


#if _OPENMP >= 200203 // at least OpenMP version 2.0

#include <omp.h>

#include "muopt.h"


#define OPENMP
#define OPENMP_YEAR (_OPENMP / 100)
#define OPENMP_MONTH (_OPENMP % 100)


namespace omp
{
    inline static void* malloc(size_t) __attribute__((alloc_size(1), malloc));

    inline static void*
    malloc(size_t size)
    {
#if defined(__ICC)
#define OMP_MALLOC_FUNCTIONS "kmp_malloc/kmp_delete"
        return kmp_malloc(size);
#elif defined(HAVE_TCMALLOC)
#define OMP_MALLOC_FUNCTIONS "tc_malloc/tc_delete"
        return tc_malloc(size);
#else
#define OMP_MALLOC_FUNCTIONS "malloc/free"
        return ::malloc(size);
#endif
    }

    inline static void
    free(void* pointer)
    {
#if defined(__ICC)
        kmp_free(pointer);
#elif defined(HAVE_TCMALLOC)
        tc_free(pointer);
#else
        ::free(pointer);
#endif
    }
} // namespace omp

#else

#undef OPENMP
#define OPENMP_YEAR 0
#define OPENMP_MONTH 0

inline void omp_set_num_threads(int) {}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_max_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
inline int omp_get_num_procs() {return 1;}
inline void omp_set_dynamic(int) {}
inline int omp_get_dynamic() {return 0;}
inline int omp_in_parallel() {return 0;}
inline void omp_set_nested(int) {}
inline int omp_get_nested() {return 0;}


namespace omp
{
    inline static void* malloc(size_t) __attribute__((alloc_size(1), malloc));

    inline static void*
    malloc(size_t size)
    {
        return std::malloc(size);
    }

    inline static void
    free(void* pointer)
    {
        std::free(pointer);
    }
} // namespace omp

#endif // _OPENMP >= 200505


// Answer whether the underlying OpenMP implementation really (thinks
// that it) supports nested parallelism.
inline static bool
have_openmp_nested()
{
    const bool openmp_nested = omp_get_nested();
    omp_set_nested(true);
    const bool result = omp_get_nested() != 0;
    omp_set_nested(openmp_nested);
    return result;
}


// Answer whether the underlying OpenMP implementation really (thinks
// that it) supports dynamic adjustment of the number of threads.
inline static bool
have_openmp_dynamic()
{
    const bool openmp_dynamic = omp_get_dynamic();
    omp_set_dynamic(true);
    const bool result = omp_get_dynamic() != 0;
    omp_set_dynamic(openmp_dynamic);
    return result;
}


namespace omp
{
    typedef sig_atomic_t atomic_t;


    class scoped_nested
    {
    public:
        scoped_nested() : level_(omp_get_nested()) {}

        explicit scoped_nested(bool allow_nested) : level_(omp_get_nested())
        {
            omp_set_nested(static_cast<int>(allow_nested));
        }

        ~scoped_nested() {omp_set_nested(level_);}

    private:
        const int level_;
    };
} // namespace omp


#endif // OPENMP_DEF_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
