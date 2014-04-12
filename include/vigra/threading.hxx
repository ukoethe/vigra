/************************************************************************/
/*                                                                      */
/*               Copyright 2013-2014 by Ullrich Koethe                  */
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

#ifndef VIGRA_THREADING_HXX
#define VIGRA_THREADING_HXX

/* Compatibility header to import threading-related functionality from boost
   when the compiler doesn't yet support C++11.
*/


#if defined(__clang__)
# if (!__has_include(<thread>) || !__has_include(<mutex>) || !__has_include(<atomic>))
#    define VIGRA_NO_STD_THREADING
# endif
#else
#  if defined(__GNUC__) && (!defined(_GLIBCXX_HAS_GTHREADS) || !defined(_GLIBCXX_USE_C99_STDINT_TR1))
#    define VIGRA_NO_STD_THREADING
#  endif
#endif

#if defined(_MSC_VER) && _MSC_VER <= 1600
#  define VIGRA_NO_STD_THREADING
#endif

#ifdef VIGRA_NO_STD_THREADING
#  include <boost/thread/thread.hpp>
#  include <boost/thread/mutex.hpp>
// #  include <boost/thread/shared_mutex.hpp>
#  include <boost/thread/locks.hpp>
#  if BOOST_VERSION >= 105300
#    include <boost/atomic.hpp>
#    define VIGRA_HAS_ATOMIC 1
#  endif
#  define VIGRA_THREADING_NAMESPACE boost
#else
#  include <thread>
#  include <mutex>
// #  include <shared_mutex>  // C++14
#  include <atomic>
#  define VIGRA_HAS_ATOMIC 1
#  define VIGRA_THREADING_NAMESPACE std
#endif

#if defined(_MSC_VER) && !defined(VIGRA_HAS_ATOMIC)
#  include "windows.h"
#endif

namespace vigra { namespace threading {

// contents of <thread>

using VIGRA_THREADING_NAMESPACE::thread;

namespace this_thread {

using VIGRA_THREADING_NAMESPACE::this_thread::yield;
using VIGRA_THREADING_NAMESPACE::this_thread::get_id;
using VIGRA_THREADING_NAMESPACE::this_thread::sleep_for;
using VIGRA_THREADING_NAMESPACE::this_thread::sleep_until;

} // namespace this_thread

// contents of <mutex>

using VIGRA_THREADING_NAMESPACE::mutex;
using VIGRA_THREADING_NAMESPACE::timed_mutex;
using VIGRA_THREADING_NAMESPACE::recursive_mutex;
using VIGRA_THREADING_NAMESPACE::recursive_timed_mutex;

using VIGRA_THREADING_NAMESPACE::lock_guard;
using VIGRA_THREADING_NAMESPACE::unique_lock;

using VIGRA_THREADING_NAMESPACE::defer_lock_t;
using VIGRA_THREADING_NAMESPACE::try_to_lock_t;
using VIGRA_THREADING_NAMESPACE::adopt_lock_t;

using VIGRA_THREADING_NAMESPACE::defer_lock;
using VIGRA_THREADING_NAMESPACE::try_to_lock;
using VIGRA_THREADING_NAMESPACE::adopt_lock;

using VIGRA_THREADING_NAMESPACE::try_lock;
using VIGRA_THREADING_NAMESPACE::lock;

using VIGRA_THREADING_NAMESPACE::once_flag;
using VIGRA_THREADING_NAMESPACE::call_once;

// contents of <shared_mutex>

// using VIGRA_THREADING_NAMESPACE::shared_mutex;   // C++14
// using VIGRA_THREADING_NAMESPACE::shared_lock;  // C++14

#ifdef VIGRA_HAS_ATOMIC

// contents of <atomic>

using VIGRA_THREADING_NAMESPACE::atomic_flag;
using VIGRA_THREADING_NAMESPACE::atomic;

using VIGRA_THREADING_NAMESPACE::atomic_char;
using VIGRA_THREADING_NAMESPACE::atomic_schar;
using VIGRA_THREADING_NAMESPACE::atomic_uchar;
using VIGRA_THREADING_NAMESPACE::atomic_short;
using VIGRA_THREADING_NAMESPACE::atomic_ushort;
using VIGRA_THREADING_NAMESPACE::atomic_int;
using VIGRA_THREADING_NAMESPACE::atomic_uint;
using VIGRA_THREADING_NAMESPACE::atomic_long;
using VIGRA_THREADING_NAMESPACE::atomic_ulong;
using VIGRA_THREADING_NAMESPACE::atomic_llong;
using VIGRA_THREADING_NAMESPACE::atomic_ullong;
using VIGRA_THREADING_NAMESPACE::atomic_char16_t;
using VIGRA_THREADING_NAMESPACE::atomic_char32_t;
using VIGRA_THREADING_NAMESPACE::atomic_wchar_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_least8_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_least8_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_least16_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_least16_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_least32_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_least32_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_least64_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_least64_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_fast8_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_fast8_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_fast16_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_fast16_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_fast32_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_fast32_t;
using VIGRA_THREADING_NAMESPACE::atomic_int_fast64_t;
using VIGRA_THREADING_NAMESPACE::atomic_uint_fast64_t;
using VIGRA_THREADING_NAMESPACE::atomic_intptr_t;
using VIGRA_THREADING_NAMESPACE::atomic_uintptr_t;
using VIGRA_THREADING_NAMESPACE::atomic_size_t;
using VIGRA_THREADING_NAMESPACE::atomic_ptrdiff_t;
using VIGRA_THREADING_NAMESPACE::atomic_intmax_t;
using VIGRA_THREADING_NAMESPACE::atomic_uintmax_t;

using VIGRA_THREADING_NAMESPACE::memory_order;
using VIGRA_THREADING_NAMESPACE::memory_order_relaxed;
using VIGRA_THREADING_NAMESPACE::memory_order_release;
using VIGRA_THREADING_NAMESPACE::memory_order_acquire;
using VIGRA_THREADING_NAMESPACE::memory_order_consume;
using VIGRA_THREADING_NAMESPACE::memory_order_acq_rel;
using VIGRA_THREADING_NAMESPACE::memory_order_seq_cst;

using VIGRA_THREADING_NAMESPACE::atomic_thread_fence;
using VIGRA_THREADING_NAMESPACE::atomic_signal_fence;

#else  // VIGRA_HAS_ATOMIC not defined

enum memory_order {
    memory_order_relaxed,
    memory_order_release,
    memory_order_acquire,
    memory_order_consume,
    memory_order_acq_rel,
    memory_order_seq_cst
};

#ifdef _MSC_VER

template <int SIZE=4>
struct atomic_long_impl
{
    typedef LONG value_type;
    
    static long load(value_type const & val)
    {
        long res = val;
        MemoryBarrier();
        return res;
    }
    
    static void store(value_type & dest, long val)
    {
        MemoryBarrier();
        dest = val;
    }
    
    static long add(value_type & dest, long val)
    {
        return InterlockedExchangeAdd(&dest, val);
    }
    
    static long sub(value_type & dest, long val)
    {
        return InterlockedExchangeAdd(&dest, -val);
    }
    
    static bool compare_exchange(value_type & dest, long & old_val, long new_val)
    {
        long check_val = old_val;
        old_val = InterlockedCompareExchange(&dest, new_val, old_val);
        return check_val == old_val;
    }
};

template <>
struct atomic_long_impl<8>
{
    typedef LONGLONG value_type;
    
    static long load(value_type const & val)
    {
        long res = val;
        MemoryBarrier();
        return res;
    }
    
    static void store(value_type & dest, long val)
    {
        MemoryBarrier();
        dest = val;
    }
    
    static long add(value_type & dest, long val)
    {
        return InterlockedExchangeAdd64(&dest, val);
    }
    
    static long sub(value_type & dest, long val)
    {
        return InterlockedExchangeAdd64(&dest, -val);
    }
    
    static bool compare_exchange(value_type & dest, long & old_val, long new_val)
    {
        long check_val = old_val;
        old_val = InterlockedCompareExchange64(&dest, new_val, old_val);
        return check_val == old_val;
    }
};

#else

template <int SIZE=4>
struct atomic_long_impl
{
    typedef long value_type;
    
    static long load(value_type const & val)
    {
        long res = val;
        __sync_synchronize();
        return res;
    }
    
    static void store(value_type & dest, long val)
    {
        __sync_synchronize();
        dest = val;
    }
    
    static long add(value_type & dest, long val)
    {
        return __sync_fetch_and_add(&dest, val);
    }
    
    static long sub(value_type & dest, long val)
    {
        return __sync_fetch_and_sub(&dest, val);
    }
    
    static bool compare_exchange(value_type & dest, long & old_val, long new_val)
    {
        long check_val = old_val;
        old_val = __sync_val_compare_and_swap(&dest, old_val, new_val);
        return check_val == old_val;
    }
};

#endif // _MSC_VER

struct atomic_long
{
    typedef atomic_long_impl<sizeof(long)>::value_type value_type;
    
    atomic_long(long v = 0)
    : value_(v)
    {}
    
    atomic_long & operator=(long val)
    {
        store(val);
        return *this;
    }
    
    long load(memory_order = memory_order_seq_cst) const
    {
        return atomic_long_impl<sizeof(long)>::load(value_);
    }
    
    void store(long v, memory_order = memory_order_seq_cst)
    {
        atomic_long_impl<sizeof(long)>::store(value_, v);
    }
    
    long fetch_add(long v, memory_order = memory_order_seq_cst)
    {
        return atomic_long_impl<sizeof(long)>::add(value_, v);
    }
    
    long fetch_sub(long v, memory_order = memory_order_seq_cst)
    {
        return atomic_long_impl<sizeof(long)>::sub(value_, v);
    }
    
    bool compare_exchange_strong(long & old_val, long new_val, memory_order = memory_order_seq_cst)
    {
        return atomic_long_impl<sizeof(long)>::compare_exchange(value_, old_val, new_val);
    }
    
    bool compare_exchange_weak(long & old_val, long new_val, memory_order = memory_order_seq_cst)
    {
        return atomic_long_impl<sizeof(long)>::compare_exchange(value_, old_val, new_val);
    }

    value_type value_;
};

#endif // VIGRA_HAS_ATOMIC

// not in boost (?)
// using VIGRA_THREADING_NAMESPACE::atomic_is_lock_free;
// using VIGRA_THREADING_NAMESPACE::atomic_storeatomic_store_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_loadatomic_load_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_exchangeatomic_exchange_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_compare_exchange_weak;
// using VIGRA_THREADING_NAMESPACE::atomic_compare_exchange_weak_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_compare_exchange_strong;
// using VIGRA_THREADING_NAMESPACE::atomic_compare_exchange_strong_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_fetch_addatomic_fetch_add_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_fetch_subatomic_fetch_sub_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_fetch_andatomic_fetch_and_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_fetch_oratomic_fetch_or_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_fetch_xoratomic_fetch_xor_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_flag_test_and_setatomic_flag_test_and_set_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_flag_clearatomic_flag_clear_explicit;
// using VIGRA_THREADING_NAMESPACE::atomic_init;
// using VIGRA_THREADING_NAMESPACE::kill_dependency;
	
}} // namespace vigra::threading

#undef VIGRA_THREADING_NAMESPACE

#endif // VIGRA_THREADING_HXX
