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

#ifndef VIGRA_ATOMIC_HXX
#define VIGRA_ATOMIC_HXX

/* Compatibility header to import atomic operation support from boost::atomic
   when std::atomic is not yet available.
*/

#ifdef VIGRA_NO_STD_ATOMIC
# include <boost/atomic.hpp>
# define VIGRA_ATOMIC_NAMESPACE boost
#else
# include <atomic>
# define VIGRA_ATOMIC_NAMESPACE std
#endif

namespace vigra { namespace atomic {

using VIGRA_ATOMIC_NAMESPACE::atomic_flag;
using VIGRA_ATOMIC_NAMESPACE::atomic;

using VIGRA_ATOMIC_NAMESPACE::atomic_char;
using VIGRA_ATOMIC_NAMESPACE::atomic_schar;
using VIGRA_ATOMIC_NAMESPACE::atomic_uchar;
using VIGRA_ATOMIC_NAMESPACE::atomic_short;
using VIGRA_ATOMIC_NAMESPACE::atomic_ushort;
using VIGRA_ATOMIC_NAMESPACE::atomic_int;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint;
using VIGRA_ATOMIC_NAMESPACE::atomic_long;
using VIGRA_ATOMIC_NAMESPACE::atomic_ulong;
using VIGRA_ATOMIC_NAMESPACE::atomic_llong;
using VIGRA_ATOMIC_NAMESPACE::atomic_ullong;
using VIGRA_ATOMIC_NAMESPACE::atomic_char16_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_char32_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_wchar_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_least8_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_least8_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_least16_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_least16_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_least32_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_least32_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_least64_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_least64_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_fast8_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_fast8_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_fast16_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_fast16_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_fast32_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_fast32_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_int_fast64_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uint_fast64_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_intptr_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uintptr_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_size_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_ptrdiff_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_intmax_t;
using VIGRA_ATOMIC_NAMESPACE::atomic_uintmax_t;

using VIGRA_ATOMIC_NAMESPACE::memory_order_relaxed;
using VIGRA_ATOMIC_NAMESPACE::memory_order_release;
using VIGRA_ATOMIC_NAMESPACE::memory_order_acquire;
using VIGRA_ATOMIC_NAMESPACE::memory_order_consume;
using VIGRA_ATOMIC_NAMESPACE::memory_order_acq_rel;
using VIGRA_ATOMIC_NAMESPACE::memory_order_seq_cst;

using VIGRA_ATOMIC_NAMESPACE::atomic_thread_fence;
using VIGRA_ATOMIC_NAMESPACE::atomic_signal_fence;
	
}} // namespace vigra::atomic

#undef VIGRA_ATOMIC_NAMESPACE

#endif // VIGRA_ATOMIC_HXX
