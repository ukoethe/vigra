/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe, Christoph L. Spiel */
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

#ifndef __MUOPT_H__
#define __MUOPT_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef __GNUC__

#define PREFETCH(m_addr) __builtin_prefetch(m_addr)
#define HINTED_PREFETCH(m_addr, m_rw_hint, m_temporal_locality_hint) \
    __builtin_prefetch((m_addr), (m_rw_hint), (m_temporal_locality_hint))

#define EXPECT_RESULT(m_condition, m_expected_result) \
    __builtin_expect((m_condition), static_cast<int>(m_expected_result))

#else

#define PREFETCH(m_addr)
#define HINTED_PREFETCH(m_addr, m_rw_hint, m_temporal_locality_hint)

#define EXPECT_RESULT(m_condition, m_expected_result) (m_condition)

#endif


typedef enum {
    PREPARE_FOR_READ,
    PREPARE_FOR_WRITE
} rw_hint;


// If data is only touched once, or if the dataset is smaller than the
// cache, prefer the non-temporal version; otherwise use one of the
// temporal versions.
typedef enum {
    // Fetch data into the first way of the L1/L2 cache, minimizing cache pollution.
    NO_TEMPORAL_LOCALITY,

    // Fetch data into the least-recently-used way of the ...
    LOW_TEMPORAL_LOCALITY,      // ... L3 cache?
    MEDIUM_TEMPORAL_LOCALITY,   // ... L2/L3 cache?
    HIGH_TEMPORAL_LOCALITY      // ... L1/L2/L3 cache just as a normal load would do.
} temporal_locality_hint;


#endif /* __MUOPT_H__ */

// Local Variables:
// mode: c++
// End:
