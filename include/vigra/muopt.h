/*
 * Copyright (C) 2012-2014 Dr. Christoph L. Spiel
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
