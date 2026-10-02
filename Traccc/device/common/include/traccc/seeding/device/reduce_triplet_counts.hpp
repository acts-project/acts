// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/triplet_counter.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {

/// Function used for performing a reduction on the triplet counters
///
/// The @c count_triplets function will only count the number of triplets per
/// each middle spacepoint independently. Thus, we need to perform a sum over
/// all middle spacepoints, as well as correctly claim the positions in which
/// the triplets will be filled in @c find:triplets
///
/// @param[in] globalIndex    The index of the current thread
/// @param[in] dc_view        Collection of doublet counters
/// @param[in,out] spM_tc_view Collection of triplet counters per middle
/// spacepoint
/// @param[out] num_triplets  The total number of triplets
TRACCC_HOST_DEVICE
inline void reduce_triplet_counts(
    global_index_t globalIndex,
    const doublet_counter_collection_types::const_view& dc_view,
    triplet_counter_spM_collection_types::view spM_tc_view,
    unsigned int& num_triplets);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/reduce_triplet_counts.ipp"
