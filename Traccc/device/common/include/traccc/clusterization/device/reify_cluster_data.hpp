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

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Fill a cluster collection with the cell indices
///
/// @param thread_id The thread identifier
/// @param disjoint_set_view The cluster/measurement index of each cell
/// @param cluster_view The collection to fill
///
TRACCC_HOST_DEVICE inline void reify_cluster_data(
    global_index_t thread_id,
    vecmem::data::vector_view<const unsigned int> disjoint_set_view,
    vecmem::data::vector_view<const unsigned int> permutation_map_view,
    traccc::edm::silicon_cluster_collection::view cluster_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/clusterization/device/impl/reify_cluster_data.ipp"
