// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vecmem/containers/data/vector_view.hpp>

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/global_index.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE inline void find_strip_pairs(
    global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos_view,
    const strip_pairing_rule_collection_types::const_view& rules_view,
    const point3& beam_spot, unsigned int& opposite_position,
    unsigned int& overlap_position,
    strip_pair_collection_types::view opposite_pairs_view,
    strip_pair_collection_types::view overlap_pairs_view);

}  // namespace traccc::device

#include "traccc/seeding/device/impl/find_strip_pairs.ipp"
