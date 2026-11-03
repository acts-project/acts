// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../sanity/contiguous_on.hpp"
#include "../sanity/ordered_on.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/edm/device/identity_projector.hpp"
#include "traccc/finding/device/geo_id_surface_comparator.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {

bool combinatorial_kalman_filter_algorithm::input_is_valid(
    const edm::measurement_collection::const_view& measurements) const {
  static constexpr std::size_t GEOMID_INDEX = 6u;
  // The measurements must be grouped by surface, and also sorted by it for
  // the CKF to find the measurement ranges of each surface.
  return is_contiguous_on<
             vecmem::device_vector<const detray::geometry::identifier>>(
             device::identity_projector{}, mr().main, copy(), stream(),
             measurements.template get<GEOMID_INDEX>()) &&
         is_ordered_on<
             vecmem::device_vector<const detray::geometry::identifier>>(
             device::geo_id_order_relation{}, mr().main, copy(), stream(),
             measurements.template get<GEOMID_INDEX>());
}

}  // namespace traccc::hip
