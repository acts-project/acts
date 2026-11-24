// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

namespace traccc::device {

/// Helper struct to compare surface descriptors by their geo IDs
///
/// This is used in the CFK to find the measurement ranges for each surface
/// efficiently.
///
struct geo_id_surface_comparator {
  template <typename sf_descriptor_t>
  TRACCC_HOST_DEVICE bool operator()(
      const sf_descriptor_t sf_desc,
      const detray::geometry::identifier& geo_id) {
    return sf_desc.identifier() < geo_id;
  }

  template <typename sf_descriptor_t>
  TRACCC_HOST_DEVICE bool operator()(const detray::geometry::identifier& bc,
                                     const sf_descriptor_t sf_desc) {
    return bc < sf_desc.identifier();
  }

};  // struct geo_id_surface_comparator

/// Relation checking that two consecutive measurements are ordered by surface
///
/// The CKF requires the measurements to be sorted in the same order that
/// @c geo_id_surface_comparator uses, to be able to find the measurement
/// ranges for each surface. Being grouped by surface is not enough.
///
struct geo_id_order_relation {
  TRACCC_HOST_DEVICE bool operator()(
      const detray::geometry::identifier& lhs,
      const detray::geometry::identifier& rhs) const {
    return (lhs <= rhs);
  }
};  // struct geo_id_order_relation

}  // namespace traccc::device
