// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include <cmath>

#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_helpers.hpp"

// Detray include(s).
#include <detray/geometry/tracking_surface.hpp>

namespace traccc::details {

/// Intersect strip lines with beam-spot planes, with explicit endpoint
/// allowances.
TRACCC_HOST_DEVICE inline bool make_strip_spacepoint(
    point3& spacepoint, const point3& first_center,
    const vector3& first_direction, const vector3& second_direction,
    const vector3& first_trajectory, const vector3& second_trajectory,
    const vector3& first_normal, const vector3& second_normal,
    const scalar first_half_length, const scalar second_half_length,
    const scalar strip_length_gap_tolerance,
    const scalar strip_length_tolerance) {
  const scalar strip_length_limit = scalar{1} + strip_length_tolerance;
  const scalar first_denominator =
      detray::vector::dot(first_direction, second_normal);
  const scalar second_denominator =
      detray::vector::dot(second_direction, first_normal);
  if ((first_half_length <= 0.f) || (second_half_length <= 0.f) ||
      (std::abs(first_denominator) < 1e-12f) ||
      (std::abs(second_denominator) < 1e-12f)) {
    return false;
  }

  const scalar a = -detray::vector::dot(first_trajectory, second_normal);
  const scalar c = -detray::vector::dot(second_trajectory, first_normal);
  const scalar first_one_over_strip = 0.5f / first_half_length;
  const scalar second_one_over_strip = 0.5f / second_half_length;
  const scalar first_pre_cut =
      strip_length_limit + first_one_over_strip * strip_length_gap_tolerance;
  const scalar second_pre_cut =
      strip_length_limit + second_one_over_strip * strip_length_gap_tolerance;

  if ((std::abs(a) > std::abs(first_denominator) * first_pre_cut) ||
      (std::abs(c) > std::abs(second_denominator) * second_pre_cut)) {
    return false;
  }

  scalar m = a / first_denominator;
  scalar n = c / second_denominator;

  if (strip_length_gap_tolerance != 0.f) {
    const scalar cs = detray::vector::dot(first_direction, second_direction) *
                      first_one_over_strip * first_one_over_strip;
    if (std::abs(cs) < 1e-12f) {
      return false;
    }
    if ((m > strip_length_limit) || (n > strip_length_limit)) {
      scalar dm = m - 1.f;
      const scalar dmn = (n - 1.f) * cs;
      if (dmn > dm) {
        dm = dmn;
      }
      m -= dm;
      n -= dm / cs;
    } else if ((m < -strip_length_limit) || (n < -strip_length_limit)) {
      scalar dm = -(1.f + m);
      const scalar dmn = -(1.f + n) * cs;
      if (dmn > dm) {
        dm = dmn;
      }
      m += dm;
      n += dm / cs;
    }

    if ((std::abs(m) > strip_length_limit) ||
        (std::abs(n) > strip_length_limit)) {
      return false;
    }
  }

  spacepoint = first_center + (0.5f * m) * first_direction;
  return true;
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_measurement(
    const edm::measurement<measurement_backend_t>& meas) {
  // We use 2D (pixel) measurements only for spacepoint creation
  return (meas.dimensions() == 2u);
}

template <typename spacepoint_backend_t, detray::concepts::detector detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_pixel_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& meas,
    const typename detector_t::geometry_context gctx) {
  // Get the global position of this silicon pixel measurement.
  const detray::tracking_surface sf{det, meas.surface_link()};
  const detray::dpoint3D<typename detector_t::algebra_type> global =
      sf.local_to_global(
          gctx,
          edm::get_measurement_local<typename detector_t::algebra_type>(meas),
          {});

  // Fill the spacepoint with the global position and the measurement.
  sp.x() = static_cast<float>(getter::element(global, 0u));
  sp.y() = static_cast<float>(getter::element(global, 1u));
  sp.z() = static_cast<float>(getter::element(global, 2u));
  sp.radius_variance() = 0.f;
  sp.z_variance() = 0.f;
}

}  // namespace traccc::details
