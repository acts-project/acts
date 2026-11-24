// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <detray/core/concepts.hpp>
#include <detray/definitions/track_parametrization.hpp>

#include "traccc/definitions/primitives.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"

namespace traccc::details {

/// Function helping with checking a measurement object for spacepoint creation
///
/// @param[in]  measurement The input measurement
template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_measurement(
    const edm::measurement<measurement_backend_t>& meas);

/// Fill a spacepoint object with the information from a measurement
///
/// @param[out] sp          The spacepoint to fill
/// @param[in]  det         The tracking geometry
/// @param[in]  measurement The measurement to create the spacepoint out of
/// @param[in]  gctx        The current geometry context
///
template <typename spacepoint_backend_t, detray::concepts::detector detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_pixel_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& meas,
    const typename detector_t::geometry_context gctx = {});

/// Intersect strip lines with beam-spot planes, with explicit endpoint
/// allowances.
TRACCC_HOST_DEVICE inline bool make_strip_spacepoint(
    point3& spacepoint, const point3& first_center,
    const vector3& first_direction, const vector3& second_direction,
    const vector3& first_trajectory, const vector3& second_trajectory,
    const vector3& first_normal, const vector3& second_normal,
    const scalar first_half_length, const scalar second_half_length,
    const scalar strip_length_gap_tolerance,
    const scalar strip_length_tolerance);

}  // namespace traccc::details

// Include the implementation.
#include "traccc/seeding/impl/spacepoint_formation.ipp"
