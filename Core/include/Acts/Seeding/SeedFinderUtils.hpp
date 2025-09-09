// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointMutableData.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"

namespace Acts {

/// @brief A partial description of a circle in u-v space.
struct LinCircle {
  LinCircle() = default;
  /// Constructor with circle parameters
  /// @param ct Cotangent of polar angle theta
  /// @param idr Inverse delta radius
  /// @param er Error on radius
  /// @param u U coordinate
  /// @param v V coordinate
  /// @param X X coordinate
  /// @param Y Y coordinate
  LinCircle(float ct, float idr, float er, float u, float v, float X, float Y)
      : cotTheta(ct), iDeltaR(idr), Er(er), U(u), V(v), x(X), y(Y) {}

  /// Cotangent of polar angle theta in coordinate transformation
  float cotTheta{0.};
  /// Inverse delta radius between space points
  float iDeltaR{0.};
  /// Error term in circle fitting for u-v transformation
  float Er{0.};
  /// U coordinate in transformed coordinate system
  float U{0.};
  /// V coordinate in transformed coordinate system
  float V{0.};
  /// X coordinate in local coordinate system
  float x{0.};
  /// Y coordinate in local coordinate system
  float y{0.};
};

/// Transform a single spacepoint to u-v space coordinates
/// @tparam external_spacepoint_t The external spacepoint type
/// @tparam callable_t The callable type for coordinate extraction
/// @param mutableData Container for mutable variables used in seeding
/// @param sp The spacepoint to transform
/// @param spM The middle reference spacepoint
/// @param bottom Whether this is a bottom spacepoint
/// @param extractFunction Function to extract coordinates from spacepoints
/// @return LinCircle representing the transformed coordinates
template <typename external_spacepoint_t, typename callable_t>
LinCircle transformCoordinates(Acts::SpacePointMutableData& mutableData,
                               const external_spacepoint_t& sp,
                               const external_spacepoint_t& spM, bool bottom,
                               callable_t&& extractFunction);

/// @brief Transform a vector of spacepoints to u-v space circles with respect
/// to a given middle spacepoint.
///
/// @tparam external_spacepoint_t The external spacepoint type.
///
/// @param mutableData Container for mutable variables used in the seeding
/// @param[in] vec The list of bottom or top spacepoints
/// @param[in] spM The middle spacepoint.
/// @param[in] bottom Should be true if vec are bottom spacepoints.
/// @param[out] linCircleVec The output vector to write to.
template <typename external_spacepoint_t>
void transformCoordinates(Acts::SpacePointMutableData& mutableData,
                          const std::vector<const external_spacepoint_t*>& vec,
                          const external_spacepoint_t& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec);

/// @brief Check the compatibility of spacepoint coordinates in xyz assuming the Bottom-Middle direction with the strip meassument details
///
/// @tparam external_spacepoint_t The external spacepoint type.
///
/// @param[in] config SeedFinder config containing the delegates to the strip measurement details.
/// @param[in] sp Input space point used in the check.
/// @param[in] spacepointPosition Spacepoint coordinates in xyz plane.
/// @param[out] outputCoordinates The output vector to write to.
/// @returns Boolean that says if spacepoint is compatible with being inside the detector element.
template <typename external_spacepoint_t>
bool xyzCoordinateCheck(
    const Acts::SeedFinderConfig<external_spacepoint_t>& config,
    const external_spacepoint_t& sp, const double* spacepointPosition,
    double* outputCoordinates);

}  // namespace Acts

#include "Acts/Seeding/SeedFinderUtils.ipp"
