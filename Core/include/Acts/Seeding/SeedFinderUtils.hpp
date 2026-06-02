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

/// Transform a single space point to u-v space coordinates
/// @tparam external_space_point_t The external space point type
/// @tparam callable_t The callable type for coordinate extraction
/// @param mutableData Container for mutable variables used in seeding
/// @param sp The space point to transform
/// @param spM The middle reference space point
/// @param bottom Whether this is a bottom space point
/// @param extractFunction Function to extract coordinates from space points
/// @return LinCircle representing the transformed coordinates
template <typename external_space_point_t, typename callable_t>
LinCircle transformCoordinates(Acts::SpacePointMutableData& mutableData,
                               const external_space_point_t& sp,
                               const external_space_point_t& spM, bool bottom,
                               callable_t&& extractFunction);

/// @brief Transform a vector of space points to u-v space circles with respect
/// to a given middle space point.
///
/// @tparam external_space_point_t The external space point type.
///
/// @param mutableData Container for mutable variables used in the seeding
/// @param[in] vec The list of bottom or top space points
/// @param[in] spM The middle space point.
/// @param[in] bottom Should be true if vec are bottom space points.
/// @param[out] linCircleVec The output vector to write to.
template <typename external_space_point_t>
void transformCoordinates(Acts::SpacePointMutableData& mutableData,
                          const std::vector<const external_space_point_t*>& vec,
                          const external_space_point_t& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec);

/// @brief Check the compatibility of space point coordinates in xyz assuming the Bottom-Middle direction with the strip meassument details
///
/// @tparam external_space_point_t The external space point type.
///
/// @param[in] config SeedFinder config containing the delegates to the strip measurement details.
/// @param[in] sp Input space point used in the check.
/// @param[in] spacePointPosition Space point coordinates in xyz plane.
/// @param[out] outputCoordinates The output vector to write to.
/// @returns Boolean that says if space point is compatible with being inside the detector element.
template <typename external_space_point_t>
bool xyzCoordinateCheck(
    const Acts::SeedFinderConfig<external_space_point_t>& config,
    const external_space_point_t& sp, const double* spacePointPosition,
    double* outputCoordinates);

}  // namespace Acts

#include "Acts/Seeding/SeedFinderUtils.ipp"
