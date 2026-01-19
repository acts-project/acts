// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/IndexGridNavigationPolicy.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {

class TrackingVolume;

namespace IndexGridNavigationJsonConverter {

/// Configure the conversion through that struct
struct Options {
  /// Whether to write the surfaces assigned to the grid
  bool writeSurfaces = true;
  /// Write projected surfaces (i.e. in grid frame for plotting)
  bool writeProjectedSurfaces = false;
  /// Number of polyhedron points to use for projected surface representation
  std::size_t numPolyhedronPoints = 20;
};

/// Convert an IndexGridNavigationPolicy for regular cylinder to Json
/// @param policy The policy to convert
/// @return The resulting json object
nlohmann::json toJson(const RegularCylinderIndexGridNavigationPolicy& policy);

/// Convert an IndexGridNavigationPolicy for regular cylinder to Json
/// @param gctx The geometry context
/// @param policy The policy to convert
/// @param volume The tracking volume holding the surfaces
/// @param options The options for the conversion
/// @return The resulting json object
nlohmann::json toJson(const GeometryContext& gctx,
                      const RegularCylinderIndexGridNavigationPolicy& policy,
                      const TrackingVolume& volume,
                      const Options& options = Options());

/// Convert an IndexGridNavigationPolicy for regular plane to Json
/// @param policy The policy to convert
/// @return The resulting json object
nlohmann::json toJson(const RegularPlaneIndexGridNavigationPolicy& policy);

/// Convert an IndexGridNavigationPolicy for regular plane to Json
/// @param gctx The geometry context
/// @param policy The policy to convert
/// @param volume The tracking volume holding the surfaces
/// @param options The options for the conversion
/// @return The resulting json object
nlohmann::json toJson(const GeometryContext& gctx,
                      const RegularPlaneIndexGridNavigationPolicy& policy,
                      const TrackingVolume& volume,
                      const Options& options = Options());

/// Convert an IndexGridNavigationPolicy for regular disc to Json
/// @param policy The policy to convert
/// @return The resulting json object
nlohmann::json toJson(const RegularDiscIndexGridNavigationPolicy& policy);

/// Convert an IndexGridNavigationPolicy for regular disc to Json
/// @param gctx The geometry context
/// @param policy The policy to convert
/// @param volume The tracking volume holding the surfaces
/// @param options The options for the conversion
/// @return The resulting json object
nlohmann::json toJson(const GeometryContext& gctx,
                      const RegularDiscIndexGridNavigationPolicy& policy,
                      const TrackingVolume& volume,
                      const Options& options = Options());

/// Convert an IndexGridNavigationPolicy for regular ring to Json
/// @param policy The policy to convert
/// @return The resulting json object
nlohmann::json toJson(const RegularRingIndexGridNavigationPolicy& policy);

/// Convert an IndexGridNavigationPolicy for regular ring to Json
/// @param gctx The geometry context
/// @param policy The policy to convert
/// @param volume The tracking volume holding the surfaces
/// @param options The options for the conversion
/// @return The resulting json object
nlohmann::json toJson(const GeometryContext& gctx,
                      const RegularRingIndexGridNavigationPolicy& policy,
                      const TrackingVolume& volume,
                      const Options& options = Options());

}  // namespace IndexGridNavigationJsonConverter
}  // namespace Acts
