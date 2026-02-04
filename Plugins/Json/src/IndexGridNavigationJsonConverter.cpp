// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/IndexGridNavigationJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"

#include <limits>

namespace {

/// Helper to write a projection
template <typename IndexPolicyType>
void writeSurfacesAndProjections(
    nlohmann::json& jPolicy, const Acts::GeometryContext& gctx,
    const IndexPolicyType& policy, const Acts::TrackingVolume& volume,
    const Acts::IndexGridNavigationJsonConverter::Options& options) {
  // In right indexed order - surfaces and projected surfaces (as vertices)
  nlohmann::json jSurfaces = nlohmann::json::array();
  nlohmann::json jProjectedSurfaces = nlohmann::json::array();
  const auto& indexGrid = policy.indexGrid();

  bool writeReferenceRange = false;
  std::array<double, 2u> referenceRange = {
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::lowest()};

  if (options.writeSurfaces || options.writeProjectedSurfaces) {
    for (const auto& surface : volume.surfaces()) {
      if (options.writeSurfaces) {
        nlohmann::json jSurface =
            Acts::SurfaceJsonConverter::toJson(gctx, surface);
        jSurfaces.push_back(jSurface);
      }

      if (options.writeProjectedSurfaces) {
        nlohmann::json jProjectedSurface;
        auto polyhedron =
            surface.polyhedronRepresentation(gctx, options.numPolyhedronPoints);
        for (const auto& vertex : polyhedron.vertices) {
          // Special fixing for ring policies
          if (indexGrid.casts.size() == 1u &&
              indexGrid.casts[0] == Acts::AxisDirection::AxisPhi) {
            std::array<Acts::AxisDirection, 2u> rphi = {
                Acts::AxisDirection::AxisR, Acts::AxisDirection::AxisPhi};
            auto pVertex = Acts::GridAccessHelpers::castPosition<
                Acts::RegularDiscIndexGrid>(indexGrid.transform * vertex, rphi);
            // Update reference range
            referenceRange[0] = std::min(referenceRange[0], pVertex[0]);
            referenceRange[1] = std::max(referenceRange[1], pVertex[0]);
            writeReferenceRange = true;
            // Write the projected vertices in R-Phi
            jProjectedSurface.push_back(pVertex);
          } else {
            // Write the projected vertices
            jProjectedSurface.push_back(
                Acts::GridAccessHelpers::castPosition<decltype(indexGrid.grid)>(
                    indexGrid.transform * vertex, indexGrid.casts));
          }
        }
        jProjectedSurfaces.push_back(jProjectedSurface);
      }
    }
  }
  // Write them if they have content
  if (!jSurfaces.empty()) {
    jPolicy["surfaces"] = jSurfaces;
  }
  if (!jProjectedSurfaces.empty()) {
    jPolicy["projectedSurfaces"] = jProjectedSurfaces;
  }
  if (writeReferenceRange) {
    jPolicy["projectedReferenceRange"] = referenceRange;
  }
}

}  // namespace

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const RegularCylinderIndexGridNavigationPolicy& policy) {
  nlohmann::json jPolicy;

  jPolicy["type"] = "RegularCylinderIndexGridNavigationPolicy";
  jPolicy["grid"] = Acts::GridJsonConverter::toJson(policy.indexGrid().grid);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const GeometryContext& gctx,
    const RegularCylinderIndexGridNavigationPolicy& policy,
    const TrackingVolume& volume, const Options& options) {
  nlohmann::json jPolicy = toJson(policy);
  writeSurfacesAndProjections(jPolicy, gctx, policy, volume, options);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const RegularPlaneIndexGridNavigationPolicy& policy) {
  nlohmann::json jPolicy;

  jPolicy["type"] = "RegularPlaneIndexGridNavigationPolicy";
  jPolicy["grid"] = Acts::GridJsonConverter::toJson(policy.indexGrid().grid);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const GeometryContext& gctx,
    const RegularPlaneIndexGridNavigationPolicy& policy,
    const TrackingVolume& volume, const Options& options) {
  nlohmann::json jPolicy = toJson(policy);
  writeSurfacesAndProjections(jPolicy, gctx, policy, volume, options);
  return jPolicy;
}
nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const RegularDiscIndexGridNavigationPolicy& policy) {
  nlohmann::json jPolicy;

  jPolicy["type"] = "RegularDiscIndexGridNavigationPolicy";
  jPolicy["grid"] = Acts::GridJsonConverter::toJson(policy.indexGrid().grid);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const GeometryContext& gctx,
    const RegularDiscIndexGridNavigationPolicy& policy,
    const TrackingVolume& volume, const Options& options) {
  nlohmann::json jPolicy = toJson(policy);
  writeSurfacesAndProjections(jPolicy, gctx, policy, volume, options);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const RegularRingIndexGridNavigationPolicy& policy) {
  nlohmann::json jPolicy;

  jPolicy["type"] = "RegularRingIndexGridNavigationPolicy";
  jPolicy["grid"] = Acts::GridJsonConverter::toJson(policy.indexGrid().grid);
  return jPolicy;
}

nlohmann::json Acts::IndexGridNavigationJsonConverter::toJson(
    const GeometryContext& gctx,
    const RegularRingIndexGridNavigationPolicy& policy,
    const TrackingVolume& volume, const Options& options) {
  nlohmann::json jPolicy = toJson(policy);
  writeSurfacesAndProjections(jPolicy, gctx, policy, volume, options);
  return jPolicy;
}
