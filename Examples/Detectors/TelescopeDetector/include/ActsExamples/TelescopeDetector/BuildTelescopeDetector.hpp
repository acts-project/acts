// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <array>
#include <memory>
#include <vector>

namespace ActsExamples {

/// The telescope detector surface type
enum class TelescopeSurfaceType {
  Plane = 0,
  Disc = 1,
};

/// Global method to build the telescope tracking geometry
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param positions are the positions of different layers in the longitudinal
///                  direction
/// @param stereoAngles are the stereo angles of different layers, which are
///                     rotation angles around the longitudinal (normal)
///                     direction
/// @param offsets is the offset (u, v) of the layers in the transverse plane
/// @param bounds is the surface bound values, i.e. halfX and halfY if plane
///               surface, and minR and maxR if disc surface
/// @param thickness is the material thickness of each layer
/// @param surfaceType is the detector surface type
/// @param rotDirection indicates which axis the detector surface normals are
/// parallel to
std::unique_ptr<const Acts::TrackingGeometry> buildTelescopeDetector(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<const Acts::SurfacePlacementBase>>&
        detectorStore,
    const std::vector<double>& positions,
    const std::vector<double>& stereoAngles,
    const std::array<double, 2>& offsets, const std::array<double, 2>& bounds,
    double thickness, TelescopeSurfaceType surfaceType,
    Acts::AxisDirection rotDirection = Acts::AxisDirection::AxisZ);

/// Global method to build the telescope Gen3 tracking geometry
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param positions are the positions of different layers in the longitudinal
///                  direction
/// @param stereoAngles are the stereo angles of different layers, which are
///                     rotation angles around the longitudinal (normal)
///                     direction
/// @param offsets is the offset (u, v) of the layers in the transverse plane
/// @param bounds is the surface bound values, i.e. halfX and halfY if plane
///               surface, and minR and maxR if disc surface
/// @param thickness is the material thickness of each layer
/// @param surfaceType is the detector surface type
/// @param rotDirection indicates which axis the detector surface normals are
///                     parallel to
/// @param envelope_x is the world size along the x-axis
/// @param envelope_y is the world size along the y-axis
/// @param envelope_z is the world size along the z-axis
/// @param logger is the Logger that (with different level of detail) describes
//                what happens in the building if the tracking geometry
std::unique_ptr<const Acts::TrackingGeometry> buildTelescopeDetectorGen3(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<const Acts::SurfacePlacementBase>>&
        detectorStore,
    const std::vector<double>& positions,
    const std::vector<double>& stereoAngles,
    const std::array<double, 2>& offsets, const std::array<double, 2>& bounds,
    double thickness, TelescopeSurfaceType surfaceType,
    Acts::AxisDirection rotDirection, const std::array<double, 2>& envelope_x,
    const std::array<double, 2>& envelope_y,
    const std::array<double, 2>& envelope_z, const Acts::Logger& logger);

}  // namespace ActsExamples
