// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <map>
#include <memory>

namespace Acts {

class ISurfaceMaterial;
class IVolumeMaterial;

namespace Experimental {
namespace DetectorMaterialHelper {

using SurfaceMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
using VolumeMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
using DetectorMaterialMaps = std::pair<SurfaceMaterialMaps, VolumeMaterialMaps>;

// A Surface material assigner
struct SurfaceMaterialAssigner {
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
      surfaceMaterialMap;

  /// @brief Call operator compaticle bith the visitMutableSurface interface
  /// @param surface that gets material assigned (if found in the map)
  void operator()(Surface* surface) {
    auto it = surfaceMaterialMap.find(surface->geometryId());
    if (it != surfaceMaterialMap.end()) {
      surface->assignSurfaceMaterial(it->second);
    }
  }
};

struct VolumeMaterialAssigner {
  std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>
      volumeMaterialMap;

  /// @brief Call operator compaticle bith the visitMutableSurface interface
  /// @param volume that gets material assigned (if found in the map)
  void operator()(DetectorVolume* volume) {
    auto it = volumeMaterialMap.find(volume->geometryId());
    if (it != volumeMaterialMap.end()) {
      volume->assignVolumeMaterial(it->second);
    }
  }
};

/// @brief Load all material maps to the detector
///
/// @param detector the detector to be loaded
/// @param materialMaps the material maps to be loaded
///
void assignMaterial(Detector& detector,
                    const DetectorMaterialMaps& materialMaps);

}  // namespace DetectorMaterialHelper
}  // namespace Experimental
}  // namespace Acts
