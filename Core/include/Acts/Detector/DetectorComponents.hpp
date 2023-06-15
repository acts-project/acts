// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"

#include <map>
#include <memory>
#include <vector>

namespace Acts {

class VolumeBounds;
class Surface;

namespace Experimental {

class DetectorVolume;
class Portal;

/// @brief The currently built detector components
/// including the constructed volumes and the current
/// shell/coating, i.e. portals ordered in a map
///
/// @note the map indices of the shell coating represent
/// their respective index in the portal vector of the
/// VolumeBounds derivative that is described by the given
/// component.
struct DetectorComponent {
  using PortalContainer = std::map<unsigned int, std::shared_ptr<Portal>>;
  /// The vector of construced volume(s)
  std::vector<std::shared_ptr<DetectorVolume>> volumes = {};
  /// The current map of outside portals
  PortalContainer portals;
};

/// @brief Holder struct for the external structure components
/// required to construct a detectorVolume
struct ExternalStructure {
  /// A 3D transform where the volume should be positioned
  Transform3 transform = Transform3::Identity();
  /// A shape definition
  std::unique_ptr<VolumeBounds> bounds = nullptr;
  /// And a portal generator
  PortalGenerator portalGenerator;
};

/// @brief Holder struct for the internal structure components of a DetectorVolume
///
/// @note the surface surfacesUpdator needs to handle also portal providing
/// of contained volumes.
struct InternalStructure {
  /// Contained surfaces of this volume, handled by the surfacesUpdator
  std::vector<std::shared_ptr<Surface>> surfaces = {};
  /// Contained volumes of this volume, handled by the volumeUpdator
  std::vector<std::shared_ptr<DetectorVolume>> volumes = {};
  /// Navigation delegate for surfaces
  SurfaceCandidatesUpdator surfacesUpdator;
  // Navigation delegate for volumes
  DetectorVolumeUpdator volumeUpdator;
};

/// @brief Container to collect root volumes for the
/// construction of a Detector
///
/// @note each root volume is expected to contain a full
/// search tree for eventually contained sub volumes
///
/// This struct collects all root volumes that will then
/// be provided to the Detector object
struct RootDetectorVolumes {
  std::vector<std::shared_ptr<DetectorVolume>> volumes = {};
};

}  // namespace Experimental
}  // namespace Acts
