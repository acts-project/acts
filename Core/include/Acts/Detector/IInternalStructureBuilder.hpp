// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
class Surface;
}

namespace Acts {
namespace Experimental {

class DetectorVolume;

/// Holder struct for the internal structure components of a DetectorVolume
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
  /// Navigaiton delegate for voluems
  DetectorVolumeUpdator volumeUpdator;
};

/// @brief This is the interface definition of internal structure
/// builders for DetectorVolume construction.
///
/// It is assumed that each builder returns a consistent set of
/// DetectorVolume internals, which in turn can be directly provided
/// to a DetectorVolume constructor.
class IInternalStructureBuilder {
 public:
  virtual ~IInternalStructureBuilder() = default;
  /// The interface definition for internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  virtual InternalStructure create(const GeometryContext& gctx) const = 0;
};

}  // namespace Experimental
}  // namespace Acts
