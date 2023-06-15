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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationState.hpp"

#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;

/// Helper extractors: all portals
struct AllPortalsExtractor {
  /// Extract the portals from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  ///
  /// @return a vector of raw Portal pointers
  inline const std::vector<const Portal*>& extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "AllPortalsExtractor: no detector volume given.");
    }
    return nState.currentVolume->portals();
  }
};

/// Helper extractors: all surfaces
struct AllSurfacesExtractor {
  /// Extract the surfaces from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices is an ignored index vector
  ///
  /// @return a vector of raw Surface pointers
  inline const std::vector<const Surface*>& extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      [[maybe_unused]] const std::vector<std::size_t>& indices = {}) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "AllSurfacesExtractor: no detector volume given.");
    }
    return nState.currentVolume->surfaces();
  }
};

/// Helper extractors: indexed surfaces
struct IndexedSurfacesExtractor {
  /// Extract the surfaces from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the surfaces store
  ///
  /// @note no out of boudns checking is done
  ///
  /// @return a vector of raw Surface pointers
  inline std::vector<const Surface*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      const std::vector<std::size_t>& indices) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "IndexedSurfacesExtractor: no detector volume given.");
    }
    // Get the surface container
    const auto& surfaces = nState.currentVolume->surfaces();
    // The extracted surfaces
    std::vector<const Surface*> eSurfaces;
    eSurfaces.reserve(indices.size());
    std::for_each(indices.begin(), indices.end(),
                  [&](const auto& i) { eSurfaces.push_back(surfaces[i]); });
    return eSurfaces;
  }
};

/// Helper extractors: all sub volumes of a volume
struct AllSubVolumesExtractor {
  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the volume store (ignored)
  ///
  /// @return a vector of raw DetectorVolume pointers
  inline const std::vector<const DetectorVolume*>& extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      [[maybe_unused]] const std::vector<std::size_t>& indices = {}) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "AllSubVolumesExtractor: no detector volume given.");
    }
    return nState.currentVolume->volumes();
  }
};

/// Helper extractors: indexed sub volume of a volume
struct IndexedSubVolumesExtractor {
  /// Extract the sub volumes from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param indices are access indices into the volume store
  ///
  /// @return a vector of raw DetectorVolume pointers
  inline std::vector<const DetectorVolume*> extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState,
      const std::vector<std::size_t>& indices) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "AllSubVolumesExtractor: no detector volume given.");
    }
    // Get the sub volumes container
    const auto& volumes = nState.currentVolume->volumes();
    // The extracted volumes
    std::vector<const DetectorVolume*> eVolumes;
    eVolumes.reserve(indices.size());
    std::for_each(indices.begin(), indices.end(),
                  [&](const auto& i) { eVolumes.push_back(volumes[i]); });
    return eVolumes;
  }
};

/// @brief A helper struct that allows to extrace a volume
/// from the detector by its index
struct IndexedDetectorVolumeExtractor {
  /// Extract the surfaces from the volume
  ///
  /// @param gctx the geometry contextfor this extraction call
  /// @param nState is the current navigation state
  /// @param index is the index in the global detector volume store
  ///
  /// @return a raw DetectorVolume pointer
  inline static const DetectorVolume* extract(
      [[maybe_unused]] const GeometryContext& gctx,
      const NavigationState& nState, std::size_t index) noexcept(false) {
    if (nState.currentDetector == nullptr) {
      throw std::runtime_error("IndexedVolumeExtractor: no detector given.");
    }
    // Get the volume container from the detector
    const auto& volumes = nState.currentDetector->volumes();
    return volumes[index];
  }
};

/// @brief This holds and extracts a collection of detector
/// volumes. Alternatively an extractor could also use the
/// detector and its associated detector volume container for
/// volume extraction.
///
struct DetectorVolumesCollection {
  /// The volumes held by this collection
  std::vector<const DetectorVolume*> dVolumes = {};

  /// Extract a voume from a collection
  ///
  /// @note that geometry context and navigation state are ignored here
  /// @param index are access indices into the surfaces store
  ///
  /// @return the extracted volume
  inline const DetectorVolume* extract(const GeometryContext& /*gctx*/,
                                       const NavigationState& /*nState*/,
                                       const std::size_t& index) const {
    return dVolumes[index];
  }
};

}  // namespace Experimental
}  // namespace Acts
