// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <exception>

namespace Acts {
namespace Experimental {

struct NoopImpl : public INavigationDelegate {
  inline void update(const GeometryContext& /*gctx*/,
                     NavigationState& /*nState*/) const {}
};

/// @brief The end of world sets the volume pointer of the
/// navigation state to nullptr, usually indicates the end of
/// the known world, hence the name
struct TrialAndErrorImpl : public INavigationDelegate {
  /// @brief a null volume link - explicitely
  ///
  /// @param gctx the geometry context for this call
  /// @param nState the navigation state into which the volume is set
  inline void update(const GeometryContext& gctx,
                     NavigationState& nState) const {
    if (nState.currentDetector == nullptr) {
      throw std::runtime_error(
          "DetectorVolumeFinders: no detector set to navigation state.");
    }

    const auto& volumes = nState.currentDetector->volumes();
    for (const auto v : volumes) {
      if (v->inside(gctx, nState.position)) {
        nState.currentVolume = v;
        v->detectorVolumeUpdator()(gctx, nState);
        return;
      }
    }
    nState.currentVolume = nullptr;
  }
};

struct TrialAndErrorVolumeImpl : public INavigationDelegate {
  inline void update(const GeometryContext& gctx,
                     NavigationState& nState) const {
    if (nState.currentVolume == nullptr) {
      throw std::runtime_error(
          "DetectorVolumeFinders: no detector set to navigation state.");
    }

    const auto& volumes = nState.currentVolume->volumes();
    for (const auto v : volumes) {
      if (v->inside(gctx, nState.position)) {
        nState.currentVolume = v;
        v->detectorVolumeUpdator()(gctx, nState);
        return;
      }
    }
  }
};

/// Generate a delegate to try all volume
inline static DetectorVolumeUpdator tryAllVolumes() {
  DetectorVolumeUpdator vFinder;
  auto tae = std::make_unique<const TrialAndErrorImpl>();
  vFinder.connect<&TrialAndErrorImpl::update>(std::move(tae));
  return vFinder;
}

/// Generate a delegate to try all volume
inline static DetectorVolumeUpdator tryAllSubVolumes() {
  DetectorVolumeUpdator vFinder;
  auto tae = std::make_unique<const TrialAndErrorVolumeImpl>();
  vFinder.connect<&TrialAndErrorVolumeImpl::update>(std::move(tae));
  return vFinder;
}

/// Generate a delegate to try no volume
inline static DetectorVolumeUpdator tryNoVolumes() {
  DetectorVolumeUpdator vFinder;
  auto noop = std::make_unique<const NoopImpl>();
  vFinder.connect<&NoopImpl::update>(std::move(noop));
  return vFinder;
}

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
      const NavigationState& nState, size_t index) noexcept(false) {
    if (nState.currentDetector == nullptr) {
      throw std::runtime_error("IndexedVolumeExtractor: no detector given.");
    }
    // Get the volume container from the detector
    const auto& volumes = nState.currentDetector->volumes();
    return volumes[index];
  }
};

/// @brief  An indexed volume implementation access
///
/// @tparam grid_type is the grid type used for this
template <typename grid_type>
using IndexedDetectorVolumeImpl =
    IndexedUpdatorImpl<grid_type, IndexedDetectorVolumeExtractor,
                       DetectorVolumeFiller>;

}  // namespace Experimental
}  // namespace Acts
