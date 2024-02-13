// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Definitions/Algebra.hpp"

#include <sstream>

namespace Acts {

/// Simple struct to select volumes
struct VolumeSelector {
  bool selectMaterial = true;
  bool selectLayer = false;
  bool selectPassive = false;

  /// VolumeSelector with options
  ///
  /// @param sMaterial is the directive to select material volumes
  /// @param sLayer is the directive to select volumes with layers
  /// @param sPassive is the directive to select passive volumes
  VolumeSelector(bool sMaterial = true, bool sLayer = false,
                 bool sPassive = false)
      : selectMaterial(sMaterial),
        selectLayer(sLayer),
        selectPassive(sPassive) {}

  /// Call operator to check if a volume should be selected
  ///
  /// @param volume is the test volume
  bool operator()(const Acts::TrackingVolume& volume) const {
    if (selectMaterial && volume.volumeMaterial() != nullptr) {
      return true;
    }
    if (selectLayer && volume.confinedLayers() != nullptr) {
      return true;
    }
    if (selectPassive) {
      return true;
    }
    return false;
  }
};

/// The information to be writtern out per hit volume
struct HitVolume : public InteractionVolume {
  /// Empty constructor
  HitVolume() = delete;

  /// Constructor from a tracking volume
  ///
  /// @param tv The tracking volume
  /// @param pos The position of the hit
  /// @param dir The direction of the hit
  HitVolume(const TrackingVolume* tv, 
  Vector3 pos, 
  Vector3 dir) : InteractionVolume(tv), position(pos), direction(dir) {};

  /// Constructor from a detector volume
  ///
  /// @param dv The detector volume
  /// @param pos The position of the hit
  /// @param dir The direction of the hit
  HitVolume(const Experimental::DetectorVolume* dv,
  Vector3 pos,
  Vector3 dir) : InteractionVolume(dv), position(pos), direction(dir) {};

  bool inside(const GeometryContext& gctx, 
    const Acts::Vector3& gpos,
    double tol = 0.0) const {
      if (trackingVolume != nullptr) {
        return trackingVolume->inside(gpos, tol);
      } else {
        return detectorVolume->inside(gctx, gpos);
      }
  };

  Vector3 position;
  Vector3 direction;
};

/// A Volume Collector struct
/// templated with a Selector type
///
/// Whenever a volume is passed in the propagation
/// that satisfies the selector, it is recorded
/// for further usage in the flow.
// template <typename volume_type, typename Selector = VolumeSelector>
template <typename Selector = VolumeSelector>
struct VolumeCollector {
  /// The selector used for this volume
  Selector selector;

  /// Simple result struct to be returned
  /// It has all the HitVolume objects that
  /// are collected (and thus have been selected)
  struct this_result {
    // std::vector<HitVolume<volume_type>> collected;
    std::vector<HitVolume> collected;
  };

  using result_type = this_result;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the propagator state has a current volume,
  /// in which case the action is performed:
  /// - it records the volume given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper used for the propagation
  /// @tparam navigator_t Type of the navigator used for the propagation
  ///
  /// @param [in,out] state is the mutable stepper state object
  /// @param [in] stepper The stepper in use
  /// @param [in] navigator The navigator in use
  /// @param [in,out] result is the mutable result object
  /// @param logger the logger object
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, result_type& result,
                  const Logger& logger) const {
    auto currentVolume = navigator.currentVolume(state.navigation);

    // The current volume has been assigned by the navigator
    if (currentVolume && selector(*currentVolume)) {
      // Create for recording
      HitVolume volume_hit(currentVolume,
                           stepper.position(state.stepping),
                           stepper.direction(state.stepping));
      bool save = true;
      // Check if the Volume ws already encountered
      for (auto const& res : result.collected) {
        if (res.geometryId() == volume_hit.geometryId()) {
          save = false;
          break;
        }
      }
      // Save if in the result if it does not already exist
      if (save) {
        result.collected.push_back(volume_hit);
        // Screen output
        ACTS_VERBOSE("Collect volume  " << currentVolume->geometryId());
      }
    }
  }
};

}  // namespace Acts
