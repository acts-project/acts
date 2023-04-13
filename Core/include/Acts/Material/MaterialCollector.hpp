// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <sstream>

namespace Acts {

/// The information to be writtern out per hit surface
struct MaterialHit {
  const Surface* surface = nullptr;
  Vector3 position;
  Vector3 direction;
  Material material;
  double pathLength;
};

/// A Material Collector struct
struct MaterialCollector {
  /// In the detailed collection mode the material
  /// per surface is collected, otherwise only the total
  /// pathlength in X0 or L0 are recorded
  bool detailedCollection = false;

  /// Simple result struct to be returned
  ///
  /// Result of the material collection process
  /// It collects the overall X0 and L0 path lengths,
  /// and optionally a detailed per-material breakdown
  struct this_result {
    std::vector<MaterialHit> collected;
    double materialInX0 = 0.;
    double materialInL0 = 0.;
  };

  using result_type = this_result;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the state has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper of the propagation
  /// @tparam navigator_t Type of the navigator of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param navigator The navigator in use
  /// @param result is the result object to be filled
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, result_type& result,
                  const Logger& logger) const {
    if (navigator.currentSurface(state.navigation)) {
      if (navigator.currentSurface(state.navigation) ==
              navigator.targetSurface(state.navigation) and
          not navigator.targetReached(state.navigation)) {
        return;
      }

      ACTS_VERBOSE("Material check on surface "
                   << navigator.currentSurface(state.navigation)->geometryId());

      if (navigator.currentSurface(state.navigation)->surfaceMaterial()) {
        // get the material propertices and only continue
        const MaterialSlab* mProperties =
            navigator.currentSurface(state.navigation)
                ->surfaceMaterial()
                ->material(stepper.position(state.stepping));
        if (mProperties) {
          // pre/post/full update
          double prepofu = 1.;
          if (navigator.startSurface(state.navigation) ==
              navigator.currentSurface(state.navigation)) {
            ACTS_VERBOSE("Update on start surface: post-update mode.");
            prepofu = navigator.currentSurface(state.navigation)
                          ->surfaceMaterial()
                          ->factor(state.stepping.navDir,
                                   MaterialUpdateStage::PostUpdate);
          } else if (navigator.targetSurface(state.navigation) ==
                     navigator.currentSurface(state.navigation)) {
            ACTS_VERBOSE("Update on target surface: pre-update mode");
            prepofu = navigator.currentSurface(state.navigation)
                          ->surfaceMaterial()
                          ->factor(state.stepping.navDir,
                                   MaterialUpdateStage::PreUpdate);
          } else {
            ACTS_VERBOSE("Update while pass through: full mode.");
          }

          // the pre/post factor has been applied
          // now check if there's still something to do
          if (prepofu == 0.) {
            ACTS_VERBOSE("Pre/Post factor set material to zero.");
            return;
          }
          // more debugging output to the screen
          ACTS_VERBOSE("Material properties found for this surface.");

          // the path correction from the surface intersection
          double pCorrection =
              prepofu * navigator.currentSurface(state.navigation)
                            ->pathCorrection(stepper.position(state.stepping),
                                             stepper.direction(state.stepping));
          // the full material
          result.materialInX0 += pCorrection * mProperties->thicknessInX0();
          result.materialInL0 += pCorrection * mProperties->thicknessInL0();

          ACTS_VERBOSE("t/X0 (t/L0) increased to "
                       << result.materialInX0 << " (" << result.materialInL0
                       << " )");

          // if configured, record the individual material hits
          if (detailedCollection) {
            // create for recording
            MaterialHit mHit;
            mHit.surface = navigator.currentSurface(state.navigation);
            mHit.position = stepper.position(state.stepping);
            mHit.direction = stepper.direction(state.stepping);
            // get the material & path length
            mHit.material = mProperties->material();
            mHit.pathLength = pCorrection * mProperties->thickness();
            // save if in the result
            result.collected.push_back(mHit);
          }
        }
      }
    }
  }
};  // namespace Acts
}  // namespace Acts
