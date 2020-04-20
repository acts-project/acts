// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// The information to be writtern out per hit surface
struct MaterialHit {
  const Surface* surface = nullptr;
  Vector3D position;
  Vector3D direction;
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
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the result object to be filled
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    if (state.navigation.currentSurface) {
      if (state.navigation.currentSurface == state.navigation.targetSurface and
          not state.navigation.targetReached) {
        return;
      }

      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Material check on surface ";
        dstream << state.navigation.currentSurface->geoID();
        return dstream.str();
      });

      if (state.navigation.currentSurface->surfaceMaterial()) {
        // get the material propertices and only continue
        const MaterialProperties* mProperties =
            state.navigation.currentSurface->surfaceMaterial()->material(
                stepper.position(state.stepping));
        if (mProperties) {
          // pre/post/full update
          double prepofu = 1.;
          if (state.navigation.startSurface ==
              state.navigation.currentSurface) {
            debugLog(state, [&] {
              return std::string("Update on start surface: post-update mode.");
            });
            prepofu =
                state.navigation.currentSurface->surfaceMaterial()->factor(
                    state.stepping.navDir, postUpdate);
          } else if (state.navigation.targetSurface ==
                     state.navigation.currentSurface) {
            debugLog(state, [&] {
              return std::string("Update on target surface: pre-update mode");
            });
            prepofu =
                state.navigation.currentSurface->surfaceMaterial()->factor(
                    state.stepping.navDir, preUpdate);
          } else {
            debugLog(state, [&] {
              return std::string("Update while pass through: full mode.");
            });
          }

          // the pre/post factor has been applied
          // now check if there's still something to do
          if (prepofu == 0.) {
            debugLog(state, [&] {
              return std::string("Pre/Post factor set material to zero.");
            });
            return;
          }
          // more debugging output to the screen
          debugLog(state, [&] {
            return std::string("Material properties found for this surface.");
          });

          // the path correction from the surface intersection
          double pCorrection =
              prepofu * state.navigation.currentSurface->pathCorrection(
                            stepper.position(state.stepping),
                            stepper.direction(state.stepping));
          // the full material
          result.materialInX0 += pCorrection * mProperties->thicknessInX0();
          result.materialInL0 += pCorrection * mProperties->thicknessInL0();

          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "t/X0 (t/L0) increased to ";
            dstream << result.materialInX0 << " (";
            dstream << result.materialInL0 << " )";
            return dstream.str();
          });

          // if configured, record the individual material hits
          if (detailedCollection) {
            // create for recording
            MaterialHit mHit;
            mHit.surface = state.navigation.currentSurface;
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

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /*state*/) const {}

 private:
  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param state the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a streamable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(state.options.debugPfxWidth);
      dstream << "material collector"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};
}  // namespace Acts
