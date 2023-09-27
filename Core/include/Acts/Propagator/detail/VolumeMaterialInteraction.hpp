// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace detail {

/// @brief Struct to handle volume material interaction
struct VolumeMaterialInteraction {
  /// Data from the propagation state
  const TrackingVolume* volume;
  /// The particle current position
  const Vector3 pos;
  /// The particle current time
  const double time;
  /// The particle current direction
  const Vector3 dir;
  /// The particle current momentum
  const double momentum;
  /// The particle current charge
  const double q;
  /// The particle current q/p
  const double qOverP;
  /// The particle mass
  const double mass;
  /// The particle pdg
  const int pdg;
  /// The covariance transport decision at the interaction
  const bool performCovarianceTransport;
  /// The navigation direction
  const NavigationDirection nav;

  /// Data evaluated within this struct
  MaterialSlab slab;
  /// The path correction factor due to non-zero incidence on the surface.
  double pathCorrection = 0;

  /// @brief Constructor
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] vVolume The current volume
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  VolumeMaterialInteraction(const TrackingVolume* vVolume,
                            const propagator_state_t& state,
                            const stepper_t& stepper)
      : volume(vVolume),
        pos(stepper.position(state.stepping)),
        time(stepper.time(state.stepping)),
        dir(stepper.direction(state.stepping)),
        momentum(stepper.momentum(state.stepping)),
        q(stepper.charge(state.stepping)),
        qOverP(q / momentum),
        mass(state.options.mass),
        pdg(state.options.absPdgCode),
        performCovarianceTransport(state.stepping.covTransport),
        nav(state.stepping.navDir) {}

  /// @brief This function evaluates the material properties to interact with
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam navigator_t Type of the propagator state
  ///
  /// @param [in] state State of the propagation
  /// @param [in] navigator Navigator of the propagation
  ///
  /// @return Boolean statement whether the material is valid
  template <typename propagator_state_t, typename navigator_t>
  bool evaluateMaterialSlab(const propagator_state_t& state,
                            const navigator_t& navigator) {
    pathCorrection = 0;
    if (navigator.currentVolume(state.navigation) != nullptr &&
        navigator.currentVolume(state.navigation)->volumeMaterial() !=
            nullptr) {
      slab = MaterialSlab(navigator.currentVolume(state.navigation)
                              ->volumeMaterial()
                              ->material(pos),
                          1);  // state.stepping.StepSize
    } else {
      slab = MaterialSlab();
    }
    return slab;
  }
};

}  // namespace detail
}  // end of namespace Acts
