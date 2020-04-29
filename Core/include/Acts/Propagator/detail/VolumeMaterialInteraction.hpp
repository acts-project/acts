// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace detail {
/// @brief Struct to handle volume material interaction
struct VolumeMaterialInteraction {
  /// Data from the propagation state
  const TrackingVolume* volume;

  const Vector3D pos;
  const double time;
  const Vector3D dir;
  const double momentum;
  const double q;
  const double qOverP;

  const double mass;
  const int pdg;
  const bool performCovarianceTransport;
  const NavigationDirection nav;

  /// Data evaluated within this struct
  MaterialProperties slab;
  double pathCorrection;

  double variancePhi = 0.;
  double varianceTheta = 0.;
  double varianceQoverP = 0.;

  double Eloss = 0.;
  double nextP;

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
        nav(state.stepping.navDir),
        variancePhi(state.stepping.cov(ePHI, ePHI)),
        varianceTheta(state.stepping.cov(eTHETA, eTHETA)),
        varianceQoverP(state.stepping.cov(eQOP, eQOP)) {}

  /// @brief This function evaluates the material properties to interact with
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param [in] state State of the propagation
  /// @param [in] updateStage The stage of the material update
  ///
  /// @return Boolean statement whether the material is valid
  template <typename propagator_state_t>
  bool evaluateMaterialProperties(const propagator_state_t& state) {
    pathCorrection = 0;
    slab = MaterialProperties(
        state.navigation.currentVolume->volumeMaterial()->material(pos),
        1);  // state.stepping.StepSize
    return slab;
  }
};
}  // namespace detail
}  // end of namespace Acts
