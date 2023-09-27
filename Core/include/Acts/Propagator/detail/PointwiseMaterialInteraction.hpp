// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace detail {
/// @brief Struct to handle pointwise material interaction
struct PointwiseMaterialInteraction {
  /// Data from the propagation state
  const Surface* surface;

  /// The particle position at the interaction.
  const Vector3 pos = Vector3(0., 0., 0);
  /// The particle time at the interaction.
  const double time = 0.0;
  /// The particle direction at the interaction.
  const Vector3 dir = Vector3(0., 0., 0);
  /// The particle momentum at the interaction
  const double momentum;
  /// The particle charge
  const double q;
  /// The particle q/p at the interaction
  const double qOverP;
  /// The particle mass
  const double mass;
  /// The particle pdg
  const int pdg;
  /// The covariance transport decision at the interaction
  const bool performCovarianceTransport;
  /// The navigation direction
  const NavigationDirection nav;

  /// The effective, passed material properties including the path correction.
  MaterialSlab slab;
  /// The path correction factor due to non-zero incidence on the surface.
  double pathCorrection = 0.;
  /// Expected phi variance due to the interactions.
  double variancePhi = 0.;
  /// Expected theta variance due to the interactions.
  double varianceTheta = 0.;
  /// Expected q/p variance due to the interactions.
  double varianceQoverP = 0.;
  /// The energy change due to the interaction.
  double Eloss = 0.;
  /// The momentum after the interaction
  double nextP = 0.;

  /// @brief Contructor
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] sSurface The current surface
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  PointwiseMaterialInteraction(const Surface* sSurface,
                               const propagator_state_t& state,
                               const stepper_t& stepper)
      : surface(sSurface),
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
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in] state State of the propagation
  /// @param [in] navigator Navigator of the propagation
  /// @param [in] updateStage The stage of the material update
  ///
  /// @return Boolean statement whether the material is valid
  template <typename propagator_state_t, typename navigator_t>
  bool evaluateMaterialSlab(
      const propagator_state_t& state, const navigator_t& navigator,
      MaterialUpdateStage updateStage = MaterialUpdateStage::FullUpdate) {
    // We are at the start surface
    if (surface == navigator.startSurface(state.navigation)) {
      updateStage = MaterialUpdateStage::PostUpdate;
      // Or is it the target surface ?
    } else if (surface == navigator.targetSurface(state.navigation)) {
      updateStage = MaterialUpdateStage::PreUpdate;
    }

    // Retrieve the material properties
    slab = navigator.currentSurface(state.navigation)
               ->surfaceMaterial()
               ->materialSlab(pos, nav, updateStage);

    // Correct the material properties for non-zero incidence
    pathCorrection = surface->pathCorrection(state.geoContext, pos, dir);
    slab.scaleThickness(pathCorrection);

    // Get the surface material & properties from them
    return slab;
  }

  /// @brief This function evaluate the material effects
  ///
  /// @param [in] multipleScattering Boolean to indiciate the application of
  /// multiple scattering
  /// @param [in] energyLoss Boolean to indiciate the application of energy loss
  void evaluatePointwiseMaterialInteraction(bool multipleScattering,
                                            bool energyLoss);

  /// @brief Update the state
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper in use
  /// @param [in] updateMode The noise update mode (in default: add noise)
  template <typename propagator_state_t, typename stepper_t>
  void updateState(propagator_state_t& state, const stepper_t& stepper,
                   NoiseUpdateMode updateMode = addNoise) {
    // in forward(backward) propagation, energy decreases(increases) and
    // variances increase(decrease)
    const auto nextE =
        std::sqrt(mass * mass + momentum * momentum) -
        std::copysign(
            Eloss,
            static_cast<std::underlying_type_t<NavigationDirection>>(nav));
    // put particle at rest if energy loss is too large
    nextP = (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;
    // minimum momentum below which we will not push particles via material
    // update
    static constexpr double minP = 10 * Acts::UnitConstants::MeV;
    nextP = std::max(minP, nextP);
    // update track parameters and covariance
    stepper.update(state.stepping, pos, dir, nextP, time);
    state.stepping.cov(eBoundPhi, eBoundPhi) = updateVariance(
        state.stepping.cov(eBoundPhi, eBoundPhi), variancePhi, updateMode);
    state.stepping.cov(eBoundTheta, eBoundTheta) =
        updateVariance(state.stepping.cov(eBoundTheta, eBoundTheta),
                       varianceTheta, updateMode);
    state.stepping.cov(eBoundQOverP, eBoundQOverP) =
        updateVariance(state.stepping.cov(eBoundQOverP, eBoundQOverP),
                       varianceQoverP, updateMode);
  }

 private:
  /// @brief Evaluates the contributions to the covariance matrix
  ///
  /// @param [in] multipleScattering Boolean to indiciate the application of
  /// multiple scattering
  /// @param [in] energyLoss Boolean to indiciate the application of energy loss
  void covarianceContributions(bool multipleScattering, bool energyLoss);

  /// @brief Convenience method for better readability
  ///
  /// @param [in] variance A diagonal entry of the covariance matrix
  /// @param [in] change The change that may be applied to it
  /// @param [in] updateMode The noise update mode (in default: add noise)
  ///
  /// @return The updated variance
  double updateVariance(double variance, double change,
                        NoiseUpdateMode updateMode = addNoise) const;
};
}  // namespace detail
}  // end of namespace Acts
