// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::detail {

/// @brief Struct to handle pointwise material interaction
struct PointwiseMaterialInteraction {
  /// The current surface at the interaction.
  const Surface* surface = nullptr;
  /// The start surface of the propagation.
  const Surface* startSurface = nullptr;
  /// The target surface of the propagation.
  const Surface* targetSurface = nullptr;
  /// The path correction factor due to non-zero incidence on the surface.
  double pathCorrection = 0.;

  /// The particle position at the interaction.
  Vector3 pos = Vector3::Zero();
  /// The particle time at the interaction.
  double time = 0.0;
  /// The particle direction at the interaction.
  Vector3 dir = Vector3::Zero();
  /// The particle q/p at the interaction
  float qOverP = 0.0;
  /// The absolute particle charge
  float absQ = 0.0;
  /// The particle momentum at the interaction
  float momentum = 0.0;
  /// The particle mass
  float mass = 0.0;
  /// The particle absolute pdg
  PdgParticle absPdg = PdgParticle::eInvalid;
  /// The covariance transport decision at the interaction
  bool performCovarianceTransport = false;
  /// The propagation direction
  Direction propDir = Direction::Forward();

  /// The final material update mode
  MaterialUpdateMode updateMode = MaterialUpdateMode::NoUpdate;
  /// The effective, passed material properties including the path correction.
  MaterialSlab slab = MaterialSlab::Nothing();

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

  /// @brief Constructor
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper in use
  /// @param [in] navigator Navigator in use
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  PointwiseMaterialInteraction(const propagator_state_t& state,
                               const stepper_t& stepper,
                               const navigator_t& navigator)
      : surface(navigator.currentSurface(state.navigation)),
        startSurface(navigator.startSurface(state.navigation)),
        targetSurface(navigator.targetSurface(state.navigation)),
        pos(stepper.position(state.stepping)),
        time(stepper.time(state.stepping)),
        dir(stepper.direction(state.stepping)),
        qOverP(stepper.qOverP(state.stepping)),
        absQ(stepper.particleHypothesis(state.stepping).absoluteCharge()),
        momentum(stepper.absoluteMomentum(state.stepping)),
        mass(stepper.particleHypothesis(state.stepping).mass()),
        absPdg(stepper.particleHypothesis(state.stepping).absolutePdg()),
        performCovarianceTransport(state.stepping.covTransport),
        propDir(state.options.direction) {
    pathCorrection =
        surface->pathCorrection(state.options.geoContext, pos, dir);
  }

  /// @brief This function evaluates the material properties to interact with
  /// This updates the slab and then returns, if the resulting slab is valid
  ///
  /// @param [in] requestedMode The requested material update mode
  ///
  /// @return Boolean statement whether the material is valid
  bool evaluateMaterialSlab(MaterialUpdateMode requestedMode);

  /// @brief This function evaluate the material effects
  ///
  /// @param [in] multipleScattering Boolean to indicate the application of
  /// multiple scattering
  /// @param [in] energyLoss Boolean to indicate the application of energy loss
  void evaluatePointwiseMaterialInteraction(bool multipleScattering,
                                            bool energyLoss);

  /// @brief Update the state
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper in use
  /// @param [in] noiseUpdateMode The noise update mode
  template <typename propagator_state_t, typename stepper_t>
  void updateState(propagator_state_t& state, const stepper_t& stepper,
                   NoiseUpdateMode noiseUpdateMode) {
    const ParticleHypothesis& particleHypothesis =
        stepper.particleHypothesis(state.stepping);
    // in forward(backward) propagation, energy decreases(increases) and
    // variances increase(decrease)
    const double nextE = fastHypot(mass, momentum) - Eloss * propDir;
    // put particle at rest if energy loss is too large
    nextP = (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;
    // minimum momentum below which we will not push particles via material
    // update
    // TODO 10 MeV might be quite low and we should make this configurable
    static constexpr double minP = 10 * Acts::UnitConstants::MeV;
    nextP = std::max(minP, nextP);
    const double nextQOverP =
        std::copysign(particleHypothesis.qOverP(nextP, absQ), qOverP);
    // update track parameters and covariance
    stepper.update(state.stepping, pos, dir, nextQOverP, time);
    state.stepping.cov(eBoundPhi, eBoundPhi) = updateVariance(
        state.stepping.cov(eBoundPhi, eBoundPhi), variancePhi, noiseUpdateMode);
    state.stepping.cov(eBoundTheta, eBoundTheta) =
        updateVariance(state.stepping.cov(eBoundTheta, eBoundTheta),
                       varianceTheta, noiseUpdateMode);
    state.stepping.cov(eBoundQOverP, eBoundQOverP) =
        updateVariance(state.stepping.cov(eBoundQOverP, eBoundQOverP),
                       varianceQoverP, noiseUpdateMode);
  }

 private:
  /// @brief Evaluates the contributions to the covariance matrix
  ///
  /// @param [in] multipleScattering Boolean to indicate the application of
  /// multiple scattering
  /// @param [in] energyLoss Boolean to indicate the application of energy loss
  void evaluateCovarianceContributions(bool multipleScattering,
                                       bool energyLoss);

  /// @brief Convenience method for better readability
  ///
  /// @param [in] variance A diagonal entry of the covariance matrix
  /// @param [in] change The change that may be applied to it
  /// @param [in] noiseUpdateMode The noise update mode
  ///
  /// @return The updated variance
  static double updateVariance(double variance, double change,
                               NoiseUpdateMode noiseUpdateMode);
};

}  // namespace Acts::detail
