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
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialInteractionError.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::detail {

/// Determine the material update mode to be used for a given surface. Depending
/// on whether the surface is the start or target surface, the update mode is
/// restricted to pre- or post-update only. This is necessary to avoid double
/// counting of material effects at the start and target surfaces.
/// @param surface The current surface
/// @param startSurface The starting surface of the propagation
/// @param targetSurface The target surface of the propagation
/// @param requestedMode The requested material update mode
/// @return The determined material update mode
MaterialUpdateMode determineMaterialUpdateMode(
    const Surface& surface, const Surface* startSurface,
    const Surface* targetSurface, MaterialUpdateMode requestedMode);

/// Determine the material update mode to be used for the current surface given
/// the propagation state. Depending on whether the surface is the start or
/// target surface, the update mode is restricted to pre- or post-update only.
/// This is necessary to avoid double counting of material effects at the start
/// and target surfaces.
/// @tparam propagator_state_t The type of the propagator state
/// @tparam navigator_t The type of the navigator
/// @param state The current propagator state
/// @param navigator The navigator used for the propagation
/// @param requestedMode The requested material update mode
/// @return The determined material update mode
template <typename propagator_state_t, typename navigator_t>
MaterialUpdateMode determineMaterialUpdateMode(
    const propagator_state_t& state, const navigator_t& navigator,
    MaterialUpdateMode requestedMode) {
  const Surface* surface = navigator.currentSurface(state.navigation);
  const Surface* startSurface = navigator.startSurface(state.navigation);
  const Surface* targetSurface = navigator.targetSurface(state.navigation);

  if (surface == nullptr) {
    return MaterialUpdateMode::NoUpdate;
  }

  return determineMaterialUpdateMode(*surface, startSurface, targetSurface,
                                     requestedMode);
}

/// Evaluate the material slab at a given surface, position, and direction
/// @param geoContext The geometry context
/// @param surface The surface at which to evaluate the material slab
/// @param propagationDirection The propagation direction
/// @param position The position at which to evaluate the material slab
/// @param direction The direction at which to evaluate the material slab
/// @param updateMode The material update mode
/// @return The evaluated material slab
Result<MaterialSlab> evaluateMaterialSlab(const GeometryContext& geoContext,
                                          const Surface& surface,
                                          Direction propagationDirection,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          MaterialUpdateMode updateMode);

/// Evaluate the material slab at the propagation state and surface
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param surface The surface at which to evaluate the material slab
/// @param updateMode The material update mode
/// @return The evaluated material slab
template <typename propagator_state_t, typename stepper_t>
Result<MaterialSlab> evaluateMaterialSlab(const propagator_state_t& state,
                                          const stepper_t& stepper,
                                          const Surface& surface,
                                          MaterialUpdateMode updateMode) {
  const GeometryContext& geoContext = state.options.geoContext;
  const Direction propagationDirection = state.options.direction;
  const Vector3 position = stepper.position(state.stepping);
  const Vector3 direction = stepper.direction(state.stepping);

  return evaluateMaterialSlab(geoContext, surface, propagationDirection,
                              position, direction, updateMode);
}

/// Struct to hold the material effects computed at a pointwise interaction
struct PointwiseMaterialEffects {
  double eLoss = 0;
  /// Projected physical scattering variance, independent of bound coordinates.
  double scatteringVariance = 0;
  double varianceQoverP = 0;
};

/// Compute the material effects given a material slab and particle properties
/// @param slab The material slab
/// @param particleHypothesis The particle hypothesis
/// @param qOverP The charge over momentum of the particle
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @param covTransport Whether to compute covariance transport effects
/// @return The computed material effects
PointwiseMaterialEffects computeMaterialEffects(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    float qOverP, bool multipleScattering, bool energyLoss, bool covTransport);

/// Express isotropic scattering noise in a bound direction basis.
/// Only the noise increment is projected; the existing covariance is untouched.
/// @param parametrization Interpretation of the covariance direction coordinates
/// @param direction Global unit direction
/// @param scatteringVariance Projected physical scattering variance
/// @pre Direction is normalized; scattering variance is finite and nonnegative.
/// @return Direction noise block, or a singular-direction error
inline Result<Matrix<2, 2>> scatteringNoise(
    const BoundParametrization& parametrization, const Vector3& direction,
    double scatteringVariance) {
  assert(direction.allFinite() &&
         std::abs(direction.squaredNorm() - 1.) < 1e-6);
  assert(std::isfinite(scatteringVariance) && scatteringVariance >= 0.);
  if (scatteringVariance == 0.) {
    return Result<Matrix<2, 2>>::success(Matrix<2, 2>::Zero());
  }
  if (parametrization.isGlobalSpherical()) {
    const double transverseSquared = direction.head<2>().squaredNorm();
    if (transverseSquared == 0.) {
      return Result<Matrix<2, 2>>::failure(
          BoundParametrizationError::SingularDirection);
    }
    const double variance0 = scatteringVariance / transverseSquared;
    if (!std::isfinite(variance0)) {
      return Result<Matrix<2, 2>>::failure(
          BoundParametrizationError::SingularDirection);
    }
    const Vector2 variances(variance0, scatteringVariance);
    return Result<Matrix<2, 2>>::success(Matrix<2, 2>(variances.asDiagonal()));
  }
  const auto jacobian =
      parametrization.directionCoordinatesJacobianFromUnitDirection(direction,
                                                                    0.);
  if (!jacobian.ok()) {
    return Result<Matrix<2, 2>>::failure(jacobian.error());
  }
  // The encoding derivative already includes the tangent projection I - d*d^T.
  return Result<Matrix<2, 2>>::success(scatteringVariance * (*jacobian) *
                                       jacobian->transpose());
}

/// Internal checked subtraction; see updateMaterialCovariance for the contract.
/// @param covariance Covariance, modified only on success
/// @param directionNoise Direction noise to subtract
/// @param varianceQoverP Q/p variance to subtract
/// @return Success or InvalidCovarianceAfterNoiseRemoval
Result<void> removeMaterialCovarianceNoise(BoundMatrix& covariance,
                                           const Matrix<2, 2>& directionNoise,
                                           double varianceQoverP);

/// Add or remove material noise in the current bound covariance basis.
/// Removal checks the whole covariance, including correlations outside the
/// direction block. It does not clamp variances or partially remove noise.
/// @param covariance Covariance, modified only on success
/// @param directionNoise Symmetric positive-semidefinite direction noise block
/// @param varianceQoverP Nonnegative q/p variance increment
/// @param mode Whether to add or remove noise
/// @pre Inputs are finite; covariance is symmetric and positive semidefinite.
/// @return Success or InvalidCovarianceAfterNoiseRemoval. The removal check
/// uses diagonal scaling and a relative eigenvalue tolerance of 64 machine
/// epsilons.
inline Result<void> updateMaterialCovariance(BoundMatrix& covariance,
                                             const Matrix<2, 2>& directionNoise,
                                             double varianceQoverP,
                                             NoiseUpdateMode mode) {
  assert(covariance.allFinite() && directionNoise.allFinite());
  assert(std::isfinite(varianceQoverP) && varianceQoverP >= 0.);
  if (mode == NoiseUpdateMode::addNoise) {
    covariance.block<2, 2>(eBoundDir0, eBoundDir0) += directionNoise;
    covariance(eBoundQOverP, eBoundQOverP) += varianceQoverP;
    return Result<void>::success();
  }
  return removeMaterialCovarianceNoise(covariance, directionNoise,
                                       varianceQoverP);
}

/// Global angular variances for legacy material diagnostics.
/// @param direction Global unit direction
/// @param scatteringVariance Projected physical scattering variance
/// @pre Direction is normalized; scattering variance is finite and nonnegative.
/// @return Phi/theta variances; phi is infinite at a pole for nonzero scattering
Vector2 globalAngularScatteringVariances(const Vector3& direction,
                                         double scatteringVariance);

/// Compute the material effects given the propagation state and material slab
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param slab The material slab
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @return The computed material effects
template <typename propagator_state_t, typename stepper_t>
PointwiseMaterialEffects computeMaterialEffects(const propagator_state_t& state,
                                                const stepper_t& stepper,
                                                const MaterialSlab& slab,
                                                bool multipleScattering,
                                                bool energyLoss) {
  const bool covTransport = stepper.hasCovariance(state.stepping);
  const float qOverP = stepper.qOverP(state.stepping);
  const ParticleHypothesis& particleHypothesis =
      stepper.particleHypothesis(state.stepping);

  return computeMaterialEffects(slab, particleHypothesis, qOverP,
                                multipleScattering, energyLoss, covTransport);
}

/// Perform the material interaction given the propagation state and material
/// slab
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param slab The material slab
/// @param parametrization Interpretation of the active covariance anchor
/// @param noiseUpdateMode The noise update mode
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @return The computed effects, or a direction-noise mapping error
template <typename propagator_state_t, typename stepper_t>
Result<PointwiseMaterialEffects> performMaterialInteraction(
    propagator_state_t& state, const stepper_t& stepper,
    const MaterialSlab& slab, const BoundParametrization& parametrization,
    NoiseUpdateMode noiseUpdateMode, bool multipleScattering, bool energyLoss) {
  if (slab.isVacuum()) {
    return PointwiseMaterialEffects{};
  }

  const PointwiseMaterialEffects effects = computeMaterialEffects(
      state, stepper, slab, multipleScattering, energyLoss);

  const auto noise =
      scatteringNoise(parametrization, stepper.direction(state.stepping),
                      effects.scatteringVariance);
  if (!noise.ok()) {
    return noise.error();
  }

  const Direction propDir = state.options.direction;
  const ParticleHypothesis& particleHypothesis =
      stepper.particleHypothesis(state.stepping);
  const double mass = particleHypothesis.mass();
  const double absQ = particleHypothesis.absoluteCharge();
  const Vector3 position = stepper.position(state.stepping);
  const double time = stepper.time(state.stepping);
  const Vector3 direction = stepper.direction(state.stepping);
  const float qOverP = stepper.qOverP(state.stepping);
  const double momentum = stepper.absoluteMomentum(state.stepping);

  //! [energy loss update]
  // in forward(backward) propagation, energy decreases(increases) and
  // variances increase(decrease)
  const double nextE = fastHypot(mass, momentum) - effects.eLoss * propDir;
  // put particle at rest if energy loss is too large
  double nextP = (mass < nextE) ? fastCathetus(nextE, mass) : 0;

  // minimum momentum below which we will not push particles via material
  // update
  // TODO 10 MeV might be quite low and we should make this configurable
  static constexpr double minP = 10 * Acts::UnitConstants::MeV;
  nextP = std::max(minP, nextP);
  const double nextQOverP =
      particleHypothesis.qOverP(nextP, std::copysign(absQ, qOverP));

  // Validate noise removal before committing either momentum or covariance.
  auto covariance = stepper.covariance(state.stepping);
  if (covariance.has_value()) {
    const auto updateResult = updateMaterialCovariance(
        *covariance, *noise, effects.varianceQoverP, noiseUpdateMode);
    if (!updateResult.ok()) {
      return updateResult.error();
    }
  }

  // update track parameters
  stepper.update(state.stepping, position, direction, nextQOverP, time);
  //! [energy loss update]

  //! [covariance update]
  if (covariance.has_value()) {
    stepper.setCovariance(state.stepping, *covariance);
  }
  //! [covariance update]

  return effects;
}

/// Perform the material interaction at the current surface given the
/// propagation state
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param surface The surface at which to perform the material interaction
/// @param updateMode The material update mode
/// @param parametrization Interpretation of the active covariance anchor
/// @param noiseUpdateMode The noise update mode
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @param logger The logger to use for verbose output
/// @return The computed effects, or a direction-noise mapping error
template <typename propagator_state_t, typename stepper_t>
Result<PointwiseMaterialEffects> performMaterialInteraction(
    propagator_state_t& state, const stepper_t& stepper, const Surface& surface,
    const BoundParametrization& parametrization, MaterialUpdateMode updateMode,
    NoiseUpdateMode noiseUpdateMode, bool multipleScattering, bool energyLoss,
    const Logger& logger) {
  const Result<MaterialSlab> slabResult =
      evaluateMaterialSlab(state, stepper, surface, updateMode);
  if (!slabResult.ok()) {
    return slabResult.error();
  }
  const MaterialSlab& slab = slabResult.value();
  if (slab.isVacuum()) {
    ACTS_VERBOSE("No material effects on surface: " << surface.geometryId()
                                                    << " with update mode: "
                                                    << updateMode);
  }

  const auto effectsResult = performMaterialInteraction(
      state, stepper, slab, parametrization, noiseUpdateMode,
      multipleScattering, energyLoss);
  if (!effectsResult.ok()) {
    return effectsResult.error();
  }
  const PointwiseMaterialEffects& effects = *effectsResult;

  const Direction propDir = state.options.direction;

  ACTS_VERBOSE("Material effects on surface: "
               << surface.geometryId() << " with update mode: " << updateMode
               << " are :");
  ACTS_VERBOSE("eLoss = " << effects.eLoss * propDir << ", "
                          << "scatteringVariance = "
                          << effects.scatteringVariance << ", "
                          << "varianceQoverP = " << effects.varianceQoverP);

  return effects;
}

}  // namespace Acts::detail
