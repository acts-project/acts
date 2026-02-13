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
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"

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
MaterialSlab evaluateMaterialSlab(const GeometryContext& geoContext,
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
MaterialSlab evaluateMaterialSlab(const propagator_state_t& state,
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

/// Evaluate the material slab at the current surface given the propagation
/// state
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @tparam navigator_t The type of the navigator
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param navigator The navigator used for the propagation
/// @param updateMode The material update mode
/// @return The evaluated material slab
template <typename propagator_state_t, typename stepper_t, typename navigator_t>
MaterialSlab evaluateMaterialSlab(const propagator_state_t& state,
                                  const stepper_t& stepper,
                                  const navigator_t& navigator,
                                  MaterialUpdateMode updateMode) {
  const Surface* surface = navigator.currentSurface(state.navigation);
  if (surface == nullptr) {
    return MaterialSlab::Nothing();
  }

  const GeometryContext& geoContext = state.options.geoContext;
  const Direction propDir = state.options.direction;
  const Vector3 pos = stepper.position(state.stepping);
  const Vector3 dir = stepper.direction(state.stepping);

  return evaluateMaterialSlab(geoContext, *surface, propDir, pos, dir,
                              updateMode);
}

/// Struct to hold the material effects computed at a pointwise interaction
struct PointwiseMaterialEffects {
  double eLoss = 0;
  double variancePhi = 0;
  double varianceTheta = 0;
  double varianceQoverP = 0;
};

/// Compute the material effects given a material slab and particle properties
/// @param slab The material slab
/// @param particleHypothesis The particle hypothesis
/// @param direction The direction of the particle
/// @param qOverP The charge over momentum of the particle
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @param covTransport Whether to compute covariance transport effects
/// @return The computed material effects
PointwiseMaterialEffects computeMaterialEffects(
    const MaterialSlab& slab, const ParticleHypothesis& particleHypothesis,
    const Vector3& direction, float qOverP, bool multipleScattering,
    bool energyLoss, bool covTransport);

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
  const bool covTransport = state.stepping.covTransport;
  const Vector3 direction = stepper.direction(state.stepping);
  const float qOverP = stepper.qOverP(state.stepping);
  const ParticleHypothesis& particleHypothesis =
      stepper.particleHypothesis(state.stepping);

  return computeMaterialEffects(slab, particleHypothesis, direction, qOverP,
                                multipleScattering, energyLoss, covTransport);
}

/// Perform the material interaction given the propagation state and material
/// slab
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param slab The material slab
/// @param noiseUpdateMode The noise update mode
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @return The computed material effects
template <typename propagator_state_t, typename stepper_t>
PointwiseMaterialEffects performMaterialInteraction(
    propagator_state_t& state, const stepper_t& stepper,
    const MaterialSlab& slab, NoiseUpdateMode noiseUpdateMode,
    bool multipleScattering, bool energyLoss) {
  if (slab.isVacuum()) {
    return {};
  }

  const PointwiseMaterialEffects effects = computeMaterialEffects(
      state, stepper, slab, multipleScattering, energyLoss);

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

  // in forward(backward) propagation, energy decreases(increases) and
  // variances increase(decrease)
  const double nextE = fastHypot(mass, momentum) - effects.eLoss * propDir;
  // put particle at rest if energy loss is too large
  double nextP = (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;

  // minimum momentum below which we will not push particles via material
  // update
  // TODO 10 MeV might be quite low and we should make this configurable
  static constexpr double minP = 10 * Acts::UnitConstants::MeV;
  nextP = std::max(minP, nextP);
  const double nextQOverP =
      particleHypothesis.qOverP(nextP, std::copysign(absQ, qOverP));

  // update track parameters
  stepper.update(state.stepping, position, direction, nextQOverP, time);

  // Convenience method to update a variance given a change and noise update
  // mode
  const auto updateVariance = [](double variance, double change,
                                 NoiseUpdateMode updateMode) {
    // Add/Subtract the change
    // Protect the variance against becoming negative
    return std::max(0.,
                    variance + std::copysign(change, toUnderlying(updateMode)));
  };

  // update covariance matrix
  state.stepping.cov(eBoundPhi, eBoundPhi) =
      updateVariance(state.stepping.cov(eBoundPhi, eBoundPhi),
                     effects.variancePhi, noiseUpdateMode);
  state.stepping.cov(eBoundTheta, eBoundTheta) =
      updateVariance(state.stepping.cov(eBoundTheta, eBoundTheta),
                     effects.varianceTheta, noiseUpdateMode);
  state.stepping.cov(eBoundQOverP, eBoundQOverP) =
      updateVariance(state.stepping.cov(eBoundQOverP, eBoundQOverP),
                     effects.varianceQoverP, noiseUpdateMode);

  return effects;
}

/// Perform the material interaction at the current surface given the
/// propagation state
/// @tparam propagator_state_t The type of the propagator state
/// @tparam stepper_t The type of the stepper
/// @tparam navigator_t The type of the navigator
/// @param state The current propagator state
/// @param stepper The stepper used for the propagation
/// @param navigator The navigator used for the propagation
/// @param updateMode The material update mode
/// @param noiseUpdateMode The noise update mode
/// @param multipleScattering Whether to compute multiple scattering effects
/// @param energyLoss Whether to compute energy loss effects
/// @param logger The logger to use for verbose output
/// @return The computed material effects
template <typename propagator_state_t, typename stepper_t, typename navigator_t>
PointwiseMaterialEffects performMaterialInteraction(
    propagator_state_t& state, const stepper_t& stepper,
    const navigator_t& navigator, MaterialUpdateMode updateMode,
    NoiseUpdateMode noiseUpdateMode, bool multipleScattering, bool energyLoss,
    const Logger& logger) {
  const Surface* surface = navigator.currentSurface(state.navigation);
  if (surface == nullptr) {
    ACTS_VERBOSE("No current surface, skip material interaction");
    return {};
  }

  const MaterialSlab slab =
      evaluateMaterialSlab(state, stepper, navigator, updateMode);
  if (slab.isVacuum()) {
    ACTS_VERBOSE("No material effects on surface: " << surface->geometryId()
                                                    << " with update mode: "
                                                    << updateMode);
  }

  const PointwiseMaterialEffects effects = performMaterialInteraction(
      state, stepper, slab, noiseUpdateMode, multipleScattering, energyLoss);

  const Direction propDir = state.options.direction;

  ACTS_VERBOSE("Material effects on surface: "
               << surface->geometryId() << " with update mode: " << updateMode
               << " are :");
  ACTS_VERBOSE("eLoss = " << effects.eLoss * propDir << ", "
                          << "variancePhi = " << effects.variancePhi << ", "
                          << "varianceTheta = " << effects.varianceTheta << ", "
                          << "varianceQoverP = " << effects.varianceQoverP);

  return effects;
}

}  // namespace Acts::detail
