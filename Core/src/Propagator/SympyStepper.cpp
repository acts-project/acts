// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyBoundToFreeScaling.hpp"
#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"

#include <cmath>
#include <span>

#include "detail/SympyStepperStep.hpp"

namespace Acts {

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

SympyStepper::SympyStepper(const Config& config) : m_bField(config.bField) {}

SympyStepper::State SympyStepper::makeState(const Options& options) const {
  State state{options, m_bField->makeCache(options.magFieldContext)};
  return state;
}

void SympyStepper::initialize(State& state, const BoundParameters& par) const {
  return initialize(state, par.parameters(), par.covariance(),
                    par.particleHypothesis(), par.referenceSurface());
}

void SympyStepper::initialize(State& state, const BoundVector& boundParams,
                              const std::optional<BoundMatrix>& cov,
                              ParticleHypothesis particleHypothesis,
                              const Surface& surface) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  state.particleHypothesis = particleHypothesis;

  state.pathAccumulated = 0;
  state.nSteps = 0;
  state.stepSize = ConstrainedStep();
  state.stepSize.setAccuracy(state.options.initialStepSize);
  state.stepSize.setUser(state.options.maxStepSize);
  state.previousStepSize = 0;
  state.statistics = StepperStatistics();

  state.pars = freeParams;
  state.field.reset();
  state.dtds = detail::sympyDtds(state);

  // Init the jacobian matrix if needed
  state.covTransport = cov.has_value();
  if (state.covTransport) {
    // set the covariance transport flag to true and copy
    state.cov = *cov;
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.segment<3>(eFreePos0),
        freeParams.segment<3>(eFreeDir0));
    detail::sympy::toScaledBoundToFree(state.jacToGlobal,
                                       freeParams[eFreeQOverP]);
    state.derivative = FreeVector::Zero();
  }
}

Result<SympyStepper::BoundParameters> SympyStepper::boundParameters(
    const State& state, const Surface& surface) const {
  return detail::boundParameters(
      state.options.geoContext, surface, state.pars,
      state.covTransport ? std::optional(state.cov) : std::nullopt,
      state.particleHypothesis);
}

bool SympyStepper::prepareCurvilinearState(State& state) const {
  // TODO implement like in EigenStepper
  static_cast<void>(state);
  return true;
}

SympyStepper::BoundParameters SympyStepper::curvilinearParameters(
    const State& state) const {
  return detail::curvilinearParameters(
      state.pars, state.covTransport ? std::optional(state.cov) : std::nullopt,
      state.particleHypothesis);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.field.reset();
  state.dtds = detail::sympyDtds(state);
  state.cov = covariance;
  if (state.covTransport) {
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.template segment<3>(eFreePos0),
        freeParams.template segment<3>(eFreeDir0));
    detail::sympy::toScaledBoundToFree(state.jacToGlobal,
                                       freeParams[eFreeQOverP]);
    state.derivative = FreeVector::Zero();
  }
  state.materialEffectsAccumulator.reset();
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qOverP,
                          double time) const {
  if (state.covTransport) {
    detail::sympy::rescaleBoundToFree(state.jacToGlobal,
                                      state.pars[eFreeQOverP], qOverP);
  }
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qOverP;
  state.dtds = detail::sympyDtds(state);
  state.field.reset();
}

SympyStepper::Jacobian SympyStepper::transportToCurvilinear(
    State& state) const {
  Jacobian jacobian = Jacobian::Identity();
  if (!state.covTransport) {
    state.materialEffectsAccumulator.reset();
    return jacobian;
  }
  const std::optional<FreeMatrix> additionalFreeCovariance =
      state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
          direction(state));
  state.materialEffectsAccumulator.reset();
  detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, jacobian, state.derivative, state.jacToGlobal,
      additionalFreeCovariance, state.pars.template segment<3>(eFreeDir0));
  detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
  return jacobian;
}

Result<SympyStepper::Jacobian> SympyStepper::transportToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  if (!surface.isOnSurface(state.options.geoContext, position(state),
                           direction(state), BoundaryTolerance::Infinite())) {
    return Result<Jacobian>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  Jacobian jacobian = Jacobian::Identity();
  if (!state.covTransport) {
    state.materialEffectsAccumulator.reset();
    return Result<Jacobian>::success(jacobian);
  }
  const std::optional<FreeMatrix> additionalFreeCovariance =
      state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
          direction(state));
  state.materialEffectsAccumulator.reset();
  detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  detail::sympy::transportCovarianceToBound(
      state.options.geoContext, surface, state.cov, jacobian, state.derivative,
      state.jacToGlobal, additionalFreeCovariance, state.pars,
      freeToBoundCorrection);
  detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
  return Result<Jacobian>::success(jacobian);
}

Result<double> SympyStepper::step(State& state, Direction propDir,
                                  const IVolumeMaterial* material) const {
  if (state.options.doDense &&
      (material != nullptr || !state.materialEffectsAccumulator.isVacuum())) {
    if (state.covTransport) {
      return detail::sympyStep<detail::SympyStepMode::Dense, true>(
          *this, state, propDir, material);
    }
    return detail::sympyStep<detail::SympyStepMode::Dense, false>(
        *this, state, propDir, material);
  }
  if (state.covTransport) {
    return detail::sympyStep<detail::SympyStepMode::Vacuum, true>(
        *this, state, propDir, nullptr);
  }
  return detail::sympyStep<detail::SympyStepMode::Vacuum, false>(
      *this, state, propDir, nullptr);
}

}  // namespace Acts
