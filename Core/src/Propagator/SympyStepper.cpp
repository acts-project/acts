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
#include "Acts/Propagator/detail/SympyBoundToFreeScaling.hpp"
#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"
#include "Acts/Propagator/detail/SympyStepperDenseStep.hpp"
#include "Acts/Propagator/detail/SympyStepperVacuumStepDecl.hpp"

#include <cmath>
#include <span>

namespace Acts {

namespace {
/// dt/ds from the current q/p.  The vacuum kernel is handed this instead of
/// forming it, so it has to be refreshed wherever q/p moves.
double computeDtds(const SympyStepper::State& state) {
  const double m = state.particleHypothesis.mass();
  const double p =
      state.particleHypothesis.extractMomentum(state.pars[eFreeQOverP]);
  return std::sqrt(1 + m * m / (p * p));
}
}  // namespace

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
  state.nStepTrials = 0;
  state.stepSize = ConstrainedStep();
  state.stepSize.setAccuracy(state.options.initialStepSize);
  state.stepSize.setUser(state.options.maxStepSize);
  state.previousStepSize = 0;
  state.statistics = StepperStatistics();

  state.pars = freeParams;
  state.field.reset();
  // after state.pars, which it reads
  state.dtds = computeDtds(state);

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
    state.jacobian = BoundMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }
}

Result<std::tuple<SympyStepper::BoundParameters, BoundMatrix, double>>
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  const bool transport = state.covTransport && transportCov;
  std::optional<FreeMatrix> additionalFreeCovariance;
  if (transport) {
    additionalFreeCovariance =
        state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
            direction(state));
  }
  state.materialEffectsAccumulator.reset();
  // The engine reads the jacobian and reinitializes it for the new surface, so
  // convert to plain before the call and back to scaled after. Without
  // transport the engine does not touch the jacobian.
  if (transport) {
    detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  }
  auto result = detail::sympy::boundState(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.derivative, state.jacToGlobal, additionalFreeCovariance, state.pars,
      state.particleHypothesis, transport, state.pathAccumulated,
      freeToBoundCorrection);
  if (transport) {
    detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
  }
  return result;
}

bool SympyStepper::prepareCurvilinearState(State& state) const {
  // TODO implement like in EigenStepper
  static_cast<void>(state);
  return true;
}

std::tuple<SympyStepper::BoundParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  const bool transport = state.covTransport && transportCov;
  std::optional<FreeMatrix> additionalFreeCovariance;
  if (transport) {
    additionalFreeCovariance =
        state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
            direction(state));
  }
  state.materialEffectsAccumulator.reset();
  if (transport) {
    detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  }
  auto result = detail::sympy::curvilinearState(
      state.cov, state.jacobian, state.derivative, state.jacToGlobal,
      additionalFreeCovariance, state.pars, state.particleHypothesis, transport,
      state.pathAccumulated);
  if (transport) {
    detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
  }
  return result;
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.field.reset();
  state.dtds = computeDtds(state);
  state.cov = covariance;
  if (state.covTransport) {
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.template segment<3>(eFreePos0),
        freeParams.template segment<3>(eFreeDir0));
    detail::sympy::toScaledBoundToFree(state.jacToGlobal,
                                       freeParams[eFreeQOverP]);
  }
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
  state.dtds = computeDtds(state);
  state.field.reset();
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  if (!state.covTransport) {
    return;
  }
  detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.derivative, state.jacToGlobal,
      std::nullopt, state.pars.template segment<3>(eFreeDir0));
  detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  if (!state.covTransport) {
    return;
  }
  detail::sympy::fromScaledBoundToFree(state.jacToGlobal, qOverP(state));
  detail::sympy::transportCovarianceToBound(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.derivative, state.jacToGlobal, std::nullopt, state.pars,
      freeToBoundCorrection);
  detail::sympy::toScaledBoundToFree(state.jacToGlobal, qOverP(state));
}

Result<double> SympyStepper::step(State& state, Direction propDir,
                                  const IVolumeMaterial* material) const {
  // Everything about material lives in the other translation unit, including
  // the branch: what stayed behind here used to share a stack frame and a
  // register allocation with the vacuum path, which is the one that runs
  // almost everywhere.
  if (state.options.doDense &&
      (material != nullptr || !state.materialEffectsAccumulator.isVacuum())) {
    if (state.covTransport) {
      return detail::sympyDenseStepFull<true>(*this, state, propDir, material);
    }
    return detail::sympyDenseStepFull<false>(*this, state, propDir, material);
  }

  // Specialised on covariance transport, each specialisation in a translation
  // unit of its own; see SympyStepperVacuumStepDecl.hpp for why.
  if (state.covTransport) {
    return detail::sympyVacuumStep<true>(*this, state, propDir);
  }
  return detail::sympyVacuumStep<false>(*this, state, propDir);
}

}  // namespace Acts
