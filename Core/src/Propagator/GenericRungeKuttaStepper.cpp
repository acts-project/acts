// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/GenericRungeKuttaStepper.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Acts {

namespace {

using FieldAndGradient = GenericRungeKuttaStepper::FieldAndGradient;

/// Cross product matrix W with W * v = u x v
SquareMatrix3 crossMatrix(const Vector3& u) {
  SquareMatrix3 w;
  // clang-format off
  w <<     0., -u.z(),  u.y(),
        u.z(),     0., -u.x(),
       -u.y(),  u.x(),     0.;
  // clang-format on
  return w;
}

/// dt/ds and its derivative by q/p
struct TimeDerivative {
  double dtds = 0;
  double dtdsDqop = 0;
};

TimeDerivative timeDerivative(const ParticleHypothesis& particleHypothesis,
                              double qop) {
  const double mass = particleHypothesis.mass();
  const double absMom = particleHypothesis.extractMomentum(qop);
  const double charge = particleHypothesis.extractCharge(qop);
  TimeDerivative result;
  // time propagates along distance as 1/b = sqrt(1 + m²/p²)
  result.dtds = fastHypot(1., mass / absMom);
  // A neutral particle has a fixed momentum hypothesis, so its time does not
  // depend on q/p.
  result.dtdsDqop =
      charge != 0. ? mass * mass * qop / (charge * charge * result.dtds) : 0.;
  return result;
}

/// The equation of motion in vacuum, dy/ds = f(y)
FreeVector motion(const FreeVector& pars, const Vector3& field, double dtds) {
  const Vector3 dir = pars.segment<3>(eFreeDir0);
  FreeVector derivative = FreeVector::Zero();
  derivative.segment<3>(eFreePos0) = dir;
  derivative[eFreeTime] = dtds;
  derivative.segment<3>(eFreeDir0) = pars[eFreeQOverP] * dir.cross(field);
  return derivative;
}

/// The jacobian of the equation of motion, df/dy
FreeMatrix motionJacobian(const FreeVector& pars, const FieldAndGradient& field,
                          double dtdsDqop) {
  const Vector3 dir = pars.segment<3>(eFreeDir0);
  const double qop = pars[eFreeQOverP];
  FreeMatrix jac = FreeMatrix::Zero();
  jac.block<3, 3>(eFreePos0, eFreeDir0).setIdentity();
  jac(eFreeTime, eFreeQOverP) = dtdsDqop;
  // d(T x B)/dr = T x (dB/dr)
  jac.block<3, 3>(eFreeDir0, eFreePos0) =
      qop * crossMatrix(dir) * field.gradient;
  // d(T x B)/dT = -B x
  jac.block<3, 3>(eFreeDir0, eFreeDir0) = -qop * crossMatrix(field.field);
  jac.block<3, 1>(eFreeDir0, eFreeQOverP) = dir.cross(field.field);
  return jac;
}

}  // namespace

GenericRungeKuttaStepper::GenericRungeKuttaStepper(
    std::shared_ptr<const MagneticFieldProvider> bField)
    : GenericRungeKuttaStepper(Config{std::move(bField)}) {}

GenericRungeKuttaStepper::GenericRungeKuttaStepper(const Config& config)
    : m_bField(config.bField), m_tableau(config.tableau) {
  if (m_tableau == nullptr) {
    throw std::invalid_argument(
        "GenericRungeKuttaStepper: the tableau is not set");
  }
}

GenericRungeKuttaStepper::State GenericRungeKuttaStepper::makeState(
    const Options& options) const {
  return State{options, m_bField->makeCache(options.magFieldContext)};
}

void GenericRungeKuttaStepper::initialize(State& state,
                                          const BoundParameters& par) const {
  initialize(state, par.parameters(), par.covariance(),
             par.particleHypothesis(), par.referenceSurface());
}

void GenericRungeKuttaStepper::initialize(State& state,
                                          const BoundVector& boundParams,
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

  state.cov = cov;
  if (state.cov.has_value()) {
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.segment<3>(eFreePos0),
        freeParams.segment<3>(eFreeDir0));
    state.jacTransport = FreeMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }
}

Result<GenericRungeKuttaStepper::BoundParameters>
GenericRungeKuttaStepper::boundParameters(const State& state,
                                          const Surface& surface) const {
  return detail::boundParameters(state.options.geoContext, surface, state.pars,
                                 state.cov, state.particleHypothesis);
}

bool GenericRungeKuttaStepper::prepareCurvilinearState(State& state) const {
  // The derivatives are only missing if no step was executed yet.
  if (state.pathAccumulated != 0) {
    return true;
  }

  if (!state.field.has_value()) {
    auto fieldRes = getFieldAndGradient(state, position(state), false);
    if (!fieldRes.ok()) {
      return false;
    }
    state.field = *fieldRes;
    state.fieldHasGradient = false;
  }

  const TimeDerivative time =
      timeDerivative(state.particleHypothesis, qOverP(state));
  state.derivative = motion(state.pars, state.field->field, time.dtds);
  return true;
}

GenericRungeKuttaStepper::BoundParameters
GenericRungeKuttaStepper::curvilinearParameters(const State& state) const {
  return detail::curvilinearParameters(state.pars, state.cov,
                                       state.particleHypothesis);
}

void GenericRungeKuttaStepper::update(State& state,
                                      const FreeVector& freeParams,
                                      const BoundVector& /*boundParams*/,
                                      const Covariance& covariance,
                                      const Surface& surface) const {
  state.pars = freeParams;
  state.field.reset();
  if (state.cov.has_value()) {
    state.cov = covariance;
  }
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

void GenericRungeKuttaStepper::update(State& state, const Vector3& uposition,
                                      const Vector3& udirection, double qop,
                                      double time) const {
  // Material interactions keep the position, so the cached field stays valid.
  if (uposition != position(state)) {
    state.field.reset();
  }
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

GenericRungeKuttaStepper::Jacobian
GenericRungeKuttaStepper::transportToCurvilinear(State& state) const {
  Jacobian jacobian = Jacobian::Identity();
  if (!state.cov.has_value()) {
    return jacobian;
  }
  detail::transportCovarianceToCurvilinear(
      *state.cov, jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, std::nullopt, direction(state));
  return jacobian;
}

Result<GenericRungeKuttaStepper::Jacobian>
GenericRungeKuttaStepper::transportToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  if (!surface.isOnSurface(state.options.geoContext, position(state),
                           direction(state), BoundaryTolerance::Infinite())) {
    return Result<Jacobian>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  Jacobian jacobian = Jacobian::Identity();
  if (!state.cov.has_value()) {
    return Result<Jacobian>::success(jacobian);
  }
  detail::transportCovarianceToBound(
      state.options.geoContext, surface, *state.cov, jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, std::nullopt,
      state.pars, freeToBoundCorrection);
  return Result<Jacobian>::success(jacobian);
}

Result<FieldAndGradient> GenericRungeKuttaStepper::getFieldAndGradient(
    State& state, const Vector3& pos, bool withGradient) const {
  FieldAndGradient result;
  Result<Vector3> field = getField(state, pos);
  if (!field.ok()) {
    return Result<FieldAndGradient>::failure(field.error());
  }
  result.field = *field;
  if (!withGradient) {
    return Result<FieldAndGradient>::success(result);
  }

  // Central finite differences, with an error proportional to epsilon^2
  const double epsilon = state.options.fieldGradientEpsilon;
  for (std::size_t j = 0; j < 3; ++j) {
    Vector3 delta = Vector3::Zero();
    delta[j] = epsilon;
    Result<Vector3> plus = getField(state, pos + delta);
    if (!plus.ok()) {
      return Result<FieldAndGradient>::failure(plus.error());
    }
    Result<Vector3> minus = getField(state, pos - delta);
    if (!minus.ok()) {
      return Result<FieldAndGradient>::failure(minus.error());
    }
    result.gradient.col(j) = (*plus - *minus) / (2. * epsilon);
  }
  return Result<FieldAndGradient>::success(result);
}

Result<double> GenericRungeKuttaStepper::step(
    State& state, Direction propDir,
    const IVolumeMaterial* /*material*/) const {
  const detail::ButcherTableau& tableau = *m_tableau;
  const std::size_t nStages = tableau.stages();
  const bool withJacobian = state.cov.has_value();
  const bool withGradient = withJacobian && state.options.includeFieldGradient;
  const bool adaptive = state.options.adaptiveStepSize && tableau.hasEmbedded();

  if (!state.field.has_value() || (withGradient && !state.fieldHasGradient)) {
    auto fieldRes = getFieldAndGradient(state, position(state), withGradient);
    if (!fieldRes.ok()) {
      return fieldRes.error();
    }
    state.field = *fieldRes;
    state.fieldHasGradient = withGradient;
  }
  const FieldAndGradient startField = *state.field;
  const FreeVector start = state.pars;
  // q/p does not change in vacuum, so neither does dt/ds
  const TimeDerivative time =
      timeDerivative(state.particleHypothesis, start[eFreeQOverP]);

  const double stepTolerance = state.options.stepTolerance;
  const auto calcStepSizeScaling = [&](const double errorEstimate) -> double {
    constexpr double safety = 0.9;
    constexpr double lower = 0.2;
    constexpr double upper = 5.0;
    const double exponent = 1. / (tableau.embeddedOrder() + 1);
    const double x = safety * std::pow(stepTolerance / errorEstimate, exponent);
    return std::clamp(x, lower, upper);
  };

  const double initialH = state.stepSize.value() * propDir;
  double h = initialH;

  std::vector<FreeVector> k(nStages);
  std::vector<FreeMatrix> dk(withJacobian ? nStages : 0);
  FreeVector end;
  double errorEstimate = 0.;
  std::size_t nStepTrials = 0;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    for (std::size_t i = 0; i < nStages; ++i) {
      FreeVector stage = start;
      FreeMatrix dStage;
      if (withJacobian) {
        dStage.setIdentity();
      }
      for (std::size_t j = 0; j < i; ++j) {
        const double aij = tableau.a(i, j);
        if (aij == 0.) {
          continue;
        }
        stage += h * aij * k[j];
        if (withJacobian) {
          dStage += h * aij * dk[j];
        }
      }

      FieldAndGradient field = startField;
      if (i > 0) {
        auto fieldRes = getFieldAndGradient(state, stage.segment<3>(eFreePos0),
                                            withGradient);
        if (!fieldRes.ok()) {
          return fieldRes.error();
        }
        field = *fieldRes;
      }

      k[i] = motion(stage, field.field, time.dtds);
      if (withJacobian) {
        dk[i] = motionJacobian(stage, field, time.dtdsDqop) * dStage;
      }
    }

    end = start;
    for (std::size_t i = 0; i < nStages; ++i) {
      end += h * tableau.b(i) * k[i];
    }

    if (!adaptive) {
      break;
    }

    FreeVector diff = FreeVector::Zero();
    for (std::size_t i = 0; i < nStages; ++i) {
      diff += h * (tableau.b(i) - tableau.bEmbedded(i)) * k[i];
    }
    errorEstimate =
        std::max({diff.segment<3>(eFreePos0).cwiseAbs().maxCoeff(),
                  std::abs(diff[eFreeTime]),
                  diff.segment<3>(eFreeDir0).cwiseAbs().maxCoeff()});
    if (withJacobian && state.options.jacobianInErrorEstimate) {
      FreeMatrix jacDiff = FreeMatrix::Zero();
      for (std::size_t i = 0; i < nStages; ++i) {
        jacDiff += h * (tableau.b(i) - tableau.bEmbedded(i)) * dk[i];
      }
      errorEstimate = std::max(errorEstimate, jacDiff.cwiseAbs().maxCoeff());
    }
    errorEstimate = std::max(errorEstimate, std::numeric_limits<double>::min());

    if (errorEstimate <= stepTolerance) {
      break;
    }

    ++state.statistics.nRejectedSteps;

    h *= calcStepSizeScaling(errorEstimate);

    if (std::abs(h) < std::abs(state.options.stepSizeCutOff)) {
      return EigenStepperError::StepSizeStalled;
    }

    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  // The field at the end is the first stage of the next step. Look it up
  // before the state changes, so that a failure leaves the state unchanged.
  auto endField =
      getFieldAndGradient(state, end.segment<3>(eFreePos0), withGradient);
  if (!endField.ok()) {
    return endField.error();
  }

  const Vector3 endDir = end.segment<3>(eFreeDir0);
  const double endDirNorm = endDir.norm();

  if (withJacobian) {
    FreeMatrix stepJacobian = FreeMatrix::Identity();
    for (std::size_t i = 0; i < nStages; ++i) {
      stepJacobian += h * tableau.b(i) * dk[i];
    }
    // The direction is normalised after the step, with the derivative
    // d(T/|T|)/dT = (I - T T^T / |T|^2) / |T|.
    const Vector3 unitDir = endDir / endDirNorm;
    const SquareMatrix3 normalisation =
        (SquareMatrix3::Identity() - unitDir * unitDir.transpose()) /
        endDirNorm;
    stepJacobian.middleRows<3>(eFreeDir0) =
        (normalisation * stepJacobian.middleRows<3>(eFreeDir0)).eval();
    state.jacTransport = (stepJacobian * state.jacTransport).eval();
  }

  end.segment<3>(eFreeDir0) /= endDirNorm;
  state.pars = end;
  state.field = *endField;
  state.fieldHasGradient = withGradient;

  if (withJacobian) {
    state.derivative = motion(end, state.field->field, time.dtds);
  }

  state.pathAccumulated += h;
  ++state.nSteps;

  ++state.statistics.nSuccessfulSteps;
  if (propDir != Direction::fromScalarZeroAsPositive(initialH)) {
    ++state.statistics.nReverseSteps;
  }
  state.statistics.pathLength += h;
  state.statistics.absolutePathLength += std::abs(h);

  if (adaptive) {
    const double nextAccuracy =
        std::abs(h * calcStepSizeScaling(errorEstimate));
    const double previousAccuracy = std::abs(state.stepSize.accuracy());
    const double initialStepLength = std::abs(initialH);
    if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
      state.stepSize.setAccuracy(nextAccuracy);
    }
  }

  return h;
}

}  // namespace Acts
